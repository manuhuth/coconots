# include <numeric>
# include <algorithm>
# include <vector>
# include <cmath>
# include <unordered_map>
# include <limits>
# include <utility>
# include <atomic>
# include <thread>
# include <mutex>
# include <string>
# include <stdexcept>
# include <Rcpp.h>
# include "coco_ad.h"

// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp ;

namespace {

// Factorials up to 170! are exactly representable in double; 171! overflows.
// Counts above 170 are outside the supported range of the model (all
// densities are built from these factorials), so a factorial of an argument
// above 170 raises an error instead of returning Inf. The error is thrown as
// a plain C++ exception (no R API call), so it is safe inside the worker
// threads of parallel_for, which re-raise its message on the main thread.
const int MAX_COUNT = 170;

struct FactTable {
  static const int size = MAX_COUNT + 1;
  double v[size];
  FactTable() {
    v[0] = 1.0;
    v[1] = 1.0;
    for (int i = 2; i < size; ++i) v[i] = i * v[i - 1];
  }
};
const FactTable FACT_TABLE;

[[noreturn]] void count_limit_error(int i) {
  throw std::range_error(
    "coconots: counts above " + std::to_string(MAX_COUNT) +
    " are not supported (factorial overflow in the C++ backend); required count = " +
    std::to_string(i));
}

inline double factorial(int i) {
  if (i < 2) return 1.0;
  if (i <= MAX_COUNT) return FACT_TABLE.v[i];
  count_limit_error(i);
}

enum LinkFunction { LINK_LOG, LINK_IDENTITY, LINK_RELU, LINK_SOFTPLUS };

LinkFunction parseLinkFunction(const std::string& link_function) {
  if (link_function == "log") return LINK_LOG;
  if (link_function == "identity") return LINK_IDENTITY;
  if (link_function == "relu") return LINK_RELU;
  if (link_function == "softplus") return LINK_SOFTPLUS;
  stop("Unknown link function");
}
// The likelihood building blocks below are templates on the scalar type S:
// S = double is the plain evaluation (the same expressions as before the
// templating, hence the same arithmetic), S = coco_ad::Dual<...> evaluates
// derivatives by forward-mode automatic differentiation (see coco_ad.h).
// Branches on S compare values; the derivative is that of the branch taken.

template <class S>
inline S applyLink(const S& mathelp, LinkFunction link) {
  switch (link) {
  case LINK_LOG:      return exp(mathelp);
  case LINK_IDENTITY: return mathelp;
  case LINK_RELU:     if (mathelp > 0) return mathelp;
                      return S(1e-10);
  default:            if (mathelp > 0) return mathelp + log(1.0 + exp(-mathelp));
                      return log(1.0 + exp(mathelp));
  }
}

// Accumulates in the same order as std::inner_product over a copied row,
// but without allocating the row.
template <class S, class V>
inline S rowInnerProductT(const NumericMatrix& xreg, int row, const V& lambdas) {
  S acc = 0.0;
  const int n = xreg.ncol();
  for (int j = 0; j < n; ++j) acc = acc + xreg(row, j) * lambdas[j];
  return acc;
}

inline double rowInnerProduct(const NumericMatrix& xreg, int row,
                              const NumericVector& lambdas) {
  return rowInnerProductT<double>(xreg, row, lambdas);
}

// Largest admissible lookup-table index. Count data anywhere near this bound
// makes the model degenerate (factorials overflow at 171), so refusing is
// safer than attempting the corresponding gigabyte-sized allocations.
const int MAX_TABLE_INDEX = 1000000;

inline void check_table_index(int n) {
  if (n > MAX_TABLE_INDEX)
    stop("Count values are too large for the (Generalized) Poisson recursions.");
}

// Simulation needs a finite, positive innovation mean lambda_t = link(x_t' beta)
// (the identity link gives lambda_t <= 0 for x_t' beta <= 0).
inline void check_sim_lambda(double lam) {
  if (!(lam > 0) || !R_FINITE(lam))
    stop("The innovation mean lambda_t = link(x_t' beta) must be finite and positive "
         "for every t, but lambda_t = %g occurs. With the identity link choose "
         "covariates and beta such that x_t' beta > 0.", lam);
}

inline int resolve_threads(int n, int cores) {
  int nthreads = std::max(1, cores);
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc > 0) nthreads = std::min(nthreads, static_cast<int>(hc));
  return std::max(1, std::min(nthreads, n));
}

// Runs f(thread_id, i) for i = 0..n-1 on `nthreads` threads with dynamic
// scheduling (work per item varies widely across (y, z) histories).
// Workers must not touch the R API; any C++ exception is captured and
// re-raised as an R error on the main thread after all workers joined.
template <class F>
void parallel_for(int n, int nthreads, F f) {
  if (nthreads <= 1 || n <= 1) {
    for (int i = 0; i < n; ++i) f(0, i);
    return;
  }
  std::atomic<int> next(0);
  std::atomic<bool> failed(false);
  std::mutex msg_mutex;
  std::string msg;
  auto work = [&](int tid) {
    try {
      int i;
      while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n) {
        if (failed.load(std::memory_order_relaxed)) return;
        f(tid, i);
      }
    } catch (const std::exception& e) {
      std::lock_guard<std::mutex> lock(msg_mutex);
      if (msg.empty()) msg = e.what();
      failed.store(true);
    } catch (...) {
      failed.store(true);
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nthreads - 1);
  for (int k = 1; k < nthreads; ++k) pool.emplace_back(work, k);
  work(0);
  for (std::size_t k = 0; k < pool.size(); ++k) pool[k].join();
  if (failed.load()) {
    if (!msg.empty()) stop(msg);
    stop("coconots: a worker thread failed (out of memory?)");
  }
}

// True iff x and (for dual numbers) all its partials are exactly zero; for
// eta this identifies the Poisson model in every AD instantiation.
inline bool is_exact_zero(double x) { return x == 0.0; }

template <class T, int N>
inline bool is_exact_zero(const coco_ad::Dual<T, N>& x) {
  if (!is_exact_zero(x.v)) return false;
  for (int i = 0; i < N; ++i) {
    if (!is_exact_zero(x.d[i])) return false;
  }
  return true;
}

// True iff the value of x (for dual numbers: the innermost primal value) is
// exactly zero, whatever its partials.
inline bool is_zero_value(double x) { return x == 0.0; }

template <class T, int N>
inline bool is_zero_value(const coco_ad::Dual<T, N>& x) {
  return is_zero_value(x.v);
}

// ------------------------------------------------------------------ order 1
//
// Quasi-binomial thinning factors at alpha on a boundary of [0, 1].
// g(r | y) = y! / (y-r)! / r! * alpha (1-alpha) pow(alpha + psi r, r-1)
//            * pow(1 - alpha + psi (y-r), y-r-1) / pow(1 + psi y, y-1).
// At alpha = 0 the r = 0 factor alpha * pow(alpha, -1) is 0 * Inf = NaN; at
// alpha = 1 (psi = 0) the r = y factor (1-alpha) * pow(1-alpha, -1) is
// 0 * Inf = NaN. Both factors equal 1 identically in the parameters for
// alpha in (0, 1) (psi * 0 is exactly 0), so their value at the boundary is
// the limit 1 with zero partials, which dual numbers then carry exactly.
// The entry points branch on gp1_alpha_on_bound (exact-zero value of alpha
// or 1 - alpha) and keep their former arithmetic otherwise.
template <class S>
inline bool gp1_alpha_on_bound(const S& alpha) {
  return is_zero_value(alpha) || is_zero_value(1 - alpha);
}

// alpha * pow(alpha + psi r, r-1), the limit 1 at r = 0 and alpha = 0.
template <class S>
inline S gp1_thin_a(const S& alpha, const S& psi, int r) {
  if (r == 0 && is_zero_value(alpha)) return S(1.0);
  return alpha * pow(alpha + psi * r, r - 1);
}

// (1-alpha) * pow(1 - alpha + psi k, k-1), the limit 1 at k = 0 and alpha = 1.
template <class S>
inline S gp1_thin_b(const S& alpha, const S& psi, int k) {
  if (k == 0 && is_zero_value(1 - alpha)) return S(1.0);
  return (1 - alpha) * pow(1 - alpha + psi * k, k - 1);
}

// y! / (y-r)! / r! * alpha (1-alpha) pow(alpha + psi r, r-1)
// * pow(1 - alpha + psi (y-r), y-r-1) for alpha on a boundary (g(r | y)
// without the normalizer pow(1 + psi y, y-1)).
template <class S>
inline S gp1_thin_bound(int r, int y, const S& alpha, const S& psi) {
  return factorial(y) / factorial(y - r) / factorial(r) * gp1_thin_a(alpha, psi, r) *
    gp1_thin_b(alpha, psi, y - r);
}

// P(x | y) = sum_r  QuasiBinomial(r; y) * GP(x - r; lambda, eta)
//
// The r-invariant factors pow(alpha+psi*r, r-1), pow(1-alpha+psi*k, k-1),
// pow(1+psi*y, y-1) and the GP innovation terms are tabulated once per
// parameter value; each table entry equals the pow()/exp() call it replaces.

template <class S>
struct GP1TablesT {
  S lambda, alpha, eta, psi;
  std::vector<S> powA, powB, denom, powL, expL;
  int cap;

  void init(const S& lambda_, const S& alpha_, const S& eta_, const S& psi_) {
    lambda = lambda_; alpha = alpha_; eta = eta_; psi = psi_;
    powA.clear(); powB.clear(); denom.clear(); powL.clear(); expL.clear();
    cap = -1;
  }

  void ensure(int n) {
    if (n <= cap) return;
    check_table_index(n);
    powA.resize(n + 1); powB.resize(n + 1); denom.resize(n + 1);
    powL.resize(n + 1); expL.resize(n + 1);
    for (int i = cap + 1; i <= n; ++i) {
      powA[i]  = pow(alpha + psi * i, i - 1);
      powB[i]  = pow(1 - alpha + psi * i, i - 1);
      denom[i] = pow(1 + psi * i, i - 1);
      powL[i]  = pow(lambda + i * eta, i - 1);
      expL[i]  = exp(-lambda - i * eta);
    }
    cap = n;
  }
};
typedef GP1TablesT<double> GP1Tables;

// gp1_density for alpha on a boundary (gp1_thin_bound instead of the
// tabulated thinning factors, whose index-0 entries are pow(0, -1) there).
template <class S>
inline S gp1_density_bound(int x, int y, const GP1TablesT<S>& tb) {
  const int minxy = std::min(x, y);
  S nain = 0.0;
  for (int r = 0; r <= minxy; ++r) {
    nain = nain + gp1_thin_bound(r, y, tb.alpha, tb.psi) / tb.denom[y] *
      tb.lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
  }
  return nain;
}

// Transition density P(x | y); term order and operation chain match the
// original inline expression exactly. With alpha on a boundary,
// gp1_density_bound.
template <class S>
inline S gp1_density(int x, int y, const GP1TablesT<S>& tb) {
  if (gp1_alpha_on_bound(tb.alpha)) return gp1_density_bound(x, y, tb);
  const S& alpha = tb.alpha;
  const S& lambda = tb.lambda;
  const int minxy = std::min(x, y);
  S nain = 0.0;
  for (int r = 0; r <= minxy; ++r) {
    nain = nain + factorial(y) / factorial(y - r) / factorial(r) * alpha * (1 - alpha) *
      tb.powA[r] * tb.powB[y - r] / tb.denom[y] *
      lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
  }
  return nain;
}

// ------------------------------------------------------------------ order 2
//
// P(x | y, z) = sum_r g(r | y, z) * GP(x - r; lambda, eta) with
// g(r | y, z) = zwerg(r) / pyz, where zwerg is a triple convolution over
// (s, v, w) and pyz the bivariate GP normalizer. pyz does not depend on r.

template <class S>
struct GP2ParamsT {
  S lambda, alpha1, alpha2, alpha3, eta;
  S U, beta1, beta2, beta3, zeta, lu13, lu13a;

  void init(const S& lambda_, const S& alpha1_, const S& alpha2_, const S& alpha3_,
            const S& eta_) {
    lambda = lambda_; alpha1 = alpha1_; alpha2 = alpha2_; alpha3 = alpha3_; eta = eta_;
    U = 1 / (1 - alpha1 - alpha2 - alpha3);
    beta1 = lambda * U * alpha1;
    beta2 = lambda * U * alpha2;
    beta3 = lambda * U * alpha3;
    zeta = lambda * U * (1 - 2 * alpha1 - alpha3);
    lu13 = lambda * U * (1 - alpha1 - alpha3);
    lu13a = lambda * U * (alpha1 + alpha3);
  }

  // Same fields as init(), with the parameter-only factors passed in:
  // U_ = 1 / (1 - alpha1 - alpha2 - alpha3), fz = 1 - 2 * alpha1 - alpha3,
  // f13 = 1 - alpha1 - alpha3, f13a = alpha1 + alpha3. Since lambda * U * a is
  // evaluated left to right as (lambda * U) * a, every field is bitwise equal
  // to the one init() computes.
  void init_hoisted(const S& lambda_, const S& alpha1_, const S& alpha2_,
                    const S& alpha3_, const S& eta_, const S& U_, const S& fz,
                    const S& f13, const S& f13a) {
    lambda = lambda_; alpha1 = alpha1_; alpha2 = alpha2_; alpha3 = alpha3_; eta = eta_;
    U = U_;
    const S lU = lambda * U;
    beta1 = lU * alpha1;
    beta2 = lU * alpha2;
    beta3 = lU * alpha3;
    zeta = lU * fz;
    lu13 = lU * f13;
    lu13a = lU * f13a;
  }
};
typedef GP2ParamsT<double> GP2Params;

template <class S>
struct GP2TablesT {
  GP2ParamsT<S> p;
  std::vector<S> powB1, expB1, powB2, expB2, powB3, expB3, powL, expL, powZ, expZ;
  int cap;

  void init(const GP2ParamsT<S>& p_) {
    p = p_;
    cap = -1;
  }

  void ensure(int n) {
    if (n <= cap) return;
    check_table_index(n);
    powB1.resize(n + 1); expB1.resize(n + 1);
    powB2.resize(n + 1); expB2.resize(n + 1);
    powB3.resize(n + 1); expB3.resize(n + 1);
    powL.resize(n + 1);  expL.resize(n + 1);
    powZ.resize(n + 1);  expZ.resize(n + 1);
    for (int i = cap + 1; i <= n; ++i) {
      powB1[i] = pow(p.beta1 + i * p.eta, i - 1);
      expB1[i] = exp(-p.beta1 - i * p.eta);
      powB2[i] = pow(p.beta2 + i * p.eta, i - 1);
      expB2[i] = exp(-p.beta2 - i * p.eta);
      powB3[i] = pow(p.beta3 + i * p.eta, i - 1);
      expB3[i] = exp(-p.beta3 - i * p.eta);
      powL[i]  = pow(p.lambda + i * p.eta, i - 1);
      expL[i]  = exp(-p.lambda - i * p.eta);
      powZ[i]  = pow(p.zeta + i * p.eta, i - 1);
      expZ[i]  = exp(-p.zeta - i * p.eta);
    }
    cap = n;
  }
};
typedef GP2TablesT<double> GP2Tables;

// gp2_pyz when lu13a = lambda U (alpha1 + alpha3) is exactly zero, i.e.
// alpha1 = alpha3 = 0 (the Julia backend returns boundary estimates exactly
// at 0). The j = 0 summand of gp2_pyz carries the factor
// (beta1 + beta3) * pow(lu13a, -1) of the GP pmf of rate lu13a at 0, which is
// 0 * Inf = NaN there. Since beta1 + beta3 = lambda U alpha1 + lambda U alpha3
// equals lu13a identically in the parameters, that factor is 1 (its limit, as
// gp2_rate_factor), so the j = 0 summand is evaluated without it. The
// summands j >= 1 are those of gp2_pyz (finite at lu13a = 0, the exponent of
// pow(lu13a + eta j, j - 1) is non-negative). Both parts are the analytic
// expression of P(y, z) at every parameter value, so dual numbers carry the
// exact partials (coco_ad::pow has zero partials for the exponent 0).
template <class S>
inline S gp2_pyz_zero_rate(int y, int z, const GP2ParamsT<S>& p) {
  const S scale = pow(p.lu13, 2) *
    exp(-(p.beta1 + p.beta3) - 2 * p.lu13 - y * p.eta - z * p.eta);
  const S head = pow(p.lu13 + p.eta * y, y - 1) / factorial(y) *
    pow(p.lu13 + p.eta * z, z - 1) / factorial(z);
  const int minyz = std::min(y, z);
  S tail = 0.0;
  for (int j = 1; j <= minyz; ++j) {
    tail = tail + pow(p.lu13 + p.eta * (y - j), y - j - 1) / factorial(y - j) *
      pow(p.lu13 + p.eta * (z - j), z - j - 1) / factorial(z - j) *
      pow(p.lu13a + p.eta * j, j - 1) / factorial(j) * exp(j * p.eta);
  }
  return (head + (p.beta1 + p.beta3) * tail) * scale;
}

// Bivariate GP normalizer P(y, z); r-invariant. With lu13a exactly zero,
// gp2_pyz_zero_rate (the arithmetic for lu13a != 0 is unchanged).
template <class S>
inline S gp2_pyz(int y, int z, const GP2ParamsT<S>& p) {
  if (is_zero_value(p.lu13a)) return gp2_pyz_zero_rate(y, z, p);
  const S help = (p.beta1 + p.beta3) * pow(p.lu13, 2) *
    exp(-(p.beta1 + p.beta3) - 2 * p.lu13 - y * p.eta - z * p.eta);
  const int minyz = std::min(y, z);
  S dwarf = 0.0;
  for (int j = 0; j <= minyz; ++j) {
    dwarf = dwarf + pow(p.lu13 + p.eta * (y - j), y - j - 1) / factorial(y - j) *
      pow(p.lu13 + p.eta * (z - j), z - j - 1) / factorial(z - j) *
      pow(p.lu13a + p.eta * j, j - 1) / factorial(j) * exp(j * p.eta);
  }
  return dwarf * help;
}

// Univariate GP pmf factor rate * pow(rate + i*eta, i-1) / i! * exp(-rate - i*eta)
// read from the tables of gp2_zwerg. At i = 0 the expression is
// rate * pow(rate, -1) * exp(-rate), which is 0 * Inf = NaN at a zero rate;
// there the factor is its limit exp(-rate) = 1 (as in Coconots.jl, whose GP
// pmf at 0 is exp(-rate) for every rate).
template <class S>
inline S gp2_rate_factor(const S& rate, const std::vector<S>& powT,
                         const std::vector<S>& expT, int i) {
  if (i == 0 && is_zero_value(rate)) return expT[0];
  return rate * powT[i] / factorial(i) * expT[i];
}

// gp2_zwerg when one of the rates beta1, beta2, beta3, zeta is exactly zero,
// e.g. at a fitted alpha on its boundary 0 (the Julia backend returns
// boundary estimates exactly at 0): same feasible index set and factors, with
// the zero-rate factor at index 0 replaced by its limit (gp2_rate_factor).
template <class S>
inline S gp2_zwerg_zero_rate(int r, int y, int z, const GP2ParamsT<S>& p,
                             const GP2TablesT<S>& tb) {
  S zwerg = 0.0;
  const int smax = std::min(y, r);
  for (int s = 0; s <= smax; ++s) {
    const int vlo = std::max(0, r - z);
    const int vhi = std::min(r - s, y - s);
    for (int v = vlo; v <= vhi; ++v) {
      const int k = r - s - v;
      const int wmax = std::min(z - r + v, y - s - v);
      const S pre = gp2_rate_factor(p.beta3, tb.powB3, tb.expB3, s) *
        gp2_rate_factor(p.beta1, tb.powB1, tb.expB1, v) *
        gp2_rate_factor(p.beta2, tb.powB2, tb.expB2, k);
      for (int w = 0; w <= wmax; ++w) {
        zwerg = zwerg + pre * gp2_rate_factor(p.beta1, tb.powB1, tb.expB1, w) *
          gp2_rate_factor(p.lambda, tb.powL, tb.expL, z - r + v - w) *
          gp2_rate_factor(p.zeta, tb.powZ, tb.expZ, y - s - v - w);
      }
    }
  }
  return zwerg;
}

// Triple convolution over (s, v, w). The loop bounds enumerate exactly the
// index triples that satisfy r-s-v >= 0, z-r+v-w >= 0 and y-s-v-w >= 0, in
// the same lexicographic order the former full 0..y cube visited them.
// Tables must cover indices up to max(y, z) (all six factor indices are
// bounded by y or z on the feasible set).
// With a zero rate beta1, beta2, beta3 or zeta, gp2_zwerg_zero_rate.
template <class S>
inline S gp2_zwerg(int r, int y, int z, const GP2ParamsT<S>& p, const GP2TablesT<S>& tb) {
  if (is_zero_value(p.beta1) || is_zero_value(p.beta2) || is_zero_value(p.beta3) ||
      is_zero_value(p.zeta)) {
    return gp2_zwerg_zero_rate(r, y, z, p, tb);
  }
  S zwerg = 0.0;
  const int smax = std::min(y, r);
  for (int s = 0; s <= smax; ++s) {
    const int vlo = std::max(0, r - z);
    const int vhi = std::min(r - s, y - s);
    for (int v = vlo; v <= vhi; ++v) {
      const int k = r - s - v;
      const int wmax = std::min(z - r + v, y - s - v);
      // leading subchain of the product; hoisting a prefix of a
      // left-associative chain does not change the arithmetic
      const S pre = p.beta3 * tb.powB3[s] / factorial(s) * tb.expB3[s] *
        p.beta1 * tb.powB1[v] / factorial(v) * tb.expB1[v] * p.beta1;
      for (int w = 0; w <= wmax; ++w) {
        zwerg = zwerg + pre * tb.powB1[w] / factorial(w) * tb.expB1[w] *
          p.beta2 * tb.powB2[k] / factorial(k) * tb.expB2[k] *
          p.lambda * tb.powL[z - r + v - w] / factorial(z - r + v - w) * tb.expL[z - r + v - w] *
          p.zeta * tb.powZ[y - s - v - w] / factorial(y - s - v - w) * tb.expZ[y - s - v - w];
      }
    }
  }
  return zwerg;
}

// g(r | y, z) for r = 0..y+z plus the value used beyond the support
// (0 * pow(pyz, -1), which is NaN when pyz underflows to 0 -- kept for
// exact behavioural equivalence of the simulation while-loops).
template <class S>
struct GVecT {
  std::vector<S> g;
  S tail;
};
typedef GVecT<double> GVec;

// Tables must already cover max(y, z, 0); with const tables this is safe to
// run concurrently for distinct `out` targets.
template <class S>
void gp2_build_g(int y, int z, const GP2ParamsT<S>& p, const GP2TablesT<S>& tb, GVecT<S>& out) {
  const S pyz = gp2_pyz(y, z, p);
  const int rmax = y + z;
  out.g.resize(rmax + 1);
  for (int r = 0; r <= rmax; ++r) {
    out.g[r] = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
  }
  out.tail = 0.0 * pow(pyz, -1);
}

// Factored evaluation of gp2_zwerg (same math, summation reordered).
//
// Each term of the triple convolution is the product of six univariate GP
// pmf values  q3[s] q1[v] q2[k] q1[w] qL[m - w] qZ[n - w]  with k = r-s-v,
// m = z-s-k and n = y-s-v, where qX[i] = betaX * powX[i] / i! * expX[i]
// (the same subchains as in gp2_zwerg). The w-sum depends on (m, n) only:
//   W(m, n) = sum_{w=0}^{min(m,n)} q1[w] qL[m-w] qZ[n-w],
// so it is tabulated once and zwerg(r) becomes a double sum. The feasible set
// of gp2_zwerg is exactly s in [0, min(y,z)], v in [0, y-s], k in [0, z-s]
// with r = s+v+k, which is enumerated directly. Only the association and
// summation order of the (non-negative, for admissible parameters) terms
// differ from gp2_zwerg, i.e. results agree up to floating-point rounding.
template <class S>
struct GP2ConvT {
  std::vector<S> q1, q2, q3, qL, qZ, W, zw, A, H;
  int mlo, mhi, nlo, nhi;
  std::size_t stride;

  // W(m, n) for m in [mlo_, mhi_], n in [nlo_, nhi_], with the products
  // q1[w] * qZ[n - w] formed once per n (array A) instead of once per (m, n).
  // Pmf tables must cover q1 and qZ up to nhi_ and qL up to mhi_.
  void build_W_factored(int mlo_, int mhi_, int nlo_, int nhi_) {
    mlo = mlo_; mhi = mhi_; nlo = nlo_; nhi = nhi_;
    stride = static_cast<std::size_t>(nhi - nlo + 1);
    W.resize(static_cast<std::size_t>(mhi - mlo + 1) * stride);
    A.resize(static_cast<std::size_t>(std::min(mhi, nhi) + 1));
    for (int n = nlo; n <= nhi; ++n) {
      const int wtop = std::min(mhi, n);
      for (int w = 0; w <= wtop; ++w) A[w] = q1[w] * qZ[n - w];
      for (int m = mlo; m <= mhi; ++m) {
        const int wmax = std::min(m, n);
        S acc = 0.0;
        for (int w = 0; w <= wmax; ++w) acc = acc + A[w] * qL[m - w];
        W[static_cast<std::size_t>(m - mlo) * stride + static_cast<std::size_t>(n - nlo)] = acc;
      }
    }
  }

  inline const S& Wat(int m, int n) const {
    return W[static_cast<std::size_t>(m - mlo) * stride + static_cast<std::size_t>(n - nlo)];
  }

  // out[r] = zwerg(r | y, z) for r = 0..rcap, with y, z >= 0, 0 <= rcap <= y+z.
  // W must cover m in [max(0, z-rcap), z] and n in [max(0, y-rcap), y].
  void zwerg_upto(int y, int z, int rcap, std::vector<S>& out) const {
    out.assign(rcap + 1, S(0.0));
    const int smax = std::min(std::min(y, z), rcap);
    for (int s = 0; s <= smax; ++s) {
      const S a3 = q3[s];
      const int vmax = std::min(y - s, rcap - s);
      for (int v = 0; v <= vmax; ++v) {
        const S a31 = a3 * q1[v];
        const int n = y - s - v;
        const int kmax = std::min(z - s, rcap - s - v);
        S* o = &out[s + v];
        for (int k = 0; k <= kmax; ++k) {
          o[k] = o[k] + a31 * q2[k] * Wat(z - s - k, n);
        }
      }
    }
  }

  // Staircase variant of W: for every n in [0, slo.size()) only the m in
  // [slo[n], shi[n]] are tabulated (none when slo[n] > shi[n]); the caller sets
  // slo and shi. Each entry is computed exactly as in build_W_factored (same
  // products A[w] = q1[w] * qZ[n - w], same w-order), so it is bitwise equal to
  // the build_W_factored entry. q1 and qZ must cover every n with a non-empty
  // range, qL every shi[n].
  std::vector<int> slo, shi;
  std::vector<std::size_t> soff;

  void build_W_stair() {
    const int ncol = static_cast<int>(slo.size());
    soff.assign(static_cast<std::size_t>(ncol) + 1, 0);
    std::size_t tot = 0;
    int atop = -1;
    for (int n = 0; n < ncol; ++n) {
      soff[n] = tot;
      if (slo[n] <= shi[n]) {
        tot += static_cast<std::size_t>(shi[n] - slo[n] + 1);
        atop = std::max(atop, std::min(shi[n], n));
      }
    }
    soff[ncol] = tot;
    W.resize(tot);
    A.resize(static_cast<std::size_t>(atop + 1));
    for (int n = 0; n < ncol; ++n) {
      const int lo = slo[n], hi = shi[n];
      if (lo > hi) continue;
      const int wtop = std::min(hi, n);
      for (int w = 0; w <= wtop; ++w) A[w] = q1[w] * qZ[n - w];
      S* col = &W[soff[n]];
      for (int m = lo; m <= hi; ++m) {
        const int wmax = std::min(m, n);
        S acc = 0.0;
        for (int w = 0; w <= wmax; ++w) acc = acc + A[w] * qL[m - w];
        col[m - lo] = acc;
      }
    }
  }

  inline const S& Wst(int m, int n) const {
    return W[soff[n] + static_cast<std::size_t>(m - slo[n])];
  }

  // zwerg_upto on the staircase W (identical loops and summation order).
  // The staircase must cover m in [max(0, z-rcap), z] for every n in
  // [max(0, y-rcap), y] (the (m, n) read for r <= rcap).
  void zwerg_upto_stair(int y, int z, int rcap, std::vector<S>& out) const {
    out.assign(rcap + 1, S(0.0));
    const int smax = std::min(std::min(y, z), rcap);
    for (int s = 0; s <= smax; ++s) {
      const S a3 = q3[s];
      const int vmax = std::min(y - s, rcap - s);
      for (int v = 0; v <= vmax; ++v) {
        const S a31 = a3 * q1[v];
        const int n = y - s - v;
        const int kmax = std::min(z - s, rcap - s - v);
        S* o = &out[s + v];
        for (int k = 0; k <= kmax; ++k) {
          o[k] = o[k] + a31 * q2[k] * Wst(z - s - k, n);
        }
      }
    }
  }
};

// Upper bound on the number of W entries tabulated at once; beyond it the
// direct gp2_zwerg evaluation is used.
const double GP2_CONV_MAX_ENTRIES = 4e6;

// Table-size limits are stated in double entries; for the wider dual-number
// scalars they are scaled down to the same number of bytes.
template <class S>
inline double entries_limit(double n_double) {
  return n_double * sizeof(double) / sizeof(S);
}

// ---------------------------------------------- per-t tables (covariates)
//
// With covariates every rate depends on lambda_t, so the pmf tables are
// rebuilt for every t. The helpers below build them with fewer
// transcendental calls (same math; values agree up to rounding).

// Largest index for which the recurrences below are used; beyond it the
// former per-entry formula is kept (factorial() raises the count-limit error
// there, since 171! overflows).
const int GP_REC_MAX = 170;

// q[i] = rate * pow(rate + i*eta, i-1) / i! * exp(-rate - i*eta), i = 0..n
// (qX[i] in the GP2ConvT notation). For 1 <= i <= GP_REC_MAX the exponential
// is exp(-rate) * eneg[i] with eneg[i] = exp(-i*eta) tabulated once per
// likelihood call (one exp per table instead of one per entry); for the
// Poisson model (eta exactly zero) q[i] = q[i-1] * rate / i (no powers).
// q[0] keeps the former formula (rate / rate == rate * pow(rate, -1) up to
// rounding) except at a zero rate, where that is NaN and q[0] is its limit
// exp(-rate) (as gp2_rate_factor). eneg must cover min(n, GP_REC_MAX).
template <class S>
inline void gp_pmf_fill(std::vector<S>& q, int n, const S& rate, const S& eta,
                        bool poisson, const std::vector<S>& eneg) {
  q.resize(static_cast<std::size_t>(n) + 1);
  const S e0 = exp(-rate);
  q[0] = is_zero_value(rate) ? e0 : rate / rate * e0;
  const int nrec = std::min(n, GP_REC_MAX);
  if (poisson) {
    for (int i = 1; i <= nrec; ++i) q[i] = q[i - 1] * rate / static_cast<double>(i);
  } else {
    for (int i = 1; i <= nrec; ++i) {
      q[i] = rate * pow(rate + i * eta, i - 1) / factorial(i) * (e0 * eneg[i]);
    }
  }
  for (int i = nrec + 1; i <= n; ++i) {
    q[i] = rate * pow(rate + i * eta, i - 1) / factorial(i) * exp(-rate - i * eta);
  }
}

// gp2_pyz for 0 <= y, z <= GP_REC_MAX with exp(j*eta) read from
// epos[j] = exp(j*eta) (bitwise the value gp2_pyz computes). The factor
// pow(lu13 + eta*i, i-1) is shared by the y- and z-parts when their index
// ranges overlap (tabulated once in hbuf; same arguments, same values). For
// the Poisson model the summands follow the ratio
//   term(j+1) / term(j) = lu13a (y-j)(z-j) / (lu13^2 (j+1))
// starting from the j = 0 summand of gp2_pyz; this path is used only when
// term(0) is a normal number well above underflow (the summands are
// log-concave in j, so later underflowing terms are below rounding relative
// to the sum) and the sum is finite, otherwise the generic loop runs.
// With lu13a exactly zero, gp2_pyz_zero_rate (as gp2_pyz).
template <class S>
inline S gp2_pyz_cov(int y, int z, const GP2ParamsT<S>& p, bool poisson,
                     const std::vector<S>& epos, std::vector<S>& hbuf) {
  if (is_zero_value(p.lu13a)) return gp2_pyz_zero_rate(y, z, p);
  const S help = (p.beta1 + p.beta3) * pow(p.lu13, 2) *
    exp(-(p.beta1 + p.beta3) - 2 * p.lu13 - y * p.eta - z * p.eta);
  const int minyz = std::min(y, z);
  if (poisson && coco_ad::value_of(p.lu13) > 0 && coco_ad::value_of(p.lu13a) > 0) {
    S term = pow(p.lu13 + p.eta * y, y - 1) / factorial(y) *
      pow(p.lu13 + p.eta * z, z - 1) / factorial(z) *
      pow(p.lu13a, -1) / factorial(0) * epos[0];
    const double t0 = coco_ad::value_of(term);
    if (std::isfinite(t0) && t0 >= 1e-290) {
      const S ratio = p.lu13a / (p.lu13 * p.lu13);
      S dwarf = term;
      for (int j = 0; j < minyz; ++j) {
        term = term * ratio *
          (static_cast<double>(y - j) * static_cast<double>(z - j) / static_cast<double>(j + 1));
        dwarf = dwarf + term;
      }
      if (std::isfinite(coco_ad::value_of(dwarf))) return dwarf * help;
    }
  }
  S dwarf = 0.0;
  const int dyz = (y > z) ? (y - z) : (z - y);
  if (dyz <= minyz) {
    // {y-j} and {z-j} (j = 0..minyz) together cover 0..max(y, z)
    const int top = std::max(y, z);
    hbuf.resize(static_cast<std::size_t>(top) + 1);
    for (int i = 0; i <= top; ++i) hbuf[i] = pow(p.lu13 + p.eta * i, i - 1);
    for (int j = 0; j <= minyz; ++j) {
      dwarf = dwarf + hbuf[y - j] / factorial(y - j) *
        hbuf[z - j] / factorial(z - j) *
        pow(p.lu13a + p.eta * j, j - 1) / factorial(j) * epos[j];
    }
  } else {
    for (int j = 0; j <= minyz; ++j) {
      dwarf = dwarf + pow(p.lu13 + p.eta * (y - j), y - j - 1) / factorial(y - j) *
        pow(p.lu13 + p.eta * (z - j), z - j - 1) / factorial(z - j) *
        pow(p.lu13a + p.eta * j, j - 1) / factorial(j) * epos[j];
    }
  }
  return dwarf * help;
}

template <class S>
using GCacheT = std::unordered_map<long long, GVecT<S> >;
typedef GCacheT<double> GCache;

inline long long g_key(int y, int z) {
  return (static_cast<long long>(y) << 32) | static_cast<unsigned int>(z);
}

template <class S>
inline GVecT<S>& g_lookup(GCacheT<S>& cache, int y, int z, const GP2ParamsT<S>& p,
                          GP2TablesT<S>& tb) {
  GVecT<S>& gv = cache[g_key(y, z)];
  if (gv.g.empty()) {
    tb.ensure(std::max(std::max(y, z), 0));
    gp2_build_g(y, z, p, tb, gv);
  }
  return gv;
}

// ------------------------------------------------------- log-likelihoods
//
// Negative log-likelihoods; the exported double entry points below and the
// automatic-differentiation entry point likelihoodAD() share these bodies.

// One observation of the order-2 likelihood with a valid history (y, z >= 0).
struct GP2Obs { int y, z, x; };

inline bool gp2obs_less(const GP2Obs& a, const GP2Obs& b) {
  if (a.y != b.y) return a.y < b.y;
  if (a.z != b.z) return a.z < b.z;
  return a.x < b.x;
}

// Order-2 negative log-likelihood without covariates. lambda is constant, so
// P(x | y, z) depends on t only through the triple (x, y, z):
// - observations are sorted by (y, z, x); identical triples contribute
//   count * log P(x | y, z) (same terms, summed in a different order);
// - for each distinct history (y, z), g(r | y, z) is built only for
//   r <= rcap = min(max x observed after (y, z), y + z), the r the
//   accumulation reads (g[r] for r <= rcap is bitwise the value of the full
//   vector: zwerg_upto adds the same contributions to out[r] in the same order);
// - W(m, n) is tabulated only over the staircase of (m, n) these r reach;
// - pmf tables come from gp_pmf_fill (one exp per table via exp(-i*eta),
//   Poisson ratio recurrence when eta is exactly zero), sized to the indices
//   read; the normalizer from gp2_pyz_cov. P(x | y, z) = sum_r g[r] qL[x-r].
// All of this is the same math as before; values agree up to rounding.
template <class S>
S likelihoodGP2_impl(const S& lambda, const S& alpha1, const S& alpha2,
                     const S& alpha3, const S& eta, int T, int seas1, int seas2,
                     const std::vector<int>& data, int cores)
{
  GP2ParamsT<S> p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);

  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  if (maxval >= 0) check_table_index(maxval);

  std::vector<GP2Obs> obs;
  std::vector<int> degen;  // t with a negative (e.g. NA) y or z
  if (T > seas2) obs.reserve(static_cast<std::size_t>(T - seas2));
  for (int t = seas2 + 1; t <= T; ++t) {
    const int x = data[t - 1];
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];
    if (y >= 0 && z >= 0) {
      GP2Obs o;
      o.y = y; o.z = z; o.x = x;
      obs.push_back(o);
    } else {
      degen.push_back(t);
    }
  }
  std::sort(obs.begin(), obs.end(), gp2obs_less);

  // distinct histories: obs[begin, end) share (y, z); rcap < 0: no g needed
  struct Hist { int y, z, rcap; std::size_t begin, end; };
  std::vector<Hist> hist;
  int maxy = 0, maxz = 0;
  for (std::size_t i = 0; i < obs.size(); ) {
    std::size_t j = i + 1;
    while (j < obs.size() && obs[j].y == obs[i].y && obs[j].z == obs[i].z) ++j;
    Hist h;
    h.y = obs[i].y; h.z = obs[i].z;
    h.rcap = std::min(obs[j - 1].x, h.y + h.z);
    h.begin = i; h.end = j;
    hist.push_back(h);
    maxy = std::max(maxy, h.y);
    maxz = std::max(maxz, h.z);
    i = j;
  }

  const bool use_conv = !hist.empty() &&
    static_cast<double>(maxz + 1) * (maxy + 1) <= entries_limit<S>(GP2_CONV_MAX_ENTRIES);

  std::vector<std::vector<S> > gv(hist.size());
  GP2ConvT<S> cv;
  GP2TablesT<S> tb;
  tb.init(p);

  if (use_conv) {
    // table extents over the histories that need g and the x read
    int ny = -1, nl = -1, n2 = -1, n3 = -1;
    for (std::size_t h = 0; h < hist.size(); ++h) {
      const Hist& hh = hist[h];
      if (hh.rcap < 0) continue;
      ny = std::max(ny, hh.y);
      nl = std::max(nl, hh.z);
      n2 = std::max(n2, std::min(hh.z, hh.rcap));
      n3 = std::max(n3, std::min(std::min(hh.y, hh.z), hh.rcap));
      nl = std::max(nl, obs[hh.end - 1].x);
    }
    if (ny >= 0) {
      const bool poisson = is_exact_zero(eta);
      const int nrec = std::max(0, std::min(maxval, GP_REC_MAX));
      std::vector<S> eneg, epos;
      if (poisson) {
        // exp(+-(i * 0)) with zero partials: exactly 1
        epos.assign(static_cast<std::size_t>(nrec) + 1, S(1.0));
      } else {
        eneg.resize(static_cast<std::size_t>(nrec) + 1);
        epos.resize(static_cast<std::size_t>(nrec) + 1);
        for (int i = 0; i <= nrec; ++i) {
          eneg[i] = exp(-(i * eta));
          epos[i] = exp(i * eta);
        }
      }
      gp_pmf_fill(cv.q1, ny, p.beta1, eta, poisson, eneg);
      gp_pmf_fill(cv.qZ, ny, p.zeta, eta, poisson, eneg);
      gp_pmf_fill(cv.q2, n2, p.beta2, eta, poisson, eneg);
      gp_pmf_fill(cv.q3, n3, p.beta3, eta, poisson, eneg);
      gp_pmf_fill(cv.qL, nl, lambda, eta, poisson, eneg);

      // staircase: for each n the m-range read by any history
      cv.slo.assign(static_cast<std::size_t>(ny) + 1, std::numeric_limits<int>::max());
      cv.shi.assign(static_cast<std::size_t>(ny) + 1, -1);
      for (std::size_t h = 0; h < hist.size(); ++h) {
        const Hist& hh = hist[h];
        if (hh.rcap < 0) continue;
        const int mlo = std::max(0, hh.z - hh.rcap);
        for (int n = std::max(0, hh.y - hh.rcap); n <= hh.y; ++n) {
          cv.slo[n] = std::min(cv.slo[n], mlo);
          cv.shi[n] = std::max(cv.shi[n], hh.z);
        }
      }
      cv.build_W_stair();

      const int nthreads = resolve_threads(static_cast<int>(hist.size()), cores);
      std::vector<std::vector<S> > hbufs(nthreads);
      parallel_for(static_cast<int>(hist.size()), nthreads, [&](int tid, int i) {
        const Hist& hh = hist[i];
        if (hh.rcap < 0) return;
        const S pyz = (std::max(hh.y, hh.z) <= GP_REC_MAX) ?
          gp2_pyz_cov(hh.y, hh.z, p, poisson, epos, hbufs[tid]) : gp2_pyz(hh.y, hh.z, p);
        std::vector<S>& g = gv[i];
        cv.zwerg_upto_stair(hh.y, hh.z, hh.rcap, g);
        const S ipyz = pow(pyz, -1);
        for (std::size_t r = 0; r < g.size(); ++r) g[r] = g[r] * ipyz;
      });
    }
  } else if (!hist.empty()) {
    if (maxval >= 0) tb.ensure(maxval);
    const int nthreads = resolve_threads(static_cast<int>(hist.size()), cores);
    parallel_for(static_cast<int>(hist.size()), nthreads, [&](int, int i) {
      const Hist& hh = hist[i];
      if (hh.rcap < 0) return;
      GVecT<S> tmp;
      gp2_build_g(hh.y, hh.z, p, tb, tmp);
      gv[i].swap(tmp.g);
    });
  }

  S mlef = 0.0;
  for (std::size_t h = 0; h < hist.size(); ++h) {
    const Hist& hh = hist[h];
    const std::vector<S>& g = gv[h];
    for (std::size_t i = hh.begin; i < hh.end; ) {
      const int x = obs[i].x;
      std::size_t j = i + 1;
      while (j < hh.end && obs[j].x == x) ++j;
      const int minxyz = std::min(x, hh.y + hh.z);
      S nain = 0.0;
      if (use_conv) {
        for (int r = 0; r <= minxyz; ++r) nain = nain + g[r] * cv.qL[x - r];
      } else {
        for (int r = 0; r <= minxyz; ++r) {
          nain = nain + g[r] * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
        }
      }
      const S ln = log(nain);
      const std::size_t cnt = j - i;
      mlef = (cnt == 1) ? mlef + ln : mlef + static_cast<double>(cnt) * ln;
      i = j;
    }
  }

  if (!degen.empty() && maxval >= 0) tb.ensure(maxval);
  for (std::size_t i = 0; i < degen.size(); ++i) {
    const int t = degen[i];
    const int x = data[t - 1];
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];
    // degenerate inputs (e.g. NA counts): zwerg/pyz vanish just like in the
    // guarded cube loops, propagating the same 0 * Inf arithmetic
    const S pyz = gp2_pyz(y, z, p);
    const int minxyz = std::min(x, y + z);
    S nain = 0.0;
    for (int r = 0; r <= minxyz; ++r) {
      const S gryz = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
      nain = nain + gryz * lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
    }
    mlef = mlef + log(nain);
  }

  return (-mlef);
}

// Order-1 negative log-likelihood. P(x | y) depends on t only through (x, y):
// observations are sorted by (y, x) and identical pairs contribute
// count * log P(x | y). For each distinct y the r-dependent factor
//   K[r] = y! / (y-r)! / r! * alpha (1-alpha) pow(alpha + psi r, r-1)
//          * pow(1 - alpha + psi (y-r), y-r-1)
// is formed once, and P(x | y) = (sum_r K[r] L[x-r]) / pow(1 + psi y, y-1)
// with the GP innovation pmf L (as in gp1_density, same factors, different
// association; values agree up to rounding). Poisson (eta exactly zero): L by
// the ratio recurrence of gp_pmf_fill.
template <class S>
S likelihoodGP1_impl(const S& lambda, const S& alpha, const S& eta, int T, int seas,
                     const std::vector<int>& data)
{
  const S psi = eta * (1 - alpha) / lambda;

  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  if (maxval >= 0) check_table_index(maxval);

  std::vector<std::pair<int, int> > obs;  // (y, x)
  if (T > seas) obs.reserve(static_cast<std::size_t>(T - seas));
  for (int t = seas + 1; t <= T; ++t) {
    obs.push_back(std::make_pair(data[t - seas - 1], data[t - 1]));
  }
  std::sort(obs.begin(), obs.end());

  // table extents over the valid pairs (x, y >= 0)
  int maxy = -1, maxx = -1, maxr = -1;
  for (std::size_t i = 0; i < obs.size(); ++i) {
    const int y = obs[i].first, x = obs[i].second;
    if (x < 0 || y < 0) continue;
    maxy = std::max(maxy, y);
    maxx = std::max(maxx, x);
    maxr = std::max(maxr, std::min(x, y));
  }

  const S a1a = alpha * (1 - alpha);
  std::vector<S> powA(static_cast<std::size_t>(maxr + 1));
  std::vector<S> powB(static_cast<std::size_t>(maxy + 1));
  std::vector<S> L;
  if (gp1_alpha_on_bound(alpha)) {
    // K[r] = y! / (y-r)! / r! * powA[r] * powB[y-r] with the factors of
    // gp1_thin_bound (limit 1 instead of 0 * Inf at index 0)
    for (int i = 0; i <= maxr; ++i) powA[i] = gp1_thin_a(alpha, psi, i);
    for (int i = 0; i <= maxy; ++i) powB[i] = gp1_thin_b(alpha, psi, i);
  } else {
    for (int i = 0; i <= maxr; ++i) powA[i] = a1a * pow(alpha + psi * i, i - 1);
    for (int i = 0; i <= maxy; ++i) powB[i] = pow(1 - alpha + psi * i, i - 1);
  }
  if (maxx >= 0) {
    if (is_exact_zero(eta)) {
      const std::vector<S> none;
      gp_pmf_fill(L, maxx, lambda, eta, true, none);
    } else {
      L.resize(static_cast<std::size_t>(maxx) + 1);
      for (int i = 0; i <= maxx; ++i) {
        L[i] = lambda * pow(lambda + i * eta, i - 1) / factorial(i) * exp(-lambda - i * eta);
      }
    }
  }

  std::vector<S> K;
  S mlef = 0.0;
  for (std::size_t i = 0; i < obs.size(); ) {
    const int y = obs[i].first;
    std::size_t e = i + 1;
    while (e < obs.size() && obs[e].first == y) ++e;
    // obs[i, e) share y, sorted by x; K over r <= min(y, max x)
    const int rk = std::min(y, obs[e - 1].second);
    S denom = 1.0;
    if (y >= 0 && rk >= 0) {
      K.resize(static_cast<std::size_t>(rk) + 1);
      for (int r = 0; r <= rk; ++r) {
        K[r] = factorial(y) / factorial(y - r) / factorial(r) * powA[r] * powB[y - r];
      }
      denom = pow(1 + psi * y, y - 1);
    }
    while (i < e) {
      const int x = obs[i].second;
      std::size_t j = i + 1;
      while (j < e && obs[j].second == x) ++j;
      const int minxy = std::min(x, y);
      S nain = 0.0;
      if (minxy >= 0) {
        for (int r = 0; r <= minxy; ++r) nain = nain + K[r] * L[x - r];
        nain = nain / denom;
      }
      const S ln = log(nain);
      const std::size_t cnt = j - i;
      mlef = (cnt == 1) ? mlef + ln : mlef + static_cast<double>(cnt) * ln;
      i = j;
    }
  }

  return (-mlef);
}

template <class S, class V>
S likelihoodGP2cov_impl(const S& alpha1, const S& alpha2, const S& alpha3,
                        const S& eta, const V& lambdas, int T, int seas1, int seas2,
                        const std::vector<int>& data, const NumericMatrix& xreg,
                        LinkFunction link, int cores)
{
  // lambda_t varies with t, so nothing can be cached across observations --
  // but each observation's contribution is independent. Compute them
  // concurrently with thread-private tables and accumulate serially in t
  // order afterwards, keeping results bit-identical to the serial loop.
  // Bound the table sizes once here so workers never raise an R error.
  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  check_table_index(maxval);

  const int nterms = T - seas2;
  if (nterms <= 0) return S(0.0);
  std::vector<S> logn(nterms);

  // parameter-only quantities, once per call (bitwise the values that
  // GP2ParamsT::init and gp2_pyz compute per t)
  const S U = 1 / (1 - alpha1 - alpha2 - alpha3);
  const S fz = 1 - 2 * alpha1 - alpha3;
  const S f13 = 1 - alpha1 - alpha3;
  const S f13a = alpha1 + alpha3;
  const bool poisson = is_exact_zero(eta);
  const int nrec = std::max(0, std::min(maxval, GP_REC_MAX));
  std::vector<S> eneg(static_cast<std::size_t>(nrec) + 1);
  std::vector<S> epos(static_cast<std::size_t>(nrec) + 1);
  for (int i = 0; i <= nrec; ++i) {
    eneg[i] = exp(-(i * eta));
    epos[i] = exp(i * eta);
  }

  const int nthreads = resolve_threads(nterms, cores);
  std::vector<GP2ParamsT<S> > pv(nthreads);
  std::vector<GP2TablesT<S> > tbv(nthreads);
  std::vector<GP2ConvT<S> > cvv(nthreads);

  parallel_for(nterms, nthreads, [&](int tid, int i) {
    const int t = seas2 + 1 + i;
    const int x = data[t - 1];
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];

    const S mathelp = rowInnerProductT<S>(xreg, t - 1, lambdas);
    const S lambda = applyLink(mathelp, link);
    GP2ParamsT<S>& p = pv[tid];
    p.init_hoisted(lambda, alpha1, alpha2, alpha3, eta, U, fz, f13, f13a);

    const int minxyz = std::min(x, y + z);
    S nain = 0.0;
    if (x >= 0 && y >= 0 && z >= 0 &&
        static_cast<double>(y + 1) * (z + 1) <= entries_limit<S>(GP2_CONV_MAX_ENTRIES)) {
      // Every pmf table is built only over the indices it is read at for
      // this (x, y, z), r <= minxyz: q1, qZ on 0..y; q2 on 0..min(z, minxyz);
      // q3 on 0..min(y, z, minxyz); qL on 0..max(x, z) (W and the final
      // r-sum). W is restricted to the (m, n) that these r reach.
      GP2ConvT<S>& cv = cvv[tid];
      gp_pmf_fill(cv.q1, y, p.beta1, eta, poisson, eneg);
      gp_pmf_fill(cv.qZ, y, p.zeta, eta, poisson, eneg);
      gp_pmf_fill(cv.q2, std::min(z, minxyz), p.beta2, eta, poisson, eneg);
      gp_pmf_fill(cv.q3, std::min(std::min(y, z), minxyz), p.beta3, eta, poisson, eneg);
      gp_pmf_fill(cv.qL, std::max(x, z), lambda, eta, poisson, eneg);
      const S pyz = (std::max(y, z) <= GP_REC_MAX) ?
        gp2_pyz_cov(y, z, p, poisson, epos, cv.H) : gp2_pyz(y, z, p);
      cv.build_W_factored(std::max(0, z - minxyz), z, std::max(0, y - minxyz), y);
      cv.zwerg_upto(y, z, minxyz, cv.zw);
      const S ipyz = pow(pyz, -1);
      for (int r = 0; r <= minxyz; ++r) {
        const S gryz = cv.zw[r] * ipyz;
        nain = nain + gryz * cv.qL[x - r];
      }
    } else {
      GP2TablesT<S>& tb = tbv[tid];
      tb.init(p);
      tb.ensure(std::max(std::max(std::max(x, y), z), 0));
      const S pyz = gp2_pyz(y, z, p);
      for (int r = 0; r <= minxyz; ++r) {
        const S gryz = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
        nain = nain + gryz * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
      }
    }

    // guard branch: constant, i.e. zero derivative
    if (nain <= 0) {
      nain = 1e-12;
    }
    logn[i] = log(nain);
  });

  S mlef = 0.0;
  for (int i = 0; i < nterms; ++i) mlef = mlef + logn[i];

  return (-mlef);
}

template <class S, class V>
S likelihoodGP1cov_impl(const S& alpha, const S& eta, const V& lambdas, int T, int seas,
                        const std::vector<int>& data, const NumericMatrix& xreg,
                        LinkFunction link)
{
  const bool bound = gp1_alpha_on_bound(alpha);
  S mlef = 0.0;
  for (int t = seas + 1; t <= T; ++t) {
    const int x = data[t - 1];
    const int y = data[t - seas - 1];

    const S mathelp = rowInnerProductT<S>(xreg, t - 1, lambdas);
    const S lambda = applyLink(mathelp, link);
    const S psi = eta * (1 - alpha) / lambda;
    const S denom = (y >= 0) ? pow(1 + psi * y, y - 1) : S(0.0);

    const int minxy = std::min(x, y);
    S nain = 0.0;
    for (int r = 0; r <= minxy; ++r) {
      if (bound) {
        nain = nain + gp1_thin_bound(r, y, alpha, psi) / denom *
          lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
        continue;
      }
      nain = nain + factorial(y) / factorial(y - r) / factorial(r) * alpha * (1 - alpha) *
        pow(alpha + psi * r, r - 1) * pow(1 - alpha + psi * (y - r), y - r - 1) / denom *
        lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
    }
    mlef = mlef + log(nain);
  }

  return (-mlef);
}

// ------------------------------------------------ automatic differentiation
//
// Model of the R objective `mlf`: the free parameters `th` in the layout of
// coco$par (lambda first without covariates; alpha(s)[, eta], betas with
// covariates); eta = 0 for Poisson innovations.
struct ADModel {
  int order;
  bool gp;
  bool cov;
  int seas1, seas2;
  const std::vector<int>& data;
  const NumericMatrix& xreg;
  LinkFunction link;

  int n_dynamic() const { return (order == 1 ? 1 : 3) + (gp ? 1 : 0); }

  template <class S>
  S operator()(const std::vector<S>& th) const {
    const int T = static_cast<int>(data.size());
    const S zero = 0.0;
    if (!cov) {
      if (order == 1) {
        return likelihoodGP1_impl<S>(th[0], th[1], gp ? th[2] : zero, T, seas1, data);
      }
      return likelihoodGP2_impl<S>(th[0], th[1], th[2], th[3], gp ? th[4] : zero,
                                   T, seas1, seas2, data, 1);
    }
    const int k = n_dynamic();
    const std::vector<S> beta(th.begin() + k, th.end());
    if (order == 1) {
      return likelihoodGP1cov_impl<S>(th[0], gp ? th[1] : zero, beta, T, seas1,
                                      data, xreg, link);
    }
    return likelihoodGP2cov_impl<S>(th[0], th[1], th[2], gp ? th[3] : zero, beta,
                                    T, seas1, seas2, data, xreg, link, 1);
  }
};

// Value and gradient, C directions per evaluation (ForwardDiff "chunks").
template <int C, class F>
void ad_value_grad(const F& f, const std::vector<double>& th, double& val,
                   std::vector<double>& g) {
  typedef coco_ad::Dual<double, C> D;
  const int n = static_cast<int>(th.size());
  g.assign(n, 0.0);
  val = NA_REAL;
  for (int b0 = 0; b0 < n; b0 += C) {
    std::vector<D> x(n);
    for (int k = 0; k < n; ++k) {
      x[k] = D(th[k]);
      const int a = k - b0;
      if (a >= 0 && a < C) x[k].d[a] = 1.0;
    }
    const D y = f(x);
    val = y.v;
    for (int a = 0; a < C && b0 + a < n; ++a) g[b0 + a] = y.d[a];
  }
}

// Value, gradient and Hessian with nested duals: the inner partials seed the
// directions of block I, the outer partials those of block J; the outer
// partial b of the inner partial a is d^2 f / (d th[I+a] d th[J+b]).
template <int C, class F>
void ad_value_grad_hess(const F& f, const std::vector<double>& th, double& val,
                        std::vector<double>& g, std::vector<double>& H) {
  typedef coco_ad::Dual<double, C> D1;
  typedef coco_ad::Dual<D1, C> D2;
  const int n = static_cast<int>(th.size());
  g.assign(n, 0.0);
  H.assign(static_cast<std::size_t>(n) * n, 0.0);
  val = NA_REAL;
  for (int bi = 0; bi < n; bi += C) {
    for (int bj = bi; bj < n; bj += C) {
      std::vector<D2> x(n);
      for (int k = 0; k < n; ++k) {
        x[k] = D2(th[k]);
        const int a = k - bi;
        if (a >= 0 && a < C) x[k].v.d[a] = 1.0;
        const int b = k - bj;
        if (b >= 0 && b < C) x[k].d[b].v = 1.0;
      }
      const D2 y = f(x);
      val = y.v.v;
      for (int a = 0; a < C && bi + a < n; ++a) g[bi + a] = y.v.d[a];
      for (int b = 0; b < C && bj + b < n; ++b) g[bj + b] = y.d[b].v;
      for (int a = 0; a < C && bi + a < n; ++a) {
        for (int b = 0; b < C && bj + b < n; ++b) {
          const double h = y.d[b].d[a];
          H[static_cast<std::size_t>(bi + a) * n + (bj + b)] = h;
          if (bi != bj) H[static_cast<std::size_t>(bj + b) * n + (bi + a)] = h;
        }
      }
    }
  }
  // diagonal blocks: average the two (mathematically equal) mixed partials
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double m = 0.5 * (H[static_cast<std::size_t>(i) * n + j] +
                              H[static_cast<std::size_t>(j) * n + i]);
      H[static_cast<std::size_t>(i) * n + j] = m;
      H[static_cast<std::size_t>(j) * n + i] = m;
    }
  }
}

// Chunk size: the smallest instantiated size >= n (at most 8; larger
// parameter vectors are processed in blocks of 8).
template <class F>
void ad_eval(const F& f, const std::vector<double>& th, int deriv, double& val,
             std::vector<double>& g, std::vector<double>& H) {
  const int n = static_cast<int>(th.size());
  if (deriv >= 2) {
    if (n <= 2)      ad_value_grad_hess<2>(f, th, val, g, H);
    else if (n <= 3) ad_value_grad_hess<3>(f, th, val, g, H);
    else if (n <= 4) ad_value_grad_hess<4>(f, th, val, g, H);
    else if (n <= 5) ad_value_grad_hess<5>(f, th, val, g, H);
    else if (n <= 6) ad_value_grad_hess<6>(f, th, val, g, H);
    else             ad_value_grad_hess<8>(f, th, val, g, H);
  } else {
    if (n <= 2)      ad_value_grad<2>(f, th, val, g);
    else if (n <= 3) ad_value_grad<3>(f, th, val, g);
    else if (n <= 4) ad_value_grad<4>(f, th, val, g);
    else if (n <= 5) ad_value_grad<5>(f, th, val, g);
    else if (n <= 6) ad_value_grad<6>(f, th, val, g);
    else             ad_value_grad<8>(f, th, val, g);
  }
}

} // namespace


// [[Rcpp::export]]
double fact(int i)
{
  return factorial(i);
}


// [[Rcpp::export]]
double applyLinkFunction(double mathelp, std::string link_function) {
  return applyLink(mathelp, parseLinkFunction(link_function));
}


// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


// [[Rcpp::export]]
double likelihoodGP2(double lambda, double alpha1, double alpha2,
                     double alpha3, double eta, int T, int seas1, int seas2,
                     const std::vector<int>& data, int cores = 1)
{
  return likelihoodGP2_impl<double>(lambda, alpha1, alpha2, alpha3, eta, T, seas1, seas2,
                                    data, cores);
}


// [[Rcpp::export]]
double likelihoodGP1(double lambda, double alpha, double eta, int T, int seas,
                     const std::vector<int>& data)
{
  return likelihoodGP1_impl<double>(lambda, alpha, eta, T, seas, data);
}


// [[Rcpp::export]]
std::vector<int> simGP2(double lambda, double alpha1, double alpha2,
                        double alpha3, double eta, int T, int N, int seas1, int seas2,
                        std::vector<int> data, const std::vector<double>& uniform,
                        const std::vector<int>& innovations)
{
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);
  GCache cache;

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];

    const double unif = uniform[t - 1];

    // the support of the thinning part is 0..y + z; if the accumulated
    // probability stays <= unif up to y + z (rounding), y + z is drawn
    const int rmax = y + z;
    int r = 0;
    double nain = 0.0;
    if (y >= 0 && z >= 0) {
      const GVec& gv = g_lookup(cache, y, z, p, tb);
      while (nain <= unif && r <= rmax) {
        nain = nain + gv.g[r];
        r = r + 1;
      }
    } else {
      const double pyz = gp2_pyz(y, z, p);
      while (nain <= unif && r <= rmax) {
        nain = nain + gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
        r = r + 1;
      }
    }

    data[t - 1] = r - 1 + innovations[t - 1];
  }

  return (data);
}


// [[Rcpp::export]]
std::vector<int> simGP1(double lambda, double alpha, double eta, int T, int N,
                        int seas, std::vector<int> data,
                        const std::vector<double>& uniform, const std::vector<int>& innovations)
{
  const double psi = eta * (1 - alpha) / lambda;
  GP1Tables tb;
  tb.init(lambda, alpha, eta, psi);
  const bool bound = gp1_alpha_on_bound(alpha);

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas - 1];
    if (y >= 0) tb.ensure(y);

    const double unif = uniform[t - 1];

    // support 0..y; if the accumulated probability stays <= unif up to y
    // (rounding), y is drawn
    int r = 0;
    double nain = 0.0;
    while (nain <= unif && r <= y) {
      if (bound) {
        nain = nain + gp1_thin_bound(r, y, tb.alpha, tb.psi) / tb.denom[y];
      } else {
        nain = nain + factorial(y) / factorial(y - r) / factorial(r) * tb.alpha * (1 - tb.alpha) *
          tb.powA[r] * tb.powB[y - r] / tb.denom[y];
      }
      r = r + 1;
    }
    data[t - 1] = r - 1 + innovations[t - 1];
  }

  return (data);
}


// [[Rcpp::export]]
double dGP2h(int x, int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);
  tb.ensure(std::max(std::max(std::max(x, y), z), 0));

  const double pyz = gp2_pyz(y, z, p);
  const int minxyz = std::min(x, y + z);
  double nain = 0.0;
  for (int r = 0; r <= minxyz; ++r) {
    const double gryz = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
    nain = nain + gryz * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
  }

  return (nain);
}


// Transition densities P(x | y, z) for x = 0..xmax in one call; the
// (y, z)-dependent convolution work is shared across all x.
// [[Rcpp::export]]
NumericVector dGP2h_vec(int xmax, int y, int z, double lambda, double alpha1,
                        double alpha2, double alpha3, double eta)
{
  if (xmax < 0) return NumericVector(0);
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);
  tb.ensure(std::max(std::max(std::max(xmax, y), z), 0));

  NumericVector out(xmax + 1);
  if (y >= 0 && z >= 0) {
    GVec gv;
    gp2_build_g(y, z, p, tb, gv);
    for (int x = 0; x <= xmax; ++x) {
      const int minxyz = std::min(x, y + z);
      double nain = 0.0;
      for (int r = 0; r <= minxyz; ++r) {
        nain = nain + gv.g[r] * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
      }
      out[x] = nain;
    }
  } else {
    for (int x = 0; x <= xmax; ++x) out[x] = dGP2h(x, y, z, lambda, alpha1, alpha2, alpha3, eta);
  }
  return out;
}


// [[Rcpp::export]]
double dGP1h(int x, int y, double lambda, double alpha, double eta)
{
  const double psi = eta * (1 - alpha) / lambda;
  // only referenced for y >= 0 (the r loop is empty otherwise)
  const double denom = (y >= 0) ? pow(1 + psi * y, y - 1) : 0.0;
  const int minxy = std::min(x, y);

  const bool bound = gp1_alpha_on_bound(alpha);

  double nain = 0.0;
  for (int r = 0; r <= minxy; ++r) {
    if (bound) {
      nain = nain + gp1_thin_bound(r, y, alpha, psi) / denom *
        lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
      continue;
    }
    nain = nain + factorial(y) / factorial(y - r) / factorial(r) * alpha * (1 - alpha) *
      pow(alpha + psi * r, r - 1) * pow(1 - alpha + psi * (y - r), y - r - 1) / denom *
      lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
  }

  return (nain);
}


// Transition densities P(x | y) for x = 0..xmax in one call.
// [[Rcpp::export]]
NumericVector dGP1h_vec(int xmax, int y, double lambda, double alpha, double eta)
{
  if (xmax < 0) return NumericVector(0);
  const double psi = eta * (1 - alpha) / lambda;
  GP1Tables tb;
  tb.init(lambda, alpha, eta, psi);
  tb.ensure(std::max(std::max(xmax, y), 0));

  NumericVector out(xmax + 1);
  for (int x = 0; x <= xmax; ++x) out[x] = gp1_density(x, y, tb);
  return out;
}


// [[Rcpp::export]]
double dR2(int r, int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);
  tb.ensure(std::max(std::max(y, z), 0));

  const double pyz = gp2_pyz(y, z, p);
  return (gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1));
}


// P(R = r | y, z) for r = 0, 1, ... until the cumulative probability
// reaches 1 - val_num, mirroring the accumulation loop in cocoResid(). The
// support of R given (y, z) is 0..y+z, so the loop stops there at the latest:
// rounding can keep the cumulative mass below 1 - val_num for tiny val_num,
// and the former unbounded loop then never terminated (the terms beyond
// y + z are zero).
// [[Rcpp::export]]
NumericVector dR2_vec(int y, int z, double lambda, double alpha1, double alpha2,
                      double alpha3, double eta, double val_num)
{
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);
  tb.ensure(std::max(std::max(y, z), 0));

  const double pyz = gp2_pyz(y, z, p);
  const int rmax = std::max(y, 0) + std::max(z, 0);
  std::vector<double> probs;
  double acc = 0.0;
  int r = 0;
  while (acc < 1 - val_num && r <= rmax) {
    const double prob = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
    probs.push_back(prob);
    acc = acc + prob;
    r = r + 1;
  }
  return wrap(probs);
}


// [[Rcpp::export]]
double Pyz(int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  return gp2_pyz(y, z, p);
}


// [[Rcpp::export]]
std::vector<int> simGP2cov(double alpha1, double alpha2, double alpha3,
                           double eta, NumericVector lambdas, int T, int N, int seas1, int seas2,
                           std::vector<int> data, NumericMatrix xreg,
                           const std::vector<double>& uniform, const std::vector<int>& innovations,
                           std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);
  GP2Params p;
  GP2Tables tb;

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];

    const double mathelp = rowInnerProduct(xreg, t - N - 1, lambdas);
    const double lambda = applyLink(mathelp, link);
    check_sim_lambda(lambda);
    p.init(lambda, alpha1, alpha2, alpha3, eta);
    tb.init(p);
    tb.ensure(std::max(std::max(y, z), 0));

    const double pyz = gp2_pyz(y, z, p);
    const double unif = uniform[t - N - 1];

    // support 0..y + z; if the accumulated probability stays <= unif up to
    // y + z (rounding), y + z is drawn
    const int rmax = y + z;
    int r = 0;
    double nain = 0.0;
    while (nain <= unif && r <= rmax) {
      nain = nain + gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
      r = r + 1;
    }

    data[t - 1] = r - 1 + innovations[t - N - 1];
  }

  return (data);
}


// [[Rcpp::export]]
std::vector<int> simGP1cov(double alpha, double eta,
                           NumericVector lambdas, int T, int N, int seas, std::vector<int> data,
                           NumericMatrix xreg, const std::vector<double>& uniform,
                           const std::vector<int>& innovations, std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);
  const bool bound = gp1_alpha_on_bound(alpha);

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas - 1];

    const double mathelp = rowInnerProduct(xreg, t - 1 - N, lambdas);
    const double lambda = applyLink(mathelp, link);
    check_sim_lambda(lambda);
    const double psi = eta * (1 - alpha) / lambda;
    const double denom = (y >= 0) ? pow(1 + psi * y, y - 1) : 0.0;

    const double unif = uniform[t - 1];
    // support 0..y; if the accumulated probability stays <= unif up to y
    // (rounding), y is drawn
    int r = 0;
    double nain = 0.0;
    while (nain <= unif && r <= y) {
      if (bound) {
        nain = nain + gp1_thin_bound(r, y, alpha, psi) / denom;
      } else {
        nain = nain + factorial(y) / factorial(y - r) / factorial(r) * alpha * (1 - alpha) *
          pow(alpha + psi * r, r - 1) * pow(1 - alpha + psi * (y - r), y - r - 1) / denom;
      }
      r = r + 1;
    }
    data[t - 1] = r - 1 + innovations[t - N - 1];
  }

  return (data);
}

namespace {

// Generalized Poisson pmf P(X = x; lambda1, lambda2), reproducing the
// arithmetic of HMMpa::dgenpois (term order, the x < 2 closed form, the
// product loop for x >= 2 and its replacement of NaN/Inf by the smallest
// denormal), so the inversion below picks the same indices as the former
// R implementation of rgenpois().
inline double hmmpa_dgenpois(int x, double l1, double l2) {
  if (x < 2) {
    const double f = l1 + x * l2;
    // factorial(0) = factorial(1) = 1; std::pow matches R's ^ for the
    // exponents -1 and 0 used here (including the Inf/NaN/0 bases)
    return l1 * pow(f, x - 1.0) * exp(-f);
  }
  const double f1 = l1 + x * l2;
  const double g1 = exp(-l2);
  double e = 1.0;
  for (int i = 2; i <= x; ++i) e = e * ((f1 * g1) / i);
  const double d1 = l1 * exp(-l1) * g1;
  double b = e * d1;
  if (std::isnan(b) || std::isinf(b)) b = 4.94065645841247e-324;
  return b;
}

// Inversion of the Generalized Poisson CDF with the cumulative table kept for
// one lambda: the smallest x with F(x) >= u, F accumulated from x = 0 upwards
// with the strict `F(x) < u` continuation rule of the former R loop.
//
// Termination. In floating point the accumulated F can stay below u forever:
// its increments round to zero once they fall below half an ulp of F, and for
// large lambda the HMMpa recursion returns only denormals. The former loop
// then never ended. Two rules bound the inversion:
//  - when the next value F(x) is still < u, the increment vanished
//    (F(x) == F(x - 1)) and 1 - F(x - 1) <= GP_INV_SAT_TOL, F has saturated:
//    beyond the bulk the pmf decreases, so F keeps this value for all larger
//    x, the remaining mass is at rounding level, and x is returned (the exact
//    quantile is >= x and is needed only with probability <= 1 - F(x - 1));
//  - a draw above MAX_COUNT = 170 is outside the supported range of counts
//    (factorials overflow at 171, see factorial()): once F(170) < u the
//    inversion stops with the count-limit error. This also ends the loop
//    when F cannot be evaluated accurately (large lambda, denormal pmf).
// Draws <= 170 that the former loop returned are unchanged, and no
// additional random numbers are used.
const double GP_INV_SAT_TOL = 1e-8;

struct GPInvCdf {
  double lam, eta;
  std::vector<double> cdf;
  long long iter;
  int sat;    // x returned once F saturated (-1: not saturated)

  GPInvCdf() : lam(0.0), eta(0.0), iter(0), sat(-1) {}

  void reset(double lam_, double eta_) {
    lam = lam_;
    eta = eta_;
    cdf.clear();
    cdf.push_back(hmmpa_dgenpois(0, lam, eta));
    sat = -1;
  }

  int draw(double u) {
    int idx = 0;
    while (true) {
      if (idx >= static_cast<int>(cdf.size())) {
        if (sat >= 0) return sat;
        // here F(idx - 1) < u, so the draw is >= idx
        if (idx > MAX_COUNT)
          stop("coconots: counts above %d are not supported (factorial overflow in "
               "the C++ backend); a simulated Generalized Poisson innovation exceeds %d "
               "(lambda = %g, eta = %g)", MAX_COUNT, MAX_COUNT, lam, eta);
        const double prev = cdf.back();
        const double next = prev + hmmpa_dgenpois(idx, lam, eta);
        if (next < u && !(next > prev) && 1.0 - prev <= GP_INV_SAT_TOL) {
          // saturated (false if next or u is NaN: error below)
          sat = idx;
          return idx;
        }
        cdf.push_back(next);
      }
      const double last = cdf[idx];
      if (std::isnan(last) || std::isnan(u))
        stop("missing value where TRUE/FALSE needed");
      if (!(last < u)) return idx;
      idx = idx + 1;
      if ((++iter & 0xFFFFF) == 0) checkUserInterrupt();
    }
  }
};

} // namespace


// Inversion sampler of the Generalized Poisson distribution: for each
// uniform u[i], the smallest x with cumulative pmf (accumulated from x = 0
// upwards, as in the former R loop) >= u[i], at lambda = lambdas[i]
// (GPInvCdf, which also documents the termination rule).
// The cumulative table is reused while consecutive lambdas are equal.
// [[Rcpp::export]]
NumericVector rgenpoisInv(NumericVector u, NumericVector lambdas, double eta)
{
  const int n = u.size();
  if (lambdas.size() != n) stop("u and lambdas must have the same length");
  NumericVector out(n);
  GPInvCdf inv;
  bool have = false;

  for (int i = 0; i < n; ++i) {
    const double lam = lambdas[i];
    if (!have || !(lam == inv.lam)) {
      inv.reset(lam, eta);
      have = true;
    }
    out[i] = inv.draw(u[i]);
  }
  return out;
}


// ------------------------------------------------------------ cocoBoot
//
// Parametric bootstrap of the sample ACF (cocoBoot, RCPP path) in one C++
// call: all replicates are simulated here and their ACFs computed, instead
// of 1000 R-level cocoSim() + forecast::Acf() calls.
//
// Per replicate the random numbers are drawn in the same order and with the
// same R generators as cocoSim (initial values, innovations, thinning
// uniforms), and every thinning step uses the same arithmetic as
// simGP1/simGP2/simGP1cov/simGP2cov, so each replicate equals the series
// cocoSim would return for the same RNG state.
//
// Speed comes from memoizing the cumulative thinning distribution: without
// covariates it depends on the history (y, z) only, so it is shared by all
// steps and replicates; with covariates lambda_t is shared by all
// replicates, so the replicates are simulated time-major in blocks and the
// distribution is shared per (t, y, z). Cache entries are extended lazily
// (only as far as some uniform required), and the cache size is capped;
// beyond the cap the distribution is evaluated directly, as before.

namespace {

// Cumulative distribution c[j] = sum_{i <= j} term(i), accumulated from
// i = 0 upwards exactly like the `nain` loops of the sim functions.
struct LazyCdf {
  std::vector<double> c;
  double aux;   // order 2: pow(pyz, -1); order 1 with covariates: denom
  bool ready;
  LazyCdf() : aux(0.0), ready(false) {}
};

// Soft cap on the doubles stored in one cache (~80 MB).
const std::size_t BOOT_CACHE_MAX_DOUBLES = 10000000;

struct CdfStore {
  std::unordered_map<long long, LazyCdf> map;
  std::size_t stored;
  LazyCdf scratch;

  CdfStore() : stored(0) {}

  void clear() {
    map.clear();
    stored = 0;
  }

  LazyCdf* get(long long key) {
    std::unordered_map<long long, LazyCdf>::iterator it = map.find(key);
    if (it != map.end()) return &it->second;
    if (stored < BOOT_CACHE_MAX_DOUBLES) {
      stored += 8;  // node overhead
      return &map[key];
    }
    scratch.c.clear();
    scratch.ready = false;
    return &scratch;
  }

  void account(const LazyCdf* e, std::size_t before) {
    if (e != &scratch) stored += e->c.size() - before;
  }
};

// Direct-indexed pointer cache in front of CdfStore::get for histories with
// 0 <= y, z < DIM (order 1: z = 0), which avoids the hash lookup per
// thinning draw. Elements of std::unordered_map keep their address until
// they are erased, so a cached pointer stays the LazyCdf that
// store.get(g_key(y, z)) returns until store.clear(); reset() must be called
// together with every store.clear(). The scratch entry (returned once the
// store is full) is never cached, and store.get is called exactly when the
// key has not been looked up since the last reset, so the store's size
// accounting is unchanged.
struct CdfIndex {
  static const int DIM = 64;
  std::vector<LazyCdf*> ptr;
  std::vector<unsigned int> stamp;
  unsigned int gen;

  CdfIndex() : ptr(DIM * DIM, static_cast<LazyCdf*>(0)), stamp(DIM * DIM, 0u), gen(1u) {}

  void reset() {
    if (++gen == 0u) {
      std::fill(stamp.begin(), stamp.end(), 0u);
      gen = 1u;
    }
  }

  LazyCdf* get(CdfStore& store, int y, int z) {
    const bool in_range = y >= 0 && z >= 0 && y < DIM && z < DIM;
    const int idx = in_range ? y * DIM + z : 0;
    if (in_range && stamp[idx] == gen) return ptr[idx];
    LazyCdf* e = store.get(g_key(y, z));
    if (in_range && e != &store.scratch) {
      ptr[idx] = e;
      stamp[idx] = gen;
    }
    return e;
  }
};

// Thinning draw of order 2 for uniform u: the index r returned by the
// simGP2 / simGP2cov while-loop, i.e. the first r with !(c[r] <= u).
// The loops add zwerg(r) * pow(pyz, -1) also beyond r = y + z (where zwerg
// is 0); if c[y + z + 1] <= u the former loop never terminated, here the
// largest feasible value y + z is returned. Tables must cover max(y, z).
inline int gp2_cdf_draw(LazyCdf& e, int y, int z, double u,
                        const GP2Params& p, const GP2Tables& tb) {
  if (!e.ready) {
    e.aux = pow(gp2_pyz(y, z, p), -1);
    e.c.clear();
    e.ready = true;
  }
  const int rmax = y + z;
  for (int j = 0; ; ++j) {
    if (j == static_cast<int>(e.c.size())) {
      if (j > rmax + 1) return rmax;
      const double prev = (j == 0) ? 0.0 : e.c[j - 1];
      e.c.push_back(prev + gp2_zwerg(j, y, z, p, tb) * e.aux);
    }
    if (!(e.c[j] <= u)) return j;
  }
}

// Per-t tables of the order-2 thinning distribution with covariates (same
// math as gp2_zwerg / gp2_pyz, evaluated as in likelihoodGP2cov_impl):
// univariate pmf tables by gp_pmf_fill and the w-sum
//   W(m, n) = sum_{w=0}^{min(m,n)} q1[w] qL[m-w] qZ[n-w]
// tabulated once per t for all histories (y, z) of the block, so that
// zwerg(r | y, z) is a double sum over W (GP2ConvT::zwerg_upto). Tables
// are sized to the largest y and z present at t; W row n holds only
// m = 0..mz[n], with mz[n] the largest z of a history with y >= n (the
// entries the histories present can read). Values agree with gp2_zwerg /
// gp2_pyz up to floating-point rounding.
struct GP2BootTabs {
  GP2Params p;
  double eta;
  bool poisson;
  int ecap, ld;
  std::vector<double> q1, q2, q3, qL, qZ, W, A, zbuf, hbuf, eneg, epos;

  GP2BootTabs() : eta(0.0), poisson(true), ecap(-1), ld(0) {}

  void init_eta(double eta_) {
    eta = eta_;
    poisson = is_exact_zero(eta_);
    ecap = -1;
    eneg.clear();
    epos.clear();
  }

  // eneg[i] = exp(-(i * eta)), epos[i] = exp(i * eta) for i = 0..min(n, GP_REC_MAX)
  void ensure_exp(int n) {
    n = std::min(n, GP_REC_MAX);
    if (n <= ecap) return;
    eneg.resize(static_cast<std::size_t>(n) + 1);
    epos.resize(static_cast<std::size_t>(n) + 1);
    for (int i = ecap + 1; i <= n; ++i) {
      eneg[i] = exp(-(i * eta));
      epos[i] = exp(i * eta);
    }
    ecap = n;
  }

  // mz must have maxy + 1 entries (non-increasing, mz[0] = maxz).
  void build(const GP2Params& p_, int maxy, int maxz, const std::vector<int>& mz) {
    p = p_;
    ensure_exp(std::max(maxy, maxz));
    gp_pmf_fill(q1, maxy, p.beta1, eta, poisson, eneg);
    gp_pmf_fill(qZ, maxy, p.zeta, eta, poisson, eneg);
    gp_pmf_fill(q2, maxz, p.beta2, eta, poisson, eneg);
    gp_pmf_fill(q3, std::min(maxy, maxz), p.beta3, eta, poisson, eneg);
    gp_pmf_fill(qL, maxz, p.lambda, eta, poisson, eneg);
    ld = maxz + 1;
    W.resize(static_cast<std::size_t>(maxy + 1) * static_cast<std::size_t>(ld));
    A.resize(static_cast<std::size_t>(std::min(maxy, maxz)) + 1);
    for (int n = 0; n <= maxy; ++n) {
      const int mtop = mz[n];
      if (mtop < 0) continue;
      const int wtop = std::min(mtop, n);
      for (int w = 0; w <= wtop; ++w) A[w] = q1[w] * qZ[n - w];
      double* row = &W[static_cast<std::size_t>(n) * static_cast<std::size_t>(ld)];
      for (int m = 0; m <= mtop; ++m) {
        const int wmax = std::min(m, n);
        double acc = 0.0;
        for (int w = 0; w <= wmax; ++w) acc = acc + A[w] * qL[m - w];
        row[m] = acc;
      }
    }
  }

  double pyz(int y, int z) {
    return (std::max(y, z) <= GP_REC_MAX) ?
      gp2_pyz_cov(y, z, p, poisson, epos, hbuf) : gp2_pyz(y, z, p);
  }

  // zbuf[r - rlo] = zwerg(r | y, z) for r = rlo..rhi, 0 <= rlo <= rhi <= y + z.
  // Terms are those of gp2_zwerg: s in [0, min(y, z)], v in [0, y - s],
  // k in [0, z - s], r = s + v + k, read at W(z - s - k, y - s - v).
  void zwerg_range(int y, int z, int rlo, int rhi) {
    zbuf.assign(static_cast<std::size_t>(rhi - rlo) + 1, 0.0);
    const int smax = std::min(std::min(y, z), rhi);
    for (int s = 0; s <= smax; ++s) {
      const double a3 = q3[s];
      const int vmax = std::min(y - s, rhi - s);
      for (int v = 0; v <= vmax; ++v) {
        const int klo = std::max(0, rlo - s - v);
        const int khi = std::min(z - s, rhi - s - v);
        if (klo > khi) continue;
        const double a31 = a3 * q1[v];
        const double* Wn = &W[static_cast<std::size_t>(y - s - v) * static_cast<std::size_t>(ld)];
        const int base = s + v - rlo;  // may be negative; base + k >= 0
        for (int k = klo; k <= khi; ++k) {
          double& o = zbuf[static_cast<std::size_t>(base + k)];
          o = o + a31 * q2[k] * Wn[z - s - k];
        }
      }
    }
  }
};

// gp2_cdf_draw with zwerg from GP2BootTabs (built for this t, covering
// (y, z)). The cumulative table is extended in chunks of doubling length
// (r = j..min(y + z, 2j - 1), at least 8 values), each value accumulated
// as c[r] = c[r - 1] + zwerg(r) * pow(pyz, -1) like the former loop; the
// stopping rule is unchanged.
inline int gp2cov_cdf_draw(LazyCdf& e, int y, int z, double u, GP2BootTabs& bt) {
  if (!e.ready) {
    e.aux = pow(bt.pyz(y, z), -1);
    e.c.clear();
    e.ready = true;
  }
  const int rmax = y + z;
  for (int j = 0; ; ++j) {
    if (j == static_cast<int>(e.c.size())) {
      if (j > rmax + 1) return rmax;
      double prev = (j == 0) ? 0.0 : e.c[j - 1];
      if (j > rmax) {
        e.c.push_back(prev + 0.0 * e.aux);
      } else {
        const int rhi = std::min(rmax, j + std::max(j, 8) - 1);
        bt.zwerg_range(y, z, j, rhi);
        for (int r = j; r <= rhi; ++r) {
          prev = prev + bt.zbuf[static_cast<std::size_t>(r - j)] * e.aux;
          e.c.push_back(prev);
        }
      }
    }
    if (!(e.c[j] <= u)) return j;
  }
}

// Thinning draw of order 1 without covariates (simGP1 loop); terms vanish
// beyond r = y, where the former loop never terminated if c[y] <= u.
inline int gp1_cdf_draw(LazyCdf& e, int y, double u, const GP1Tables& tb) {
  e.ready = true;
  for (int j = 0; ; ++j) {
    if (j == static_cast<int>(e.c.size())) {
      if (j > y) return y;
      const double prev = (j == 0) ? 0.0 : e.c[j - 1];
      if (gp1_alpha_on_bound(tb.alpha)) {
        e.c.push_back(prev + gp1_thin_bound(j, y, tb.alpha, tb.psi) / tb.denom[y]);
      } else {
        e.c.push_back(prev + factorial(y) / factorial(y - j) / factorial(j) * tb.alpha *
                      (1 - tb.alpha) * tb.powA[j] * tb.powB[y - j] / tb.denom[y]);
      }
    }
    if (!(e.c[j] <= u)) return j;
  }
}

// Thinning draw of order 1 with covariates (simGP1cov loop) at psi(lambda_t).
inline int gp1cov_cdf_draw(LazyCdf& e, int y, double u, double alpha, double psi) {
  if (!e.ready) {
    e.aux = pow(1 + psi * y, y - 1);
    e.c.clear();
    e.ready = true;
  }
  for (int j = 0; ; ++j) {
    if (j == static_cast<int>(e.c.size())) {
      if (j > y) return y;
      const double prev = (j == 0) ? 0.0 : e.c[j - 1];
      if (gp1_alpha_on_bound(alpha)) {
        e.c.push_back(prev + gp1_thin_bound(j, y, alpha, psi) / e.aux);
      } else {
        e.c.push_back(prev + factorial(y) / factorial(y - j) / factorial(j) * alpha * (1 - alpha) *
                      pow(alpha + psi * j, j - 1) * pow(1 - alpha + psi * (y - j), y - j - 1) / e.aux);
      }
    }
    if (!(e.c[j] <= u)) return j;
  }
}

inline int count_from_double(double v) {
  if (!R_FINITE(v) || v < 0 || v > 2147483647.0)
    stop("Simulated innovation is not a finite non-negative count (invalid innovation mean?)");
  return static_cast<int>(v);
}

// Sample ACF at lags 1..nlags as stats::acf(type = "correlation",
// demean = TRUE) computes it for a complete series: mean accumulated in long
// double (colMeans), lag-k autocovariance sum / n, correlation
// acov / (se * se) with se = sqrt(acov[0]), clamped to [-1, 1]. Lags beyond
// n - 1 are NA (stats::acf truncates lag.max to n - 1).
void sample_acf(const int* x, int n, int nlags, std::vector<double>& xd, double* out) {
  for (int k = 0; k < nlags; ++k) out[k] = NA_REAL;
  if (n < 2) return;
  long double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[i];
  s /= n;
  const double mean = static_cast<double>(s);
  xd.resize(static_cast<std::size_t>(n) + nlags + 1);  // [n, n + nlags] holds the lag sums
  for (int i = 0; i < n; ++i) xd[i] = static_cast<double>(x[i]) - mean;
  double acov0 = 0.0;
  for (int i = 0; i < n; ++i) acov0 += xd[i] * xd[i];
  acov0 = acov0 / n;
  const double se = sqrt(acov0);
  const double den = se * se;
  const int lagmax = std::min(nlags, n - 1);
  // sum[lag] = sum_{i = 0}^{n - 1 - lag} xd[i + lag] * xd[i], accumulated for
  // all lags in one pass over i: each sum[lag] receives the same terms in the
  // same (increasing i) order as a separate loop per lag, but the lagmax
  // accumulators are independent, so the loop is not bound by the latency
  // of one serial addition chain and the inner loop vectorizes.
  double* __restrict__ sum = &xd[n];
  std::fill(sum, sum + lagmax + 1, 0.0);
  const double* __restrict__ xp = xd.data();
  const int nfull = n - lagmax;  // i < nfull: all lags 1..lagmax have a term
  for (int i = 0; i < nfull; ++i) {
    const double xi = xp[i];
    const double* __restrict__ xl = xp + i;
    for (int lag = 1; lag <= lagmax; ++lag) sum[lag] += xl[lag] * xi;
  }
  for (int i = std::max(nfull, 0); i < n - 1; ++i) {
    const double xi = xp[i];
    const int lm = n - 1 - i;
    for (int lag = 1; lag <= lm; ++lag) sum[lag] += xp[i + lag] * xi;
  }
  for (int lag = 1; lag <= lagmax; ++lag) {
    const double a = (sum[lag] / n) / den;
    out[lag - 1] = (a > 1.) ? 1. : ((a < -1.) ? -1. : a);
  }
}

} // namespace


// Bootstrap ACFs without covariates. Each replicate is the cocoSim_base
// series of length `size` (initial values, innovations and thinning
// uniforms drawn as there); its last n_keep values enter the ACF.
// Returns an nlags x nB matrix (column b = replicate b).
// [[Rcpp::export]]
NumericMatrix bootAcfBase(int order, bool gp, double lambda, double alpha1,
                          double alpha2, double alpha3, double eta,
                          int size, int n_keep, int nlags, int nB)
{
  if (order != 1 && order != 2) stop("order must be 1 or 2");
  const int N = order;
  if (size <= N || n_keep < 1 || n_keep > size || nlags < 1 || nB < 0)
    stop("invalid bootstrap dimensions");

  NumericMatrix out(nlags, nB);
  std::vector<int> data(size), innov(size);
  std::vector<double> unif(size), xd;

  GPInvCdf inv;
  if (gp) inv.reset(lambda, eta);

  GP1Tables tb1;
  GP2Params p2;
  GP2Tables tb2;
  if (order == 1) {
    tb1.init(lambda, alpha1, eta, eta * (1 - alpha1) / lambda);
  } else {
    p2.init(lambda, alpha1, alpha2, alpha3, eta);
    tb2.init(p2);
  }
  CdfStore store;
  CdfIndex cdf_index;

  for (int b = 0; b < nB; ++b) {
    // same draw order as cocoSim_base: initial values, innovations, uniforms
    for (int i = 0; i < N; ++i)
      data[i] = gp ? inv.draw(R::runif(0.0, 1.0)) : count_from_double(R::rpois(lambda));
    if (gp) {
      for (int i = 0; i < size; ++i) innov[i] = inv.draw(R::runif(0.0, 1.0));
    } else {
      for (int i = 0; i < size; ++i) innov[i] = count_from_double(R::rpois(lambda));
    }
    for (int i = 0; i < size; ++i) unif[i] = R::runif(0.0, 1.0);

    for (int t = N + 1; t <= size; ++t) {
      const double u = unif[t - 1];
      int r;
      if (order == 1) {
        const int y = data[t - 2];
        tb1.ensure(y);
        LazyCdf* e = cdf_index.get(store, y, 0);
        const std::size_t before = e->c.size();
        r = gp1_cdf_draw(*e, y, u, tb1);
        store.account(e, before);
      } else {
        const int y = data[t - 2];
        const int z = data[t - 3];
        tb2.ensure(std::max(y, z));
        LazyCdf* e = cdf_index.get(store, y, z);
        const std::size_t before = e->c.size();
        r = gp2_cdf_draw(*e, y, z, u, p2, tb2);
        store.account(e, before);
      }
      data[t - 1] = r + innov[t - 1];
    }

    sample_acf(&data[size - n_keep], n_keep, nlags, xd,
               out.begin() + static_cast<std::size_t>(b) * nlags);
    checkUserInterrupt();
  }
  return out;
}


// Bootstrap ACFs with covariates. Each replicate is the cocoSim_cov series
// of length `size` without burn-in: initial values from init_lambdas (the
// innovation distribution at xreg rows 1..order, GP draws with eta for
// type GP), innovations from lambdas (one per xreg row), n_unif thinning
// uniforms; observation t (1-based) uses lambda_t from xreg row t and
// innovation t, as in the likelihood, and the uniform t - order (order 2) or
// t (order 1), as cocoSim_cov passes them to simGP2cov / simGP1cov.
// Replicates are drawn in order (so the RNG stream is that of consecutive
// cocoSim calls) and then simulated time-major in blocks, sharing per-t
// tables and caches.
// [[Rcpp::export]]
NumericMatrix bootAcfCov(int order, bool gp, double alpha1, double alpha2,
                         double alpha3, double eta, NumericVector betas,
                         NumericMatrix xreg, NumericVector lambdas,
                         NumericVector init_lambdas, int size, int n_unif,
                         int nlags, int nB, std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);
  if (order != 1 && order != 2) stop("order must be 1 or 2");
  const int N = order;
  const int n_innov = lambdas.size();
  if (size <= N || nlags < 1 || nB < 0) stop("invalid bootstrap dimensions");
  if (init_lambdas.size() != N) stop("init_lambdas must have one value per initial observation");
  if (betas.size() != xreg.ncol()) stop("betas and xreg do not match");
  if (xreg.nrow() < size || n_innov < size) stop("xreg has too few rows");
  if (n_unif < ((order == 2) ? size - N : size)) stop("too few thinning uniforms");
  for (int i = 0; i < size; ++i) check_sim_lambda(lambdas[i]);
  for (int i = 0; i < N; ++i) check_sim_lambda(init_lambdas[i]);

  NumericMatrix out(nlags, nB);
  const std::size_t per_rep = static_cast<std::size_t>(size) + n_innov + n_unif + 1;
  const std::size_t blk_cap = std::max<std::size_t>(1, 4000000 / per_rep);
  const int blk = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(std::max(nB, 1)), blk_cap));

  std::vector<int> data(static_cast<std::size_t>(blk) * size);
  std::vector<double> innov(static_cast<std::size_t>(blk) * n_innov);  // Poisson: count, GP: uniform
  std::vector<double> unif(static_cast<std::size_t>(blk) * n_unif);
  std::vector<double> xd;

  GP2Params p2;
  GP2Tables tb2;
  GP2BootTabs bt2;
  bt2.init_eta(eta);
  std::vector<int> mz;
  CdfStore store;
  CdfIndex cdf_index;
  GPInvCdf inv, inv_init;

  for (int b0 = 0; b0 < nB; b0 += blk) {
    const int nb = std::min(blk, nB - b0);

    // random numbers per replicate, in the order of cocoSim_cov
    for (int k = 0; k < nb; ++k) {
      int* d = &data[static_cast<std::size_t>(k) * size];
      for (int i = 0; i < N; ++i) {
        if (gp) {
          inv_init.reset(init_lambdas[i], eta);
          d[i] = inv_init.draw(R::runif(0.0, 1.0));
        } else {
          d[i] = count_from_double(R::rpois(init_lambdas[i]));
        }
      }
      double* iv = &innov[static_cast<std::size_t>(k) * n_innov];
      if (gp) {
        for (int i = 0; i < n_innov; ++i) iv[i] = R::runif(0.0, 1.0);
      } else {
        for (int i = 0; i < n_innov; ++i) iv[i] = R::rpois(lambdas[i]);
      }
      double* uv = &unif[static_cast<std::size_t>(k) * n_unif];
      for (int i = 0; i < n_unif; ++i) uv[i] = R::runif(0.0, 1.0);
    }

    for (int t = N + 1; t <= size; ++t) {
      const int row = t - 1;
      const int ui = (order == 2) ? (t - N - 1) : (t - 1);
      const double lam = applyLink(rowInnerProduct(xreg, row, betas), link);
      check_sim_lambda(lam);
      const double psi = eta * (1 - alpha1) / lam;
      bool conv = false;
      if (order == 2) {
        p2.init(lam, alpha1, alpha2, alpha3, eta);
        // histories (y, z) = (d[t-2], d[t-3]) of the block at this t
        int maxy = 0, maxz = 0;
        for (int k = 0; k < nb; ++k) {
          const int* d = &data[static_cast<std::size_t>(k) * size];
          maxy = std::max(maxy, d[t - 2]);
          maxz = std::max(maxz, d[t - 3]);
        }
        conv = static_cast<double>(maxy + 1) * (maxz + 1) <= GP2_CONV_MAX_ENTRIES;
        if (conv) {
          check_table_index(std::max(maxy, maxz));
          mz.assign(static_cast<std::size_t>(maxy) + 1, -1);
          for (int k = 0; k < nb; ++k) {
            const int* d = &data[static_cast<std::size_t>(k) * size];
            if (d[t - 3] > mz[d[t - 2]]) mz[d[t - 2]] = d[t - 3];
          }
          for (int n = maxy - 1; n >= 0; --n) mz[n] = std::max(mz[n], mz[n + 1]);
          bt2.build(p2, maxy, maxz, mz);
        } else {
          tb2.init(p2);
        }
      }
      store.clear();
      cdf_index.reset();
      if (gp) inv.reset(lambdas[row], eta);

      for (int k = 0; k < nb; ++k) {
        int* d = &data[static_cast<std::size_t>(k) * size];
        const double u = unif[static_cast<std::size_t>(k) * n_unif + ui];
        int r;
        if (order == 1) {
          const int y = d[t - 2];
          LazyCdf* e = cdf_index.get(store, y, 0);
          const std::size_t before = e->c.size();
          r = gp1cov_cdf_draw(*e, y, u, alpha1, psi);
          store.account(e, before);
        } else {
          const int y = d[t - 2];
          const int z = d[t - 3];
          LazyCdf* e = cdf_index.get(store, y, z);
          const std::size_t before = e->c.size();
          if (conv) {
            r = gp2cov_cdf_draw(*e, y, z, u, bt2);
          } else {
            tb2.ensure(std::max(y, z));
            r = gp2_cdf_draw(*e, y, z, u, p2, tb2);
          }
          store.account(e, before);
        }
        const double iv = innov[static_cast<std::size_t>(k) * n_innov + row];
        const int x = gp ? inv.draw(iv) : count_from_double(iv);
        d[t - 1] = r + x;
      }
    }

    for (int k = 0; k < nb; ++k) {
      sample_acf(&data[static_cast<std::size_t>(k) * size], size, nlags, xd,
                 out.begin() + static_cast<std::size_t>(b0 + k) * nlags);
    }
    checkUserInterrupt();
  }
  return out;
}



// ------------------------------------------------------------ k-step forecast
//
// Simulated forecast paths of predict(k > 1) (cocoForecastKStepsRCPP) in one
// C++ call, instead of number_simulations R-level cocoSim_base /
// cocoSim_cov(size = k + order, init = <last observations>) calls.
//
// Per path the random numbers are drawn in the same order and with the same
// R generators as there: the `order` initial values (drawn and then replaced
// by init, as cocoSim does), the innovations, the thinning uniforms. The
// thinning draws use the memoized cumulative distributions of cocoBoot
// (gp1_cdf_draw, gp2_cdf_draw, gp1cov_cdf_draw, gp2cov_cdf_draw: same
// arithmetic as simGP1/simGP2/simGP1cov/simGP2cov up to rounding), so each
// path is the series the former call returned for the same RNG state.
// Returns an nsim x nsteps integer matrix (row = path, column = step).

// Without covariates (cocoSim_base): innovations and uniforms at t = 1..size,
// the thinning at step t uses uniform t - 1 and innovation t - 1.
// [[Rcpp::export]]
IntegerMatrix forecastPathsBase(int order, bool gp, double lambda, double alpha1,
                                double alpha2, double alpha3, double eta,
                                IntegerVector init, int nsteps, int nsim)
{
  if (order != 1 && order != 2) stop("order must be 1 or 2");
  const int N = order;
  if (init.size() != N) stop("init must hold one value per lag");
  if (nsteps < 1 || nsim < 0) stop("invalid forecast dimensions");
  const int size = N + nsteps;

  IntegerMatrix out(nsim, nsteps);
  std::vector<int> data(size), innov(size);
  std::vector<double> unif(size);

  GPInvCdf inv;
  if (gp) inv.reset(lambda, eta);

  GP1Tables tb1;
  GP2Params p2;
  GP2Tables tb2;
  if (order == 1) {
    tb1.init(lambda, alpha1, eta, eta * (1 - alpha1) / lambda);
  } else {
    p2.init(lambda, alpha1, alpha2, alpha3, eta);
    tb2.init(p2);
  }
  CdfStore store;

  for (int s = 0; s < nsim; ++s) {
    // initial values as drawn by cocoSim_base, then replaced by init
    for (int i = 0; i < N; ++i) {
      if (gp) {
        R::runif(0.0, 1.0);
      } else {
        R::rpois(lambda);
      }
      data[i] = init[i];
    }
    if (gp) {
      for (int i = 0; i < size; ++i) innov[i] = inv.draw(R::runif(0.0, 1.0));
    } else {
      for (int i = 0; i < size; ++i) innov[i] = count_from_double(R::rpois(lambda));
    }
    for (int i = 0; i < size; ++i) unif[i] = R::runif(0.0, 1.0);

    for (int t = N + 1; t <= size; ++t) {
      const double u = unif[t - 1];
      int r;
      if (order == 1) {
        const int y = data[t - 2];
        tb1.ensure(y);
        LazyCdf* e = store.get(g_key(y, 0));
        const std::size_t before = e->c.size();
        r = gp1_cdf_draw(*e, y, u, tb1);
        store.account(e, before);
      } else {
        const int y = data[t - 2];
        const int z = data[t - 3];
        tb2.ensure(std::max(y, z));
        LazyCdf* e = store.get(g_key(y, z));
        const std::size_t before = e->c.size();
        r = gp2_cdf_draw(*e, y, z, u, p2, tb2);
        store.account(e, before);
      }
      data[t - 1] = r + innov[t - 1];
    }

    for (int j = 0; j < nsteps; ++j) out(s, j) = data[N + j];
    if ((s & 255) == 255) checkUserInterrupt();
  }
  return out;
}


// With covariates (cocoSim_cov(init = <last observations>), xreg = xcast):
// initial values from init_lambdas (drawn, GP draws with eta for type GP,
// then replaced by init), innovations from lambdas (one per xreg row),
// n_unif thinning uniforms; forecast step j = t - order (data position t,
// 1-based, after the `order` observations of init) uses lambda_j from xcast
// row j, innovation j and the uniform j (order 2) or t (order 1), as
// cocoSim_cov(init = ...) passes them to simGP2cov / simGP1cov. Paths are
// drawn in order (so the RNG stream is that of consecutive cocoSim_cov calls)
// and then simulated time-major in blocks, sharing per-t tables and caches,
// as in bootAcfCov.
// [[Rcpp::export]]
IntegerMatrix forecastPathsCov(int order, bool gp, double alpha1, double alpha2,
                               double alpha3, double eta, NumericVector betas,
                               NumericMatrix xreg, NumericVector lambdas,
                               NumericVector init_lambdas, IntegerVector init,
                               int nsteps, int n_unif, int nsim,
                               std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);
  if (order != 1 && order != 2) stop("order must be 1 or 2");
  const int N = order;
  const int n_innov = lambdas.size();
  if (nsteps < 1 || nsim < 0) stop("invalid forecast dimensions");
  const int size = N + nsteps;
  if (init.size() != N) stop("init must hold one value per lag");
  if (init_lambdas.size() != N) stop("init_lambdas must have one value per initial observation");
  if (betas.size() != xreg.ncol()) stop("betas and xreg do not match");
  if (xreg.nrow() < nsteps || n_innov < nsteps) stop("xreg has too few rows");
  if (n_unif < ((order == 2) ? nsteps : size)) stop("too few thinning uniforms");

  IntegerMatrix out(nsim, nsteps);
  const std::size_t per_rep = static_cast<std::size_t>(size) + n_innov + n_unif + 1;
  const std::size_t blk_cap = std::max<std::size_t>(1, 4000000 / per_rep);
  const int blk = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(std::max(nsim, 1)), blk_cap));

  std::vector<int> data(static_cast<std::size_t>(blk) * size);
  std::vector<double> innov(static_cast<std::size_t>(blk) * n_innov);  // Poisson: count, GP: uniform
  std::vector<double> unif(static_cast<std::size_t>(blk) * n_unif);

  GP2Params p2;
  GP2Tables tb2;
  GP2BootTabs bt2;
  bt2.init_eta(eta);
  std::vector<int> mz;
  CdfStore store;
  GPInvCdf inv;

  for (int b0 = 0; b0 < nsim; b0 += blk) {
    const int nb = std::min(blk, nsim - b0);

    // random numbers per path, in the order of cocoSim_cov
    for (int k = 0; k < nb; ++k) {
      int* d = &data[static_cast<std::size_t>(k) * size];
      for (int i = 0; i < N; ++i) {
        // GP: rgenpois(n = 1, ...) (one uniform); Poisson: rpois
        if (gp) {
          R::runif(0.0, 1.0);
        } else {
          R::rpois(init_lambdas[i]);
        }
        d[i] = init[i];
      }
      double* iv = &innov[static_cast<std::size_t>(k) * n_innov];
      if (gp) {
        for (int i = 0; i < n_innov; ++i) iv[i] = R::runif(0.0, 1.0);
      } else {
        for (int i = 0; i < n_innov; ++i) iv[i] = R::rpois(lambdas[i]);
      }
      double* uv = &unif[static_cast<std::size_t>(k) * n_unif];
      for (int i = 0; i < n_unif; ++i) uv[i] = R::runif(0.0, 1.0);
    }

    for (int t = N + 1; t <= size; ++t) {
      const int row = t - N - 1;  // forecast step t - N uses xcast row t - N
      const int ui = (order == 2) ? (t - N - 1) : (t - 1);
      const double lam = applyLink(rowInnerProduct(xreg, row, betas), link);
      check_sim_lambda(lam);
      const double psi = eta * (1 - alpha1) / lam;
      bool conv = false;
      if (order == 2) {
        p2.init(lam, alpha1, alpha2, alpha3, eta);
        // histories (y, z) = (d[t-2], d[t-3]) of the block at this t
        int maxy = 0, maxz = 0;
        for (int k = 0; k < nb; ++k) {
          const int* d = &data[static_cast<std::size_t>(k) * size];
          maxy = std::max(maxy, d[t - 2]);
          maxz = std::max(maxz, d[t - 3]);
        }
        conv = static_cast<double>(maxy + 1) * (maxz + 1) <= GP2_CONV_MAX_ENTRIES;
        if (conv) {
          check_table_index(std::max(maxy, maxz));
          mz.assign(static_cast<std::size_t>(maxy) + 1, -1);
          for (int k = 0; k < nb; ++k) {
            const int* d = &data[static_cast<std::size_t>(k) * size];
            if (d[t - 3] > mz[d[t - 2]]) mz[d[t - 2]] = d[t - 3];
          }
          for (int n = maxy - 1; n >= 0; --n) mz[n] = std::max(mz[n], mz[n + 1]);
          bt2.build(p2, maxy, maxz, mz);
        } else {
          tb2.init(p2);
        }
      }
      store.clear();
      if (gp) inv.reset(lambdas[row], eta);

      for (int k = 0; k < nb; ++k) {
        int* d = &data[static_cast<std::size_t>(k) * size];
        const double u = unif[static_cast<std::size_t>(k) * n_unif + ui];
        int r;
        if (order == 1) {
          const int y = d[t - 2];
          LazyCdf* e = store.get(g_key(y, 0));
          const std::size_t before = e->c.size();
          r = gp1cov_cdf_draw(*e, y, u, alpha1, psi);
          store.account(e, before);
        } else {
          const int y = d[t - 2];
          const int z = d[t - 3];
          LazyCdf* e = store.get(g_key(y, z));
          const std::size_t before = e->c.size();
          if (conv) {
            r = gp2cov_cdf_draw(*e, y, z, u, bt2);
          } else {
            tb2.ensure(std::max(y, z));
            r = gp2_cdf_draw(*e, y, z, u, p2, tb2);
          }
          store.account(e, before);
        }
        const double iv = innov[static_cast<std::size_t>(k) * n_innov + row];
        const int x = gp ? inv.draw(iv) : count_from_double(iv);
        d[t - 1] = r + x;
      }
    }

    for (int k = 0; k < nb; ++k) {
      const int* d = &data[static_cast<std::size_t>(k) * size];
      for (int j = 0; j < nsteps; ++j) out(b0 + k, j) = d[N + j];
    }
    checkUserInterrupt();
  }
  return out;
}



// [[Rcpp::export]]
double likelihoodGP2cov(double alpha1, double alpha2, double alpha3,
                        double eta, NumericVector lambdas, int T, int seas1, int seas2,
                        const std::vector<int>& data, NumericMatrix xreg, std::string link_function,
                        int cores = 1)
{
  const LinkFunction link = parseLinkFunction(link_function);
  return likelihoodGP2cov_impl<double>(alpha1, alpha2, alpha3, eta, lambdas, T, seas1, seas2,
                                       data, xreg, link, cores);
}


// [[Rcpp::export]]
double likelihoodGP1cov(double alpha, double eta, NumericVector lambdas, int T, int seas,
                        const std::vector<int>& data, NumericMatrix xreg, std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);
  return likelihoodGP1cov_impl<double>(alpha, eta, lambdas, T, seas, data, xreg, link);
}


// Negative log-likelihood (the objective `mlf` of cocoReg) with its gradient
// (deriv = 1) or gradient and Hessian (deriv = 2) by forward-mode automatic
// differentiation. `par` is laid out as coco$par: lambda, alpha[, eta] /
// lambda, alpha1, alpha2, alpha3[, eta] without covariates and
// alpha[, eta], beta / alpha1, alpha2, alpha3[, eta], beta with covariates
// (xreg not NULL). Evaluated serially.
// [[Rcpp::export]]
List likelihoodAD(NumericVector par, int order, bool gp, const std::vector<int>& data,
                  int seas1, int seas2, SEXP xreg = R_NilValue,
                  std::string link_function = "log", int deriv = 2)
{
  if (order != 1 && order != 2) stop("order must be 1 or 2");
  const bool cov = !Rf_isNull(xreg);
  const NumericMatrix X = cov ? NumericMatrix(xreg) : NumericMatrix(0, 0);
  const LinkFunction link = cov ? parseLinkFunction(link_function) : LINK_LOG;

  const ADModel model{order, gp, cov, seas1, seas2, data, X, link};
  const int n_expected = cov ? model.n_dynamic() + X.ncol() : model.n_dynamic() + 1;
  if (par.size() != n_expected) stop("likelihoodAD: wrong number of parameters");
  if (cov && X.nrow() < static_cast<int>(data.size()))
    stop("likelihoodAD: xreg has fewer rows than data");

  const std::vector<double> th(par.begin(), par.end());
  double val;
  std::vector<double> g, H;
  ad_eval(model, th, deriv, val, g, H);

  const int n = static_cast<int>(th.size());
  NumericVector grad(g.begin(), g.end());
  if (deriv < 2) {
    return List::create(_["value"] = val, _["gradient"] = grad);
  }
  NumericMatrix hess(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) hess(i, j) = H[static_cast<std::size_t>(i) * n + j];
  return List::create(_["value"] = val, _["gradient"] = grad, _["hessian"] = hess);
}
