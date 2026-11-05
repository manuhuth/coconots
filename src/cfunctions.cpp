# include <numeric>
# include <algorithm>
# include <vector>
# include <cmath>
# include <unordered_map>
# include <atomic>
# include <thread>
# include <Rcpp.h>

// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp ;

namespace {

// Factorials up to 170! are exactly representable in double; larger values
// overflow to Inf just like the former recursive implementation
// (fact(i) = i * fact(i-1)), whose multiplication order the table reproduces.
struct FactTable {
  static const int size = 180;
  double v[size];
  FactTable() {
    v[0] = 1.0;
    v[1] = 1.0;
    for (int i = 2; i < size; ++i) v[i] = i * v[i - 1];
  }
};
const FactTable FACT_TABLE;

inline double factorial(int i) {
  if (i < 2) return 1.0;
  if (i < FactTable::size) return FACT_TABLE.v[i];
  double f = FACT_TABLE.v[FactTable::size - 1];
  for (int k = FactTable::size; k <= i; ++k) f = k * f;
  return f;
}

enum LinkFunction { LINK_LOG, LINK_IDENTITY, LINK_RELU, LINK_SOFTPLUS };

LinkFunction parseLinkFunction(const std::string& link_function) {
  if (link_function == "log") return LINK_LOG;
  if (link_function == "identity") return LINK_IDENTITY;
  if (link_function == "relu") return LINK_RELU;
  if (link_function == "softplus") return LINK_SOFTPLUS;
  stop("Unknown link function");
}

inline double applyLink(double mathelp, LinkFunction link) {
  switch (link) {
  case LINK_LOG:      return exp(mathelp);
  case LINK_IDENTITY: return mathelp;
  case LINK_RELU:     return (mathelp > 0) ? mathelp : 1e-10;
  default:            return (mathelp > 0) ? mathelp + log(1.0 + exp(-mathelp))
                                           : log(1.0 + exp(mathelp));
  }
}

// Accumulates in the same order as std::inner_product over a copied row,
// but without allocating the row.
inline double rowInnerProduct(const NumericMatrix& xreg, int row,
                              const NumericVector& lambdas) {
  double acc = 0.0;
  const int n = xreg.ncol();
  for (int j = 0; j < n; ++j) acc = acc + xreg(row, j) * lambdas[j];
  return acc;
}

// Largest admissible lookup-table index. Count data anywhere near this bound
// makes the model degenerate (factorials overflow at 171), so refusing is
// safer than attempting the corresponding gigabyte-sized allocations.
const int MAX_TABLE_INDEX = 1000000;

inline void check_table_index(int n) {
  if (n > MAX_TABLE_INDEX)
    stop("Count values are too large for the (Generalized) Poisson recursions.");
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
  auto work = [&](int tid) {
    try {
      int i;
      while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n) f(tid, i);
    } catch (...) {
      failed.store(true);
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nthreads - 1);
  for (int k = 1; k < nthreads; ++k) pool.emplace_back(work, k);
  work(0);
  for (std::size_t k = 0; k < pool.size(); ++k) pool[k].join();
  if (failed.load()) stop("coconots: a worker thread failed (out of memory?)");
}

// ------------------------------------------------------------------ order 1
//
// P(x | y) = sum_r  QuasiBinomial(r; y) * GP(x - r; lambda, eta)
//
// The r-invariant factors pow(alpha+psi*r, r-1), pow(1-alpha+psi*k, k-1),
// pow(1+psi*y, y-1) and the GP innovation terms are tabulated once per
// parameter value; each table entry equals the pow()/exp() call it replaces.

struct GP1Tables {
  double lambda, alpha, eta, psi;
  std::vector<double> powA, powB, denom, powL, expL;
  int cap;

  void init(double lambda_, double alpha_, double eta_, double psi_) {
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

// Transition density P(x | y); term order and operation chain match the
// original inline expression exactly.
inline double gp1_density(int x, int y, const GP1Tables& tb) {
  const double alpha = tb.alpha, lambda = tb.lambda;
  const int minxy = std::min(x, y);
  double nain = 0.0;
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

struct GP2Params {
  double lambda, alpha1, alpha2, alpha3, eta;
  double U, beta1, beta2, beta3, zeta, lu13, lu13a;

  void init(double lambda_, double alpha1_, double alpha2_, double alpha3_, double eta_) {
    lambda = lambda_; alpha1 = alpha1_; alpha2 = alpha2_; alpha3 = alpha3_; eta = eta_;
    U = 1 / (1 - alpha1 - alpha2 - alpha3);
    beta1 = lambda * U * alpha1;
    beta2 = lambda * U * alpha2;
    beta3 = lambda * U * alpha3;
    zeta = lambda * U * (1 - 2 * alpha1 - alpha3);
    lu13 = lambda * U * (1 - alpha1 - alpha3);
    lu13a = lambda * U * (alpha1 + alpha3);
  }
};

struct GP2Tables {
  GP2Params p;
  std::vector<double> powB1, expB1, powB2, expB2, powB3, expB3, powL, expL, powZ, expZ;
  int cap;

  void init(const GP2Params& p_) {
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

// Bivariate GP normalizer P(y, z); r-invariant.
inline double gp2_pyz(int y, int z, const GP2Params& p) {
  const double help = (p.beta1 + p.beta3) * pow(p.lu13, 2) *
    exp(-(p.beta1 + p.beta3) - 2 * p.lu13 - y * p.eta - z * p.eta);
  const int minyz = std::min(y, z);
  double dwarf = 0.0;
  for (int j = 0; j <= minyz; ++j) {
    dwarf = dwarf + pow(p.lu13 + p.eta * (y - j), y - j - 1) / factorial(y - j) *
      pow(p.lu13 + p.eta * (z - j), z - j - 1) / factorial(z - j) *
      pow(p.lu13a + p.eta * j, j - 1) / factorial(j) * exp(j * p.eta);
  }
  return dwarf * help;
}

// Triple convolution over (s, v, w). The loop bounds enumerate exactly the
// index triples that satisfy r-s-v >= 0, z-r+v-w >= 0 and y-s-v-w >= 0, in
// the same lexicographic order the former full 0..y cube visited them.
// Tables must cover indices up to max(y, z) (all six factor indices are
// bounded by y or z on the feasible set).
inline double gp2_zwerg(int r, int y, int z, const GP2Params& p, const GP2Tables& tb) {
  double zwerg = 0.0;
  const int smax = std::min(y, r);
  for (int s = 0; s <= smax; ++s) {
    const int vlo = std::max(0, r - z);
    const int vhi = std::min(r - s, y - s);
    for (int v = vlo; v <= vhi; ++v) {
      const int k = r - s - v;
      const int wmax = std::min(z - r + v, y - s - v);
      // leading subchain of the product; hoisting a prefix of a
      // left-associative chain does not change the arithmetic
      const double pre = p.beta3 * tb.powB3[s] / factorial(s) * tb.expB3[s] *
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
struct GVec {
  std::vector<double> g;
  double tail;
};

// Tables must already cover max(y, z, 0); with const tables this is safe to
// run concurrently for distinct `out` targets.
void gp2_build_g(int y, int z, const GP2Params& p, const GP2Tables& tb, GVec& out) {
  const double pyz = gp2_pyz(y, z, p);
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
struct GP2Conv {
  std::vector<double> q1, q2, q3, qL, qZ, W, zw;
  int mlo, mhi, nlo, nhi;
  std::size_t stride;

  // pmf tables for indices 0..n; tb must cover n.
  void build_pmf(const GP2Params& p, const GP2Tables& tb, int n) {
    q1.resize(n + 1); q2.resize(n + 1); q3.resize(n + 1);
    qL.resize(n + 1); qZ.resize(n + 1);
    for (int i = 0; i <= n; ++i) {
      const double f = factorial(i);
      q1[i] = p.beta1 * tb.powB1[i] / f * tb.expB1[i];
      q2[i] = p.beta2 * tb.powB2[i] / f * tb.expB2[i];
      q3[i] = p.beta3 * tb.powB3[i] / f * tb.expB3[i];
      qL[i] = p.lambda * tb.powL[i] / f * tb.expL[i];
      qZ[i] = p.zeta * tb.powZ[i] / f * tb.expZ[i];
    }
  }

  // W(m, n) for m in [mlo_, mhi_], n in [nlo_, nhi_]; pmf tables must cover
  // max(mhi_, nhi_).
  void build_W(int mlo_, int mhi_, int nlo_, int nhi_) {
    mlo = mlo_; mhi = mhi_; nlo = nlo_; nhi = nhi_;
    stride = static_cast<std::size_t>(nhi - nlo + 1);
    W.resize(static_cast<std::size_t>(mhi - mlo + 1) * stride);
    for (int m = mlo; m <= mhi; ++m) {
      double* row = &W[static_cast<std::size_t>(m - mlo) * stride];
      for (int n = nlo; n <= nhi; ++n) {
        const int wmax = std::min(m, n);
        double acc = 0.0;
        for (int w = 0; w <= wmax; ++w) acc = acc + q1[w] * qL[m - w] * qZ[n - w];
        row[n - nlo] = acc;
      }
    }
  }

  inline double Wat(int m, int n) const {
    return W[static_cast<std::size_t>(m - mlo) * stride + static_cast<std::size_t>(n - nlo)];
  }

  // out[r] = zwerg(r | y, z) for r = 0..rcap, with y, z >= 0, 0 <= rcap <= y+z.
  // W must cover m in [max(0, z-rcap), z] and n in [max(0, y-rcap), y].
  void zwerg_upto(int y, int z, int rcap, std::vector<double>& out) const {
    out.assign(rcap + 1, 0.0);
    const int smax = std::min(std::min(y, z), rcap);
    for (int s = 0; s <= smax; ++s) {
      const double a3 = q3[s];
      const int vmax = std::min(y - s, rcap - s);
      for (int v = 0; v <= vmax; ++v) {
        const double a31 = a3 * q1[v];
        const int n = y - s - v;
        const int kmax = std::min(z - s, rcap - s - v);
        double* o = &out[s + v];
        for (int k = 0; k <= kmax; ++k) {
          o[k] = o[k] + a31 * q2[k] * Wat(z - s - k, n);
        }
      }
    }
  }
};

// Upper bound on the number of W entries tabulated at once; beyond it the
// direct gp2_zwerg evaluation is used.
const double GP2_CONV_MAX_ENTRIES = 4e6;

// gp2_build_g via GP2Conv; cv must hold pmf tables and W covering
// m in [0, z], n in [0, y]. Safe to run concurrently (cv is read-only).
void gp2_build_g_conv(int y, int z, const GP2Params& p, const GP2Conv& cv, GVec& out) {
  const double pyz = gp2_pyz(y, z, p);
  cv.zwerg_upto(y, z, y + z, out.g);
  const double ipyz = pow(pyz, -1);
  for (std::size_t r = 0; r < out.g.size(); ++r) out.g[r] = out.g[r] * ipyz;
  out.tail = 0.0 * pow(pyz, -1);
}

typedef std::unordered_map<long long, GVec> GCache;

inline long long g_key(int y, int z) {
  return (static_cast<long long>(y) << 32) | static_cast<unsigned int>(z);
}

inline GVec& g_lookup(GCache& cache, int y, int z, const GP2Params& p, GP2Tables& tb) {
  GVec& gv = cache[g_key(y, z)];
  if (gv.g.empty()) {
    tb.ensure(std::max(std::max(y, z), 0));
    gp2_build_g(y, z, p, tb, gv);
  }
  return gv;
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
  GP2Params p;
  p.init(lambda, alpha1, alpha2, alpha3, eta);
  GP2Tables tb;
  tb.init(p);

  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  if (maxval >= 0) tb.ensure(maxval);

  GCache cache;

  // Register the distinct (y, z) histories first, then build their g-vectors
  // concurrently: each build is independent and the tables are read-only once
  // ensured up to maxval. Pointers into the node-based unordered_map stay
  // valid across later insertions. The serial accumulation below is unchanged,
  // so results are bit-identical to the single-threaded evaluation.
  struct Job { int y, z; GVec* gv; };
  std::vector<Job> jobs;
  int maxy = 0, maxz = 0;
  for (int t = seas2 + 1; t <= T; ++t) {
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];
    if (y >= 0 && z >= 0) {
      std::pair<GCache::iterator, bool> ins = cache.emplace(g_key(y, z), GVec());
      if (ins.second) {
        jobs.push_back(Job{y, z, &ins.first->second});
        maxy = std::max(maxy, y);
        maxz = std::max(maxz, z);
      }
    }
  }

  // lambda is constant, so the pmf tables and W(m, n) are shared by all
  // histories; build them once (serially) before the read-only g builds.
  const bool use_conv = !jobs.empty() &&
    static_cast<double>(maxz + 1) * (maxy + 1) <= GP2_CONV_MAX_ENTRIES;
  GP2Conv cv;
  if (use_conv) {
    cv.build_pmf(p, tb, maxval);
    cv.build_W(0, maxz, 0, maxy);
  }

  const int nthreads = resolve_threads(static_cast<int>(jobs.size()), cores);
  parallel_for(static_cast<int>(jobs.size()), nthreads, [&](int, int i) {
    if (use_conv) {
      gp2_build_g_conv(jobs[i].y, jobs[i].z, p, cv, *jobs[i].gv);
    } else {
      gp2_build_g(jobs[i].y, jobs[i].z, p, tb, *jobs[i].gv);
    }
  });

  double mlef = 0.0;
  for (int t = seas2 + 1; t <= T; ++t) {
    const int x = data[t - 1];
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];

    double nain = 0.0;
    if (y >= 0 && z >= 0) {
      const GVec& gv = g_lookup(cache, y, z, p, tb);
      const int minxyz = std::min(x, y + z);
      for (int r = 0; r <= minxyz; ++r) {
        nain = nain + gv.g[r] * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
      }
    } else {
      // degenerate inputs (e.g. NA counts): zwerg/pyz vanish just like in the
      // guarded cube loops, propagating the same 0 * Inf arithmetic
      const double pyz = gp2_pyz(y, z, p);
      const int minxyz = std::min(x, y + z);
      for (int r = 0; r <= minxyz; ++r) {
        const double gryz = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
        nain = nain + gryz * lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
      }
    }
    mlef = mlef + log(nain);
  }

  return (-mlef);
}


// [[Rcpp::export]]
double likelihoodGP1(double lambda, double alpha, double eta, int T, int seas,
                     const std::vector<int>& data)
{
  const double psi = eta * (1 - alpha) / lambda;
  GP1Tables tb;
  tb.init(lambda, alpha, eta, psi);

  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  if (maxval >= 0) tb.ensure(maxval);

  // conditional densities depend on t only through (x, y); cache them
  const bool use_memo = maxval >= 0 &&
    static_cast<double>(maxval + 1) * (maxval + 1) <= 4e6;
  std::vector<double> memo;
  if (use_memo) memo.assign(static_cast<std::size_t>(maxval + 1) * (maxval + 1), -1.0);

  double mlef = 0.0;
  for (int t = seas + 1; t <= T; ++t) {
    const int x = data[t - 1];
    const int y = data[t - seas - 1];

    double p;
    if (use_memo && x >= 0 && y >= 0) {
      const std::size_t key = static_cast<std::size_t>(x) * (maxval + 1) + y;
      p = memo[key];
      if (p < 0.0) {
        p = gp1_density(x, y, tb);
        memo[key] = p;
      }
    } else {
      p = gp1_density(x, y, tb);
    }
    mlef = mlef + log(p);
  }

  return (-mlef);
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

    int r = 0;
    double nain = 0.0;
    if (y >= 0 && z >= 0) {
      const GVec& gv = g_lookup(cache, y, z, p, tb);
      const int rmax = y + z;
      while (nain <= unif) {
        nain = nain + ((r <= rmax) ? gv.g[r] : gv.tail);
        r = r + 1;
      }
    } else {
      const double pyz = gp2_pyz(y, z, p);
      while (nain <= unif) {
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

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas - 1];
    if (y >= 0) tb.ensure(y);

    const double unif = uniform[t - 1];

    int r = 0;
    double nain = 0.0;
    while (nain <= unif) {
      if (y >= r) {
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

  double nain = 0.0;
  for (int r = 0; r <= minxy; ++r) {
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
// reaches 1 - val_num, mirroring the accumulation loop in cocoResid().
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
  std::vector<double> probs;
  double acc = 0.0;
  int r = 0;
  while (acc < 1 - val_num) {
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
    p.init(lambda, alpha1, alpha2, alpha3, eta);
    tb.init(p);
    tb.ensure(std::max(std::max(y, z), 0));

    const double pyz = gp2_pyz(y, z, p);
    const double unif = uniform[t - N - 1];

    int r = 0;
    double nain = 0.0;
    while (nain <= unif) {
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

  for (int t = N + 1; t <= T; ++t) {
    const int y = data[t - seas - 1];

    const double mathelp = rowInnerProduct(xreg, t - 1 - N, lambdas);
    const double lambda = applyLink(mathelp, link);
    const double psi = eta * (1 - alpha) / lambda;
    const double denom = (y >= 0) ? pow(1 + psi * y, y - 1) : 0.0;

    const double unif = uniform[t - 1];
    int r = 0;
    double nain = 0.0;
    while (nain <= unif) {
      if (y >= r) {
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

} // namespace


// Inversion sampler of the Generalized Poisson distribution: for each
// uniform u[i], the smallest x with cumulative pmf (accumulated from x = 0
// upwards, as in the former R loop) >= u[i], at lambda = lambdas[i].
// The cumulative table is reused while consecutive lambdas are equal.
// [[Rcpp::export]]
NumericVector rgenpoisInv(NumericVector u, NumericVector lambdas, double eta)
{
  const int n = u.size();
  if (lambdas.size() != n) stop("u and lambdas must have the same length");
  NumericVector out(n);
  std::vector<double> cdf;
  double cur = 0.0;
  bool have = false;
  long long iter = 0;

  for (int i = 0; i < n; ++i) {
    const double lam = lambdas[i];
    if (!have || !(lam == cur)) {
      cdf.clear();
      cdf.push_back(hmmpa_dgenpois(0, lam, eta));
      cur = lam;
      have = true;
    }
    const double ui = u[i];
    int idx = 0;
    while (true) {
      if (idx >= static_cast<int>(cdf.size())) {
        const double next = cdf.back() + hmmpa_dgenpois(idx, cur, eta);
        cdf.push_back(next);
      }
      const double last = cdf[idx];
      if (std::isnan(last) || std::isnan(ui))
        stop("missing value where TRUE/FALSE needed");
      if (!(last < ui)) break;
      idx = idx + 1;
      if ((++iter & 0xFFFFF) == 0) checkUserInterrupt();
    }
    out[i] = idx;
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

// Thinning draw of order 1 without covariates (simGP1 loop); terms vanish
// beyond r = y, where the former loop never terminated if c[y] <= u.
inline int gp1_cdf_draw(LazyCdf& e, int y, double u, const GP1Tables& tb) {
  e.ready = true;
  for (int j = 0; ; ++j) {
    if (j == static_cast<int>(e.c.size())) {
      if (j > y) return y;
      const double prev = (j == 0) ? 0.0 : e.c[j - 1];
      e.c.push_back(prev + factorial(y) / factorial(y - j) / factorial(j) * tb.alpha *
                    (1 - tb.alpha) * tb.powA[j] * tb.powB[y - j] / tb.denom[y]);
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
      e.c.push_back(prev + factorial(y) / factorial(y - j) / factorial(j) * alpha * (1 - alpha) *
                    pow(alpha + psi * j, j - 1) * pow(1 - alpha + psi * (y - j), y - j - 1) / e.aux);
    }
    if (!(e.c[j] <= u)) return j;
  }
}

// Inversion of the Generalized Poisson CDF exactly as in rgenpoisInv, with
// the cumulative table kept for one lambda.
struct GPInvCdf {
  double lam, eta;
  std::vector<double> cdf;
  long long iter;

  GPInvCdf() : lam(0.0), eta(0.0), iter(0) {}

  void reset(double lam_, double eta_) {
    lam = lam_;
    eta = eta_;
    cdf.clear();
    cdf.push_back(hmmpa_dgenpois(0, lam, eta));
  }

  int draw(double u) {
    int idx = 0;
    while (true) {
      if (idx >= static_cast<int>(cdf.size())) {
        const double next = cdf.back() + hmmpa_dgenpois(idx, lam, eta);
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
  xd.resize(n);
  for (int i = 0; i < n; ++i) xd[i] = static_cast<double>(x[i]) - mean;
  double acov0 = 0.0;
  for (int i = 0; i < n; ++i) acov0 += xd[i] * xd[i];
  acov0 = acov0 / n;
  const double se = sqrt(acov0);
  const double den = se * se;
  const int lagmax = std::min(nlags, n - 1);
  for (int lag = 1; lag <= lagmax; ++lag) {
    double sum = 0.0;
    for (int i = 0; i < n - lag; ++i) sum += xd[i + lag] * xd[i];
    const double a = (sum / n) / den;
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

    sample_acf(&data[size - n_keep], n_keep, nlags, xd,
               out.begin() + static_cast<std::size_t>(b) * nlags);
    checkUserInterrupt();
  }
  return out;
}


// Bootstrap ACFs with covariates. Each replicate is the cocoSim_cov series
// of length `size` without burn-in: initial values from init_lambdas,
// innovations from lambdas (one per xreg row), n_unif thinning uniforms,
// and the thinning at step t uses lambda from xreg row t - order - 1, the
// innovation t - order - 1 and the uniform t - order - 1 (order 2) or
// t - 1 (order 1), as in simGP2cov / simGP1cov. Replicates are drawn in
// order (so the RNG stream is that of consecutive cocoSim calls) and then
// simulated time-major in blocks, sharing per-t tables and caches.
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
  if (xreg.nrow() < size - N || n_innov < size - N) stop("xreg has too few rows");
  if (n_unif < ((order == 2) ? size - N : size)) stop("too few thinning uniforms");

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
  CdfStore store;
  GPInvCdf inv, inv_init;

  for (int b0 = 0; b0 < nB; b0 += blk) {
    const int nb = std::min(blk, nB - b0);

    // random numbers per replicate, in the order of cocoSim_cov
    for (int k = 0; k < nb; ++k) {
      int* d = &data[static_cast<std::size_t>(k) * size];
      for (int i = 0; i < N; ++i) {
        if (order == 1 && gp) {
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
      const int row = t - N - 1;
      const int ui = (order == 2) ? (t - N - 1) : (t - 1);
      const double lam = applyLink(rowInnerProduct(xreg, row, betas), link);
      const double psi = eta * (1 - alpha1) / lam;
      if (order == 2) {
        p2.init(lam, alpha1, alpha2, alpha3, eta);
        tb2.init(p2);
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
          tb2.ensure(std::max(y, z));
          LazyCdf* e = store.get(g_key(y, z));
          const std::size_t before = e->c.size();
          r = gp2_cdf_draw(*e, y, z, u, p2, tb2);
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



// [[Rcpp::export]]
double likelihoodGP2cov(double alpha1, double alpha2, double alpha3,
                        double eta, NumericVector lambdas, int T, int seas1, int seas2,
                        const std::vector<int>& data, NumericMatrix xreg, std::string link_function,
                        int cores = 1)
{
  const LinkFunction link = parseLinkFunction(link_function);

  // lambda_t varies with t, so nothing can be cached across observations --
  // but each observation's contribution is independent. Compute them
  // concurrently with thread-private tables and accumulate serially in t
  // order afterwards, keeping results bit-identical to the serial loop.
  // Bound the table sizes once here so workers never raise an R error.
  int maxval = -1;
  for (std::size_t i = 0; i < data.size(); ++i) maxval = std::max(maxval, data[i]);
  check_table_index(maxval);

  const int nterms = T - seas2;
  if (nterms <= 0) return 0.0;
  std::vector<double> logn(nterms);

  const int nthreads = resolve_threads(nterms, cores);
  std::vector<GP2Params> pv(nthreads);
  std::vector<GP2Tables> tbv(nthreads);
  std::vector<GP2Conv> cvv(nthreads);

  parallel_for(nterms, nthreads, [&](int tid, int i) {
    const int t = seas2 + 1 + i;
    const int x = data[t - 1];
    const int y = data[t - seas1 - 1];
    const int z = data[t - seas2 - 1];

    const double mathelp = rowInnerProduct(xreg, t - 1, lambdas);
    const double lambda = applyLink(mathelp, link);
    GP2Params& p = pv[tid];
    GP2Tables& tb = tbv[tid];
    p.init(lambda, alpha1, alpha2, alpha3, eta);
    tb.init(p);
    tb.ensure(std::max(std::max(std::max(x, y), z), 0));

    const double pyz = gp2_pyz(y, z, p);
    const int minxyz = std::min(x, y + z);
    double nain = 0.0;
    if (x >= 0 && y >= 0 && z >= 0 &&
        static_cast<double>(y + 1) * (z + 1) <= GP2_CONV_MAX_ENTRIES) {
      // only r <= minxyz is needed, so W is restricted to the (m, n) that
      // these r reach
      GP2Conv& cv = cvv[tid];
      cv.build_pmf(p, tb, std::max(y, z));
      cv.build_W(std::max(0, z - minxyz), z, std::max(0, y - minxyz), y);
      cv.zwerg_upto(y, z, minxyz, cv.zw);
      for (int r = 0; r <= minxyz; ++r) {
        const double gryz = cv.zw[r] * pow(pyz, -1);
        nain = nain + gryz * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
      }
    } else {
      for (int r = 0; r <= minxyz; ++r) {
        const double gryz = gp2_zwerg(r, y, z, p, tb) * pow(pyz, -1);
        nain = nain + gryz * lambda * tb.powL[x - r] * tb.expL[x - r] / factorial(x - r);
      }
    }

    if (nain <= 0) {
      nain = 1e-12;
    }
    logn[i] = log(nain);
  });

  double mlef = 0.0;
  for (int i = 0; i < nterms; ++i) mlef = mlef + logn[i];

  return (-mlef);
}


// [[Rcpp::export]]
double likelihoodGP1cov(double alpha, double eta, NumericVector lambdas, int T, int seas,
                        const std::vector<int>& data, NumericMatrix xreg, std::string link_function)
{
  const LinkFunction link = parseLinkFunction(link_function);

  double mlef = 0.0;
  for (int t = seas + 1; t <= T; ++t) {
    const int x = data[t - 1];
    const int y = data[t - seas - 1];

    const double mathelp = rowInnerProduct(xreg, t - 1, lambdas);
    const double lambda = applyLink(mathelp, link);
    const double psi = eta * (1 - alpha) / lambda;
    const double denom = (y >= 0) ? pow(1 + psi * y, y - 1) : 0.0;

    const int minxy = std::min(x, y);
    double nain = 0.0;
    for (int r = 0; r <= minxy; ++r) {
      nain = nain + factorial(y) / factorial(y - r) / factorial(r) * alpha * (1 - alpha) *
        pow(alpha + psi * r, r - 1) * pow(1 - alpha + psi * (y - r), y - r - 1) / denom *
        lambda * pow(lambda + (x - r) * eta, x - r - 1) * exp(-lambda - (x - r) * eta) / factorial(x - r);
    }
    mlef = mlef + log(nain);
  }

  return (-mlef);
}
