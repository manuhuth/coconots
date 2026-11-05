// Forward-mode automatic differentiation with dual numbers.
//
// Dual<T, N> carries a value v and N partial derivatives d[0..N-1]; T is
// double for first derivatives and Dual<double, N> for second derivatives
// (nested duals, the same construction as ForwardDiff.jl's Hessian).
// Only the operations used by the coconots likelihoods are provided.

#ifndef COCONOTS_COCO_AD_H
#define COCONOTS_COCO_AD_H

#include <cmath>
#include <type_traits>

namespace coco_ad {

template <class T, int N>
struct Dual {
  T v;
  T d[N];

  // d() value-initializes the partials (zeros)
  Dual() : v(0.0), d() {}
  // constants (zero partials)
  template <class A, typename std::enable_if<std::is_arithmetic<A>::value, int>::type = 0>
  Dual(A x) : v(static_cast<double>(x)), d() {}
};

inline double value_of(double x) { return x; }

template <class T, int N>
inline double value_of(const Dual<T, N>& x) { return value_of(x.v); }

// ------------------------------------------------------------- arithmetic

template <class T, int N>
inline Dual<T, N> operator+(const Dual<T, N>& a) { return a; }

template <class T, int N>
inline Dual<T, N> operator-(const Dual<T, N>& a) {
  Dual<T, N> r;
  r.v = -a.v;
  for (int i = 0; i < N; ++i) r.d[i] = -a.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> operator+(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> r;
  r.v = a.v + b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] + b.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> operator+(const Dual<T, N>& a, double b) {
  Dual<T, N> r = a;
  r.v = a.v + b;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator+(double a, const Dual<T, N>& b) {
  Dual<T, N> r = b;
  r.v = a + b.v;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator-(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> r;
  r.v = a.v - b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] - b.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> operator-(const Dual<T, N>& a, double b) {
  Dual<T, N> r = a;
  r.v = a.v - b;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator-(double a, const Dual<T, N>& b) {
  Dual<T, N> r;
  r.v = a - b.v;
  for (int i = 0; i < N; ++i) r.d[i] = -b.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> operator*(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> r;
  r.v = a.v * b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.v * b.d[i] + a.d[i] * b.v;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator*(const Dual<T, N>& a, double b) {
  Dual<T, N> r;
  r.v = a.v * b;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * b;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator*(double a, const Dual<T, N>& b) {
  Dual<T, N> r;
  r.v = a * b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a * b.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> operator/(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> r;
  const T inv = 1.0 / b.v;
  r.v = a.v * inv;
  for (int i = 0; i < N; ++i) r.d[i] = (a.d[i] - r.v * b.d[i]) * inv;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator/(const Dual<T, N>& a, double b) {
  Dual<T, N> r;
  r.v = a.v / b;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] / b;
  return r;
}

template <class T, int N>
inline Dual<T, N> operator/(double a, const Dual<T, N>& b) {
  Dual<T, N> r;
  const T inv = 1.0 / b.v;
  r.v = a * inv;
  const T k = -(r.v * inv);
  for (int i = 0; i < N; ++i) r.d[i] = k * b.d[i];
  return r;
}

// ------------------------------------------------------------- comparisons
// (on the value; a branch contributes the derivative of the branch taken)

template <class T, int N>
inline bool operator<(const Dual<T, N>& a, double b) { return value_of(a) < b; }
template <class T, int N>
inline bool operator<=(const Dual<T, N>& a, double b) { return value_of(a) <= b; }
template <class T, int N>
inline bool operator>(const Dual<T, N>& a, double b) { return value_of(a) > b; }
template <class T, int N>
inline bool operator>=(const Dual<T, N>& a, double b) { return value_of(a) >= b; }

// ------------------------------------------------------------- functions

template <class T, int N>
inline Dual<T, N> exp(const Dual<T, N>& a) {
  using std::exp;
  Dual<T, N> r;
  r.v = exp(a.v);
  for (int i = 0; i < N; ++i) r.d[i] = r.v * a.d[i];
  return r;
}

template <class T, int N>
inline Dual<T, N> log(const Dual<T, N>& a) {
  using std::log;
  Dual<T, N> r;
  r.v = log(a.v);
  const T inv = 1.0 / a.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * inv;
  return r;
}

// a^n for an integer constant n; d/da a^0 = 0 also where a = 0.
template <class T, int N>
inline Dual<T, N> pow(const Dual<T, N>& a, int n) {
  using std::pow;
  Dual<T, N> r;
  if (n == 0) {
    r.v = pow(a.v, 0);
    return r;
  }
  r.v = pow(a.v, n);
  const T k = static_cast<double>(n) * pow(a.v, n - 1);
  for (int i = 0; i < N; ++i) r.d[i] = k * a.d[i];
  return r;
}

// a^e for a real constant e.
template <class T, int N>
inline Dual<T, N> pow(const Dual<T, N>& a, double e) {
  using std::pow;
  Dual<T, N> r;
  if (e == 0.0) {
    r.v = pow(a.v, 0.0);
    return r;
  }
  r.v = pow(a.v, e);
  const T k = e * pow(a.v, e - 1.0);
  for (int i = 0; i < N; ++i) r.d[i] = k * a.d[i];
  return r;
}

} // namespace coco_ad

#endif
