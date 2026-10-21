#pragma once
// Double-double arithmetic (about 106 significant bits) for accumulations that R carries in long double.
// `long double` is 80-bit on x86-64, 128-bit on aarch64 Linux and plain double on macOS arm64 and MSVC, so it would
// make results depend on the platform; double-double gives the same bits everywhere. It relies on IEEE double
// rounding and no fused multiply-add contraction (-ffp-contract=off).
// Error-free transformations: Knuth TwoSum, Dekker's split product (Dekker 1971); algorithms as in Bailey's QD
// library (Hida, Li & Bailey 2001).
#include <cmath>

namespace rcli {

struct DD {
  double hi = 0.0, lo = 0.0;
  DD() = default;
  DD(double h) : hi(h) {}  // NOLINT: implicit widening, like double -> long double
  DD(double h, double l) : hi(h), lo(l) {}
  explicit operator double() const { return hi + lo; }
};

namespace dd_detail {
inline DD fast_two_sum(double a, double b) {  // |a| >= |b|
  const double s = a + b;
  return {s, b - (s - a)};
}
inline DD two_sum(double a, double b) {
  const double s = a + b;
  const double bb = s - a;
  return {s, (a - (s - bb)) + (b - bb)};
}
inline DD two_prod(double a, double b) {
  const double p = a * b;
  if (!std::isfinite(p)) return {p, 0.0};
  constexpr double kSplit = 134217729.0;  // 2^27 + 1
  const double ta = kSplit * a, tb = kSplit * b;
  if (!std::isfinite(ta) || !std::isfinite(tb)) return {p, 0.0};  // too large to split: plain product
  const double ah = ta - (ta - a), al = a - ah;
  const double bh = tb - (tb - b), bl = b - bh;
  return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
}
inline DD norm(double s, double e) { return std::isfinite(s) ? fast_two_sum(s, e) : DD{s, 0.0}; }
}  // namespace dd_detail

inline DD operator+(const DD& a, const DD& b) {
  const DD s = dd_detail::two_sum(a.hi, b.hi);
  if (!std::isfinite(s.hi)) return {s.hi, 0.0};
  const DD t = dd_detail::two_sum(a.lo, b.lo);
  const DD u = dd_detail::fast_two_sum(s.hi, s.lo + t.hi);
  return dd_detail::norm(u.hi, u.lo + t.lo);
}
inline DD operator-(const DD& a) { return {-a.hi, -a.lo}; }
inline DD operator-(const DD& a, const DD& b) { return a + (-b); }
inline DD operator*(const DD& a, const DD& b) {
  const DD p = dd_detail::two_prod(a.hi, b.hi);
  if (!std::isfinite(p.hi)) return {p.hi, 0.0};
  return dd_detail::norm(p.hi, p.lo + (a.hi * b.lo + a.lo * b.hi));
}
inline DD operator/(const DD& a, const DD& b) {
  const double q1 = a.hi / b.hi;
  if (!std::isfinite(q1)) return {q1, 0.0};
  const DD r = a - b * DD(q1);
  const double q2 = r.hi / b.hi;
  const DD r2 = r - b * DD(q2);
  const double q3 = r2.hi / b.hi;
  return dd_detail::fast_two_sum(q1, q2) + DD(q3);
}
inline DD& operator+=(DD& a, const DD& b) { return a = a + b; }
inline bool operator==(const DD& a, const DD& b) { return a.hi == b.hi && a.lo == b.lo; }

inline DD sqrt(const DD& a) {
  if (!(a.hi > 0.0) || !std::isfinite(a.hi)) return {std::sqrt(a.hi), 0.0};
  const double x = std::sqrt(a.hi);
  const DD r = a - dd_detail::two_prod(x, x);
  return dd_detail::fast_two_sum(x, r.hi / (2.0 * x));
}
inline DD fabs(const DD& a) { return a.hi < 0.0 ? -a : a; }

}  // namespace rcli
