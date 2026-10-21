// Interpolating cubic splines and piecewise linear/constant interpolation (clean-room).
// References:
//   Forsythe, G. E., Malcolm, M. A. and Moler, C. B. (1977) Computer Methods for Mathematical Computations, Wiley,
//     sec. 4.4: the end conditions make the spline's third derivative on the first (last) interval equal to that of
//     the cubic through the first (last) four points.
//   de Boor, C. (1978) A Practical Guide to Splines, Springer, ch. IV (second-derivative form, natural and periodic
//     end conditions).
//   Fritsch, F. N. and Carlson, R. E. (1980) Monotone piecewise cubic interpolation. SIAM J. Numer. Anal. 17, 238-246.
// Representation: on [x_i, x_{i+1}) the spline is y_i + b_i t + c_i t^2 + d_i t^3 with t = x - x_i; the coefficients
// stored for the last knot describe the extrapolation to the right.
#include "rcli/interp.hpp"

#include "rcli/dd.hpp"

#include <algorithm>
#include <numeric>

namespace rcli {
namespace {

// Sort by x and collapse tied x by the mean of their y (R's ties = mean).
void collapse_ties(const Vec& x, const Vec& y, Vec& ux, Vec& uy, const char* fn) {
  if (x.size() != y.size()) throw Error(std::string(fn) + ": x and y lengths differ");
  for (size_t i = 0; i < x.size(); ++i)
    if (std::isnan(x[i]) || std::isnan(y[i])) throw Error(std::string(fn) + ": missing values (NaN) are not allowed");
  std::vector<size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  ux.clear();
  uy.clear();
  for (size_t i = 0; i < idx.size();) {
    size_t j = i;
    DD s;  // double-double where R's mean() uses long double (same bits on every platform)
    while (j < idx.size() && x[idx[j]] == x[idx[i]]) s += y[idx[j++]];
    const DD cnt = static_cast<double>(j - i);
    DD m = s / cnt;
    if (j - i > 1) {  // second pass as R's mean()
      DD t;
      for (size_t k = i; k < j; ++k) t += DD(y[idx[k]]) - m;
      m += t / cnt;
    }
    ux.push_back(x[idx[i]]);
    uy.push_back(static_cast<double>(m));
    i = j;
  }
}

// Solve a tridiagonal system (sub a, diag b, super c, rhs r; a[0], c[n-1] unused) by elimination without pivoting;
// callers only pass diagonally dominant systems.
Vec thomas(Vec a, Vec b, Vec c, Vec r) {
  const size_t n = b.size();
  for (size_t i = 1; i < n; ++i) {
    const double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    r[i] -= m * r[i - 1];
  }
  Vec x(n);
  x[n - 1] = r[n - 1] / b[n - 1];
  for (size_t i = n - 1; i-- > 0;) x[i] = (r[i] - c[i] * x[i + 1]) / b[i];
  return x;
}

// Cyclic tridiagonal: row i couples x[i-1], x[i], x[i+1] with wrap-around (Sherman-Morrison for n >= 3).
Vec cyclic(const Vec& a, const Vec& b, const Vec& c, const Vec& r) {
  const size_t n = b.size();
  if (n == 1) return {r[0] / (a[0] + b[0] + c[0])};
  if (n == 2) {  // [b0, a0 + c0; a1 + c1, b1]
    const double p = b[0], q = a[0] + c[0], s = a[1] + c[1], t = b[1];
    const double det = p * t - q * s;
    return {(r[0] * t - q * r[1]) / det, (p * r[1] - s * r[0]) / det};
  }
  const double alpha = c[n - 1];  // A[n-1][0]
  const double beta = a[0];       // A[0][n-1]
  const double gamma = -b[0];
  Vec bb(b);
  bb[0] = b[0] - gamma;
  bb[n - 1] = b[n - 1] - alpha * beta / gamma;
  const Vec x = thomas(a, bb, c, r);
  Vec u(n, 0.0);
  u[0] = gamma;
  u[n - 1] = alpha;
  const Vec z = thomas(a, bb, c, u);
  const double fact = (x[0] + beta * x[n - 1] / gamma) / (1.0 + z[0] + beta * z[n - 1] / gamma);
  Vec out(n);
  for (size_t i = 0; i < n; ++i) out[i] = x[i] - fact * z[i];
  return out;
}

}  // namespace

CubicSpline::CubicSpline(const Vec& x, const Vec& y, const std::string& method) : method_(method) {
  if (method != "fmm" && method != "natural" && method != "periodic" && method != "monoH.FC")
    throw Error("spline: unknown method '" + method + "' (use fmm, natural, periodic or monoH.FC)");
  collapse_ties(x, y, x_, y_, "spline");
  const size_t n = x_.size();
  if (n < 2) throw Error("spline: need at least two distinct x values");
  if (method == "periodic") y_[n - 1] = y_[0];  // R: "first and last y values differ - using y[1] for both"
  Vec h(n - 1), del(n - 1);
  for (size_t i = 0; i + 1 < n; ++i) {
    h[i] = x_[i + 1] - x_[i];
    del[i] = (y_[i + 1] - y_[i]) / h[i];
  }
  b_.assign(n, 0.0);
  c_.assign(n, 0.0);
  d_.assign(n, 0.0);

  if (method == "monoH.FC") {
    // Fritsch & Carlson (1980): three-point slopes, then per-interval restriction to the monotonicity region.
    Vec m(n);
    m[0] = del[0];
    m[n - 1] = del[n - 2];
    for (size_t i = 1; i + 1 < n; ++i) m[i] = (del[i - 1] + del[i]) / 2.0;
    for (size_t k = 0; k + 1 < n; ++k) {
      if (del[k] == 0.0) {
        m[k] = m[k + 1] = 0.0;
        continue;
      }
      const double alpha = m[k] / del[k], beta = m[k + 1] / del[k];
      const double a2b2 = alpha * alpha + beta * beta;
      if (alpha >= 0.0 && beta >= 0.0 && a2b2 > 9.0) {  // only locally monotone pieces are restricted (as R)
        const double tau = 3.0 / std::sqrt(a2b2);
        m[k] = tau * alpha * del[k];
        m[k + 1] = tau * beta * del[k];
      }
    }
    // Cubic Hermite on each interval in power form.
    for (size_t i = 0; i + 1 < n; ++i) {
      b_[i] = m[i];
      c_[i] = (3.0 * del[i] - 2.0 * m[i] - m[i + 1]) / h[i];
      d_[i] = (m[i] + m[i + 1] - 2.0 * del[i]) / (h[i] * h[i]);
    }
    b_[n - 1] = m[n - 1];  // linear beyond the last knot
    return;
  }

  // Second derivatives M at the knots.
  Vec M(n, 0.0);
  if (n == 2) {
    // a straight line for every method (periodic: y[1] == y[0], so a constant)
  } else if (method == "periodic") {
    const size_t mm = n - 1;  // unknowns M_0..M_{n-2}; M_{n-1} = M_0
    Vec a(mm), b(mm), c(mm), r(mm);
    for (size_t i = 0; i < mm; ++i) {
      const size_t im = (i + mm - 1) % mm;  // interval to the left of knot i (wrapping)
      a[i] = h[im];
      c[i] = h[i];
      b[i] = 2.0 * (h[im] + h[i]);
      r[i] = 6.0 * (del[i] - del[im]);
    }
    const Vec s = cyclic(a, b, c, r);
    for (size_t i = 0; i < mm; ++i) M[i] = s[i];
    M[n - 1] = M[0];
  } else {
    // Interior equations for M_1..M_{n-2}.
    const size_t mm = n - 2;
    Vec a(mm), b(mm), c(mm), r(mm);
    for (size_t k = 0; k < mm; ++k) {
      const size_t i = k + 1;
      a[k] = h[i - 1];
      c[k] = h[i];
      b[k] = 2.0 * (h[i - 1] + h[i]);
      r[k] = 6.0 * (del[i] - del[i - 1]);
    }
    if (method == "fmm") {
      // s''' on the end intervals equals 6 * the third divided difference of the four end points (0 when n == 3):
      // M_0 = M_1 - 6 h_0 D_first, M_{n-1} = M_{n-2} + 6 h_{n-2} D_last; substitute into the first/last rows.
      double dfirst = 0.0, dlast = 0.0;
      if (n >= 4) {
        auto third_dd = [&](size_t i) {  // f[x_i, ..., x_{i+3}]
          const double d2a = (del[i + 1] - del[i]) / (x_[i + 2] - x_[i]);
          const double d2b = (del[i + 2] - del[i + 1]) / (x_[i + 3] - x_[i + 1]);
          return (d2b - d2a) / (x_[i + 3] - x_[i]);
        };
        dfirst = third_dd(0);
        dlast = third_dd(n - 4);
      }
      b[0] += h[0];
      r[0] += 6.0 * h[0] * h[0] * dfirst;
      b[mm - 1] += h[n - 2];
      r[mm - 1] -= 6.0 * h[n - 2] * h[n - 2] * dlast;
      const Vec s = thomas(a, b, c, r);
      for (size_t k = 0; k < mm; ++k) M[k + 1] = s[k];
      M[0] = M[1] - 6.0 * h[0] * dfirst;
      M[n - 1] = M[n - 2] + 6.0 * h[n - 2] * dlast;
    } else {  // natural: M_0 = M_{n-1} = 0
      const Vec s = thomas(a, b, c, r);
      for (size_t k = 0; k < mm; ++k) M[k + 1] = s[k];
    }
  }
  for (size_t i = 0; i + 1 < n; ++i) {
    b_[i] = del[i] - h[i] * (2.0 * M[i] + M[i + 1]) / 6.0;
    c_[i] = M[i] / 2.0;
    d_[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
  }
  const size_t L = n - 2;
  if (method == "periodic") {
    b_[n - 1] = b_[0];
    c_[n - 1] = c_[0];
    d_[n - 1] = d_[0];
  } else {
    b_[n - 1] = b_[L] + h[L] * (2.0 * c_[L] + 3.0 * d_[L] * h[L]);
    c_[n - 1] = M[n - 1] / 2.0;
    d_[n - 1] = method == "natural" ? 0.0 : d_[L];  // natural: linear to the right; fmm: continue the end cubic
  }
}

double CubicSpline::deriv(double x, int order) const {
  if (order < 0 || order > 3) throw Error("spline: deriv must be between 0 and 3");
  if (std::isnan(x)) return kNaN;
  const size_t n = x_.size();
  double u = x;
  if (method_ == "periodic" && (u < x_[0] || u > x_[n - 1])) {
    const double period = x_[n - 1] - x_[0];
    double r = std::fmod(u - x_[0], period);
    if (r < 0) r += period;
    u = x_[0] + r;
  }
  // knot i with x_i <= u < x_{i+1}; knot 0 to the left of the data, knot n-1 at and beyond the last x
  size_t i = static_cast<size_t>(std::upper_bound(x_.begin(), x_.end(), u) - x_.begin());
  i = i == 0 ? 0 : i - 1;
  double b = b_[i], c = c_[i], d = d_[i];
  if (method_ == "natural" && u <= x_[0]) d = 0.0;  // c_0 = 0: linear to the left
  if (method_ == "monoH.FC" && u < x_[0]) c = d = 0.0;
  const double t = u - x_[i];
  switch (order) {
    case 0: return y_[i] + t * (b + t * (c + t * d));
    case 1: return b + t * (2.0 * c + 3.0 * d * t);
    case 2: return 2.0 * c + 6.0 * d * t;
    default: return 6.0 * d;
  }
}

double CubicSpline::operator()(double x) const { return deriv(x, 0); }

Approx::Approx(const Vec& x, const Vec& y, const std::string& method, int rule_left, int rule_right, double f)
    : rule_left_(rule_left), rule_right_(rule_right), f_(f) {
  if (method == "linear") constant_ = false;
  else if (method == "constant") constant_ = true;
  else throw Error("approx: unknown method '" + method + "' (use linear or constant)");
  if ((rule_left != 1 && rule_left != 2) || (rule_right != 1 && rule_right != 2))
    throw Error("approx: rule must be 1 or 2");
  if (constant_ && !(f >= 0.0 && f <= 1.0)) throw Error("approx: f must be in [0, 1]");
  collapse_ties(x, y, x_, y_, "approx");
  if (x_.size() < (constant_ ? 1u : 2u))
    throw Error(constant_ ? "approx: need at least one point" : "approx: need at least two distinct x values");
}

double Approx::operator()(double v) const {
  if (std::isnan(v)) return kNaN;
  const size_t n = x_.size();
  if (v < x_[0]) return rule_left_ == 2 ? y_[0] : kNaN;
  if (v > x_[n - 1]) return rule_right_ == 2 ? y_[n - 1] : kNaN;
  // j: first knot with x_j >= v
  const size_t j = static_cast<size_t>(std::lower_bound(x_.begin(), x_.end(), v) - x_.begin());
  if (x_[j] == v) return y_[j];
  const size_t i = j - 1;
  if (constant_) {
    if (f_ == 0.0) return y_[i];
    if (f_ == 1.0) return y_[j];
    return (1.0 - f_) * y_[i] + f_ * y_[j];
  }
  return y_[i] + (y_[j] - y_[i]) * ((v - x_[i]) / (x_[j] - x_[i]));
}

}  // namespace rcli
