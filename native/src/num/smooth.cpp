// Scatterplot smoothers (clean-room): penalised cubic smoothing spline and loess.
// References:
//   Reinsch, C. H. (1967) Smoothing by spline functions. Numer. Math. 10, 177-183.
//   Green, P. J. and Silverman, B. W. (1994) Nonparametric Regression and Generalized Linear Models, Chapman & Hall,
//     sections 2.3 (Reinsch algorithm: Q, R matrices, g = y - lambda W^-1 Q gamma) and 3.2 (GCV, CV, leverages).
//   Hutchinson, M. F. and de Hoog, F. R. (1985) Smoothing noisy data with spline functions. Numer. Math. 47, 99-106
//     (band of the inverse of a banded matrix from its LDL' factor, giving the leverages in O(n)).
//   Craven, P. and Wahba, G. (1979) Smoothing noisy data with spline functions. Numer. Math. 31, 377-403 (GCV).
//   Kiefer, J. (1953) Sequential minimax search for a maximum. Proc. AMS 4, 502-506 (golden-section search).
//   Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. JASA 74, 829-836.
//   Cleveland, W. S., Grosse, E. and Shyu, W. M. (1992) Local regression models. Ch. 8 of Statistical Models in S
//     (span > 1, delta1/delta2 statistics, equivalent number of parameters).
//
// Smoothing spline conventions. The criterion is
//     sum_i w_i (ybar_i - g(x_i))^2 + lambda * integral g''(t)^2 dt
// over unique x_i (w_i = tie counts, ybar_i = tie means) with x rescaled to t = (x - min x) / (max x - min x) in
// [0, 1]. This is the scale R's smooth.spline uses for its 'lambda' when all distinct x are knots (nx <= 49), so an
// explicit lambda gives the same fit as R's smooth.spline(x, y, lambda = ) there. R's 'spar' is not supported.
#include "rcli/smooth.hpp"

#include <algorithm>
#include <numeric>

#include "rcli/summary.hpp"

namespace rcli {
namespace {

// ---------------------------------------------------------------------------------------------------------------
// Smoothing spline (Reinsch form)

struct SplineProblem {
  Vec t;     // knots on [0, 1]
  Vec h;     // t[i+1] - t[i]
  Vec w;     // weights at knots
  Vec ybar;  // weighted mean response at knots
  // Band of R (tridiagonal) and of Q' W^-1 Q (pentadiagonal), each of order m = n - 2.
  Vec r0, r1, p0, p1, p2;
};

struct SplineSolution {
  Vec g;      // fitted values at the knots
  Vec gamma;  // second derivatives (t scale) at the knots, zero at the ends
  Vec lev;
  double trace = 0.0;
};

// Q has columns j = 0..m-1 (interior knot j+1) with entries at rows j, j+1, j+2.
inline void q_column(const Vec& h, size_t j, double q[3]) {
  q[0] = 1.0 / h[j];
  q[1] = -1.0 / h[j] - 1.0 / h[j + 1];
  q[2] = 1.0 / h[j + 1];
}

void setup(SplineProblem& P) {
  const size_t n = P.t.size(), m = n - 2;
  P.h.resize(n - 1);
  for (size_t i = 0; i + 1 < n; ++i) P.h[i] = P.t[i + 1] - P.t[i];
  P.r0.assign(m, 0.0);
  P.r1.assign(m, 0.0);
  P.p0.assign(m, 0.0);
  P.p1.assign(m, 0.0);
  P.p2.assign(m, 0.0);
  for (size_t j = 0; j < m; ++j) {
    P.r0[j] = (P.h[j] + P.h[j + 1]) / 3.0;
    if (j + 1 < m) P.r1[j] = P.h[j + 1] / 6.0;
  }
  // (Q' W^-1 Q)_{jl} = sum_i Q_ij Q_il / w_i; column j touches rows j..j+2.
  for (size_t j = 0; j < m; ++j) {
    double a[3], b[3];
    q_column(P.h, j, a);
    for (int k = 0; k < 3; ++k) P.p0[j] += a[k] * a[k] / P.w[j + k];
    if (j + 1 < m) {
      q_column(P.h, j + 1, b);  // rows j+1..j+3
      P.p1[j] = a[1] * b[0] / P.w[j + 1] + a[2] * b[1] / P.w[j + 2];
    }
    if (j + 2 < m) {
      q_column(P.h, j + 2, b);  // rows j+2..j+4
      P.p2[j] = a[2] * b[0] / P.w[j + 2];
    }
  }
}

SplineSolution solve(const SplineProblem& P, double lambda, bool want_lev) {
  const size_t n = P.t.size(), m = n - 2;
  // B = R + lambda Q'W^-1 Q = L D L' (unit lower triangular L with two sub-diagonals l1, l2).
  Vec D(m), l1(m, 0.0), l2(m, 0.0);
  for (size_t i = 0; i < m; ++i) {
    const double b0 = P.r0[i] + lambda * P.p0[i];
    double d = b0;
    if (i >= 1) d -= l1[i - 1] * l1[i - 1] * D[i - 1];
    if (i >= 2) d -= l2[i - 2] * l2[i - 2] * D[i - 2];
    D[i] = d;
    if (i + 1 < m) {
      double b1 = P.r1[i] + lambda * P.p1[i];
      if (i >= 1) b1 -= l2[i - 1] * l1[i - 1] * D[i - 1];
      l1[i] = b1 / d;
    }
    if (i + 2 < m) l2[i] = lambda * P.p2[i] / d;
  }
  // Right-hand side Q' ybar.
  Vec z(m);
  for (size_t j = 0; j < m; ++j) {
    double q[3];
    q_column(P.h, j, q);
    z[j] = q[0] * P.ybar[j] + q[1] * P.ybar[j + 1] + q[2] * P.ybar[j + 2];
  }
  for (size_t i = 0; i < m; ++i) {
    if (i >= 1) z[i] -= l1[i - 1] * z[i - 1];
    if (i >= 2) z[i] -= l2[i - 2] * z[i - 2];
  }
  Vec gam(m);
  for (size_t i = m; i-- > 0;) {
    double v = z[i] / D[i];
    if (i + 1 < m) v -= l1[i] * gam[i + 1];
    if (i + 2 < m) v -= l2[i] * gam[i + 2];
    gam[i] = v;
  }
  SplineSolution S;
  S.g = P.ybar;
  for (size_t j = 0; j < m; ++j) {
    double q[3];
    q_column(P.h, j, q);
    for (int k = 0; k < 3; ++k) S.g[j + k] -= lambda * q[k] * gam[j] / P.w[j + k];
  }
  S.gamma.assign(n, 0.0);
  for (size_t j = 0; j < m; ++j) S.gamma[j + 1] = gam[j];
  if (!want_lev) return S;
  // Band (|i-j| <= 2) of Sigma = B^-1 (Hutchinson & de Hoog): Sigma = D^-1 L^-1 + (I - L') Sigma, from the bottom.
  Vec s0(m), s1(m, 0.0), s2(m, 0.0);  // Sigma_ii, Sigma_{i,i+1}, Sigma_{i,i+2}
  for (size_t i = m; i-- > 0;) {
    const double a1 = i + 1 < m ? l1[i] : 0.0, a2 = i + 2 < m ? l2[i] : 0.0;
    if (i + 2 < m) s2[i] = -a1 * s1[i + 1] - a2 * s0[i + 2];
    if (i + 1 < m) s1[i] = -a1 * s0[i + 1] - (i + 2 < m ? a2 * s1[i + 1] : 0.0);
    s0[i] = 1.0 / D[i] - a1 * s1[i] - a2 * s2[i];
  }
  auto sigma = [&](size_t a, size_t b) {
    if (a > b) std::swap(a, b);
    const size_t d = b - a;
    return d == 0 ? s0[a] : d == 1 ? s1[a] : d == 2 ? s2[a] : 0.0;
  };
  S.lev.assign(n, 0.0);
  S.trace = 0.0;
  for (size_t i = 0; i < n; ++i) {
    // Q row i has entries in columns j = i-2, i-1, i (those in 0..m-1): Q_{i,j} = q_column(j)[i - j].
    size_t cols[3];
    double vals[3];
    int k = 0;
    for (size_t j = (i >= 2 ? i - 2 : 0); j <= i && j < m; ++j) {
      double q[3];
      q_column(P.h, j, q);
      cols[k] = j;
      vals[k] = q[i - j];
      ++k;
    }
    double qbq = 0.0;
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) qbq += vals[a] * vals[b] * sigma(cols[a], cols[b]);
    S.lev[i] = 1.0 - lambda * qbq / P.w[i];
    S.trace += S.lev[i];
  }
  return S;
}

double trace_at(const SplineProblem& P, double rho) { return solve(P, std::exp(rho), true).trace; }

}  // namespace

SmoothSplineFit smooth_spline(const Vec& x, const Vec& y, const SmoothSplineOptions& opt) {
  if (x.size() != y.size()) throw Error("smooth.spline: x and y lengths differ");
  for (size_t i = 0; i < x.size(); ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
      throw Error("smooth.spline: missing or infinite values are not allowed");
  const size_t nobs = x.size();
  if (nobs < 4) throw Error("smooth.spline: need at least four unique x values");
  // Group x values closer than tol = 1e-6 * IQR(x) (R's default), each group represented by its smallest x.
  const double tol = 1e-6 * (summary::quantile(x, 0.75) - summary::quantile(x, 0.25));
  std::vector<size_t> idx(nobs);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  const double xm = summary::mean(x);
  auto key = [&](size_t i) { return tol > 0 ? std::round((x[i] - xm) / tol) : x[i]; };
  SmoothSplineFit fit;
  Vec ysum;
  std::vector<std::vector<size_t>> groups;
  for (size_t k = 0; k < nobs;) {
    size_t l = k;
    const double kk = key(idx[k]);
    groups.emplace_back();
    while (l < nobs && key(idx[l]) == kk) groups.back().push_back(idx[l++]);
    fit.x.push_back(x[idx[k]]);
    k = l;
  }
  const size_t nx = fit.x.size();
  if (nx < 4) throw Error("smooth.spline: need at least four unique x values");
  SplineProblem P;
  const double xmin = fit.x.front(), range = fit.x.back() - fit.x.front();
  P.t.resize(nx);
  P.w.resize(nx);
  P.ybar.resize(nx);
  double within_ss = 0.0;  // residual sum of squares inside tie groups (part of the GCV numerator)
  for (size_t i = 0; i < nx; ++i) {
    P.t[i] = (fit.x[i] - xmin) / range;
    Vec yy;
    for (size_t j : groups[i]) yy.push_back(y[j]);
    P.w[i] = static_cast<double>(yy.size());
    P.ybar[i] = summary::mean(yy);
    for (double v : yy) within_ss += (v - P.ybar[i]) * (v - P.ybar[i]);
  }
  setup(P);

  // Search interval for rho = log(lambda): from (nearly) interpolation to (nearly) a straight line.
  // With n points on [0,1] the equivalent kernel half-width is about (lambda / n)^(1/4).
  const double sumw = static_cast<double>(nobs);
  double hmin = kInf;
  for (double hh : P.h) hmin = std::min(hmin, hh);
  double rho_lo = std::log(sumw * std::pow(hmin, 3.0)) - 4.0;
  double rho_hi = std::log(sumw) + 12.0;
  for (int it = 0; it < 40 && trace_at(P, rho_lo) < static_cast<double>(nx) - 0.05; ++it) rho_lo -= 3.0;
  for (int it = 0; it < 40 && trace_at(P, rho_hi) > 2.05; ++it) rho_hi += 3.0;

  auto crit = [&](const SplineSolution& S) {
    double c = 0.0;
    if (opt.cv) {
      for (size_t i = 0; i < nx; ++i) {
        const double r = (P.ybar[i] - S.g[i]) / (1.0 - S.lev[i]);
        c += P.w[i] * r * r;
      }
      return c / sumw;
    }
    double rss = within_ss;
    for (size_t i = 0; i < nx; ++i) rss += P.w[i] * (P.ybar[i] - S.g[i]) * (P.ybar[i] - S.g[i]);
    const double den = 1.0 - S.trace / sumw;
    return (rss / sumw) / (den * den);
  };

  double lambda;
  if (!std::isnan(opt.lambda)) {
    if (!(opt.lambda >= 0.0) || !std::isfinite(opt.lambda)) throw Error("smooth.spline: lambda must be >= 0 and finite");
    lambda = opt.lambda;
  } else if (!std::isnan(opt.df)) {
    if (!(opt.df > 1.0 && opt.df <= static_cast<double>(nx)))
      throw Error("smooth.spline: df must be in (1, nx], nx = number of unique x");
    // trace(S_lambda) decreases monotonically in rho; bisection to machine precision in rho.
    double lo = rho_lo, hi = rho_hi;
    if (trace_at(P, lo) <= opt.df) hi = lo;
    else if (trace_at(P, hi) >= opt.df) lo = hi;
    for (int it = 0; it < 200 && hi - lo > 1e-13 * std::max(1.0, std::fabs(lo)); ++it) {
      const double mid = 0.5 * (lo + hi);
      if (trace_at(P, mid) > opt.df) lo = mid;
      else hi = mid;
    }
    lambda = std::exp(0.5 * (lo + hi));
  } else {
    // GCV / CV: grid over rho, then golden-section search on the bracket around the best grid point.
    auto score = [&](double rho) { return crit(solve(P, std::exp(rho), true)); };
    const int ngrid = 80;
    double best = kInf;
    int kbest = 0;
    Vec grid(ngrid + 1);
    for (int k = 0; k <= ngrid; ++k) {
      grid[k] = rho_lo + (rho_hi - rho_lo) * k / ngrid;
      const double s = score(grid[k]);
      if (s < best) { best = s; kbest = k; }
    }
    double a = grid[std::max(kbest - 1, 0)], b = grid[std::min(kbest + 1, ngrid)];
    const double g = (std::sqrt(5.0) - 1.0) / 2.0;
    double c = b - g * (b - a), d = a + g * (b - a);
    double fc = score(c), fd = score(d);
    while (b - a > 1e-9) {
      if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = score(c); }
      else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = score(d); }
    }
    double rho = 0.5 * (a + b);
    if (best < score(rho)) rho = grid[kbest];
    lambda = std::exp(rho);
  }

  const SplineSolution S = solve(P, lambda, true);
  fit.y = S.g;
  fit.w = P.w;
  fit.lev = S.lev;
  fit.df = S.trace;
  fit.lambda = lambda;
  fit.cv_crit = crit(S);
  fit.pen_crit = 0.0;
  for (size_t i = 0; i < nx; ++i) fit.pen_crit += P.w[i] * (P.ybar[i] - S.g[i]) * (P.ybar[i] - S.g[i]);
  // Piecewise cubic in x units: second derivatives scale by 1/range^2.
  fit.coef_b.assign(nx, 0.0);
  fit.coef_c.assign(nx, 0.0);
  fit.coef_d.assign(nx, 0.0);
  const double s2 = 1.0 / (range * range);
  for (size_t i = 0; i + 1 < nx; ++i) {
    const double hx = fit.x[i + 1] - fit.x[i];
    const double M0 = S.gamma[i] * s2, M1 = S.gamma[i + 1] * s2;
    fit.coef_b[i] = (fit.y[i + 1] - fit.y[i]) / hx - hx * (2.0 * M0 + M1) / 6.0;
    fit.coef_c[i] = M0 / 2.0;
    fit.coef_d[i] = (M1 - M0) / (6.0 * hx);
  }
  const size_t L = nx - 2;
  const double hl = fit.x[nx - 1] - fit.x[L];
  fit.coef_b[nx - 1] = fit.coef_b[L] + hl * (2.0 * fit.coef_c[L] + 3.0 * fit.coef_d[L] * hl);
  return fit;
}

double SmoothSplineFit::predict(double x0, int deriv) const {
  if (deriv < 0 || deriv > 3) throw Error("predict.smooth.spline: deriv must be between 0 and 3");
  if (std::isnan(x0)) return kNaN;
  const size_t n = x.size();
  if (x0 < x.front() || x0 > x.back()) {  // linear extrapolation
    const size_t e = x0 < x.front() ? 0 : n - 1;
    if (deriv == 0) return y[e] + coef_b[e] * (x0 - x[e]);
    return deriv == 1 ? coef_b[e] : 0.0;
  }
  size_t i = static_cast<size_t>(std::upper_bound(x.begin(), x.end(), x0) - x.begin());
  i = i == 0 ? 0 : std::min(i - 1, n - 2);
  const double t = x0 - x[i], b = coef_b[i], c = coef_c[i], d = coef_d[i];
  switch (deriv) {
    case 0: return y[i] + t * (b + t * (c + t * d));
    case 1: return b + t * (2.0 * c + 3.0 * d * t);
    case 2: return 2.0 * c + 6.0 * d * t;
    default: return 6.0 * d;
  }
}

// ---------------------------------------------------------------------------------------------------------------
// loess

namespace {

// Hat-matrix row of the local fit at x0: fitted value = sum_j row[j] * y[lo + j], support xs[lo..hi).
struct LocalRow {
  size_t lo = 0, hi = 0;
  Vec row;
};

// Local weighted polynomial regression at x0 on sorted xs with prior weights pw (robustness weights).
// Returns false if no point has positive weight.
bool local_row(const Vec& xs, const Vec& pw, double x0, double span, int degree, LocalRow& out) {
  const size_t n = xs.size();
  Vec dist(n);
  for (size_t i = 0; i < n; ++i) dist[i] = std::fabs(xs[i] - x0);
  size_t q = static_cast<size_t>(std::floor(static_cast<double>(n) * span));
  q = std::clamp<size_t>(q, 1, n);
  double h;
  if (span <= 1.0) {
    Vec tmp(dist);
    std::nth_element(tmp.begin(), tmp.begin() + (q - 1), tmp.end());
    h = tmp[q - 1];
  } else {
    // All points; the maximum distance is enlarged. ?loess documents alpha^(1/p) (p = 1 here), but R 4.3.3 uses
    // sqrt(alpha) for one predictor (black-box: loess(y ~ x, span = 1.5) on 1:10); we follow R.
    h = *std::max_element(dist.begin(), dist.end()) * std::sqrt(span);
  }
  // support: points strictly within h (tricube weight positive); all points at distance 0 if h == 0
  std::vector<size_t> sup;
  Vec wts;
  for (size_t i = 0; i < n; ++i) {
    double wi;
    if (h > 0) {
      const double u = dist[i] / h;
      if (u >= 1.0) continue;
      const double v = 1.0 - u * u * u;
      wi = v * v * v;
    } else {
      if (dist[i] > 0) continue;
      wi = 1.0;
    }
    wi *= pw[i];
    if (wi <= 0) continue;
    sup.push_back(i);
    wts.push_back(wi);
  }
  if (sup.empty()) return false;
  const size_t k = sup.size();
  // Design columns 1, u, u^2 with u = (x - x0)/scale, rows scaled by sqrt(w); Householder QR.
  const double scale = h > 0 ? h : 1.0;
  int p = degree + 1;
  std::vector<Vec> A(p, Vec(k));
  for (size_t r = 0; r < k; ++r) {
    const double sw = std::sqrt(wts[r]), u = (xs[sup[r]] - x0) / scale;
    double pw_ = sw;
    for (int c = 0; c < p; ++c) { A[c][r] = pw_; pw_ *= u; }
  }
  // Modified Gram-Schmidt (twice) with rank-deficiency dropping: columns whose residual norm is negligible
  // relative to their original norm are removed (the local fit then uses a lower degree).
  std::vector<Vec> Qc;
  std::vector<Vec> Rcol;  // R entries per kept column
  std::vector<int> kept;
  for (int c = 0; c < p; ++c) {
    Vec v = A[c];
    double norm0 = 0;
    for (double t : v) norm0 += t * t;
    norm0 = std::sqrt(norm0);
    Vec rc(Qc.size() + 1, 0.0);
    for (int pass = 0; pass < 2; ++pass)
      for (size_t j = 0; j < Qc.size(); ++j) {
        double dot = 0;
        for (size_t r = 0; r < k; ++r) dot += Qc[j][r] * v[r];
        for (size_t r = 0; r < k; ++r) v[r] -= dot * Qc[j][r];
        rc[j] += dot;
      }
    double nv = 0;
    for (double t : v) nv += t * t;
    nv = std::sqrt(nv);
    if (!(nv > 1e-7 * norm0) || norm0 == 0) continue;
    for (double& t : v) t /= nv;
    rc.back() = nv;
    Qc.push_back(v);
    Rcol.push_back(rc);
    kept.push_back(c);
  }
  // Fitted value at x0 is the intercept coefficient: beta = R^-1 Q' sqrt(W) y; the intercept column is kept first
  // (A[0] = sqrt(w) > 0), so beta_0 = e1' R^-1 Q' sqrt(W) y. Solve R' v = e1, row = sqrt(w) * (Q v).
  const size_t r = Qc.size();
  Vec v(r, 0.0);
  for (size_t i = 0; i < r; ++i) {  // forward substitution with R' (lower triangular): R'(i, j) = Rcol[i][j]
    double s = i == 0 ? 1.0 : 0.0;
    for (size_t j = 0; j < i; ++j) s -= Rcol[i][j] * v[j];
    v[i] = s / Rcol[i][i];
  }
  out.lo = sup.front();
  out.hi = sup.back() + 1;
  out.row.assign(out.hi - out.lo, 0.0);
  for (size_t rr = 0; rr < k; ++rr) {
    double s = 0;
    for (size_t i = 0; i < r; ++i) s += Qc[i][rr] * v[i];
    out.row[sup[rr] - out.lo] = std::sqrt(wts[rr]) * s;
  }
  return true;
}

double bisquare(double u) {
  const double a = std::fabs(u);
  if (a <= 1e-3) return 1.0;  // R's robustness weights are exactly 1 / 0 within 0.1% of 0 / 1 (black-box)
  if (a >= 0.999) return 0.0;
  const double v = 1.0 - u * u;
  return v * v;
}

}  // namespace

LoessFit loess(const Vec& x, const Vec& y, const LoessOptions& opt) {
  if (x.size() != y.size()) throw Error("loess: x and y lengths differ");
  for (size_t i = 0; i < x.size(); ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) throw Error("loess: missing or infinite values are not allowed");
  if (opt.degree < 0 || opt.degree > 2) throw Error("loess: degree must be 0, 1 or 2");
  if (opt.family != "gaussian" && opt.family != "symmetric")
    throw Error("loess: family must be 'gaussian' or 'symmetric'");
  if (!(opt.span > 0) || !std::isfinite(opt.span)) throw Error("loess: span must be positive");
  const size_t n = x.size();
  if (n < static_cast<size_t>(opt.degree) + 1 || n < 2) throw Error("loess: too few data points");
  if (static_cast<double>(n) * opt.span < static_cast<double>(opt.degree) + 1.0)
    throw Error("loess: span too small (fewer than degree + 1 points in each neighbourhood)");

  LoessFit fit;
  fit.x = x;
  fit.y = y;
  fit.opt = opt;
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  Vec xs(n), ys(n);
  for (size_t i = 0; i < n; ++i) { xs[i] = x[ord[i]]; ys[i] = y[ord[i]]; }

  auto fit_all = [&](const Vec& pw, std::vector<LocalRow>& rows, Vec& fv) {
    rows.assign(n, LocalRow());
    fv.assign(n, kNaN);
    for (size_t i = 0; i < n; ++i) {
      if (!local_row(xs, pw, xs[i], opt.span, opt.degree, rows[i])) continue;
      double s = 0;
      for (size_t j = rows[i].lo; j < rows[i].hi; ++j) s += rows[i].row[j - rows[i].lo] * ys[j];
      fv[i] = s;
    }
  };

  Vec rw(n, 1.0);
  std::vector<LocalRow> rows0;
  Vec fv;
  fit_all(rw, rows0, fv);
  // Statistics from the (non-robust) operator L: trace, delta1 = tr (I-L)'(I-L), delta2 = tr ((I-L)'(I-L))^2.
  {
    double tr = 0, ll = 0;
    for (size_t i = 0; i < n; ++i) {
      const LocalRow& R = rows0[i];
      for (size_t j = R.lo; j < R.hi; ++j) {
        const double v = R.row[j - R.lo];
        ll += v * v;
        if (j == i) tr += v;
      }
    }
    fit.trace_hat = tr;
    fit.enp = ll;  // R's enp is tr(L'L)
    // A = I - L as dense rows over their support [lo, hi) (the diagonal is always inside the support).
    // M = A A' (same nonzero spectrum as A'A): delta1 = tr M, delta2 = ||M||_F^2.
    std::vector<size_t> lo(n), hi(n);
    std::vector<Vec> A(n);
    for (size_t i = 0; i < n; ++i) {
      const LocalRow& R = rows0[i];
      lo[i] = std::min(R.lo, i);
      hi[i] = std::max(R.hi, i + 1);
      A[i].assign(hi[i] - lo[i], 0.0);
      for (size_t j = R.lo; j < R.hi; ++j) A[i][j - lo[i]] = -R.row[j - R.lo];
      A[i][i - lo[i]] += 1.0;
    }
    // delta2 = sum_{i,k} M_ik^2 over the upper triangle (M symmetric): O(n^2 q) work, so exact statistics cost
    // O(n^3) for a fixed span, as with R's statistics = "exact". Dot products of contiguous windows use four
    // partial sums.
    double d1 = 0, d2 = 0;
    for (size_t i = 0; i < n; ++i) {
      double row2 = 0;  // sum over k > i of M_ik^2
      for (size_t k = i; k < n; ++k) {
        const size_t a = std::max(lo[i], lo[k]), b = std::min(hi[i], hi[k]);
        if (a >= b) continue;
        const double* p = A[i].data() + (a - lo[i]);
        const double* q = A[k].data() + (a - lo[k]);
        const size_t len = b - a;
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        size_t j = 0;
        for (; j + 4 <= len; j += 4) {
          s0 += p[j] * q[j];
          s1 += p[j + 1] * q[j + 1];
          s2 += p[j + 2] * q[j + 2];
          s3 += p[j + 3] * q[j + 3];
        }
        for (; j < len; ++j) s0 += p[j] * q[j];
        const double sum = (s0 + s1) + (s2 + s3);
        if (k == i) { d1 += sum; d2 += sum * sum; }
        else row2 += sum * sum;
      }
      d2 += 2 * row2;
    }
    fit.one_delta = d1;
    fit.two_delta = d2;
  }
  if (opt.family == "symmetric") {
    // R's iterations = 4 means four fits in all: the least-squares fit and three bisquare reweightings
    // (black-box: matches R to 1e-10 only with three reweightings). The robustness weights reported are those used in
    // the last fit.
    std::vector<LocalRow> rows;
    for (int iter = 0; iter < 3; ++iter) {
      Vec ares(n);
      for (size_t i = 0; i < n; ++i) ares[i] = std::fabs(ys[i] - fv[i]);
      const double s6 = 6.0 * summary::median(ares);
      for (size_t i = 0; i < n; ++i) rw[i] = s6 > 0 ? bisquare((ys[i] - fv[i]) / s6) : (ares[i] == 0 ? 1.0 : 0.0);
      fit_all(rw, rows, fv);
    }
  }
  fit.fitted.assign(n, kNaN);
  fit.residuals.assign(n, kNaN);
  fit.robustness_weights.assign(n, 1.0);
  double rss = 0;
  for (size_t i = 0; i < n; ++i) {
    fit.fitted[ord[i]] = fv[i];
    fit.residuals[ord[i]] = ys[i] - fv[i];
    fit.robustness_weights[ord[i]] = rw[i];
    rss += (ys[i] - fv[i]) * (ys[i] - fv[i]);
  }
  fit.s = std::sqrt(rss / fit.one_delta);
  if (opt.family == "symmetric") {
    // Robust scale from pseudovalues (Cleveland, Grosse & Shyu 1992, sec. 8.4): with u = r / (6 median|r|),
    // ytilde = yhat + 6 median|r| * psi(u) / mean(psi'(u)) for the bisquare psi(u) = u (1 - u^2)^2; s is the residual
    // scale of the least-squares smooth of ytilde. This approximates R's s for family = "symmetric" (see tests).
    Vec ares(n);
    for (size_t i = 0; i < n; ++i) ares[i] = std::fabs(ys[i] - fv[i]);
    const double S = 6.0 * summary::median(ares);
    double mdpsi = 0;
    Vec u(n);
    for (size_t i = 0; i < n; ++i) {
      u[i] = S > 0 ? (ys[i] - fv[i]) / S : 0.0;
      if (std::fabs(u[i]) < 1) mdpsi += (1 - u[i] * u[i]) * (1 - 5 * u[i] * u[i]);
    }
    mdpsi /= static_cast<double>(n);
    if (S > 0 && mdpsi > 0) {
      Vec yp(n);
      for (size_t i = 0; i < n; ++i) {
        const double v = std::fabs(u[i]) < 1 ? 1 - u[i] * u[i] : 0.0;
        yp[i] = fv[i] + S * u[i] * v * v / mdpsi;
      }
      double ss = 0;
      for (size_t i = 0; i < n; ++i) {
        const LocalRow& R = rows0[i];
        double sm = 0;
        for (size_t j = R.lo; j < R.hi; ++j) sm += R.row[j - R.lo] * yp[j];
        ss += (yp[i] - sm) * (yp[i] - sm);
      }
      fit.s = std::sqrt(ss / fit.one_delta);
    }
  }
  return fit;
}

double LoessFit::predict(double x0) const {
  if (std::isnan(x0) || x.empty()) return kNaN;
  const auto mm = std::minmax_element(x.begin(), x.end());
  if (x0 < *mm.first || x0 > *mm.second) return kNaN;
  const size_t n = x.size();
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  Vec xs(n), ys(n), pw(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    xs[i] = x[ord[i]];
    ys[i] = y[ord[i]];
    if (!robustness_weights.empty()) pw[i] = robustness_weights[ord[i]];
  }
  LocalRow R;
  if (!local_row(xs, pw, x0, opt.span, opt.degree, R)) return kNaN;
  double s = 0;
  for (size_t j = R.lo; j < R.hi; ++j) s += R.row[j - R.lo] * ys[j];
  return s;
}

}  // namespace rcli
