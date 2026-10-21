// Noncentral t, chi-square and F distributions and the studentized range (clean-room; written from published
// mathematics only).
//
// Sources, per algorithm:
// - Noncentral chi-square: the Poisson mixture of central chi-squares, F(x) = sum_j Pois(j; ncp/2) P(chisq_{df+2j} <= x)
//   (Johnson, Kotz & Balakrishnan 1995, Continuous Univariate Distributions vol. 2, eq. 29.2; the same series as
//   Ding 1992, AS 275, Applied Statistics 41, 478-482). As in Ding's algorithm the summation starts at the largest
//   term and runs outward in both directions; each tail is summed on its own (lower: P(a, x), upper: Q(a, x)), so both
//   tails keep full relative accuracy instead of one being 1 minus the other. Density: the same mixture of chi-square
//   densities (JKB eq. 29.4).
// - Noncentral F: Poisson mixture of incomplete beta functions, F(x) = sum_j Pois(j; ncp/2) I_y(df1/2 + j, df2/2) with
//   y = df1 x / (df1 x + df2) (JKB eq. 30.6; the noncentral beta series of Lenth 1987, AS 226, and Frick 1990, AS R84);
//   density as the mixture of beta densities. df2 = Inf reduces to chisq(df1, ncp) / df1.
// - Noncentral t: R. V. Lenth (1989), Algorithm AS 243, Applied Statistics 38, 185-189. For t >= 0, with
//   x = t^2 / (t^2 + df), P(T <= t) = Phi(-d) + 1/2 sum_j [p_j I_x(j + 1/2, df/2) + q_j I_x(j + 1, df/2)],
//   p_j = e^{-d^2/2} (d^2/2)^j / j!, q_j = d e^{-d^2/2} (d^2/2)^j / (sqrt(2) Gamma(j + 3/2)). Since the same sum at
//   x = 1 equals 1 - Phi(-d), the upper tail is 1/2 sum_j [p_j (1 - I_x(j + 1/2, .)) + q_j (1 - I_x(j + 1, .))],
//   which for d > 0 is a sum of positive terms. Negative t uses P(T <= t; d) = P(T >= -t; -d). When d < 0 makes the
//   upper-tail series cancel, the tail is computed from the mixture integral
//   P(T > t) = int_0^Inf Phi(-(t sqrt(s/df) - d)) f_{chisq_df}(s) ds (JKB eq. 31.4), whose integrand is positive.
//   Density: f(t) = int_0^Inf phi(t v - d) v g(v) dv with v = sqrt(s/df) (the definition T = (Z + d) / V).
// - Studentized range: P(Q <= q) = int_0^Inf f_S(s) W(q s)^nranges ds, S = sqrt(chisq_df / df), with the range
//   distribution of k standard normals W(w) = k int phi(z) [Phi(z) - Phi(z - w)]^{k-1} dz (Hartley 1942;
//   Copenhaver & Holland 1988, J. Statist. Comput. Simul. 30, 1-15, who evaluate the same double integral by
//   Gauss-Legendre quadrature). Here both integrals are evaluated by adaptive Gauss-Kronrod quadrature (QUADPACK
//   QAGS, quad.cpp) in log-scaled form around the located mode of the integrand.
// - Quantiles: inversion of the smaller tail on the log scale by the Illinois variant of regula falsi (Dowell &
//   Jarratt 1971) after a doubling search for a bracket, in log x for the positive distributions.
//
// Integration variables: for S^2 = chisq_df / df we integrate over w = log(S^2), where the log density
// -(df/2)(e^w - 1 - w) + C(df) is concave and stays accurate for huge df (no lgamma(df/2) cancellation).
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "rcli/dist.hpp"
#include "rcli/optim.hpp"
#include "rcli/special.hpp"

namespace rcli {

namespace {

namespace sp = rcli::special;

constexpr double kLn2 = 0.693147180559945309417232121458;
constexpr double kLnSqrt2Pi = 0.918938533204672741780329736406;

inline double log1mexp(double x) { return x > -kLn2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x)); }
inline double d0(bool lg) { return lg ? -kInf : 0.0; }
inline double from_log(double l, bool lg) { return lg ? l : std::exp(l); }
inline double pstep(bool below, bool lower, bool lg) { return (below == lower) ? d0(lg) : (lg ? 0.0 : 1.0); }
inline double logadd(double a, double b) {
  if (a == -kInf) return b;
  if (b == -kInf) return a;
  const double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::fabs(a - b)));
}
inline double lphi_lower(double z) {
  double l, u;
  sp::pnorm_both(z, &l, &u, true);
  return l;
}
inline double lphi_upper(double z) {
  double l, u;
  sp::pnorm_both(z, &l, &u, true);
  return u;
}
inline double ldnorm(double z) { return -0.5 * z * z - kLnSqrt2Pi; }

// e^w - 1 - w without cancellation for small |w|.
double expm1mx(double w) {
  if (std::fabs(w) > 0.5) return std::expm1(w) - w;
  double term = 0.5 * w * w, sum = term;
  for (int k = 3; k < 40; ++k) {
    term *= w / k;
    sum += term;
    if (std::fabs(term) <= 1e-17 * sum) break;
  }
  return sum;
}

// For W = log(chisq_df / df): log density -(df/2) (e^w - 1 - w) + C(df), with
// C(df) = (a) log a - a - lgamma(a) = 1/2 log a - log sqrt(2 pi) - stirlerr(a), a = df/2.
double lchisq_w_const(double df) {
  const double a = 0.5 * df;
  return 0.5 * std::log(a) - kLnSqrt2Pi - sp::stirlerr(a);
}
inline double lchisq_w(double w, double df, double c) { return -0.5 * df * expm1mx(w) + c; }

// --------------------------------------------------------------------------------------------------------------
// Sum over integers j >= 0 of exp(lt(j)) for a log-concave (hence unimodal) sequence; returns the log of the sum.
// The largest term is located by a doubling search from j0 followed by an integer ternary search; the sum then runs
// outward and stops in each direction when the geometric bound t r / (1 - r) on the remainder (valid because the term
// ratios r of a log-concave sequence decrease away from the mode) is below 1e-17 of the sum.
template <class F>
double log_series(F&& lt, double j0) {
  double j = std::max(0.0, std::floor(j0));
  double lj = lt(j);
  for (int g = 0; lj == -kInf && g < 8; ++g) lj = lt(++j);
  double lo = j, hi = j;
  const double lup = lt(j + 1);
  if (lup > lj) {
    double a = j, b = j + 1, lb = lup, step = 1;
    for (;;) {
      step *= 2;
      const double c = b + step, lc = lt(c);
      if (!(lc > lb)) { lo = a; hi = c; break; }
      a = b; b = c; lb = lc;
      if (step > 1e15) { lo = a; hi = c; break; }
    }
  } else if (j > 0) {
    const double ldn = lt(j - 1);
    if (ldn > lj) {
      double a = j, b = j - 1, lb = ldn, step = 1;
      for (;;) {
        step *= 2;
        const double c = std::max(0.0, b - step), lc = lt(c);
        if (!(lc > lb) || c == 0) { lo = c; hi = a; break; }
        a = b; b = c; lb = lc;
      }
    }
  }
  while (hi - lo > 2) {
    const double third = std::floor((hi - lo) / 3);
    const double m1 = lo + third, m2 = hi - third;
    if (lt(m1) < lt(m2)) lo = m1;
    else hi = m2;
  }
  double peak = lo, lmax = lt(lo);
  for (double i = lo + 1; i <= hi; i += 1) {
    const double li = lt(i);
    if (li > lmax) { lmax = li; peak = i; }
  }
  if (lmax == -kInf || std::isnan(lmax)) return lmax;
  double s = 1, prev = 1;
  for (double i = peak + 1; i - peak < 1e7; i += 1) {
    const double t = std::exp(lt(i) - lmax);
    s += t;
    if (!(t > 0)) break;
    const double r = t / prev;
    prev = t;
    if (r < 1 && t * r / (1 - r) <= 1e-17 * s) break;
  }
  prev = 1;
  for (double i = peak - 1; i >= 0; i -= 1) {
    const double t = std::exp(lt(i) - lmax);
    s += t;
    if (!(t > 0)) break;
    const double r = t / prev;
    prev = t;
    if (r < 1 && t * r / (1 - r) <= 1e-17 * s) break;
  }
  return lmax + std::log(s);
}

// --------------------------------------------------------------------------------------------------------------
// log of int exp(h(u)) du over the real line for a unimodal h (log-concave in every use here). The mode is bracketed
// from u0 with steps of `scale` and refined by golden-section search (to 1e-4 scale: the mode only splits the range
// and scales the integrand); the range is cut where h has dropped by `drop` and integrated by QAGS to relative
// tolerance `rel_tol` on both sides of the mode.
template <class H>
double log_integral(H&& h, double u0, double scale, double rel_tol = 1e-14, double drop = 50) {
  auto hv = [&](double u) {
    const double v = h(u);
    return std::isnan(v) ? -kInf : v;
  };
  double a = u0 - scale, b = u0, c = u0 + scale;
  double fa = hv(a), fb = hv(b), fc = hv(c);
  for (int it = 0; it < 200 && !(fb >= fa && fb >= fc); ++it) {
    if (fc > fb || (fa == -kInf && fb == -kInf && fc == -kInf && it % 2 == 0)) {
      const double step = 2 * (c - b);
      a = b; fa = fb; b = c; fb = fc; c = b + step; fc = hv(c);
    } else {
      const double step = 2 * (b - a);
      c = b; fc = fb; b = a; fb = fa; a = b - step; fa = hv(a);
    }
  }
  constexpr double kG = 0.381966011250105151795413165634;  // 2 - golden ratio
  for (int it = 0; it < 200 && (c - a) > 1e-4 * scale; ++it) {
    double x;
    if (b - a > c - b) x = b - kG * (b - a);
    else x = b + kG * (c - b);
    const double fx = hv(x);
    if (fx > fb) {
      if (x < b) { c = b; fc = fb; } else { a = b; fa = fb; }
      b = x; fb = fx;
    } else {
      if (x < b) { a = x; fa = fx; } else { c = x; fc = fx; }
    }
  }
  const double mode = b, hmax = fb;
  if (hmax == -kInf || !std::isfinite(hmax)) return hmax == kInf ? kInf : -kInf;
  auto edge = [&](double dir) {
    double d = scale;
    for (int it = 0; it < 80; ++it) {
      if (hv(mode + dir * d) < hmax - drop) return mode + dir * d;
      d *= 2;
    }
    return mode + dir * d;
  };
  const double lo = edge(-1), hi = edge(1);
  IntegrateOptions opt;
  opt.subdivisions = 500;
  opt.rel_tol = rel_tol;
  opt.abs_tol = 1e-300;
  const Fn1 f = [&](double u) { return std::exp(hv(u) - hmax); };
  const double v = integrate(f, lo, mode, opt).value + integrate(f, mode, hi, opt).value;
  return hmax + std::log(v);
}

// --------------------------------------------------------------------------------------------------------------
// Noncentral chi-square

// log P(X <= x) (lower) or log P(X > x), 0 < x < Inf.
double ncchisq_ltail(double x, double k, double lam, bool lower) {
  const double mu = 0.5 * lam, y = 0.5 * x;
  auto lt = [&](double j) {
    const double a = 0.5 * k + j;
    const double g = a == 0 ? (lower ? 0.0 : -kInf) : sp::gamma_inc(a, y, lower, true);
    return sp::dpois_raw(j, mu, true) + g;
  };
  return std::min(0.0, log_series(lt, mu));
}

double ncchisq_p(double x, double k, double lam, bool lower, bool lg) {
  if (x < 0) return pstep(true, lower, lg);
  if (x == 0) {
    if (k > 0) return pstep(true, lower, lg);
    const double l0 = -0.5 * lam;  // point mass at 0 when df = 0
    return from_log(lower ? l0 : log1mexp(l0), lg);
  }
  if (x == kInf) return pstep(false, lower, lg);
  const double l = ncchisq_ltail(x, k, lam, lower);
  if (lg && l > -kLn2) return log1mexp(ncchisq_ltail(x, k, lam, !lower));
  return from_log(l, lg);
}

double ncchisq_d(double x, double k, double lam, bool lg) {
  if (x < 0 || x == kInf) return d0(lg);
  const double mu = 0.5 * lam;
  if (x == 0) {
    if (k < 2) return kInf;
    if (k == 2) return from_log(-mu - kLn2, lg);
    return d0(lg);
  }
  auto lt = [&](double j) {
    const double nu = k + 2 * j;
    const double ld = nu == 0 ? -kInf : Dist::chisq(nu).d(x, true);
    return sp::dpois_raw(j, mu, true) + ld;
  };
  return from_log(log_series(lt, std::max(mu, k == 0 ? 1.0 : 0.0)), lg);
}

// --------------------------------------------------------------------------------------------------------------
// Noncentral F

// y = m x / (m x + n), yc = n / (m x + n) and their logs, all accurate.
void f_split(double x, double m, double n, double* y, double* yc, double* ly, double* lyc) {
  const double mx = m * x;
  if (mx < n) {
    const double t = mx / n;
    *y = t / (1 + t);
    *yc = 1 / (1 + t);
    *lyc = -std::log1p(t);
    *ly = std::log(m) + std::log(x) - std::log(n) + *lyc;
  } else {
    const double t = n / mx;
    *y = 1 / (1 + t);
    *yc = t / (1 + t);
    *ly = -std::log1p(t);
    *lyc = std::log(n) - std::log(m) - std::log(x) + *ly;
  }
}

double ncf_ltail(double x, double m, double n, double lam, bool lower) {
  double y, yc, ly, lyc;
  f_split(x, m, n, &y, &yc, &ly, &lyc);
  const double mu = 0.5 * lam;
  auto lt = [&](double j) { return sp::dpois_raw(j, mu, true) + sp::beta_inc_tail(0.5 * m + j, 0.5 * n, y, yc, lower, true); };
  return std::min(0.0, log_series(lt, mu));
}

double ncf_p(double x, double m, double n, double lam, bool lower, bool lg) {
  if (std::isinf(n)) return ncchisq_p(x * m, m, lam, lower, lg);
  if (x <= 0) return pstep(true, lower, lg);
  if (x == kInf) return pstep(false, lower, lg);
  const double l = ncf_ltail(x, m, n, lam, lower);
  if (lg && l > -kLn2) return log1mexp(ncf_ltail(x, m, n, lam, !lower));
  return from_log(l, lg);
}

// log of the beta(a, b) density at y (yc = 1 - y), Loader's binomial form for a, b > 2.
double ldbeta(double y, double yc, double ly, double lyc, double a, double b) {
  if (a <= 2 || b <= 2) return (a - 1) * ly + (b - 1) * lyc - sp::lbeta(a, b);
  return std::log(a + b - 1) + sp::dbinom_raw(a - 1, a + b - 2, y, yc, true);
}

double ncf_d(double x, double m, double n, double lam, bool lg) {
  if (std::isinf(n)) {
    const double v = ncchisq_d(x * m, m, lam, true) + std::log(m);
    return from_log(v, lg);
  }
  if (x < 0 || x == kInf) return d0(lg);
  const double mu = 0.5 * lam;
  if (x == 0) {
    if (m < 2) return kInf;
    if (m == 2) return from_log(-mu, lg);
    return d0(lg);
  }
  double y, yc, ly, lyc;
  f_split(x, m, n, &y, &yc, &ly, &lyc);
  const double ljac = ly + lyc - std::log(x);  // dy/dx = y (1 - y) / x
  auto lt = [&](double j) { return sp::dpois_raw(j, mu, true) + ldbeta(y, yc, ly, lyc, 0.5 * m + j, 0.5 * n); };
  return from_log(log_series(lt, mu) + ljac, lg);
}

// --------------------------------------------------------------------------------------------------------------
// Noncentral t

// log P(T > t) for t > 0, finite df, from the chi-square mixture integral (positive integrand).
double nct_upper_integral(double t, double n, double del) {
  const double c = lchisq_w_const(n);
  auto h = [&](double w) { return lchisq_w(w, n, c) + lphi_upper(t * std::exp(0.5 * w) - del); };
  return std::min(0.0, log_integral(h, 0.0, std::min(2.0, std::sqrt(2 / n))));
}

// log P(T <= t) (lower) or log P(T > t) for t > 0 (finite), finite df, d != 0: Lenth's series, with the integral
// where the series cancels (d < 0, upper tail). Only the requested tail is summed.
double nct_tail_pos(double t, double n, double del, bool lower) {
  double x, xc;
  if (t > 1e150) { xc = n / t / t; x = 1; }
  else { const double t2 = t * t; x = t2 / (t2 + n); xc = n / (t2 + n); }
  const double mu = 0.5 * del * del;
  auto lp = [&](double j, bool lw) {
    return sp::dpois_raw(j, mu, true) + sp::beta_inc_tail(j + 0.5, 0.5 * n, x, xc, lw, true);
  };
  auto lq = [&](double j, bool lw) {
    // |q_j| = |d| e^{-mu} mu^j / (sqrt(2) Gamma(j + 3/2)) = e^{-mu} mu^{j+1/2} / Gamma(j + 3/2), mu = d^2/2:
    // a Poisson density at the real argument j + 1/2, evaluated in Loader's saddle-point form.
    return sp::dpois_raw(j + 0.5, mu, true) + sp::beta_inc_tail(j + 1, 0.5 * n, x, xc, lw, true);
  };
  auto series = [&](bool lw, double* sp_, double* sq_) {
    *sp_ = log_series([&](double j) { return lp(j, lw); }, mu) - kLn2;
    *sq_ = log_series([&](double j) { return lq(j, lw); }, mu) - kLn2;
  };
  double a, b;
  if (del > 0) {
    series(lower, &a, &b);
    const double s = logadd(a, b);
    return std::min(0.0, lower ? logadd(lphi_lower(-del), s) : s);
  }
  // d < 0: P(T > t) <= Phi(d) < 1/2; the series is a difference.
  series(false, &a, &b);
  double up;
  if (b - a <= std::log(0.999)) up = a + log1mexp(b - a);
  else up = nct_upper_integral(t, n, del);
  return lower ? log1mexp(up) : up;
}

double nct_tail(double t, double n, double del, bool lower) {
  if (std::isinf(n)) {
    double l, u;
    sp::pnorm_both(t - del, &l, &u, true);
    return lower ? l : u;
  }
  if (t == kInf) return lower ? 0.0 : -kInf;
  if (t == -kInf) return lower ? -kInf : 0.0;
  if (t == 0) return lower ? lphi_lower(-del) : lphi_upper(-del);
  if (t > 0) return nct_tail_pos(t, n, del, lower);
  return nct_tail_pos(-t, n, -del, !lower);  // P(T <= t; d) = P(T >= -t; -d)
}

double nct_p(double t, double n, double del, bool lower, bool lg) {
  double l = nct_tail(t, n, del, lower);
  if (lg && l > -kLn2) l = log1mexp(nct_tail(t, n, del, !lower));
  return from_log(l, lg);
}

double nct_d(double x, double n, double del, bool lg) {
  if (!std::isfinite(x)) return d0(lg);
  if (std::isinf(n)) return from_log(ldnorm(x - del), lg);
  if (x == 0) return from_log(Dist::t(n).d(0, true) - 0.5 * del * del, lg);
  const double c = lchisq_w_const(n);
  auto h = [&](double w) {
    const double v = std::exp(0.5 * w);
    return lchisq_w(w, n, c) + ldnorm(x * v - del) + 0.5 * w;
  };
  return from_log(log_integral(h, 0.0, std::min(2.0, std::sqrt(2 / n))), lg);
}

// --------------------------------------------------------------------------------------------------------------
// Studentized range

// log(Phi(z) - Phi(z - w)), w > 0, computed in the tail that avoids cancellation.
double ldiff_phi(double z, double w) {
  if (z > 0.5 * w) {
    const double l1 = lphi_upper(z - w), l2 = lphi_upper(z);
    return l1 + log1mexp(l2 - l1);
  }
  const double l1 = lphi_lower(z), l2 = lphi_lower(z - w);
  return l1 + log1mexp(l2 - l1);
}

// log W(w): distribution function of the range of k independent standard normals.
double lrange_cdf(double w, double k) {
  if (w <= 0) return -kInf;
  if (w == kInf) return 0;
  const double lk = std::log(k);
  auto h = [&](double z) { return lk + ldnorm(z) + (k - 1) * ldiff_phi(z, w); };
  return std::min(0.0, log_integral(h, 0.5 * w, 1 / std::sqrt(k), 1e-12, 36));
}

double ptukey_lower_log(double q, double k, double df, double r) {
  if (q <= 0) return -kInf;
  if (q == kInf) return 0;
  if (std::isinf(df)) return r * lrange_cdf(q, k);
  const double c = lchisq_w_const(df);
  auto h = [&](double w) { return lchisq_w(w, df, c) + r * lrange_cdf(q * std::exp(0.5 * w), k); };
  return std::min(0.0, log_integral(h, 0.0, std::min(2.0, std::sqrt(2 / df)), 1e-12, 36));
}

double ptukey(double q, double k, double df, double r, bool lower, bool lg) {
  const double l = ptukey_lower_log(q, k, df, r);
  if (lower) return from_log(l, lg);
  return lg ? log1mexp(l) : -std::expm1(l);
}

// --------------------------------------------------------------------------------------------------------------
// Quantiles

// Root of an increasing function g on [umin, umax] by bracketing from u0 (doubling steps) and the Illinois method.
// Returns -Inf/+Inf when the root lies below umin / above umax.
template <class G>
double solve_increasing(G&& g, double u0, double step, double umin, double umax) {
  double a = std::min(std::max(u0, umin), umax), ga = g(a);
  if (ga == 0) return a;
  double b, gb;
  if (ga < 0) {
    for (;;) {
      b = a + step;
      if (b >= umax) { b = umax; gb = g(b); if (gb < 0) return kInf; break; }
      gb = g(b);
      if (gb >= 0) break;
      a = b; ga = gb; step *= 2;
    }
  } else {
    for (;;) {
      b = a - step;
      if (b <= umin) { b = umin; gb = g(b); if (gb > 0) return -kInf; break; }
      gb = g(b);
      if (gb <= 0) break;
      a = b; ga = gb; step *= 2;
    }
  }
  if (gb == 0) return b;
  // a and b bracket the root (g(a), g(b) of opposite signs).
  int side = 0;
  for (int it = 0; it < 400; ++it) {
    const double lo = std::min(a, b), hi = std::max(a, b);
    if (hi - lo <= 4 * DBL_EPSILON * std::max(1.0, std::max(std::fabs(lo), std::fabs(hi)))) break;
    double c = (std::isfinite(ga) && std::isfinite(gb) && ga != gb) ? b - gb * (b - a) / (gb - ga) : 0.5 * (a + b);
    if (!(c > lo && c < hi) || it % 8 == 7) c = 0.5 * (a + b);
    const double gc = g(c);
    if (gc == 0) return c;
    if ((gc < 0) == (gb < 0)) {
      // Illinois step: the retained end point's value is halved when the same side is replaced twice.
      if (side == 1) ga *= 0.5;
      side = 1;
    } else {
      a = b; ga = gb;
      side = 0;
    }
    b = c; gb = gc;
  }
  return std::fabs(ga) < std::fabs(gb) ? a : b;
}

// x with log P_tail(x) = target (tail = lower if use_lower), for a positive distribution given its log tail; solved in
// u = log x from the log-normal start (mlog, slog).
template <class P>
double quantile_positive(P&& ltail, bool use_lower, double target, double mlog, double slog) {
  const double s = use_lower ? 1.0 : -1.0;
  auto g = [&](double u) {
    const double v = s * (ltail(std::exp(u), use_lower) - target);
    return std::isnan(v) ? -1.0 : v;
  };
  const double u = solve_increasing(g, mlog, std::max(0.05, slog), std::log(DBL_TRUE_MIN), std::log(DBL_MAX));
  return std::exp(u);
}

double z_lower(double ll) {
  // Standard normal quantile of exp(ll) via the Dist module.
  return Dist::norm(0, 1).q(ll, true, true);
}

}  // namespace

// ==============================================================================================================
// Dispatch

namespace ncdist_detail {

double d(DistKind kind, const double* par, double x, bool lg) {
  const double a = par[0], b = par[1], c = par[2];
  switch (kind) {
    case DistKind::NcT: return nct_d(x, a, b, lg);
    case DistKind::NcChisq: return ncchisq_d(x, a, b, lg);
    case DistKind::NcF: return ncf_d(x, a, b, c, lg);
    case DistKind::Tukey: throw Error("dtukey is not available (R has no studentized range density)");
    default: break;
  }
  return kNaN;
}

double p(DistKind kind, const double* par, double x, bool lower, bool lg) {
  const double a = par[0], b = par[1], c = par[2];
  switch (kind) {
    case DistKind::NcT: return nct_p(x, a, b, lower, lg);
    case DistKind::NcChisq: return ncchisq_p(x, a, b, lower, lg);
    case DistKind::NcF: return ncf_p(x, a, b, c, lower, lg);
    case DistKind::Tukey: return ptukey(x, a, b, c, lower, lg);
    default: break;
  }
  return kNaN;
}

double q(DistKind kind, const double* par, double p, bool lower, bool lg) {
  const double a = par[0], b = par[1], c = par[2];
  double ll, lu;
  if (std::isnan(p)) return kNaN;
  if (lg) {
    if (p > 0) return kNaN;
    ll = lower ? p : log1mexp(p);
    lu = lower ? log1mexp(p) : p;
  } else {
    if (p < 0 || p > 1) return kNaN;
    ll = lower ? std::log(p) : std::log1p(-p);
    lu = lower ? std::log1p(-p) : std::log(p);
  }
  const bool use_lower = ll <= lu;
  const double target = use_lower ? ll : lu;
  const double z = z_lower(ll);
  switch (kind) {
    case DistKind::NcT: {
      const double n = a, del = b;
      if (ll == -kInf) return -kInf;
      if (lu == -kInf) return kInf;
      if (std::isinf(n)) return del + z;
      const double s = use_lower ? 1.0 : -1.0;
      auto g = [&](double x) {
        const double v = s * (nct_tail(x, n, del, use_lower) - target);
        return std::isnan(v) ? -1.0 : v;
      };
      const double sd = std::sqrt(1 + 0.5 * del * del / n) * (n > 2 ? std::sqrt(n / (n - 2)) : 2.0);
      return solve_increasing(g, del + z * sd, sd, -DBL_MAX, DBL_MAX);
    }
    case DistKind::NcChisq: {
      const double k = a, lam = b;
      if (lu == -kInf) return kInf;
      if (ll == -kInf) return 0;
      if (k == 0 && ll <= -0.5 * lam) return 0;  // inside the point mass at 0
      const double m = k + lam, v = std::log1p(2 * (k + 2 * lam) / (m * m));
      auto lt = [&](double x, bool lw) { return ncchisq_p(x, k, lam, lw, true); };
      return quantile_positive(lt, use_lower, target, std::log(m) - 0.5 * v + z * std::sqrt(v), std::sqrt(v));
    }
    case DistKind::NcF: {
      const double m = a, n = b, lam = c;
      if (lu == -kInf) return kInf;
      if (ll == -kInf) return 0;
      const double mm = m + lam;
      double v = std::log1p(2 * (m + 2 * lam) / (mm * mm));
      if (std::isfinite(n)) v += sp::trigamma(0.5 * n);
      const double mlog = std::log(mm / m) + (std::isfinite(n) && n > 2 ? std::log(n / (n - 2)) : 0.0);
      auto lt = [&](double x, bool lw) { return ncf_p(x, m, n, lam, lw, true); };
      return quantile_positive(lt, use_lower, target, mlog - 0.5 * v + z * std::sqrt(v), std::sqrt(v));
    }
    case DistKind::Tukey: {
      if (lu == -kInf) return kInf;
      if (ll == -kInf) return 0;
      auto lt = [&](double x, bool lw) { return ptukey(x, a, b, c, lw, true); };
      // Start: log of the median range of a normal sample, shifted by z times a rough log-scale spread.
      const double mlog = std::log(std::max(1.0, 2 * std::sqrt(std::log(a + c)) + 0.3)), slog = 0.3;
      return quantile_positive(lt, use_lower, target, mlog + z * slog, slog);
    }
    default: break;
  }
  return kNaN;
}

}  // namespace ncdist_detail

}  // namespace rcli
