// Probability distributions with R's d/p/q semantics (clean-room; written from published mathematics only).
//
// Sources, per algorithm:
// - Normal quantile: M. J. Wichura (1988), Algorithm AS 241 (PPND16), Applied Statistics 37, 477-484, with a Newton
//   refinement on the log scale in the extreme tails.
// - Normal, gamma, beta tails: see special.cpp (continued fractions, Loader's saddle-point densities).
// - Densities of binomial, Poisson, negative binomial, hypergeometric, gamma, beta, F: C. Loader (2000), "Fast and
//   accurate computation of binomial probabilities".
// - Student t density: ratio Gamma((n+1)/2)/Gamma(n/2) through the Stirling remainder (DLMF 5.11.1); distribution
//   function through the incomplete beta function (A&S 26.7.1).
// - Hypergeometric distribution function: summation of the term ratio d(x-1)/d(x) downward from x (direct series),
//   using the smaller tail.
// - Continuous quantiles without closed form (t, chisq, gamma, beta, F): safeguarded Newton iteration on the log of
//   the smaller tail probability, in a transformed variable (log x or logit x), bracketed with bisection fallback.
//   Starting values: Cornish-Fisher / Wilson-Hilferty normal approximations and the log-normal approximation of
//   log X with mean digamma(a) and variance trigamma(a).
// - Discrete quantiles: smallest x with P(X <= x) >= p, found by an exponential then binary search from a
//   Cornish-Fisher start; a relative fuzz of 8 eps (log scale 2.5 eps) as R applies, so q(p(k)) == k.
// - Noncentral t, chi-square, F and the studentized range: see ncdist.cpp (this file only dispatches to it).
#include "rcli/dist.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "rcli/special.hpp"

namespace rcli {

namespace {

namespace sp = rcli::special;

constexpr double kPi = 3.14159265358979323846264338328;
constexpr double kLn2 = 0.693147180559945309417232121458;
constexpr double kLnSqrt2Pi = 0.918938533204672741780329736406;
constexpr double kLnPi = 1.14472988584940017414342735135;

inline double log1mexp(double x) { return x > -kLn2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x)); }
inline double d0(bool lg) { return lg ? -kInf : 0.0; }
inline double d1(bool lg) { return lg ? 0.0 : 1.0; }
inline double dval(double v, bool lg) { return lg ? std::log(v) : v; }
inline double from_log(double l, bool lg) { return lg ? l : std::exp(l); }
// Distribution function at x when x lies entirely below (below = true) or at/above all of the mass.
inline double pstep(bool below, bool lower, bool lg) { return (below == lower) ? d0(lg) : d1(lg); }
inline bool nonint(double x) { return std::fabs(x - std::nearbyint(x)) > 1e-7 * std::max(1.0, std::fabs(x)); }

[[noreturn]] void fail(const std::string& fn, const std::string& what) { throw Error("dist: " + fn + ": " + what); }

// Splits a probability argument into log lower-tail and log upper-tail values (one exact, the other complemented).
// Returns false if p is outside its domain.
bool tails(double p, bool lower, bool lg, double* ll, double* lu) {
  if (std::isnan(p)) return false;
  if (lg) {
    if (p > 0) return false;
    const double other = log1mexp(p);
    *ll = lower ? p : other;
    *lu = lower ? other : p;
  } else {
    if (p < 0 || p > 1) return false;
    const double a = std::log(p), b = std::log1p(-p);
    *ll = lower ? a : b;
    *lu = lower ? b : a;
  }
  return true;
}

// ------------------------------------------------------------------------------------------------------------
// Normal

double dnorm_std(double z, bool lg) {
  if (!std::isfinite(z)) return d0(lg);
  const double az = std::fabs(z);
  if (az < 5) {
    const double l = -0.5 * z * z;
    return lg ? l - kLnSqrt2Pi : 0.398942280401432677939946059934 * std::exp(l);
  }
  const double hi = std::ldexp(std::trunc(std::ldexp(az, 16)), -16), lo = az - hi;
  const double e1 = -0.5 * hi * hi, e2 = -0.5 * lo * (az + hi);
  if (lg) return e1 + e2 - kLnSqrt2Pi;
  return 0.398942280401432677939946059934 * std::exp(e1) * std::exp(e2);
}

// Wichura AS 241 (PPND16).
double qnorm_std(double p, bool lower, bool lg) {
  double ll, lu;
  if (!tails(p, lower, lg, &ll, &lu)) return kNaN;
  if (ll == -kInf) return -kInf;
  if (lu == -kInf) return kInf;
  double pl;  // lower-tail probability (non-log)
  if (lg) pl = lower ? std::exp(p) : -std::expm1(p);
  else pl = lower ? p : (0.5 - p) + 0.5;
  const double q = pl - 0.5;
  if (std::fabs(q) <= 0.425) {
    const double r = 0.180625 - q * q;
    return q *
           (((((((2509.0809287301226727 * r + 33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
             133.14166789178437745) * r + 3.387132872796366608) /
           (((((((5226.495278852545925 * r + 28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
             42.313330701600911252) * r + 1.0);
  }
  const double lt = q < 0 ? ll : lu;
  double r = std::sqrt(-lt), val;
  if (r <= 5) {
    r -= 1.6;
    val = (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r + 0.24178072517745061177) * r +
              1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
            4.6303378461565452959) * r + 1.42343711074968357734) /
          (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
              0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r +
            2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5;
    val = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
              0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r +
            5.4637849111641143699) * r + 6.6579046435011037772) /
          (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
              7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r +
            0.59983220655588793769) * r + 1.0);
    // Newton refinement on the log tail (the rational approximation is tuned for |log p| <~ 700).
    for (int it = 0; it < 4 && std::isfinite(val); ++it) {
      double a, b;
      sp::pnorm_both(-val, &a, &b, true);  // a = log Phi(-val) = log tail
      const double ld = dnorm_std(val, true);
      const double step = (a - lt) / std::exp(ld - a);
      val += step;
      if (std::fabs(step) < 1e-15 * val) break;
    }
  }
  return q < 0 ? -val : val;
}

// ------------------------------------------------------------------------------------------------------------
// Loader-form densities

double dbinom_k(double x, double n, double p, bool lg) {
  if (nonint(x) || x < 0 || !std::isfinite(x)) return d0(lg);
  x = std::nearbyint(x);
  return sp::dbinom_raw(x, n, p, 1 - p, lg);
}

// dbeta with x and y = 1 - x supplied (and optionally log x, log y).
double dbeta_xy(double x, double y, double a, double b, bool lg, const double* lxp = nullptr,
                const double* lyp = nullptr) {
  double l;
  if (a <= 2 || b <= 2) {
    const double lx = lxp ? *lxp : (x < 0.5 ? std::log(x) : std::log1p(-y));
    const double ly = lyp ? *lyp : (y < 0.5 ? std::log(y) : std::log1p(-x));
    l = (a - 1) * lx + (b - 1) * ly - sp::lbeta(a, b);
  } else {
    l = std::log(a + b - 1) + sp::dbinom_raw(a - 1, a + b - 2, x, y, true);
  }
  return from_log(l, lg);
}

double dgamma_r(double x, double shape, double rate, bool lg) {
  if (x < 0) return d0(lg);
  if (shape == 0) return x == 0 ? kInf : d0(lg);
  if (x == 0) {
    if (shape < 1) return kInf;
    if (shape > 1) return d0(lg);
    return dval(rate, lg);
  }
  if (x == kInf) return d0(lg);
  if (shape < 1) {
    const double l = sp::dpois_raw(shape, x * rate, true) + std::log(shape / x);
    return from_log(l, lg);
  }
  const double l = sp::dpois_raw(shape - 1, x * rate, true) + std::log(rate);
  return from_log(l, lg);
}

double pgamma_r(double x, double shape, double rate, bool lower, bool lg) {
  if (shape == 0) return pstep(x < 0, lower, lg);  // point mass at 0
  if (x <= 0) return pstep(true, lower, lg);
  return sp::gamma_inc(shape, x * rate, lower, lg);
}

// ------------------------------------------------------------------------------------------------------------
// Generic continuous quantile solver.

enum class Domain { Positive, Unit };

// Find x in the domain with log P_tail(x) = target (tail lower if use_lower), starting at u0 in the transformed
// variable (log x or logit x).
double solve_quantile(const Dist& D, bool use_lower, double target, Domain dom, double u0) {
  auto X = [&](double u) {
    if (dom == Domain::Positive) return std::exp(u);
    return u >= 0 ? 1 / (1 + std::exp(-u)) : std::exp(u) / (1 + std::exp(u));
  };
  const double s = use_lower ? 1.0 : -1.0;
  if (dom == Domain::Positive) {
    // Quantiles beyond the double range: R returns Inf (or 0).
    if (s * (D.p(DBL_MAX, use_lower, true) - target) < 0) return kInf;
    if (s * (D.p(DBL_TRUE_MIN, use_lower, true) - target) > 0) return 0;
  }
  double lo = -kInf, hi = kInf, u = std::isfinite(u0) ? u0 : 0.0, step = 1;
  double width_prev = kInf, width_prev2 = kInf;
  for (int it = 0; it < 2000; ++it) {
    const double x = X(u);
    const double lP = D.p(x, use_lower, true);
    const double h = s * (lP - target);
    if (h == 0) return x;
    if (h < 0 || std::isnan(h)) lo = std::max(lo, u);
    else hi = std::min(hi, u);
    const bool bracketed = std::isfinite(lo) && std::isfinite(hi);
    if (bracketed && hi - lo <= 2 * DBL_EPSILON * std::max(1.0, std::fabs(u))) break;
    const double dxdu = dom == Domain::Positive ? x : x * (1 - x);
    const double deriv = std::exp(D.d(x, true) - lP) * dxdu;
    double un = u - h / deriv;
    if (!bracketed) {
      // Expand geometrically until the root is bracketed; Newton only if it moves farther.
      const double dir = (h < 0 || std::isnan(h)) ? 1.0 : -1.0;
      const double nstep = std::isfinite(un) ? (un - u) * dir : -1;
      un = u + dir * std::max(step, nstep);
      step *= 2;
    } else {
      const double width = hi - lo;
      const bool slow = width > 0.5 * width_prev2;
      if (!std::isfinite(un) || un <= lo || un >= hi || slow) un = 0.5 * (lo + hi);
      width_prev2 = width_prev;
      width_prev = width;
    }
    if (un == u) break;
    const double du = std::fabs(un - u);
    u = un;
    if (du <= 1e-15 * std::max(1.0, std::fabs(u)) && bracketed) return X(u);
  }
  return X(u);
}

// ------------------------------------------------------------------------------------------------------------
// Generic discrete quantile: smallest integer x in [lo, hi] with P(X <= x) >= p (with fuzz).
double discrete_quantile(const Dist& D, double p, bool lower, bool lg, double lo, double hi, double mean, double sd,
                         double skew) {
  double ll, lu;
  if (!tails(p, lower, lg, &ll, &lu)) return kNaN;
  if (ll == -kInf) return lo;
  if (lu == -kInf) return hi;
  // Each request is compared in the tail it names (R's semantics, so q(p(k)) == k). For an upper-tail p above 1/2 the
  // fuzz is applied to the small complement 1 - p (exact there): applied to p itself, p * (1 + fuzz) >= 1 near 1 let
  // every x qualify, so qbinom(1 - 1.3e-15, 533, 0.885, lower.tail = FALSE) returned 0 instead of R's 407.
  const double upper_cap = lg ? std::log1p(-std::exp(ll) * (1 - 8 * DBL_EPSILON)) : 1 - (1 - p) * (1 - 8 * DBL_EPSILON);
  auto ok = [&](double x) {
    if (x >= hi) return true;
    if (x < lo) return false;
    if (!lg) {
      if (lower) return D.p(x, true, false) >= p * (1 - 8 * DBL_EPSILON);
      if (p > 0.5) return D.p(x, false, false) <= upper_cap;
      return D.p(x, false, false) <= p * (1 + 8 * DBL_EPSILON);
    }
    if (lower) return D.p(x, true, true) >= p * (1 + 2.5 * DBL_EPSILON);
    if (p > -M_LN2) return D.p(x, false, true) <= upper_cap;
    return D.p(x, false, true) <= p * (1 - 2.5 * DBL_EPSILON);
  };
  // Cornish-Fisher start.
  const double z = qnorm_std(ll, true, true);
  double g = std::floor(mean + sd * (z + skew * (z * z - 1) / 6) + 0.5);
  if (!std::isfinite(g)) g = mean;
  g = std::min(std::max(g, lo), std::isfinite(hi) ? hi : g);
  double a, b;  // invariant: !ok(a), ok(b)
  if (ok(g)) {
    b = g;
    double stepd = 1;
    for (;;) {
      a = b - stepd;
      if (a < lo) { a = lo - 1; break; }
      if (!ok(a)) break;
      b = a;
      stepd *= 2;
    }
  } else {
    a = g;
    double stepd = 1;
    for (;;) {
      b = a + stepd;
      if (b >= hi) { b = hi; break; }
      if (ok(b)) break;
      a = b;
      stepd *= 2;
    }
  }
  while (b - a > 1) {
    const double m = std::floor(a + (b - a) / 2);
    if (ok(m)) b = m;
    else a = m;
  }
  return b;
}

// ------------------------------------------------------------------------------------------------------------
// t distribution helpers

double dt_std(double x, double n, bool lg) {
  if (std::isinf(n)) return dnorm_std(x, lg);
  if (!std::isfinite(x)) return d0(lg);
  const double a = 0.5 * n;
  // lr = lgamma(a + 1/2) - lgamma(a) - 0.5 log a; for large a through the Stirling remainder
  // w(z) = lgamma(z) - [(z - 1/2) log z - z + log sqrt(2 pi)] (= stirlerr(z)):
  // lr = a log1pmx(1/(2a)) + w(a + 1/2) - w(a).
  const double lr = (a < 10) ? sp::lgamma(a + 0.5) - sp::lgamma(a) - 0.5 * std::log(a)
                             : a * sp::log1pmx(0.5 / a) + sp::stirlerr(a + 0.5) - sp::stirlerr(a);
  const double t = x * x / n;
  double l1p;
  if (!std::isfinite(t) || t > 1e300) l1p = 2 * std::log(std::fabs(x)) - std::log(n);
  else l1p = std::log1p(t);
  const double l = lr - kLnSqrt2Pi - 0.5 * (n + 1) * l1p;
  return from_log(l, lg);
}

// log P(T <= x) and log P(T > x).
void pt_both(double x, double n, double* ll, double* lu) {
  if (std::isinf(n)) {
    sp::pnorm_both(x, ll, lu, true);
    return;
  }
  if (std::isnan(x)) { *ll = *lu = x; return; }
  if (x == 0) { *ll = *lu = -kLn2; return; }
  const double ax = std::fabs(x);
  double lI, lIc;  // log P(|T| > ax), log P(|T| < ax)
  if (ax == kInf) { lI = -kInf; lIc = 0; }
  else if (ax > 1e150) {
    lI = 0.5 * n * (std::log(n) - 2 * std::log(ax)) - std::log(0.5 * n) - sp::lbeta(0.5 * n, 0.5);
    lIc = log1mexp(lI);
  } else {
    double xb, yb;
    const double t = ax * ax / n;
    if (t < 1) { xb = 1 / (1 + t); yb = t / (1 + t); }
    else { const double r = n / (ax * ax); xb = r / (1 + r); yb = 1 / (1 + r); }
    lI = sp::beta_inc_tail(0.5 * n, 0.5, xb, yb, true, true);
    lIc = sp::beta_inc_tail(0.5 * n, 0.5, xb, yb, false, true);
  }
  const double tail = lI - kLn2;                       // one-sided tail beyond ax
  const double big = log1mexp(tail);  // 1 - tail
  (void)lIc;
  if (x < 0) { *ll = tail; *lu = big; }
  else { *lu = tail; *ll = big; }
}

double cf_start(double z, double n) {
  // Cornish-Fisher expansion of the t quantile in terms of the normal quantile z (A&S 26.7.5).
  const double z2 = z * z;
  return z + z * (z2 + 1) / (4 * n) + z * ((5 * z2 + 16) * z2 + 3) / (96 * n * n) +
         z * (((3 * z2 + 19) * z2 + 17) * z2 - 15) / (384 * n * n * n);
}

// F helper: x_b = m x/(n + m x), y_b = n/(n + m x), both accurate, and optionally their logs.
void f_split(double x, double m, double n, double* xb, double* yb, double* lxb = nullptr, double* lyb = nullptr) {
  const double mx = m * x;
  double lx, ly;
  if (mx < n) {
    const double t = mx / n;
    *xb = t / (1 + t);
    *yb = 1 / (1 + t);
    ly = -std::log1p(t);
    lx = std::log(m) + std::log(x) - std::log(n) + ly;
  } else {
    const double t = n / mx;
    *xb = 1 / (1 + t);
    *yb = t / (1 + t);
    lx = -std::log1p(t);
    ly = std::log(n) - std::log(m) - std::log(x) + lx;
  }
  if (lxb) *lxb = lx;
  if (lyb) *lyb = ly;
}

// hypergeometric: log P(X <= x) via the downward series, x <= mean preferably.
double phyper_lower_log(double x, double m, double n, double k) {
  // d(x-1)/d(x) = x (n - k + x) / ((m - x + 1)(k - x + 1))
  double sum = 0, term = 1;
  double xx = x;
  while (xx > 0 && term >= DBL_EPSILON * sum) {
    term *= xx * (n - k + xx) / ((k + 1 - xx) * (m + 1 - xx));
    sum += term;
    xx -= 1;
  }
  const double pp = k / (m + n), qq = (m + n - k) / (m + n);
  const double ld = sp::dbinom_raw(x, m, pp, qq, true) + sp::dbinom_raw(k - x, n, pp, qq, true) -
                    sp::dbinom_raw(k, m + n, pp, qq, true);
  return ld + std::log1p(sum);
}

}  // namespace

// ==============================================================================================================
// Constructors

Dist Dist::norm(double mean, double sd) {
  if (!std::isfinite(mean) || std::isnan(sd) || sd < 0 || std::isinf(sd)) fail("norm", "need finite mean and sd >= 0");
  return Dist(DistKind::Norm, mean, sd, kNaN);
}
Dist Dist::lnorm(double meanlog, double sdlog) {
  if (!std::isfinite(meanlog) || std::isnan(sdlog) || sdlog < 0 || std::isinf(sdlog))
    fail("lnorm", "need finite meanlog and sdlog >= 0");
  return Dist(DistKind::Lnorm, meanlog, sdlog, kNaN);
}
Dist Dist::t(double df) {
  if (std::isnan(df) || df <= 0) fail("t", "df must be > 0");
  return Dist(DistKind::T, df, kNaN, kNaN);
}
Dist Dist::chisq(double df) {
  if (std::isnan(df) || df < 0 || std::isinf(df)) fail("chisq", "df must be finite and >= 0");
  return Dist(DistKind::Chisq, df, kNaN, kNaN);
}
Dist Dist::f(double df1, double df2) {
  if (std::isnan(df1) || std::isnan(df2) || df1 <= 0 || df2 <= 0) fail("f", "df1 and df2 must be > 0");
  return Dist(DistKind::F, df1, df2, kNaN);
}
Dist Dist::gamma(double shape, double rate) {
  if (std::isnan(shape) || std::isnan(rate) || shape < 0 || std::isinf(shape) || rate <= 0 || std::isinf(rate))
    fail("gamma", "need finite shape >= 0 and finite rate > 0");
  return Dist(DistKind::Gamma, shape, rate, kNaN);
}
Dist Dist::beta(double shape1, double shape2) {
  if (std::isnan(shape1) || std::isnan(shape2) || shape1 < 0 || shape2 < 0) fail("beta", "shapes must be >= 0");
  return Dist(DistKind::Beta, shape1, shape2, kNaN);
}
Dist Dist::exp(double rate) {
  if (std::isnan(rate) || rate <= 0 || std::isinf(rate)) fail("exp", "rate must be finite and > 0");
  return Dist(DistKind::Exp, rate, kNaN, kNaN);
}
Dist Dist::unif(double min, double max) {
  if (!std::isfinite(min) || !std::isfinite(max) || min > max) fail("unif", "need finite min <= max");
  return Dist(DistKind::Unif, min, max, kNaN);
}
Dist Dist::weibull(double shape, double scale) {
  if (std::isnan(shape) || std::isnan(scale) || shape <= 0 || scale <= 0 || std::isinf(shape) || std::isinf(scale))
    fail("weibull", "shape and scale must be finite and > 0");
  return Dist(DistKind::Weibull, shape, scale, kNaN);
}
Dist Dist::cauchy(double location, double scale) {
  if (!std::isfinite(location) || std::isnan(scale) || scale <= 0 || std::isinf(scale))
    fail("cauchy", "need finite location and finite scale > 0");
  return Dist(DistKind::Cauchy, location, scale, kNaN);
}
Dist Dist::logis(double location, double scale) {
  if (!std::isfinite(location) || std::isnan(scale) || scale <= 0 || std::isinf(scale))
    fail("logis", "need finite location and finite scale > 0");
  return Dist(DistKind::Logis, location, scale, kNaN);
}
Dist Dist::binom(double size, double prob) {
  if (std::isnan(size) || size < 0 || !std::isfinite(size) || nonint(size)) fail("binom", "size must be a non-negative integer");
  if (std::isnan(prob) || prob < 0 || prob > 1) fail("binom", "prob must be in [0, 1]");
  return Dist(DistKind::Binom, std::nearbyint(size), prob, kNaN);
}
Dist Dist::pois(double lambda) {
  if (std::isnan(lambda) || lambda < 0 || std::isinf(lambda)) fail("pois", "lambda must be finite and >= 0");
  return Dist(DistKind::Pois, lambda, kNaN, kNaN);
}
Dist Dist::geom(double prob) {
  if (std::isnan(prob) || prob <= 0 || prob > 1) fail("geom", "prob must be in (0, 1]");
  return Dist(DistKind::Geom, prob, kNaN, kNaN);
}
Dist Dist::nbinom(double size, double prob) {
  if (std::isnan(size) || size < 0 || std::isinf(size)) fail("nbinom", "size must be finite and >= 0");
  if (std::isnan(prob) || prob <= 0 || prob > 1) fail("nbinom", "prob must be in (0, 1]");
  return Dist(DistKind::Nbinom, size, prob, kNaN);
}
Dist Dist::hyper(double m, double n, double k) {
  for (double v : {m, n, k})
    if (std::isnan(v) || v < 0 || !std::isfinite(v) || nonint(v)) fail("hyper", "m, n and k must be non-negative integers");
  m = std::nearbyint(m);
  n = std::nearbyint(n);
  k = std::nearbyint(k);
  if (k > m + n) fail("hyper", "k must be <= m + n");
  return Dist(DistKind::Hyper, m, n, k);
}

Dist Dist::t(double df, double ncp) {
  if (std::isnan(df) || df <= 0) fail("t", "df must be > 0");
  if (!std::isfinite(ncp)) fail("t", "ncp must be finite");
  if (ncp == 0) return t(df);
  return Dist(DistKind::NcT, df, ncp, kNaN);
}
Dist Dist::chisq(double df, double ncp) {
  if (std::isnan(df) || df < 0 || std::isinf(df)) fail("chisq", "df must be finite and >= 0");
  if (std::isnan(ncp) || ncp < 0 || std::isinf(ncp)) fail("chisq", "ncp must be finite and >= 0");
  if (ncp == 0) return chisq(df);
  return Dist(DistKind::NcChisq, df, ncp, kNaN);
}
Dist Dist::f(double df1, double df2, double ncp) {
  if (std::isnan(df1) || std::isnan(df2) || df1 <= 0 || df2 <= 0) fail("f", "df1 and df2 must be > 0");
  if (std::isnan(ncp) || ncp < 0 || std::isinf(ncp)) fail("f", "ncp must be finite and >= 0");
  if (ncp == 0) return f(df1, df2);
  if (std::isinf(df1)) fail("f", "df1 must be finite when ncp > 0");
  return Dist(DistKind::NcF, df1, df2, ncp);
}
Dist Dist::tukey(double nmeans, double df, double nranges) {
  if (std::isnan(nmeans) || nmeans < 2 || std::isinf(nmeans)) fail("tukey", "nmeans must be finite and >= 2");
  if (std::isnan(df) || df < 2) fail("tukey", "df must be >= 2");
  if (std::isnan(nranges) || nranges < 1 || std::isinf(nranges)) fail("tukey", "nranges must be finite and >= 1");
  return Dist(DistKind::Tukey, nmeans, df, nranges);
}

Dist Dist::by_name(const std::string& name, double p1, double p2, double p3) {
  auto need = [&](int k) {
    const double v[3] = {p1, p2, p3};
    for (int i = 0; i < k; ++i)
      if (std::isnan(v[i])) fail(name, "missing parameter " + std::to_string(i + 1));
  };
  if (name == "norm") { if (std::isnan(p1)) p1 = 0; if (std::isnan(p2)) p2 = 1; return norm(p1, p2); }
  if (name == "lnorm") { if (std::isnan(p1)) p1 = 0; if (std::isnan(p2)) p2 = 1; return lnorm(p1, p2); }
  if (name == "t") { need(1); return std::isnan(p2) ? t(p1) : t(p1, p2); }
  if (name == "chisq") { need(1); return std::isnan(p2) ? chisq(p1) : chisq(p1, p2); }
  if (name == "f") { need(2); return std::isnan(p3) ? f(p1, p2) : f(p1, p2, p3); }
  if (name == "tukey") { need(2); return tukey(p1, p2, std::isnan(p3) ? 1.0 : p3); }
  if (name == "gamma") { need(1); if (std::isnan(p2)) p2 = 1; return gamma(p1, p2); }
  if (name == "beta") { need(2); return beta(p1, p2); }
  if (name == "exp") { if (std::isnan(p1)) p1 = 1; return exp(p1); }
  if (name == "unif") { if (std::isnan(p1)) p1 = 0; if (std::isnan(p2)) p2 = 1; return unif(p1, p2); }
  if (name == "weibull") { need(1); if (std::isnan(p2)) p2 = 1; return weibull(p1, p2); }
  if (name == "cauchy") { if (std::isnan(p1)) p1 = 0; if (std::isnan(p2)) p2 = 1; return cauchy(p1, p2); }
  if (name == "logis") { if (std::isnan(p1)) p1 = 0; if (std::isnan(p2)) p2 = 1; return logis(p1, p2); }
  if (name == "binom") { need(2); return binom(p1, p2); }
  if (name == "pois") { need(1); return pois(p1); }
  if (name == "geom") { need(1); return geom(p1); }
  if (name == "nbinom") { need(2); return nbinom(p1, p2); }
  if (name == "hyper") { need(3); return hyper(p1, p2, p3); }
  throw Error("dist: unknown distribution '" + name + "'");
}

std::string Dist::name() const {
  switch (kind_) {
    case DistKind::Norm: return "norm";
    case DistKind::Lnorm: return "lnorm";
    case DistKind::T: return "t";
    case DistKind::Chisq: return "chisq";
    case DistKind::F: return "f";
    case DistKind::Gamma: return "gamma";
    case DistKind::Beta: return "beta";
    case DistKind::Exp: return "exp";
    case DistKind::Unif: return "unif";
    case DistKind::Weibull: return "weibull";
    case DistKind::Cauchy: return "cauchy";
    case DistKind::Logis: return "logis";
    case DistKind::Binom: return "binom";
    case DistKind::Pois: return "pois";
    case DistKind::Geom: return "geom";
    case DistKind::Nbinom: return "nbinom";
    case DistKind::Hyper: return "hyper";
    case DistKind::NcT: return "t";
    case DistKind::NcChisq: return "chisq";
    case DistKind::NcF: return "f";
    case DistKind::Tukey: return "tukey";
  }
  return "?";
}

bool Dist::discrete() const {
  switch (kind_) {
    case DistKind::Binom: case DistKind::Pois: case DistKind::Geom: case DistKind::Nbinom: case DistKind::Hyper:
      return true;
    default:
      return false;
  }
}

// ==============================================================================================================
// Densities

double Dist::d(double x, bool lg) const {
  if (std::isnan(x)) return x;
  const double a = par_[0], b = par_[1], c = par_[2];
  switch (kind_) {
    case DistKind::NcT: case DistKind::NcChisq: case DistKind::NcF: case DistKind::Tukey:
      return ncdist_detail::d(kind_, par_, x, lg);
    case DistKind::Norm: {
      if (b == 0) return x == a ? kInf : d0(lg);
      const double v = dnorm_std((x - a) / b, lg);
      return lg ? v - std::log(b) : v / b;
    }
    case DistKind::Lnorm: {
      if (x <= 0) return d0(lg);
      if (b == 0) return std::log(x) == a ? kInf : d0(lg);
      const double lx = std::log(x);
      const double v = dnorm_std((lx - a) / b, true) - std::log(b) - lx;
      return from_log(v, lg);
    }
    case DistKind::T: return dt_std(x, a, lg);
    case DistKind::Chisq: return dgamma_r(x, 0.5 * a, 0.5, lg);
    case DistKind::Gamma: return dgamma_r(x, a, b, lg);
    case DistKind::F: {
      const double m = a, n = b;
      if (x < 0) return d0(lg);
      if (x == 0) return m > 2 ? d0(lg) : (m == 2 ? d1(lg) : kInf);
      if (std::isinf(m) && std::isinf(n)) return x == 1 ? kInf : d0(lg);
      if (std::isinf(n)) return dgamma_r(x, 0.5 * m, 0.5 * m, lg);
      if (std::isinf(m)) {
        const double v = dgamma_r(1 / x, 0.5 * n, 0.5 * n, true) - 2 * std::log(x);
        return from_log(v, lg);
      }
      if (std::isinf(x)) return d0(lg);
      double xb, yb, lxb, lyb;
      f_split(x, m, n, &xb, &yb, &lxb, &lyb);
      // dF(x) = dbeta(xb; m/2, n/2) * xb * yb / x, with xb / x = m / (n + m x)
      const double lfac = std::log(m / n) + lyb;
      const double v = dbeta_xy(xb, yb, 0.5 * m, 0.5 * n, true, &lxb, &lyb) + lfac + lyb;
      return from_log(v, lg);
    }
    case DistKind::Beta: {
      if (x < 0 || x > 1) return d0(lg);
      if (a == 0 || b == 0 || std::isinf(a) || std::isinf(b)) {
        if (a == 0 && b == 0) return (x == 0 || x == 1) ? kInf : d0(lg);
        if (a == 0 || std::isinf(b)) return x == 0 ? kInf : d0(lg);
        if (b == 0 || std::isinf(a)) return x == 1 ? kInf : d0(lg);
        return x == 0.5 ? kInf : d0(lg);
      }
      if (x == 0) {
        if (a > 1) return d0(lg);
        if (a < 1) return kInf;
        return dval(b, lg);
      }
      if (x == 1) {
        if (b > 1) return d0(lg);
        if (b < 1) return kInf;
        return dval(a, lg);
      }
      return dbeta_xy(x, 1 - x, a, b, lg);
    }
    case DistKind::Exp: return x < 0 ? d0(lg) : (lg ? std::log(a) - a * x : a * std::exp(-a * x));
    case DistKind::Unif:
      if (x < a || x > b) return d0(lg);
      if (a == b) return kInf;
      return lg ? -std::log(b - a) : 1 / (b - a);
    case DistKind::Weibull: {
      if (x < 0 || std::isinf(x)) return d0(lg);
      const double z = x / b;
      if (x == 0) return a < 1 ? kInf : (a > 1 ? d0(lg) : dval(1 / b, lg));
      const double tk = std::pow(z, a);
      const double v = std::log(a / b) + (a - 1) * std::log(z) - tk;
      return lg ? v : (a / b) * std::pow(z, a - 1) * std::exp(-tk);
    }
    case DistKind::Cauchy: {
      const double z = (x - a) / b;
      if (std::isinf(z)) return d0(lg);
      if (lg) {
        const double az = std::fabs(z);
        const double l1 = az > 1e150 ? 2 * std::log(az) : std::log1p(z * z);
        return -(kLnPi + std::log(b) + l1);
      }
      return 1 / (kPi * b * (1 + z * z));
    }
    case DistKind::Logis: {
      const double z = (x - a) / b;
      if (std::isinf(z)) return d0(lg);
      const double az = std::fabs(z), e = std::exp(-az);
      if (lg) return -az - std::log(b) - 2 * std::log1p(e);
      const double f = 1 + e;
      return e / (b * f * f);
    }
    case DistKind::Binom: return dbinom_k(x, a, b, lg);
    case DistKind::Pois: {
      if (nonint(x) || x < 0 || !std::isfinite(x)) return d0(lg);
      return sp::dpois_raw(std::nearbyint(x), a, lg);
    }
    case DistKind::Geom: {
      if (nonint(x) || x < 0 || !std::isfinite(x)) return d0(lg);
      x = std::nearbyint(x);
      if (a == 1) return x == 0 ? d1(lg) : d0(lg);
      const double l = std::log(a) + x * std::log1p(-a);
      return from_log(l, lg);
    }
    case DistKind::Nbinom: {
      if (nonint(x) || x < 0 || !std::isfinite(x)) return d0(lg);
      x = std::nearbyint(x);
      if (x == 0 && a == 0) return d1(lg);
      if (a == 0) return d0(lg);
      if (b == 1) return x == 0 ? d1(lg) : d0(lg);
      const double l = sp::dbinom_raw(a, x + a, b, 1 - b, true) + std::log(a / (a + x));
      return from_log(l, lg);
    }
    case DistKind::Hyper: {
      const double m = a, n = b, k = c;
      if (nonint(x) || x < 0 || !std::isfinite(x)) return d0(lg);
      x = std::nearbyint(x);
      if (x < std::max(0.0, k - n) || x > std::min(k, m)) return d0(lg);
      if (k == 0) return x == 0 ? d1(lg) : d0(lg);
      const double pp = k / (m + n), qq = (m + n - k) / (m + n);
      const double l = sp::dbinom_raw(x, m, pp, qq, true) + sp::dbinom_raw(k - x, n, pp, qq, true) -
                       sp::dbinom_raw(k, m + n, pp, qq, true);
      return from_log(l, lg);
    }
  }
  return kNaN;
}

// ==============================================================================================================
// Distribution functions

double Dist::p(double x, bool lower, bool lg) const {
  if (std::isnan(x)) return x;
  const double a = par_[0], b = par_[1], c = par_[2];
  switch (kind_) {
    case DistKind::NcT: case DistKind::NcChisq: case DistKind::NcF: case DistKind::Tukey:
      return ncdist_detail::p(kind_, par_, x, lower, lg);
    case DistKind::Norm: {
      if (b == 0) return pstep(x < a, lower, lg);
      const double z = (x - a) / b;
      double l, u;
      sp::pnorm_both(z, &l, &u, lg);
      return lower ? l : u;
    }
    case DistKind::Lnorm: {
      if (x <= 0) return pstep(true, lower, lg);
      if (b == 0) return pstep(std::log(x) < a, lower, lg);
      double l, u;
      sp::pnorm_both((std::log(x) - a) / b, &l, &u, lg);
      return lower ? l : u;
    }
    case DistKind::T: {
      double l, u;
      pt_both(x, a, &l, &u);
      return from_log(lower ? l : u, lg);
    }
    case DistKind::Chisq: return pgamma_r(x, 0.5 * a, 0.5, lower, lg);
    case DistKind::Gamma: return pgamma_r(x, a, b, lower, lg);
    case DistKind::F: {
      const double m = a, n = b;
      if (x <= 0) return pstep(true, lower, lg);
      if (x == kInf) return pstep(false, lower, lg);
      if (std::isinf(m) && std::isinf(n)) return pstep(x < 1, lower, lg);
      if (std::isinf(n)) return pgamma_r(x * m, 0.5 * m, 0.5, lower, lg);
      if (std::isinf(m)) return pgamma_r(n / x, 0.5 * n, 0.5, !lower, lg);
      double xb, yb;
      f_split(x, m, n, &xb, &yb);
      return sp::beta_inc_tail(0.5 * m, 0.5 * n, xb, yb, lower, lg);
    }
    case DistKind::Beta: {
      if (x <= 0) return pstep(true, lower, lg);
      if (x >= 1) return pstep(false, lower, lg);
      return sp::beta_inc_tail(a, b, x, 1 - x, lower, lg);
    }
    case DistKind::Exp: {
      if (x <= 0) return pstep(true, lower, lg);
      const double lu = -a * x;
      if (lower) return lg ? log1mexp(lu) : -std::expm1(lu);
      return lg ? lu : std::exp(lu);
    }
    case DistKind::Unif: {
      if (x < a) return pstep(true, lower, lg);
      if (x >= b) return pstep(false, lower, lg);
      const double v = lower ? (x - a) / (b - a) : (b - x) / (b - a);
      return dval(v, lg);
    }
    case DistKind::Weibull: {
      if (x <= 0) return pstep(true, lower, lg);
      const double lu = -std::pow(x / b, a);
      if (lower) return lg ? log1mexp(lu) : -std::expm1(lu);
      return lg ? lu : std::exp(lu);
    }
    case DistKind::Cauchy: {
      double z = (x - a) / b;
      if (!lower) z = -z;  // P(X > x) = P(-Z <= -z)
      // lower tail of standard Cauchy at z
      if (std::isinf(z)) return pstep(z < 0, true, lg);
      if (z < -1) return dval(std::atan(-1 / z) / kPi, lg);
      if (z > 1) {
        const double t = std::atan(1 / z) / kPi;  // upper tail
        return lg ? std::log1p(-t) : 1 - t;
      }
      return dval(0.5 + std::atan(z) / kPi, lg);
    }
    case DistKind::Logis: {
      double z = (x - a) / b;
      if (!lower) z = -z;
      if (std::isinf(z)) return pstep(z < 0, true, lg);
      if (lg) return z < 0 ? z - std::log1p(std::exp(z)) : -std::log1p(std::exp(-z));
      if (z < 0) { const double e = std::exp(z); return e / (1 + e); }
      return 1 / (1 + std::exp(-z));
    }
    case DistKind::Binom: {
      if (x < 0) return pstep(true, lower, lg);
      const double k = std::floor(x + 1e-7);
      if (k >= a) return pstep(false, lower, lg);
      if (b == 0) return pstep(false, lower, lg);
      if (b == 1) return pstep(true, lower, lg);
      // P(X <= k) = 1 - I_p(k + 1, n - k)
      return sp::beta_inc_tail(k + 1, a - k, b, 1 - b, !lower, lg);
    }
    case DistKind::Pois: {
      if (x < 0) return pstep(true, lower, lg);
      if (a == 0 || std::isinf(x)) return pstep(false, lower, lg);
      const double k = std::floor(x + 1e-7);
      return sp::gamma_inc(k + 1, a, !lower, lg);
    }
    case DistKind::Geom: {
      if (x < 0) return pstep(true, lower, lg);
      if (a == 1 || std::isinf(x)) return pstep(false, lower, lg);
      const double k = std::floor(x + 1e-7);
      const double lu = (k + 1) * std::log1p(-a);
      if (lower) return lg ? log1mexp(lu) : -std::expm1(lu);
      return lg ? lu : std::exp(lu);
    }
    case DistKind::Nbinom: {
      if (x < 0) return pstep(true, lower, lg);
      if (a == 0 || b == 1 || std::isinf(x)) return pstep(false, lower, lg);
      const double k = std::floor(x + 1e-7);
      return sp::beta_inc_tail(a, k + 1, b, 1 - b, lower, lg);
    }
    case DistKind::Hyper: {
      const double m = a, n = b, k = c;
      const double xx = std::floor(x + 1e-7);
      if (xx < std::max(0.0, k - n)) return pstep(true, lower, lg);
      if (xx >= std::min(k, m)) return pstep(false, lower, lg);
      // Use the smaller tail: if x is above the mean, swap colours (X' = k - X ~ hyper(n, m, k)).
      if (xx * (m + n) > k * m) {
        const double lu = phyper_lower_log(k - xx - 1, n, m, k);  // P(X > x) = P(X' <= k - x - 1)
        return from_log(lower ? log1mexp(lu) : lu, lg);
      }
      const double ll = phyper_lower_log(xx, m, n, k);
      return from_log(lower ? ll : log1mexp(ll), lg);
    }
  }
  return kNaN;
}

// ==============================================================================================================
// Quantiles

double Dist::q(double p, bool lower, bool lg) const {
  const double a = par_[0], b = par_[1], c = par_[2];
  double ll, lu;
  if (!tails(p, lower, lg, &ll, &lu)) return kNaN;
  switch (kind_) {
    case DistKind::NcT: case DistKind::NcChisq: case DistKind::NcF: case DistKind::Tukey:
      return ncdist_detail::q(kind_, par_, p, lower, lg);
    case DistKind::Norm:
      if (b == 0) return a;
      return a + b * qnorm_std(p, lower, lg);
    case DistKind::Lnorm:
      if (b == 0) return std::exp(a);
      return std::exp(a + b * qnorm_std(p, lower, lg));
    case DistKind::Exp:
      if (ll == -kInf) return 0;
      return -lu / a;
    case DistKind::Unif: {
      if (lg || !lower) {
        const double pl = std::exp(ll);
        return a + pl * (b - a);
      }
      return a + p * (b - a);
    }
    case DistKind::Weibull:
      if (ll == -kInf) return 0;
      return b * std::pow(-lu, 1 / a);
    case DistKind::Cauchy: {
      if (ll == -kInf) return -kInf;
      if (lu == -kInf) return kInf;
      const bool neg = ll < lu;
      const double t = std::exp(neg ? ll : lu);  // smaller tail, <= 0.5
      if (t == 0.5) return a;
      // x = 1 / tan(pi t) on the side of the smaller tail
      const double s = std::sin(kPi * t), co = std::cos(kPi * t);
      const double v = co / s;
      return a + b * (neg ? -v : v);
    }
    case DistKind::Logis:
      if (ll == -kInf) return -kInf;
      if (lu == -kInf) return kInf;
      return a + b * (ll - lu);
    case DistKind::T: {
      const double n = a;
      if (ll == -kInf) return -kInf;
      if (lu == -kInf) return kInf;
      if (std::isinf(n)) return qnorm_std(p, lower, lg);
      if (ll == lu) return 0;
      const bool neg = ll < lu;
      const double lt = neg ? ll : lu;  // log of the smaller tail
      if (n == 1) {
        const double t = std::exp(lt);
        const double v = std::cos(kPi * t) / std::sin(kPi * t);
        return neg ? -v : v;
      }
      const double z = -qnorm_std(lt, true, true);  // >= 0
      double y0 = cf_start(z, n);
      const double lpow = ((0.5 * n - 1) * std::log(n) - sp::lbeta(0.5 * n, 0.5) - lt) / n;
      double u0 = (y0 > 0) ? std::log(y0) : lpow;
      if (lpow > u0) u0 = lpow;
      if (!std::isfinite(u0)) u0 = 0;
      const double y = solve_quantile(*this, false, lt, Domain::Positive, u0);
      return neg ? -y : y;
    }
    case DistKind::Chisq:
    case DistKind::Gamma: {
      const double shape = kind_ == DistKind::Chisq ? 0.5 * a : a;
      const double rate = kind_ == DistKind::Chisq ? 0.5 : b;
      if (ll == -kInf || shape == 0) return 0;
      if (lu == -kInf) return kInf;
      const bool use_lower = ll <= lu;
      const double z = qnorm_std(ll, true, true);
      double u0;
      // Wilson-Hilferty when it is positive, else the log-normal approximation of log X.
      const double cc = 1 / (9 * shape);
      const double w = 1 - cc + z * std::sqrt(cc);
      if (shape > 0.5 && w > 0.05) u0 = std::log(shape) + 3 * std::log(w);
      else if (shape < 1 && (ll + sp::lgamma1p(shape)) / shape < 0) u0 = (ll + sp::lgamma1p(shape)) / shape;
      else {
        u0 = sp::digamma(shape) + z * std::sqrt(sp::trigamma(shape));
        // small-x tail: P(a, x) ~ x^a / Gamma(a + 1)
        const double us = (ll + sp::lgamma1p(shape)) / shape;
        if (use_lower && us < u0) u0 = us;
      }
      // Solve on the standard (rate 1) scale.
      const Dist g1 = Dist::gamma(shape, 1.0);
      const double x = solve_quantile(g1, use_lower, use_lower ? ll : lu, Domain::Positive, u0);
      return x / rate;
    }
    case DistKind::Beta: {
      if (ll == -kInf) return 0;
      if (lu == -kInf) return 1;
      if (a == 0 || b == 0 || std::isinf(a) || std::isinf(b)) {
        if (a == 0 && b == 0) return ll < -kLn2 ? 0.0 : (lu < -kLn2 ? 1.0 : 0.5);
        if (a == 0 || std::isinf(b)) return 0;
        if (b == 0 || std::isinf(a)) return 1;
        return 0.5;
      }
      const bool use_lower = ll <= lu;
      const double z = qnorm_std(ll, true, true);
      double u0 = sp::digamma(a) - sp::digamma(b) + z * std::sqrt(sp::trigamma(a) + sp::trigamma(b));
      // tail starts: I_x(a,b) ~ x^a / (a B(a,b)) near 0; 1 - I ~ (1-x)^b / (b B(a,b)) near 1
      const double lB = sp::lbeta(a, b);
      if (use_lower) {
        const double lx = (ll + std::log(a) + lB) / a;
        if (lx < -1 && lx < u0) u0 = lx;
      } else {
        const double ly = (lu + std::log(b) + lB) / b;
        if (ly < -1 && -ly > u0) u0 = -ly;
      }
      return solve_quantile(*this, use_lower, use_lower ? ll : lu, Domain::Unit, u0);
    }
    case DistKind::F: {
      const double m = a, n = b;
      if (ll == -kInf) return 0;
      if (lu == -kInf) return kInf;
      if (std::isinf(m) && std::isinf(n)) return 1;
      if (std::isinf(n)) return Dist::chisq(m).q(p, lower, lg) / m;
      if (std::isinf(m)) return n / Dist::chisq(n).q(p, !lower, lg);
      const bool use_lower = ll <= lu;
      const double z = qnorm_std(ll, true, true);
      const double am = 0.5 * m, an = 0.5 * n;
      double u0 = std::log(n / m) + sp::digamma(am) - sp::digamma(an) + z * std::sqrt(sp::trigamma(am) + sp::trigamma(an));
      return solve_quantile(*this, use_lower, use_lower ? ll : lu, Domain::Positive, u0);
    }
    case DistKind::Binom: {
      const double n = a, pr = b;
      const double sd = std::sqrt(n * pr * (1 - pr));
      const double sk = sd > 0 ? (1 - 2 * pr) / sd : 0;
      return discrete_quantile(*this, p, lower, lg, 0, n, n * pr, sd, sk);
    }
    case DistKind::Pois: {
      if (a == 0) return 0;
      const double sd = std::sqrt(a);
      return discrete_quantile(*this, p, lower, lg, 0, kInf, a, sd, 1 / sd);
    }
    case DistKind::Geom: {
      if (a == 1) return 0;
      if (ll == -kInf) return 0;
      if (lu == -kInf) return kInf;
      return std::max(0.0, std::ceil(lu / std::log1p(-a) - 1 - 1e-12));
    }
    case DistKind::Nbinom: {
      if (a == 0 || b == 1) return 0;
      const double Q = 1 / b, P = (1 - b) / b;
      const double mean = a * P, sd = std::sqrt(a * P * Q);
      return discrete_quantile(*this, p, lower, lg, 0, kInf, mean, sd, (Q + P) / sd);
    }
    case DistKind::Hyper: {
      const double m = a, n = b, k = c, N = m + n;
      const double lo = std::max(0.0, k - n), hi = std::min(k, m);
      const double mean = k * m / N;
      const double var = N > 1 ? k * (m / N) * (n / N) * (N - k) / (N - 1) : 0;
      return discrete_quantile(*this, p, lower, lg, lo, hi, mean, std::sqrt(var), 0);
    }
  }
  return kNaN;
}

// ==============================================================================================================
// Moments

double Dist::mean() const {
  const double a = par_[0], b = par_[1], c = par_[2];
  switch (kind_) {
    case DistKind::Norm: return a;
    case DistKind::Lnorm: return std::exp(a + 0.5 * b * b);
    case DistKind::T: return a > 1 ? 0.0 : kNaN;
    case DistKind::Chisq: return a;
    case DistKind::F: return b > 2 ? b / (b - 2) : kInf;
    case DistKind::Gamma: return a / b;
    case DistKind::Beta: return a / (a + b);
    case DistKind::Exp: return 1 / a;
    case DistKind::Unif: return 0.5 * (a + b);
    case DistKind::Weibull: return b * sp::gamma(1 + 1 / a);
    case DistKind::Cauchy: return kNaN;
    case DistKind::Logis: return a;
    case DistKind::Binom: return a * b;
    case DistKind::Pois: return a;
    case DistKind::Geom: return (1 - a) / a;
    case DistKind::Nbinom: return a * (1 - b) / b;
    case DistKind::Hyper: return c * a / (a + b);
    case DistKind::NcT: {  // E T = ncp sqrt(df/2) Gamma((df-1)/2) / Gamma(df/2)  (Johnson, Kotz & Balakrishnan 31.5)
      if (!(a > 1)) return kNaN;
      if (std::isinf(a)) return b;
      return b * std::sqrt(0.5 * a) * std::exp(sp::lgamma(0.5 * (a - 1)) - sp::lgamma(0.5 * a));
    }
    case DistKind::NcChisq: return a + b;
    case DistKind::NcF: {
      if (std::isinf(b)) return (a + c) / a;
      return b > 2 ? b * (a + c) / (a * (b - 2)) : kInf;
    }
    case DistKind::Tukey: return kNaN;
  }
  return kNaN;
}

double Dist::var() const {
  const double a = par_[0], b = par_[1], c = par_[2];
  switch (kind_) {
    case DistKind::Norm: return b * b;
    case DistKind::Lnorm: return std::expm1(b * b) * std::exp(2 * a + b * b);
    case DistKind::T: return a > 2 ? (std::isinf(a) ? 1.0 : a / (a - 2)) : (a > 1 ? kInf : kNaN);
    case DistKind::Chisq: return 2 * a;
    case DistKind::F: {
      const double m = a, n = b;
      if (n <= 2) return kNaN;
      if (n <= 4) return kInf;
      if (std::isinf(m) && std::isinf(n)) return 0;
      if (std::isinf(n)) return 2 / m;
      if (std::isinf(m)) return 2 * n * n / ((n - 2) * (n - 2) * (n - 4));
      return 2 * n * n * (m + n - 2) / (m * (n - 2) * (n - 2) * (n - 4));
    }
    case DistKind::Gamma: return a / (b * b);
    case DistKind::Beta: { const double s = a + b; return a * b / (s * s * (s + 1)); }
    case DistKind::Exp: return 1 / (a * a);
    case DistKind::Unif: return (b - a) * (b - a) / 12;
    case DistKind::Weibull: {
      const double g1 = sp::gamma(1 + 1 / a), g2 = sp::gamma(1 + 2 / a);
      return b * b * (g2 - g1 * g1);
    }
    case DistKind::Cauchy: return kNaN;
    case DistKind::Logis: return b * b * kPi * kPi / 3;
    case DistKind::Binom: return a * b * (1 - b);
    case DistKind::Pois: return a;
    case DistKind::Geom: return (1 - a) / (a * a);
    case DistKind::Nbinom: return a * (1 - b) / (b * b);
    case DistKind::Hyper: {
      const double N = a + b;
      if (N <= 1) return 0;
      return c * (a / N) * (b / N) * (N - c) / (N - 1);
    }
    case DistKind::NcT: {  // Var T = df (1 + ncp^2) / (df - 2) - (E T)^2
      if (!(a > 1)) return kNaN;
      if (a <= 2) return kInf;
      if (std::isinf(a)) return 1.0;
      const double m = mean();
      return a * (1 + b * b) / (a - 2) - m * m;
    }
    case DistKind::NcChisq: return 2 * (a + 2 * b);
    case DistKind::NcF: {
      const double m = a, n = b;
      if (std::isinf(n)) return 2 * (m + 2 * c) / (m * m);
      if (n <= 2) return kNaN;
      if (n <= 4) return kInf;
      return 2 * (n / m) * (n / m) * ((m + c) * (m + c) + (m + 2 * c) * (n - 2)) / ((n - 2) * (n - 2) * (n - 4));
    }
    case DistKind::Tukey: return kNaN;
  }
  return kNaN;
}

}  // namespace rcli
