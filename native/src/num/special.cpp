// Special functions for the native core (clean-room; written from published mathematics only).
//
// Sources, per algorithm:
// - log Gamma near 1 and 2: the Taylor series ln Gamma(2 + a) = (1 - gamma) a + sum_{k>=2} (-1)^k (zeta(k) - 1) a^k / k
//   (Abramowitz & Stegun 6.1.33 rearranged; DLMF 5.7.3). zeta(k) - 1 by Euler-Maclaurin summation (DLMF 2.10.1).
// - Stirling series with Bernoulli numbers (DLMF 5.11.1) for large arguments; reflection formula (DLMF 5.5.3).
// - digamma / trigamma: recurrence, reflection and asymptotic series (DLMF 5.5, 5.11.2, 5.15.8).
// - stirlerr / bd0 / dbinom_raw / dpois_raw: C. Loader (2000), "Fast and accurate computation of binomial
//   probabilities" (saddle-point deviance form).
// - Regularized incomplete gamma: power series (DLMF 8.11.4 / 8.7.1), Legendre continued fraction (DLMF 8.9.2)
//   evaluated by the modified Lentz method (Lentz 1976; Thompson & Barnett 1986), and the small-a, small-x
//   representation of DiDonato & Morris (1986, ACM TOMS 12) for the upper tail.
// - Regularized incomplete beta: continued fraction DLMF 8.17.22 with the symmetry I_x(a,b) = 1 - I_{1-x}(b,a),
//   prefactor in Loader's deviance form.
// - erf / erfc / normal tails: Maclaurin series of erf (DLMF 7.6.1) and erfc(x) = Q(1/2, x^2) through the Legendre
//   continued fraction; exp(-x^2) evaluated with an exact split of x (a standard technique, cf. Cody 1969).
// - Bessel functions: see the section header further below.
#include "rcli/special.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <complex>

namespace rcli::special {

namespace {

constexpr double kPi = 3.14159265358979323846264338328;
constexpr double kLnSqrt2Pi = 0.918938533204672741780329736406;  // log(sqrt(2 pi))
constexpr double kLn2Pi = 1.83787706640934548356065947281;       // log(2 pi)
constexpr double kLnPi = 1.14472988584940017414342735135;        // log(pi)
constexpr double kLnSqrtPi = 0.572364942924700087071713675677;   // log(sqrt(pi))
constexpr double kEulerGamma = 0.577215664901532860606512090082;
constexpr double kLn2 = 0.693147180559945309417232121458;

// log(1 - exp(x)) for x <= 0, accurate everywhere.
inline double log1mexp(double x) { return x > -kLn2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x)); }

// ---------------------------------------------------------------------------------------------------------------
// zeta(k) - 1, k = 2..kZetaMax, by Euler-Maclaurin: sum_{n=2}^{N-1} n^-k + N^{1-k}/(k-1) + N^-k/2
//   + sum_j B_{2j}/(2j)! k(k+1)...(k+2j-2) N^{-k-2j+1}.
constexpr int kZetaMax = 64;
const std::array<double, kZetaMax + 1>& zeta_m1() {
  static const std::array<double, kZetaMax + 1> table = [] {
    std::array<double, kZetaMax + 1> z{};
    const int N = 20;
    const double B[] = {1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6};
    for (int k = 2; k <= kZetaMax; ++k) {
      double s = 0;
      for (int n = N - 1; n >= 2; --n) s += std::pow(static_cast<double>(n), -k);
      const double Nk = std::pow(static_cast<double>(N), -k);
      double tail = N * Nk / (k - 1) + 0.5 * Nk;
      double fact = 1, rising = 1, Np = Nk / N;
      for (int j = 1; j <= 7; ++j) {
        fact *= (2.0 * j - 1) * (2.0 * j);
        rising = (j == 1) ? k : rising * (k + 2.0 * j - 3) * (k + 2.0 * j - 2);
        tail += B[j - 1] / fact * rising * Np;
        Np /= static_cast<double>(N) * N;
      }
      z[k] = s + tail;
    }
    return z;
  }();
  return table;
}

// lgamma(2 + a), |a| <= 1 (used for |a| <= 0.5).
double lgamma2p(double a) {
  const auto& z = zeta_m1();
  double s = 0, p = -a;
  for (int k = 2; k <= kZetaMax; ++k) {
    p *= -a;  // (-1)^k a^k
    const double t = z[k] * p / k;
    s += t;
    if (std::fabs(t) < 1e-18 * std::fabs(s) + 1e-300) break;
  }
  return a * (1 - kEulerGamma) + s;
}

// psi(2 + a), |a| <= 0.5: derivative of lgamma2p.
double digamma2p(double a) {
  const auto& z = zeta_m1();
  double s = 0, p = a;  // p = (-1)^k a^{k-1}
  for (int k = 2; k <= kZetaMax; ++k) {
    if (k > 2) p *= -a;
    const double t = z[k] * p;
    s += t;
    if (std::fabs(t) < 1e-18 * std::fabs(s) + 1e-300) break;
  }
  return (1 - kEulerGamma) + s;
}

// sin(pi x), cos(pi x) with exact argument reduction.
double sinpi(double x) {
  if (!std::isfinite(x)) return kNaN;
  double r = std::fmod(x, 2.0);  // (-2, 2)
  if (r <= -1) r += 2;
  else if (r > 1) r -= 2;  // (-1, 1]
  if (r == 0 || r == 1) return 0.0;
  if (r == 0.5) return 1.0;
  if (r == -0.5) return -1.0;
  if (r > 0.5) r = 1 - r;
  else if (r < -0.5) r = -1 - r;
  return std::sin(kPi * r);
}
double cospi(double x) {
  if (!std::isfinite(x)) return kNaN;
  double r = std::fmod(std::fabs(x), 2.0);  // [0, 2)
  if (r > 1) r = 2 - r;                     // cos(pi r) symmetric about 1
  return sinpi(0.5 - r);
}

// Stirling series remainder w(x) = lgamma(x) - [(x - 1/2) log x - x + log sqrt(2 pi)], x >= 10.
double stirling_w(double x) {
  const double c[] = {1.0 / 12, -1.0 / 360, 1.0 / 1260, -1.0 / 1680, 1.0 / 1188, -691.0 / 360360, 1.0 / 156,
                      -3617.0 / 122400};
  const double x2 = 1 / (x * x);
  double s = 0;
  for (int k = 7; k >= 0; --k) s = s * x2 + c[k];
  return s / x;
}

// Gamma(x) for 0.5 <= x <= 2.5 via the zeta series, relative accuracy ~1e-16.
double gamma_core(double x) {
  if (x < 1.5) return std::exp(lgamma2p(x - 1) - std::log1p(x - 1));
  return std::exp(lgamma2p(x - 2));
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// Gamma family

double lgamma1p(double a) {
  if (std::isnan(a)) return a;
  if (std::fabs(a) <= 0.5) return lgamma2p(a) - std::log1p(a);
  return lgamma(1 + a);
}

double lgamma(double x) {
  if (std::isnan(x)) return x;
  if (std::isinf(x)) return kInf;
  if (x <= 0 && x == std::floor(x)) return kInf;
  if (x < 0) {
    if (x < -4503599627370496.0) return kInf;  // all such doubles are integers
    const double s = std::fabs(sinpi(x));
    return kLnPi - std::log(s) - lgamma(1 - x);
  }
  if (x < 0.5) return lgamma1p(x) - std::log(x);
  if (x < 1.5) return lgamma1p(x - 1);
  if (x < 2.5) return lgamma2p(x - 2);
  if (x < 15) {
    const double n = std::floor(x - 1.5);
    const double y = x - n;  // [1.5, 2.5)
    double prod = 1;
    for (double t = y; t < x - 0.25; t += 1) prod *= t;
    return std::log(prod) + lgamma2p(y - 2);
  }
  if (x > 1e305) return kInf;
  return (x - 0.5) * std::log(x) - x + kLnSqrt2Pi + stirling_w(x);
}

double gamma(double x) {
  if (std::isnan(x)) return x;
  if (x == kInf) return kInf;
  if (x == -kInf) return kNaN;
  if (x <= 0 && x == std::floor(x)) return kNaN;
  if (x < 0) {
    const double s = sinpi(x);
    if (1 - x > 170) {
      const double l = kLnPi - std::log(std::fabs(s)) - lgamma(1 - x);
      return (s < 0 ? -1.0 : 1.0) * std::exp(l);
    }
    return kPi / (s * gamma(1 - x));
  }
  if (x > 171.61447887182298) return kInf;
  if (x < 0.5) return std::exp(lgamma1p(x)) / x;
  if (x <= 2.5) return gamma_core(x);
  // Shift down to [1.5, 2.5) and multiply back with a compensated (fma-based) product.
  const double n = std::floor(x - 1.5);
  const double y = x - n;
  double hi = 1, lo = 0;
  for (double t = y; t < x - 0.25; t += 1) {
    const double p = hi * t;
    const double e = std::fma(hi, t, -p);
    lo = lo * t + e;
    hi = p;
  }
  return gamma_core(y) * (hi + lo);
}

double beta(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return a + b;
  if (a < 0 || b < 0) return kNaN;
  if (a == 0 || b == 0) return kInf;
  if (std::isinf(a) || std::isinf(b)) return 0;
  if (a + b < 171.61447887182298) {
    const double p = std::min(a, b), q = std::max(a, b);
    return gamma(q) / gamma(p + q) * gamma(p);
  }
  return std::exp(lbeta(a, b));
}

double lbeta(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return a + b;
  const double p = std::min(a, b), q = std::max(a, b);
  if (p < 0) return kNaN;
  if (p == 0) return kInf;
  if (std::isinf(q)) return -kInf;
  if (p >= 10) {
    const double corr = stirling_w(p) + stirling_w(q) - stirling_w(p + q);
    return -0.5 * std::log(q) + kLnSqrt2Pi + corr + (p - 0.5) * std::log(p / (p + q)) + q * std::log1p(-p / (p + q));
  }
  if (q >= 10) {
    const double corr = stirling_w(q) - stirling_w(p + q);
    return lgamma(p) + corr + p - p * std::log(p + q) + (q - 0.5) * std::log1p(-p / (p + q));
  }
  if (p + q < 171.61447887182298) return std::log(gamma(p) * (gamma(q) / gamma(p + q)));
  return lgamma(p) + lgamma(q) - lgamma(p + q);
}

double digamma(double x) {
  if (std::isnan(x)) return x;
  if (x == kInf) return kInf;
  if (x == -kInf) return kNaN;
  if (x <= 0 && x == std::floor(x)) return kNaN;
  if (x < 0) return digamma(1 - x) - kPi * cospi(x) / sinpi(x);
  if (x < 0.5) return digamma2p(x) - 1 / x - 1 / (x + 1);  // psi(x) = psi(x + 2) - 1/(x + 1) - 1/x
  if (x < 1.5) return digamma2p(x - 1) - 1 / x;
  if (x < 2.5) return digamma2p(x - 2);
  if (x < 10) {
    const double n = std::floor(x - 1.5);
    const double y = x - n;
    double s = 0;
    for (double t = y; t < x - 0.25; t += 1) s += 1 / t;
    return digamma2p(y - 2) + s;
  }
  // DLMF 5.11.2
  const double c[] = {1.0 / 12, -1.0 / 120, 1.0 / 252, -1.0 / 240, 1.0 / 132, -691.0 / 32760, 1.0 / 12,
                      -3617.0 / 8160};
  const double x2 = 1 / (x * x);
  double s = 0;
  for (int k = 7; k >= 0; --k) s = s * x2 + c[k];
  return std::log(x) - 0.5 / x - s * x2;
}

double trigamma(double x) {
  if (std::isnan(x)) return x;
  if (x == kInf) return 0;
  if (x == -kInf) return kNaN;
  if (x <= 0 && x == std::floor(x)) return kInf;
  if (x < 0) {
    const double s = sinpi(x);
    return -trigamma(1 - x) + kPi * kPi / (s * s);
  }
  double acc = 0;
  while (x < 10) {
    acc += 1 / (x * x);
    x += 1;
  }
  // DLMF 5.15.8: 1/x + 1/(2x^2) + sum B_2k / x^{2k+1}
  const double c[] = {1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6, -3617.0 / 510};
  const double x2 = 1 / (x * x);
  double s = 0;
  for (int k = 7; k >= 0; --k) s = s * x2 + c[k];
  return acc + 1 / x + 0.5 * x2 + s * x2 / x;
}

namespace {
inline bool is_int(double x) { return std::fabs(x - std::nearbyint(x)) <= 1e-7 * std::max(1.0, std::fabs(x)); }
double lfastchoose(double n, double k) { return -std::log(n + 1) - lbeta(n - k + 1, k + 1); }
// log|choose(n, k)| for n < k - 1 (non-integer n), with the sign of Gamma(n - k + 1).
double lfastchoose2(double n, double k, int* sign) {
  const double r = lgamma(n - k + 1);
  const double v = n - k + 1;
  *sign = (static_cast<long long>(std::floor(v)) % 2 != 0) ? -1 : 1;
  return lgamma(n + 1) - lgamma(k + 1) - r;
}
}  // namespace

double choose(double n, double k) {
  if (std::isnan(n) || std::isnan(k)) return n + k;
  k = std::nearbyint(k);
  const double kSmall = 30;
  if (k < kSmall) {
    if (n - k < k && n >= 0 && is_int(n)) k = std::nearbyint(n - k);
    if (k < 0) return 0;
    if (k == 0) return 1;
    double r = n;
    for (int j = 2; j <= static_cast<int>(k); ++j) r *= (n - j + 1) / j;
    return is_int(n) ? std::nearbyint(r) : r;
  }
  if (n < 0) {
    const double r = choose(-n + k - 1, k);
    return std::fmod(k, 2.0) != 0 ? -r : r;
  }
  if (is_int(n)) {
    n = std::nearbyint(n);
    if (n < k) return 0;
    if (n - k < kSmall) return choose(n, n - k);
    return std::nearbyint(std::exp(lfastchoose(n, k)));
  }
  if (n < k - 1) {
    int s;
    const double r = lfastchoose2(n, k, &s);
    return s * std::exp(r);
  }
  return std::exp(lfastchoose(n, k));
}

double lchoose(double n, double k) {
  if (std::isnan(n) || std::isnan(k)) return n + k;
  k = std::nearbyint(k);
  if (k < 2) {
    if (k < 0) return -kInf;
    if (k == 0) return 0;
    return std::log(std::fabs(n));
  }
  if (n < 0) return lchoose(-n + k - 1, k);
  if (is_int(n)) {
    n = std::nearbyint(n);
    if (n < k) return -kInf;
    if (n - k < 2) return lchoose(n, n - k);
    return lfastchoose(n, k);
  }
  if (n < k - 1) {
    int s;
    return lfastchoose2(n, k, &s);
  }
  return lfastchoose(n, k);
}

double factorial(double x) { return gamma(x + 1); }
double lfactorial(double x) { return lgamma(x + 1); }

double log1pmx(double x) {
  if (std::isnan(x)) return x;
  if (x <= -1) return x == -1 ? -kInf : kNaN;
  if (std::fabs(x) > 0.5) return std::log1p(x) - x;
  // log(1+x) = 2 atanh(r), r = x / (2 + x); log(1+x) - x = -r x + 2 r (r^2/3 + r^4/5 + ...)
  const double r = x / (2 + x), r2 = r * r;
  double s = 0, p = r2;
  for (int j = 1; j < 200; ++j) {
    const double t = p / (2 * j + 1);
    s += t;
    if (t < 1e-18 * s) break;
    p *= r2;
  }
  return -r * x + 2 * r * s;
}

// ---------------------------------------------------------------------------------------------------------------
// Loader's saddle-point kernels

double stirlerr(double n) {
  if (std::isnan(n)) return n;
  if (n <= 0) return n == 0 ? kInf : kNaN;
  if (n >= 15) return stirling_w(n);
  if (n < 0.5) return lgamma1p(n) - (n + 0.5) * std::log(n) + n - kLnSqrt2Pi;
  // s(z) = s(z + 1) + sum_{j>=1} u^{2j} / (2j + 1), u = 1 / (2z + 1).
  double acc = 0, z = n;
  while (z < 15) {
    const double u = 1 / (2 * z + 1), u2 = u * u;
    double p = u2, s = 0;
    for (int j = 1; j < 200; ++j) {
      const double t = p / (2 * j + 1);
      s += t;
      if (t < 1e-18 * s) break;
      p *= u2;
    }
    acc += s;
    z += 1;
  }
  return acc + stirling_w(z);
}

double bd0(double x, double np) {
  if (std::isnan(x) || std::isnan(np)) return x + np;
  if (np == 0) return x == 0 ? 0 : kInf;
  if (x == 0) return np;
  if (std::fabs(x - np) < 0.1 * (x + np)) {
    const double v = (x - np) / (x + np);
    double s = (x - np) * v;
    if (std::fabs(s) < DBL_MIN) return s;
    double ej = 2 * x * v;
    const double v2 = v * v;
    for (int j = 1; j < 1000; ++j) {
      ej *= v2;
      const double s1 = s + ej / (2 * j + 1);
      if (s1 == s) return s1;
      s = s1;
    }
    return s;
  }
  double l = std::log(x / np);
  if (!std::isfinite(l)) l = std::log(x) - std::log(np);
  return x * l + np - x;
}

namespace {
// log[ Gamma(a+b+1) / (Gamma(a+1) Gamma(b+1)) x^a y^b ], a, b > 0, y = 1 - x.
double log_binom_kernel(double a, double b, double x, double y) {
  const double n = a + b;
  const double lc = stirlerr(n) - stirlerr(a) - stirlerr(b) - bd0(a, n * x) - bd0(b, n * y);
  return lc + 0.5 * (std::log(n) - kLn2Pi - std::log(a) - std::log(b));
}
}  // namespace

double dbinom_raw(double x, double n, double p, double q, bool give_log) {
  const double zero = give_log ? -kInf : 0.0, one = give_log ? 0.0 : 1.0;
  if (p == 0) return x == 0 ? one : zero;
  if (q == 0) return x == n ? one : zero;
  double l;
  if (x == 0) {
    if (n == 0) return one;
    l = n * (p < 0.5 ? std::log1p(-p) : std::log(q));
  } else if (x == n) {
    l = n * (q < 0.5 ? std::log1p(-q) : std::log(p));
  } else if (x < 0 || x > n) {
    return zero;
  } else {
    l = log_binom_kernel(x, n - x, p, q);
  }
  return give_log ? l : std::exp(l);
}

double dpois_raw(double x, double lambda, bool give_log) {
  const double zero = give_log ? -kInf : 0.0, one = give_log ? 0.0 : 1.0;
  if (lambda == 0) return x == 0 ? one : zero;
  if (!std::isfinite(lambda) || x < 0) return zero;
  double l;
  if (x <= lambda * DBL_MIN) l = -lambda;
  else if (lambda < x * DBL_MIN) {
    if (!std::isfinite(x)) return zero;
    l = -lambda + x * std::log(lambda) - lgamma(x + 1);
  } else if (x < 1) {
    l = x * std::log(lambda) - lambda - lgamma1p(x);
  } else {
    l = -stirlerr(x) - bd0(x, lambda) - 0.5 * (kLn2Pi + std::log(x));
  }
  return give_log ? l : std::exp(l);
}

// ---------------------------------------------------------------------------------------------------------------
// Incomplete gamma

namespace {
constexpr double kTiny = 1e-300;

// Legendre continued fraction F with Q(a, x) = x^a e^-x / Gamma(a) / F (DLMF 8.9.2), modified Lentz.
double gamma_cf(double a, double x) {
  double f = x + 1 - a;
  if (std::fabs(f) < kTiny) f = kTiny;
  double C = f, D = 0;
  for (long n = 1; n < 100000000L; ++n) {
    const double an = -n * (n - a);
    const double bn = x + 2.0 * n + 1 - a;
    D = bn + an * D;
    if (std::fabs(D) < kTiny) D = kTiny;
    C = bn + an / C;
    if (std::fabs(C) < kTiny) C = kTiny;
    D = 1 / D;
    const double delta = C * D;
    f *= delta;
    if (std::fabs(delta - 1) < 4e-16) break;
  }
  return f;
}

// log P and log Q.
void inc_gamma_log(double a, double x, double* lp, double* lq) {
  if (a < 1 && x < 1) {
    const double lpre = a * std::log(x) - lgamma1p(a);  // log(x^a / Gamma(a+1))
    // P: e^-x x^a / Gamma(a+1) * sum_n x^n / ((a+1)...(a+n))
    double sum = 1, term = 1;
    for (int n = 1; n < 10000; ++n) {
      term *= x / (a + n);
      sum += term;
      if (term < 1e-17 * sum) break;
    }
    *lp = lpre - x + std::log(sum);
    // Q: 1 - x^a/Gamma(a+1) - x^a/Gamma(a+1) * a * sum_{n>=1} (-x)^n / (n! (a+n))
    double S = 0, c = 1;
    for (int n = 1; n < 10000; ++n) {
      c *= -x / n;
      const double t = c / (a + n);
      S += t;
      if (std::fabs(t) < 1e-17 * std::fabs(S)) break;
    }
    const double q = -std::expm1(lpre) - std::exp(lpre) * a * S;
    if (*lp < -kLn2) *lq = log1mexp(*lp);
    else { *lq = std::log(q); *lp = std::log1p(-q); }
    return;
  }
  if (x < a) {
    double sum = 1, term = 1;
    for (long n = 1; n < 100000000L; ++n) {
      term *= x / (a + n);
      sum += term;
      if (term < 1e-17 * sum) break;
    }
    *lp = dpois_raw(a, x, true) + std::log(sum);
    *lq = log1mexp(*lp);
    return;
  }
  const double lpre = (a < 1) ? std::log(a) + a * std::log(x) - x - lgamma1p(a) : std::log(x) + dpois_raw(a - 1, x, true);
  *lq = lpre - std::log(gamma_cf(a, x));
  *lp = log1mexp(*lq);
}
}  // namespace

double gamma_inc(double a, double x, bool lower, bool log_p) {
  if (std::isnan(a) || std::isnan(x)) return a + x;
  if (a < 0) return kNaN;
  double lp, lq;
  if (x <= 0) { lp = -kInf; lq = 0; }
  else if (a == 0 || x == kInf) { lp = 0; lq = -kInf; }
  else if (a == kInf) { lp = -kInf; lq = 0; }
  else inc_gamma_log(a, x, &lp, &lq);
  const double r = lower ? lp : lq;
  return log_p ? r : std::exp(r);
}

double gamma_p(double a, double x) { return gamma_inc(a, x, true, false); }
double gamma_q(double a, double x) { return gamma_inc(a, x, false, false); }

// ---------------------------------------------------------------------------------------------------------------
// Incomplete beta

namespace {
// lgamma(b + a) - lgamma(b) for small a without cancellation: shift b up to >= 10 with log1p terms, then the
// Stirling form (z - 1/2) log1p(a/z) + a log(z + a) - a + w(z + a) - w(z).
double lgamma_shift(double b, double a) {
  double acc = 0, z = b;
  while (z < 10) {
    acc -= std::log1p(a / z);
    z += 1;
  }
  return acc + (z - 0.5) * std::log1p(a / z) + a * std::log(z + a) - a + stirling_w(z + a) - stirling_w(z);
}

// log of x^a y^b / (a B(a, b)).
double beta_log_prefactor(double a, double b, double x, double y) {
  if (a >= 1 && b >= 1) return std::log(b / (a + b)) + log_binom_kernel(a, b, x, y);
  const double lx = x < 0.5 ? std::log(x) : std::log1p(-y);
  const double ly = y < 0.5 ? std::log(y) : std::log1p(-x);
  return a * lx + b * ly - std::log(a) - lbeta(a, b);
}

// Continued fraction DLMF 8.17.22: I_x(a,b) = prefactor / F, F = 1 + d1/(1 + d2/(1 + ...)).
// F is assembled as (1 + d1 + H)/(1 + H), H = d2/(1 + d3/(1 + ...)), with 1 + d1 formed without cancellation;
// this keeps relative accuracy when F is small (x near 1 with b < 1).
double beta_cf(double a, double b, double x, double y) {
  auto dk = [&](long k) {
    const long m = k / 2;
    if (k % 2 == 1) return -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1));
    return m * (b - m) * x / ((a + 2.0 * m - 1) * (a + 2.0 * m));
  };
  double g = 1, C = 1, D = 0;
  for (long k = 3; k < 100000000L; ++k) {
    const double d = dk(k);
    D = 1 + d * D;
    if (std::fabs(D) < kTiny) D = kTiny;
    C = 1 + d / C;
    if (std::fabs(C) < kTiny) C = kTiny;
    D = 1 / D;
    const double delta = C * D;
    g *= delta;
    if (std::fabs(delta - 1) < 4e-16) break;
  }
  const double H = dk(2) / g;
  const double r = (a + b) * x / (a + 1);
  const double one_d1 = r < 0.5 ? 1 - r : (1 - b + (a + b) * y) / (a + 1);
  return (one_d1 + H) / (1 + H);
}

// Large a, small b: with t = e^{-s}, I_x(a,b) = (1/B) int_u^inf e^{-a s} s^{b-1} f(s) ds, u = -log x,
// f(s) = ((1 - e^{-s})/s)^{b-1} = sum_k c_k s^k (radius 2 pi), giving
// I_x(a,b) = R sum_k (b)_k a^{-k} c_k Q(b + k, a u),  1 - I_x(a,b) = R sum_k (b)_k a^{-k} c_k P(b + k, a u),
// R = Gamma(a+b) / (Gamma(a) a^b). (The idea of expanding in incomplete gamma functions is that of
// DiDonato & Morris 1992, ACM TOMS 18, "BGRAT"; this form is derived directly from the integral.)
// The c_k come from the J.C.P. Miller recurrence for powers of a power series.
bool beta_large_a(double a, double b, double x, double y, double* li, double* lic) {
  constexpr int K = 40;
  double g[K + 1], c[K + 1];
  double fact = 1;
  for (int j = 0; j <= K; ++j) {
    fact *= (j + 1);
    g[j] = ((j % 2) ? -1.0 : 1.0) / fact;  // (1 - e^{-s})/s = sum (-1)^j s^j / (j+1)!
  }
  const double beta = b - 1;
  c[0] = 1;
  for (int k = 1; k <= K; ++k) {
    double t = 0;
    for (int j = 1; j <= k; ++j) t += ((beta + 1) * j - k) * g[j] * c[k - j];
    c[k] = t / k;
  }
  if (x < 0.5) return false;  // need u = -log x well inside the radius 2 pi
  const double u = -std::log1p(-y);
  const double au = a * u;
  const double lR = lgamma(b) - lbeta(a, b) - b * std::log(a);
  // sums relative to scale exp(lR)
  double sq = 0, sp = 0, lq0 = 0, lp0 = 0, coef = 1;  // coef = (b)_k a^{-k} c_k
  bool first = true;
  double prev = kInf;
  for (int k = 0; k <= K; ++k) {
    if (k > 0) coef *= (b + k - 1) / a;
    const double ck = coef * c[k];
    if (ck == 0) continue;
    const double lq = gamma_inc(b + k, au, false, true), lpv = gamma_inc(b + k, au, true, true);
    if (first) { lq0 = lq; lp0 = lpv; first = false; }
    const double tq = ck * std::exp(lq - lq0), tp = ck * std::exp(lpv - lp0);
    sq += tq;
    sp += tp;
    const double mag = std::fabs(ck);
    if (mag > prev && k > 3) return false;  // diverging: not in the asymptotic regime
    prev = mag;
    if (std::fabs(tq) <= 1e-17 * std::fabs(sq) && std::fabs(tp) <= 1e-17 * std::fabs(sp) && k > 1) break;
    if (k == K) return false;
  }
  if (!(sq > 0) || !(sp > 0)) return false;
  *li = lR + lq0 + std::log(sq);
  *lic = lR + lp0 + std::log(sp);
  if (*li < *lic) *lic = log1mexp(*li);
  else *li = log1mexp(*lic);
  return true;
}

// log I and log (1 - I), a, b > 0 finite, 0 < x < 1.
void inc_beta_log(double a, double b, double x, double y, double* li, double* lic) {
  // Large-a expansion when the other parameter is small: b^2 << a (terms decay like b^2 / (2 a k)).
  if (a >= 30 && 4 * std::max(b, 1.0) * std::max(b, 1.0) <= a && beta_large_a(a, b, x, y, li, lic)) return;
  if (b >= 30 && 4 * std::max(a, 1.0) * std::max(a, 1.0) <= b && beta_large_a(b, a, y, x, lic, li)) return;
  const bool swap = x > (a + 1) / (a + b + 2);
  if (swap) {
    std::swap(a, b);
    std::swap(x, y);
  }
  double l = beta_log_prefactor(a, b, x, y) - std::log(beta_cf(a, b, x, y));
  double lc = log1mexp(l);
  // Small a with I close to 1: the complement 1 - I from the hypergeometric series (DLMF 8.17.7)
  // I = K [1 + a sum_{n>=1} (1-b)_n x^n / (n! (a+n))], K = x^a Gamma(a+b) / (Gamma(a+1) Gamma(b)).
  if (a < 1 && l > -kLn2 && x <= 0.5 && b * x <= 2) {
    const double lK = a * std::log(x) + lgamma_shift(b, a) - lgamma1p(a);
    double S = 0, c = 1;
    for (int n = 1; n < 2000; ++n) {
      c *= (n - b) * x / n;
      const double t = c / (a + n);
      S += t;
      if (std::fabs(t) <= 1e-17 * std::fabs(S) || c == 0) break;
    }
    const double comp = -std::expm1(lK) - std::exp(lK) * a * S;
    if (comp > 0) {
      lc = std::log(comp);
      l = std::log1p(-comp);
    }
  }
  if (swap) { *lic = l; *li = lc; }
  else { *li = l; *lic = lc; }
}
}  // namespace

double beta_inc_tail(double a, double b, double x, double y, bool lower, bool log_p) {
  if (std::isnan(a) || std::isnan(b) || std::isnan(x) || std::isnan(y)) return a + b + x + y;
  if (a < 0 || b < 0) return kNaN;
  double li, lic;
  if (x <= 0) { li = -kInf; lic = 0; }
  else if (y <= 0) { li = 0; lic = -kInf; }
  else if ((a == 0 || std::isinf(a)) && (b == 0 || std::isinf(b))) {
    if (a == 0 && b == 0) { li = lic = -kLn2; }
    else if (a == 0) { li = 0; lic = -kInf; }
    else if (b == 0) { li = -kInf; lic = 0; }
    else if (x < 0.5) { li = -kInf; lic = 0; }  // both infinite: point mass at 1/2
    else { li = 0; lic = -kInf; }
  } else if (a == 0 || std::isinf(b)) { li = 0; lic = -kInf; }
  else if (b == 0 || std::isinf(a)) { li = -kInf; lic = 0; }
  else inc_beta_log(a, b, x, y, &li, &lic);
  const double r = lower ? li : lic;
  return log_p ? r : std::exp(r);
}

double beta_inc(double a, double b, double x) { return beta_inc_tail(a, b, x, 1 - x, true, false); }
double beta_inc_c(double a, double b, double x) { return beta_inc_tail(a, b, x, 1 - x, false, false); }

// ---------------------------------------------------------------------------------------------------------------
// erf / erfc

namespace {
// erf(x) for |x| < 1 by the Maclaurin series 2/sqrt(pi) sum (-1)^n x^{2n+1} / (n! (2n+1)).
double erf_series(double x) {
  const double x2 = x * x;
  double s = 0, p = x;
  for (int n = 0; n < 60; ++n) {
    const double t = p / (2 * n + 1);
    s += t;
    if (std::fabs(t) < 1e-18 * std::fabs(s)) break;
    p *= -x2 / (n + 1);
  }
  return 1.1283791670955125738961589031215 * s;  // 2 / sqrt(pi)
}
// log(exp(-x^2)) = -x^2 computed with an exact split of x: x = hi + lo, hi with 16 fractional bits.
inline double neg_sq(double x) {
  const double hi = std::ldexp(std::trunc(std::ldexp(x, 16)), -16);
  const double lo = x - hi;
  return -hi * hi - lo * (x + hi);
}
}  // namespace

double erfc(double x) {
  if (std::isnan(x)) return x;
  if (std::fabs(x) < 1) return 1 - erf_series(x);
  if (x < 0) return 2 - erfc(-x);
  if (x > 27.3) return 0;
  // erfc(x) = Q(1/2, x^2) = exp(-x^2) x / sqrt(pi) / F(x^2)
  return std::exp(neg_sq(x)) * x / (1.7724538509055160272981674833411 * gamma_cf(0.5, x * x));
}

double erf(double x) {
  if (std::isnan(x)) return x;
  if (std::fabs(x) < 1) return erf_series(x);
  return x > 0 ? 1 - erfc(x) : erfc(-x) - 1;
}


// ---------------------------------------------------------------------------------------------------------------
// Bessel functions of real order nu >= 0 (negative orders by reflection).
// - I_nu: ascending series (DLMF 10.25.2), all terms positive, rescaled against overflow; Hankel's asymptotic
//   expansion (DLMF 10.40.1) for large x.
// - K_nu: K_nu(x) = (1/2) int_{-inf}^{inf} exp(-x cosh t + nu t) dt (DLMF 10.32.9), trapezoidal rule centred at the
//   saddle point t* = asinh(nu/x); the integrand is entire, so the rule converges geometrically in 1/h
//   (Trefethen & Weideman, SIAM Review 2014). The step is halved until two successive sums agree.
// - J_nu for x^2/4 < nu + 1: ascending series (DLMF 10.2.2), terms decrease monotonically.
// - J_nu, Y_nu otherwise: Steed's method (Barnett, Feng, Steed & Seaton 1974; Thompson & Barnett 1986):
//   CF1 for J'_nu/J_nu with the sign of J_nu from the denominators, downward recurrence to |mu| <= 1/2,
//   complex CF2 for (J' + iY')/(J + iY) and the Wronskian J Y' - J' Y = 2/(pi x) (DLMF 10.5.2); for x < 2,
//   Temme's (1976, J. Comput. Phys. 21) power series for Y_mu, Y_{mu+1}; forward recurrence for Y (DLMF 10.6.1).
//   Hankel's expansion (DLMF 10.17.3/10.17.4) for very large x.
// - Negative order: DLMF 10.4.7/10.4.8 (J, Y; as R, A&S 9.1.2), 10.27.2 (I; A&S 9.6.2), K is even in nu.

namespace {

constexpr double kLn10 = 2.30258509299404568401799145468;

double bessel_i_series(double x, double nu, bool scaled) {
  const double lpre = nu * std::log(0.5 * x) - lgamma(nu + 1) - (scaled ? x : 0.0);
  const double q = 0.25 * x * x;
  double t = 1, s = 1, lsc = 0;
  for (long k = 1; k < 100000000L; ++k) {
    t *= q / (k * (nu + k));
    s += t;
    if (t < 1e-17 * s) break;
    if (s > 1e280) {
      s *= 1e-280;
      t *= 1e-280;
      lsc += 280 * kLn10;
    }
  }
  return std::exp(lpre + lsc + std::log(s));
}

// Hankel expansion of exp(-x) I_nu(x), large x with 4 nu^2 << 8 x.
double bessel_i_scaled_asym(double x, double nu) {
  const double m = 4 * nu * nu;
  double s = 1, t = 1;
  for (int k = 1; k < 200; ++k) {
    const double nt = -t * (m - (2.0 * k - 1) * (2.0 * k - 1)) / (8.0 * k * x);
    if (std::fabs(nt) >= std::fabs(t) && k > 1) break;
    t = nt;
    s += t;
    if (std::fabs(t) < 1e-17 * std::fabs(s)) break;
  }
  return s / std::sqrt(2 * kPi * x);
}

double i_nonneg(double x, double nu, bool scaled) {
  if (x == 0) return nu == 0 ? 1.0 : 0.0;
  if (x > 1000 && nu * nu < 0.25 * x) {
    const double v = bessel_i_scaled_asym(x, nu);
    return scaled ? v : (x > 720 ? kInf : v * std::exp(x));
  }
  if (!scaled && x > 720) {
    // I grows like e^x / sqrt(2 pi x): overflows
    const double v = bessel_i_series(x, nu, true);
    return v == 0 ? 0 : std::exp(std::log(v) + x);
  }
  return bessel_i_series(x, nu, scaled);
}

// sinh(s) - s, accurate for small s.
inline double sinh_minus(double s) {
  if (std::fabs(s) < 1) {
    const double s2 = s * s;
    double t = s * s2 / 6, r = t;
    for (int k = 2; k < 30; ++k) {
      t *= s2 / ((2.0 * k) * (2.0 * k + 1));
      r += t;
      if (std::fabs(t) < 1e-18 * std::fabs(r)) break;
    }
    return r;
  }
  return std::sinh(s) - s;
}

double k_nonneg(double x, double nu, bool scaled) {
  if (x == 0) return kInf;
  const double R = std::hypot(x, nu);
  const double tstar = std::asinh(nu / x);
  const double lpeak = scaled ? (-nu * nu / (R + x) + nu * tstar) : (-R + nu * tstar);
  // f(t* + s) - f(t*) = -R (cosh s - 1) - nu (sinh s - s) = -nu (e^s - 1 - s) - (R - nu)(cosh s - 1),
  // written with non-negative terms only (no cancellation when R ~ nu, i.e. x << nu).
  const double rmn = x * x / (R + nu);
  auto g = [&](double s) {
    const double sh = std::sinh(0.5 * s);
    const double em = (std::fabs(s) < 1) ? 2 * sh * sh + sinh_minus(s) : std::expm1(s) - s;
    return std::exp(-nu * em - rmn * 2 * sh * sh);
  };
  // Sum of g over the lattice offset + j h, j in Z, truncated where terms are negligible.
  auto lattice = [&](double offset, double h) {
    double sum = 0;
    for (int dir = 0; dir < 2; ++dir) {
      for (long j = (dir == 0 ? 0 : -1);; j += (dir == 0 ? 1 : -1)) {
        const double s = offset + j * h;
        const double v = g(s);
        sum += v;
        if (v < 1e-19 * sum && std::fabs(s) > 0.5 * h) break;
        if (std::labs(j) > 100000000L) break;
      }
    }
    return sum;
  };
  double h = std::min(0.5, 1.5 / std::sqrt(R));
  double S = lattice(0, h);
  double T = h * S;
  for (int it = 0; it < 12; ++it) {
    const double S2 = S + lattice(0.5 * h, h);
    const double T2 = 0.5 * h * S2;
    const bool done = std::fabs(T2 - T) <= 1e-12 * T2;
    h *= 0.5;
    S = S2;
    T = T2;
    if (done) break;
  }
  return std::exp(lpeak + std::log(0.5 * T));
}

double j_series(double x, double nu) {
  const double lpre = nu * std::log(0.5 * x) - lgamma(nu + 1);
  const double q = -0.25 * x * x;
  double t = 1, s = 1;
  for (int k = 1; k < 100000; ++k) {
    t *= q / (k * (nu + k));
    s += t;
    if (std::fabs(t) < 1e-17 * std::fabs(s)) break;
  }
  return std::exp(lpre + std::log(s));
}

// Temme's series for Y_mu(x), Y_{mu+1}(x), |mu| <= 1/2, 0 < x < 2.
void temme_y(double x, double mu, double* ymu, double* ymu1) {
  const double x2 = 0.5 * x, ln2x = -std::log(x2);
  const double sigma = mu * ln2x;
  const double rp = std::exp(-lgamma1p(mu));   // 1/Gamma(1+mu)
  const double rm = std::exp(-lgamma1p(-mu));  // 1/Gamma(1-mu)
  const double g2 = 0.5 * (rm + rp);
  // g1 = (rm - rp)/(2 mu) without cancellation: rm - rp = rp expm1(D), D = L(mu) - L(-mu), L = lgamma1p;
  // D/(2 mu) = -atanh(mu)/mu + (1 - gamma) - sum_{k odd >= 3} (zeta(k) - 1) mu^{k-1} / k.
  const auto& z = zeta_m1();
  double dh = (mu == 0 ? -1.0 : -std::atanh(mu) / mu) + (1 - kEulerGamma);
  {
    double p = mu * mu;
    for (int k = 3; k <= kZetaMax; k += 2) {
      const double t = z[k] * p / k;
      dh -= t;
      if (std::fabs(t) < 1e-18) break;
      p *= mu * mu;
    }
  }
  const double D = 2 * mu * dh;
  const double em = (D == 0) ? 1.0 : std::expm1(D) / D;
  const double g1 = rp * em * dh;
  const double fact = (mu == 0) ? 1.0 : mu * kPi / sinpi(mu);
  const double shs = (sigma == 0) ? 1.0 : std::sinh(sigma) / sigma;
  double f = (2 / kPi) * fact * (std::cosh(sigma) * g1 + shs * ln2x * g2);
  const double e = std::exp(sigma);
  double p = e / (kPi * rp);
  double q = 1 / (e * kPi * rm);
  const double hp = 0.5 * kPi * mu;
  const double s3 = (mu == 0) ? 1.0 : std::sin(hp) / hp;
  const double r = kPi * hp * s3 * s3;
  double c = 1;
  const double d = -x2 * x2;
  double sum = f + r * q, sum1 = p;
  for (int k = 1; k < 10000; ++k) {
    f = (k * f + p + q) / (k * k - mu * mu);
    c *= d / k;
    p /= (k - mu);
    q /= (k + mu);
    const double del = c * (f + r * q);
    sum += del;
    const double del1 = c * p - k * del;
    sum1 += del1;
    if (std::fabs(del) < 1e-17 * (1 + std::fabs(sum)) && std::fabs(del1) < 1e-17 * (1 + std::fabs(sum1))) break;
  }
  *ymu = -sum;
  *ymu1 = -sum1 * 2 / x;
}

// Hankel expansion: J and Y for large x.
void jy_hankel(double x, double nu, double* J, double* Y) {
  const double m = 4 * nu * nu;
  double P = 1, Q = 0, t = 1;
  for (int k = 1; k < 200; ++k) {
    const double nt = t * (m - (2.0 * k - 1) * (2.0 * k - 1)) / (8.0 * k * x);
    if (std::fabs(nt) >= std::fabs(t) && k > 2) break;
    t = nt;
    // a_k / x^k with sign (-1)^{floor(k/2)} pattern: P gets even k, Q odd k
    const int r = k % 4;
    if (r == 1) Q += t;
    else if (r == 2) P -= t;
    else if (r == 3) Q -= t;
    else P += t;
    if (std::fabs(t) < 1e-17) break;
  }
  const double cphi = cospi(0.5 * nu + 0.25), sphi = sinpi(0.5 * nu + 0.25);
  const double cx = std::cos(x), sx = std::sin(x);
  const double cchi = cx * cphi + sx * sphi, schi = sx * cphi - cx * sphi;
  const double a = std::sqrt(2 / (kPi * x));
  *J = a * (P * cchi - Q * schi);
  *Y = a * (P * schi + Q * cchi);
}

// J_nu and Y_nu for x > 0, nu >= 0. Y only if want_y.
void jy_nonneg(double x, double nu, double* J, double* Y, bool want_j, bool want_y) {
  if (want_j && !want_y && 0.25 * x * x < nu + 1) {
    *J = j_series(x, nu);
    return;
  }
  if (x > 1000 && 4 * nu * nu < 0.5 * x) {
    jy_hankel(x, nu, J, Y);
    return;
  }
  const double n = std::floor(nu + 0.5);
  const double mu = nu - n;  // [-1/2, 1/2)
  if (x < 2) {
    // J from the series (always applicable here since x^2/4 < 1 <= nu + 1).
    if (want_j) *J = j_series(x, nu);
    if (want_y) {
      double y0, y1;
      temme_y(x, mu, &y0, &y1);
      double m = mu + 1;
      for (double i = 0; i < n; i += 1) {
        const double y2 = (2 * m / x) * y1 - y0;
        y0 = y1;
        y1 = y2;
        m += 1;
        if (!std::isfinite(y0)) break;
      }
      *Y = std::isfinite(y0) ? y0 : -kInf;
    }
    return;
  }
  // CF1 at order nu: f = J'_nu / J_nu = nu/x - 1/(b1 - 1/(b2 - ...)), b_k = 2(nu + k)/x.
  double f = nu / x;
  if (f == 0) f = kTiny;
  double C = f, D = 0;
  int isign = 1;
  for (long k = 1; k < 100000000L; ++k) {
    const double b = 2 * (nu + k) / x;
    D = b - D;
    if (std::fabs(D) < kTiny) D = kTiny;
    C = b - 1 / C;
    if (std::fabs(C) < kTiny) C = kTiny;
    D = 1 / D;
    const double delta = C * D;
    f *= delta;
    if (D < 0) isign = -isign;
    if (std::fabs(delta - 1) < 4e-16) break;
  }
  // Downward recurrence from nu to mu with trial J_nu = isign.
  double jl = isign, jp = f * jl, lsc = 0;
  double m = nu;
  for (double i = 0; i < n; i += 1) {
    const double jm1 = (m / x) * jl + jp;
    jp = ((m - 1) / x) * jm1 - jl;
    jl = jm1;
    m -= 1;
    if (std::fabs(jl) > 1e250) {
      jl *= 1e-250;
      jp *= 1e-250;
      lsc += 250 * kLn10;
    }
  }
  if (jl == 0) jl = kTiny;
  const double fmu = jp / jl;
  // CF2 at mu: p + iq = -1/(2x) + i + (i/x) a1 / (b1 + a2/(b2 + ...)), a_k = (k - 1/2)^2 - mu^2, b_k = 2(x + k i).
  using cd = std::complex<double>;
  const double a1 = 0.25 - mu * mu;
  cd pq(-0.5 / x, 1.0);
  if (a1 != 0) {
    cd W(2 * x, 2.0);
    cd Cc = W, Dc = 0;
    for (long k = 2; k < 100000000L; ++k) {
      const double ak = (k - 0.5) * (k - 0.5) - mu * mu;
      const cd bk(2 * x, 2.0 * k);
      Dc = bk + ak * Dc;
      if (std::abs(Dc) < kTiny) Dc = kTiny;
      Cc = bk + ak / Cc;
      if (std::abs(Cc) < kTiny) Cc = kTiny;
      Dc = 1.0 / Dc;
      const cd delta = Cc * Dc;
      W *= delta;
      if (std::abs(delta - 1.0) < 4e-16) break;
    }
    pq += cd(0, 1.0 / x) * (a1 / W);
  }
  const double p = pq.real(), q = pq.imag();
  const double gam = (p - fmu) / q;
  const double mag = std::sqrt(2 / (kPi * x) / (q + gam * (p - fmu)));
  const double jmu = jl < 0 ? -mag : mag;
  if (want_j) {
    const double l = std::log(mag) - std::log(std::fabs(jl)) - lsc;
    *J = isign * std::exp(l);
  }
  if (want_y) {
    double y0 = gam * jmu;
    const double yp = q * jmu + p * y0;
    double y1 = (mu / x) * y0 - yp;
    double mm = mu + 1;
    for (double i = 0; i < n; i += 1) {
      const double y2 = (2 * mm / x) * y1 - y0;
      y0 = y1;
      y1 = y2;
      mm += 1;
      if (!std::isfinite(y0)) break;
    }
    *Y = std::isfinite(y0) ? y0 : -kInf;
  }
}

}  // namespace

double bessel_j(double x, double nu) {
  if (std::isnan(x) || std::isnan(nu)) return x + nu;
  if (x < 0) return kNaN;
  if (x == kInf) return 0;
  if (nu < 0) {
    const double s = sinpi(nu), c = cospi(nu);
    if (x == 0) return (s == 0) ? (nu == 0 ? 1.0 : 0.0) : kNaN;
    double J = 0, Y = 0;
    jy_nonneg(x, -nu, &J, &Y, true, s != 0);
    return s == 0 ? c * J : c * J + s * Y;  // J_{-v} = cos(v pi) J_v - sin(v pi) Y_v with v = -nu
  }
  if (x == 0) return nu == 0 ? 1.0 : 0.0;
  double J = 0, Y = 0;
  jy_nonneg(x, nu, &J, &Y, true, false);
  return J;
}

double bessel_y(double x, double nu) {
  if (std::isnan(x) || std::isnan(nu)) return x + nu;
  if (x < 0) return kNaN;
  if (x == 0) return -kInf;
  if (x == kInf) return 0;
  if (nu < 0) {
    const double s = sinpi(nu), c = cospi(nu);
    double J = 0, Y = 0;
    jy_nonneg(x, -nu, &J, &Y, s != 0, true);
    return s == 0 ? c * Y : c * Y - s * J;  // Y_{-v} = sin(v pi) J_v + cos(v pi) Y_v with v = -nu
  }
  double J = 0, Y = 0;
  jy_nonneg(x, nu, &J, &Y, false, true);
  return Y;
}

double bessel_i(double x, double nu, bool expon_scaled) {
  if (std::isnan(x) || std::isnan(nu)) return x + nu;
  if (x < 0) return kNaN;
  if (x == kInf) return expon_scaled ? 0.0 : kInf;
  if (nu < 0) {
    const double v = -nu;
    const double s = sinpi(v);
    const double iv = i_nonneg(x, v, expon_scaled);
    if (s == 0) return iv;
    if (x == 0) return kInf;
    const double kv = k_nonneg(x, v, expon_scaled);
    return iv + (2 / kPi) * s * (expon_scaled ? kv * std::exp(-2 * x) : kv);
  }
  return i_nonneg(x, nu, expon_scaled);
}

double bessel_k(double x, double nu, bool expon_scaled) {
  if (std::isnan(x) || std::isnan(nu)) return x + nu;
  if (x < 0) return kNaN;
  if (x == kInf) return 0;
  return k_nonneg(x, std::fabs(nu), expon_scaled);
}

// ---------------------------------------------------------------------------------------------------------------
// Normal distribution kernels (exported to dist.cpp through these non-header functions).

// Phi(z) and 1 - Phi(z), both accurate; log scale if log_p.
void pnorm_both(double z, double* lower, double* upper, bool log_p) {
  if (std::isnan(z)) { *lower = *upper = z; return; }
  if (std::fabs(z) < 1.4142135623730950488) {
    const double e = erf_series(z * 0.70710678118654752440);
    *lower = 0.5 + 0.5 * e;
    *upper = 0.5 - 0.5 * e;
    if (log_p) { *lower = std::log(*lower); *upper = std::log(*upper); }
    return;
  }
  const double az = std::fabs(z);
  double tail, other;  // in the requested scale
  if (az == kInf) {
    tail = log_p ? -kInf : 0.0;
    other = log_p ? 0.0 : 1.0;
  } else if (az > 1e8) {
    const double lt = -0.5 * az * az - std::log(az) - kLnSqrt2Pi + std::log1p(-1 / (az * az));
    tail = log_p ? lt : std::exp(lt);
    other = log_p ? log1mexp(lt) : 1.0;
  } else {
    // Q_N(z) = z phi(z) / (2 F(z^2/2)), phi(z) = exp(-z^2/2) / sqrt(2 pi); exp(-z^2/2) via an exact split of z.
    const double hi = std::ldexp(std::trunc(std::ldexp(az, 16)), -16);
    const double lo = az - hi;
    const double e1 = -0.5 * hi * hi, e2 = -0.5 * lo * (az + hi);
    const double F = gamma_cf(0.5, 0.5 * az * az);
    if (log_p) {
      tail = e1 + e2 + std::log(az / (2 * F)) - kLnSqrt2Pi;
      other = log1mexp(tail);
    } else {
      tail = std::exp(e1) * std::exp(e2) * (az / (2 * F)) * 0.398942280401432677939946059934;  // 1/sqrt(2 pi)
      other = 1 - tail;
    }
  }
  if (z > 0) { *upper = tail; *lower = other; }
  else { *lower = tail; *upper = other; }
}

}  // namespace rcli::special
