// Classical hypothesis tests with R's htest semantics (clean-room; see htest.hpp for the contract).
// Every formula below is implemented from the published literature or from first principles; R was used only as a
// black-box oracle to confirm method strings, branch rules and numeric results.
//
// References:
//   Student (1908); Welch, B. L. (1947) The generalization of Student's problem. Biometrika 34, 28-35.
//   Wilcoxon, F. (1945) Individual comparisons by ranking methods. Biometrics Bulletin 1, 80-83.
//   Mann, H. B. and Whitney, D. R. (1947) Ann. Math. Statist. 18, 50-60 (rank-sum distribution recurrence).
//   Hollander, M. and Wolfe, D. A. (1973) Nonparametric Statistical Methods. Wiley (Walsh averages, Hodges-Lehmann
//     estimators, tie-corrected variances, distribution-free confidence intervals).
//   Bauer, D. F. (1972) Constructing confidence sets using rank statistics. JASA 67, 687-690.
//   Kendall, M. G. (1970) Rank Correlation Methods, 4th ed. Griffin (tau, tie-corrected variance of S, the exact
//     null distribution via the inversion-count recurrence).
//   Best, D. J. and Roberts, D. E. (1975) Algorithm AS 89: The upper tail probabilities of Spearman's rho.
//     Applied Statistics 24, 377-379.
//   Royston, P. (1995) Remark AS R94: A remark on Algorithm AS 181: The W-test for normality.
//     Applied Statistics 44, 547-551.
//   Marsaglia, G., Tsang, W. W. and Wang, J. (2003) Evaluating Kolmogorov's distribution. J. Stat. Softw. 8(18).
//   Birnbaum, Z. W. and Tingey, F. H. (1951) One-sided confidence contours for probability distribution functions.
//     Ann. Math. Statist. 22, 592-596.
//   Schroer, G. and Trenkler, D. (1995) Exact and randomization distributions of Kolmogorov-Smirnov tests for two or
//     three samples. CSDA 20, 185-202 (lattice-path counting, conditional on ties).
//   Feller, W. (1948) Ann. Math. Statist. 19, 177-189 (limiting Kolmogorov / Smirnov distributions).
//   Fisher, R. A. (1935) The logic of inductive inference. JRSS 98, 39-82 (noncentral hypergeometric).
//   Clopper, C. J. and Pearson, E. S. (1934) Biometrika 26, 404-413.
//   Wilson, E. B. (1927) Probable inference ... JASA 22, 209-212; Newcombe, R. G. (1998) Statistics in Medicine 17,
//     857-872 and 873-890 (continuity-corrected score intervals).
//   Yates, F. (1934) Contingency tables involving small numbers and the chi-squared test. JRSS Suppl. 1, 217-235.
//   Brent, R. P. (1973) Algorithms for Minimization without Derivatives, ch. 4 (procedure zero).
#include "rcli/htest.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>

#include "rcli/dist.hpp"
#include "rcli/summary.hpp"

namespace rcli {
namespace {

enum class Alt { Two, Less, Greater };

Alt parse_alt(const std::string& s, const std::string& test) {
  // R matches the argument partially ("t", "l", "g").
  static const char* names[] = {"two.sided", "less", "greater"};
  if (!s.empty())
    for (int i = 0; i < 3; ++i)
      if (std::string(names[i]).compare(0, s.size(), s) == 0) return static_cast<Alt>(i);
  throw Error(test + ": 'alternative' must be one of \"two.sided\", \"less\", \"greater\"");
}

const char* alt_name(Alt a) {
  return a == Alt::Two ? "two.sided" : a == Alt::Less ? "less" : "greater";
}

void check_conf_level(double cl, const std::string& test) {
  if (!(cl > 0.0 && cl < 1.0)) throw Error(test + ": 'conf.level' must be a single number between 0 and 1");
}

Vec drop_nan(const Vec& x) {
  Vec r;
  r.reserve(x.size());
  for (double v : x)
    if (!std::isnan(v)) r.push_back(v);
  return r;
}

Vec drop_nonfinite(const Vec& x) {
  Vec r;
  r.reserve(x.size());
  for (double v : x)
    if (std::isfinite(v)) r.push_back(v);
  return r;
}

const Dist& stdnorm() {
  static const Dist d = Dist::norm(0.0, 1.0);
  return d;
}
double pnorm(double z, bool lower = true) { return stdnorm().p(z, lower); }
double qnorm(double p, bool lower = true) { return stdnorm().q(p, lower); }

// Two-sided / one-sided p-value for a statistic with a continuous symmetric-or-not reference distribution.
double p_cont(const Dist& d, double stat, Alt a) {
  if (a == Alt::Less) return d.p(stat, true);
  if (a == Alt::Greater) return d.p(stat, false);
  return 2.0 * std::min(d.p(stat, true), d.p(stat, false));
}

// Sizes of tie groups in x (only groups of size > 1 matter for corrections, but all are returned).
std::vector<double> tie_sizes(Vec x) {
  std::sort(x.begin(), x.end());
  std::vector<double> t;
  for (size_t i = 0; i < x.size();) {
    size_t j = i + 1;
    while (j < x.size() && x[j] == x[i]) ++j;
    t.push_back(static_cast<double>(j - i));
    i = j;
  }
  return t;
}

bool has_ties(const Vec& x) {
  Vec s = x;
  std::sort(s.begin(), s.end());
  return std::adjacent_find(s.begin(), s.end()) != s.end();
}

double sum_t3_t(const std::vector<double>& t) {
  double s = 0.0;
  for (double v : t) s += v * v * v - v;
  return s;
}

// Round to two significant digits (decimal, ties to even on the exact binary value, like R's signif()).
double signif2(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.1e", v);
  return std::strtod(buf, nullptr);
}

// Median of a vector (sorted in place).
double median_inplace(Vec& v) {
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  if (n == 0) return kNaN;
  return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
}

// ---------------------------------------------------------------------------------------------------------------
// Root finding: Brent's method (Brent 1973, chapter 4: bisection safeguarding secant / inverse quadratic
// interpolation) for f(a), f(b) of opposite signs, absolute tolerance tol (R's uniroot convention: the step
// tolerance is 2 eps |b| + tol / 2). The interpolation step is accepted only when it is smaller than half the
// previous step; this coding reproduces R's uniroot iterates (checked black-box on smooth and step functions).
double brent_zero(const Fn1& f, double a, double b, double fa, double fb, double tol, int maxit = 1000) {
  double c = a, fc = fa;
  for (int it = 0; it < maxit; ++it) {
    const double prev = b - a;
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    const double tol_act = 2.0 * kEps * std::fabs(b) + tol / 2.0;
    double step = (c - b) / 2.0;
    if (std::fabs(step) <= tol_act || fb == 0.0) return b;
    if (std::fabs(prev) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
      const double cb = c - b;
      double p, q;
      if (a == c) {  // secant
        const double t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {  // inverse quadratic interpolation
        q = fa / fc;
        const double t1 = fb / fc, t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0.0) q = -q;
      else p = -p;
      if (p < 0.75 * cb * q - std::fabs(tol_act * q) / 2.0 && p < std::fabs(prev * q / 2.0)) step = p / q;
    }
    if (std::fabs(step) < tol_act) step = step > 0 ? tol_act : -tol_act;
    a = b;
    fa = fb;
    b += step;
    fb = f(b);
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
  return b;
}

// ---------------------------------------------------------------------------------------------------------------
// Exact null distributions (probabilities, accumulated with positive terms only).

// Signed-rank statistic V for n: P(V = s), s = 0..n(n+1)/2. Recurrence on adding rank k with probability 1/2.
Vec signrank_pmf(int n) {
  const size_t M = static_cast<size_t>(n) * (n + 1) / 2;
  Vec p(M + 1, 0.0);
  p[0] = 1.0;
  size_t top = 0;
  for (int k = 1; k <= n; ++k) {
    top += k;
    for (size_t s = top; s >= static_cast<size_t>(k); --s) p[s] = 0.5 * (p[s] + p[s - k]);
    for (size_t s = std::min(top, static_cast<size_t>(k) - 1) + 1; s-- > 0;) p[s] *= 0.5;
  }
  return p;
}

// Rank-sum (Mann-Whitney U) statistic for sizes m, n: P(U = u), u = 0..mn. Recurrence
// p_{a,b}(u) = a/(a+b) p_{a-1,b}(u-b) + b/(a+b) p_{a,b-1}(u) (the largest observation is from x or from y).
Vec wilcox_pmf(int m, int n) {
  const double cost = static_cast<double>(m) * n * (static_cast<double>(m) * n + 1.0);
  if (cost > 4e9) throw Error("wilcox.test: samples too large for the exact distribution");
  // prev[b] = pmf for (a-1, b); cur[b] = pmf for (a, b).
  std::vector<Vec> prev(n + 1), cur(n + 1);
  for (int b = 0; b <= n; ++b) prev[b] = Vec{1.0};  // a = 0
  for (int a = 1; a <= m; ++a) {
    cur[0] = Vec{1.0};
    for (int b = 1; b <= n; ++b) {
      const size_t len = static_cast<size_t>(a) * b + 1;
      Vec v(len, 0.0);
      const double wa = static_cast<double>(a) / (a + b), wb = static_cast<double>(b) / (a + b);
      const Vec& pa = prev[b];   // (a-1, b), support 0..(a-1)b, shifted by b
      const Vec& pb = cur[b - 1];  // (a, b-1), support 0..a(b-1)
      for (size_t u = 0; u < pa.size(); ++u) v[u + b] += wa * pa[u];
      for (size_t u = 0; u < pb.size(); ++u) v[u] += wb * pb[u];
      cur[b] = std::move(v);
    }
    std::swap(prev, cur);
  }
  return prev[n];
}

// Kendall: number of inversions of a random permutation of n, P(I = k), k = 0..n(n-1)/2.
// Recurrence P_n(k) = (1/n) sum_{j=0}^{min(k, n-1)} P_{n-1}(k - j) (insert element n at any of n positions).
Vec kendall_pmf(int n) {
  Vec p{1.0};
  for (int m = 2; m <= n; ++m) {
    const size_t len = p.size() + m - 1;
    Vec q(len, 0.0);
    double window = 0.0;
    for (size_t k = 0; k < len; ++k) {
      if (k < p.size()) window += p[k];
      if (k >= static_cast<size_t>(m) && k - m < p.size()) window -= p[k - m];
      q[k] = window / m;
    }
    p = std::move(q);
  }
  return p;
}

// Lower cumulative sum P(X <= q) and upper P(X >= q) of a pmf on 0..N (q integer-valued).
double cdf_lower(const Vec& pmf, double q) {
  if (q < 0) return 0.0;
  const size_t top = std::min(pmf.size() - 1, static_cast<size_t>(std::floor(q + 1e-7)));
  double s = 0.0;
  for (size_t i = 0; i <= top; ++i) s += pmf[i];
  return std::min(1.0, s);
}
double cdf_upper(const Vec& pmf, double q) {  // P(X >= q)
  const double qc = std::ceil(q - 1e-7);
  if (qc <= 0) return 1.0;
  if (qc > static_cast<double>(pmf.size() - 1)) return 0.0;
  double s = 0.0;
  for (size_t i = static_cast<size_t>(qc); i < pmf.size(); ++i) s += pmf[i];
  return std::min(1.0, s);
}
// Smallest q with P(X <= q) >= p (R's discrete quantile convention).
int pmf_quantile(const Vec& pmf, double p) {
  double s = 0.0;
  for (size_t i = 0; i < pmf.size(); ++i) {
    s += pmf[i];
    if (s >= p * (1.0 - 64.0 * kEps)) return static_cast<int>(i);
  }
  return static_cast<int>(pmf.size() - 1);
}

// ---------------------------------------------------------------------------------------------------------------
// Wilcoxon helpers

double rank_sum_positive(const Vec& d) {  // V = sum of ranks of |d| for positive d (d contains no zeros)
  Vec a(d.size());
  for (size_t i = 0; i < d.size(); ++i) a[i] = std::fabs(d[i]);
  const Vec r = summary::rank(a);
  double v = 0.0;
  for (size_t i = 0; i < d.size(); ++i)
    if (d[i] > 0) v += r[i];
  return v;
}

// Standardised signed-rank statistic of x - mu (zeros dropped), with optional continuity correction.
double signrank_z(const Vec& x, double mu, Alt a, bool correct) {
  Vec d;
  for (double v : x) {
    const double w = v - mu;
    if (w != 0.0) d.push_back(w);
  }
  const double n = static_cast<double>(d.size());
  if (n == 0) return kNaN;
  Vec ad(d.size());
  for (size_t i = 0; i < d.size(); ++i) ad[i] = std::fabs(d[i]);
  const double v = rank_sum_positive(d);
  const double z = v - n * (n + 1) / 4;
  const double sigma = std::sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum_t3_t(tie_sizes(ad)) / 48);
  double corr = 0.0;
  if (correct) corr = a == Alt::Two ? (z > 0 ? 0.5 : z < 0 ? -0.5 : 0.0) : a == Alt::Greater ? 0.5 : -0.5;
  return (z - corr) / sigma;
}

// Standardised rank-sum statistic of (x - mu, y), with optional continuity correction.
double ranksum_z(const Vec& x, const Vec& y, double mu, Alt a, bool correct) {
  const double nx = static_cast<double>(x.size()), ny = static_cast<double>(y.size());
  Vec all;
  all.reserve(x.size() + y.size());
  for (double v : x) all.push_back(v - mu);
  for (double v : y) all.push_back(v);
  const Vec r = summary::rank(all);
  double w = 0.0;
  for (size_t i = 0; i < x.size(); ++i) w += r[i];
  w -= nx * (nx + 1) / 2;
  const double z = w - nx * ny / 2;
  const double N = nx + ny;
  const double sigma = std::sqrt(nx * ny / 12 * ((N + 1) - sum_t3_t(tie_sizes(all)) / (N * (N - 1))));
  double corr = 0.0;
  if (correct) corr = a == Alt::Two ? (z > 0 ? 0.5 : z < 0 ? -0.5 : 0.0) : a == Alt::Greater ? 0.5 : -0.5;
  return (z - corr) / sigma;
}

// Asymptotic interval / estimate by root finding on a decreasing standardised statistic f(mu) - zq over [lo, hi].
double asymptotic_root(const Fn1& f, double zq, double lo, double hi) {
  auto g = [&](double mu) { return f(mu) - zq; };
  const double glo = g(lo), ghi = g(hi);
  if (std::isnan(glo) || std::isnan(ghi)) return kNaN;
  if ((glo > 0 && ghi > 0) || (glo < 0 && ghi < 0)) return kNaN;  // not bracketed: R returns NaN bounds
  // The statistic is a step function of mu; R reports the iterate at which Brent's method with absolute
  // tolerance 1e-4 stops (verified black-box against wilcox.test). Reproducing that iteration exactly keeps the
  // estimate identical to R's even when the root lies on a flat stretch of the step function.
  return brent_zero(g, lo, hi, glo, ghi, 1e-4);
}

void set_ci(HtestResult& r, double lo, double hi, double level) {
  r.has_conf_int = true;
  r.conf_low = lo;
  r.conf_high = hi;
  r.conf_level = level;
}

// ---------------------------------------------------------------------------------------------------------------
// Kolmogorov-Smirnov helpers

// Upper tail of the limiting Kolmogorov distribution, P(K > x), K = sup |Brownian bridge|.
double kolmogorov_upper(double x) {
  if (!(x > 0)) return 1.0;
  if (x < 1.0) {
    // P(K <= x) = sqrt(2 pi)/x sum_{k>=1} exp(-(2k-1)^2 pi^2 / (8 x^2))
    const double pi = 3.14159265358979323846;
    const double w = -pi * pi / (8.0 * x * x);
    double s = 0.0;
    for (int k = 1; k < 100; ++k) {
      const double term = std::exp((2.0 * k - 1) * (2.0 * k - 1) * w);
      s += term;
      if (term < 1e-300 || term < s * 1e-17) break;
    }
    return 1.0 - std::sqrt(2.0 * pi) / x * s;
  }
  // P(K > x) = 2 sum_{k>=1} (-1)^{k-1} exp(-2 k^2 x^2)
  double s = 0.0;
  for (int k = 1; k < 100; ++k) {
    const double term = std::exp(-2.0 * k * k * x * x);
    s += (k % 2 ? 1.0 : -1.0) * term;
    if (term < 1e-300 || term < std::fabs(s) * 1e-17) break;
  }
  return std::min(1.0, 2.0 * s);
}

// Marsaglia, Tsang & Wang (2003): P(D_n < d) for the two-sided one-sample statistic, exact, via the n-th power
// of an m x m matrix (m = 2k - 1, k = floor(n d) + 1), with binary exponent tracking to avoid underflow.
double mtw_cdf(int n, double d) {
  if (d <= 0) return 0.0;
  if (d >= 1) return 1.0;
  const int k = static_cast<int>(n * d) + 1;
  const int m = 2 * k - 1;
  const double h = k - n * d;
  Matrix H(m, m, 0.0);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j) H(i, j) = (i - j + 1 >= 0) ? 1.0 : 0.0;
  for (int i = 0; i < m; ++i) {
    H(i, 0) -= std::pow(h, i + 1);
    H(m - 1, i) -= std::pow(h, m - i);
  }
  H(m - 1, 0) += (2 * h - 1 > 0) ? std::pow(2 * h - 1, m) : 0.0;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      if (i - j + 1 > 0)
        for (int g = 1; g <= i - j + 1; ++g) H(i, j) /= g;
  // R = H^n by repeated squaring; matrices carry a shared binary exponent.
  auto mul = [m](const Matrix& A, const Matrix& B) {
    Matrix C(m, m, 0.0);
    for (int i = 0; i < m; ++i)
      for (int l = 0; l < m; ++l) {
        const double a = A(i, l);
        if (a == 0.0) continue;
        for (int j = 0; j < m; ++j) C(i, j) += a * B(l, j);
      }
    return C;
  };
  auto renorm = [](Matrix& A, long& e) {
    double mx = 0.0;
    for (double v : A.a) mx = std::max(mx, std::fabs(v));
    if (mx == 0.0) return;
    int ex;
    std::frexp(mx, &ex);
    for (double& v : A.a) v = std::ldexp(v, -ex);
    e += ex;
  };
  Matrix R(m, m, 0.0);
  for (int i = 0; i < m; ++i) R(i, i) = 1.0;
  long eR = 0, eP = 0;
  Matrix P = H;
  for (int p = n; p > 0; p >>= 1) {
    if (p & 1) {
      R = mul(R, P);
      eR += eP;
      renorm(R, eR);
    }
    if (p > 1) {
      P = mul(P, P);
      eP *= 2;
      renorm(P, eP);
    }
  }
  // P(D_n < d) = n!/n^n * R[k-1][k-1]
  double s = R(k - 1, k - 1);
  long e = eR;
  for (int i = 1; i <= n; ++i) {
    s = s * i / n;
    int ex;
    s = std::frexp(s, &ex);
    e += ex;
  }
  return std::ldexp(s, static_cast<int>(std::max<long>(std::min<long>(e, 4000), -4000)));
}

// Birnbaum & Tingey (1951): P(D_n^+ >= d) exactly.
double birnbaum_tingey(int n, double d) {
  if (d <= 0) return 1.0;
  if (d >= 1) return 0.0;
  const int jmax = static_cast<int>(std::floor(n * (1 - d) + 1e-12));
  double s = 0.0;
  for (int j = 0; j <= jmax; ++j) {
    const double a = 1.0 - d - static_cast<double>(j) / n;
    const double b = d + static_cast<double>(j) / n;
    double lt = std::lgamma(n + 1.0) - std::lgamma(j + 1.0) - std::lgamma(n - j + 1.0);
    if (n - j > 0) {
      if (a <= 0) continue;
      lt += (n - j) * std::log(a);
    }
    lt += (j - 1) * std::log(b);
    s += std::exp(lt);
  }
  return std::min(1.0, std::max(0.0, d * s));
}

// Exact conditional Smirnov distribution (Schroer & Trenkler 1995): random interleaving of m x's and n y's in the
// pooled order, the statistic only examined at the ends of tie groups of the pooled sample. Returns
// P(D >= K / (m n)), where m n D is |i n - j m| (side 0), i n - j m (side 1, D^+) or j m - i n (side -1, D^-).
// tie_end[k] (k = 1..m+n) tells whether position k closes a tie group.
double smirnov_upper(int m, int n, long K, int side, const std::vector<char>& tie_end) {
  // P[j] = probability of reaching (i, j) without a violation; Q[j] = with a violation already.
  Vec P(n + 1, 0.0), Q(n + 1, 0.0), P2(n + 1), Q2(n + 1);
  P[0] = 1.0;
  const double N = static_cast<double>(m) + n;
  for (int k = 1; k <= m + n; ++k) {
    std::fill(P2.begin(), P2.end(), 0.0);
    std::fill(Q2.begin(), Q2.end(), 0.0);
    const double left = N - (k - 1);
    for (int j = std::max(0, k - m); j <= std::min(k, n); ++j) {
      const int i = k - j;
      double p = 0.0, q = 0.0;
      if (i > 0) {  // previous (i-1, j), step is an x with probability (m - i + 1)/left
        const double w = (m - i + 1) / left;
        p += P[j] * w;
        q += Q[j] * w;
      }
      if (j > 0) {  // previous (i, j-1), step is a y
        const double w = (n - j + 1) / left;
        p += P[j - 1] * w;
        q += Q[j - 1] * w;
      }
      if (tie_end[k]) {
        const long diff = static_cast<long>(i) * n - static_cast<long>(j) * m;
        const bool viol = side == 0 ? std::labs(diff) >= K : side * diff >= K;
        if (viol) {
          q += p;
          p = 0.0;
        }
      }
      P2[j] = p;
      Q2[j] = q;
    }
    std::swap(P, P2);
    std::swap(Q, Q2);
  }
  return std::min(1.0, Q[n]);
}

// ---------------------------------------------------------------------------------------------------------------
// Spearman (AS 89 and exact enumeration)

// Exact distribution of S = sum (i - p_i)^2 over all permutations of n (n <= 9): counts by S/2.
Vec spearman_exact_pmf(int n) {
  std::vector<int> p(n);
  std::iota(p.begin(), p.end(), 0);
  const int smax = n * (n * n - 1) / 3;  // maximum S
  Vec c(smax / 2 + 1, 0.0);
  double total = 0;
  do {
    int s = 0;
    for (int i = 0; i < n; ++i) s += (i - p[i]) * (i - p[i]);
    c[s / 2] += 1.0;
    total += 1.0;
  } while (std::next_permutation(p.begin(), p.end()));
  for (double& v : c) v /= total;
  return c;
}

// Best & Roberts (1975), AS 89: P(S >= js) for n >= 2 by the Edgeworth series (exact enumeration for small n is
// handled by the caller).
double as89_upper_edgeworth(double js, int n) {
  const double c1 = 0.2274, c2 = 0.2531, c3 = 0.1745, c4 = 0.0758, c5 = 0.1033, c6 = 0.3932, c7 = 0.0879,
               c8 = 0.0151, c9 = 0.0072, c10 = 0.0831, c11 = 0.0131, c12 = 4.6e-4;
  const double b = 1.0 / n;
  const double x = (6.0 * (js - 1) * b / (static_cast<double>(n) * n - 1) - 1.0) * std::sqrt(1.0 / b - 1.0);
  double y = x * x;
  const double u =
      x * b *
      (c1 + b * (c2 + c3 * b) +
       y * (-c4 + b * (c5 + c6 * b) - y * b * (c7 + c8 * b - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))));
  y = u / std::exp(y / 2.0);
  const double p = y + pnorm(x, false);
  return std::min(1.0, std::max(0.0, p));
}

// ---------------------------------------------------------------------------------------------------------------
// Fisher 2 x 2 helpers: the conditional distribution of x = table[0][0] given margins is noncentral
// hypergeometric with log weights logdc[i] (central) + i log(psi) over the support lo..hi.
struct NcHyper {
  int lo = 0, hi = 0;
  Vec logdc;
  Vec dens(double psi) const {
    const size_t L = logdc.size();
    Vec d(L);
    if (psi == 0.0 || std::isinf(psi)) {
      std::fill(d.begin(), d.end(), 0.0);
      d[psi == 0.0 ? 0 : L - 1] = 1.0;
      return d;
    }
    const double lp = std::log(psi);
    double mx = -kInf;
    for (size_t i = 0; i < L; ++i) {
      d[i] = logdc[i] + lp * (lo + static_cast<double>(i));
      mx = std::max(mx, d[i]);
    }
    double s = 0.0;
    for (size_t i = 0; i < L; ++i) {
      d[i] = std::exp(d[i] - mx);
      s += d[i];
    }
    for (double& v : d) v /= s;
    return d;
  }
  double mean(double psi) const {
    const Vec d = dens(psi);
    double s = 0.0;
    for (size_t i = 0; i < d.size(); ++i) s += d[i] * (lo + static_cast<double>(i));
    return s;
  }
  double p_le(int x, double psi) const {  // P(X <= x)
    const Vec d = dens(psi);
    double s = 0.0;
    for (int i = lo; i <= x; ++i) s += d[i - lo];
    return std::min(1.0, s);
  }
  double p_ge(int x, double psi) const {  // P(X >= x)
    const Vec d = dens(psi);
    double s = 0.0;
    for (int i = x; i <= hi; ++i) s += d[i - lo];
    return std::min(1.0, s);
  }
};

// Root of g(psi) = target for g monotone in psi. As in R (verified black-box against fisher.test), the search is
// over psi in (0, 1) when the root lies below 1 and over 1/psi in (eps, 1) otherwise, by Brent's method with
// uniroot's default tolerance eps^(1/4); reproducing this iteration keeps estimates and limits identical to R's
// (the limits are therefore only accurate to about 1e-4 relative, like R's).
double solve_psi(const Fn1& g, double target, bool increasing) {
  const double g1 = g(1.0) - target;
  if (g1 == 0.0) return 1.0;
  const bool below = increasing ? g1 > 0 : g1 < 0;  // the root is in (0, 1)
  if (below) {
    auto f = [&](double t) { return g(t) - target; };
    return brent_zero(f, 0.0, 1.0, f(0.0), g1, kEps25);
  }
  auto f = [&](double t) { return g(1.0 / t) - target; };
  return 1.0 / brent_zero(f, kEps, 1.0, f(kEps), g1, kEps25);
}

int as_count(double v, const std::string& test, const char* what) {
  if (!std::isfinite(v) || v < 0) throw Error(test + ": " + what + " must be nonnegative and finite");
  const double r = std::nearbyint(v);
  if (std::fabs(r - v) > 1e-7 * std::max(1.0, std::fabs(v))) throw Error(test + ": " + what + " must be integers");
  if (r > 2147483647.0) throw Error(test + ": " + what + " too large");
  return static_cast<int>(r);
}

}  // namespace

// =================================================================================================================
// t.test

HtestResult t_test(const Vec& x0, const TestOptions& o) {
  const std::string T = "t.test";
  const Alt alt = parse_alt(o.alternative, T);
  check_conf_level(o.conf_level, T);
  if (!std::isfinite(o.mu)) throw Error("t.test: 'mu' must be a single finite number");
  const Vec x = drop_nan(x0);
  const double nx = static_cast<double>(x.size());
  if (nx < 2) throw Error("t.test: not enough 'x' observations");
  for (double v : x)
    if (std::isinf(v)) throw Error("t.test: infinite values are not allowed");  // R fails with an NA error
  const double mx = summary::mean(x), vx = summary::var(x);
  const double df = nx - 1;
  const double stderr_ = std::sqrt(vx / nx);
  if (stderr_ < 10 * kEps * std::fabs(mx)) throw Error("t.test: data are essentially constant");
  const double t = (mx - o.mu) / stderr_;
  const Dist td = Dist::t(df);
  HtestResult r;
  r.method = "One Sample t-test";
  r.alternative = alt_name(alt);
  r.statistic_name = "t";
  r.statistic = t;
  r.parameters = {{"df", df}};
  r.p_value = p_cont(td, t, alt);
  const double cl = o.conf_level;
  if (alt == Alt::Less) set_ci(r, -kInf, o.mu + (t + td.q(cl)) * stderr_, cl);
  else if (alt == Alt::Greater) set_ci(r, o.mu + (t - td.q(cl)) * stderr_, kInf, cl);
  else {
    const double q = td.q(1 - (1 - cl) / 2);
    set_ci(r, o.mu + (t - q) * stderr_, o.mu + (t + q) * stderr_, cl);
  }
  r.estimates = {{"mean of x", mx}};
  r.null_name = "mean";
  r.null_value = o.mu;
  return r;
}

HtestResult t_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  const std::string T = "t.test";
  const Alt alt = parse_alt(o.alternative, T);
  check_conf_level(o.conf_level, T);
  if (!std::isfinite(o.mu)) throw Error("t.test: 'mu' must be a single finite number");
  if (o.paired) {
    if (x0.size() != y0.size()) throw Error("t.test: 'x' and 'y' must have the same length when paired");
    Vec d;
    for (size_t i = 0; i < x0.size(); ++i)
      if (!std::isnan(x0[i]) && !std::isnan(y0[i])) d.push_back(x0[i] - y0[i]);
    if (d.size() < 2) throw Error("t.test: not enough 'x' observations");
    HtestResult r = t_test(d, o);
    r.method = "Paired t-test";
    r.estimates = {{"mean difference", r.estimates[0].second}};
    r.null_name = "mean difference";
    return r;
  }
  const Vec x = drop_nan(x0), y = drop_nan(y0);
  const double nx = static_cast<double>(x.size()), ny = static_cast<double>(y.size());
  if (nx < 1 || (!o.var_equal && nx < 2)) throw Error("t.test: not enough 'x' observations");
  if (ny < 1 || (!o.var_equal && ny < 2)) throw Error("t.test: not enough 'y' observations");
  if (o.var_equal && nx + ny < 3) throw Error("t.test: not enough observations");
  for (const Vec* v : {&x, &y})
    for (double w : *v)
      if (std::isinf(w)) throw Error("t.test: infinite values are not allowed");
  const double mx = summary::mean(x), my = summary::mean(y);
  const double vx = nx > 1 ? summary::var(x) : 0.0, vy = ny > 1 ? summary::var(y) : 0.0;
  double df, se;
  if (o.var_equal) {
    df = nx + ny - 2;
    double v = 0.0;
    if (nx > 1) v += (nx - 1) * vx;
    if (ny > 1) v += (ny - 1) * vy;
    v /= df;
    se = std::sqrt(v * (1 / nx + 1 / ny));
  } else {
    const double sx = vx / nx, sy = vy / ny;
    se = std::sqrt(sx + sy);
    df = (sx + sy) * (sx + sy) / (sx * sx / (nx - 1) + sy * sy / (ny - 1));
  }
  if (se < 10 * kEps * std::max(std::fabs(mx), std::fabs(my))) throw Error("t.test: data are essentially constant");
  const double t = (mx - my - o.mu) / se;
  const Dist td = Dist::t(df);
  HtestResult r;
  r.method = o.var_equal ? " Two Sample t-test" : "Welch Two Sample t-test";
  r.alternative = alt_name(alt);
  r.statistic_name = "t";
  r.statistic = t;
  r.parameters = {{"df", df}};
  r.p_value = p_cont(td, t, alt);
  const double cl = o.conf_level;
  if (alt == Alt::Less) set_ci(r, -kInf, o.mu + (t + td.q(cl)) * se, cl);
  else if (alt == Alt::Greater) set_ci(r, o.mu + (t - td.q(cl)) * se, kInf, cl);
  else {
    const double q = td.q(1 - (1 - cl) / 2);
    set_ci(r, o.mu + (t - q) * se, o.mu + (t + q) * se, cl);
  }
  r.estimates = {{"mean of x", mx}, {"mean of y", my}};
  r.null_name = "difference in means";
  r.null_value = o.mu;
  return r;
}

// =================================================================================================================
// wilcox.test

namespace {

HtestResult signed_rank(const Vec& x, const TestOptions& o, const std::string& null_name) {
  const Alt alt = parse_alt(o.alternative, "wilcox.test");
  const double mu = o.mu;
  // d = x - mu without zeros
  Vec d;
  bool zeroes = false;
  for (double v : x) {
    const double w = v - mu;
    if (w == 0.0) zeroes = true;
    else d.push_back(w);
  }
  const int n = static_cast<int>(d.size());
  const bool exact = o.exact < 0 ? n < 50 : o.exact == 1;
  Vec ad(d.size());
  for (size_t i = 0; i < d.size(); ++i) ad[i] = std::fabs(d[i]);
  const bool ties = has_ties(ad);
  HtestResult r;
  r.alternative = alt_name(alt);
  r.statistic_name = "V";
  r.null_name = null_name;
  r.null_value = mu;
  const double cl = o.conf_level;
  if (n == 0) {
    // Every observation equals mu: V = 0 and the normal approximation is undefined (R returns NaN).
    r.method = std::string("Wilcoxon signed rank test") + (o.correct ? " with continuity correction" : "");
    r.statistic = 0.0;
    r.p_value = kNaN;
    if (o.conf_int) {
      set_ci(r, kNaN, kNaN, 0.0);  // R reports conf.level 0 here
      r.estimates = {{"(pseudo)median", kNaN}};
    }
    return r;
  }
  const double V = rank_sum_positive(d);
  r.statistic = V;
  if (exact && !ties && !zeroes) {
    r.method = "Wilcoxon signed rank exact test";
    const Vec pmf = signrank_pmf(n);
    const double M = static_cast<double>(n) * (n + 1) / 2;
    double p;
    if (alt == Alt::Two) {
      p = V > M / 2 ? cdf_upper(pmf, V) : cdf_lower(pmf, V);
      p = std::min(1.0, 2 * p);
    } else if (alt == Alt::Greater) {
      p = cdf_upper(pmf, V);
    } else {
      p = cdf_lower(pmf, V);
    }
    r.p_value = p;
    if (o.conf_int) {
      // Walsh averages (x_i + x_j)/2, i <= j, of the original data (Hollander & Wolfe; Bauer 1972).
      Vec xs = x;
      Vec w;
      w.reserve(static_cast<size_t>(n) * (n + 1) / 2);
      for (size_t i = 0; i < xs.size(); ++i)
        for (size_t j = i; j < xs.size(); ++j) w.push_back((xs[i] + xs[j]) / 2);
      std::sort(w.begin(), w.end());
      const long Mi = static_cast<long>(w.size());
      double achieved = cl;
      auto lower_q = [&](double a) {
        int q = pmf_quantile(pmf, a);
        if (q == 0) q = 1;
        return q;
      };
      if (alt == Alt::Two) {
        const int q = lower_q((1 - cl) / 2);
        achieved = 1 - 2 * cdf_lower(pmf, q - 1);
        set_ci(r, w[q - 1], w[Mi - q], cl);
      } else if (alt == Alt::Greater) {
        const int q = lower_q(1 - cl);
        achieved = 1 - cdf_lower(pmf, q - 1);
        set_ci(r, w[q - 1], kInf, cl);
      } else {
        const int q = lower_q(1 - cl);
        achieved = 1 - cdf_lower(pmf, q - 1);
        set_ci(r, -kInf, w[Mi - q], cl);
      }
      // R substitutes the achievable level (with a warning), reported as 1 - alpha with alpha rounded to two
      // significant digits in the signed-rank case (oracle: n = 4, conf.level 0.99 -> 0.88, n = 6 -> 0.984).
      if (achieved < cl) r.conf_level = 1.0 - signif2(1.0 - achieved);
      r.estimates = {{"(pseudo)median", median_inplace(w)}};
    }
    return r;
  }
  r.method = std::string("Wilcoxon signed rank test") + (o.correct ? " with continuity correction" : "");
  const double z = signrank_z(x, mu, alt, o.correct);
  r.p_value = p_cont(stdnorm(), z, alt);
  if (o.conf_int) {
    // Observations equal to mu are dropped before the interval is inverted, as for the test itself.
    Vec xs;
    for (double v : d) xs.push_back(v + mu);
    const double lo = *std::min_element(xs.begin(), xs.end()), hi = *std::max_element(xs.begin(), xs.end());
    auto fc = [&](double m) { return signrank_z(xs, m, alt, o.correct); };
    auto f0 = [&](double m) { return signrank_z(xs, m, alt, false); };
    if (alt == Alt::Two) {
      const double zq = qnorm((1 - cl) / 2, false);
      set_ci(r, asymptotic_root(fc, zq, lo, hi), asymptotic_root(fc, -zq, lo, hi), cl);
    } else if (alt == Alt::Greater) {
      set_ci(r, asymptotic_root(fc, qnorm(1 - cl, false), lo, hi), kInf, cl);
    } else {
      set_ci(r, -kInf, asymptotic_root(fc, qnorm(1 - cl, true), lo, hi), cl);
    }
    r.estimates = {{"(pseudo)median", asymptotic_root(f0, 0.0, lo, hi)}};
  }
  return r;
}

}  // namespace

HtestResult wilcox_test(const Vec& x0, const TestOptions& o) {
  parse_alt(o.alternative, "wilcox.test");
  if (!std::isfinite(o.mu)) throw Error("wilcox.test: 'mu' must be a single finite number");
  if (o.conf_int) check_conf_level(o.conf_level, "wilcox.test");
  const Vec x = drop_nan(x0);  // R keeps infinite values (they take the extreme ranks)
  if (x.empty()) throw Error("wilcox.test: not enough (non-missing) 'x' observations");
  return signed_rank(x, o, "location");
}

HtestResult wilcox_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  const std::string T = "wilcox.test";
  const Alt alt = parse_alt(o.alternative, T);
  if (!std::isfinite(o.mu)) throw Error("wilcox.test: 'mu' must be a single finite number");
  if (o.conf_int) check_conf_level(o.conf_level, T);
  if (o.paired) {
    if (x0.size() != y0.size()) throw Error("wilcox.test: 'x' and 'y' must have the same length");
    Vec d;
    for (size_t i = 0; i < x0.size(); ++i) {
      const double v = x0[i] - y0[i];
      if (!std::isnan(v)) d.push_back(v);
    }
    if (d.empty()) throw Error("wilcox.test: not enough (non-missing) 'x' observations");
    return signed_rank(d, o, "location shift");
  }
  const Vec x = drop_nan(x0), y = drop_nan(y0);
  if (x.empty()) throw Error("wilcox.test: not enough (non-missing) 'x' observations");
  if (y.empty()) throw Error("wilcox.test: not enough 'y' observations");
  const int nx = static_cast<int>(x.size()), ny = static_cast<int>(y.size());
  Vec all;
  for (double v : x) all.push_back(v - o.mu);
  for (double v : y) all.push_back(v);
  const Vec rk = summary::rank(all);
  double W = 0.0;
  for (int i = 0; i < nx; ++i) W += rk[i];
  W -= static_cast<double>(nx) * (nx + 1) / 2;
  const bool exact = o.exact < 0 ? (nx < 50 && ny < 50) : o.exact == 1;
  const bool ties = has_ties(all);
  HtestResult r;
  r.alternative = alt_name(alt);
  r.statistic_name = "W";
  r.statistic = W;
  r.null_name = "location shift";
  r.null_value = o.mu;
  const double cl = o.conf_level;
  if (exact && !ties) {
    r.method = "Wilcoxon rank sum exact test";
    const Vec pmf = wilcox_pmf(nx, ny);
    double p;
    if (alt == Alt::Two) {
      p = W > static_cast<double>(nx) * ny / 2 ? cdf_upper(pmf, W) : cdf_lower(pmf, W);
      p = std::min(1.0, 2 * p);
    } else if (alt == Alt::Greater) {
      p = cdf_upper(pmf, W);
    } else {
      p = cdf_lower(pmf, W);
    }
    r.p_value = p;
    if (o.conf_int) {
      Vec dif;
      dif.reserve(static_cast<size_t>(nx) * ny);
      for (double a : x)
        for (double b : y) dif.push_back(a - b);
      std::sort(dif.begin(), dif.end());
      const long Mi = static_cast<long>(dif.size());
      double achieved = cl;
      auto lower_q = [&](double a) {
        int q = pmf_quantile(pmf, a);
        if (q == 0) q = 1;
        return q;
      };
      if (alt == Alt::Two) {
        const int q = lower_q((1 - cl) / 2);
        achieved = 1 - 2 * cdf_lower(pmf, q - 1);
        set_ci(r, dif[q - 1], dif[Mi - q], cl);
      } else if (alt == Alt::Greater) {
        const int q = lower_q(1 - cl);
        achieved = 1 - cdf_lower(pmf, q - 1);
        set_ci(r, dif[q - 1], kInf, cl);
      } else {
        const int q = lower_q(1 - cl);
        achieved = 1 - cdf_lower(pmf, q - 1);
        set_ci(r, -kInf, dif[Mi - q], cl);
      }
      if (achieved < cl) r.conf_level = achieved;
      r.estimates = {{"difference in location", median_inplace(dif)}};
    }
    return r;
  }
  r.method = std::string("Wilcoxon rank sum test") + (o.correct ? " with continuity correction" : "");
  const double z = ranksum_z(x, y, o.mu, alt, o.correct);
  r.p_value = p_cont(stdnorm(), z, alt);
  if (o.conf_int) {
    const double lo = *std::min_element(x.begin(), x.end()) - *std::max_element(y.begin(), y.end());
    const double hi = *std::max_element(x.begin(), x.end()) - *std::min_element(y.begin(), y.end());
    auto fc = [&](double m) { return ranksum_z(x, y, m, alt, o.correct); };
    auto f0 = [&](double m) { return ranksum_z(x, y, m, alt, false); };
    if (alt == Alt::Two) {
      const double zq = qnorm((1 - cl) / 2, false);
      set_ci(r, asymptotic_root(fc, zq, lo, hi), asymptotic_root(fc, -zq, lo, hi), cl);
    } else if (alt == Alt::Greater) {
      set_ci(r, asymptotic_root(fc, qnorm(1 - cl, false), lo, hi), kInf, cl);
    } else {
      set_ci(r, -kInf, asymptotic_root(fc, qnorm(1 - cl, true), lo, hi), cl);
    }
    r.estimates = {{"difference in location", asymptotic_root(f0, 0.0, lo, hi)}};
  }
  return r;
}

// =================================================================================================================
// var.test

HtestResult var_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  const std::string T = "var.test";
  const Alt alt = parse_alt(o.alternative, T);
  check_conf_level(o.conf_level, T);
  if (!(o.ratio > 0 && std::isfinite(o.ratio))) throw Error("var.test: 'ratio' must be a single positive number");
  const Vec x = drop_nonfinite(x0), y = drop_nonfinite(y0);
  if (x.size() < 2) throw Error("var.test: not enough 'x' observations");
  if (y.size() < 2) throw Error("var.test: not enough 'y' observations");
  const double dfx = x.size() - 1.0, dfy = y.size() - 1.0;
  const double est = summary::var(x) / summary::var(y);
  const double F = est / o.ratio;
  const Dist fd = Dist::f(dfx, dfy);
  HtestResult r;
  r.method = "F test to compare two variances";
  r.alternative = alt_name(alt);
  r.statistic_name = "F";
  r.statistic = F;
  r.parameters = {{"num df", dfx}, {"denom df", dfy}};
  const double pv = fd.p(F);
  r.p_value = alt == Alt::Less ? pv : alt == Alt::Greater ? fd.p(F, false) : 2 * std::min(pv, fd.p(F, false));
  const double cl = o.conf_level;
  if (alt == Alt::Two) {
    const double b = (1 - cl) / 2;
    set_ci(r, est / fd.q(1 - b), est / fd.q(b), cl);
  } else if (alt == Alt::Greater) {
    set_ci(r, est / fd.q(cl), kInf, cl);
  } else {
    set_ci(r, 0.0, est / fd.q(1 - cl), cl);
  }
  r.estimates = {{"ratio of variances", est}};
  r.null_name = "ratio of variances";
  r.null_value = o.ratio;
  return r;
}

// =================================================================================================================
// cor.test

HtestResult cor_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  const std::string T = "cor.test";
  const Alt alt = parse_alt(o.alternative, T);
  if (x0.size() != y0.size()) throw Error("cor.test: 'x' and 'y' must have the same length");
  Vec x, y;
  for (size_t i = 0; i < x0.size(); ++i)
    if (!std::isnan(x0[i]) && !std::isnan(y0[i])) {
      x.push_back(x0[i]);
      y.push_back(y0[i]);
    }
  const int n = static_cast<int>(x.size());
  HtestResult r;
  r.alternative = alt_name(alt);
  r.null_value = 0.0;
  if (o.method == "pearson") {
    check_conf_level(o.conf_level, T);
    if (n < 3) throw Error("cor.test: not enough finite observations");
    const double est = summary::cor(x, y, "pearson");
    const double df = n - 2.0;
    const double t = std::sqrt(df) * est / std::sqrt(1 - est * est);
    r.method = "Pearson's product-moment correlation";
    r.statistic_name = "t";
    r.statistic = t;
    r.parameters = {{"df", df}};
    r.p_value = p_cont(Dist::t(df), t, alt);
    if (n > 3) {
      const double z = std::atanh(est), sigma = 1 / std::sqrt(n - 3.0), cl = o.conf_level;
      if (alt == Alt::Less) set_ci(r, -1.0, std::tanh(z + sigma * qnorm(cl)), cl);
      else if (alt == Alt::Greater) set_ci(r, std::tanh(z - sigma * qnorm(cl)), 1.0, cl);
      else {
        const double q = qnorm((1 + cl) / 2);
        set_ci(r, std::tanh(z - sigma * q), std::tanh(z + sigma * q), cl);
      }
    }
    r.estimates = {{"cor", est}};
    r.null_name = "correlation";
    return r;
  }
  if (n < 2) throw Error("cor.test: not enough finite observations");
  const bool ties = has_ties(x) || has_ties(y);
  if (o.method == "kendall") {
    r.method = "Kendall's rank correlation tau";
    r.null_name = "tau";
    // S = concordant - discordant pairs; tau-b = S / sqrt((T0 - Tx)(T0 - Ty)).
    double S = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        const double a = (x[j] > x[i]) - (x[j] < x[i]), b = (y[j] > y[i]) - (y[j] < y[i]);
        S += a * b;
      }
    const double dn = n, T0 = dn * (dn - 1) / 2;
    const std::vector<double> tx = tie_sizes(x), ty = tie_sizes(y);
    double sx = 0, sy = 0;
    for (double t : tx) sx += t * (t - 1) / 2;
    for (double t : ty) sy += t * (t - 1) / 2;
    const double tau = S / std::sqrt((T0 - sx) * (T0 - sy));
    r.estimates = {{"tau", tau}};
    const bool exact = o.exact < 0 ? n < 50 : o.exact == 1;
    if (exact && !ties) {
      const double q = std::round((S + T0) / 2);  // number of concordant pairs
      r.statistic_name = "T";
      r.statistic = q;
      const Vec pmf = kendall_pmf(n);
      double p;
      if (alt == Alt::Two) {
        p = q > T0 / 2 ? cdf_upper(pmf, q) : cdf_lower(pmf, q);
        p = std::min(1.0, 2 * p);
      } else if (alt == Alt::Greater) {
        p = cdf_upper(pmf, q);
      } else {
        p = cdf_lower(pmf, q);
      }
      r.p_value = p;
    } else {
      double vt = 0, vu = 0, v1a = 0, v1b = 0, v2a = 0, v2b = 0;
      for (double t : tx) {
        vt += t * (t - 1) * (2 * t + 5);
        v1a += t * (t - 1);
        v2a += t * (t - 1) * (t - 2);
      }
      for (double t : ty) {
        vu += t * (t - 1) * (2 * t + 5);
        v1b += t * (t - 1);
        v2b += t * (t - 1) * (t - 2);
      }
      const double v0 = dn * (dn - 1) * (2 * dn + 5);
      const double var = (v0 - vt - vu) / 18 + v1a * v1b / (2 * dn * (dn - 1)) +
                         (n > 2 ? v2a * v2b / (9 * dn * (dn - 1) * (dn - 2)) : 0.0);
      const double z = S / std::sqrt(var);
      r.statistic_name = "z";
      r.statistic = z;
      r.p_value = p_cont(stdnorm(), z, alt);
    }
    return r;
  }
  if (o.method == "spearman") {
    r.method = "Spearman's rank correlation rho";
    r.null_name = "rho";
    const double rho = summary::cor(summary::rank(x), summary::rank(y), "pearson");
    const double dn = n;
    const double S = (dn * dn * dn - dn) * (1 - rho) / 6;
    r.statistic_name = "S";
    r.statistic = S;
    r.estimates = {{"rho", rho}};
    const bool exact = o.exact < 0 ? n < 1290 : o.exact == 1;
    if (exact && !ties && n < 1290) {
      // Exact for n <= 9 by enumeration; AS 89 Edgeworth series otherwise. S takes even values only.
      double lower, upper;  // P(S <= s), P(S >= s)
      const double s = std::round(S);
      if (n <= 9) {
        const Vec pmf = spearman_exact_pmf(n);
        lower = cdf_lower(pmf, s / 2);
        upper = cdf_upper(pmf, s / 2);
      } else {
        upper = as89_upper_edgeworth(s, n);
        lower = 1 - as89_upper_edgeworth(s + 2, n);
      }
      r.p_value = alt == Alt::Greater ? lower : alt == Alt::Less ? upper : std::min(1.0, 2 * std::min(lower, upper));
    } else {
      const double t = rho / std::sqrt((1 - rho * rho) / (dn - 2));
      r.p_value = p_cont(Dist::t(dn - 2), t, alt);
    }
    return r;
  }
  throw Error("cor.test: 'method' must be one of \"pearson\", \"kendall\", \"spearman\"");
}

// =================================================================================================================
// shapiro.test (Royston 1995, AS R94)

HtestResult shapiro_test(const Vec& x0) {
  Vec x = drop_nan(x0);
  const int n = static_cast<int>(x.size());
  if (n < 3 || n > 5000) throw Error("shapiro.test: sample size must be between 3 and 5000");
  std::sort(x.begin(), x.end());
  const double range = x.back() - x.front();
  if (!std::isfinite(range)) throw Error("shapiro.test: infinite values are not allowed");
  if (range < 1e-10) throw Error("shapiro.test: all 'x' values are identical");
  const int nn2 = n / 2;
  // a[i], i = 0..nn2-1: coefficient for the i-th largest (positive); the smallest gets -a[i].
  Vec a(nn2);
  auto poly = [](const double* c, int nc, double v) {
    double s = c[nc - 1];
    for (int i = nc - 2; i >= 0; --i) s = s * v + c[i];
    return s;
  };
  if (n == 3) {
    a[0] = std::sqrt(0.5);
  } else {
    // m_i = Phi^{-1}((i - 3/8) / (n + 1/4)); by symmetry use the upper half, m[k] for the (k+1)-th largest.
    Vec m(nn2);
    double summ2 = 0.0;
    for (int k = 0; k < nn2; ++k) {
      m[k] = qnorm((n - k - 0.375) / (n + 0.25));
      summ2 += 2 * m[k] * m[k];
    }
    const double ssumm2 = std::sqrt(summ2);
    const double u = 1.0 / std::sqrt(static_cast<double>(n));
    static const double c1[6] = {0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056};
    static const double c2[6] = {0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633};
    a[0] = m[0] / ssumm2 + poly(c1, 6, u);
    int first;
    double phi;
    if (n > 5) {
      a[1] = m[1] / ssumm2 + poly(c2, 6, u);
      phi = (summ2 - 2 * m[0] * m[0] - 2 * m[1] * m[1]) / (1 - 2 * a[0] * a[0] - 2 * a[1] * a[1]);
      first = 2;
    } else {
      phi = (summ2 - 2 * m[0] * m[0]) / (1 - 2 * a[0] * a[0]);
      first = 1;
    }
    const double sp = std::sqrt(phi);
    for (int k = first; k < nn2; ++k) a[k] = m[k] / sp;
  }
  // W = (sum a_i x_(i))^2 / sum (x - mean)^2, data centred and scaled by the range for stability.
  const double mean = summary::mean(x);
  double num = 0.0, ss = 0.0;
  for (int k = 0; k < nn2; ++k) num += a[k] * ((x[n - 1 - k] - x[k]) / range);
  for (double v : x) {
    const double d = (v - mean) / range;
    ss += d * d;
  }
  double W = num * num / ss;
  if (W > 1.0) W = 1.0;
  double p;
  if (n == 3) {
    const double pi = 3.14159265358979323846;
    p = 6.0 / pi * (std::asin(std::sqrt(W)) - std::asin(std::sqrt(0.75)));
    p = std::max(0.0, std::min(1.0, p));
  } else {
    const double w1 = std::log1p(-W);
    const double dn = n;
    if (n <= 11) {
      const double gamma = -2.273 + 0.459 * dn;
      if (w1 >= gamma) {
        p = 0.0;
      } else {
        const double y = -std::log(gamma - w1);
        static const double c3[4] = {0.5440, -0.39978, 0.025054, -6.714e-4};
        static const double c4[4] = {1.3822, -0.77857, 0.062767, -0.0020322};
        const double mm = poly(c3, 4, dn), s = std::exp(poly(c4, 4, dn));
        p = pnorm((y - mm) / s, false);
      }
    } else {
      const double ln = std::log(dn);
      static const double c5[4] = {-1.5861, -0.31082, -0.083751, 0.0038915};
      static const double c6[3] = {-0.4803, -0.082676, 0.0030302};
      const double mm = poly(c5, 4, ln), s = std::exp(poly(c6, 3, ln));
      p = pnorm((w1 - mm) / s, false);
    }
  }
  HtestResult r;
  r.method = "Shapiro-Wilk normality test";
  r.alternative = "";
  r.statistic_name = "W";
  r.statistic = W;
  r.p_value = p;
  return r;
}

// =================================================================================================================
// ks.test

HtestResult ks_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  const std::string T = "ks.test";
  const Alt alt = parse_alt(o.alternative, T);
  const Vec x = drop_nan(x0), y = drop_nan(y0);
  if (x.empty()) throw Error("ks.test: not enough 'x' data");
  if (y.empty()) throw Error("ks.test: not enough 'y' data");
  const int m = static_cast<int>(x.size()), n = static_cast<int>(y.size());
  // Pooled order: 1 for x, 0 for y; ties grouped.
  std::vector<std::pair<double, int>> w;
  for (double v : x) w.push_back({v, 1});
  for (double v : y) w.push_back({v, 0});
  std::sort(w.begin(), w.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  std::vector<char> tie_end(m + n + 1, 0);
  long i = 0, j = 0, dmax_pos = 0, dmax_neg = 0;  // in units of 1/(m n): i n - j m
  for (int k = 0; k < m + n; ++k) {
    if (w[k].second) ++i;
    else ++j;
    if (k + 1 == m + n || w[k + 1].first != w[k].first) {
      tie_end[k + 1] = 1;
      const long diff = i * n - j * static_cast<long>(m);
      dmax_pos = std::max(dmax_pos, diff);
      dmax_neg = std::max(dmax_neg, -diff);
    }
  }
  long K = alt == Alt::Two ? std::max(dmax_pos, dmax_neg) : alt == Alt::Greater ? dmax_pos : dmax_neg;
  const double mn = static_cast<double>(m) * n;
  const double D = K / mn;
  const bool exact = o.exact < 0 ? mn < 10000 : o.exact == 1;
  HtestResult r;
  r.statistic_name = alt == Alt::Two ? "D" : alt == Alt::Greater ? "D^+" : "D^-";
  r.statistic = D;
  r.alternative = alt == Alt::Two ? "two-sided"
                  : alt == Alt::Less ? "the CDF of x lies below that of y"
                                     : "the CDF of x lies above that of y";
  if (exact) {
    r.method = "Exact two-sample Kolmogorov-Smirnov test";
    if (K == 0) {
      r.p_value = 1.0;
    } else {
      r.p_value = smirnov_upper(m, n, K, alt == Alt::Two ? 0 : alt == Alt::Greater ? 1 : -1, tie_end);
    }
  } else {
    r.method = "Asymptotic two-sample Kolmogorov-Smirnov test";
    const double z = std::sqrt(mn / (m + n)) * D;
    // R reports the complement of the limiting distribution function, 1 - P(K <= z), so tails below the
    // double-precision resolution of 1 - p are rounded the same way (e.g. 2.2e-16 or 0).
    r.p_value = 1.0 - (1.0 - (alt == Alt::Two ? kolmogorov_upper(z) : std::exp(-2 * z * z)));
  }
  return r;
}

HtestResult ks_test(const Vec& x0, const Dist& d, const TestOptions& o) {
  const std::string T = "ks.test";
  const Alt alt = parse_alt(o.alternative, T);
  if (d.discrete()) throw Error("ks.test: only continuous distributions are supported");
  Vec x = drop_nan(x0);
  const int n = static_cast<int>(x.size());
  if (n < 1) throw Error("ks.test: not enough 'x' data");
  std::sort(x.begin(), x.end());
  double dplus = -kInf, dminus = -kInf;
  for (int i = 0; i < n; ++i) {
    const double F = d.p(x[i]);
    dplus = std::max(dplus, (i + 1.0) / n - F);
    dminus = std::max(dminus, F - static_cast<double>(i) / n);
  }
  const double D = alt == Alt::Two ? std::max(dplus, dminus) : alt == Alt::Greater ? dplus : dminus;
  const bool ties = has_ties(x);
  const bool exact = (o.exact < 0 ? n < 100 : o.exact == 1) && !ties;
  HtestResult r;
  r.statistic_name = alt == Alt::Two ? "D" : alt == Alt::Greater ? "D^+" : "D^-";
  r.statistic = D;
  r.alternative = alt == Alt::Two ? "two-sided"
                  : alt == Alt::Less ? "the CDF of x lies below the null hypothesis"
                                     : "the CDF of x lies above the null hypothesis";
  if (exact) {
    r.method = "Exact one-sample Kolmogorov-Smirnov test";
    r.p_value = alt == Alt::Two ? std::min(1.0, std::max(0.0, 1.0 - mtw_cdf(n, D))) : birnbaum_tingey(n, D);
  } else {
    r.method = "Asymptotic one-sample Kolmogorov-Smirnov test";
    const double z = std::sqrt(static_cast<double>(n)) * D;
    r.p_value = 1.0 - (1.0 - (alt == Alt::Two ? kolmogorov_upper(z) : std::exp(-2 * z * z)));
  }
  return r;
}

// =================================================================================================================
// chisq.test

HtestResult chisq_test_gof(const Vec& counts, const Vec& p0) {
  const std::string T = "chisq.test";
  const size_t k = counts.size();
  if (k < 1) throw Error("chisq.test: 'x' must have at least one element");
  for (double v : counts)
    if (!(v >= 0) || !std::isfinite(v)) throw Error("chisq.test: all entries of 'x' must be nonnegative and finite");
  Vec p = p0.empty() ? Vec(k, 1.0) : p0;
  if (p.size() != k) throw Error("chisq.test: 'x' and 'p' lengths differ");
  double sp = 0.0, n = 0.0;
  for (double v : p) {
    if (!(v >= 0) || !std::isfinite(v)) throw Error("chisq.test: probabilities must be nonnegative and finite");
    sp += v;
  }
  if (!(sp > 0)) throw Error("chisq.test: probabilities must not all be zero");
  for (double v : counts) n += v;
  if (!(n > 0)) throw Error("chisq.test: at least one entry of 'x' must be positive");
  HtestResult r;
  r.method = "Chi-squared test for given probabilities";
  r.alternative = "";
  r.statistic_name = "X-squared";
  r.expected.resize(k);
  double x2 = 0.0;
  for (size_t i = 0; i < k; ++i) {
    const double E = n * (p[i] / sp);
    r.expected[i] = E;
    x2 += (counts[i] - E) * (counts[i] - E) / E;
  }
  r.statistic = x2;
  const double df = static_cast<double>(k) - 1;
  r.parameters = {{"df", df}};
  r.p_value = Dist::chisq(df).p(x2, false);
  return r;
}

HtestResult chisq_test(const Matrix& t, const TestOptions& o) {
  if (t.rows == 0 || t.cols == 0) throw Error("chisq.test: empty table");
  if (t.rows == 1 || t.cols == 1) return chisq_test_gof(t.a, {});
  for (double v : t.a)
    if (!(v >= 0) || !std::isfinite(v)) throw Error("chisq.test: all entries of 'x' must be nonnegative and finite");
  const size_t R = t.rows, C = t.cols;
  Vec sr(R, 0.0), sc(C, 0.0);
  double n = 0.0;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j) {
      sr[i] += t(i, j);
      sc[j] += t(i, j);
      n += t(i, j);
    }
  if (!(n > 0)) throw Error("chisq.test: at least one entry of 'x' must be positive");
  const bool yates = o.correct && R == 2 && C == 2;
  HtestResult r;
  r.method = yates ? "Pearson's Chi-squared test with Yates' continuity correction" : "Pearson's Chi-squared test";
  r.alternative = "";
  r.statistic_name = "X-squared";
  r.expected.resize(R * C);
  double corr = 0.0;
  if (yates) {
    corr = 0.5;
    for (size_t i = 0; i < R; ++i)
      for (size_t j = 0; j < C; ++j) corr = std::min(corr, std::fabs(t(i, j) - sr[i] * sc[j] / n));
  }
  double x2 = 0.0;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j) {
      const double E = sr[i] * sc[j] / n;
      r.expected[i * C + j] = E;
      const double dev = std::fabs(t(i, j) - E) - corr;
      x2 += dev * dev / E;
    }
  r.statistic = x2;
  const double df = (R - 1.0) * (C - 1.0);
  r.parameters = {{"df", df}};
  r.p_value = Dist::chisq(df).p(x2, false);
  return r;
}

// =================================================================================================================
// fisher.test

HtestResult fisher_test(const Matrix& t, const TestOptions& o) {
  const std::string T = "fisher.test";
  if (t.rows < 2 || t.cols < 2) throw Error("fisher.test: 'x' must have at least 2 rows and columns");
  std::vector<std::vector<int>> c(t.rows, std::vector<int>(t.cols));  // counts, row-major
  for (size_t i = 0; i < t.rows; ++i)
    for (size_t j = 0; j < t.cols; ++j) c[i][j] = as_count(t(i, j), T, "all entries of 'x'");
  HtestResult r;
  r.method = "Fisher's Exact Test for Count Data";
  const double relErr = 1 + 1e-7;
  if (t.rows == 2 && t.cols == 2) {
    const Alt alt = parse_alt(o.alternative, T);
    check_conf_level(o.conf_level, T);
    if (!(o.or_ >= 0) || std::isnan(o.or_)) throw Error("fisher.test: 'or' must be a single number between 0 and Inf");
    const int x = c[0][0];
    const int m = c[0][0] + c[1][0];  // first column total
    const int nn = c[0][1] + c[1][1];  // second column total
    const int k = c[0][0] + c[0][1];  // first row total
    NcHyper h;
    h.lo = std::max(0, k - nn);
    h.hi = std::min(k, m);
    const Dist hy = Dist::hyper(m, nn, k);
    for (int s = h.lo; s <= h.hi; ++s) h.logdc.push_back(hy.d(s, true));
    r.alternative = alt_name(alt);
    // p-value
    double pv;
    const double OR = o.or_;
    if (alt == Alt::Less) {
      pv = h.p_le(x, OR);
    } else if (alt == Alt::Greater) {
      pv = h.p_ge(x, OR);
    } else if (OR == 0) {
      pv = x == h.lo ? 1.0 : 0.0;
    } else if (std::isinf(OR)) {
      pv = x == h.hi ? 1.0 : 0.0;
    } else {
      const Vec d = h.dens(OR);
      const double dx = d[x - h.lo] * relErr;
      pv = 0.0;
      for (double v : d)
        if (v <= dx) pv += v;
      pv = std::min(1.0, pv);
    }
    r.p_value = pv;
    // conditional MLE
    double est;
    if (h.lo == h.hi) est = kNaN;
    else if (x == h.lo) est = 0.0;
    else if (x == h.hi) est = kInf;
    else {
      est = solve_psi([&](double psi) { return h.mean(psi); }, x, true);
    }
    // exact confidence limits
    auto ncp_upper = [&](double alpha) {  // largest psi with P(X <= x; psi) >= alpha
      if (x == h.hi) return kInf;
      return solve_psi([&](double psi) { return h.p_le(x, psi); }, alpha, false);
    };
    auto ncp_lower = [&](double alpha) {  // smallest psi with P(X >= x; psi) >= alpha
      if (x == h.lo) return 0.0;
      return solve_psi([&](double psi) { return h.p_ge(x, psi); }, alpha, true);
    };
    const double cl = o.conf_level;
    if (alt == Alt::Less) set_ci(r, 0.0, ncp_upper(1 - cl), cl);
    else if (alt == Alt::Greater) set_ci(r, ncp_lower(1 - cl), kInf, cl);
    else {
      const double a = (1 - cl) / 2;
      set_ci(r, ncp_lower(a), ncp_upper(a), cl);
    }
    r.estimates = {{"odds ratio", est}};
    r.null_name = "odds ratio";
    r.null_value = OR;
    return r;
  }
  // r x c: sum the probabilities of all tables with the observed margins that are no more probable than the
  // observed one (relative tolerance 1e-7). Tables are generated column by column; as in the network algorithm of
  // Mehta & Patel (1983), a partial table whose completions are all at or below the threshold contributes their
  // total probability in closed form, and one whose completions are all above it is skipped. The bounds on the
  // completions come from relaxing the row-sum constraints (each remaining column filled independently, as evenly
  // as possible for the most probable completion and as concentrated as possible for the least probable one).
  r.alternative = "two.sided";
  size_t R = t.rows, C = t.cols;
  if (R > C) {  // fewer rows than columns: shorter inner compositions
    std::vector<std::vector<int>> tc(C, std::vector<int>(R));
    for (size_t i = 0; i < R; ++i)
      for (size_t j = 0; j < C; ++j) tc[j][i] = c[i][j];
    c.swap(tc);
    std::swap(R, C);
  }
  std::vector<int> rs(R, 0), cs(C, 0);
  int N = 0;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j) {
      rs[i] += c[i][j];
      cs[j] += c[i][j];
      N += c[i][j];
    }
  // Larger columns first: the bounds become tight sooner.
  std::sort(cs.begin(), cs.end(), std::greater<int>());
  Vec lf(N + 1, 0.0);
  for (int i = 1; i <= N; ++i) lf[i] = lf[i - 1] + std::log(static_cast<double>(i));
  // log P(table) = base - sum lf(n_ij), base = sum lf(rs) + sum lf(cs) - lf(N)
  double base = -lf[N];
  for (int v : rs) base += lf[v];
  for (int v : cs) base += lf[v];
  double lobs = base;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j) lobs -= lf[c[i][j]];
  const double thresh = lobs + std::log(relErr);
  Vec suffix_lf(C + 1, 0.0);  // sum of lf(cs[k]) for k >= j
  std::vector<int> suffix_n(C + 1, 0);
  for (size_t j = C; j-- > 0;) {
    suffix_lf[j] = suffix_lf[j + 1] + lf[cs[j]];
    suffix_n[j] = suffix_n[j + 1] + cs[j];
  }
  double pv = 0.0;
  long long work = 0;
  const long long limit = 20000000LL;
  std::vector<int> rem = rs, caps(R);
  // Bounds on -sum lf(n) over the completions of columns j.. given remaining row sums rem.
  auto bounds = [&](size_t j, double& most, double& least) {
    caps = rem;
    std::sort(caps.begin(), caps.end());
    most = 0.0;
    least = 0.0;
    for (size_t k = j; k < C; ++k) {
      // most probable: water-filling (caps ascending)
      int left = cs[k];
      size_t rows = R;
      for (size_t i = 0; i < R; ++i, --rows) {
        const int share = left / static_cast<int>(rows);
        if (caps[i] <= share) {
          most -= lf[caps[i]];
          left -= caps[i];
        } else {
          const int q = left / static_cast<int>(rows), extra = left % static_cast<int>(rows);
          most -= extra * lf[q + 1] + (static_cast<double>(rows) - extra) * lf[q];
          break;
        }
      }
      // least probable: fill the largest capacities first (caps descending)
      left = cs[k];
      for (size_t i = R; i-- > 0 && left > 0;) {
        const int v = std::min(left, caps[i]);
        least -= lf[v];
        left -= v;
      }
    }
  };
  std::function<void(size_t, double)> node;
  std::function<void(size_t, size_t, int, double)> fill = [&](size_t j, size_t i, int left, double lp) {
    if (i == R - 1) {
      if (left > rem[i]) return;
      rem[i] -= left;
      node(j + 1, lp - lf[left]);
      rem[i] += left;
      return;
    }
    int below = 0;
    for (size_t ii = i + 1; ii < R; ++ii) below += rem[ii];
    const int lo = std::max(0, left - below), hi = std::min(left, rem[i]);
    for (int v = lo; v <= hi; ++v) {
      rem[i] -= v;
      fill(j, i + 1, left - v, lp - lf[v]);
      rem[i] += v;
    }
  };
  node = [&](size_t j, double lp) {
    if (++work > limit) throw Error("fisher.test: table too large for exact enumeration (more than 2e7 partial tables)");
    if (j == C - 1) {  // last column is forced
      double l = lp;
      for (size_t i = 0; i < R; ++i) l -= lf[rem[i]];
      if (l <= thresh) pv += std::exp(l);
      return;
    }
    if (j > 0) {
      double most, least;
      bounds(j, most, least);
      if (lp + least > thresh) return;  // every completion is more probable than the observed table
      if (lp + most <= thresh) {      // every completion counts: closed-form total
        double lsum = lf[suffix_n[j]] - suffix_lf[j];
        for (size_t i = 0; i < R; ++i) lsum -= lf[rem[i]];
        pv += std::exp(lp + lsum);
        return;
      }
    }
    fill(j, 0, cs[j], lp);
  };
  node(0, base);
  r.p_value = std::min(1.0, pv);
  return r;
}

// =================================================================================================================
// binom.test

HtestResult binom_test(double xd, double nd, const TestOptions& o) {
  const std::string T = "binom.test";
  const Alt alt = parse_alt(o.alternative, T);
  check_conf_level(o.conf_level, T);
  // R's order of checks (black box): x first, then n against x.
  if (!(xd >= 0) || !std::isfinite(xd) || std::fabs(xd - std::nearbyint(xd)) > 1e-7 * std::max(1.0, xd))
    throw Error("binom.test: 'x' must be nonnegative and integer");
  if (!(nd >= 1) || !std::isfinite(nd) || std::fabs(nd - std::nearbyint(nd)) > 1e-7 * nd || nd < xd)
    throw Error("binom.test: 'n' must be a positive integer >= 'x'");
  const double x = std::nearbyint(xd), n = std::nearbyint(nd), p = o.p;
  if (!(p >= 0 && p <= 1)) throw Error("binom.test: 'p' must be a single number between 0 and 1");
  const Dist b = Dist::binom(n, p);
  double pv;
  // Tail sums are only requested at nonnegative quantiles (P(X <= -1) = 0 and P(X > -1) = 1 are explicit).
  auto lower_tail = [&](double q) { return q < 0 ? 0.0 : b.p(q); };
  auto upper_tail = [&](double q) { return q < 0 ? 1.0 : b.p(q, false); };  // P(X > q)
  if (alt == Alt::Less) pv = lower_tail(x);
  else if (alt == Alt::Greater) pv = upper_tail(x - 1);
  else {
    // Sum of the probabilities of all outcomes no more likely than x (relative tolerance 1e-7).
    const double relErr = 1 + 1e-7;
    if (p == 0) pv = x == 0 ? 1.0 : 0.0;
    else if (p == 1) pv = x == n ? 1.0 : 0.0;
    else {
      const double d = b.d(x);
      const double m = n * p;
      if (x == m) pv = 1.0;
      else if (x < m) {
        // count outcomes in the upper tail (ceil(m)..n) with probability <= d; they form a tail n-y+1..n
        double y = 0;
        for (double i = std::ceil(m); i <= n; ++i)
          if (b.d(i) <= d * relErr) ++y;
        pv = lower_tail(x) + upper_tail(n - y);
      } else {
        double y = 0;
        for (double i = 0; i <= std::floor(m); ++i)
          if (b.d(i) <= d * relErr) ++y;
        pv = lower_tail(y - 1) + upper_tail(x - 1);
      }
      pv = std::min(1.0, pv);
    }
  }
  // Clopper-Pearson interval.
  const double cl = o.conf_level;
  auto lower = [&](double a) { return x == 0 ? 0.0 : Dist::beta(x, n - x + 1).q(a); };
  auto upper = [&](double a) { return x == n ? 1.0 : Dist::beta(x + 1, n - x).q(1 - a); };
  HtestResult r;
  r.method = "Exact binomial test";
  r.alternative = alt_name(alt);
  r.statistic_name = "number of successes";
  r.statistic = x;
  r.parameters = {{"number of trials", n}};
  r.p_value = pv;
  if (alt == Alt::Less) set_ci(r, 0.0, upper(1 - cl), cl);
  else if (alt == Alt::Greater) set_ci(r, lower(1 - cl), 1.0, cl);
  else set_ci(r, lower((1 - cl) / 2), upper((1 - cl) / 2), cl);
  r.estimates = {{"probability of success", x / n}};
  r.null_name = "probability of success";
  r.null_value = p;
  return r;
}

// =================================================================================================================
// prop.test

HtestResult prop_test(const Vec& x, const Vec& n, const Vec& p, const TestOptions& o) {
  const std::string T = "prop.test";
  Alt alt = parse_alt(o.alternative, T);
  const size_t k = x.size();
  if (k == 0) throw Error("prop.test: not enough data");
  if (n.size() != k) throw Error("prop.test: 'x' and 'n' must have the same length");
  for (size_t i = 0; i < k; ++i) {
    if (!(n[i] > 0) || !std::isfinite(n[i])) throw Error("prop.test: elements of 'n' must be positive");
    if (!(x[i] >= 0) || !std::isfinite(x[i])) throw Error("prop.test: elements of 'x' must be nonnegative");
    if (x[i] > n[i]) throw Error("prop.test: elements of 'x' must not be greater than those of 'n'");
  }
  Vec pp = p;
  if (k == 1 && pp.empty()) pp = {0.5};
  if (!pp.empty()) {
    if (pp.size() != k) throw Error("prop.test: 'p' must have the same length as 'x' and 'n'");
    for (double v : pp)
      if (!(v > 0 && v < 1)) throw Error("prop.test: elements of 'p' must be in (0,1)");
  }
  if (k == 1 || (k == 2 && pp.empty())) check_conf_level(o.conf_level, T);
  if (k > 2 || (k == 2 && !pp.empty())) alt = Alt::Two;
  const bool use_corr = o.correct && k <= 2;
  Vec est(k);
  for (size_t i = 0; i < k; ++i) est[i] = x[i] / n[i];
  HtestResult r;
  r.alternative = alt_name(alt);
  r.statistic_name = "X-squared";
  r.expected.assign(2 * k, 0.0);
  const std::string corr_txt = use_corr ? "with" : "without";
  // Expected counts E (k x 2: successes, failures) and the chi-squared statistic with Yates' correction.
  double x2 = 0.0, df;
  if (pp.empty()) {
    double sx = 0, sn = 0;
    for (size_t i = 0; i < k; ++i) {
      sx += x[i];
      sn += n[i];
    }
    const double ph = sx / sn;
    double yates = use_corr ? 0.5 : 0.0;
    for (size_t i = 0; i < k; ++i) {
      r.expected[2 * i] = n[i] * ph;
      r.expected[2 * i + 1] = n[i] * (1 - ph);
      if (use_corr) {
        yates = std::min(yates, std::fabs(x[i] - n[i] * ph));
        yates = std::min(yates, std::fabs((n[i] - x[i]) - n[i] * (1 - ph)));
      }
    }
    for (size_t i = 0; i < k; ++i) {
      const double o1 = x[i], o2 = n[i] - x[i], e1 = r.expected[2 * i], e2 = r.expected[2 * i + 1];
      const double d1 = std::fabs(o1 - e1) - yates, d2 = std::fabs(o2 - e2) - yates;
      x2 += d1 * d1 / e1 + d2 * d2 / e2;
    }
    df = static_cast<double>(k) - 1;
    r.method = std::to_string(k) + "-sample test for equality of proportions " + corr_txt + " continuity correction";
  } else {
    double yates = use_corr ? 0.5 : 0.0;
    for (size_t i = 0; i < k; ++i) {
      r.expected[2 * i] = n[i] * pp[i];
      r.expected[2 * i + 1] = n[i] * (1 - pp[i]);
      // R caps the correction at |x - n p| only in the one-sample case (oracle: prop.test(c(6, 9), c(20, 30),
      // p = c(0.31, 0.2)) uses the full 1/2 although |6 - 6.2| < 1/2).
      if (use_corr && k == 1) yates = std::min(yates, std::fabs(x[i] - n[i] * pp[i]));
    }
    for (size_t i = 0; i < k; ++i) {
      const double d1 = std::fabs(x[i] - r.expected[2 * i]) - yates;
      const double d2 = std::fabs((n[i] - x[i]) - r.expected[2 * i + 1]) - yates;
      x2 += d1 * d1 / r.expected[2 * i] + d2 * d2 / r.expected[2 * i + 1];
    }
    df = static_cast<double>(k);
    r.method = k == 1 ? "1-sample proportions test " + corr_txt + " continuity correction"
                      : std::to_string(k) + "-sample test for given proportions " + corr_txt + " continuity correction";
  }
  r.statistic = x2;
  r.parameters = {{"df", df}};
  if (alt == Alt::Two) {
    r.p_value = Dist::chisq(df).p(x2, false);
  } else {
    const double delta = k == 1 ? est[0] - pp[0] : est[0] - est[1];
    const double z = (delta > 0 ? 1.0 : delta < 0 ? -1.0 : 0.0) * std::sqrt(x2);
    r.p_value = pnorm(z, alt == Alt::Less);
  }
  const double cl = o.conf_level;
  if (k == 1) {
    // Wilson score interval inverted from the (continuity-corrected) score test.
    const double nn = n[0], e = est[0];
    const double yates = use_corr ? std::min(0.5, std::fabs(x[0] - nn * pp[0])) : 0.0;
    const double z = qnorm(alt == Alt::Two ? (1 + cl) / 2 : cl);
    const double z22n = z * z / (2 * nn);
    auto upper = [&]() {
      const double pc = e + yates / nn;
      if (pc >= 1) return 1.0;
      return (pc + z22n + z * std::sqrt(pc * (1 - pc) / nn + z22n / (2 * nn))) / (1 + 2 * z22n);
    };
    auto lower = [&]() {
      const double pc = e - yates / nn;
      if (pc <= 0) return 0.0;
      return (pc + z22n - z * std::sqrt(pc * (1 - pc) / nn + z22n / (2 * nn))) / (1 + 2 * z22n);
    };
    if (alt == Alt::Two) set_ci(r, std::max(lower(), 0.0), std::min(upper(), 1.0), cl);
    else if (alt == Alt::Less) set_ci(r, 0.0, std::min(upper(), 1.0), cl);
    else set_ci(r, std::max(lower(), 0.0), 1.0, cl);
    r.estimates = {{"p", e}};
    r.null_name = "p";
    r.null_value = pp[0];
  } else {
    for (size_t i = 0; i < k; ++i) r.estimates.push_back({"prop " + std::to_string(i + 1), est[i]});
    if (k == 2 && pp.empty()) {
      const double delta = est[0] - est[1];
      const double sinv = 1 / n[0] + 1 / n[1];
      const double yates = use_corr ? std::min(0.5, std::fabs(delta) / sinv) : 0.0;
      const double z = qnorm(alt == Alt::Two ? (1 + cl) / 2 : cl);
      const double width =
          z * std::sqrt(est[0] * (1 - est[0]) / n[0] + est[1] * (1 - est[1]) / n[1]) + yates * sinv;
      if (alt == Alt::Two) set_ci(r, std::max(delta - width, -1.0), std::min(delta + width, 1.0), cl);
      else if (alt == Alt::Less) set_ci(r, -1.0, std::min(delta + width, 1.0), cl);
      else set_ci(r, std::max(delta - width, -1.0), 1.0, cl);
    }
  }
  return r;
}

}  // namespace rcli
