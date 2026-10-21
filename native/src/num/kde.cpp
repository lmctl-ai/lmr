// Kernel density estimation, bandwidth selection, ECDF and histogram helpers (clean-room; contract in kde.hpp).
// References:
//   Silverman, B. W. (1986) Density Estimation for Statistics and Data Analysis. Chapman & Hall (eq. 3.31, rule of
//     thumb 0.9 min(sd, IQR/1.34) n^-1/5).
//   Scott, D. W. (1979) On optimal and data-based histograms. Biometrika 66, 605-610 (h = 3.5 sd n^-1/3).
//   Scott, D. W. (1992) Multivariate Density Estimation. Wiley (1.06 sd n^-1/5; UCV/BCV with the Gaussian kernel).
//   Freedman, D. and Diaconis, P. (1981) On the histogram as a density estimator: L2 theory. Z. Wahrsch. 57, 453-476.
//   Sturges, H. A. (1926) The choice of a class interval. JASA 21, 65-66.
//   Rudemo, M. (1982) Empirical choice of histograms and kernel density estimators. Scand. J. Statist. 9, 65-78, and
//   Bowman, A. W. (1984) An alternative method of cross-validation for the smoothing of density estimates.
//     Biometrika 71, 353-360 (least-squares / unbiased cross-validation).
//   Scott, D. W. and Terrell, G. R. (1987) Biased and unbiased cross-validation in density estimation. JASA 82,
//     1131-1146 (UCV and BCV criteria for the Gaussian kernel).
//   Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth selection method for kernel density
//     estimation. JRSS B 53, 683-690 (solve-the-equation and direct plug-in; pilot bandwidths a = 0.920 lambda n^-1/7,
//     b = 0.912 lambda n^-1/9 with lambda the IQR, alpha2(h) = 1.357 (S(a)/T(b))^1/7 h^5/7).
//   Wand, M. P. and Jones, M. C. (1995) Kernel Smoothing. Chapman & Hall (section 3.6: two-stage direct plug-in,
//     g = (2.394 / (n T))^1/7 with 2.394 = 6/sqrt(2 pi); kernel functionals R(K), mu2(K)).
//   Brent, R. P. (1973) Algorithms for Minimization without Derivatives (fmin via rcli::optimize; zero as in the
//     public-domain netlib C zeroin, see zeroin_netlib below).
//
// Details found by black-box testing against R 4.3.3 (no R source was read):
//   - ucv/bcv/SJ bin the data: d = 1.01 * (max - min) / nb, bin index trunc(x / d) (toward zero, so the result depends
//     slightly on translation, as ?bw.ucv says), all pairs (including i < j in the same bin, at distance 0) at
//     distance |k| d. This reproduces R's criteria to ~1e-12.
//   - ucv/bcv minimise on [0.1 hmax, hmax] with tol 0.1 * lower, hmax = 1.144 * sd * n^-1/5; bw.SJ uses
//     hmax = 1.144 * min(sd, IQR/1.349) * n^-1/5 instead.
//   - SJ: scale = min(sd, IQR/1.349), a = 1.24 scale n^-1/7, b = 1.23 scale n^-1/9 (Sheather & Jones's constants
//     multiplied by 1.349); the estimates of int f''^2 and int f'''^2 include the i = j terms and divide by n(n-1).
//     "ste": when f(lower) f(upper) > 0 the interval is enlarged alternately (upper * 1.2, then lower / 1.2), and the
//     root is found with tol = 0.1 * (final) lower by the netlib-style zeroin below (same iterates as R).
//   - nrd0's fallbacks when IQR = 0: sd, then |x[0]|, then 1.
//   - nclass.FD's fallback quantiles: a = 1/8, 1/16, ..., 1/512 with the spread scaled by 0.5 / (1 - 2a).
//   - hist's boundary fuzz: 1e-7 times median(diff(breaks)) for >= 5 cells, the smallest positive diff(breaks) for 3-4
//     cells, and diff(range(x)) for 1-2 cells.
//   - pretty(): see the comments in pretty() below.
#include "rcli/kde.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <numeric>

#include "rcli/dd.hpp"
#include "rcli/optim.hpp"
#include "rcli/summary.hpp"

namespace rcli::kde {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kSqrt2Pi = 2.50662827463100050242;  // sqrt(2 pi)

void require_finite(const Vec& x, const char* fn) {
  for (double v : x)
    if (!std::isfinite(v)) throw Error(std::string(fn) + ": non-finite values (NaN/Inf) are not allowed");
}

void require_n2(const Vec& x, const char* fn) {
  require_finite(x, fn);
  if (x.size() < 2) throw Error(std::string(fn) + ": need at least 2 data points");
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

double iqr(const Vec& x) { return summary::quantile(x, 0.75, 7) - summary::quantile(x, 0.25, 7); }

double dnorm0(double u) { return std::exp(-0.5 * u * u) / kSqrt2Pi; }

// ---------------------------------------------------------------------------------------------------------------
// Pairwise distances, binned as R does (nb > 0) or exact (nb == 0). Each entry is (distance, number of pairs i < j).
struct Pairs {
  Vec dist, count;
  size_t n = 0;
};

Pairs make_pairs(const Vec& x, int nb, const char* fn) {
  Pairs p;
  p.n = x.size();
  if (nb < 0) throw Error(std::string(fn) + ": 'nb' must be >= 0");
  if (nb == 0) {
    p.dist.reserve(p.n * (p.n - 1) / 2);
    for (size_t i = 0; i < p.n; ++i)
      for (size_t j = i + 1; j < p.n; ++j) p.dist.push_back(std::fabs(x[i] - x[j]));
    p.count.assign(p.dist.size(), 1.0);
    return p;
  }
  const auto mm = std::minmax_element(x.begin(), x.end());
  const double d = 1.01 * (*mm.second - *mm.first) / nb;
  if (!(d > 0)) {  // all values equal: every pair at distance 0
    p.dist.push_back(0.0);
    p.count.push_back(0.5 * static_cast<double>(p.n) * static_cast<double>(p.n - 1));
    return p;
  }
  std::vector<long long> idx(p.n);
  for (size_t i = 0; i < p.n; ++i) idx[i] = static_cast<long long>(std::trunc(x[i] / d));
  const long long lo = *std::min_element(idx.begin(), idx.end());
  const long long hi = *std::max_element(idx.begin(), idx.end());
  const size_t nbin = static_cast<size_t>(hi - lo + 1);
  std::vector<double> cnt(nbin, 0.0);
  for (long long v : idx) cnt[static_cast<size_t>(v - lo)] += 1.0;
  p.dist.resize(nbin);
  p.count.assign(nbin, 0.0);
  for (size_t k = 0; k < nbin; ++k) p.dist[k] = static_cast<double>(k) * d;
  for (size_t b = 0; b < nbin; ++b) {
    if (cnt[b] == 0.0) continue;
    p.count[0] += cnt[b] * (cnt[b] - 1.0) / 2.0;
    for (size_t c = b + 1; c < nbin; ++c) p.count[c - b] += cnt[b] * cnt[c];
  }
  return p;
}

double hmax_of(const Vec& x) {
  return 1.144 * summary::sd(x) * std::pow(static_cast<double>(x.size()), -0.2);
}

// UCV(h) for the Gaussian kernel (Scott & Terrell 1987 with the n^-2 normalisation of both pair terms):
//   1/(2 n h sqrt(pi)) + 1/(n^2 h sqrt(pi)) sum_{i<j} [exp(-D^2/4) - sqrt(8) exp(-D^2/2)],  D = (x_i - x_j)/h.
double ucv_crit(const Pairs& p, double h) {
  const double n = static_cast<double>(p.n);
  DD s;
  for (size_t k = 0; k < p.dist.size(); ++k) {
    if (p.count[k] == 0.0) continue;
    const double dd = (p.dist[k] / h) * (p.dist[k] / h);
    s += p.count[k] * (std::exp(-dd / 4.0) - std::sqrt(8.0) * std::exp(-dd / 2.0));
  }
  return 1.0 / (2.0 * n * h * std::sqrt(kPi)) + static_cast<double>(s) / (n * n * h * std::sqrt(kPi));
}

// BCV(h) for the Gaussian kernel (Scott & Terrell 1987):
//   1/(2 n h sqrt(pi)) + 1/(64 n^2 h sqrt(pi)) sum_{i<j} (D^4 - 12 D^2 + 12) exp(-D^2/4).
double bcv_crit(const Pairs& p, double h) {
  const double n = static_cast<double>(p.n);
  DD s;
  for (size_t k = 0; k < p.dist.size(); ++k) {
    if (p.count[k] == 0.0) continue;
    const double dd = (p.dist[k] / h) * (p.dist[k] / h);
    s += p.count[k] * ((dd * dd - 12.0 * dd + 12.0) * std::exp(-dd / 4.0));
  }
  return 1.0 / (2.0 * n * h * std::sqrt(kPi)) + static_cast<double>(s) / (64.0 * n * n * h * std::sqrt(kPi));
}

// Kernel estimates of int f''^2 (phi4) and int f'''^2 (phi6) with pilot bandwidth g (Sheather & Jones 1991, eq. 5-6):
//   S(g) =  1/(n(n-1) g^5) sum_{i,j} phi^(4)((x_i - x_j)/g),  T(g) = -1/(n(n-1) g^7) sum_{i,j} phi^(6)(...),
// the double sums running over all i, j (the diagonal included).
double sj_S(const Pairs& p, double g) {
  const double n = static_cast<double>(p.n);
  DD s = n * 3.0 / kSqrt2Pi;  // diagonal: phi4(0) = 3 phi(0)
  for (size_t k = 0; k < p.dist.size(); ++k) {
    if (p.count[k] == 0.0) continue;
    const double u = p.dist[k] / g, u2 = u * u;
    s += 2.0 * p.count[k] * ((u2 * u2 - 6.0 * u2 + 3.0) * dnorm0(u));
  }
  return static_cast<double>(s) / (n * (n - 1.0) * std::pow(g, 5.0));
}
double sj_T(const Pairs& p, double g) {
  const double n = static_cast<double>(p.n);
  DD s = n * -15.0 / kSqrt2Pi;  // phi6(0) = -15 phi(0)
  for (size_t k = 0; k < p.dist.size(); ++k) {
    if (p.count[k] == 0.0) continue;
    const double u = p.dist[k] / g, u2 = u * u;
    s += 2.0 * p.count[k] * ((u2 * u2 * u2 - 15.0 * u2 * u2 + 45.0 * u2 - 15.0) * dnorm0(u));
  }
  return -static_cast<double>(s) / (n * (n - 1.0) * std::pow(g, 7.0));
}

// Brent's zero finder (Brent 1973, ch. 4) in the formulation of the public-domain netlib C zeroin: the interpolated
// step is accepted only if it is less than half of the most recent step (rather than of the step before it, as in
// Brent's ALGOL procedure), with tol_act = 2 eps |b| + tol / 2. R's uniroot() follows the same iterates (verified by
// tracing evaluation points), which matters here because bw.SJ's default tolerance (0.1 * lower) is coarse.
double zeroin_netlib(const Fn1& f, double a, double b, double fa, double fb, double tol, int maxit) {
  double c = a, fc = fa;
  if (fa == 0.0) return a;
  if (fb == 0.0) return b;
  for (int it = 0; it < maxit; ++it) {
    const double prev_step = b - a;
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    const double tol_act = 2.0 * kEps * std::fabs(b) + tol / 2.0;
    double new_step = (c - b) / 2.0;
    if (std::fabs(new_step) <= tol_act || fb == 0.0) return b;
    if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
      const double cb = c - b;
      double p, q;
      if (a == c) {  // secant
        const double t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {  // inverse quadratic interpolation
        const double qa = fa / fc, t1 = fb / fc, t2 = fb / fa;
        p = t2 * (cb * qa * (qa - t1) - (b - a) * (t1 - 1.0));
        q = (qa - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0) q = -q;
      else p = -p;
      if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2.0) && p < std::fabs(prev_step * q / 2.0)) new_step = p / q;
    }
    if (std::fabs(new_step) < tol_act) new_step = new_step > 0 ? tol_act : -tol_act;
    a = b;
    fa = fb;
    b += new_step;
    fb = f(b);
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
  return b;  // not converged in maxit steps (R warns); the best approximation
}

double cv_select(const Vec& x, const BwOptions& opt, bool unbiased, const char* fn) {
  require_n2(x, fn);
  const Pairs p = make_pairs(x, opt.nb, fn);
  const double hmax = hmax_of(x);
  const double lo = std::isnan(opt.lower) ? 0.1 * hmax : opt.lower;
  const double up = std::isnan(opt.upper) ? hmax : opt.upper;
  if (!(lo > 0 && up > lo)) throw Error(std::string(fn) + ": invalid 'lower'/'upper' (need 0 < lower < upper)");
  const double tol = std::isnan(opt.tol) ? 0.1 * lo : opt.tol;
  const auto r = unbiased ? optimize([&](double h) { return ucv_crit(p, h); }, lo, up, false, tol)
                          : optimize([&](double h) { return bcv_crit(p, h); }, lo, up, false, tol);
  return r.minimum;
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// Rules of thumb

double bw_nrd0(const Vec& x) {
  require_n2(x, "bw_nrd0");
  const double hi = summary::sd(x);
  double lo = std::min(hi, iqr(x) / 1.34);
  if (!(lo > 0)) {
    lo = hi;
    if (!(lo > 0)) lo = std::fabs(x[0]);
    if (!(lo > 0)) lo = 1.0;
  }
  return 0.9 * lo * std::pow(static_cast<double>(x.size()), -0.2);
}

double bw_nrd(const Vec& x) {
  require_n2(x, "bw_nrd");
  const double s = std::min(summary::sd(x), iqr(x) / 1.34);
  return 1.06 * s * std::pow(static_cast<double>(x.size()), -0.2);
}

double bw_ucv(const Vec& x, const BwOptions& opt) { return cv_select(x, opt, true, "bw_ucv"); }
double bw_bcv(const Vec& x, const BwOptions& opt) { return cv_select(x, opt, false, "bw_bcv"); }

double bw_sj(const Vec& x, const std::string& method, const BwOptions& opt) {
  const char* fn = "bw_sj";
  require_n2(x, fn);
  const std::string m = lower(method);
  if (m != "ste" && m != "dpi") throw Error("bw_sj: 'method' must be \"ste\" or \"dpi\"");
  const Pairs p = make_pairs(x, opt.nb, fn);
  const double n = static_cast<double>(x.size());
  const double scale = std::min(summary::sd(x), iqr(x) / 1.349);
  const double a = 1.24 * scale * std::pow(n, -1.0 / 7.0);
  const double b = 1.23 * scale * std::pow(n, -1.0 / 9.0);
  const double c1 = 1.0 / (2.0 * std::sqrt(kPi) * n);
  const double TDb = sj_T(p, b);
  if (!std::isfinite(TDb) || TDb <= 0) throw Error("bw_sj: sample is too sparse to find TD");
  if (m == "dpi") {
    const double g = std::pow(2.394 / (n * TDb), 1.0 / 7.0);
    return std::pow(c1 / sj_S(p, g), 1.0 / 5.0);
  }
  const double SDa = sj_S(p, a);
  const double alph = 1.357 * std::pow(SDa / TDb, 1.0 / 7.0);
  if (!std::isfinite(alph)) throw Error("bw_sj: sample is too sparse to find alph2");
  auto fSD = [&](double h) { return std::pow(c1 / sj_S(p, alph * std::pow(h, 5.0 / 7.0)), 1.0 / 5.0) - h; };
  const double hmax = 1.144 * scale * std::pow(n, -0.2);  // SJ's hmax uses the robust scale (black-box)
  const bool given = !std::isnan(opt.lower) || !std::isnan(opt.upper);
  double lo = std::isnan(opt.lower) ? 0.1 * hmax : opt.lower;
  double up = std::isnan(opt.upper) ? hmax : opt.upper;
  if (!(lo > 0 && up > lo)) throw Error("bw_sj: invalid 'lower'/'upper' (need 0 < lower < upper)");
  double flo = fSD(lo), fup = fSD(up);
  for (int itry = 1; flo * fup > 0; ++itry) {
    if (given || itry > 99) throw Error("bw_sj: no solution in the specified range of bandwidths");
    if (itry % 2) { up *= 1.2; fup = fSD(up); }
    else { lo /= 1.2; flo = fSD(lo); }
  }
  if (std::isnan(flo) || std::isnan(fup)) throw Error("bw_sj: non-finite equation value");
  const double tol = std::isnan(opt.tol) ? 0.1 * lo : opt.tol;
  return zeroin_netlib(fSD, lo, up, flo, fup, tol, 1000);
}

double bandwidth(const Vec& x, const std::string& rule) {
  const std::string r = lower(rule);
  if (r == "nrd0") return bw_nrd0(x);
  if (r == "nrd") return bw_nrd(x);
  if (r == "ucv") return bw_ucv(x);
  if (r == "bcv") return bw_bcv(x);
  if (r == "sj" || r == "sj-ste") return bw_sj(x, "ste");
  if (r == "sj-dpi") return bw_sj(x, "dpi");
  throw Error("bandwidth: unknown bandwidth rule '" + rule + "'");
}

// ---------------------------------------------------------------------------------------------------------------
// Kernels with standard deviation bw. For a kernel K0 on [-1, 1] with variance s0^2 the support half-width is
// a = bw / s0: epanechnikov s0^2 = 1/5, rectangular 1/3, triangular 1/6, biweight 1/7, cosine (1 + cos(pi t))/2:
// 1/3 - 2/pi^2, optcosine (pi/4) cos(pi t/2): 1 - 8/pi^2.

namespace {
const char* const kKernels[] = {"gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine",
                                "optcosine"};

int kernel_id(const std::string& kernel) {
  if (kernel.empty()) throw Error("density: 'kernel' must be one of gaussian, epanechnikov, rectangular, triangular, "
                                  "biweight, cosine, optcosine");
  int found = -1;
  for (int i = 0; i < 7; ++i) {
    const std::string k = kKernels[i];
    if (k == kernel) return i;
    if (k.compare(0, kernel.size(), kernel) == 0) {
      if (found >= 0) throw Error("density: ambiguous kernel '" + kernel + "'");
      found = i;
    }
  }
  if (found < 0)
    throw Error("density: unknown kernel '" + kernel + "' (gaussian, epanechnikov, rectangular, triangular, biweight, "
                "cosine, optcosine)");
  return found;
}

// Half-width multiplier: a = bw * width_factor (infinite for the Gaussian).
double width_factor(int id) {
  switch (id) {
    case 1: return std::sqrt(5.0);
    case 2: return std::sqrt(3.0);
    case 3: return std::sqrt(6.0);
    case 4: return std::sqrt(7.0);
    case 5: return 1.0 / std::sqrt(1.0 / 3.0 - 2.0 / (kPi * kPi));
    case 6: return 1.0 / std::sqrt(1.0 - 8.0 / (kPi * kPi));
    default: return kInf;
  }
}

double kval(int id, double u, double bw) {
  if (id == 0) {
    const double z = u / bw;
    return std::exp(-0.5 * z * z) / (kSqrt2Pi * bw);
  }
  const double a = bw * width_factor(id);
  const double au = std::fabs(u);
  if (!(au < a)) return 0.0;
  const double t = au / a;
  switch (id) {
    case 1: return 0.75 * (1.0 - t * t) / a;
    case 2: return 0.5 / a;
    case 3: return (1.0 - t) / a;
    case 4: { const double s = 1.0 - t * t; return 15.0 / 16.0 * s * s / a; }
    case 5: return (1.0 + std::cos(kPi * t)) / (2.0 * a);
    default: return kPi / 4.0 * std::cos(kPi * t / 2.0) / a;
  }
}

// Sum over sorted data within the kernel's support (|u| <= 40 bw for the Gaussian, where exp(-800) is already 0).
double kde_sorted(const Vec& xs, const Vec& w, double x0, double bw, int id) {
  const double reach = id == 0 ? 40.0 * bw : bw * width_factor(id);
  const auto first = std::lower_bound(xs.begin(), xs.end(), x0 - reach);
  const auto last = std::upper_bound(first, xs.end(), x0 + reach);
  DD s;
  for (auto it = first; it != last; ++it) {
    const size_t i = static_cast<size_t>(it - xs.begin());
    s += w[i] * kval(id, x0 - *it, bw);
  }
  return static_cast<double>(s);
}

// Sort x (and its weights alongside); empty weights become 1/n.
void sorted_with_weights(const Vec& x, const Vec& weights, Vec& xs, Vec& ws, const char* fn) {
  require_finite(x, fn);
  if (x.empty()) throw Error(std::string(fn) + ": need at least 1 data point");
  const size_t n = x.size();
  if (!weights.empty()) {
    if (weights.size() != n) throw Error(std::string(fn) + ": 'x' and 'weights' have unequal length");
    for (double v : weights)
      if (!std::isfinite(v) || v < 0) throw Error(std::string(fn) + ": 'weights' must all be finite and >= 0");
  }
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](size_t i, size_t j) { return x[i] < x[j]; });
  xs.resize(n);
  ws.resize(n);
  for (size_t k = 0; k < n; ++k) {
    xs[k] = x[ord[k]];
    ws[k] = weights.empty() ? 1.0 / static_cast<double>(n) : weights[ord[k]];
  }
}
}  // namespace

std::string kernel_name(const std::string& kernel) { return kKernels[kernel_id(kernel)]; }

double kernel_value(const std::string& kernel, double u, double bw) {
  if (!(bw > 0)) throw Error("kernel_value: 'bw' must be positive");
  return kval(kernel_id(kernel), u, bw);
}

double kernel_rk(const std::string& kernel) {
  const int id = kernel_id(kernel);
  if (id == 0) return 1.0 / (2.0 * std::sqrt(kPi));
  const double a = width_factor(id);
  switch (id) {
    case 1: return 3.0 / (5.0 * a);          // (9/16) * (16/15) / a
    case 2: return 1.0 / (2.0 * a);
    case 3: return 2.0 / (3.0 * a);
    case 4: return 5.0 / (7.0 * a);          // (225/256) * (256/315) / a
    case 5: return 3.0 / (4.0 * a);
    default: return kPi * kPi / (16.0 * a);
  }
}

double kde_at(const Vec& x, double x0, double bw, const std::string& kernel, const Vec& weights) {
  if (!(bw > 0) || !std::isfinite(bw)) throw Error("kde_at: 'bw' must be positive and finite");
  Vec xs, ws;
  sorted_with_weights(x, weights, xs, ws, "kde_at");
  return kde_sorted(xs, ws, x0, bw, kernel_id(kernel));
}

double Density::at(double x0) const { return kde_sorted(data, w, x0, bw, kernel_id(kernel)); }

Density density(const Vec& x, const DensityOptions& opt) {
  Density d;
  const int id = kernel_id(opt.kernel);
  d.kernel = kKernels[id];
  sorted_with_weights(x, opt.weights, d.data, d.w, "density");
  d.n = x.size();
  if (opt.n < 1) throw Error("density: 'n' must be >= 1");
  if (!(opt.adjust > 0) || !std::isfinite(opt.adjust)) throw Error("density: 'adjust' must be positive and finite");
  double bw = opt.bw;
  if (std::isnan(bw)) bw = bandwidth(x, opt.bw_rule);
  if (!std::isfinite(bw)) throw Error("density: non-finite 'bw'");
  bw *= opt.adjust;
  if (!(bw > 0)) throw Error("density: 'bw' is not positive.");
  d.bw = bw;
  const double from = std::isnan(opt.from) ? d.data.front() - opt.cut * bw : opt.from;
  const double to = std::isnan(opt.to) ? d.data.back() + opt.cut * bw : opt.to;
  if (!std::isfinite(from) || !std::isfinite(to)) throw Error("density: non-finite 'from' or 'to'");
  // Grid as R's seq.int(from, to, length.out = n): from + i * by with by = (to - from) / (n - 1).
  const int n = opt.n;
  d.x.resize(static_cast<size_t>(n));
  d.y.resize(static_cast<size_t>(n));
  const double by = n > 1 ? (to - from) / (n - 1) : 0.0;
  for (int i = 0; i < n; ++i) d.x[static_cast<size_t>(i)] = from + i * by;
  if (n > 1) d.x.back() = to;
  for (size_t i = 0; i < d.x.size(); ++i) d.y[i] = kde_sorted(d.data, d.w, d.x[i], bw, id);
  return d;
}

// ---------------------------------------------------------------------------------------------------------------
// ECDF

Ecdf ecdf(const Vec& x) {
  for (double v : x)
    if (std::isnan(v)) throw Error("ecdf: missing values (NaN) are not allowed");
  if (x.empty()) throw Error("ecdf: 'x' must have 1 or more non-missing values");
  Ecdf e;
  e.sorted = x;
  std::sort(e.sorted.begin(), e.sorted.end());
  const double n = static_cast<double>(e.sorted.size());
  for (size_t i = 0; i < e.sorted.size(); ++i) {
    if (i + 1 < e.sorted.size() && e.sorted[i + 1] == e.sorted[i]) continue;
    e.knots.push_back(e.sorted[i]);
    e.cumprob.push_back(static_cast<double>(i + 1) / n);
  }
  return e;
}

double Ecdf::operator()(double t) const {
  if (std::isnan(t)) return kNaN;
  const auto it = std::upper_bound(sorted.begin(), sorted.end(), t);
  return static_cast<double>(it - sorted.begin()) / static_cast<double>(sorted.size());
}

Vec Ecdf::operator()(const Vec& t) const {
  Vec r(t.size());
  for (size_t i = 0; i < t.size(); ++i) r[i] = (*this)(t[i]);
  return r;
}

double Ecdf::quantile(double p, int type) const { return summary::quantile(sorted, p, type); }

// ---------------------------------------------------------------------------------------------------------------
// Histograms

Vec hist_counts(const Vec& x, const Vec& breaks, bool right, bool include_lowest, double fuzz) {
  require_finite(x, "hist_counts");
  require_finite(breaks, "hist_counts");
  const size_t nb = breaks.size();
  if (nb < 2) throw Error("hist_counts: need at least 2 breaks");
  for (size_t i = 1; i < nb; ++i)
    if (!(breaks[i] > breaks[i - 1])) throw Error("hist_counts: 'breaks' are not strictly increasing");
  if (!(fuzz >= 0)) throw Error("hist_counts: 'fuzz' must be >= 0");
  // Fuzz scale (black-box): >= 5 cells: median diff; 3-4 cells: smallest diff; 1-2 cells: range of x.
  Vec h(nb - 1);
  for (size_t i = 0; i + 1 < nb; ++i) h[i] = breaks[i + 1] - breaks[i];
  double scale;
  if (nb > 5) scale = summary::median(h);
  else if (nb > 3) scale = *std::min_element(h.begin(), h.end());
  else if (x.empty()) scale = 0.0;
  else {
    const auto mm = std::minmax_element(x.begin(), x.end());
    scale = *mm.second - *mm.first;
  }
  const double diddle = fuzz * scale;
  Vec fb(breaks);
  if (right) {
    for (size_t i = 0; i < nb; ++i) fb[i] += diddle;
    if (include_lowest) fb[0] = breaks[0] - diddle;
  } else {
    for (size_t i = 0; i < nb; ++i) fb[i] -= diddle;
    if (include_lowest) fb[nb - 1] = breaks[nb - 1] + diddle;
  }
  Vec counts(nb - 1, 0.0);
  size_t counted = 0;
  for (double v : x) {
    size_t cell;
    if (right) {  // (fb[i], fb[i+1]], first closed when include_lowest
      const auto it = std::lower_bound(fb.begin(), fb.end(), v);  // first fb >= v
      if (it == fb.end()) continue;
      size_t j = static_cast<size_t>(it - fb.begin());
      if (j == 0) {
        if (!(include_lowest && v == fb[0])) continue;
        j = 1;
      }
      cell = j - 1;
    } else {  // [fb[i], fb[i+1]), last closed when include_lowest
      const auto it = std::upper_bound(fb.begin(), fb.end(), v);  // first fb > v
      size_t j = static_cast<size_t>(it - fb.begin());
      if (j == 0) continue;
      if (j == nb) {
        if (!(include_lowest && v == fb[nb - 1])) continue;
        j = nb - 1;
      }
      cell = j - 1;
    }
    counts[cell] += 1.0;
    ++counted;
  }
  if (counted != x.size()) throw Error("hist_counts: some 'x' not counted; maybe 'breaks' do not span range of 'x'");
  return counts;
}

int nclass_sturges(const Vec& x) {
  if (x.empty()) throw Error("nclass_sturges: 'x' is empty");
  return static_cast<int>(std::ceil(std::log2(static_cast<double>(x.size())) + 1.0));
}

int nclass_scott(const Vec& x) {
  require_finite(x, "nclass_scott");
  if (x.size() < 2) return 1;
  const double h = 3.5 * std::sqrt(summary::var(x)) * std::pow(static_cast<double>(x.size()), -1.0 / 3.0);
  if (!(h > 0)) return 1;
  const auto mm = std::minmax_element(x.begin(), x.end());
  return std::max(1, static_cast<int>(std::ceil((*mm.second - *mm.first) / h)));
}

namespace {
// R's signif(v, digits) via decimal formatting (round half to even on the exact binary value).
double signif(double v, int digits) {
  if (v == 0 || !std::isfinite(v)) return v;
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*e", digits - 1, v);
  return std::strtod(buf, nullptr);
}
}  // namespace

int nclass_fd(const Vec& x, int digits) {
  require_finite(x, "nclass_fd");
  if (x.size() < 2) return 1;
  Vec xs(x.size());
  for (size_t i = 0; i < x.size(); ++i) xs[i] = signif(x[i], digits);
  double h = iqr(xs);
  // Fallback when the IQR is 0: quantile spreads at a = 1/8, ..., 1/512, scaled by 0.5 / (1 - 2a).
  for (double a = 0.125; h == 0 && a >= 1.0 / 512.0; a /= 2) {
    const double s = summary::quantile(xs, 1.0 - a, 7) - summary::quantile(xs, a, 7);
    h = s * 0.5 / (1.0 - 2.0 * a);
  }
  if (!(h > 0)) return nclass_scott(x);
  h = 2.0 * h * std::pow(static_cast<double>(x.size()), -1.0 / 3.0);
  const auto mm = std::minmax_element(x.begin(), x.end());
  return std::max(1, static_cast<int>(std::ceil((*mm.second - *mm.first) / h)));
}

// pretty(), characterised by black-box tests of R 4.3.3 (?pretty describes only the first steps):
//   - dx = hi - lo is "small" when dx < cell * 3 * max(1, n) * (1 + 1/(1 + h)) * eps, cell = max(|lo|, |hi|)
//     (measured thresholds 21.5 eps (n = 5, h = 1.5), 25.5 (h = 0.5), 18.5 (h = 3), proportional to n).
//     Then c = cell * shrink / max(1, min.n) as ?pretty documents (cell = 1 when it is 0). For |x| > 10 R uses a
//     smaller cell than this; that undocumented detail is deliberately not reproduced (clean room), so pretty() of a
//     (near-)constant range above 10 can differ from R's.
//     Otherwise c = dx / n.
//   - b = 10^floor(log10 c); the unit starts at b and moves to 2b when 2b - c < h (c - unit), then to 5b when
//     5b - c < f (c - unit), then to 10b when 10b - c < h (c - unit) (switch points c/b = 1.4, 2.8, 7 at the defaults).
//   - ns = floor(lo/unit + 1e-10), nu = ceil(hi/unit - 1e-10); the ends are widened to lo/hi themselves when
//     ns*unit > lo or nu*unit < hi (so the result always covers [lo, hi]).
//   - if nu - ns < min.n, k = min.n - (nu - ns) cells are added: all above when lo = 0 < hi, all below when
//     lo < 0 = hi, otherwise k/2 (rounded down) above and the rest below when ns >= 0, the other way round when ns < 0.
//   - values: seq(l, u, length.out = nu - ns + 1), with the value where the grid crosses 0 set to exactly 0.
Vec pretty(double lo, double hi, int n, int min_n, double shrink_sml, double high_u_bias, double u5_bias) {
  if (!std::isfinite(lo) || !std::isfinite(hi)) throw Error("pretty: non-finite range");
  if (lo > hi) std::swap(lo, hi);
  if (n < 0) throw Error("pretty: invalid 'n' argument");
  if (min_n < 0) min_n = n / 3;
  if (min_n > n) throw Error("pretty: invalid 'min.n' argument");
  if (!(shrink_sml > 0)) throw Error("pretty: 'shrink.sml' must be positive");
  if (!(high_u_bias >= 0)) throw Error("pretty: 'high.u.bias' must be non-negative");
  if (std::isnan(u5_bias)) u5_bias = 0.5 + 1.5 * high_u_bias;
  if (!(u5_bias >= 0)) throw Error("pretty: 'u5.bias' must be non-negative");
  const double h = high_u_bias, f = u5_bias;
  const double dx = hi - lo;
  double cell = std::max(std::fabs(lo), std::fabs(hi));
  const double ndiv = std::max(1, n);
  const bool small = dx == 0 || dx < cell * 3.0 * ndiv * (1.0 + 1.0 / (1.0 + h)) * kEps;
  double c;
  if (small) {
    if (cell == 0) cell = 1.0;
    c = cell * shrink_sml / std::max(1, min_n);
  } else {
    c = dx / ndiv;
  }
  const double b = std::pow(10.0, std::floor(std::log10(c)));
  double unit = b;
  if (2 * b - c < h * (c - unit)) unit = 2 * b;
  if (5 * b - c < f * (c - unit)) unit = 5 * b;
  if (10 * b - c < h * (c - unit)) unit = 10 * b;
  double ns = std::floor(lo / unit + 1e-10);
  double nu = std::ceil(hi / unit - 1e-10);
  if (nu - ns < min_n) {
    const double k = min_n - (nu - ns);
    const double half = std::floor(k / 2), rest = k - half;
    if (lo == 0 && hi > 0) nu += k;
    else if (hi == 0 && lo < 0) ns -= k;
    else if (ns >= 0) { nu += half; ns -= rest; }
    else { ns -= half; nu += rest; }
  }
  double l = ns * unit, u = nu * unit;
  if (l > lo) l = lo;
  if (u < hi) u = hi;
  const int cells = static_cast<int>(nu - ns);
  Vec out(static_cast<size_t>(cells + 1));
  if (cells == 0) { out[0] = l; return out; }
  const double by = (u - l) / cells;
  for (int i = 0; i <= cells; ++i) {
    const double v = l + i * by;
    out[static_cast<size_t>(i)] = std::fabs(v) < 1e-10 * by ? 0.0 : v;  // exact 0 where the grid crosses zero
  }
  out[0] = l;
  out.back() = u;
  return out;
}

Vec hist_breaks(const Vec& x, const std::string& rule) {
  require_finite(x, "hist_breaks");
  if (x.empty()) throw Error("hist_breaks: 'x' is empty");
  const std::string r = lower(rule);
  int nc;
  if (r == "sturges") nc = nclass_sturges(x);
  else if (r == "scott") nc = nclass_scott(x);
  else if (r == "fd" || r == "freedman-diaconis") nc = nclass_fd(x);
  else throw Error("hist_breaks: unknown breaks algorithm '" + rule + "'");
  nc = std::min(nc, 1000000);
  const auto mm = std::minmax_element(x.begin(), x.end());
  return pretty(*mm.first, *mm.second, nc, 1);
}

Histogram hist(const Vec& x, const Vec& breaks, bool right, bool include_lowest) {
  Histogram hg;
  hg.breaks = breaks;
  hg.counts = hist_counts(x, breaks, right, include_lowest);
  const size_t nc = hg.counts.size();
  const double n = static_cast<double>(x.size());
  hg.density.resize(nc);
  hg.mids.resize(nc);
  Vec w(nc);
  for (size_t i = 0; i < nc; ++i) {
    w[i] = breaks[i + 1] - breaks[i];
    hg.density[i] = hg.counts[i] / (n * w[i]);
    hg.mids[i] = 0.5 * (breaks[i] + breaks[i + 1]);
  }
  // R: diff(range(diff(breaks))) < 1e-7 * mean(diff(breaks))
  const auto mm = std::minmax_element(w.begin(), w.end());
  hg.equidist = *mm.second - *mm.first < 1e-7 * summary::mean(w);
  return hg;
}

Histogram hist(const Vec& x, const std::string& rule) { return hist(x, hist_breaks(x, rule)); }

}  // namespace rcli::kde
