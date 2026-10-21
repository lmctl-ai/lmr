// Descriptive statistics (clean-room; see summary.hpp for the contract).
// References:
//   Neumaier, A. (1974) Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen. ZAMM 54, 39-51.
//   Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages. American Statistician 50, 361-365.
//   Kendall, M. G. (1945) The treatment of ties in ranking problems. Biometrika 33, 239-251 (tau-b).
//   Knight, W. R. (1966) A computer method for calculating Kendall's tau with ungrouped data. JASA 61, 436-439.
#include "rcli/summary.hpp"

#include "rcli/dd.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>

namespace rcli::summary {
namespace {

void require_finite_input(const Vec& x, const char* fn) {
  for (double v : x)
    if (std::isnan(v)) throw Error(std::string(fn) + ": missing values (NaN) are not allowed");
}

// Sums in double-double (rcli/dd.hpp), where R accumulates in long double: the same result on every platform, and at
// least as accurate as R's 80-bit x86-64 accumulator.
DD ksum(const Vec& x) {
  DD s;
  for (double v : x) s += v;
  return s;
}

// Two-pass mean as R's: m = sum/n, then m += sum(x - m)/n.
DD mean_dd(const Vec& x) {
  const DD n = static_cast<double>(x.size());
  DD m = ksum(x) / n;
  if (std::isfinite(m.hi)) {
    DD t;
    for (double v : x) t += DD(v) - m;
    m += t / n;
  }
  return m;
}

double pearson(const Vec& x, const Vec& y) {
  const size_t n = x.size();
  if (n < 2) return kNaN;
  const DD mx = mean_dd(x), my = mean_dd(y);
  DD sxx, syy, sxy;
  for (size_t i = 0; i < n; ++i) {
    const DD dx = DD(x[i]) - mx, dy = DD(y[i]) - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  if (sxx.hi == 0.0 || syy.hi == 0.0) return kNaN;  // R: "the standard deviation is zero" -> NA
  double r = static_cast<double>(sxy / (sqrt(sxx) * sqrt(syy)));
  return std::clamp(r, -1.0, 1.0);
}

// Merge sort that counts strict inversions (a[i] > a[j] for i < j).
std::int64_t merge_count(std::vector<double>& a, std::vector<double>& tmp, size_t lo, size_t hi) {
  if (hi - lo < 2) return 0;
  const size_t mid = lo + (hi - lo) / 2;
  std::int64_t inv = merge_count(a, tmp, lo, mid) + merge_count(a, tmp, mid, hi);
  size_t i = lo, j = mid, k = lo;
  while (i < mid && j < hi) {
    if (a[j] < a[i]) {
      inv += static_cast<std::int64_t>(mid - i);
      tmp[k++] = a[j++];
    } else {
      tmp[k++] = a[i++];
    }
  }
  while (i < mid) tmp[k++] = a[i++];
  while (j < hi) tmp[k++] = a[j++];
  std::copy(tmp.begin() + lo, tmp.begin() + hi, a.begin() + lo);
  return inv;
}

// Kendall's tau-b by Knight's O(n log n) algorithm.
double kendall(const Vec& x, const Vec& y) {
  const size_t n = x.size();
  if (n < 2) return kNaN;
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b] || (x[a] == x[b] && y[a] < y[b]); });
  auto pairs = [](std::int64_t t) { return t * (t - 1) / 2; };
  const std::int64_t n0 = pairs(static_cast<std::int64_t>(n));
  std::int64_t n1 = 0, n3 = 0;
  for (size_t i = 0; i < n;) {
    size_t j = i;
    while (j < n && x[idx[j]] == x[idx[i]]) ++j;
    n1 += pairs(static_cast<std::int64_t>(j - i));
    for (size_t k = i; k < j;) {  // joint ties within the x group
      size_t l = k;
      while (l < j && y[idx[l]] == y[idx[k]]) ++l;
      n3 += pairs(static_cast<std::int64_t>(l - k));
      k = l;
    }
    i = j;
  }
  std::vector<double> ys(n), tmp(n);
  for (size_t i = 0; i < n; ++i) ys[i] = y[idx[i]];
  const std::int64_t swaps = merge_count(ys, tmp, 0, n);
  std::int64_t n2 = 0;
  for (size_t i = 0; i < n;) {
    size_t j = i;
    while (j < n && ys[j] == ys[i]) ++j;
    n2 += pairs(static_cast<std::int64_t>(j - i));
    i = j;
  }
  const std::int64_t s = n0 - n1 - n2 + n3 - 2 * swaps;
  // The pair counts are exact in double up to 2^53; the square roots and quotient go through double-double.
  const DD den = sqrt(DD(static_cast<double>(n0 - n1))) * sqrt(DD(static_cast<double>(n0 - n2)));
  if (den.hi == 0.0) return kNaN;
  return std::clamp(static_cast<double>(DD(static_cast<double>(s)) / den), -1.0, 1.0);
}

double quantile_sorted(const Vec& s, double p, int type) {
  const size_t n = s.size();
  if (std::isnan(p)) return kNaN;
  if (n == 0) return kNaN;
  const double fuzz = 4.0 * kEps;
  const double nd = static_cast<double>(n);
  // x[j] with 1-based j, clamped to the sample (x[0] = x[1], x[n+1] = x[n]).
  auto at = [&](double j) {
    if (j < 1) return s.front();
    if (j > nd) return s.back();
    return s[static_cast<size_t>(j) - 1];
  };
  auto mix = [&](double j, double h) {
    const double lo = at(j), hi = at(j + 1);
    if (h == 0.0 || lo == hi) return lo;
    if (h == 1.0) return hi;
    return (1.0 - h) * lo + h * hi;
  };
  if (type <= 3) {
    const double nppm = type == 3 ? nd * p - 0.5 : nd * p;
    const double j = std::floor(nppm + fuzz);
    const bool frac = nppm > j + fuzz;  // g > 0 beyond the fuzz
    double h;
    if (type == 1) h = frac ? 1.0 : 0.0;
    else if (type == 2) h = frac ? 1.0 : 0.5;
    else h = (frac || std::fmod(j, 2.0) != 0.0) ? 1.0 : 0.0;
    return mix(j, h);
  }
  if (type == 7) {
    const double index = 1.0 + (nd - 1.0) * p;
    const double lo = std::floor(index);
    return mix(lo, index - lo);
  }
  double a, b;
  switch (type) {
    case 4: a = 0.0; b = 1.0; break;
    case 5: a = 0.5; b = 0.5; break;
    case 6: a = 0.0; b = 0.0; break;
    case 8: a = 1.0 / 3.0; b = 1.0 / 3.0; break;
    default: a = 3.0 / 8.0; b = 3.0 / 8.0; break;  // type 9
  }
  const double nppm = a + p * (nd + 1.0 - a - b);
  const double j = std::floor(nppm + fuzz);
  double h = nppm - j;
  if (std::fabs(h) < fuzz) h = 0.0;
  return mix(j, h);
}

Vec sorted_copy(const Vec& x, const char* fn) {
  require_finite_input(x, fn);
  Vec s(x);
  std::sort(s.begin(), s.end());
  return s;
}

}  // namespace

double sum(const Vec& x) {
  require_finite_input(x, "sum");
  return static_cast<double>(ksum(x));
}

double mean(const Vec& x) {
  require_finite_input(x, "mean");
  if (x.empty()) return kNaN;
  return static_cast<double>(mean_dd(x));
}

double var(const Vec& x) {
  require_finite_input(x, "var");
  const size_t n = x.size();
  if (n < 2) return kNaN;
  const DD m = mean_dd(x);
  DD ss;
  for (double v : x) {
    const DD d = DD(v) - m;
    ss += d * d;
  }
  return static_cast<double>(ss / DD(static_cast<double>(n - 1)));
}

double sd(const Vec& x) { return std::sqrt(var(x)); }

double median(const Vec& x) {
  require_finite_input(x, "median");
  const size_t n = x.size();
  if (n == 0) return kNaN;
  Vec s(x);
  const size_t half = (n + 1) / 2;  // 1-based
  std::nth_element(s.begin(), s.begin() + (half - 1), s.end());
  const double a = s[half - 1];
  if (n % 2 == 1) return a;
  const double b = *std::min_element(s.begin() + half, s.end());
  return static_cast<double>((DD(a) + DD(b)) / DD(2.0));
}

double mad(const Vec& x, double constant) {
  require_finite_input(x, "mad");
  if (x.empty()) return kNaN;
  const double c = median(x);
  Vec d(x.size());
  for (size_t i = 0; i < x.size(); ++i) d[i] = std::fabs(x[i] - c);
  return constant * median(d);
}

double cov(const Vec& x, const Vec& y) {
  if (x.size() != y.size()) throw Error("cov: x and y have different lengths");
  require_finite_input(x, "cov");
  require_finite_input(y, "cov");
  const size_t n = x.size();
  if (n < 2) return kNaN;
  const DD mx = mean_dd(x), my = mean_dd(y);
  DD s;
  for (size_t i = 0; i < n; ++i) s += (DD(x[i]) - mx) * (DD(y[i]) - my);
  return static_cast<double>(s / DD(static_cast<double>(n - 1)));
}

double cor(const Vec& x, const Vec& y, const std::string& method) {
  if (x.size() != y.size()) throw Error("cor: x and y have different lengths");
  require_finite_input(x, "cor");
  require_finite_input(y, "cor");
  if (method == "pearson") return pearson(x, y);
  if (method == "spearman") return pearson(rank(x), rank(y));
  if (method == "kendall") return kendall(x, y);
  throw Error("cor: unknown method '" + method + "' (use pearson, spearman or kendall)");
}

Vec rank(const Vec& x) {
  require_finite_input(x, "rank");
  const size_t n = x.size();
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  Vec r(n);
  for (size_t i = 0; i < n;) {
    size_t j = i;
    while (j < n && x[idx[j]] == x[idx[i]]) ++j;
    const double avg = (static_cast<double>(i + 1) + static_cast<double>(j)) / 2.0;
    for (size_t k = i; k < j; ++k) r[idx[k]] = avg;
    i = j;
  }
  return r;
}

double quantile(const Vec& x, double prob, int type) { return quantiles(x, Vec{prob}, type)[0]; }

Vec quantiles(const Vec& x, const Vec& probs, int type) {
  if (type < 1 || type > 9) throw Error("quantile: type must be an integer between 1 and 9");
  const double pfuzz = 2e-14;  // R accepts probabilities up to 2e-14 outside [0, 1]
  Vec p(probs);
  for (double& v : p) {
    if (std::isnan(v)) continue;
    if (v < -pfuzz || v > 1.0 + pfuzz) throw Error("quantile: 'probs' outside [0,1]");
    v = std::clamp(v, 0.0, 1.0);
  }
  const Vec s = sorted_copy(x, "quantile");
  Vec out(p.size());
  for (size_t i = 0; i < p.size(); ++i) out[i] = quantile_sorted(s, p[i], type);
  return out;
}

}  // namespace rcli::summary
