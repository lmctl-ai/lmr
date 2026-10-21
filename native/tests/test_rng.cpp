// Tests for rcli::Rng: determinism (recorded first values), moments of large samples within 4 standard errors,
// discrete frequencies against the pmf, sample() properties.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>

#include "check.hpp"
#include "rcli/rng.hpp"

using rcli::Dist;
using rcli::Rng;
using rcli::Vec;

namespace {

// Mean and variance of n draws must be within 4 standard errors of the theory (variance SE via the 4th moment
// estimated from the sample).
void moments(Rng& r, const Dist& d, size_t n, const std::string& what) {
  Vec x = r.draw(d, n);
  double m = 0;
  for (double v : x) m += v;
  m /= n;
  double s2 = 0, m4 = 0;
  for (double v : x) {
    const double e = v - m;
    s2 += e * e;
    m4 += e * e * e * e;
  }
  s2 /= (n - 1);
  m4 /= n;
  const double mu = d.mean(), var = d.var();
  check::ok(std::fabs(m - mu) <= 4 * std::sqrt(var / n), what + " mean " + std::to_string(m) + " vs " + std::to_string(mu));
  const double se_var = std::sqrt(std::max(m4 - s2 * s2, 0.0) / n);
  check::ok(std::fabs(s2 - var) <= 4 * se_var + 1e-12, what + " var " + std::to_string(s2) + " vs " + std::to_string(var));
}

// Chi-square goodness of fit of a discrete distribution on its cells [lo, hi] (tails pooled); require the statistic
// below the 0.9999 quantile.
void chisq_fit(Rng& r, const Dist& d, size_t n, double lo, double hi, const std::string& what) {
  const int K = static_cast<int>(hi - lo) + 1;
  std::vector<double> obs(K + 2, 0.0);
  bool integers = true;
  for (size_t i = 0; i < n; ++i) {
    const double v = r.draw(d);
    integers = integers && v == std::floor(v);
    if (v < lo) obs[0] += 1;
    else if (v > hi) obs[K + 1] += 1;
    else obs[1 + static_cast<int>(v - lo)] += 1;
  }
  check::ok(integers, what + " integer values");
  double stat = 0;
  int cells = 0;
  for (int c = 0; c < K + 2; ++c) {
    double p;
    if (c == 0) p = lo > 0 ? d.p(lo - 1) : 0;
    else if (c == K + 1) p = d.p(hi, false);
    else p = d.d(lo + c - 1);
    const double e = p * n;
    if (e < 5) continue;  // tiny cells are skipped (their counts are checked for zero-probability below)
    stat += (obs[c] - e) * (obs[c] - e) / e;
    ++cells;
  }
  const double crit = Dist::chisq(std::max(cells - 1, 1)).q(0.9999);
  check::ok(stat < crit, what + " chisq " + std::to_string(stat) + " < " + std::to_string(crit));
}

}  // namespace

int main() {
  // Determinism: recorded first values of xoshiro256** seeded through splitmix64(42).
  {
    Rng r(42);
    check::ok(r.next_u64() == 1546998764402558742ULL, "u64[0] seed 42");
    check::ok(r.next_u64() == 6990951692964543102ULL, "u64[1] seed 42");
    check::ok(r.next_u64() == 12544586762248559009ULL, "u64[2] seed 42");
    Rng s(42);
    check::close(s.unif(), 0.083862971059882219, 0, "unif[0] seed 42");
    check::close(s.unif(), 0.37898025066266866, 0, "unif[1] seed 42");
    Rng t(1);
    check::close(t.norm(), 0.53282267825248775, 1e-15, "norm seed 1");
    check::close(t.exp(), 0.65308716599008509, 1e-15, "exp seed 1");
    check::close(t.gamma(2.5), 2.4534830787154482, 1e-15, "gamma seed 1");
    check::ok(t.below(10) == 6, "below seed 1");
    // Same seed, same stream; different seeds differ.
    Rng a(123), b(123), c(124);
    bool same = true, diff = false;
    for (int i = 0; i < 1000; ++i) {
      const uint64_t x = a.next_u64();
      same = same && x == b.next_u64();
      diff = diff || x != c.next_u64();
    }
    check::ok(same, "same seed same stream");
    check::ok(diff, "different seed different stream");
  }
  // unif in the open interval.
  {
    Rng r(7);
    bool inside = true;
    for (int i = 0; i < 100000; ++i) {
      const double u = r.unif();
      inside = inside && u > 0 && u < 1;
    }
    check::ok(inside, "unif in (0,1)");
    check::throws([&] { r.unif(2, 1); }, "unif bad range");
    check::throws([&] { r.gamma(-1); }, "gamma bad shape");
    check::throws([&] { r.below(0); }, "below(0)");
  }
  // Moments of continuous families.
  {
    Rng r(2024);
    const size_t N = 200000;
    moments(r, Dist::norm(1, 2), N, "norm");
    moments(r, Dist::exp(3), N, "exp");
    moments(r, Dist::unif(-1, 5), N, "unif");
    moments(r, Dist::gamma(0.3, 2), N, "gamma 0.3");
    moments(r, Dist::gamma(1, 1), N, "gamma 1");
    moments(r, Dist::gamma(7.5, 0.5), N, "gamma 7.5");
    moments(r, Dist::beta(0.5, 0.5), N, "beta .5 .5");
    moments(r, Dist::beta(2, 5), N, "beta 2 5");
    moments(r, Dist::beta(0.05, 3), N, "beta .05 3");
    moments(r, Dist::chisq(4), N, "chisq 4");
    moments(r, Dist::chisq(0.7), N, "chisq 0.7");
    moments(r, Dist::t(10), N, "t 10");
    moments(r, Dist::f(5, 12), N, "f 5 12");
    moments(r, Dist::lnorm(0, 0.5), N, "lnorm");
    moments(r, Dist::weibull(1.5, 2), N, "weibull");
    moments(r, Dist::logis(1, 2), N, "logis");
    // Cauchy: median and quartiles instead of moments.
    Vec x = r.draw(Dist::cauchy(1, 2), N);
    std::sort(x.begin(), x.end());
    check::ok(std::fabs(x[N / 2] - 1) < 0.05, "cauchy median");
    check::ok(std::fabs(x[3 * N / 4] - 3) < 0.1, "cauchy upper quartile");
    // Normal tail frequency: P(Z > 3) = 0.0013499
    size_t tail = 0;
    for (int i = 0; i < 400000; ++i) tail += r.norm() > 3;
    const double ph = tail / 400000.0, p0 = 0.0013498980316301;
    check::ok(std::fabs(ph - p0) < 4 * std::sqrt(p0 * (1 - p0) / 400000), "norm tail frequency");
  }
  // Discrete families: moments and chi-square goodness of fit.
  {
    Rng r(99);
    const size_t N = 200000;
    moments(r, Dist::binom(20, 0.3), N, "binom inversion");
    moments(r, Dist::binom(1000, 0.4), N, "binom BTRS");
    moments(r, Dist::binom(1000, 0.97), N, "binom BTRS flipped");
    moments(r, Dist::pois(3.5), N, "pois inversion");
    moments(r, Dist::pois(250), N, "pois PTRS");
    moments(r, Dist::geom(0.2), N, "geom");
    moments(r, Dist::nbinom(3.5, 0.4), N, "nbinom");
    moments(r, Dist::hyper(30, 50, 20), N, "hyper");
    chisq_fit(r, Dist::binom(20, 0.3), N, 0, 14, "binom 20 .3");
    chisq_fit(r, Dist::binom(200, 0.5), N, 70, 130, "binom 200 .5 (BTRS)");
    chisq_fit(r, Dist::binom(5000, 0.02), N, 70, 130, "binom 5000 .02 (BTRS)");
    chisq_fit(r, Dist::pois(4), N, 0, 15, "pois 4");
    chisq_fit(r, Dist::pois(30), N, 10, 55, "pois 30 (PTRS)");
    chisq_fit(r, Dist::pois(1000), N, 900, 1100, "pois 1000 (PTRS)");
    chisq_fit(r, Dist::geom(0.3), N, 0, 20, "geom");
    chisq_fit(r, Dist::nbinom(2, 0.3), N, 0, 25, "nbinom");
    chisq_fit(r, Dist::hyper(10, 7, 8), N, 1, 8, "hyper small");
    chisq_fit(r, Dist::hyper(500, 800, 300), N, 90, 140, "hyper large");
    // Degenerate parameters.
    check::ok(r.draw(Dist::binom(10, 0)) == 0 && r.draw(Dist::binom(10, 1)) == 10, "binom degenerate");
    check::ok(r.draw(Dist::pois(0)) == 0, "pois 0");
    check::ok(r.draw(Dist::geom(1)) == 0, "geom 1");
    check::ok(r.draw(Dist::hyper(5, 0, 3)) == 3, "hyper degenerate");
  }
  // below(): uniform frequencies, range.
  {
    Rng r(5);
    std::vector<double> cnt(7, 0);
    bool range = true;
    for (int i = 0; i < 70000; ++i) {
      const uint64_t v = r.below(7);
      range = range && v < 7;
      if (v < 7) cnt[v] += 1;
    }
    check::ok(range, "below range");
    double stat = 0;
    for (double c : cnt) stat += (c - 10000) * (c - 10000) / 10000;
    check::ok(stat < Dist::chisq(6).q(0.9999), "below uniform");
    check::ok(r.below(1) == 0, "below(1)");
    const uint64_t big = 0xFFFFFFFFFFFFFFF1ULL;
    check::ok(r.below(big) < big, "below(big)");
  }
  // sample()
  {
    Rng r(11);
    auto s = r.sample(100, 100, false);
    std::set<size_t> u(s.begin(), s.end());
    check::ok(u.size() == 100 && *u.rbegin() == 99, "sample permutation");
    auto s2 = r.sample(1000000, 50, false);
    std::set<size_t> u2(s2.begin(), s2.end());
    check::ok(u2.size() == 50, "sample distinct");
    check::throws([&] { r.sample(5, 6, false); }, "sample too large");
    // uniform without replacement: first position frequencies
    std::vector<double> first(5, 0);
    for (int i = 0; i < 50000; ++i) first[r.sample(5, 3, false)[0]] += 1;
    double stat = 0;
    for (double c : first) stat += (c - 10000) * (c - 10000) / 10000;
    check::ok(stat < Dist::chisq(4).q(0.9999), "sample uniform first position");
    // weighted with replacement (alias)
    Vec w = {1, 2, 3, 0, 4};
    auto ws = r.sample(5, 100000, true, &w);
    std::vector<double> c(5, 0);
    for (size_t v : ws) c[v] += 1;
    check::ok(c[3] == 0, "alias zero weight never drawn");
    stat = 0;
    for (int k = 0; k < 5; ++k) {
      if (w[k] == 0) continue;
      const double e = 100000 * w[k] / 10;
      stat += (c[k] - e) * (c[k] - e) / e;
    }
    check::ok(stat < Dist::chisq(3).q(0.9999), "alias frequencies");
    // weighted without replacement: distinct, zero weight excluded, first draw proportional to weights
    std::vector<double> f(5, 0);
    bool distinct = true, nozero = true;
    for (int i = 0; i < 40000; ++i) {
      auto t = r.sample(5, 4, false, &w);
      std::set<size_t> tu(t.begin(), t.end());
      distinct = distinct && tu.size() == 4;
      nozero = nozero && !tu.count(3);
      f[t[0]] += 1;
    }
    check::ok(distinct, "weighted no-replace distinct");
    check::ok(nozero, "weighted no-replace excludes zero weight");
    stat = 0;
    for (int k = 0; k < 5; ++k) {
      if (w[k] == 0) continue;
      const double e = 40000 * w[k] / 10;
      stat += (f[k] - e) * (f[k] - e) / e;
    }
    check::ok(stat < Dist::chisq(3).q(0.9999), "weighted no-replace first draw");
    check::throws([&] { r.sample(5, 5, false, &w); }, "too few positive probabilities");
    Vec bad = {1, -1};
    check::throws([&] { r.sample(2, 1, true, &bad); }, "negative weight");
    check::ok(r.sample(3, 0, false).empty(), "sample size 0");
  }
  return check::report("test_rng");
}
