// Seeded random numbers (clean-room; written from the published algorithms).
//
// - Generator: xoshiro256** (D. Blackman & S. Vigna, "Scrambled linear pseudorandom number generators", 2018,
//   public domain), state seeded by splitmix64 (G. Steele, D. Lea & C. Flood 2014; public domain) from the seed.
// - unif(): one 64-bit draw, top 53 bits, (k + 1/2) * 2^-53: open interval (0, 1), resolution 2^-53.
// - norm(): inversion, Dist::norm(0, 1).q(unif()) (Wichura AS 241), one 53-bit uniform per variate.
// - exp(): inversion, -log(unif()).
// - gamma(): G. Marsaglia & W. W. Tsang (2000), "A simple method for generating gamma variables", ACM TOMS 26,
//   with the boost Gamma(a) = Gamma(a + 1) U^(1/a) for a < 1.
// - Binomial: inversion (sequential search) when n min(p, 1-p) < 30, else BTRS, W. Hormann (1993), "The generation
//   of binomial random variates", J. Statist. Comput. Simul. 46.
// - Poisson: inversion when lambda < 30, else PTRS, W. Hormann (1993), "The transformed rejection method for
//   generating Poisson random variables", Insurance: Mathematics and Economics 12.
// - Geometric / negative binomial: gamma-Poisson mixture. Hypergeometric: inversion by chop-down search from the mode.
// - below(n): D. Lemire (2019), "Fast random integer generation in an interval", ACM TOMACS 29 (unbiased).
// - sample(): partial Fisher-Yates (Durstenfeld) without replacement; Walker's alias method (Vose's construction) for
//   weighted sampling with replacement; sequential weighted draws without replacement.
//
// Determinism: integer and IEEE double arithmetic (+, -, *, /, sqrt) are exact-rounded and identical on every
// platform. Values that pass through libm (log, exp, pow, the normal quantile) can differ in the last ulp between
// libm implementations; accept/reject decisions that involve libm (Marsaglia-Tsang squeeze, BTRS/PTRS final test)
// could in principle flip, with probability of order 1e-16 per draw. Build with -ffp-contract=off (implied by
// -std=c++17 under GCC; needed explicitly with -std=gnu++17 on FMA targets such as arm64) so that no fused
// multiply-adds change the rounding.
#include "rcli/rng.hpp"

#include <algorithm>
#include <cmath>

#include "rcli/special.hpp"

namespace rcli {

namespace {

inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

uint64_t splitmix64(uint64_t& x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// High and low 64 bits of a 64 x 64 bit product.
inline void mul128(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
  const uint64_t a0 = a & 0xFFFFFFFFULL, a1 = a >> 32, b0 = b & 0xFFFFFFFFULL, b1 = b >> 32;
  const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
  *lo = (mid << 32) | (p00 & 0xFFFFFFFFULL);
  *hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
}

[[noreturn]] void fail(const std::string& what) { throw Error("rng: " + what); }

// x^n for integer n >= 0 by binary powering (exact IEEE operations only).
double ipow(double x, double n) {
  double r = 1;
  while (n > 0) {
    const double h = std::floor(n / 2);
    if (n - 2 * h == 1) r *= x;
    x *= x;
    n = h;
  }
  return r;
}

}  // namespace

Rng::Rng(uint64_t seed) {
  uint64_t x = seed;
  for (auto& s : s_) s = splitmix64(x);
}

uint64_t Rng::next_u64() {
  const uint64_t result = rotl(s_[1] * 5, 7) * 9;
  const uint64_t t = s_[1] << 17;
  s_[2] ^= s_[0];
  s_[3] ^= s_[1];
  s_[1] ^= s_[2];
  s_[0] ^= s_[3];
  s_[2] ^= t;
  s_[3] = rotl(s_[3], 45);
  return result;
}

double Rng::unif() { return (static_cast<double>(next_u64() >> 11) + 0.5) * 0x1.0p-53; }

double Rng::unif(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b) || a > b) fail("unif: need finite min <= max");
  return a + (b - a) * unif();
}

double Rng::norm() {
  static const Dist std_norm = Dist::norm(0, 1);
  return std_norm.q(unif());
}

double Rng::exp() { return -std::log(unif()); }

double Rng::gamma(double a) {
  if (std::isnan(a) || a < 0 || std::isinf(a)) fail("gamma: shape must be finite and >= 0");
  if (a == 0) return 0;
  if (a < 1) {
    const double g = gamma(a + 1);
    return g * std::pow(unif(), 1 / a);
  }
  const double d = a - 1.0 / 3, c = 1 / std::sqrt(9 * d);
  for (;;) {
    double x, v;
    do {
      x = norm();
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;
    const double u = unif();
    const double x2 = x * x;
    if (u < 1 - 0.0331 * x2 * x2) return d * v;
    if (std::log(u) < 0.5 * x2 + d * (1 - v + std::log(v))) return d * v;
  }
}

namespace {

// log of a Gamma(a) variate, safe for tiny a (no underflow).
double log_gamma_variate(Rng& r, double a) {
  if (a < 1) return std::log(r.gamma(a + 1)) + std::log(r.unif()) / a;
  return std::log(r.gamma(a));
}

double rbinom(Rng& r, double n, double p) {
  if (n == 0 || p == 0) return 0;
  if (p == 1) return n;
  const bool flip = p > 0.5;
  const double pp = flip ? 1 - p : p, q = 1 - pp;
  double k;
  if (n * pp < 30) {
    const double f0 = ipow(q, n), ratio = pp / q;
    for (;;) {
      double u = r.unif(), f = f0, x = 0;
      bool done = false;
      while (x <= n) {
        if (u <= f) { done = true; break; }
        u -= f;
        x += 1;
        f *= (n - x + 1) / x * ratio;
      }
      if (done) { k = x; break; }
    }
  } else {
    // BTRS (Hormann 1993).
    const double spq = std::sqrt(n * pp * q);
    const double b = 1.15 + 2.53 * spq, a = -0.0873 + 0.0248 * b + 0.01 * pp, c = n * pp + 0.5;
    const double vr = 0.92 - 4.2 / b, alpha = (2.83 + 5.1 / b) * spq, lpq = std::log(pp / q);
    const double m = std::floor((n + 1) * pp);
    const double h = special::lgamma(m + 1) + special::lgamma(n - m + 1);
    for (;;) {
      const double u = r.unif() - 0.5;
      double v = r.unif();
      const double us = 0.5 - std::fabs(u);
      const double kk = std::floor((2 * a / us + b) * u + c);
      if (kk < 0 || kk > n) continue;
      if (us >= 0.07 && v <= vr) { k = kk; break; }
      v = std::log(v * alpha / (a / (us * us) + b));
      if (v <= h - special::lgamma(kk + 1) - special::lgamma(n - kk + 1) + (kk - m) * lpq) { k = kk; break; }
    }
  }
  return flip ? n - k : k;
}

double rpois(Rng& r, double lambda) {
  if (lambda == 0) return 0;
  if (lambda < 30) {
    const double f0 = std::exp(-lambda);
    for (;;) {
      double u = r.unif(), f = f0, x = 0;
      while (u > f && x < 1000) {
        u -= f;
        x += 1;
        f *= lambda / x;
      }
      if (x < 1000) return x;
    }
  }
  // PTRS (Hormann 1993).
  const double slam = std::sqrt(lambda), loglam = std::log(lambda);
  const double b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
  const double invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2);
  for (;;) {
    const double U = r.unif() - 0.5, V = r.unif();
    const double us = 0.5 - std::fabs(U);
    const double k = std::floor((2 * a / us + b) * U + lambda + 0.43);
    if (us >= 0.07 && V <= vr) return k;
    if (k < 0 || (us < 0.013 && V > us)) continue;
    if (std::log(V) + std::log(invalpha) - std::log(a / (us * us) + b) <= -lambda + k * loglam - special::lgamma(k + 1))
      return k;
  }
}

double rhyper(Rng& r, const Dist& d, double m, double n, double k) {
  const double L = std::max(0.0, k - n), H = std::min(k, m);
  if (L == H) return L;
  double mode = std::floor((k + 1) * (m + 1) / (m + n + 2));
  mode = std::min(std::max(mode, L), H);
  const double fm = d.d(mode);
  for (;;) {
    double u = r.unif() - fm;
    if (u <= 0) return mode;
    double xl = mode, xh = mode, fl = fm, fh = fm;
    while (xl > L || xh < H) {
      if (xl > L) {
        fl *= xl * (n - k + xl) / ((m - xl + 1) * (k - xl + 1));
        xl -= 1;
        u -= fl;
        if (u <= 0) return xl;
      }
      if (xh < H) {
        fh *= (m - xh) * (k - xh) / ((xh + 1) * (n - k + xh + 1));
        xh += 1;
        u -= fh;
        if (u <= 0) return xh;
      }
    }
    // rounding left a sliver of probability: draw again
  }
}

}  // namespace

double Rng::draw(const Dist& d) {
  const double a = d.param(0), b = d.param(1), c = d.param(2);
  switch (d.kind()) {
    case DistKind::Norm: return a + b * norm();
    case DistKind::Lnorm: return std::exp(a + b * norm());
    case DistKind::T: {
      if (std::isinf(a)) return norm();
      const double z = norm();
      return z / std::sqrt(2 * gamma(0.5 * a) / a);
    }
    case DistKind::Chisq: return 2 * gamma(0.5 * a);
    case DistKind::F: {
      const double x1 = std::isinf(a) ? 1.0 : 2 * gamma(0.5 * a) / a;
      const double x2 = std::isinf(b) ? 1.0 : 2 * gamma(0.5 * b) / b;
      return x1 / x2;
    }
    case DistKind::Gamma: return gamma(a) / b;
    case DistKind::Beta: {
      if (a == 0 && b == 0) return unif() < 0.5 ? 0.0 : 1.0;
      if (a == 0 || std::isinf(b)) return std::isinf(a) ? 0.5 : 0.0;
      if (b == 0 || std::isinf(a)) return 1.0;
      const double lx = log_gamma_variate(*this, a), ly = log_gamma_variate(*this, b);
      return 1 / (1 + std::exp(ly - lx));
    }
    case DistKind::Exp: return exp() / a;
    case DistKind::Unif: return a + (b - a) * unif();
    case DistKind::Weibull: return b * std::pow(exp(), 1 / a);
    case DistKind::Cauchy: {
      const double u = unif();
      return a + b * std::tan(3.14159265358979323846 * (u - 0.5));
    }
    case DistKind::Logis: {
      const double u = unif();
      return a + b * (std::log(u) - std::log1p(-u));
    }
    case DistKind::Binom: return rbinom(*this, a, b);
    case DistKind::Pois: return rpois(*this, a);
    case DistKind::Geom: return a == 1 ? 0.0 : rpois(*this, exp() * (1 - a) / a);
    case DistKind::Nbinom:
      if (a == 0 || b == 1) return 0;
      return rpois(*this, gamma(a) * (1 - b) / b);
    case DistKind::Hyper: return rhyper(*this, d, a, b, c);
  }
  return kNaN;
}

Vec Rng::draw(const Dist& d, size_t n) {
  Vec out(n);
  for (size_t i = 0; i < n; ++i) out[i] = draw(d);
  return out;
}

uint64_t Rng::below(uint64_t n) {
  if (n == 0) fail("below: n must be > 0");
  uint64_t hi, lo;
  mul128(next_u64(), n, &hi, &lo);
  if (lo < n) {
    const uint64_t t = (0 - n) % n;
    while (lo < t) mul128(next_u64(), n, &hi, &lo);
  }
  return hi;
}

std::vector<size_t> Rng::sample(size_t n, size_t size, bool replace, const Vec* prob) {
  std::vector<size_t> out;
  out.reserve(size);
  if (size == 0) return out;
  if (n == 0) fail("sample: cannot take a sample from an empty population");
  if (!prob) {
    if (replace) {
      for (size_t i = 0; i < size; ++i) out.push_back(static_cast<size_t>(below(n)));
      return out;
    }
    if (size > n) fail("sample: cannot take a sample larger than the population when replace = false");
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    for (size_t i = 0; i < size; ++i) {
      const size_t j = i + static_cast<size_t>(below(n - i));
      std::swap(idx[i], idx[j]);
      out.push_back(idx[i]);
    }
    return out;
  }
  const Vec& w = *prob;
  if (w.size() != n) fail("sample: prob must have length n");
  double total = 0;
  size_t npos = 0;
  for (double v : w) {
    if (!std::isfinite(v) || v < 0) fail("sample: prob must be finite and non-negative");
    total += v;
    if (v > 0) ++npos;
  }
  if (npos == 0) fail("sample: no positive probabilities");
  if (replace) {
    // Walker alias table, Vose's construction.
    std::vector<double> pr(n);
    std::vector<size_t> alias(n, 0), small, large;
    for (size_t i = 0; i < n; ++i) {
      pr[i] = w[i] * static_cast<double>(n) / total;
      (pr[i] < 1 ? small : large).push_back(i);
    }
    while (!small.empty() && !large.empty()) {
      const size_t s = small.back(), l = large.back();
      small.pop_back();
      alias[s] = l;
      pr[l] = (pr[l] + pr[s]) - 1;
      if (pr[l] < 1) {
        large.pop_back();
        small.push_back(l);
      }
    }
    for (size_t i : large) pr[i] = 1;
    for (size_t i : small) pr[i] = 1;  // rounding leftovers
    for (size_t i = 0; i < size; ++i) {
      const size_t j = static_cast<size_t>(below(n));
      out.push_back(unif() < pr[j] ? j : alias[j]);
    }
    return out;
  }
  if (size > npos) fail("sample: too few positive probabilities for sampling without replacement");
  Vec rem(w);
  for (size_t i = 0; i < size; ++i) {
    double t = 0;
    for (double v : rem) t += v;
    const double target = unif() * t;
    double acc = 0;
    size_t pick = n;
    for (size_t j = 0; j < n; ++j) {
      if (rem[j] <= 0) continue;
      acc += rem[j];
      pick = j;
      if (target < acc) break;
    }
    out.push_back(pick);
    rem[pick] = 0;
  }
  return out;
}

}  // namespace rcli
