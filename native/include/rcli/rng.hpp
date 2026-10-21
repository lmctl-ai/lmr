// Seeded random numbers. Streams are NOT R-identical (an accepted decision): the generator is xoshiro256**
// seeded through splitmix64; variates use published exact methods (inversion for the normal via qnorm,
// Marsaglia-Tsang for gamma, BTPE or inversion for binomial, PTRS or inversion for Poisson, ...).
// Guarantee: the same seed and the same call sequence give the same values on every platform (no libm-dependent
// rejection paths that differ by ulps are allowed to change the accepted/rejected decision in practice; document any).
#pragma once
#include <cstdint>
#include "rcli/common.hpp"
#include "rcli/dist.hpp"

namespace rcli {

class Rng {
 public:
  explicit Rng(uint64_t seed);
  uint64_t next_u64();
  double unif();                       // uniform on the open interval (0, 1)
  double unif(double a, double b);
  double norm();                       // standard normal
  double norm(double mean, double sd) { return mean + sd * norm(); }
  double exp();                        // standard exponential
  double gamma(double shape);          // shape > 0, scale 1
  double draw(const Dist& d);          // one variate from any Dist
  Vec draw(const Dist& d, size_t n);
  uint64_t below(uint64_t n);          // uniform integer in [0, n), unbiased
  // Indices in [0, n): `size` of them, with or without replacement, optionally weighted (weights need not sum to 1).
  std::vector<size_t> sample(size_t n, size_t size, bool replace, const Vec* prob = nullptr);

 private:
  uint64_t s_[4];
};

}  // namespace rcli
