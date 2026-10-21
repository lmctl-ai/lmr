// Probability distributions with R's d/p/q semantics (dnorm, pnorm(lower.tail, log.p), qnorm, ...).
// Parameterisations follow R exactly: gamma(shape, rate), weibull(shape, scale), nbinom(size, prob), hyper(m, n, k), etc.
// Constructors validate parameters and throw rcli::Error (R would return NaN with a warning).
// Accuracy target: d/p 1e-13 relative (including far tails via log/upper tails), q 1e-10 relative.
#pragma once
#include "rcli/common.hpp"

namespace rcli {

enum class DistKind {
  Norm, Lnorm, T, Chisq, F, Gamma, Beta, Exp, Unif, Weibull, Cauchy, Logis,  // continuous
  Binom, Pois, Geom, Nbinom, Hyper                                             // discrete
};

class Dist {
 public:
  static Dist norm(double mean, double sd);
  static Dist lnorm(double meanlog, double sdlog);
  static Dist t(double df);
  static Dist chisq(double df);
  static Dist f(double df1, double df2);
  static Dist gamma(double shape, double rate);
  static Dist beta(double shape1, double shape2);
  static Dist exp(double rate);
  static Dist unif(double min, double max);
  static Dist weibull(double shape, double scale);
  static Dist cauchy(double location, double scale);
  static Dist logis(double location, double scale);
  static Dist binom(double size, double prob);
  static Dist pois(double lambda);
  static Dist geom(double prob);
  static Dist nbinom(double size, double prob);
  static Dist hyper(double m, double n, double k);  // m white, n black, k drawn (R dhyper(x, m, n, k))
  // By R name ("norm", "t", "chisq", "f", "gamma", "beta", "exp", "unif", "lnorm", "weibull", "cauchy", "logis",
  // "binom", "pois", "geom", "nbinom", "hyper") with up to three parameters in R's positional order.
  static Dist by_name(const std::string& name, double p1, double p2 = kNaN, double p3 = kNaN);

  DistKind kind() const { return kind_; }
  std::string name() const;                 // R name, e.g. "norm"
  bool discrete() const;
  double param(int i) const { return par_[i]; }

  double d(double x, bool give_log = false) const;
  double p(double q, bool lower_tail = true, bool log_p = false) const;
  double q(double p, bool lower_tail = true, bool log_p = false) const;  // discrete: smallest x with P(X <= x) >= p
  double mean() const;   // NaN or Inf where undefined, as in theory
  double var() const;

 private:
  Dist(DistKind k, double a, double b, double c) : kind_(k), par_{a, b, c} {}
  DistKind kind_;
  double par_[3];
};

}  // namespace rcli
