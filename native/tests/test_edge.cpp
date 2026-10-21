// Edge cases across modules: degenerate inputs must give a value or an rcli::Error (never a crash, hang or UB),
// and callback-driven algorithms must terminate in bounded evaluations for hostile callbacks.
// Run under -fsanitize=address,undefined as well (cmake -DLMR_SANITIZE=address,undefined).
#include <functional>
#include <string>

#include "check.hpp"
#include "rcli/cluster.hpp"
#include "rcli/dist.hpp"
#include "rcli/fft.hpp"
#include "rcli/htest.hpp"
#include "rcli/interp.hpp"
#include "rcli/models.hpp"
#include "rcli/optim.hpp"
#include "rcli/rng.hpp"
#include "rcli/smooth.hpp"
#include "rcli/summary.hpp"
#include "rcli/ts.hpp"

using namespace rcli;

namespace {

// Either outcome is acceptable; anything else (another exception type, a crash) fails the test run.
template <class F>
void value_or_error(F&& f, const std::string& what) {
  try {
    f();
    check::ok(true, what);
  } catch (const Error&) {
    check::ok(true, what);
  } catch (const std::exception& e) {
    check::ok(false, what + ": non-rcli exception " + e.what());
  }
}

struct Counter {
  int n = 0;
  int limit;
  explicit Counter(int l) : limit(l) {}
  void tick(const std::string& what) {
    if (++n > limit) throw std::runtime_error(what + ": evaluation budget exceeded (runaway loop)");
  }
};

GlmOptions poisson_opts() { GlmOptions o; o.family = "poisson"; return o; }
GlmOptions binomial_opts() { GlmOptions o; o.family = "binomial"; return o; }

}  // namespace

int main() {
  const double nan = kNaN, inf = kInf;

  // ---- hostile callbacks: bounded evaluations or a clean error ----
  const std::vector<std::pair<std::string, std::function<double(double)>>> hostile = {
      {"NaN", [nan](double) { return nan; }},
      {"+Inf", [inf](double) { return inf; }},
      {"-Inf", [inf](double) { return -inf; }},
      {"constant", [](double) { return 1.0; }},
      {"huge", [](double x) { return 1e308 * (x > 0 ? 1 : -1); }},
      {"alternating", [](double x) { return std::fmod(std::floor(x * 1e6), 2.0) == 0 ? 1.0 : -1.0; }},
  };
  for (const auto& h : hostile) {
    for (const char* method : {"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"}) {
      Counter c(200000);
      Rng rng(1);
      OptimOptions o;
      o.method = method;
      o.rng = &rng;
      o.maxit = 200;
      value_or_error([&] { optim([&](const Vec& p) { c.tick("optim"); return h.second(p[0] + p[1]); }, {0.5, 0.5}, o); },
                     std::string("optim ") + method + " with a " + h.first + " objective terminates");
    }
    {
      Counter c(100000);
      value_or_error([&] { optimize([&](double x) { c.tick("optimize"); return h.second(x); }, 0.0, 1.0); },
                     "optimize with a " + h.first + " objective terminates");
    }
    for (const char* ext : {"no", "yes", "upX", "downX"}) {
      Counter c(100000);
      RootOptions ro;
      ro.extend_int = ext;
      value_or_error([&] { uniroot([&](double x) { c.tick("uniroot"); return h.second(x); }, 0.0, 1.0, ro); },
                     std::string("uniroot extend=") + ext + " with a " + h.first + " function terminates");
    }
    for (auto range : {std::pair<double, double>{0, 1}, {-inf, inf}, {0, inf}, {-inf, 0}}) {
      Counter c(200000);
      value_or_error([&] { integrate([&](double x) { c.tick("integrate"); return h.second(x); }, range.first, range.second); },
                     "integrate with a " + h.first + " integrand terminates");
    }
    {
      Counter c(100000);
      Vec y = {1, 2, 3, 4, 5};
      for (const char* alg : {"default", "port"}) {
        NlsOptions no;
        no.algorithm = alg;
        value_or_error([&] {
          nls([&](const Vec& p, Vec& out) { c.tick("nls"); for (size_t i = 0; i < out.size(); ++i) out[i] = h.second(p[0] + double(i)); },
              y, {"a"}, {1.0}, no);
        }, std::string("nls ") + alg + " with a " + h.first + " model terminates");
      }
    }
  }
  // An exception thrown from the callback reaches the caller unchanged at every evaluation index.
  for (int k = 1; k <= 12; ++k) {
    int seen = 0;
    struct Stop {};
    auto f1 = [&](double x) { if (++seen == k) throw Stop{}; return x * x - 0.5; };
    bool caught = false;
    try { integrate(f1, 0.0, 1.0); } catch (const Stop&) { caught = true; }
    seen = 0;
    try { uniroot(f1, 0.0, 1.0); } catch (const Stop&) { caught = caught && true; }
    check::ok(caught, "a callback exception at evaluation " + std::to_string(k) + " propagates");
  }

  // ---- degenerate data ----
  value_or_error([] { summary::mean({}); }, "mean of an empty vector");
  value_or_error([] { summary::var({1.0}); }, "var of one value");
  value_or_error([] { summary::quantile({}, 0.5); }, "quantile of an empty vector");
  value_or_error([] { summary::quantile({1, 2}, 1.5); }, "quantile with p > 1");
  value_or_error([] { summary::cor({1, 1, 1}, {1, 2, 3}); }, "cor with zero variance");
  value_or_error([] { summary::cor({1, 2}, {1, 2, 3}); }, "cor with unequal lengths");
  check::ok(fft({}).empty(), "fft of nothing is nothing");
  value_or_error([] { CubicSpline({1}, {1}); }, "spline with one point");
  value_or_error([] { CubicSpline({1, 1, 1}, {1, 2, 3}); }, "spline with one distinct x");
  value_or_error([nan] { CubicSpline(Vec{0, 1, 2}, Vec{0, nan, 2}); }, "spline with NaN y");
  value_or_error([] { Approx({}, {}); }, "approx with no points");
  value_or_error([] { smooth_spline({1, 2, 3}, {1, 2, 3}); }, "smooth.spline with 3 points");
  value_or_error([] { loess({1, 2, 3}, {1, 2, 3}, LoessOptions{0.1, 2, "gaussian"}); }, "loess with a tiny span");
  value_or_error([] { loess({1, 1, 1, 1, 1}, {1, 2, 3, 4, 5}); }, "loess with constant x");
  value_or_error([] { t_test({1.0}); }, "t-test with one value");
  value_or_error([] { t_test({2, 2, 2}); }, "t-test with zero variance");
  value_or_error([] { t_test(Vec{1, 2, 3}, Vec{}); }, "two-sample t with an empty y");
  value_or_error([] { wilcox_test({}); }, "wilcox with nothing");
  value_or_error([] { shapiro_test({1, 1, 1, 1}); }, "shapiro with identical values");
  value_or_error([] { shapiro_test({1, 2}); }, "shapiro with 2 values");
  value_or_error([] { var_test({1, 2}, {3, 3}); }, "var.test with zero variance");
  value_or_error([] { cor_test({1, 2, 3}, {1, 1, 1}); }, "cor.test with constant y");
  value_or_error([] { chisq_test(Matrix(2, 2, 0.0)); }, "chisq of an all-zero table");
  value_or_error([] { fisher_test(Matrix(2, 2, 0.0)); }, "fisher of an all-zero table");
  value_or_error([] { binom_test(5, 0); }, "binom.test with n = 0");
  value_or_error([] { binom_test(5, 3); }, "binom.test with x > n");
  value_or_error([] { prop_test({0}, {0}); }, "prop.test with no trials");
  value_or_error([] { ks_test({}, Dist::norm(0, 1)); }, "ks with no data");
  {
    Data d;
    d.add("y", {1, 1, 1});
    d.add("x", {1, 2, 3});
    value_or_error([&] { lm(Formula::parse("y ~ x"), d); }, "lm with a constant response");
    value_or_error([&] { glm(Formula::parse("y ~ x"), d, poisson_opts()); }, "poisson glm, all equal counts");
    Data z;
    z.add("y", {0, 0, 0, 0});
    z.add("x", {1, 2, 3, 4});
    value_or_error([&] { glm(Formula::parse("y ~ x"), z, poisson_opts()); }, "poisson glm, all zero counts");
    Data sep;
    sep.add("y", {0, 0, 0, 1, 1, 1});
    sep.add("x", {1, 2, 3, 4, 5, 6});
    value_or_error([&] { glm(Formula::parse("y ~ x"), sep, binomial_opts()); }, "binomial glm, perfect separation");
    Data wide;
    wide.add("y", {1, 2});
    wide.add("a", {1, 3});
    wide.add("b", {2, 1});
    wide.add("c", {5, 4});
    value_or_error([&] { lm(Formula::parse("y ~ a + b + c"), wide); }, "lm with more columns than rows");
    value_or_error([&] { Formula::parse("y ~ (a"); }, "unbalanced formula");
    value_or_error([&] { Formula::parse("~ a"); }, "formula without response");
    value_or_error([&] { lm(Formula::parse("y ~ q"), wide); }, "unknown column");
    Data empty;
    empty.add("y", {});
    empty.add("x", {});
    value_or_error([&] { lm(Formula::parse("y ~ x"), empty); }, "lm on no rows");
  }
  for (size_t n : {size_t(1), size_t(2), size_t(3), size_t(5)}) {
    Vec x(n, 1.0);
    value_or_error([&] { acf(x); }, "acf of a constant series, n = " + std::to_string(n));
    value_or_error([&] { pacf(x); }, "pacf of a constant series");
    ArimaSpec s;
    s.p = 1;
    value_or_error([&] { arima(x, s); }, "arima on a short constant series");
  }
  {
    Vec x(30);
    for (size_t i = 0; i < x.size(); ++i) x[i] = std::sin(double(i));
    ArimaSpec big;
    big.p = 10;
    big.q = 10;
    value_or_error([&] { arima(x, big); }, "arima with more parameters than the data supports");
    value_or_error([&] { stl(x, 20); }, "stl with fewer than two periods");
    Vec withnan = x;
    withnan[5] = nan;
    value_or_error([&] { stl(withnan, 6); }, "stl with NaN");
    ArimaSpec s;
    s.p = 1;
    value_or_error([&] { arima(x, s).forecast(1000); }, "forecast 1000 steps");
  }
  {
    Rng rng(1);
    Matrix one(3, 1);
    one(0, 0) = 1; one(1, 0) = 1; one(2, 0) = 1;
    value_or_error([&] { kmeans(one, 2, rng); }, "kmeans with identical rows");
    value_or_error([&] { kmeans(one, 3, rng); }, "kmeans with k = n");
    value_or_error([&] { hclust(Vec{0.0}, 2, "average"); }, "hclust with n = 2");
    value_or_error([&] { hclust(Vec{1, 1, 1}, 3, "ward.D2"); }, "hclust with tied distances");
    value_or_error([&] { hclust(Vec{inf, 1, 1}, 3, "complete"); }, "hclust with an infinite distance");
    value_or_error([&] { hclust(Vec{1, 1}, 3, "single"); }, "hclust with a short distance vector");
  }
  {
    // Distributions at boundaries: value or error, never a hang.
    for (double df : {1e-300, 1e-10, 1e10, 1e300}) {
      value_or_error([&] { Dist::t(df).q(0.975); }, "qt at extreme df");
      value_or_error([&] { Dist::chisq(df).q(0.5); }, "qchisq at extreme df");
    }
    value_or_error([] { Dist::binom(1e15, 0.5).q(0.5); }, "qbinom with a huge size");
    value_or_error([] { Dist::pois(1e15).q(0.5); }, "qpois with a huge lambda");
    value_or_error([] { Dist::nbinom(1e-10, 0.5).q(0.999); }, "qnbinom with a tiny size");
    value_or_error([] { Dist::beta(1e-300, 1e-300).p(0.5); }, "pbeta with tiny shapes");
    value_or_error([] { Dist::gamma(1e300, 1).q(0.5); }, "qgamma with a huge shape");
    value_or_error([nan] { Dist::norm(0, 1).q(nan); }, "qnorm(NaN)");
    value_or_error([] { Dist::norm(0, 1).q(1.5); }, "qnorm(1.5)");
    Rng rng(2);
    value_or_error([&] { rng.draw(Dist::binom(1e15, 0.3)); }, "rbinom with a huge size");
    value_or_error([&] { rng.draw(Dist::pois(1e15)); }, "rpois with a huge lambda");
    value_or_error([&] { rng.sample(3, 5, false); }, "sample larger than the population without replacement");
    value_or_error([&] { Vec p = {0, 0, 0}; rng.sample(3, 2, true, &p); }, "sample with all-zero weights");
  }
  // Degenerate inputs that R answers with a specific value or error (black-box Rscript, R 4.3.3).
  auto message = [](std::function<void()> f) -> std::string {
    try { f(); } catch (const Error& e) { return e.what(); }
    return "";
  };
  check::ok(message([] { binom_test(5, 3); }) == "binom.test: 'n' must be a positive integer >= 'x'", "binom.test(5, 3) message");
  check::ok(message([] { binom_test(-1, 3); }) == "binom.test: 'x' must be nonnegative and integer", "binom.test(-1, 3) message");
  check::ok(message([] { binom_test(2.5, 3); }) == "binom.test: 'x' must be nonnegative and integer", "binom.test(2.5, 3) message");
  check::ok(message([] { binom_test(0, 0); }) == "binom.test: 'n' must be a positive integer >= 'x'", "binom.test(0, 0) message");
  check::ok(message([] { t_test(Vec{2, 2, 2}); }) == "t.test: data are essentially constant", "t.test(c(2, 2, 2))");
  check::ok(message([] { shapiro_test({1, 1, 1, 1}); }) == "shapiro.test: all 'x' values are identical", "shapiro.test(rep(1, 4))");
  check::ok(message([] { chisq_test(Matrix(2, 2, 0.0)); }) == "chisq.test: at least one entry of 'x' must be positive", "chisq.test(matrix(0, 2, 2))");
  check::ok(message([] { summary::quantile({1, 2}, 1.5); }) == "quantile: 'probs' outside [0,1]", "quantile(1:2, 1.5)");
  check::close(var_test({1, 2}, {3, 3}).p_value, 0.0, 0.0, "var.test(1:2, c(3, 3))$p.value is 0");
  check::close(fisher_test(Matrix(2, 2, 0.0)).p_value, 1.0, 0.0, "fisher.test(matrix(0, 2, 2))$p.value is 1");
  check::ok(std::isnan(summary::var({1})) && std::isnan(summary::mean({})) && std::isnan(summary::cor({1, 1, 1}, {1, 2, 3})),
            "var(1), mean(numeric(0)), cor with constant x are NA/NaN");
  check::ok(std::isnan(Dist::norm(0, 1).q(1.5)) && Dist::t(1e-300).q(0.975) == kInf, "qnorm(1.5) NaN, qt(.975, 1e-300) Inf");
  check::close(Dist::binom(1e15, 0.5).q(0.5), 5e14, 0.0, "qbinom(.5, 1e15, .5)");
  {
    Data z;
    z.add("y", {0, 0, 0, 0});
    z.add("x", {1, 2, 3, 4});
    GlmOptions o;
    o.family = "poisson";
    check::close(glm(Formula::parse("y ~ x"), z, o).coef.estimate[0], -24.30258509, 1e-9,
                 "coef(glm(c(0,0,0,0) ~ I(1:4), poisson))[1]");
  }
  // Upper-tail discrete quantiles with p near 1 (found by tests/oracle/fuzz_r.py): the fuzz must not push the target
  // to 1. R values: qbinom(0.9999999999999987, 533, 0.8853877112555562, lower.tail = FALSE) = 407, and the same rule
  // for qhyper(1.923934826718089e-15, 132, 71, 155, lower.tail = FALSE) gives 124 by R's own phyper.
  check::close(Dist::binom(533, 0.8853877112555562).q(0.9999999999999987, false), 407, 0, "qbinom upper tail near 1");
  check::close(Dist::binom(533, 0.885).q(1 - std::ldexp(1.0, -30), false), 423, 0, "qbinom(1 - 2^-30, upper)");
  check::close(Dist::binom(533, 0.885).q(0.7, false), 468, 0, "qbinom(0.7, upper)");
  check::close(Dist::hyper(132, 71, 155).q(1.923934826718089e-15, false), 124, 0, "qhyper tiny upper tail");
  check::close(Dist::binom(533, 0.885).q(std::log(0.9999999999999987), false, true), 407, 0, "qbinom upper, log scale");
  return check::report("edge");
}
