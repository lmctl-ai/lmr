// Classical hypothesis tests with R's htest semantics (statistic, parameters, p-value, interval, estimates).
// Exact p-values use published recurrences (Wilcoxon / signed-rank distributions by counting, the exact Kolmogorov
// distribution of Marsaglia, Tsang & Wang 2003, Fisher's noncentral hypergeometric conditional MLE, Royston 1995
// AS R94 for Shapiro-Wilk). Defaults (exact vs normal approximation, continuity corrections) follow R's documented
// rules so that the same inputs pick the same branch.
#pragma once
#include <string>
#include <utility>
#include "rcli/common.hpp"

namespace rcli {

class Dist;

struct HtestResult {
  std::string method;                 // R's method string, e.g. "Welch Two Sample t-test"
  std::string alternative = "two.sided";
  std::string statistic_name;         // "t", "W", "V", "X-squared", "D", "p" ...; empty when the test has none
  double statistic = kNaN;
  std::vector<std::pair<std::string, double>> parameters;  // e.g. {"df", 17.8}; {"num df", 9}, {"denom df", 9}
  double p_value = kNaN;
  bool has_conf_int = false;
  double conf_low = kNaN, conf_high = kNaN, conf_level = kNaN;
  std::vector<std::pair<std::string, double>> estimates;   // e.g. {"mean of x", ..}, {"odds ratio", ..}
  std::string null_name;              // e.g. "difference in means"
  double null_value = kNaN;
  // chisq.test / chisq_test_gof: expected counts under the null (row-major for tables, as R's $expected);
  // prop.test: the k x 2 table of expected successes / failures, row-major. Empty for other tests.
  Vec expected;
};

struct TestOptions {
  std::string alternative = "two.sided";  // "two.sided", "less", "greater"
  double mu = 0.0;                         // t / wilcox location; the hypothesised difference for two samples
  bool paired = false;
  bool var_equal = false;
  double conf_level = 0.95;
  int exact = -1;                          // -1 auto (R's NULL), 0 false, 1 true
  bool correct = true;                     // continuity correction (wilcox, chisq 2x2, prop)
  bool conf_int = false;                   // wilcox: compute the Hodges-Lehmann estimate and interval
  double ratio = 1.0;                      // var.test hypothesised ratio
  double p = 0.5;                          // binom.test success probability
  std::string method = "pearson";          // cor.test: "pearson", "kendall", "spearman"
  double or_ = 1.0;                        // fisher.test null odds ratio
};

HtestResult t_test(const Vec& x, const TestOptions& o = {});
HtestResult t_test(const Vec& x, const Vec& y, const TestOptions& o = {});
HtestResult wilcox_test(const Vec& x, const TestOptions& o = {});
HtestResult wilcox_test(const Vec& x, const Vec& y, const TestOptions& o = {});
HtestResult var_test(const Vec& x, const Vec& y, const TestOptions& o = {});
HtestResult cor_test(const Vec& x, const Vec& y, const TestOptions& o = {});
HtestResult shapiro_test(const Vec& x);                         // 3 <= n <= 5000
HtestResult ks_test(const Vec& x, const Vec& y, const TestOptions& o = {});      // two-sample
HtestResult ks_test(const Vec& x, const Dist& d, const TestOptions& o = {});     // one-sample against a distribution
// Goodness of fit: observed counts against probabilities p (rescaled to sum 1 when they do not).
HtestResult chisq_test_gof(const Vec& counts, const Vec& p);
// Contingency table r x c (row-major), Pearson chi-squared; Yates correction for 2 x 2 when o.correct.
HtestResult chisq_test(const Matrix& table, const TestOptions& o = {});
// 2 x 2: conditional MLE odds ratio and exact interval; r x c: exact p-value by enumeration (refuses large tables).
HtestResult fisher_test(const Matrix& table, const TestOptions& o = {});
HtestResult binom_test(double x, double n, const TestOptions& o = {});  // Clopper-Pearson interval
// One or more proportions: successes x[i] out of n[i]. p empty: equality of k >= 2 proportions; else the given p.
HtestResult prop_test(const Vec& x, const Vec& n, const Vec& p = {}, const TestOptions& o = {});

// k-sample, block and symmetry tests (htest2.cpp). Groups / blocks are numeric labels (any distinct values), ordered
// by label value as R's factor of numbers. For the k-sample tests a pair with NaN in x or g is dropped (R drops NA);
// Inf is kept. For the block designs NaN in groups/blocks is an error and a block with a NaN response is removed.
HtestResult kruskal_test(const Vec& x, const Vec& g);   // Kruskal-Wallis with tie correction
HtestResult bartlett_test(const Vec& x, const Vec& g);  // each group needs >= 2 observations
HtestResult fligner_test(const Vec& x, const Vec& g);   // Fligner-Killeen:med
HtestResult oneway_test(const Vec& x, const Vec& g, bool var_equal = false);  // F or Welch (1951)
// Unreplicated complete block designs: exactly one y per (group, block); at least 2 groups and 2 complete blocks.
HtestResult friedman_test(const Vec& y, const Vec& groups, const Vec& blocks);
HtestResult quade_test(const Vec& y, const Vec& groups, const Vec& blocks);
// k x k table: McNemar (k = 2) / Bowker symmetry; continuity correction only for 2 x 2 with n12 != n21.
HtestResult mcnemar_test(const Matrix& table, bool correct = true);
// Ansari-Bradley; uses o.alternative and o.exact (auto: exact when both samples have < 50 finite values and there
// are no ties). NaN is dropped, Inf kept. o.conf_int is not supported (throws).
HtestResult ansari_test(const Vec& x, const Vec& y, const TestOptions& o = {});

}  // namespace rcli
