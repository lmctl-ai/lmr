// k-sample, block-design, paired-table and scale tests with R's htest semantics (clean-room; see htest.hpp).
// Every formula is implemented from the published literature cited below; R was used only as a black-box oracle to
// confirm method strings, branch rules, degenerate-case behaviour and numeric results.
//
// References:
//   Kruskal, W. H. and Wallis, W. A. (1952) Use of ranks in one-criterion variance analysis. JASA 47, 583-621;
//     Hollander, M. and Wolfe, D. A. (1973) Nonparametric Statistical Methods. Wiley, pp. 115-120 (tie correction).
//   Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proc. R. Soc. Lond. A 160, 268-282.
//   Fligner, M. A. and Killeen, T. J. (1976) Distribution-free two-sample tests for scale. JASA 71, 210-213;
//     Conover, W. J., Johnson, M. E. and Johnson, M. M. (1981) A comparative study of tests for homogeneity of
//     variances. Technometrics 23, 351-361 (the "F-K:med X^2" statistic: median-centred samples, scores
//     a(i) = Phi^{-1}((1 + i/(N+1)) / 2), chi-squared approximation of the linear rank statistic).
//   Fisher, R. A. (1925) Statistical Methods for Research Workers (one-way analysis of variance);
//     Welch, B. L. (1951) On the comparison of several mean values: an alternative approach. Biometrika 38, 330-336.
//   Friedman, M. (1937) The use of ranks to avoid the assumption of normality implicit in the analysis of variance.
//     JASA 32, 675-701; Hollander & Wolfe (1973) pp. 139-146 (tie-corrected statistic).
//   McNemar, Q. (1947) Psychometrika 12, 153-157; Bowker, A. H. (1948) A test for symmetry in contingency tables.
//     JASA 43, 572-574 (k x k); Edwards, A. L. (1948) Psychometrika 13, 185-187 (continuity correction);
//     Agresti, A. (1990) Categorical Data Analysis, pp. 350-354.
//   Quade, D. (1979) Using weighted rankings in the analysis of complete blocks with additive block effects.
//     JASA 74, 680-683; Conover, W. J. (1999) Practical Nonparametric Statistics, 3rd ed., pp. 373-380 (including the
//     rule p = (1/k!)^(b-1) when A = B).
//   Ansari, A. R. and Bradley, R. A. (1960) Rank-sum tests for dispersions. Ann. Math. Statist. 31, 1174-1189;
//     Hollander & Wolfe (1973) pp. 83-92 (null mean and variance, with ties); Hajek, J., Sidak, Z. and Sen, P. K.
//     (1999) Theory of Rank Tests, pp. 131ff (mid-rank scores).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>

#include "rcli/dist.hpp"
#include "rcli/htest.hpp"
#include "rcli/summary.hpp"

namespace rcli {
namespace {

// Compensated sum (Neumaier); a plain sum when an element is not finite (the compensation would turn Inf into NaN).
double ksum(const Vec& v) {
  double s = 0.0, c = 0.0;
  for (double x : v)
    if (!std::isfinite(x)) return std::accumulate(v.begin(), v.end(), 0.0);
  for (double x : v) {
    const double t = s + x;
    c += std::fabs(s) >= std::fabs(x) ? (s - t) + x : (x - t) + s;
    s = t;
  }
  return s + c;
}

double chisq_upper(double stat, double df) {
  if (std::isnan(stat) || std::isnan(df)) return kNaN;
  return Dist::chisq(df).p(stat, false);
}
double f_upper(double stat, double df1, double df2) {
  if (std::isnan(stat) || !(df1 > 0) || !(df2 > 0)) return kNaN;
  return Dist::f(df1, df2).p(stat, false);
}

// Observations split by group label. R would drop NA pairs (complete cases); here NaN in x or g drops the pair.
// Groups are ordered by label value (as R's factor of numbers).
struct Groups {
  std::vector<Vec> x;  // samples per group
  Vec all;             // all kept observations, in input order
  std::vector<size_t> gi;  // group index of each kept observation
  size_t n = 0;
};

Groups split(const Vec& x, const Vec& g, const std::string& test) {
  if (x.size() != g.size()) throw Error(test + ": 'x' and 'g' must have the same length");
  std::map<double, size_t> lab;
  for (size_t i = 0; i < x.size(); ++i)
    if (!std::isnan(x[i]) && !std::isnan(g[i])) lab.emplace(g[i] == 0.0 ? 0.0 : g[i], 0);
  if (lab.size() < 2) throw Error(test + ": all observations are in the same group");
  size_t k = 0;
  for (auto& kv : lab) kv.second = k++;
  Groups G;
  G.x.resize(k);
  for (size_t i = 0; i < x.size(); ++i) {
    if (std::isnan(x[i]) || std::isnan(g[i])) continue;
    const size_t j = lab[g[i] == 0.0 ? 0.0 : g[i]];
    G.x[j].push_back(x[i]);
    G.all.push_back(x[i]);
    G.gi.push_back(j);
  }
  G.n = G.all.size();
  return G;
}

// Sum over tie groups of t^3 - t.
double ties_t3(Vec v) {
  std::sort(v.begin(), v.end());
  double s = 0.0;
  for (size_t i = 0; i < v.size();) {
    size_t j = i + 1;
    while (j < v.size() && v[j] == v[i]) ++j;
    const double t = static_cast<double>(j - i);
    s += t * t * t - t;
    i = j;
  }
  return s;
}

double median_of(Vec v) {
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  if (n == 0) return kNaN;
  return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
}

double var_of(const Vec& v) {
  const size_t n = v.size();
  if (n < 2) return kNaN;
  const double m = summary::mean(v);
  Vec d(n);
  for (size_t i = 0; i < n; ++i) d[i] = (v[i] - m) * (v[i] - m);
  return ksum(d) / static_cast<double>(n - 1);
}

HtestResult chisq_result(const char* method, const char* stat_name, double stat, double df) {
  HtestResult r;
  r.method = method;
  r.statistic_name = stat_name;
  r.statistic = stat;
  r.parameters = {{"df", df}};
  r.p_value = chisq_upper(stat, df);
  return r;
}

// Complete block design: y[block][group], groups and blocks ordered by label. NaN in groups/blocks is an error; a
// block with a NaN response is removed (as R documents).
std::vector<Vec> block_design(const Vec& y, const Vec& groups, const Vec& blocks, const std::string& test) {
  if (y.size() != groups.size() || y.size() != blocks.size())
    throw Error(test + ": 'y', 'groups' and 'blocks' must have the same length");
  std::map<double, size_t> gl, bl;
  for (size_t i = 0; i < y.size(); ++i) {
    if (std::isnan(groups[i]) || std::isnan(blocks[i])) throw Error(test + ": NA's are not allowed in 'groups' or 'blocks'");
    gl.emplace(groups[i] == 0.0 ? 0.0 : groups[i], 0);
    bl.emplace(blocks[i] == 0.0 ? 0.0 : blocks[i], 0);
  }
  size_t k = 0, b = 0;
  for (auto& kv : gl) kv.second = k++;
  for (auto& kv : bl) kv.second = b++;
  if (k * b != y.size()) throw Error(test + ": not an unreplicated complete block design");
  std::vector<Vec> m(b, Vec(k, kNaN));
  std::vector<std::vector<char>> seen(b, std::vector<char>(k, 0));
  for (size_t i = 0; i < y.size(); ++i) {
    const size_t gi = gl[groups[i] == 0.0 ? 0.0 : groups[i]], bi = bl[blocks[i] == 0.0 ? 0.0 : blocks[i]];
    if (seen[bi][gi]) throw Error(test + ": not an unreplicated complete block design");
    seen[bi][gi] = 1;
    m[bi][gi] = y[i];
  }
  std::vector<Vec> out;
  for (auto& row : m)
    if (std::none_of(row.begin(), row.end(), [](double v) { return std::isnan(v); })) out.push_back(row);
  if (k < 2) throw Error(test + ": need at least 2 groups");
  if (out.size() < 2) throw Error(test + ": need at least 2 complete blocks");
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// Kruskal-Wallis: H = [12 / (N (N+1)) sum R_j^2 / n_j - 3 (N+1)] / [1 - sum (t^3 - t) / (N^3 - N)], df = k - 1.
HtestResult kruskal_test(const Vec& x, const Vec& g) {
  const Groups G = split(x, g, "kruskal.test");
  const double N = static_cast<double>(G.n);
  const size_t k = G.x.size();
  const Vec r = summary::rank(G.all);
  Vec rs(k, 0.0), ns(k, 0.0);
  for (size_t i = 0; i < G.n; ++i) rs[G.gi[i]] += r[i], ns[G.gi[i]] += 1.0;
  Vec terms(k);
  for (size_t j = 0; j < k; ++j) terms[j] = rs[j] * rs[j] / ns[j];
  double h = 12.0 * ksum(terms) / (N * (N + 1.0)) - 3.0 * (N + 1.0);
  h /= 1.0 - ties_t3(G.all) / (N * N * N - N);
  return chisq_result("Kruskal-Wallis rank sum test", "Kruskal-Wallis chi-squared", h, static_cast<double>(k - 1));
}

// ---------------------------------------------------------------------------------------------------------------
// Bartlett (1937): K^2 = [(N-k) ln s_p^2 - sum (n_j - 1) ln s_j^2] / [1 + (sum 1/(n_j-1) - 1/(N-k)) / (3 (k-1))].
HtestResult bartlett_test(const Vec& x, const Vec& g) {
  const Groups G = split(x, g, "bartlett.test");
  const size_t k = G.x.size();
  for (const Vec& s : G.x)
    if (s.size() < 2) throw Error("bartlett.test: there must be at least 2 observations in each group");
  const double nk = static_cast<double>(G.n - k);
  Vec pooled(k), logs(k), inv(k);
  for (size_t j = 0; j < k; ++j) {
    const double df = static_cast<double>(G.x[j].size() - 1);
    const double v = var_of(G.x[j]);
    pooled[j] = df * v;
    logs[j] = df * std::log(v);
    inv[j] = 1.0 / df;
  }
  const double sp = ksum(pooled) / nk;
  const double kk = static_cast<double>(k);
  const double stat = (nk * std::log(sp) - ksum(logs)) / (1.0 + (ksum(inv) - 1.0 / nk) / (3.0 * (kk - 1.0)));
  return chisq_result("Bartlett test of homogeneity of variances", "Bartlett's K-squared", stat, kk - 1.0);
}

// ---------------------------------------------------------------------------------------------------------------
// Fligner-Killeen:med (Conover, Johnson & Johnson 1981): centre each sample at its median, rank |x - med| over all
// N observations (mid-ranks), score a = qnorm((1 + rank/(N+1)) / 2), X^2 = sum n_j (abar_j - abar)^2 / var(a).
HtestResult fligner_test(const Vec& x, const Vec& g) {
  const Groups G = split(x, g, "fligner.test");
  const size_t k = G.x.size();
  const double N = static_cast<double>(G.n);
  Vec med(k);
  for (size_t j = 0; j < k; ++j) med[j] = median_of(G.x[j]);
  Vec dev(G.n);
  for (size_t i = 0; i < G.n; ++i) dev[i] = std::fabs(G.all[i] - med[G.gi[i]]);
  const Vec r = summary::rank(dev);
  const Dist z = Dist::norm(0.0, 1.0);
  Vec a(G.n);
  for (size_t i = 0; i < G.n; ++i) a[i] = z.q((1.0 + r[i] / (N + 1.0)) / 2.0);
  const double abar = summary::mean(a);
  std::vector<Vec> ag(k);
  for (size_t i = 0; i < G.n; ++i) ag[G.gi[i]].push_back(a[i]);
  Vec terms(k);
  for (size_t j = 0; j < k; ++j) {
    const double d = summary::mean(ag[j]) - abar;
    terms[j] = static_cast<double>(ag[j].size()) * d * d;
  }
  const double stat = ksum(terms) / var_of(a);
  return chisq_result("Fligner-Killeen test of homogeneity of variances", "Fligner-Killeen:med chi-squared", stat,
                      static_cast<double>(k - 1));
}

// ---------------------------------------------------------------------------------------------------------------
// One-way analysis of means. var_equal: the classical F = [SSB/(k-1)] / [SSW/(N-k)]. Otherwise Welch (1951):
// w_j = n_j / s_j^2, m = sum w_j xbar_j / sum w_j, A = sum w_j (xbar_j - m)^2 / (k-1),
// L = sum (1 - w_j / sum w)^2 / (n_j - 1) / (k^2 - 1), F = A / (1 + 2 (k-2) L) on (k-1, 1 / (3 L)) df.
HtestResult oneway_test(const Vec& x, const Vec& g, bool var_equal) {
  if (x.size() != g.size()) throw Error("oneway.test: 'x' and 'g' must have the same length");
  size_t present = 0;
  {
    std::map<double, int> lab;
    for (size_t i = 0; i < x.size(); ++i)
      if (!std::isnan(x[i]) && !std::isnan(g[i])) lab[g[i] == 0.0 ? 0.0 : g[i]] = 1;
    present = lab.size();
  }
  if (present < 2) throw Error("oneway.test: not enough groups");
  const Groups G = split(x, g, "oneway.test");
  const size_t k = G.x.size();
  for (const Vec& s : G.x)
    if (s.size() < 2) throw Error("oneway.test: not enough observations");
  const double kk = static_cast<double>(k), N = static_cast<double>(G.n);
  Vec m(k), v(k), n(k);
  for (size_t j = 0; j < k; ++j) {
    n[j] = static_cast<double>(G.x[j].size());
    m[j] = summary::mean(G.x[j]);
    v[j] = var_of(G.x[j]);
  }
  HtestResult r;
  r.statistic_name = "F";
  if (var_equal) {
    r.method = "One-way analysis of means";
    const double gm = summary::mean(G.all);
    Vec ssb(k), ssw(k);
    for (size_t j = 0; j < k; ++j) {
      ssb[j] = n[j] * (m[j] - gm) * (m[j] - gm);
      ssw[j] = (n[j] - 1.0) * v[j];
    }
    const double df1 = kk - 1.0, df2 = N - kk;
    r.statistic = (ksum(ssb) / df1) / (ksum(ssw) / df2);
    r.parameters = {{"num df", df1}, {"denom df", df2}};
    r.p_value = f_upper(r.statistic, df1, df2);
    return r;
  }
  r.method = "One-way analysis of means (not assuming equal variances)";
  Vec w(k), wm(k);
  for (size_t j = 0; j < k; ++j) {
    w[j] = n[j] / v[j];
    wm[j] = w[j] * m[j];
  }
  const double sw = ksum(w);
  const double mw = ksum(wm) / sw;
  Vec ta(k), tl(k);
  for (size_t j = 0; j < k; ++j) {
    ta[j] = w[j] * (m[j] - mw) * (m[j] - mw);
    const double q = 1.0 - w[j] / sw;
    tl[j] = q * q / (n[j] - 1.0);
  }
  const double lam = ksum(tl) / (kk * kk - 1.0);
  const double df1 = kk - 1.0, df2 = 1.0 / (3.0 * lam);
  r.statistic = ksum(ta) / df1 / (1.0 + 2.0 * (kk - 2.0) * lam);
  r.parameters = {{"num df", df1}, {"denom df", df2}};
  r.p_value = f_upper(r.statistic, df1, df2);
  return r;
}

// ---------------------------------------------------------------------------------------------------------------
// Friedman: rank within blocks (mid-ranks), R_j the group rank sums over b blocks,
// X^2 = 12 sum (R_j - b (k+1)/2)^2 / [b k (k+1) - sum_blocks sum (t^3 - t) / (k - 1)], df = k - 1.
HtestResult friedman_test(const Vec& y, const Vec& groups, const Vec& blocks) {
  const std::vector<Vec> m = block_design(y, groups, blocks, "friedman.test");
  const size_t b = m.size(), k = m[0].size();
  const double bb = static_cast<double>(b), kk = static_cast<double>(k);
  std::vector<Vec> cols(k, Vec(b));
  Vec tie(b);
  for (size_t i = 0; i < b; ++i) {
    const Vec r = summary::rank(m[i]);
    for (size_t j = 0; j < k; ++j) cols[j][i] = r[j];
    tie[i] = ties_t3(m[i]);
  }
  Vec terms(k);
  for (size_t j = 0; j < k; ++j) {
    const double d = ksum(cols[j]) - bb * (kk + 1.0) / 2.0;
    terms[j] = d * d;
  }
  const double stat = 12.0 * ksum(terms) / (bb * kk * (kk + 1.0) - ksum(tie) / (kk - 1.0));
  return chisq_result("Friedman rank sum test", "Friedman chi-squared", stat, kk - 1.0);
}

// ---------------------------------------------------------------------------------------------------------------
// Quade (1979), as in Conover (1999): R_ij the within-block mid-ranks, Q_i the mid-rank of block i's range,
// S_ij = Q_i (R_ij - (k+1)/2), A = sum S_ij^2, B = sum_j (sum_i S_ij)^2 / b, F = (b-1) B / (A - B) on
// (k-1, (b-1)(k-1)) df. When A = B the F statistic is undefined and Conover's exact rule p = (1/k!)^(b-1) is used
// (statistic and parameters are then NaN, as R reports them).
HtestResult quade_test(const Vec& y, const Vec& groups, const Vec& blocks) {
  const std::vector<Vec> m = block_design(y, groups, blocks, "quade.test");
  const size_t b = m.size(), k = m[0].size();
  const double bb = static_cast<double>(b), kk = static_cast<double>(k);
  Vec range(b);
  for (size_t i = 0; i < b; ++i) {
    const auto mm = std::minmax_element(m[i].begin(), m[i].end());
    range[i] = *mm.second - *mm.first;
  }
  const Vec q = summary::rank(range);
  Vec sq;
  std::vector<Vec> cols(k, Vec(b));
  for (size_t i = 0; i < b; ++i) {
    const Vec r = summary::rank(m[i]);
    for (size_t j = 0; j < k; ++j) {
      const double s = q[i] * (r[j] - (kk + 1.0) / 2.0);
      cols[j][i] = s;
      sq.push_back(s * s);
    }
  }
  const double A = ksum(sq);
  Vec bt(k);
  for (size_t j = 0; j < k; ++j) {
    const double s = ksum(cols[j]);
    bt[j] = s * s;
  }
  const double B = ksum(bt) / bb;
  HtestResult r;
  r.method = "Quade test";
  r.statistic_name = "Quade F";
  if (A == B) {
    double kf = 1.0;
    for (size_t i = 2; i <= k; ++i) kf *= static_cast<double>(i);
    r.statistic = kNaN;
    r.parameters = {{"num df", kNaN}, {"denom df", kNaN}};
    r.p_value = std::pow(1.0 / kf, bb - 1.0);
    return r;
  }
  const double df1 = kk - 1.0, df2 = (bb - 1.0) * (kk - 1.0);
  r.statistic = (bb - 1.0) * B / (A - B);
  r.parameters = {{"num df", df1}, {"denom df", df2}};
  r.p_value = f_upper(r.statistic, df1, df2);
  return r;
}

// ---------------------------------------------------------------------------------------------------------------
// McNemar / Bowker symmetry test: X^2 = sum_{i<j} (n_ij - n_ji)^2 / (n_ij + n_ji), df = k (k-1) / 2. For 2 x 2 with
// correct and n_12 != n_21, Edwards' correction (|n_12 - n_21| - 1)^2 / (n_12 + n_21).
HtestResult mcnemar_test(const Matrix& t, bool correct) {
  if (t.rows < 2 || t.rows != t.cols || t.a.size() != t.rows * t.cols)
    throw Error("mcnemar.test: 'x' must be square with at least two rows and columns");
  for (double v : t.a)
    if (!std::isfinite(v) || v < 0) throw Error("mcnemar.test: all entries of 'x' must be nonnegative and finite");
  const size_t k = t.rows;
  HtestResult r;
  r.statistic_name = "McNemar's chi-squared";
  r.parameters = {{"df", static_cast<double>(k * (k - 1) / 2)}};
  if (correct && k == 2 && t(0, 1) != t(1, 0)) {
    r.method = "McNemar's Chi-squared test with continuity correction";
    const double d = std::fabs(t(0, 1) - t(1, 0)) - 1.0;
    r.statistic = d * d / (t(0, 1) + t(1, 0));
  } else {
    r.method = "McNemar's Chi-squared test";
    Vec terms;
    for (size_t i = 0; i < k; ++i)
      for (size_t j = i + 1; j < k; ++j) {
        const double d = t(i, j) - t(j, i);
        terms.push_back(d * d / (t(i, j) + t(j, i)));
      }
    r.statistic = ksum(terms);
  }
  r.p_value = chisq_upper(r.statistic, r.parameters[0].second);
  return r;
}

// ---------------------------------------------------------------------------------------------------------------
// Ansari-Bradley: AB = sum over x of min(r, N - r + 1), r the mid-ranks in the pooled sample. Exact null distribution
// (no ties) by counting m-subsets of the scores with a given sum; otherwise the normal approximation with
// E = m (N+2)/4 (N even) or m (N+1)^2 / (4N) (N odd) and Var = m n / (N (N-1)) (sum a_i^2 - N (E/m)^2), which is
// Hollander & Wolfe's tie-corrected variance written over the mid-rank scores (no continuity correction).
// "greater" (x more dispersed) rejects for small AB.
HtestResult ansari_test(const Vec& x0, const Vec& y0, const TestOptions& o) {
  static const char* names[] = {"two.sided", "less", "greater"};
  int alt = -1;
  if (!o.alternative.empty())
    for (int i = 0; i < 3; ++i)
      if (std::string(names[i]).compare(0, o.alternative.size(), o.alternative) == 0) alt = i;
  if (alt < 0) throw Error("ansari.test: 'alternative' must be one of \"two.sided\", \"less\", \"greater\"");
  if (o.conf_int) throw Error("ansari.test: conf.int is not supported");
  // NaN is dropped (R drops NA); +-Inf is kept and ranked, as R does. The automatic exact rule counts finite values.
  Vec x, y;
  size_t fx = 0, fy = 0;
  for (double v : x0)
    if (!std::isnan(v)) x.push_back(v), fx += std::isfinite(v) ? 1 : 0;
  for (double v : y0)
    if (!std::isnan(v)) y.push_back(v), fy += std::isfinite(v) ? 1 : 0;
  if (x.empty()) throw Error("ansari.test: not enough 'x' observations");
  if (y.empty()) throw Error("ansari.test: not enough 'y' observations");
  const size_t m = x.size(), n = y.size(), N = m + n;
  const double dm = static_cast<double>(m), dn = static_cast<double>(n), dN = static_cast<double>(N);
  Vec all = x;
  all.insert(all.end(), y.begin(), y.end());
  const Vec r = summary::rank(all);
  Vec a(N);
  for (size_t i = 0; i < N; ++i) a[i] = std::min(r[i], dN - r[i] + 1.0);
  const double stat = ksum(Vec(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(m)));
  Vec sorted = all;
  std::sort(sorted.begin(), sorted.end());
  const bool ties = std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end();
  const bool exact = (o.exact == -1 ? (fx < 50 && fy < 50) : o.exact == 1) && !ties;

  HtestResult res;
  res.method = "Ansari-Bradley test";
  res.alternative = names[alt];
  res.statistic_name = "AB";
  res.statistic = stat;
  res.null_name = "ratio of scales";
  res.null_value = 1.0;
  if (exact) {
    if (static_cast<double>(m) * dN * dN * dN / 4.0 > 4e9) throw Error("ansari.test: samples too large for the exact distribution");
    // Integer scores 1, 1, 2, 2, ... (N even) or 1, 1, ..., (N+1)/2 (N odd); cnt[j][s] = number of j-subsets of the
    // scores seen so far with sum s.
    const size_t smax = static_cast<size_t>(ksum(a)) + 1;
    std::vector<Vec> cnt(m + 1, Vec(smax + 1, 0.0));
    cnt[0][0] = 1.0;
    for (size_t i = 1; i <= N; ++i) {
      const size_t sc = std::min(i, N - i + 1);
      for (size_t j = std::min(i, m); j >= 1; --j)
        for (size_t s = smax; s >= sc; --s) cnt[j][s] += cnt[j - 1][s - sc];
    }
    const Vec& f = cnt[m];
    const double total = ksum(f);
    const size_t s0 = static_cast<size_t>(stat);
    Vec le(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(s0 + 1)), ge(f.begin() + static_cast<std::ptrdiff_t>(s0), f.end());
    const double ple = std::min(1.0, ksum(le) / total), pge = std::min(1.0, ksum(ge) / total);
    res.p_value = alt == 0 ? std::min(1.0, 2.0 * std::min(ple, pge)) : alt == 1 ? pge : ple;
    return res;
  }
  const double mu = N % 2 == 0 ? (dN + 2.0) / 4.0 : (dN + 1.0) * (dN + 1.0) / (4.0 * dN);
  Vec a2(N);
  for (size_t i = 0; i < N; ++i) a2[i] = a[i] * a[i];
  const double v = dm * dn / (dN * (dN - 1.0)) * (ksum(a2) - dN * mu * mu);
  const double z = (stat - dm * mu) / std::sqrt(v);
  const Dist nd = Dist::norm(0.0, 1.0);
  const double pl = nd.p(z, true), pu = nd.p(z, false);
  res.p_value = alt == 0 ? std::min(1.0, 2.0 * std::min(pl, pu)) : alt == 1 ? pu : pl;
  return res;
}

}  // namespace rcli
