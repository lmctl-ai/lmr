// Linear, generalised linear and nonlinear least-squares models (R lm, glm, nls semantics).
// lm uses Householder QR with R's limited column pivoting (a column whose remaining norm falls below tol = 1e-7 times
// its original norm is moved to the end and reported as aliased; coefficient NaN), so rank-deficient designs report the
// same aliased terms as R. glm uses IRLS with R's convergence rule |dev - dev_old| / (|dev| + 0.1) < epsilon.
#pragma once
#include <map>
#include <string>
#include "rcli/common.hpp"

namespace rcli {

// A data table of named columns of equal length: numeric columns (names/cols) and text (categorical) columns
// (text_names/text_cols). A text column used in a model formula is a factor (R's treatment of character data).
struct Data {
  std::vector<std::string> names;                     // numeric columns
  std::vector<Vec> cols;
  std::vector<std::string> text_names;                // text (categorical) columns
  std::vector<std::vector<std::string>> text_cols;
  size_t nrow() const;                                // from the numeric columns, else the text columns
  const Vec& col(const std::string& name) const;      // numeric column; throws when absent or text
  const std::vector<std::string>& text(const std::string& name) const;  // text column; throws when absent
  bool is_text(const std::string& name) const;
  bool has(const std::string& name) const;            // numeric or text
  void add(const std::string& name, Vec v);           // throws on length mismatch
  void add_text(const std::string& name, std::vector<std::string> v);  // throws on length mismatch
  // All column names in insertion order (columns pushed directly into names/text_names come after the added ones,
  // numeric first). This is the order of R's "." expansion.
  std::vector<std::string> column_order() const;
  std::vector<std::string> order;                     // insertion order maintained by add/add_text
};

// Model formula: `y ~ a + b + a:b`, `y ~ a*b` (= a + b + a:b), `(a + b + c)^2`, `y ~ x - 1` / `y ~ 0 + x` (no
// intercept), `y ~ .` (all other columns, text columns as factors), and these transforms of one column or of the
// response: log(v), exp(v), sqrt(v), I(v^k) for a numeric literal k; poly is not supported.
// Factors: a text column, or factor(v) / as.factor(v) of a numeric (or text) column. Levels are sorted: text in plain
// byte order (R's order in the C locale; R's default locale collation can order mixed case differently), numeric
// values numerically and labelled as R's as.character prints them (15 significant digits, fixed or scientific
// notation, whichever is narrower: 0.5 -> "0.5", 1e5 -> "1e+05"). Coding is by treatment contrasts (first level is the
// baseline) or by indicator columns for every level when the term's margin is absent (Chambers & Hastie 1992 rule,
// see models.cpp). Columns are named as R's model.matrix: "suppVC", "factor(dose)1", "suppVC:factor(dose)1", "wt:am1".
struct Formula {
  std::string response;                         // the response expression, e.g. "y" or "log(y)"
  bool intercept = true;
  std::vector<std::string> term_labels;          // R term labels, e.g. "x", "log(z)", "a:b", "factor(dose)"
  // Every variable the right-hand side mentions, also those removed again by "-" (R's "variables" attribute): lm and
  // glm drop the rows where any of them is missing, as R does (y ~ . - w drops the rows where w is NA).
  std::vector<std::string> variables;
  // Factor levels by variable label (R's xlevels). Set by lm/glm from the fitted rows; when a factor's levels are
  // present here, design() uses them and a value outside them is an error ("factor x has new level q").
  std::map<std::string, std::vector<std::string>> xlevels;
  static Formula parse(const std::string& text, const Data* data = nullptr);  // data resolves "."
  Vec response_values(const Data& d) const;     // throws for a text (factor) response
  // n x p model matrix, "(Intercept)" first. assign: per column 0 for the intercept, k for term_labels[k - 1] (R's
  // "assign" attribute). Rows with a missing value (NaN) in a variable are NaN.
  Matrix design(const Data& d, std::vector<std::string>* colnames = nullptr, std::vector<int>* assign = nullptr) const;
};

struct CoefTable {
  std::vector<std::string> names;
  Vec estimate, std_error, statistic, p_value;  // NaN for aliased
  std::vector<bool> aliased;
  int index(const std::string& name) const;     // -1 when absent
};

struct LmFit {
  Formula formula;
  CoefTable coef;
  int n = 0, rank = 0, df_residual = 0;
  double sigma = kNaN, r_squared = kNaN, adj_r_squared = kNaN;
  double fstatistic = kNaN, f_numdf = kNaN, f_dendf = kNaN, f_pvalue = kNaN;
  double rss = kNaN, loglik = kNaN, aic = kNaN, bic = kNaN;
  Vec residuals, fitted;
  Matrix cov_unscaled;                          // (X'X)^-1 over the non-aliased columns, in coef order (aliased rows 0)
  std::vector<int> assign;                      // per coefficient: 0 intercept, k = formula.term_labels[k - 1]
  Vec effects;                                  // Q'y of the (weighted) least-squares fit, QR pivoted order
  std::vector<int> qr_pivot;                    // qr_pivot[k] = coefficient index at QR position k (first rank used)
  Vec predict(const Data& newdata) const;       // point predictions
  // Interval half-widths for predictions: interval "confidence" or "prediction" at level.
  Vec predict_se(const Data& newdata) const;    // standard error of the fitted mean
};
LmFit lm(const Formula& f, const Data& d, const Vec* weights = nullptr);

struct GlmOptions {
  std::string family = "gaussian";  // "gaussian", "binomial", "poisson", "Gamma"
  std::string link = "default";     // "identity", "log", "inverse", "logit", "probit", "cloglog", "cauchit", "sqrt"
  int maxit = 25;
  double epsilon = 1e-8;
  Vec weights;                      // prior weights; for binomial with proportions these are the trials
  Vec offset;
};
struct GlmFit {
  Formula formula;
  std::string family, link;
  CoefTable coef;                   // z statistics for binomial/poisson, t for gaussian/Gamma (dispersion estimated)
  int n = 0, rank = 0, df_residual = 0, df_null = 0, iter = 0;
  bool converged = false;
  double deviance = kNaN, null_deviance = kNaN, aic = kNaN, dispersion = kNaN;
  Vec fitted, linear_predictors, residuals_deviance, residuals_pearson, residuals_working, weights_working;
  Matrix cov_unscaled;
  std::vector<int> assign;          // per coefficient: 0 intercept, k = formula.term_labels[k - 1]
  Matrix model_matrix;              // the fitted rows (after NA removal), kept for anova refits
  Vec y, prior_weights, offset;     // the fitted rows
  int maxit = 25;
  double epsilon = 1e-8;
  Vec predict(const Data& newdata, const std::string& type = "link") const;  // "link" or "response"
};
GlmFit glm(const Formula& f, const Data& d, const GlmOptions& opt = {});

// An R-style analysis of variance / deviance table: named rows and named numeric columns; NaN is R's NA (blank cell).
struct AnovaTable {
  std::string heading;
  std::vector<std::string> row_names;
  std::vector<std::string> col_names;
  std::vector<Vec> cols;                        // cols[j][i] is row i of column col_names[j]
  const Vec& col(const std::string& name) const;                       // throws when absent
  double at(const std::string& row, const std::string& col) const;     // throws when absent
};
// Sequential (type I) sums of squares (R's anova.lm): one row per term with Df > 0 in term order, then "Residuals";
// columns "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)".
AnovaTable anova(const LmFit& fit);
// Sequential analysis of deviance (R's anova.glm): rows "NULL" and each term; columns "Df", "Deviance", "Resid. Df",
// "Resid. Dev", plus for test "Chisq" (or "LRT") "Pr(>Chi)" and for test "F" "F" and "Pr(>F)". test "none" (or "")
// adds no test. The dispersion of the full model is used for the tests (1 for binomial and poisson, for which the F
// test reduces to the chi-square test).
AnovaTable anova(const GlmFit& fit, const std::string& test = "none");
// Nested model comparison (R's anova(m1, m2)): rows "1", "2"; columns "Res.Df", "RSS", "Df", "Sum of Sq", "F",
// "Pr(>F)". The F test uses the residual mean square of the model with fewer residual df.
AnovaTable anova(const LmFit& smaller, const LmFit& bigger);

struct NlsOptions {
  std::string algorithm = "default";  // "default" (Gauss-Newton with step halving) or "port" (bounded)
  Vec lower, upper;                   // port only; empty = unbounded
  int maxiter = 50;
  double tol = 1e-5;                  // R's relative-offset convergence criterion (Bates & Watts)
  double min_factor = 1.0 / 1024.0;
};
struct NlsFit {
  CoefTable coef;
  int n = 0, df_residual = 0, iterations = 0;
  bool converged = false;
  std::string message;                // "converged", "step factor reduced below minFactor", ...
  double sigma = kNaN, rss = kNaN, conv_criterion = kNaN;
  Vec fitted, residuals;
  Matrix cov_unscaled;
};
// model(par, out) writes the n predictions for parameter vector par into out (out has size n on entry).
// The Jacobian is by forward differences as R's numericDeriv (delta = sqrt(eps) * max(|par|, ...)).
using NlsModel = std::function<void(const Vec& par, Vec& out)>;
NlsFit nls(const NlsModel& model, const Vec& y, const std::vector<std::string>& names, const Vec& start,
           const NlsOptions& opt = {});

}  // namespace rcli
