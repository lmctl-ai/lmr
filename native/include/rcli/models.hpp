// Linear, generalised linear and nonlinear least-squares models (R lm, glm, nls semantics).
// lm uses Householder QR with R's limited column pivoting (a column whose remaining norm falls below tol = 1e-7 times
// its original norm is moved to the end and reported as aliased; coefficient NaN), so rank-deficient designs report the
// same aliased terms as R. glm uses IRLS with R's convergence rule |dev - dev_old| / (|dev| + 0.1) < epsilon.
#pragma once
#include <map>
#include <string>
#include "rcli/common.hpp"

namespace rcli {

// A data table of named numeric columns of equal length.
struct Data {
  std::vector<std::string> names;
  std::vector<Vec> cols;
  size_t nrow() const { return cols.empty() ? 0 : cols[0].size(); }
  const Vec& col(const std::string& name) const;  // throws when absent
  void add(const std::string& name, Vec v);       // throws on length mismatch
};

// Model formula: `y ~ a + b + a:b`, `y ~ a*b` (= a + b + a:b), `y ~ x - 1` / `y ~ 0 + x` (no intercept),
// `y ~ .` (all other columns), and these transforms of one column or of the response: log(v), exp(v), sqrt(v),
// I(v^k) for a numeric literal k, poly is not supported. Terms are numeric (no factors).
struct Formula {
  std::string response;                         // the response expression, e.g. "y" or "log(y)"
  bool intercept = true;
  std::vector<std::string> term_labels;          // R term labels, e.g. "x", "log(z)", "a:b"
  static Formula parse(const std::string& text, const Data* data = nullptr);  // data resolves "."
  Vec response_values(const Data& d) const;
  Matrix design(const Data& d, std::vector<std::string>* colnames = nullptr) const;  // n x p, "(Intercept)" first
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
  Vec predict(const Data& newdata, const std::string& type = "link") const;  // "link" or "response"
};
GlmFit glm(const Formula& f, const Data& d, const GlmOptions& opt = {});

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
