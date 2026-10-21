// Time series: autocorrelation, partial autocorrelation, ARIMA estimation/forecasting, simulation and STL.
// ARIMA follows the published methods R documents for arima(): exact Gaussian likelihood through a state-space
// Kalman filter (Gardner, Harvey & Phillips 1980; Durbin & Koopman), conditional sum of squares, the CSS-ML default,
// and the Jones (1980) reparameterisation that keeps AR parts stationary during ML optimisation. Optimiser: BFGS from
// optim.hpp. Measured against R 4.3.3: coefficients to about 1e-11 and standard errors to 6e-8 for ordinary models
// (see tests/test_ts.cpp); near-non-stationary seasonal fits can differ (R starts from Gardner et al. 1980).
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

class Rng;

Vec acf(const Vec& x, int lag_max = -1);   // lags 0..lag_max; default lag_max = floor(10 * log10(n)), capped at n - 1
Vec pacf(const Vec& x, int lag_max = -1);  // lags 1..lag_max (Durbin-Levinson)

struct ArimaSpec {
  int p = 0, d = 0, q = 0;
  int P = 0, D = 0, Q = 0, period = 1;   // seasonal part; period > 1 when P, D or Q > 0
  bool include_mean = true;              // only used when d + D == 0
  std::string method = "CSS-ML";         // "CSS-ML", "ML", "CSS"
};
struct ArimaForecast { Vec mean, se; };
struct ArimaFit {
  ArimaSpec spec;
  std::vector<std::string> coef_names;   // "ar1".., "ma1".., "sar1".., "sma1".., "intercept"
  Vec coef, se;
  Matrix vcov;
  double sigma2 = kNaN, loglik = kNaN, aic = kNaN;  // aic NaN for CSS, as R
  int nobs = 0;                          // used observations (n - d - D*period)
  int convergence = 0;                   // optim code
  Vec residuals;                         // length n (innovations)
  Vec x;                                 // the data (kept for forecasting)
  ArimaForecast forecast(int n_ahead) const;
};
ArimaFit arima(const Vec& x, const ArimaSpec& spec);

// Simulate an ARMA(p, q) series (optionally integrated d times, as arima.sim with order c(p, d, q)) with Gaussian
// innovations of standard deviation sd. Burn-in n_start defaults like R: p + q + ceil(6 / log(minroots)) for AR parts.
// Requires a stationary AR part. Returns n values (n + d when d > 0, starting from 0 like arima.sim).
Vec arima_sim(Rng& rng, int n, const Vec& ar, const Vec& ma, double sd = 1.0, int d = 0, int n_start = -1);

struct StlOptions {
  int s_window = -1;          // odd >= 7 (R requires >= 3); -1 with periodic = true means "periodic"
  bool periodic = false;
  int s_degree = 0;
  int t_window = -1;          // default: next odd ceil(1.5 * period / (1 - 1.5 / s_window))
  int t_degree = 1;
  int l_window = -1;          // default: next odd >= period
  bool robust = false;        // inner 1 / outer 15 when robust, inner 2 / outer 0 otherwise (R defaults)
};
struct StlFit { Vec seasonal, trend, remainder, weights; };
// Cleveland, Cleveland, McRae & Terpenning (1990). period >= 2 and at least two full periods.
StlFit stl(const Vec& x, int period, const StlOptions& opt = {});

}  // namespace rcli
