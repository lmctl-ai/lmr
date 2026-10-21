// Sample autocorrelation and partial autocorrelation functions.
// acf: the usual biased estimator r(k) = c(k) / c(0), c(k) = (1/n) sum_{t=1}^{n-k} (x_t - xbar)(x_{t+k} - xbar)
// (Box & Jenkins 1976, sec. 2.1.4; Brockwell & Davis 1991, sec. 1.4).
// pacf: Durbin-Levinson recursion on the sample autocorrelations (Durbin 1960; Brockwell & Davis 1991, sec. 5.2).
#include <cmath>
#include "rcli/ts.hpp"

namespace rcli {

namespace {
int default_lag_max(size_t n) {
  int lag = static_cast<int>(std::floor(10.0 * std::log10(static_cast<double>(n))));
  if (lag > static_cast<int>(n) - 1) lag = static_cast<int>(n) - 1;
  return lag;
}
void check_series(const Vec& x, const char* fn) {
  if (x.size() < 2) throw Error(std::string(fn) + ": 'x' must have at least 2 observations");
  for (double v : x)
    if (!std::isfinite(v)) throw Error(std::string(fn) + ": 'x' contains missing or non-finite values");
}
}  // namespace

Vec acf(const Vec& x, int lag_max) {
  check_series(x, "acf");
  const size_t n = x.size();
  if (lag_max < 0) lag_max = default_lag_max(n);
  if (lag_max > static_cast<int>(n) - 1) lag_max = static_cast<int>(n) - 1;
  double m = 0.0;
  for (double v : x) m += v;
  m /= static_cast<double>(n);
  Vec y(n);
  for (size_t i = 0; i < n; ++i) y[i] = x[i] - m;
  Vec r(static_cast<size_t>(lag_max) + 1);
  double c0 = 0.0;
  for (size_t t = 0; t < n; ++t) c0 += y[t] * y[t];
  c0 /= static_cast<double>(n);
  for (int k = 0; k <= lag_max; ++k) {
    double s = 0.0;
    for (size_t t = 0; t + static_cast<size_t>(k) < n; ++t) s += y[t] * y[t + static_cast<size_t>(k)];
    r[static_cast<size_t>(k)] = (s / static_cast<double>(n)) / c0;
  }
  return r;
}

Vec pacf(const Vec& x, int lag_max) {
  check_series(x, "pacf");
  const size_t n = x.size();
  if (lag_max < 0) lag_max = default_lag_max(n);
  if (lag_max > static_cast<int>(n) - 1) lag_max = static_cast<int>(n) - 1;
  if (lag_max < 1) throw Error("pacf: 'lag.max' must be at least 1");
  Vec r = acf(x, lag_max);
  const size_t L = static_cast<size_t>(lag_max);
  Vec out(L), phi(L + 1, 0.0), prev(L + 1, 0.0);
  // Durbin-Levinson: phi_kk = (r_k - sum_{j<k} phi_{k-1,j} r_{k-j}) / (1 - sum_{j<k} phi_{k-1,j} r_j)
  for (size_t k = 1; k <= L; ++k) {
    double num = r[k], den = 1.0;
    for (size_t j = 1; j < k; ++j) {
      num -= prev[j] * r[k - j];
      den -= prev[j] * r[j];
    }
    double a = num / den;
    phi[k] = a;
    for (size_t j = 1; j < k; ++j) phi[j] = prev[j] - a * prev[k - j];
    out[k - 1] = a;
    prev = phi;
  }
  return out;
}

}  // namespace rcli
