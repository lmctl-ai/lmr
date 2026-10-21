// ARIMA estimation, forecasting and simulation.
// Sources (published methods, as documented in R's ?arima, ?KalmanLike, ?predict.Arima and ?arima.sim):
// - State-space form of an ARMA(p, q) process with state dimension r = max(p, q + 1) (Harvey 1993, sec. 3.3;
//   Brockwell & Davis 1996, sec. 8.3): a_{t+1} = T a_t + R e_{t+1}, w_t = a_t[0], T has the AR coefficients in its
//   first column and an identity superdiagonal, R = (1, theta_1, ..., theta_{r-1})'.
// - Differencing handled inside the state (Durbin & Koopman 2001, sec. 3.2.4 / 5.6): the state is augmented with the
//   last d + D*s observations, which receive a diffuse prior kappa * I (kappa = 1e6); observations whose prediction
//   variance is still at least 1e4 are excluded from the likelihood (?arima, "Fitting methods").
// - Initial state covariance of the stationary ARMA part: the solution of the discrete Lyapunov equation
//   P = T P T' + R R' (Gardner, Harvey & Phillips 1980), solved here directly as a linear system in vech(P).
// - Exact Gaussian likelihood by the Kalman filter prediction-error decomposition, concentrated over sigma^2:
//   0.5 log(S / nu) + 0.5 sum(log F_t) / nu, S = sum v_t^2 / F_t (Harvey 1993, sec. 3.4).
// - Conditional sum of squares (Box & Jenkins 1976, sec. 7.1): innovations from observation n.cond on, earlier
//   innovations set to zero.
// - Jones (1980) reparameterisation: AR parameters through tanh of unconstrained values, read as partial
//   autocorrelations and mapped to AR coefficients by the Durbin-Levinson recursion; applied separately to the
//   non-seasonal and seasonal AR parts during ML optimisation. MA parts are made invertible afterwards by reflecting
//   roots inside the unit circle (?arima).
// - Standard errors from the inverse of the numerically differentiated Hessian of the objective in the original
//   parameterisation, scaled by the number of used observations.
#include <algorithm>
#include <cmath>
#include <complex>
#include "rcli/optim.hpp"
#include "rcli/rng.hpp"
#include "rcli/ts.hpp"

namespace rcli {

namespace {

constexpr double kKappa = 1e6;
constexpr double kDiffuseLimit = 1e4;
const double kLog2Pi = std::log(2.0 * 3.14159265358979323846);

using cplx = std::complex<double>;

// Roots of c[0] + c[1] z + ... + c[m] z^m (c[m] != 0) by the Aberth-Ehrlich simultaneous iteration
// (Aberth 1973), polished with Newton steps.
std::vector<cplx> poly_roots(Vec c) {
  while (!c.empty() && c.back() == 0.0) c.pop_back();
  const size_t m = c.empty() ? 0 : c.size() - 1;
  std::vector<cplx> z(m);
  if (m == 0) return z;
  auto eval = [&](cplx x, cplx& d) {
    cplx p = c[m];
    d = 0.0;
    for (size_t k = m; k-- > 0;) {
      d = d * x + p;
      p = p * x + c[k];
    }
    return p;
  };
  // initial points on a circle of radius given by the Cauchy-type bound
  double rad = 0.0;
  for (size_t k = 0; k < m; ++k) rad = std::max(rad, std::pow(std::fabs(c[k] / c[m]), 1.0 / static_cast<double>(m - k)));
  if (rad == 0.0) rad = 1.0;
  for (size_t k = 0; k < m; ++k) z[k] = std::polar(rad, 2.0 * 3.14159265358979323846 * (k + 0.25) / m + 0.4);
  for (int it = 0; it < 500; ++it) {
    double maxstep = 0.0;
    for (size_t k = 0; k < m; ++k) {
      cplx d;
      cplx p = eval(z[k], d);
      if (p == 0.0) continue;
      cplx ratio = p / d;
      cplx s = 0.0;
      for (size_t j = 0; j < m; ++j)
        if (j != k) s += 1.0 / (z[k] - z[j]);
      cplx step = ratio / (1.0 - ratio * s);
      z[k] -= step;
      maxstep = std::max(maxstep, std::abs(step) / std::max(1.0, std::abs(z[k])));
    }
    if (maxstep < 1e-15) break;
  }
  for (size_t k = 0; k < m; ++k)
    for (int it = 0; it < 3; ++it) {
      cplx d;
      cplx p = eval(z[k], d);
      if (d == 0.0) break;
      z[k] -= p / d;
    }
  return z;
}

// Jones (1980): unconstrained values -> partial autocorrelations tanh(u) -> AR coefficients (Durbin-Levinson).
Vec pacf_to_ar(const Vec& u) {
  const size_t p = u.size();
  Vec phi(p), work(p);
  for (size_t k = 0; k < p; ++k) phi[k] = std::tanh(u[k]);
  for (size_t k = 1; k < p; ++k) {
    double a = phi[k];
    for (size_t j = 0; j < k; ++j) work[j] = phi[j] - a * phi[k - j - 1];
    for (size_t j = 0; j < k; ++j) phi[j] = work[j];
  }
  return phi;
}

// Inverse of pacf_to_ar without the tanh: AR coefficients -> partial autocorrelations (reverse Durbin-Levinson).
// Returns false if some |partial autocorrelation| >= 1, i.e. the AR polynomial is not stationary (Schur-Cohn).
bool ar_to_pacf(const Vec& phi_in, Vec& r) {
  const size_t p = phi_in.size();
  Vec phi = phi_in, work(p);
  r.assign(p, 0.0);
  for (size_t k = p; k-- > 0;) {
    double a = phi[k];
    r[k] = a;
    if (!(std::fabs(a) < 1.0)) return false;
    for (size_t j = 0; j < k; ++j) work[j] = (phi[j] + a * phi[k - j - 1]) / (1.0 - a * a);
    for (size_t j = 0; j < k; ++j) phi[j] = work[j];
  }
  return true;
}

bool ar_stationary(const Vec& phi) {
  if (phi.empty()) return true;
  Vec r;
  return ar_to_pacf(phi, r);
}

// Reflect roots of 1 + sum theta_j z^j lying inside the unit circle (invertible MA representation).
Vec ma_invert(const Vec& theta) {
  size_t q0 = theta.size();
  while (q0 > 0 && theta[q0 - 1] == 0.0) --q0;
  if (q0 == 0) return theta;
  Vec c(q0 + 1);
  c[0] = 1.0;
  for (size_t j = 0; j < q0; ++j) c[j + 1] = theta[j];
  std::vector<cplx> roots = poly_roots(c);
  bool any = false;
  for (auto& z : roots)
    if (std::abs(z) < 1.0) { z = 1.0 / z; any = true; }
  if (!any) return theta;
  // rebuild prod (1 - z / root_k)
  std::vector<cplx> poly(1, 1.0);
  for (const auto& z : roots) {
    std::vector<cplx> next(poly.size() + 1, 0.0);
    for (size_t k = 0; k < poly.size(); ++k) {
      next[k] += poly[k];
      next[k + 1] -= poly[k] / z;
    }
    poly = next;
  }
  Vec out = theta;
  for (size_t j = 0; j < q0; ++j) out[j] = poly[j + 1].real();
  return out;
}

struct Layout {
  int p, q, P, Q, s, ncoef, nmean;
  bool mean;
};

Layout layout_of(const ArimaSpec& sp) {
  Layout L{};
  L.p = sp.p; L.q = sp.q; L.P = sp.P; L.Q = sp.Q; L.s = sp.period;
  L.mean = sp.include_mean && sp.d + sp.D == 0;
  L.nmean = L.mean ? 1 : 0;
  L.ncoef = sp.p + sp.q + sp.P + sp.Q + L.nmean;
  return L;
}

// Expanded (multiplicative seasonal) AR and MA coefficient vectors.
void expand(const Layout& L, const Vec& c, Vec& phi, Vec& theta) {
  const int s = L.s;
  phi.assign(static_cast<size_t>(L.p + s * L.P), 0.0);
  theta.assign(static_cast<size_t>(L.q + s * L.Q), 0.0);
  const double* ar = c.data();
  const double* ma = ar + L.p;
  const double* sar = ma + L.q;
  const double* sma = sar + L.P;
  for (int i = 0; i < L.p; ++i) phi[static_cast<size_t>(i)] = ar[i];
  for (int i = 0; i < L.q; ++i) theta[static_cast<size_t>(i)] = ma[i];
  for (int j = 0; j < L.P; ++j) {
    phi[static_cast<size_t>((j + 1) * s - 1)] += sar[j];
    for (int i = 0; i < L.p; ++i) phi[static_cast<size_t>((j + 1) * s + i)] -= ar[i] * sar[j];
  }
  for (int j = 0; j < L.Q; ++j) {
    theta[static_cast<size_t>((j + 1) * s - 1)] += sma[j];
    for (int i = 0; i < L.q; ++i) theta[static_cast<size_t>((j + 1) * s + i)] += ma[i] * sma[j];
  }
}

// Coefficients of (1 - B)^d (1 - B^s)^D written as 1 - sum delta_i B^i.
Vec delta_of(int d, int D, int s) {
  Vec poly{1.0};
  auto mul = [&](int lag) {
    Vec next(poly.size() + static_cast<size_t>(lag), 0.0);
    for (size_t k = 0; k < poly.size(); ++k) {
      next[k] += poly[k];
      next[k + static_cast<size_t>(lag)] -= poly[k];
    }
    poly = next;
  };
  for (int i = 0; i < d; ++i) mul(1);
  for (int i = 0; i < D; ++i) mul(s);
  Vec delta(poly.size() - 1);
  for (size_t i = 1; i < poly.size(); ++i) delta[i - 1] = -poly[i];
  return delta;
}

// Solve the dense system A x = b by Gaussian elimination with partial pivoting. Returns false if singular.
bool solve_dense(std::vector<double>& A, Vec& b, size_t m) {
  for (size_t k = 0; k < m; ++k) {
    size_t piv = k;
    double best = std::fabs(A[k * m + k]);
    for (size_t i = k + 1; i < m; ++i)
      if (std::fabs(A[i * m + k]) > best) { best = std::fabs(A[i * m + k]); piv = i; }
    if (best == 0.0) return false;
    if (piv != k) {
      for (size_t j = 0; j < m; ++j) std::swap(A[k * m + j], A[piv * m + j]);
      std::swap(b[k], b[piv]);
    }
    const double inv = 1.0 / A[k * m + k];
    for (size_t i = k + 1; i < m; ++i) {
      double f = A[i * m + k] * inv;
      if (f == 0.0) continue;
      for (size_t j = k; j < m; ++j) A[i * m + j] -= f * A[k * m + j];
      b[i] -= f * b[k];
    }
  }
  for (size_t k = m; k-- > 0;) {
    double s = b[k];
    for (size_t j = k + 1; j < m; ++j) s -= A[k * m + j] * b[j];
    b[k] = s / A[k * m + k];
  }
  return true;
}

// Stationary covariance of the ARMA state, i.e. the solution of P = T P T' + R R' (Gardner, Harvey & Phillips 1980),
// computed from the process autocovariances (Brockwell & Davis 1991, sec. 3.3, method 2). With theta_0 = 1 the
// state elements are a[i] = sum_{j=i+1}^{r} phi_j w_{t+i-j} + sum_{j=i}^{r-1} theta_j e_{t+i-j}, so P follows from
// gamma(h) = cov(w_t, w_{t-h}), cov(w_a, e_b) = psi_{a-b} (zero for a < b) and white-noise e (unit variance).
Matrix stationary_cov(const Vec& phi, const Vec& theta, size_t r) {
  const size_t p = phi.size();
  auto th = [&](size_t j) { return j == 0 ? 1.0 : (j - 1 < theta.size() ? theta[j - 1] : 0.0); };
  auto ph = [&](size_t j) { return (j >= 1 && j <= p) ? phi[j - 1] : 0.0; };
  // psi weights up to lag 2r
  const size_t H = 2 * r + 2;
  Vec psi(H, 0.0);
  for (size_t j = 0; j < H; ++j) {
    double v = th(j);
    for (size_t i = 1; i <= std::min(j, p); ++i) v += phi[i - 1] * psi[j - i];
    psi[j] = v;
  }
  // gamma(0..H-1): solve gamma(k) - sum_j phi_j gamma(|k - j|) = sum_{j>=k} theta_j psi_{j-k} for k = 0..p
  const size_t qmax = theta.size();
  auto rhs = [&](size_t k) {
    double c = 0.0;
    for (size_t j = k; j <= qmax; ++j) c += th(j) * psi[j - k];
    return c;
  };
  Vec gamma(H, 0.0);
  {
    const size_t m = p + 1;
    std::vector<double> A(m * m, 0.0);
    Vec b(m);
    for (size_t k = 0; k <= p; ++k) {
      A[k * m + k] += 1.0;
      for (size_t j = 1; j <= p; ++j) {
        size_t lag = k >= j ? k - j : j - k;
        A[k * m + lag] -= phi[j - 1];
      }
      b[k] = rhs(k);
    }
    if (!solve_dense(A, b, m)) throw Error("arima: non-stationary AR part");
    for (size_t k = 0; k <= p; ++k) gamma[k] = b[k];
    for (size_t k = p + 1; k < H; ++k) {
      double v = rhs(k);
      for (size_t j = 1; j <= p; ++j) v += phi[j - 1] * gamma[k - j];
      gamma[k] = v;
    }
  }
  auto gam = [&](long h) { return gamma[static_cast<size_t>(h < 0 ? -h : h)]; };
  auto cwe = [&](long a, long b) { return a >= b ? psi[static_cast<size_t>(a - b)] : 0.0; };  // cov(w_a, e_b)
  Matrix P(r, r, 0.0);
  const long R = static_cast<long>(r);
  for (long i = 0; i < R; ++i)
    for (long k = i; k < R; ++k) {
      double v = 0.0;
      // w-w terms: phi_j phi_l gamma((i - j) - (k - l))
      for (long j = i + 1; j <= R; ++j) {
        double pj = ph(static_cast<size_t>(j));
        if (pj == 0.0) continue;
        for (long l = k + 1; l <= R; ++l) {
          double pl = ph(static_cast<size_t>(l));
          if (pl != 0.0) v += pj * pl * gam((i - j) - (k - l));
        }
        // w-e terms: phi_j theta_l cov(w_{i-j}, e_{k-l})
        for (long l = k; l < R; ++l) {
          double tl = th(static_cast<size_t>(l));
          if (tl != 0.0) v += pj * tl * cwe(i - j, k - l);
        }
      }
      for (long j = i; j < R; ++j) {
        double tj = th(static_cast<size_t>(j));
        if (tj == 0.0) continue;
        // e-w terms: theta_j phi_l cov(w_{k-l}, e_{i-j})
        for (long l = k + 1; l <= R; ++l) {
          double pl = ph(static_cast<size_t>(l));
          if (pl != 0.0) v += tj * pl * cwe(k - l, i - j);
        }
        // e-e terms
        long l = k - i + j;
        if (l >= k && l < R) v += tj * th(static_cast<size_t>(l));
      }
      P(static_cast<size_t>(i), static_cast<size_t>(k)) = P(static_cast<size_t>(k), static_cast<size_t>(i)) = v;
    }
  return P;
}

// Complete state-space model of the (differenced) ARMA process.
struct SSModel {
  size_t r = 0, nd = 0, rd = 0;
  Vec phi, theta, delta;
  Vec Z;         // observation vector (length rd)
  Matrix V;      // R R' padded to rd x rd
  Matrix P0;     // initial state covariance
  // sparse transition matrix: entries (row, col, value)
  std::vector<size_t> ti, tj;
  Vec tv;
};

SSModel make_model(const Vec& phi, const Vec& theta, const Vec& delta) {
  SSModel m;
  m.phi = phi; m.theta = theta; m.delta = delta;
  m.r = std::max(phi.size(), theta.size() + 1);
  m.nd = delta.size();
  m.rd = m.r + m.nd;
  const size_t r = m.r, rd = m.rd;
  m.Z.assign(rd, 0.0);
  m.Z[0] = 1.0;
  for (size_t i = 0; i < m.nd; ++i) m.Z[r + i] = delta[i];
  auto add = [&](size_t i, size_t j, double v) {
    if (v == 0.0) return;
    m.ti.push_back(i); m.tj.push_back(j); m.tv.push_back(v);
  };
  for (size_t i = 0; i < r; ++i) {
    if (i < phi.size()) add(i, 0, phi[i]);
    if (i + 1 < r) add(i, i + 1, 1.0);
  }
  if (m.nd > 0) {
    for (size_t j = 0; j < rd; ++j) add(r, j, m.Z[j]);  // new lagged observation y_t = Z a_t
    for (size_t i = 1; i < m.nd; ++i) add(r + i, r + i - 1, 1.0);
  }
  Vec R(r, 0.0);
  R[0] = 1.0;
  for (size_t i = 1; i < r; ++i) R[i] = i - 1 < theta.size() ? theta[i - 1] : 0.0;
  m.V = Matrix(rd, rd, 0.0);
  for (size_t i = 0; i < r; ++i)
    for (size_t j = 0; j < r; ++j) m.V(i, j) = R[i] * R[j];
  m.P0 = Matrix(rd, rd, 0.0);
  Matrix Q0 = stationary_cov(phi, theta, r);
  for (size_t i = 0; i < r; ++i)
    for (size_t j = 0; j < r; ++j) m.P0(i, j) = Q0(i, j);
  for (size_t i = r; i < rd; ++i) m.P0(i, i) = kKappa;
  return m;
}

struct FilterState {
  Vec a;
  Matrix P;
};

// a <- T a; P <- T P T' + V
void predict_step(const SSModel& m, FilterState& st) {
  const size_t rd = m.rd, nnz = m.tv.size();
  Vec an(rd, 0.0);
  for (size_t k = 0; k < nnz; ++k) an[m.ti[k]] += m.tv[k] * st.a[m.tj[k]];
  st.a.swap(an);
  Matrix TP(rd, rd, 0.0);
  for (size_t k = 0; k < nnz; ++k) {
    const double v = m.tv[k];
    const double* src = &st.P.a[m.tj[k] * rd];
    double* dst = &TP.a[m.ti[k] * rd];
    for (size_t j = 0; j < rd; ++j) dst[j] += v * src[j];
  }
  Matrix Pn = m.V;
  for (size_t k = 0; k < nnz; ++k) {
    const double v = m.tv[k];
    const size_t ci = m.ti[k], cj = m.tj[k];
    for (size_t i = 0; i < rd; ++i) Pn.a[i * rd + ci] += v * TP.a[i * rd + cj];
  }
  st.P = Pn;
}

struct FilterResult {
  double ssq = 0.0, sumlog = 0.0;
  int nu = 0;
  Vec resid;
  FilterState state;
};

FilterResult kalman_filter(const SSModel& m, const Vec& y, bool want_resid) {
  const size_t n = y.size(), rd = m.rd;
  FilterResult out;
  FilterState& st = out.state;
  st.a.assign(rd, 0.0);
  st.P = m.P0;
  if (want_resid) out.resid.assign(n, 0.0);
  Vec PZ(rd);
  for (size_t t = 0; t < n; ++t) {
    if (t > 0) predict_step(m, st);
    double pred = 0.0;
    for (size_t i = 0; i < rd; ++i) pred += m.Z[i] * st.a[i];
    for (size_t i = 0; i < rd; ++i) {
      double s = 0.0;
      const double* row = &st.P.a[i * rd];
      for (size_t j = 0; j < rd; ++j) s += row[j] * m.Z[j];
      PZ[i] = s;
    }
    double F = 0.0;
    for (size_t i = 0; i < rd; ++i) F += m.Z[i] * PZ[i];
    double v = y[t] - pred;
    if (!(F > 0.0)) throw Error("arima: non-positive prediction variance in the Kalman filter");
    for (size_t i = 0; i < rd; ++i) st.a[i] += PZ[i] * v / F;
    for (size_t i = 0; i < rd; ++i)
      for (size_t j = 0; j < rd; ++j) st.P.a[i * rd + j] -= PZ[i] * PZ[j] / F;
    if (F < kDiffuseLimit) {
      out.ssq += v * v / F;
      out.sumlog += std::log(F);
      ++out.nu;
    }
    if (want_resid) out.resid[t] = v / std::sqrt(F);
  }
  return out;
}

struct Problem {
  Vec x;
  ArimaSpec spec;
  Layout L;
  Vec delta;
  int n_used = 0;
  int ncond = 0;
};

Vec centered(const Problem& pr, const Vec& coef) {
  Vec y = pr.x;
  if (pr.L.mean) {
    double mu = coef[static_cast<size_t>(pr.L.ncoef - 1)];
    for (double& v : y) v -= mu;
  }
  return y;
}

// ML objective in the original parameterisation: 0.5 log(s2) + 0.5 sum log F / nu.
double ml_objective(const Problem& pr, const Vec& coef, FilterResult* keep = nullptr) {
  Vec phi, theta;
  expand(pr.L, coef, phi, theta);
  if (!ar_stationary(phi)) return kInf;
  FilterResult fr;
  try {
    SSModel m = make_model(phi, theta, pr.delta);
    fr = kalman_filter(m, centered(pr, coef), keep != nullptr);
  } catch (const Error&) {
    return kInf;  // degenerate parameter values: treated as an infeasible point by the optimiser
  }
  if (fr.nu == 0) return kInf;
  double s2 = fr.ssq / fr.nu;
  double val = 0.5 * std::log(s2) + 0.5 * fr.sumlog / fr.nu;
  if (keep) *keep = std::move(fr);
  return std::isfinite(val) ? val : kInf;
}

// Conditional sum of squares; returns 0.5 log(ssq / nu) and optionally residuals, ssq, nu.
double css_objective(const Problem& pr, const Vec& coef, Vec* resid = nullptr, double* ssq_out = nullptr,
                     int* nu_out = nullptr) {
  Vec phi, theta;
  expand(pr.L, coef, phi, theta);
  Vec w = centered(pr, coef);
  const size_t n = w.size(), nd = pr.delta.size();
  for (size_t t = n; t-- > nd;) {
    // difference from the end so that lagged values are still the undifferenced ones
    double v = w[t];
    for (size_t i = 0; i < nd; ++i) v -= pr.delta[i] * w[t - i - 1];
    w[t] = v;
  }
  const size_t nc = static_cast<size_t>(pr.ncond);
  Vec e(n, 0.0);
  double ssq = 0.0;
  int nu = 0;
  for (size_t t = nc; t < n; ++t) {
    double v = w[t];
    for (size_t i = 0; i < phi.size(); ++i) v -= phi[i] * w[t - i - 1];
    for (size_t j = 0; j < theta.size() && j < t; ++j) v -= theta[j] * e[t - j - 1];
    e[t] = v;
    ssq += v * v;
    ++nu;
  }
  if (resid) *resid = e;
  if (ssq_out) *ssq_out = ssq;
  if (nu_out) *nu_out = nu;
  double val = 0.5 * std::log(ssq / nu);
  return std::isfinite(val) ? val : kInf;
}

// Original -> optimisation parameterisation (AR and seasonal AR parts via atanh of partial autocorrelations).
Vec to_free(const Layout& L, const Vec& c) {
  Vec u = c;
  auto conv = [&](size_t off, int len) {
    if (len == 0) return;
    Vec phi(c.begin() + static_cast<long>(off), c.begin() + static_cast<long>(off) + len), r;
    if (!ar_to_pacf(phi, r)) throw Error("arima: non-stationary AR part");
    for (int k = 0; k < len; ++k) u[off + static_cast<size_t>(k)] = std::atanh(r[static_cast<size_t>(k)]);
  };
  conv(0, L.p);
  conv(static_cast<size_t>(L.p + L.q), L.P);
  return u;
}

Vec from_free(const Layout& L, const Vec& u) {
  Vec c = u;
  auto conv = [&](size_t off, int len) {
    if (len == 0) return;
    Vec v(u.begin() + static_cast<long>(off), u.begin() + static_cast<long>(off) + len);
    Vec phi = pacf_to_ar(v);
    for (int k = 0; k < len; ++k) c[off + static_cast<size_t>(k)] = phi[static_cast<size_t>(k)];
  };
  conv(0, L.p);
  conv(static_cast<size_t>(L.p + L.q), L.P);
  return c;
}

// Invert a symmetric matrix by Gauss-Jordan elimination; false if singular.
bool invert(Matrix& A) {
  const size_t m = A.rows;
  Matrix I(m, m, 0.0);
  for (size_t i = 0; i < m; ++i) I(i, i) = 1.0;
  for (size_t k = 0; k < m; ++k) {
    size_t piv = k;
    for (size_t i = k + 1; i < m; ++i)
      if (std::fabs(A(i, k)) > std::fabs(A(piv, k))) piv = i;
    if (A(piv, k) == 0.0 || !std::isfinite(A(piv, k))) return false;
    for (size_t j = 0; j < m; ++j) { std::swap(A(k, j), A(piv, j)); std::swap(I(k, j), I(piv, j)); }
    double inv = 1.0 / A(k, k);
    for (size_t j = 0; j < m; ++j) { A(k, j) *= inv; I(k, j) *= inv; }
    for (size_t i = 0; i < m; ++i) {
      if (i == k) continue;
      double f = A(i, k);
      if (f == 0.0) continue;
      for (size_t j = 0; j < m; ++j) { A(i, j) -= f * A(k, j); I(i, j) -= f * I(k, j); }
    }
  }
  A = I;
  return true;
}

void validate(const Vec& x, const ArimaSpec& sp) {
  if (sp.p < 0 || sp.d < 0 || sp.q < 0 || sp.P < 0 || sp.D < 0 || sp.Q < 0)
    throw Error("arima: 'order' and 'seasonal$order' must be non-negative");
  if (sp.period < 1) throw Error("arima: 'period' must be a positive integer");
  if ((sp.P > 0 || sp.D > 0 || sp.Q > 0) && sp.period < 2)
    throw Error("arima: a seasonal part needs 'period' > 1");
  if (sp.method != "CSS-ML" && sp.method != "ML" && sp.method != "CSS")
    throw Error("arima: 'method' must be one of \"CSS-ML\", \"ML\", \"CSS\"");
  for (double v : x)
    if (!std::isfinite(v)) throw Error("arima: 'x' contains missing or non-finite values");
}

Problem make_problem(const Vec& x, const ArimaSpec& spec) {
  Problem pr;
  pr.x = x;
  pr.spec = spec;
  pr.L = layout_of(spec);
  pr.delta = delta_of(spec.d, spec.D, spec.period);
  pr.n_used = static_cast<int>(x.size()) - static_cast<int>(pr.delta.size());
  pr.ncond = static_cast<int>(pr.delta.size()) + spec.p + spec.period * spec.P;
  return pr;
}

}  // namespace

ArimaFit arima(const Vec& x, const ArimaSpec& spec_in) {
  ArimaSpec spec = spec_in;
  if (spec.P == 0 && spec.D == 0 && spec.Q == 0 && spec.period < 1) spec.period = 1;
  validate(x, spec);
  Problem pr = make_problem(x, spec);
  const Layout& L = pr.L;
  const size_t k = static_cast<size_t>(L.ncoef);
  if (pr.n_used <= 0 || pr.n_used <= pr.ncond - static_cast<int>(pr.delta.size()) ||
      static_cast<int>(x.size()) - pr.ncond < static_cast<int>(k) + 1)
    throw Error("arima: too few non-missing observations");

  ArimaFit fit;
  fit.spec = spec;
  fit.x = x;
  fit.nobs = pr.n_used;
  for (int i = 0; i < L.p; ++i) fit.coef_names.push_back("ar" + std::to_string(i + 1));
  for (int i = 0; i < L.q; ++i) fit.coef_names.push_back("ma" + std::to_string(i + 1));
  for (int i = 0; i < L.P; ++i) fit.coef_names.push_back("sar" + std::to_string(i + 1));
  for (int i = 0; i < L.Q; ++i) fit.coef_names.push_back("sma" + std::to_string(i + 1));
  if (L.mean) fit.coef_names.push_back("intercept");

  Vec init(k, 0.0);
  if (L.mean) {
    double s = 0.0;
    for (double v : x) s += v;
    init[k - 1] = s / static_cast<double>(x.size());
  }
  OptimOptions oo;
  oo.method = "BFGS";
  // The intercept is optimised on its own scale (optim's parscale idea): mu = scale * z with scale = 10 times the
  // iid standard error of the mean, so that its finite-difference steps are commensurate with those of the ARMA
  // parameters. Hessians are taken in unscaled coordinates.
  double mscale = 1.0;
  if (L.mean) {
    double m = init[k - 1], ss = 0.0;
    for (double v : x) ss += (v - m) * (v - m);
    double sd = std::sqrt(ss / static_cast<double>(x.size() - 1));
    if (sd > 0.0 && std::isfinite(sd)) mscale = 10.0 * sd / std::sqrt(static_cast<double>(x.size()));
  }
  auto unscale = [&](Vec z) {
    if (L.mean) z[k - 1] *= mscale;
    return z;
  };
  auto rescale = [&](Vec u) {
    if (L.mean) u[k - 1] /= mscale;
    return u;
  };
  Vec coef = init;
  const std::string& method = spec.method;
  // vcov = A (n_used * H)^{-1} A', H the numerical Hessian of f at `at` and A the Jacobian of the map from the
  // optimisation parameters to the coefficients (identity unless the Jones transformation was used).
  auto hessian_vcov = [&](const FnN& f, const Vec& at, const Matrix* A) {
    fit.vcov = Matrix(k, k, kNaN);
    fit.se.assign(k, kNaN);
    if (k == 0) return;
    Vec h = numeric_hessian(f, at, oo.ndeps);
    Matrix H(k, k);
    for (size_t i = 0; i < k * k; ++i) H.a[i] = h[i] * pr.n_used;
    if (!invert(H)) return;
    if (A) {
      Matrix AH(k, k, 0.0), V(k, k, 0.0);
      for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < k; ++j)
          for (size_t l = 0; l < k; ++l) AH(i, j) += (*A)(i, l) * H(l, j);
      for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < k; ++j)
          for (size_t l = 0; l < k; ++l) V(i, j) += AH(i, l) * (*A)(j, l);
      H = V;
    }
    fit.vcov = H;
    for (size_t i = 0; i < k; ++i) fit.se[i] = H(i, i) >= 0.0 ? std::sqrt(H(i, i)) : kNaN;
  };

  if (method == "CSS" || method == "CSS-ML") {
    FnN fcss = [&](const Vec& c) { return css_objective(pr, c); };
    if (k > 0) {
      if (!std::isfinite(fcss(coef))) throw Error("arima: non-finite value supplied by optim");
      FnN fz = [&](const Vec& z) { return fcss(unscale(z)); };
      OptimResult res = optim(fz, rescale(coef), oo);
      coef = unscale(res.par);
      fit.convergence = res.convergence;
    }
    if (method == "CSS") {
      double ssq = 0.0;
      int nu = 0;
      double val = css_objective(pr, coef, &fit.residuals, &ssq, &nu);
      fit.coef = coef;
      fit.sigma2 = ssq / nu;
      fit.loglik = -0.5 * (2.0 * pr.n_used * val + pr.n_used + pr.n_used * kLog2Pi);
      fit.aic = kNaN;
      hessian_vcov(fcss, coef, nullptr);
      return fit;
    }
    Vec phi(coef.begin(), coef.begin() + L.p), sphi(coef.begin() + L.p + L.q, coef.begin() + L.p + L.q + L.P);
    if (!ar_stationary(phi) || !ar_stationary(sphi)) throw Error("arima: non-stationary AR part from CSS");
  }

  // ML with the Jones reparameterisation
  FnN fml = [&](const Vec& c) { return ml_objective(pr, c); };
  FnN ftr = [&](const Vec& u) { return ml_objective(pr, from_free(L, u)); };
  Vec uhat;
  bool inverted = false;
  if (k > 0) {
    Vec u0 = to_free(L, coef);
    if (!std::isfinite(ftr(u0))) throw Error("arima: non-finite value supplied by optim");
    FnN fz = [&](const Vec& z) { return ftr(unscale(z)); };
    OptimResult res = optim(fz, rescale(u0), oo);
    fit.convergence = res.convergence;
    uhat = unscale(res.par);
    coef = from_free(L, uhat);
    // invertible MA parts
    auto invert_part = [&](size_t off, int len) {
      if (len == 0) return;
      Vec part(coef.begin() + static_cast<long>(off), coef.begin() + static_cast<long>(off) + len);
      Vec inv = ma_invert(part);
      if (inv != part) inverted = true;
      std::copy(inv.begin(), inv.end(), coef.begin() + static_cast<long>(off));
    };
    invert_part(static_cast<size_t>(L.p), L.q);
    invert_part(static_cast<size_t>(L.p + L.q + L.P), L.Q);
  }
  FilterResult fr;
  double val = ml_objective(pr, coef, &fr);
  if (!std::isfinite(val)) throw Error("arima: non-finite likelihood at the estimate");
  fit.coef = coef;
  fit.sigma2 = fr.ssq / fr.nu;
  fit.residuals = fr.resid;
  fit.loglik = -0.5 * (2.0 * pr.n_used * val + pr.n_used + pr.n_used * kLog2Pi);
  fit.aic = -2.0 * fit.loglik + 2.0 * (static_cast<double>(k) + 1.0);
  if (k > 0 && L.p + L.P > 0 && !inverted) {
    // Hessian in the transformed parameters at the optimum, mapped back with the forward-difference Jacobian of
    // the transformation (step ndeps); this is what R's var.coef matches when transform.pars = TRUE.
    Matrix A(k, k, 0.0);
    const Vec c0 = from_free(L, uhat);
    for (size_t j = 0; j < k; ++j) {
      Vec up = uhat;
      up[j] += oo.ndeps;
      Vec cp = from_free(L, up);
      for (size_t i = 0; i < k; ++i) A(i, j) = (cp[i] - c0[i]) / oo.ndeps;
    }
    hessian_vcov(ftr, uhat, &A);
  } else {
    hessian_vcov(fml, coef, nullptr);
  }
  return fit;
}

ArimaForecast ArimaFit::forecast(int n_ahead) const {
  if (n_ahead < 1) throw Error("predict.Arima: 'n.ahead' must be at least 1");
  if (coef.size() != static_cast<size_t>(layout_of(spec).ncoef) || x.empty())
    throw Error("predict.Arima: the fit is incomplete");
  Problem pr = make_problem(x, spec);
  Vec phi, theta;
  expand(pr.L, coef, phi, theta);
  SSModel m = make_model(phi, theta, pr.delta);
  FilterResult fr = kalman_filter(m, centered(pr, coef), false);
  FilterState st = fr.state;
  const double mu = pr.L.mean ? coef.back() : 0.0;
  ArimaForecast out;
  out.mean.resize(static_cast<size_t>(n_ahead));
  out.se.resize(static_cast<size_t>(n_ahead));
  for (int h = 0; h < n_ahead; ++h) {
    predict_step(m, st);
    double pred = 0.0, F = 0.0;
    for (size_t i = 0; i < m.rd; ++i) {
      pred += m.Z[i] * st.a[i];
      double s = 0.0;
      for (size_t j = 0; j < m.rd; ++j) s += st.P(i, j) * m.Z[j];
      F += m.Z[i] * s;
    }
    out.mean[static_cast<size_t>(h)] = pred + mu;
    out.se[static_cast<size_t>(h)] = std::sqrt(F * sigma2);
  }
  return out;
}

Vec arima_sim(Rng& rng, int n, const Vec& ar, const Vec& ma, double sd, int d, int n_start) {
  if (n < 1) throw Error("arima.sim: 'n' must be a positive integer");
  if (d < 0) throw Error("arima.sim: 'd' must be non-negative");
  if (!(sd >= 0.0) || !std::isfinite(sd)) throw Error("arima.sim: 'sd' must be a non-negative number");
  for (double v : ar)
    if (!std::isfinite(v)) throw Error("arima.sim: 'ar' contains non-finite values");
  for (double v : ma)
    if (!std::isfinite(v)) throw Error("arima.sim: 'ma' contains non-finite values");
  size_t p = ar.size();
  while (p > 0 && ar[p - 1] == 0.0) --p;
  const size_t q = ma.size();
  double minroots = kInf;
  if (p > 0) {
    Vec c(p + 1);
    c[0] = 1.0;
    for (size_t i = 0; i < p; ++i) c[i + 1] = -ar[i];
    for (const auto& z : poly_roots(c)) minroots = std::min(minroots, std::abs(z));
    if (!(minroots > 1.0) || !ar_stationary(Vec(ar.begin(), ar.begin() + static_cast<long>(p))))
      throw Error("arima.sim: 'ar' part of model is not stationary");
  }
  if (n_start < 0) {
    n_start = static_cast<int>(ar.size() + q);
    if (p > 0) n_start += static_cast<int>(std::ceil(6.0 / std::log(minroots)));
  } else if (static_cast<size_t>(n_start) < ar.size() + q) {
    throw Error("arima.sim: burn-in 'n.start' must be as long as 'ar + ma'");
  }
  const size_t N = static_cast<size_t>(n_start) + static_cast<size_t>(n);
  Vec e(N);
  for (size_t i = 0; i < N; ++i) e[i] = sd * rng.norm();
  Vec x(N);
  for (size_t t = 0; t < N; ++t) {
    double v = e[t];
    if (t < q) v = 0.0;  // the MA filter is undefined for the first q values; they are set to zero
    else
      for (size_t j = 0; j < q; ++j) v += ma[j] * e[t - j - 1];
    x[t] = v;
  }
  for (size_t t = 0; t < N; ++t)
    for (size_t i = 0; i < p && i < t; ++i) x[t] += ar[i] * x[t - i - 1];
  Vec out(x.begin() + n_start, x.end());
  for (int k = 0; k < d; ++k) {
    Vec y(out.size() + 1, 0.0);
    for (size_t i = 0; i < out.size(); ++i) y[i + 1] = y[i] + out[i];
    out.swap(y);
  }
  return out;
}

}  // namespace rcli
