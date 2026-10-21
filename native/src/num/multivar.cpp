// Multivariate analysis: symmetric eigendecomposition, SVD, PCA, covariance/correlation, Cholesky, Mahalanobis
// distances, classical multidimensional scaling.
//
// Clean-room implementation from published descriptions:
//  - Golub, G. H. and Van Loan, C. F. (2013). Matrix Computations, 4th ed. Johns Hopkins.
//    Householder vectors (Alg. 5.1.1), Householder tridiagonalisation (Alg. 8.3.1), implicit symmetric QR step
//    with Wilkinson shift and the symmetric QR algorithm with deflation (Alg. 8.3.2, 8.3.3), Cholesky (Sec. 4.2).
//  - Hestenes, M. R. (1958). Inversion of matrices by biorthogonalization and related results. J. SIAM 6, 51-90.
//    Demmel, J. and Veselic, K. (1992). Jacobi's method is more accurate than QR. SIAM J. Matrix Anal. Appl. 13,
//    1204-1245. (One-sided Jacobi SVD: rotate column pairs until mutually orthogonal.)
//  - Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis.
//    Biometrika 53, 325-328; Torgerson, W. S. (1958). Theory and Methods of Scaling. (Principal coordinates.)
//  - Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979). Multivariate Analysis. (PCA, Mahalanobis distance.)
// Behavioural details (prcomp's divisor, scale() for uncentred data, error messages, cmdscale's k range and GOF) were
// taken from R's help pages and black-box comparison with R 4.3.3.
#include "rcli/multivar.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "rcli/cluster.hpp"

namespace rcli {

namespace {

void check_matrix(const Matrix& x, const char* fn, bool allow_empty = false) {
  if (x.a.size() != x.rows * x.cols) throw Error(std::string(fn) + ": matrix storage size does not match dimensions");
  if (!allow_empty && (x.rows == 0 || x.cols == 0)) throw Error(std::string(fn) + ": empty matrix");
  for (double v : x.a)
    if (!std::isfinite(v)) throw Error(std::string(fn) + ": infinite or missing values in 'x'");
}

double sgn1(double v) { return v < 0.0 ? -1.0 : 1.0; }

// Index of the component that decides a vector's sign: the largest magnitude, ties (within a relative 1e-9) going
// to the lowest index.
size_t sign_index(const double* v, size_t n, size_t stride) {
  double mx = 0.0;
  for (size_t i = 0; i < n; ++i) mx = std::max(mx, std::fabs(v[i * stride]));
  for (size_t i = 0; i < n; ++i)
    if (std::fabs(v[i * stride]) >= mx * (1.0 - 1e-9)) return i;
  return 0;
}

// Column means with a second correction pass (mean of residuals added), so a constant column has exactly its value.
Vec col_means(const Matrix& x) {
  const size_t n = x.rows, p = x.cols;
  Vec m(p, 0.0);
  for (size_t j = 0; j < p; ++j) {
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += x(i, j);
    double mu = s / static_cast<double>(n);
    if (std::isfinite(mu)) {
      double t = 0.0;
      for (size_t i = 0; i < n; ++i) t += x(i, j) - mu;
      mu += t / static_cast<double>(n);
    }
    m[j] = mu;
  }
  return m;
}

// Householder vector (GvL Alg. 5.1.1) for x[0..m): v[0] = 1 and (I - beta v v') x = mu e1, mu = ||x|| >= 0.
// Returns mu; beta = 0 when x is already a non-negative multiple of e1.
double house(const double* x, size_t m, Vec& v, double& beta) {
  v.assign(x, x + m);
  double scale = 0.0;
  for (size_t i = 0; i < m; ++i) scale = std::max(scale, std::fabs(v[i]));
  if (scale == 0.0) { beta = 0.0; return 0.0; }
  for (double& t : v) t /= scale;
  double sigma = 0.0;
  for (size_t i = 1; i < m; ++i) sigma += v[i] * v[i];
  const double x0 = v[0];
  if (sigma == 0.0) {
    if (x0 >= 0.0) { beta = 0.0; v[0] = 1.0; return x0 * scale; }
    // Reflect the sign: H = I - 2 e1 e1'.
    beta = 2.0;
    v[0] = 1.0;
    return -x0 * scale;
  }
  const double mu = std::sqrt(x0 * x0 + sigma);
  const double v0 = x0 <= 0.0 ? x0 - mu : -sigma / (x0 + mu);
  beta = 2.0 * v0 * v0 / (sigma + v0 * v0);
  for (size_t i = 1; i < m; ++i) v[i] /= v0;
  v[0] = 1.0;
  return mu * scale;
}

// Symmetric eigenproblem for a symmetric n x n matrix stored row-major in `a` (full storage, symmetric).
// Returns eigenvalues (unsorted) and q (row-major n x n, columns eigenvectors).
void sym_eigen_raw(Vec a, size_t n, Vec& d, Vec& q) {
  d.assign(n, 0.0);
  Vec e(n > 0 ? n - 1 : 0, 0.0);
  q.assign(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) q[i * n + i] = 1.0;
  if (n == 0) return;
  auto A = [&](size_t i, size_t j) -> double& { return a[i * n + j]; };
  // Householder tridiagonalisation (GvL Alg. 8.3.1): A = Q T Q'.
  std::vector<Vec> vs(n >= 2 ? n - 2 : 0);
  Vec betas(vs.size(), 0.0);
  Vec x, p(n), w(n);
  for (size_t k = 0; k + 2 < n; ++k) {
    const size_t m = n - k - 1;
    x.resize(m);
    for (size_t i = 0; i < m; ++i) x[i] = A(k + 1 + i, k);
    double beta;
    const double mu = house(x.data(), m, vs[k], beta);
    betas[k] = beta;
    const Vec& v = vs[k];
    if (beta != 0.0) {
      for (size_t i = 0; i < m; ++i) {
        double s = 0.0;
        for (size_t j = 0; j < m; ++j) s += A(k + 1 + i, k + 1 + j) * v[j];
        p[i] = beta * s;
      }
      double pv = 0.0;
      for (size_t i = 0; i < m; ++i) pv += p[i] * v[i];
      const double f = beta * pv / 2.0;
      for (size_t i = 0; i < m; ++i) w[i] = p[i] - f * v[i];
      for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < m; ++j) A(k + 1 + i, k + 1 + j) -= v[i] * w[j] + w[i] * v[j];
    }
    e[k] = mu;
    d[k] = A(k, k);
  }
  if (n >= 2) {
    d[n - 2] = A(n - 2, n - 2);
    e[n - 2] = A(n - 1, n - 2);
  }
  d[n - 1] = A(n - 1, n - 1);
  // Backward accumulation Q = H_0 H_1 ... H_{n-3}.
  for (size_t kk = vs.size(); kk-- > 0;) {
    if (betas[kk] == 0.0) continue;
    const Vec& v = vs[kk];
    const size_t off = kk + 1, m = n - off;
    for (size_t j = off; j < n; ++j) {
      double s = 0.0;
      for (size_t i = 0; i < m; ++i) s += v[i] * q[(off + i) * n + j];
      s *= betas[kk];
      for (size_t i = 0; i < m; ++i) q[(off + i) * n + j] -= s * v[i];
    }
  }
  // Symmetric QR iteration with Wilkinson shift and deflation (GvL Alg. 8.3.2, 8.3.3).
  const size_t max_steps = 60 * n + 60;
  size_t steps = 0;
  while (true) {
    for (size_t i = 0; i + 1 < n; ++i)
      if (std::fabs(e[i]) <= kEps * (std::fabs(d[i]) + std::fabs(d[i + 1]))) e[i] = 0.0;
    size_t hi = n - 1;
    while (hi > 0 && e[hi - 1] == 0.0) --hi;
    if (hi == 0) break;
    size_t lo = hi - 1;
    while (lo > 0 && e[lo - 1] != 0.0) --lo;
    if (++steps > max_steps) throw Error("eigen_sym: QR iteration did not converge");
    // Wilkinson shift from the trailing 2 x 2 block of the unreduced block [lo, hi].
    const double dd = (d[hi - 1] - d[hi]) / 2.0;
    const double eh = e[hi - 1];
    const double den = dd + sgn1(dd) * std::hypot(dd, eh);
    const double mu = d[hi] - eh * (eh / den);
    double xx = d[lo] - mu, zz = e[lo];
    for (size_t k = lo; k < hi; ++k) {
      const double r = std::hypot(xx, zz);
      double c = 1.0, s = 0.0;
      if (r != 0.0) { c = xx / r; s = zz / r; }
      if (k > lo) e[k - 1] = r;
      const double a0 = d[k], b0 = e[k], c0 = d[k + 1];
      d[k] = c * c * a0 + 2.0 * c * s * b0 + s * s * c0;
      d[k + 1] = s * s * a0 - 2.0 * c * s * b0 + c * c * c0;
      e[k] = c * s * (c0 - a0) + (c * c - s * s) * b0;
      if (k + 1 < hi) {
        zz = s * e[k + 1];
        e[k + 1] *= c;
        xx = e[k];
      }
      for (size_t i = 0; i < n; ++i) {
        const double qk = q[i * n + k], qk1 = q[i * n + k + 1];
        q[i * n + k] = c * qk + s * qk1;
        q[i * n + k + 1] = -s * qk + c * qk1;
      }
    }
  }
}

// One-sided Jacobi SVD of an n x p matrix with n >= p, given column-major columns w[j] (length n).
// On return w[j] = sigma_j u_j (unnormalised), v column-major p x p.
void jacobi_svd_raw(std::vector<Vec>& w, size_t n, std::vector<Vec>& v) {
  const size_t p = w.size();
  v.assign(p, Vec(p, 0.0));
  for (size_t j = 0; j < p; ++j) v[j][j] = 1.0;
  for (int sweep = 0;; ++sweep) {
    if (sweep >= 80) throw Error("svd: one-sided Jacobi did not converge");
    bool rotated = false;
    for (size_t i = 0; i + 1 < p; ++i) {
      for (size_t j = i + 1; j < p; ++j) {
        double alpha = 0.0, beta = 0.0, gamma = 0.0;
        const double* wi = w[i].data();
        const double* wj = w[j].data();
        for (size_t r = 0; r < n; ++r) {
          alpha += wi[r] * wi[r];
          beta += wj[r] * wj[r];
          gamma += wi[r] * wj[r];
        }
        if (gamma == 0.0 || std::fabs(gamma) <= kEps * std::sqrt(alpha) * std::sqrt(beta)) continue;
        rotated = true;
        const double zeta = (beta - alpha) / (2.0 * gamma);
        double t;
        if (std::fabs(zeta) > 1e150) t = 1.0 / (2.0 * zeta);
        else t = sgn1(zeta) / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
        const double c = 1.0 / std::sqrt(1.0 + t * t), s = c * t;
        for (size_t r = 0; r < n; ++r) {
          const double a = w[i][r], b = w[j][r];
          w[i][r] = c * a - s * b;
          w[j][r] = s * a + c * b;
        }
        for (size_t r = 0; r < p; ++r) {
          const double a = v[i][r], b = v[j][r];
          v[i][r] = c * a - s * b;
          v[j][r] = s * a + c * b;
        }
      }
    }
    if (!rotated) break;
  }
}

double dot(const Vec& a, const Vec& b) {
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Completes unit columns u[j] for which ok[j] is false so that all columns are orthonormal: the coordinate vector
// with the largest residual after (twice repeated) Gram-Schmidt against the columns already set.
void complete_basis(std::vector<Vec>& u, std::vector<bool>& ok, size_t n) {
  for (size_t j = 0; j < u.size(); ++j) {
    if (ok[j]) continue;
    Vec best;
    double best_norm = -1.0;
    for (size_t t = 0; t < n; ++t) {
      Vec r(n, 0.0);
      r[t] = 1.0;
      for (int pass = 0; pass < 2; ++pass)
        for (size_t i = 0; i < u.size(); ++i) {
          if (!ok[i]) continue;
          const double c = dot(u[i], r);
          for (size_t l = 0; l < n; ++l) r[l] -= c * u[i][l];
        }
      const double nr = std::sqrt(dot(r, r));
      if (nr > best_norm * (1.0 + 1e-9)) { best_norm = nr; best = r; }
    }
    for (double& t : best) t /= best_norm;
    u[j] = best;
    ok[j] = true;
  }
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------

EigenSym eigen_sym(const Matrix& a) {
  check_matrix(a, "eigen_sym");
  if (a.rows != a.cols) throw Error("eigen_sym: matrix is not square");
  const size_t n = a.rows;
  double scale = 0.0;
  for (double v : a.a) scale = std::max(scale, std::fabs(v));
  if (scale == 0.0) scale = 1.0;
  Vec s(n * n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= i; ++j) s[i * n + j] = s[j * n + i] = a(i, j) / scale;
  Vec d, q;
  sym_eigen_raw(std::move(s), n, d, q);
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), size_t{0});
  std::stable_sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return d[x] > d[y]; });
  EigenSym r;
  r.values.resize(n);
  r.vectors = Matrix(n, n);
  for (size_t c = 0; c < n; ++c) {
    const size_t k = ord[c];
    r.values[c] = d[k] * scale;
    double nrm = 0.0;
    for (size_t i = 0; i < n; ++i) nrm += q[i * n + k] * q[i * n + k];
    nrm = std::sqrt(nrm);
    const size_t si = sign_index(&q[k], n, n);
    const double f = (q[si * n + k] < 0.0 ? -1.0 : 1.0) / nrm;
    for (size_t i = 0; i < n; ++i) r.vectors(i, c) = q[i * n + k] * f;
  }
  return r;
}

Svd svd(const Matrix& x) {
  check_matrix(x, "svd");
  const bool tr = x.rows < x.cols;
  const size_t n = tr ? x.cols : x.rows, p = tr ? x.rows : x.cols;  // working matrix n x p, n >= p
  double scale = 0.0;
  for (double v : x.a) scale = std::max(scale, std::fabs(v));
  if (scale == 0.0) scale = 1.0;
  std::vector<Vec> w(p, Vec(n));
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < n; ++i) w[j][i] = (tr ? x(j, i) : x(i, j)) / scale;
  std::vector<Vec> v;
  jacobi_svd_raw(w, n, v);
  Vec d(p);
  std::vector<bool> ok(p, true);
  for (size_t j = 0; j < p; ++j) {
    d[j] = std::sqrt(dot(w[j], w[j]));
    if (!(d[j] >= std::numeric_limits<double>::min())) {
      d[j] = 0.0;
      ok[j] = false;
    } else {
      for (double& t : w[j]) t /= d[j];
    }
  }
  complete_basis(w, ok, n);
  std::vector<size_t> ord(p);
  std::iota(ord.begin(), ord.end(), size_t{0});
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return d[a] > d[b]; });
  // Map back: x = U d V' (not transposed) or x' = U d V' => x = V d U'.
  const size_t rn = x.rows, rp = x.cols;
  Svd r;
  r.d.resize(p);
  r.u = Matrix(rn, p);
  r.v = Matrix(rp, p);
  for (size_t c = 0; c < p; ++c) {
    const size_t k = ord[c];
    r.d[c] = d[k] * scale;
    const Vec& left = tr ? v[k] : w[k];   // length rn
    const Vec& right = tr ? w[k] : v[k];  // length rp
    const size_t si = sign_index(right.data(), rp, 1);
    const double f = right[si] < 0.0 ? -1.0 : 1.0;
    for (size_t i = 0; i < rn; ++i) r.u(i, c) = f * left[i];
    for (size_t i = 0; i < rp; ++i) r.v(i, c) = f * right[i];
  }
  return r;
}

Pca prcomp(const Matrix& x, bool center, bool scale) {
  check_matrix(x, "prcomp");
  const size_t n = x.rows, p = x.cols;
  Matrix z = x;
  Pca out;
  if (center) {
    out.center = col_means(x);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < p; ++j) z(i, j) -= out.center[j];
  }
  if (scale) {
    out.scale.assign(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
      double ss = 0.0;
      bool constant = true;
      for (size_t i = 0; i < n; ++i) {
        ss += z(i, j) * z(i, j);
        if (z(i, j) != z(0, j)) constant = false;
      }
      const double sc = n >= 2 ? std::sqrt(ss / static_cast<double>(n - 1)) : kNaN;
      if (!(sc > 0.0) || !std::isfinite(sc) || (center && constant))
        throw Error("prcomp: cannot rescale a constant/zero column to unit variance");
      out.scale[j] = sc;
      for (size_t i = 0; i < n; ++i) z(i, j) /= sc;
    }
  }
  Svd s = svd(z);
  const size_t k = s.d.size();
  const double div = std::sqrt(static_cast<double>(std::max<size_t>(1, n - 1)));
  out.sdev.resize(k);
  for (size_t c = 0; c < k; ++c) out.sdev[c] = s.d[c] / div;
  out.rotation = s.v;
  out.x = Matrix(n, k);
  for (size_t i = 0; i < n; ++i)
    for (size_t c = 0; c < k; ++c) {
      double t = 0.0;
      for (size_t j = 0; j < p; ++j) t += z(i, j) * s.v(j, c);
      out.x(i, c) = t;
    }
  double tot = 0.0;
  for (double sd : out.sdev) tot += sd * sd;
  out.prop_var.resize(k);
  out.cum_prop.resize(k);
  double cum = 0.0;
  for (size_t c = 0; c < k; ++c) {
    out.prop_var[c] = tot > 0.0 ? out.sdev[c] * out.sdev[c] / tot : kNaN;
    cum += out.sdev[c] * out.sdev[c];
    out.cum_prop[c] = tot > 0.0 ? cum / tot : kNaN;
  }
  return out;
}

Matrix cov_matrix(const Matrix& x) {
  check_matrix(x, "cov_matrix");
  const size_t n = x.rows, p = x.cols;
  Matrix c(p, p, kNaN);
  if (n < 2) return c;
  const Vec m = col_means(x);
  for (size_t j = 0; j < p; ++j)
    for (size_t l = 0; l <= j; ++l) {
      double s = 0.0;
      for (size_t i = 0; i < n; ++i) s += (x(i, j) - m[j]) * (x(i, l) - m[l]);
      c(j, l) = c(l, j) = s / static_cast<double>(n - 1);
    }
  return c;
}

Matrix cor_matrix(const Matrix& x) {
  Matrix c = cov_matrix(x);
  const size_t p = c.rows;
  if (x.rows < 2) return c;
  Vec sd(p);
  for (size_t j = 0; j < p; ++j) sd[j] = std::sqrt(c(j, j));
  for (size_t j = 0; j < p; ++j)
    for (size_t l = 0; l < p; ++l) {
      if (!(sd[j] > 0.0) || !(sd[l] > 0.0)) c(j, l) = kNaN;
      else if (j == l) c(j, l) = 1.0;
      else c(j, l) = std::min(1.0, std::max(-1.0, c(j, l) / (sd[j] * sd[l])));
    }
  return c;
}

Matrix cholesky(const Matrix& a) {
  check_matrix(a, "chol");
  if (a.rows != a.cols) throw Error("chol: matrix is not square");
  const size_t n = a.rows;
  Matrix r(n, n, 0.0);
  for (size_t j = 0; j < n; ++j) {
    double s = a(j, j);
    for (size_t k = 0; k < j; ++k) s -= r(k, j) * r(k, j);
    if (!(s > 0.0)) throw Error("chol: the leading minor of order " + std::to_string(j + 1) + " is not positive");
    const double rjj = std::sqrt(s);
    r(j, j) = rjj;
    for (size_t i = j + 1; i < n; ++i) {
      double t = a(j, i);
      for (size_t k = 0; k < j; ++k) t -= r(k, j) * r(k, i);
      r(j, i) = t / rjj;
    }
  }
  return r;
}

Matrix chol_inverse(const Matrix& a) {
  const Matrix r = cholesky(a);
  const size_t n = r.rows;
  // ri = r^{-1} (upper triangular), by back substitution column by column.
  Matrix ri(n, n, 0.0);
  for (size_t j = 0; j < n; ++j) {
    ri(j, j) = 1.0 / r(j, j);
    for (size_t ii = j; ii-- > 0;) {
      double s = 0.0;
      for (size_t k = ii + 1; k <= j; ++k) s += r(ii, k) * ri(k, j);
      ri(ii, j) = -s / r(ii, ii);
    }
  }
  Matrix inv(n, n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i; j < n; ++j) {
      double s = 0.0;
      for (size_t k = j; k < n; ++k) s += ri(i, k) * ri(j, k);
      inv(i, j) = inv(j, i) = s;
    }
  return inv;
}

Vec chol_solve(const Matrix& a, const Vec& b) {
  const Matrix r = cholesky(a);
  const size_t n = r.rows;
  if (b.size() != n) throw Error("chol_solve: 'b' has the wrong length");
  for (double v : b)
    if (!std::isfinite(v)) throw Error("chol_solve: infinite or missing values in 'b'");
  Vec y(n);
  for (size_t i = 0; i < n; ++i) {  // r' y = b
    double s = b[i];
    for (size_t k = 0; k < i; ++k) s -= r(k, i) * y[k];
    y[i] = s / r(i, i);
  }
  for (size_t ii = n; ii-- > 0;) {  // r x = y
    double s = y[ii];
    for (size_t k = ii + 1; k < n; ++k) s -= r(ii, k) * y[k];
    y[ii] = s / r(ii, ii);
  }
  return y;
}

Vec mahalanobis(const Matrix& x, const Vec& center, const Matrix& cov, bool inverted) {
  check_matrix(x, "mahalanobis", true);
  check_matrix(cov, "mahalanobis");
  const size_t p = x.cols;
  if (center.size() != p) throw Error("mahalanobis: 'center' must have one value per column of 'x'");
  if (cov.rows != p || cov.cols != p) throw Error("mahalanobis: 'cov' must be a p x p matrix, p = ncol(x)");
  for (double v : center)
    if (!std::isfinite(v)) throw Error("mahalanobis: infinite or missing values in 'center'");
  Vec out(x.rows);
  Vec diff(p), z(p);
  if (inverted) {
    for (size_t i = 0; i < x.rows; ++i) {
      for (size_t j = 0; j < p; ++j) diff[j] = x(i, j) - center[j];
      double s = 0.0;
      for (size_t j = 0; j < p; ++j) {
        double t = 0.0;
        for (size_t l = 0; l < p; ++l) t += cov(j, l) * diff[l];
        s += diff[j] * t;
      }
      out[i] = s;
    }
    return out;
  }
  Matrix r;
  try {
    r = cholesky(cov);
  } catch (const Error&) {
    throw Error("mahalanobis: covariance matrix is not positive definite");
  }
  for (size_t i = 0; i < x.rows; ++i) {
    for (size_t j = 0; j < p; ++j) diff[j] = x(i, j) - center[j];
    double s = 0.0;
    for (size_t j = 0; j < p; ++j) {  // r' z = diff, distance = z'z
      double t = diff[j];
      for (size_t k = 0; k < j; ++k) t -= r(k, j) * z[k];
      z[j] = t / r(j, j);
      s += z[j] * z[j];
    }
    out[i] = s;
  }
  return out;
}

Cmdscale cmdscale(const Matrix& d, int k) {
  if (d.a.size() != d.rows * d.cols) throw Error("cmdscale: matrix storage size does not match dimensions");
  for (double v : d.a)
    if (std::isnan(v)) throw Error("cmdscale: NA values not allowed in 'd'");
  for (double v : d.a)
    if (!std::isfinite(v)) throw Error("cmdscale: infinite values in 'd'");
  if (d.rows != d.cols) throw Error("cmdscale: 'd' must be a square (distance) matrix");
  const size_t n = d.rows;
  if (k < 1 || n < 2 || static_cast<size_t>(k) > n - 1) throw Error("cmdscale: 'k' must be in {1, 2, ..  n - 1}");
  double mx = 0.0;
  for (double v : d.a) mx = std::max(mx, std::fabs(v));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < i; ++j)
      if (std::fabs(d(i, j) - d(j, i)) > 1e-12 * mx) throw Error("cmdscale: 'd' is not symmetric");
  // B = -1/2 J D^2 J (double centring of the squared dissimilarities).
  Matrix b(n, n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) b(i, j) = d(i, j) * d(i, j);
  Vec rm(n, 0.0), cm(n, 0.0);
  double gm = 0.0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      rm[i] += b(i, j);
      cm[j] += b(i, j);
    }
  for (size_t i = 0; i < n; ++i) {
    gm += rm[i];
    rm[i] /= static_cast<double>(n);
    cm[i] /= static_cast<double>(n);
  }
  gm /= static_cast<double>(n) * static_cast<double>(n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) b(i, j) = -0.5 * (((b(i, j) - rm[i]) - cm[j]) + gm);
  EigenSym es = eigen_sym(b);
  Cmdscale out;
  out.eig = es.values;
  size_t m = 0;
  while (m < static_cast<size_t>(k) && out.eig[m] > 0.0) ++m;
  if (m == 0) throw Error("cmdscale: only 0 of the first " + std::to_string(k) + " eigenvalues are > 0");
  if (m < static_cast<size_t>(k))
    out.warning = "only " + std::to_string(m) + " of the first " + std::to_string(k) + " eigenvalues are > 0";
  out.points = Matrix(n, m);
  for (size_t c = 0; c < m; ++c) {
    const double f = std::sqrt(out.eig[c]);
    for (size_t i = 0; i < n; ++i) out.points(i, c) = es.vectors(i, c) * f;
  }
  double top = 0.0, sabs = 0.0, spos = 0.0;
  for (size_t c = 0; c < m; ++c) top += out.eig[c];
  for (double v : out.eig) {
    sabs += std::fabs(v);
    spos += std::max(v, 0.0);
  }
  out.gof[0] = top / sabs;
  out.gof[1] = top / spos;
  return out;
}

Cmdscale cmdscale(const Vec& d, size_t n, int k) {
  if (d.size() != (n < 2 ? 0 : n * (n - 1) / 2)) throw Error("cmdscale: dist vector length does not match n");
  Matrix m(n, n, 0.0);
  size_t idx = 0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) m(i, j) = m(j, i) = d[idx++];
  return cmdscale(m, k);
}

Matrix dist_matrix(const Matrix& x, const std::string& method, double p) {
  const Vec dv = dist(x, method, p);
  const size_t n = x.rows;
  Matrix m(n, n, 0.0);
  size_t idx = 0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) m(i, j) = m(j, i) = dv[idx++];
  return m;
}

}  // namespace rcli
