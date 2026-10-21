// Multivariate analysis: symmetric eigendecomposition, singular value decomposition, principal components (like R's
// prcomp), covariance/correlation matrices, Cholesky inverse, Mahalanobis distances, classical multidimensional
// scaling (like R's cmdscale) and full distance matrices. Data matrices are row-major: rows are observations,
// columns are variables. None of these functions accept NaN or infinite input (they throw rcli::Error).
//
// Sign convention (eigenvectors, singular vectors, PCA rotations, MDS axes): every returned vector is scaled to unit
// length and its sign is chosen so that its largest-magnitude component is positive; components whose magnitudes tie
// within a relative 1e-9 go to the lowest index. R's signs come from LAPACK and are platform dependent ("the signs of
// the columns of the rotation matrix are arbitrary", ?prcomp), so compare with R only up to column sign.
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

// Symmetric eigendecomposition of an n x n matrix (only the lower triangle is read, as LAPACK-based
// eigen(symmetric = TRUE) does). Householder tridiagonalisation + implicit symmetric QR with Wilkinson shifts. values sorted decreasing (as R's eigen(symmetric = TRUE)); vectors are the columns
// of an n x n orthogonal matrix, column j belonging to values[j].
struct EigenSym {
  Vec values;
  Matrix vectors;
};
EigenSym eigen_sym(const Matrix& a);

// Thin singular value decomposition of an n x p matrix (any shape): x = u diag(d) v', k = min(n, p), d decreasing
// and non-negative, u n x k and v p x k with orthonormal columns (also for rank-deficient x, where the null-space
// columns are completed deterministically). Sign convention applies to the columns of v; u follows.
struct Svd {
  Vec d;
  Matrix u, v;
};
Svd svd(const Matrix& x);

// Principal components like R's prcomp(x, center, scale.): SVD of the centred (and scaled) data.
// sdev = d / sqrt(max(1, n - 1)) (length k = min(n, p)), rotation p x k (= v), x scores n x k (centred data times
// rotation). center holds the column means (empty when center = false); scale holds the column standard deviations
// (n - 1 divisor), or the root mean square sqrt(sum x^2 / (n - 1)) when center = false, as R's scale() does
// (empty when scale = false). prop_var and cum_prop are summary(prcomp)'s importance rows, unrounded.
// scale = true with a constant (or, uncentred, all-zero) column throws
// "prcomp: cannot rescale a constant/zero column to unit variance".
struct Pca {
  Vec sdev;
  Matrix rotation;
  Matrix x;
  Vec center, scale;
  Vec prop_var, cum_prop;
};
Pca prcomp(const Matrix& x, bool center = true, bool scale = false);

// Sample covariance and Pearson correlation matrices (n - 1 divisor), p x p. n < 2 gives an all-NaN matrix (R
// returns NA). In cor_matrix, rows and columns of a zero-variance variable are NaN, diagonal included (R gives NA
// with a warning); other diagonal entries are exactly 1 and off-diagonals are clamped to [-1, 1].
Matrix cov_matrix(const Matrix& x);
Matrix cor_matrix(const Matrix& x);

// Cholesky factor of a symmetric positive definite matrix: upper-triangular r with a = r'r (R's chol()). Throws
// "chol: the leading minor of order k is not positive" (k 1-based) otherwise.
Matrix cholesky(const Matrix& a);
// Inverse of a symmetric positive definite matrix via its Cholesky factor (R's chol2inv(chol(a))), exactly symmetric.
Matrix chol_inverse(const Matrix& a);
// Solution of a x = b for symmetric positive definite a.
Vec chol_solve(const Matrix& a, const Vec& b);

// Squared Mahalanobis distances of the rows of x from center, like R's mahalanobis(x, center, cov, inverted).
// cov must be symmetric positive definite (it is inverted through its Cholesky factor); with inverted = true, cov
// is taken to be the inverse covariance already and is used as given.
Vec mahalanobis(const Matrix& x, const Vec& center, const Matrix& cov, bool inverted = false);

// Classical (metric) multidimensional scaling, principal coordinates (Gower 1966; Torgerson 1958), like R's
// cmdscale(d, k, eig = TRUE). d is a full symmetric n x n dissimilarity matrix (its diagonal is used as given and
// should be 0). eig holds all n eigenvalues of the doubly centred matrix
// B = -1/2 J D^2 J (decreasing). points is n x m where m is the number of positive eigenvalues among the first k
// (R returns fewer columns with a warning in that case); column j = eigenvector j times sqrt(eig[j]).
// gof = (sum of the first m eigenvalues) / (sum |eig|, sum max(eig, 0)), as R's GOF.
// Errors as R: "'k' must be in {1, 2, ..  n - 1}", "only 0 of the first k eigenvalues are > 0",
// "NA values not allowed in 'd'"; also a non-square or non-symmetric d.
struct Cmdscale {
  Matrix points;
  Vec eig;
  double gof[2] = {kNaN, kNaN};
  std::string warning;  // R's warning when m < k: "only m of the first k eigenvalues are > 0"; empty otherwise
};
Cmdscale cmdscale(const Matrix& d, int k = 2);
// Same from a dist() vector (R's lower-triangle order, see rcli::dist) on n objects.
Cmdscale cmdscale(const Vec& d, size_t n, int k = 2);

// Full n x n distance matrix between rows (as.matrix(dist(x, method, p = p))); methods as rcli::dist:
// "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" (with p).
Matrix dist_matrix(const Matrix& x, const std::string& method = "euclidean", double p = 2.0);

}  // namespace rcli
