// Matrices and multivariate analysis: Mat (dense, row-major), eigen/svd/chol/solve, prcomp, cmdscale, mahalanobis.
#include <algorithm>
#include <cstdio>

#include "host.hpp"
#include "rcli/multivar.hpp"

namespace rcli::host {

namespace {

struct PcaModel {
  Pca fit;
  Strs columns;  // variable names (Frame columns, or V1..Vp for a Mat)
};

Strs default_names(size_t p) {
  Strs out;
  for (size_t j = 0; j < p; ++j) out.push_back("V" + std::to_string(j + 1));
  return out;
}

Matrix numeric_matrix(const Frame& f, const char* what) {
  for (size_t j = 0; j < f.cols.size(); ++j)
    if (f.cols[j].text)
      throw dslvm::Error(std::string(what) + ": column '" + f.names[j] + "' is text; select the numeric columns first (Frame.select)");
  if (f.cols.empty() || f.nrow() == 0) throw dslvm::Error(std::string(what) + ": the frame is empty");
  Matrix m(f.nrow(), f.cols.size());
  for (size_t j = 0; j < f.cols.size(); ++j)
    for (size_t i = 0; i < f.nrow(); ++i) {
      double v = f.cols[j].num[i];
      if (std::isnan(v)) throw dslvm::Error(std::string(what) + ": column '" + f.names[j] + "' has missing values (Frame.drop_nan)");
      m(i, j) = v;
    }
  return m;
}

size_t idx(int64_t i, size_t n, const char* what) { return static_cast<size_t>(as_index(i, n, what)); }

void need_square(const Matrix& a, const char* what) {
  if (a.rows != a.cols || a.rows == 0) throw dslvm::Error(std::string(what) + ": the matrix must be square and non-empty");
}

void need_finite(const Matrix& a, const char* what) {
  for (double v : a.a)
    if (!std::isfinite(v)) throw dslvm::Error(std::string(what) + ": the matrix has NaN or infinite entries");
}

Matrix transpose(const Matrix& a) {
  Matrix t(a.cols, a.rows);
  for (size_t i = 0; i < a.rows; ++i)
    for (size_t j = 0; j < a.cols; ++j) t(j, i) = a(i, j);
  return t;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
  if (a.cols != b.rows)
    throw dslvm::Error("Mat.mul: non-conformable (" + std::to_string(a.rows) + " x " + std::to_string(a.cols) + " times " +
                       std::to_string(b.rows) + " x " + std::to_string(b.cols) + ")");
  Matrix c(a.rows, b.cols);
  for (size_t i = 0; i < a.rows; ++i)
    for (size_t k = 0; k < a.cols; ++k) {
      const double aik = a(i, k);
      for (size_t j = 0; j < b.cols; ++j) c(i, j) += aik * b(k, j);
    }
  return c;
}

// LU factorisation with partial pivoting (Golub & Van Loan, Matrix Computations, Algorithm 3.4.1). Returns the
// determinant sign through `sign`; throws on an exactly singular or numerically singular pivot (R: "system is
// exactly singular" / "computationally singular").
struct Lu {
  Matrix lu;
  std::vector<size_t> perm;
  int sign = 1;
};

Lu lu_factor(const Matrix& a, const char* what) {
  need_square(a, what);
  need_finite(a, what);
  Lu f{a, {}, 1};
  const size_t n = a.rows;
  f.perm.resize(n);
  for (size_t i = 0; i < n; ++i) f.perm[i] = i;
  double scale = 0.0;
  for (double v : a.a) scale = std::max(scale, std::fabs(v));
  for (size_t k = 0; k < n; ++k) {
    size_t p = k;
    for (size_t i = k + 1; i < n; ++i)
      if (std::fabs(f.lu(i, k)) > std::fabs(f.lu(p, k))) p = i;
    if (f.lu(p, k) == 0.0) throw dslvm::Error(std::string(what) + ": the matrix is exactly singular");
    if (std::fabs(f.lu(p, k)) < scale * static_cast<double>(n) * kEps)
      throw dslvm::Error(std::string(what) + ": the matrix is computationally singular");
    if (p != k) {
      for (size_t j = 0; j < n; ++j) std::swap(f.lu(k, j), f.lu(p, j));
      std::swap(f.perm[k], f.perm[p]);
      f.sign = -f.sign;
    }
    for (size_t i = k + 1; i < n; ++i) {
      f.lu(i, k) /= f.lu(k, k);
      const double m = f.lu(i, k);
      for (size_t j = k + 1; j < n; ++j) f.lu(i, j) -= m * f.lu(k, j);
    }
  }
  return f;
}

Vec lu_solve(const Lu& f, const Vec& b) {
  const size_t n = f.lu.rows;
  Vec x(n);
  for (size_t i = 0; i < n; ++i) {
    double s = b[f.perm[i]];
    for (size_t j = 0; j < i; ++j) s -= f.lu(i, j) * x[j];
    x[i] = s;
  }
  for (size_t i = n; i-- > 0;) {
    double s = x[i];
    for (size_t j = i + 1; j < n; ++j) s -= f.lu(i, j) * x[j];
    x[i] = s / f.lu(i, i);
  }
  return x;
}

double determinant(const Matrix& a) {
  need_square(a, "Mat.det");
  need_finite(a, "Mat.det");
  Lu f{a, {}, 1};
  const size_t n = a.rows;
  double det = 1.0;
  for (size_t k = 0; k < n; ++k) {
    size_t p = k;
    for (size_t i = k + 1; i < n; ++i)
      if (std::fabs(f.lu(i, k)) > std::fabs(f.lu(p, k))) p = i;
    if (f.lu(p, k) == 0.0) return 0.0;
    if (p != k) {
      for (size_t j = 0; j < n; ++j) std::swap(f.lu(k, j), f.lu(p, j));
      det = -det;
    }
    det *= f.lu(k, k);
    for (size_t i = k + 1; i < n; ++i) {
      const double m = f.lu(i, k) / f.lu(k, k);
      for (size_t j = k + 1; j < n; ++j) f.lu(i, j) -= m * f.lu(k, j);
    }
  }
  return det;
}

Matrix inverse(const Matrix& a) {
  Lu f = lu_factor(a, "Mat.inverse");
  const size_t n = a.rows;
  Matrix inv(n, n);
  for (size_t j = 0; j < n; ++j) {
    Vec e(n, 0.0);
    e[j] = 1.0;
    Vec col = lu_solve(f, e);
    for (size_t i = 0; i < n; ++i) inv(i, j) = col[i];
  }
  return inv;
}

std::string mat_text(const Matrix& a, const std::string& spec) {
  std::string out;
  for (size_t i = 0; i < a.rows; ++i) {
    if (i) out += "\n";
    for (size_t j = 0; j < a.cols; ++j) {
      if (j) out += " ";
      out += std::isnan(a(i, j)) ? "NA" : format_number(spec, a(i, j));
    }
  }
  return out;
}

}  // namespace

void bind_multivar(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };
  auto mat = [&h](Matrix m) { return h.hold("Mat", std::move(m)); };

  // ---- construction ----
  b.fn("matrix", {"Vec", "int", "bool"}, "Mat", [mat](A& a) {
    const Vec& v = obj<Vec>(a[0]);
    const int64_t nr = a[1].i;
    if (nr < 1 || v.empty() || v.size() % static_cast<size_t>(nr) != 0)
      throw dslvm::Error("matrix: the length " + std::to_string(v.size()) + " is not a positive multiple of nrow " + std::to_string(nr));
    const size_t r = static_cast<size_t>(nr), c = v.size() / r;
    Matrix m(r, c);
    for (size_t k = 0; k < v.size(); ++k) {
      if (a[2].b) m(k / c, k % c) = v[k];  // byrow
      else m(k % r, k / r) = v[k];         // R's default: fill by column
    }
    return mat(std::move(m));
  });
  b.fn("as_matrix", {"Frame"}, "Mat", [mat](A& a) { return mat(numeric_matrix(obj<Frame>(a[0]), "as_matrix")); });
  b.fn("diag", {"Vec"}, "Mat", [mat](A& a) {
    const Vec& v = obj<Vec>(a[0]);
    if (v.empty()) throw dslvm::Error("diag: empty vector");
    Matrix m(v.size(), v.size());
    for (size_t i = 0; i < v.size(); ++i) m(i, i) = v[i];
    return mat(std::move(m));
  });
  b.fn("identity", {"int"}, "Mat", [mat](A& a) {
    if (a[0].i < 1 || a[0].i > 10000) throw dslvm::Error("identity: n must be between 1 and 10000");
    Matrix m(static_cast<size_t>(a[0].i), static_cast<size_t>(a[0].i));
    for (size_t i = 0; i < m.rows; ++i) m(i, i) = 1.0;
    return mat(std::move(m));
  });

  // ---- access ----
  b.m("Mat", "nrow", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<Matrix>(s).rows)); });
  b.m("Mat", "ncol", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<Matrix>(s).cols)); });
  b.m("Mat", "at", {"int", "int"}, "float", [](void* s, A& a) {
    auto& m = self<Matrix>(s);
    return Value::real(m(idx(a[0].i, m.rows, "Mat.at row"), idx(a[1].i, m.cols, "Mat.at column")));
  });
  b.m("Mat", "set", {"int", "int", "float"}, "Mat", [](void* s, A& a) {
    auto& m = self<Matrix>(s);
    m(idx(a[0].i, m.rows, "Mat.set row"), idx(a[1].i, m.cols, "Mat.set column")) = a[2].f;
    return Value::object("Mat", s);
  });
  b.m("Mat", "row", {"int"}, "Vec", [vec](void* s, A& a) {
    auto& m = self<Matrix>(s);
    const size_t i = idx(a[0].i, m.rows, "Mat.row");
    return vec(Vec(m.a.begin() + static_cast<std::ptrdiff_t>(i * m.cols), m.a.begin() + static_cast<std::ptrdiff_t>((i + 1) * m.cols)));
  });
  b.m("Mat", "col", {"int"}, "Vec", [vec](void* s, A& a) {
    auto& m = self<Matrix>(s);
    const size_t j = idx(a[0].i, m.cols, "Mat.col");
    Vec out(m.rows);
    for (size_t i = 0; i < m.rows; ++i) out[i] = m(i, j);
    return vec(std::move(out));
  });
  b.m("Mat", "values", {}, "Vec", [vec](void* s, A&) {  // column-major, as R's as.vector(m)
    auto& m = self<Matrix>(s);
    Vec out;
    out.reserve(m.a.size());
    for (size_t j = 0; j < m.cols; ++j)
      for (size_t i = 0; i < m.rows; ++i) out.push_back(m(i, j));
    return vec(std::move(out));
  });
  b.m("Mat", "diag", {}, "Vec", [vec](void* s, A&) {
    auto& m = self<Matrix>(s);
    Vec out;
    for (size_t i = 0; i < std::min(m.rows, m.cols); ++i) out.push_back(m(i, i));
    return vec(std::move(out));
  });
  b.m("Mat", "to_frame", {}, "Frame", [&h](void* s, A&) {
    auto& m = self<Matrix>(s);
    Frame f;
    Strs names = default_names(m.cols);
    for (size_t j = 0; j < m.cols; ++j) {
      Vec c(m.rows);
      for (size_t i = 0; i < m.rows; ++i) c[i] = m(i, j);
      Column col;
      col.num = std::move(c);
      f.set(names[j], std::move(col));
    }
    return h.hold("Frame", std::move(f));
  });
  b.m("Mat", "format", {"string"}, "string", [](void* s, A& a) { return Value::string(mat_text(self<Matrix>(s), a[0].s)); });
  b.m("Mat", "str", {}, "string", [](void* s, A&) { return Value::string(mat_text(self<Matrix>(s), "%g")); });

  // ---- algebra ----
  b.m("Mat", "t", {}, "Mat", [mat](void* s, A&) { return mat(transpose(self<Matrix>(s))); });
  b.m("Mat", "mul", {"Mat"}, "Mat", [mat](void* s, A& a) { return mat(multiply(self<Matrix>(s), obj<Matrix>(a[0]))); });
  b.m("Mat", "mul_vec", {"Vec"}, "Vec", [vec](void* s, A& a) {
    auto& m = self<Matrix>(s);
    const Vec& v = obj<Vec>(a[0]);
    if (v.size() != m.cols) throw dslvm::Error("Mat.mul_vec: the vector has length " + std::to_string(v.size()) + ", the matrix " + std::to_string(m.cols) + " columns");
    Vec out(m.rows, 0.0);
    for (size_t i = 0; i < m.rows; ++i)
      for (size_t j = 0; j < m.cols; ++j) out[i] += m(i, j) * v[j];
    return vec(std::move(out));
  });
  auto elementwise = [mat](const char* what, double sign) {
    return [mat, what, sign](void* s, A& a) {
      auto& x = self<Matrix>(s);
      auto& y = obj<Matrix>(a[0]);
      if (x.rows != y.rows || x.cols != y.cols) throw dslvm::Error(std::string(what) + ": the matrices have different shapes");
      Matrix c(x.rows, x.cols);
      for (size_t k = 0; k < x.a.size(); ++k) c.a[k] = x.a[k] + sign * y.a[k];
      return mat(std::move(c));
    };
  };
  b.m("Mat", "add", {"Mat"}, "Mat", elementwise("Mat.add", 1.0));
  b.m("Mat", "sub", {"Mat"}, "Mat", elementwise("Mat.sub", -1.0));
  b.m("Mat", "scale", {"float"}, "Mat", [mat](void* s, A& a) {
    Matrix c = self<Matrix>(s);
    for (double& v : c.a) v *= a[0].f;
    return mat(std::move(c));
  });
  b.m("Mat", "solve", {"Vec"}, "Vec", [vec](void* s, A& a) {
    const Vec& rhs = obj<Vec>(a[0]);
    auto& m = self<Matrix>(s);
    if (rhs.size() != m.rows) throw dslvm::Error("Mat.solve: the right-hand side has length " + std::to_string(rhs.size()) + ", the matrix " + std::to_string(m.rows) + " rows");
    return vec(lu_solve(lu_factor(m, "Mat.solve"), rhs));
  });
  b.m("Mat", "inverse", {}, "Mat", [mat](void* s, A&) { return mat(inverse(self<Matrix>(s))); });
  b.m("Mat", "det", {}, "float", [](void* s, A&) { return Value::real(determinant(self<Matrix>(s))); });
  b.m("Mat", "chol", {}, "Mat", [mat](void* s, A&) { return mat(cholesky(self<Matrix>(s))); });
  b.m("Mat", "cov", {}, "Mat", [mat](void* s, A&) { return mat(cov_matrix(self<Matrix>(s))); });
  b.m("Mat", "cor", {}, "Mat", [mat](void* s, A&) { return mat(cor_matrix(self<Matrix>(s))); });
  b.m("Mat", "dist", {"string"}, "Mat", [mat](void* s, A& a) { return mat(dist_matrix(self<Matrix>(s), a[0].s)); });

  // ---- decompositions ----
  b.m("Mat", "eigen", {}, "Eigen", [&h](void* s, A&) {
    auto& m = self<Matrix>(s);
    need_square(m, "Mat.eigen");
    return h.hold("Eigen", eigen_sym(m));
  });
  b.m("Eigen", "values", {}, "Vec", [vec](void* s, A&) { return vec(self<EigenSym>(s).values); });
  b.m("Eigen", "vectors", {}, "Mat", [mat](void* s, A&) { return mat(self<EigenSym>(s).vectors); });
  b.m("Mat", "svd", {}, "Svd", [&h](void* s, A&) { return h.hold("Svd", svd(self<Matrix>(s))); });
  b.m("Svd", "d", {}, "Vec", [vec](void* s, A&) { return vec(self<Svd>(s).d); });
  b.m("Svd", "u", {}, "Mat", [mat](void* s, A&) { return mat(self<Svd>(s).u); });
  b.m("Svd", "v", {}, "Mat", [mat](void* s, A&) { return mat(self<Svd>(s).v); });

  // ---- PCA, MDS, Mahalanobis ----
  auto pca = [&h](Matrix x, Strs names, bool center, bool scale) {
    PcaModel m{prcomp(x, center, scale), std::move(names)};
    return h.hold("Pca", std::move(m));
  };
  b.fn("prcomp", {"Frame", "bool"}, "Pca", [pca](A& a) {
    return pca(numeric_matrix(obj<Frame>(a[0]), "prcomp"), obj<Frame>(a[0]).names, true, a[1].b);
  });
  b.fn("prcomp", {"Frame", "bool", "bool"}, "Pca", [pca](A& a) {
    return pca(numeric_matrix(obj<Frame>(a[0]), "prcomp"), obj<Frame>(a[0]).names, a[1].b, a[2].b);
  });
  b.fn("prcomp", {"Mat", "bool"}, "Pca", [pca](A& a) {
    const Matrix& x = obj<Matrix>(a[0]);
    return pca(x, default_names(x.cols), true, a[1].b);
  });
  b.m("Pca", "sdev", {}, "Vec", [vec](void* s, A&) { return vec(self<PcaModel>(s).fit.sdev); });
  b.m("Pca", "rotation", {}, "Mat", [mat](void* s, A&) { return mat(self<PcaModel>(s).fit.rotation); });
  b.m("Pca", "scores", {}, "Mat", [mat](void* s, A&) { return mat(self<PcaModel>(s).fit.x); });
  b.m("Pca", "center", {}, "Vec", [vec](void* s, A&) { return vec(self<PcaModel>(s).fit.center); });
  b.m("Pca", "scale", {}, "Vec", [vec](void* s, A&) { return vec(self<PcaModel>(s).fit.scale); });
  b.m("Pca", "proportion", {}, "Vec", [vec](void* s, A&) { return vec(self<PcaModel>(s).fit.prop_var); });
  b.m("Pca", "cumulative", {}, "Vec", [vec](void* s, A&) { return vec(self<PcaModel>(s).fit.cum_prop); });
  b.m("Pca", "loading", {"string", "int"}, "float", [](void* s, A& a) {
    auto& m = self<PcaModel>(s);
    for (size_t j = 0; j < m.columns.size(); ++j)
      if (m.columns[j] == a[0].s) return Value::real(m.fit.rotation(j, idx(a[1].i, m.fit.rotation.cols, "Pca.loading component")));
    throw dslvm::Error("Pca.loading: no variable '" + a[0].s + "'");
  });
  b.m("Pca", "summary", {}, "string", [](void* s, A&) {
    auto& f = self<PcaModel>(s).fit;
    std::string out = "Importance of components:\n                         ";
    char buf[64];
    for (size_t k = 0; k < f.sdev.size(); ++k) { std::snprintf(buf, sizeof buf, "%9s", ("PC" + std::to_string(k + 1)).c_str()); out += buf; }
    const char* labels[] = {"Standard deviation      ", "Proportion of Variance  ", "Cumulative Proportion   "};
    const Vec* rows[] = {&f.sdev, &f.prop_var, &f.cum_prop};
    for (int r = 0; r < 3; ++r) {
      out += "\n";
      out += labels[r];
      for (double v : *rows[r]) { std::snprintf(buf, sizeof buf, " %8.5f", v); out += buf; }
    }
    return Value::string(out);
  });

  b.fn("cmdscale", {"Mat", "int"}, "Mds", [&h](A& a) {
    if (a[1].i < 1) throw dslvm::Error("cmdscale: 'k' must be in {1, 2, ..  n - 1}");
    return h.hold("Mds", cmdscale(obj<Matrix>(a[0]), static_cast<int>(a[1].i)));
  });
  b.m("Mds", "points", {}, "Mat", [mat](void* s, A&) { return mat(self<Cmdscale>(s).points); });
  b.m("Mds", "eig", {}, "Vec", [vec](void* s, A&) { return vec(self<Cmdscale>(s).eig); });
  b.m("Mds", "gof", {}, "Vec", [vec](void* s, A&) { return vec(Vec{self<Cmdscale>(s).gof[0], self<Cmdscale>(s).gof[1]}); });
  b.m("Mds", "warning", {}, "string", [](void* s, A&) { return Value::string(self<Cmdscale>(s).warning); });

  b.fn("mahalanobis", {"Mat", "Vec", "Mat"}, "Vec", [vec](A& a) {
    return vec(mahalanobis(obj<Matrix>(a[0]), obj<Vec>(a[1]), obj<Matrix>(a[2])));
  });
}

}  // namespace rcli::host
