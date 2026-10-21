// Linear models (lm), generalised linear models (glm) and nonlinear least squares (nls) with R semantics.
//
// Clean-room implementation from published sources only:
//  - Householder QR: G. H. Golub and C. F. Van Loan, "Matrix Computations", 4th ed., sections 5.1-5.2 and 5.4.1;
//    the limited column pivoting rule (columns whose remaining 2-norm is negligible relative to their original norm,
//    tol = 1e-7, are moved to the right-hand edge) follows the description in R's ?qr.
//  - Least squares via QR, (X'X)^-1 = R^-1 R^-T: Golub and Van Loan section 5.3; Dongarra et al., LINPACK Users' Guide.
//  - Summary statistics (R^2, F, logLik): R's ?summary.lm and ?logLik documentation; standard theory (Seber & Lee,
//    "Linear Regression Analysis").
//  - glm: iteratively reweighted least squares, P. McCullagh and J. A. Nelder, "Generalized Linear Models", 2nd ed.,
//    chapter 2.5; starting values, convergence rule and family definitions as documented in ?glm, ?glm.control and
//    ?family; deviance residuals and AIC from the family log densities (McCullagh & Nelder ch. 2, 6, 8).
//  - nls: Gauss-Newton with step halving and the relative-offset convergence criterion of D. M. Bates and
//    D. G. Watts, "Nonlinear Regression Analysis and Its Applications" (Wiley 1988), sections 2.2 and 2.2.3, and
//    Bates & Watts (1981) "A relative offset orthogonality convergence criterion for nonlinear least squares",
//    Technometrics 23, 179-183. Forward-difference Jacobian with step sqrt(eps)*|theta| (sqrt(eps) at 0), established
//    by black-box comparison with R's numericDeriv. "port": a projected Levenberg-Marquardt method (Nocedal & Wright,
//    "Numerical Optimization", 2nd ed., section 10.3 and 16.7) on the box constraints.
#include <algorithm>
#include <cctype>
#include <iterator>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <set>

#include "rcli/dist.hpp"
#include "rcli/models.hpp"

namespace rcli {

// ===================================================================================================================
// Data
// ===================================================================================================================

size_t Data::nrow() const {
  if (!cols.empty()) return cols[0].size();
  if (!text_cols.empty()) return text_cols[0].size();
  return 0;
}

bool Data::is_text(const std::string& name) const {
  return std::find(text_names.begin(), text_names.end(), name) != text_names.end();
}

bool Data::has(const std::string& name) const {
  return is_text(name) || std::find(names.begin(), names.end(), name) != names.end();
}

const Vec& Data::col(const std::string& name) const {
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name) return cols[i];
  if (is_text(name)) throw Error("Data: column '" + name + "' is text, not numeric");
  throw Error("Data: object '" + name + "' not found");
}

const std::vector<std::string>& Data::text(const std::string& name) const {
  for (size_t i = 0; i < text_names.size(); ++i)
    if (text_names[i] == name) return text_cols[i];
  throw Error("Data: text column '" + name + "' not found");
}

void Data::add(const std::string& name, Vec v) {
  if ((!cols.empty() || !text_cols.empty()) && v.size() != nrow())
    throw Error("Data: column '" + name + "' has length " + std::to_string(v.size()) + ", expected " +
                std::to_string(nrow()));
  if (is_text(name)) throw Error("Data: column '" + name + "' already exists as a text column");
  if (std::find(order.begin(), order.end(), name) == order.end()) order.push_back(name);
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name) {
      cols[i] = std::move(v);
      return;
    }
  names.push_back(name);
  cols.push_back(std::move(v));
}

void Data::add_text(const std::string& name, std::vector<std::string> v) {
  if ((!cols.empty() || !text_cols.empty()) && v.size() != nrow())
    throw Error("Data: column '" + name + "' has length " + std::to_string(v.size()) + ", expected " +
                std::to_string(nrow()));
  if (std::find(names.begin(), names.end(), name) != names.end())
    throw Error("Data: column '" + name + "' already exists as a numeric column");
  if (std::find(order.begin(), order.end(), name) == order.end()) order.push_back(name);
  for (size_t i = 0; i < text_names.size(); ++i)
    if (text_names[i] == name) {
      text_cols[i] = std::move(v);
      return;
    }
  text_names.push_back(name);
  text_cols.push_back(std::move(v));
}

std::vector<std::string> Data::column_order() const {
  std::vector<std::string> out;
  for (const auto& n : order)
    if (has(n) && std::find(out.begin(), out.end(), n) == out.end()) out.push_back(n);
  for (const auto& n : names)
    if (std::find(out.begin(), out.end(), n) == out.end()) out.push_back(n);
  for (const auto& n : text_names)
    if (std::find(out.begin(), out.end(), n) == out.end()) out.push_back(n);
  return out;
}

int CoefTable::index(const std::string& name) const {
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name) return static_cast<int>(i);
  return -1;
}

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

// ===================================================================================================================
// Formula parsing
// ===================================================================================================================

struct Tok {
  enum Kind { Name, Num, Op, End } kind;
  std::string text;
};

bool name_start(char c) { return std::isalpha(static_cast<unsigned char>(c)) || c == '.' || c == '_'; }
bool name_char(char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '.' || c == '_'; }

std::vector<Tok> tokenize(const std::string& s) {
  std::vector<Tok> out;
  size_t i = 0;
  while (i < s.size()) {
    char c = s[i];
    if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
    if (std::isdigit(static_cast<unsigned char>(c)) ||
        (c == '.' && i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i + 1])))) {
      size_t j = i;
      while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
      if (j < s.size() && (s[j] == 'e' || s[j] == 'E')) {
        size_t k = j + 1;
        if (k < s.size() && (s[k] == '+' || s[k] == '-')) ++k;
        if (k < s.size() && std::isdigit(static_cast<unsigned char>(s[k]))) {
          j = k;
          while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
        }
      }
      if (j < s.size() && s[j] == 'L') ++j;  // R integer literal suffix
      out.push_back({Tok::Num, s.substr(i, j - i)});
      i = j;
    } else if (name_start(c)) {
      size_t j = i;
      while (j < s.size() && name_char(s[j])) ++j;
      out.push_back({Tok::Name, s.substr(i, j - i)});
      i = j;
    } else if (c == '`') {
      size_t j = s.find('`', i + 1);
      if (j == std::string::npos) throw Error("formula: unterminated backquote in '" + s + "'");
      out.push_back({Tok::Name, s.substr(i + 1, j - i - 1)});
      i = j + 1;
    } else if (std::string("~+-*:^()/,").find(c) != std::string::npos) {
      out.push_back({Tok::Op, std::string(1, c)});
      ++i;
    } else {
      throw Error(std::string("formula: unexpected character '") + c + "' in '" + s + "'");
    }
  }
  out.push_back({Tok::End, ""});
  return out;
}

// Deparse a numeric literal the way R prints it (15 significant digits).
std::string num_text(const std::string& lit) {
  std::string t = lit;
  if (!t.empty() && t.back() == 'L') t.pop_back();
  char* end = nullptr;
  double v = std::strtod(t.c_str(), &end);
  if (end == t.c_str() || *end != '\0') throw Error("formula: bad number '" + lit + "'");
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.15g", v);
  return buf;
}

// A set of model terms (each a sorted list of variable indices) plus an intercept instruction.
struct TermSet {
  std::vector<std::vector<int>> terms;
  int icpt = 0;  // 0: unspecified, +1: add intercept, -1: remove intercept
};

bool has_term(const TermSet& s, const std::vector<int>& t) {
  return std::find(s.terms.begin(), s.terms.end(), t) != s.terms.end();
}

TermSet ts_plus(const TermSet& a, const TermSet& b) {
  TermSet r = a;
  for (const auto& t : b.terms)
    if (!has_term(r, t)) r.terms.push_back(t);
  if (b.icpt != 0) r.icpt = b.icpt;
  return r;
}

TermSet ts_minus(const TermSet& a, const TermSet& b) {
  TermSet r;
  for (const auto& t : a.terms)
    if (!has_term(b, t)) r.terms.push_back(t);
  r.icpt = a.icpt;
  if (b.icpt == 1) r.icpt = -1;
  if (b.icpt == -1) r.icpt = 1;
  return r;
}

TermSet ts_colon(const TermSet& a, const TermSet& b) {
  TermSet r;
  for (const auto& x : a.terms)
    for (const auto& y : b.terms) {
      std::vector<int> u;
      std::set_union(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(u));
      if (!has_term(r, u)) r.terms.push_back(u);
    }
  return r;
}

TermSet ts_times(const TermSet& a, const TermSet& b) {
  TermSet ab = a, bb = b;
  ab.icpt = bb.icpt = 0;
  return ts_plus(ts_plus(ab, bb), ts_colon(ab, bb));
}

class FormulaParser {
 public:
  FormulaParser(const std::string& text, const Data* data) : text_(text), toks_(tokenize(text)), data_(data) {}

  Formula run() {
    Formula f;
    bool has_tilde = false;
    for (const auto& t : toks_)
      if (t.kind == Tok::Op && t.text == "~") has_tilde = true;
    if (!has_tilde) throw Error("formula: missing '~' in '" + text_ + "'");
    if (!is_op("~")) {
      f.response = primary_label();
      var_index(f.response);
      response_vars_ = base_columns(f.response);
    }
    expect("~");
    TermSet rhs = sum();
    if (peek().kind != Tok::End) throw Error("formula: unexpected '" + peek().text + "' in '" + text_ + "'");
    f.intercept = rhs.icpt != -1;
    for (const auto& v : vars_)
      if (v != f.response) f.variables.push_back(v);
    std::stable_sort(rhs.terms.begin(), rhs.terms.end(),
                     [](const std::vector<int>& x, const std::vector<int>& y) { return x.size() < y.size(); });
    for (const auto& t : rhs.terms) {
      std::string lab;
      for (size_t k = 0; k < t.size(); ++k) lab += (k ? ":" : "") + vars_[t[k]];
      f.term_labels.push_back(lab);
    }
    return f;
  }

 private:
  const Tok& peek() const { return toks_[pos_]; }
  bool is_op(const char* op) const { return peek().kind == Tok::Op && peek().text == op; }
  void expect(const char* op) {
    if (!is_op(op)) throw Error(std::string("formula: expected '") + op + "' in '" + text_ + "'");
    ++pos_;
  }

  int var_index(const std::string& label) {
    for (size_t i = 0; i < vars_.size(); ++i)
      if (vars_[i] == label) return static_cast<int>(i);
    vars_.push_back(label);
    return static_cast<int>(vars_.size() - 1);
  }

  // Plain column names referenced by a variable expression (used to exclude the response from ".").
  static std::vector<std::string> base_columns(const std::string& label) {
    std::vector<std::string> out;
    for (const auto& t : tokenize(label))
      if (t.kind == Tok::Name) out.push_back(t.text);
    return out;
  }

  // name | name(args): returns the R deparsed label of a variable expression.
  std::string primary_label() {
    if (peek().kind != Tok::Name) throw Error("formula: expected a variable in '" + text_ + "'");
    std::string name = peek().text;
    ++pos_;
    if (!is_op("(")) return name;
    ++pos_;
    std::string inner;
    if (name == "I") {
      inner = i_body();
    } else {
      inner = primary_label();
    }
    expect(")");
    return name + "(" + inner + ")";
  }

  // Body of I(...): v, v^k, v^-k (deparsed without spaces around ^ as R does).
  std::string i_body() {
    std::string base = primary_label();
    if (!is_op("^")) return base;
    ++pos_;
    std::string sign;
    if (is_op("-")) { sign = "-"; ++pos_; }
    if (peek().kind != Tok::Num) throw Error("formula: I() supports only v^k for a numeric literal k");
    std::string k = num_text(peek().text);
    ++pos_;
    return base + "^" + sign + k;
  }

  TermSet sum() {
    TermSet r;
    bool first = true;
    while (true) {
      bool minus = false;
      if (is_op("+") || is_op("-")) {
        minus = is_op("-");
        ++pos_;
      } else if (!first) {
        break;
      }
      TermSet t = prod();
      r = minus ? ts_minus(r, t) : ts_plus(r, t);
      first = false;
    }
    return r;
  }

  TermSet prod() {
    TermSet r = inter();
    while (is_op("*") || is_op("/")) {
      if (is_op("/")) throw Error("formula: nesting operator '/' is not supported");
      ++pos_;
      r = ts_times(r, inter());
    }
    return r;
  }

  TermSet inter() {
    TermSet r = power();
    while (is_op(":")) {
      ++pos_;
      r = ts_colon(r, power());
    }
    return r;
  }

  TermSet power() {
    TermSet r = atom();
    if (is_op("^")) {
      ++pos_;
      if (peek().kind != Tok::Num) throw Error("formula: '^' needs a positive integer power");
      double k = std::strtod(peek().text.c_str(), nullptr);
      ++pos_;
      if (!(k >= 1) || k != std::floor(k)) throw Error("formula: invalid power in formula");
      TermSet base = r;
      for (int i = 1; i < static_cast<int>(k); ++i) r = ts_times(r, base);
    }
    return r;
  }

  TermSet atom() {
    TermSet r;
    const Tok& t = peek();
    if (t.kind == Tok::Op && t.text == "(") {
      ++pos_;
      r = sum();
      expect(")");
      return r;
    }
    if (t.kind == Tok::Num) {
      double v = std::strtod(t.text.c_str(), nullptr);
      ++pos_;
      if (v == 0) r.icpt = -1;
      else if (v == 1) r.icpt = 1;
      else throw Error("formula: invalid model formula term '" + t.text + "'");
      return r;
    }
    if (t.kind == Tok::Name && t.text == "." && !(toks_[pos_ + 1].kind == Tok::Op && toks_[pos_ + 1].text == "(")) {
      ++pos_;
      if (!data_) throw Error("formula: '.' in formula and no 'data' argument");
      for (const auto& nm : data_->column_order()) {
        if (std::find(response_vars_.begin(), response_vars_.end(), nm) != response_vars_.end()) continue;
        std::vector<int> term{var_index(nm)};
        if (!has_term(r, term)) r.terms.push_back(term);
      }
      return r;
    }
    if (t.kind == Tok::Name) {
      std::string lab = primary_label();
      r.terms.push_back({var_index(lab)});
      return r;
    }
    throw Error("formula: unexpected '" + t.text + "' in '" + text_ + "'");
  }

  std::string text_;
  std::vector<Tok> toks_;
  size_t pos_ = 0;
  const Data* data_;
  std::vector<std::string> vars_;
  std::vector<std::string> response_vars_;
};

// Evaluate a variable expression label ("x", "log(x)", "I(x^2)", ...) against a data table.
Vec eval_var(const std::string& label, const Data& d) {
  for (size_t i = 0; i < d.names.size(); ++i)
    if (d.names[i] == label) return d.cols[i];
  size_t open = label.find('(');
  if (open == std::string::npos || label.back() != ')') return d.col(label);  // throws "not found"
  std::string fn = label.substr(0, open);
  std::string inner = label.substr(open + 1, label.size() - open - 2);
  if (fn == "I") {
    // inner is v or v^k (k may be negative); find the top-level '^'
    int depth = 0;
    size_t caret = std::string::npos;
    for (size_t i = 0; i < inner.size(); ++i) {
      if (inner[i] == '(') ++depth;
      else if (inner[i] == ')') --depth;
      else if (inner[i] == '^' && depth == 0) caret = i;
    }
    if (caret == std::string::npos) return eval_var(inner, d);
    Vec v = eval_var(inner.substr(0, caret), d);
    std::string ks = inner.substr(caret + 1);
    char* end = nullptr;
    double k = std::strtod(ks.c_str(), &end);
    if (end == ks.c_str() || *end != '\0') throw Error("formula: cannot evaluate '" + label + "'");
    for (double& x : v) x = std::pow(x, k);
    return v;
  }
  Vec v = eval_var(inner, d);
  if (fn == "log") for (double& x : v) x = std::log(x);
  else if (fn == "exp") for (double& x : v) x = std::exp(x);
  else if (fn == "sqrt") for (double& x : v) x = std::sqrt(x);
  else throw Error("formula: unsupported function '" + fn + "' in '" + label + "'");
  return v;
}

// Split a term label "a:log(b)" into its variable labels at top-level ':'.
std::vector<std::string> split_term(const std::string& lab) {
  std::vector<std::string> out;
  int depth = 0;
  size_t start = 0;
  for (size_t i = 0; i < lab.size(); ++i) {
    if (lab[i] == '(') ++depth;
    else if (lab[i] == ')') --depth;
    else if (lab[i] == ':' && depth == 0) {
      out.push_back(lab.substr(start, i - start));
      start = i + 1;
    }
  }
  out.push_back(lab.substr(start));
  return out;
}

// ===================================================================================================================
// Householder QR with limited column pivoting (Golub & Van Loan 5.1.2, 5.2.1; pivot rule as described in ?qr)
// ===================================================================================================================

struct QR {
  size_t n = 0, p = 0, rank = 0;
  std::vector<Vec> c;          // columns (pivoted order); upper triangle holds R after factorisation
  std::vector<Vec> v;          // Householder vectors (full length n, zeros above the diagonal)
  Vec beta;                    // H_k = I - beta_k v_k v_k'
  std::vector<size_t> pivot;   // pivot[j] = original index of the column at position j

  QR(const Matrix& X, double tol) : n(X.rows), p(X.cols), c(X.cols, Vec(X.rows)), pivot(X.cols) {
    for (size_t j = 0; j < p; ++j) {
      pivot[j] = j;
      for (size_t i = 0; i < n; ++i) c[j][i] = X(i, j);
    }
    Vec onorm(p);
    for (size_t j = 0; j < p; ++j) onorm[j] = norm2(c[j], 0);
    size_t lup = p;  // columns [lup, p) have been moved to the end as negligible
    size_t l = 0;
    for (; l < std::min(n, lup); ++l) {
      // Move columns with negligible remaining norm to the right-hand edge.
      while (l < lup) {
        double rem = norm2(c[l], l);
        if (onorm[l] > 0 && rem >= tol * onorm[l]) break;
        std::rotate(c.begin() + l, c.begin() + l + 1, c.end());
        std::rotate(pivot.begin() + l, pivot.begin() + l + 1, pivot.end());
        std::rotate(onorm.begin() + l, onorm.begin() + l + 1, onorm.end());
        --lup;
      }
      if (l >= lup) break;
      // Householder vector for c[l][l..n): x -> alpha e1, alpha = -sign(x0) ||x||.
      Vec hv(n, 0.0);
      double nx = norm2(c[l], l);
      double alpha = c[l][l] >= 0 ? -nx : nx;
      for (size_t i = l; i < n; ++i) hv[i] = c[l][i];
      hv[l] -= alpha;
      double vtv = 0;
      for (size_t i = l; i < n; ++i) vtv += hv[i] * hv[i];
      double b = vtv > 0 ? 2.0 / vtv : 0.0;
      for (size_t j = l + 1; j < p; ++j) apply(hv, b, l, c[j]);
      c[l][l] = alpha;
      for (size_t i = l + 1; i < n; ++i) c[l][i] = 0.0;
      v.push_back(std::move(hv));
      beta.push_back(b);
    }
    rank = v.size();
  }

  static double norm2(const Vec& x, size_t from) {
    // scaled to avoid overflow/underflow
    double scale = 0, ssq = 1;
    for (size_t i = from; i < x.size(); ++i) {
      if (x[i] == 0) continue;
      double a = std::fabs(x[i]);
      if (scale < a) {
        ssq = 1 + ssq * (scale / a) * (scale / a);
        scale = a;
      } else {
        ssq += (a / scale) * (a / scale);
      }
    }
    return scale * std::sqrt(ssq);
  }

  static void apply(const Vec& hv, double b, size_t l, Vec& y) {
    double s = 0;
    for (size_t i = l; i < hv.size(); ++i) s += hv[i] * y[i];
    s *= b;
    for (size_t i = l; i < hv.size(); ++i) y[i] -= s * hv[i];
  }

  Vec qty(Vec y) const {
    for (size_t k = 0; k < rank; ++k) apply(v[k], beta[k], k, y);
    return y;
  }
  Vec qy(Vec y) const {
    for (size_t k = rank; k-- > 0;) apply(v[k], beta[k], k, y);
    return y;
  }
  double r(size_t i, size_t j) const { return c[j][i]; }

  // Coefficients in pivoted order for the first `rank` columns.
  Vec solve_r(const Vec& qtyv) const {
    Vec b(rank);
    for (size_t i = rank; i-- > 0;) {
      double s = qtyv[i];
      for (size_t j = i + 1; j < rank; ++j) s -= r(i, j) * b[j];
      b[i] = s / r(i, i);
    }
    return b;
  }

  // (R'R)^-1 over the first rank columns (pivoted order), Golub & Van Loan 5.3.
  Matrix rinv_cov() const {
    size_t k = rank;
    Matrix Ri(k, k, 0.0);  // R^-1 (upper)
    for (size_t j = 0; j < k; ++j) {
      Ri(j, j) = 1.0 / r(j, j);
      for (size_t i = j; i-- > 0;) {
        double s = 0;
        for (size_t m = i + 1; m <= j; ++m) s += r(i, m) * Ri(m, j);
        Ri(i, j) = -s / r(i, i);
      }
    }
    Matrix C(k, k, 0.0);
    for (size_t i = 0; i < k; ++i)
      for (size_t j = i; j < k; ++j) {
        double s = 0;
        for (size_t m = j; m < k; ++m) s += Ri(i, m) * Ri(j, m);
        C(i, j) = C(j, i) = s;
      }
    return C;
  }
};

// Unpivot (R'R)^-1 into a p x p matrix in original column order; aliased rows/columns are 0.
Matrix unpivot_cov(const QR& qr) {
  Matrix C = qr.rinv_cov();
  Matrix out(qr.p, qr.p, 0.0);
  for (size_t i = 0; i < qr.rank; ++i)
    for (size_t j = 0; j < qr.rank; ++j) out(qr.pivot[i], qr.pivot[j]) = C(i, j);
  return out;
}

double pt_two_sided(double t, double df) {
  if (std::isnan(t) || std::isnan(df) || !(df > 0)) return kNaN;
  return 2.0 * Dist::t(df).p(std::fabs(t), false);
}
double pnorm_two_sided(double z) {
  if (std::isnan(z)) return kNaN;
  return 2.0 * Dist::norm(0, 1).p(-std::fabs(z));
}

// Fill a coefficient table from estimates (original order), aliased flags, cov_unscaled and a scale.
CoefTable make_table(const std::vector<std::string>& names, const Vec& est, const std::vector<bool>& aliased,
                     const Matrix& cov, double scale, double df, bool use_t) {
  CoefTable t;
  t.names = names;
  t.estimate = est;
  t.aliased = aliased;
  size_t p = names.size();
  t.std_error.assign(p, kNaN);
  t.statistic.assign(p, kNaN);
  t.p_value.assign(p, kNaN);
  for (size_t j = 0; j < p; ++j) {
    if (aliased[j]) continue;
    t.std_error[j] = std::sqrt(cov(j, j) * scale);
    t.statistic[j] = est[j] / t.std_error[j];
    t.p_value[j] = use_t ? pt_two_sided(t.statistic[j], df) : pnorm_two_sided(t.statistic[j]);
  }
  return t;
}

// Design rows for prediction (no response needed); NaN rows give NaN.
Vec linear_predict(const Formula& f, const CoefTable& coef, const Data& newdata) {
  std::vector<std::string> cn;
  Matrix X = f.design(newdata, &cn);
  if (cn.size() != coef.names.size()) throw Error("predict: design does not match the fitted model");
  Vec out(X.rows, 0.0);
  for (size_t i = 0; i < X.rows; ++i)
    for (size_t j = 0; j < X.cols; ++j)
      if (!coef.aliased[j]) out[i] += X(i, j) * coef.estimate[j];
  return out;
}

void check_finite_design(const char* fn, const Matrix& X, const Vec& y) {
  for (double v : X.a)
    if (std::isinf(v)) throw Error(std::string(fn) + ": NA/NaN/Inf in 'x'");
  for (double v : y)
    if (std::isinf(v)) throw Error(std::string(fn) + ": NA/NaN/Inf in 'y'");
}

Vec take(const Vec& v, const std::vector<size_t>& rows) {
  Vec out(rows.size());
  for (size_t i = 0; i < rows.size(); ++i) out[i] = v[rows[i]];
  return out;
}

// ===================================================================================================================
// Factors and the model matrix
//
// Coding follows J. M. Chambers and T. J. Hastie (eds.), "Statistical Models in S" (Wadsworth 1992), section 2.4.1:
// the columns of a term are the products of one column set per variable (the first variable varying fastest); a
// numeric variable gives itself; a factor F in term T is coded by contrasts (here treatment contrasts: indicators of
// levels 2..L, the first level being the baseline) when the margin T\F is already spanned by the terms before T, and
// by indicators of all L levels otherwise. "Spanned" was established by black-box experiments with R's model.matrix:
//  - a non-empty margin T\F is spanned when it is contained in an earlier term (y ~ z:x + x:a codes a by contrasts;
//    y ~ a:c + a:b codes c in a:c by indicators but b in a:b by contrasts);
//  - the empty margin (a factor main effect) is spanned by the intercept or, without an intercept, by an earlier
//    factor main effect (y ~ 0 + a + b: indicators for a, contrasts for b; y ~ 0 + x + a: indicators for a).
// Terms are already ordered by degree (R's terms()), so main effects precede interactions.
// ===================================================================================================================

// A double printed as R's as.character: up to 15 significant digits (the fewest that reproduce the value rounded to
// 15 digits), in fixed notation unless scientific notation is narrower (R's ?options "scipen" = 0 rule).
std::string r_num_label(double x) {
  if (x == 0) return "0";
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.14e", x);
  double target = std::strtod(buf, nullptr);
  int sig = 15;
  for (int k = 1; k <= 15; ++k) {
    std::snprintf(buf, sizeof buf, "%.*e", k - 1, x);
    if (std::strtod(buf, nullptr) == target) {
      sig = k;
      break;
    }
  }
  std::snprintf(buf, sizeof buf, "%.*e", sig - 1, x);
  std::string sci = buf;
  int ex = std::atoi(sci.c_str() + sci.find('e') + 1);
  int neg = x < 0 ? 1 : 0;
  int decimals = std::max(0, sig - 1 - ex);
  int fixed_width = ex >= 0 ? neg + ex + 1 + (decimals > 0 ? 1 + decimals : 0) : neg + 1 + (-ex) + sig;
  if (fixed_width > static_cast<int>(sci.size())) return sci;
  std::snprintf(buf, sizeof buf, "%.*f", decimals, x);
  return buf;
}

// factor(v) / as.factor(v): returns the inner expression.
bool factor_call(const std::string& label, std::string* inner) {
  for (const char* fn : {"factor(", "as.factor("}) {
    std::string pre = fn;
    if (label.size() > pre.size() && label.compare(0, pre.size(), pre) == 0 && label.back() == ')') {
      *inner = label.substr(pre.size(), label.size() - pre.size() - 1);
      return true;
    }
  }
  return false;
}

// One evaluated formula variable: numeric values, or a factor's per-row labels.
struct VarData {
  bool factor = false;
  bool numeric_levels = false;     // factor(v) of a numeric v: levels sort by value
  Vec num;                         // numeric variable values, or v for factor(v)
  std::vector<std::string> lab;    // factor labels per row
  std::vector<bool> miss;          // factor value missing (NaN in v)
};

VarData eval_variable(const std::string& label, const Data& d) {
  VarData v;
  size_t n = d.nrow();
  std::string inner;
  if (d.is_text(label)) {
    v.factor = true;
    v.lab = d.text(label);
    v.miss.assign(n, false);
  } else if (factor_call(label, &inner)) {
    v.factor = true;
    if (d.is_text(inner)) {
      v.lab = d.text(inner);
      v.miss.assign(n, false);
    } else {
      v.numeric_levels = true;
      v.num = eval_var(inner, d);
      v.lab.resize(v.num.size());
      v.miss.resize(v.num.size());
      for (size_t i = 0; i < v.num.size(); ++i) {
        v.miss[i] = std::isnan(v.num[i]);
        if (!v.miss[i]) v.lab[i] = r_num_label(v.num[i]);
      }
    }
  } else {
    v.num = eval_var(label, d);
  }
  if ((v.factor ? v.lab.size() : v.num.size()) != n)
    throw Error("formula: variable '" + label + "' has the wrong length");
  return v;
}

// Levels of a factor over the rows in mask: text in byte order (std::string comparison is unsigned-byte
// lexicographic, R's order in the C locale), numeric values in numeric order labelled by r_num_label.
std::vector<std::string> learn_levels(const VarData& v, const std::vector<bool>& mask) {
  std::vector<std::string> out;
  if (v.numeric_levels) {
    Vec vals;
    for (size_t i = 0; i < v.num.size(); ++i)
      if (mask[i] && !v.miss[i]) vals.push_back(v.num[i]);
    std::sort(vals.begin(), vals.end());
    std::set<std::string> seen;
    for (double x : vals) {
      std::string s = r_num_label(x);
      if (seen.insert(s).second) out.push_back(s);
    }
  } else {
    std::set<std::string> s;
    for (size_t i = 0; i < v.lab.size(); ++i)
      if (mask[i] && !v.miss[i]) s.insert(v.lab[i]);
    out.assign(s.begin(), s.end());
  }
  return out;
}

// The variables and terms of a formula: terms as lists of indices into the variable labels.
struct TermStructure {
  std::vector<std::string> vars;
  std::vector<std::vector<int>> terms;
};

TermStructure term_structure(const Formula& f) {
  TermStructure ts;
  for (const auto& lab : f.term_labels) {
    std::vector<int> t;
    for (const auto& var : split_term(lab)) {
      auto it = std::find(ts.vars.begin(), ts.vars.end(), var);
      if (it == ts.vars.end()) {
        ts.vars.push_back(var);
        t.push_back(static_cast<int>(ts.vars.size() - 1));
      } else {
        t.push_back(static_cast<int>(it - ts.vars.begin()));
      }
    }
    ts.terms.push_back(t);
  }
  return ts;
}

// Model matrix rows `rows` from evaluated variables and factor levels (levels[k] used when vars[k] is a factor).
Matrix build_design(const Formula& f, const TermStructure& ts, const std::vector<VarData>& vd,
                    const std::vector<std::vector<std::string>>& levels, const std::vector<size_t>& rows,
                    std::vector<std::string>* colnames, std::vector<int>* assign) {
  size_t n = rows.size();
  // Level codes per factor variable (-1: missing).
  std::vector<std::vector<int>> code(ts.vars.size());
  for (size_t k = 0; k < ts.vars.size(); ++k) {
    if (!vd[k].factor) continue;
    if (levels[k].size() < 2) throw Error("contrasts can be applied only to factors with 2 or more levels");
    std::map<std::string, int> idx;
    for (size_t l = 0; l < levels[k].size(); ++l) idx[levels[k][l]] = static_cast<int>(l);
    code[k].assign(n, -1);
    for (size_t i = 0; i < n; ++i) {
      size_t r = rows[i];
      if (vd[k].miss[r]) continue;
      auto it = idx.find(vd[k].lab[r]);
      if (it == idx.end()) throw Error("factor " + ts.vars[k] + " has new level " + vd[k].lab[r]);
      code[k][i] = it->second;
    }
  }
  std::vector<Vec> columns;
  std::vector<std::string> names;
  std::vector<int> asg;
  if (f.intercept) {
    columns.emplace_back(n, 1.0);
    names.push_back("(Intercept)");
    asg.push_back(0);
  }
  for (size_t t = 0; t < ts.terms.size(); ++t) {
    const std::vector<int>& term = ts.terms[t];
    // Per variable: the list of levels to take (factor), or {-1} for a numeric variable.
    std::vector<std::vector<int>> pieces;
    for (int v : term) {
      std::vector<int> pc;
      if (!vd[v].factor) {
        pc.push_back(-1);
      } else {
        std::vector<int> margin;
        for (int u : term)
          if (u != v) margin.push_back(u);
        bool spanned = false;
        for (size_t s = 0; s < t && !spanned; ++s) {
          const std::vector<int>& prev = ts.terms[s];
          if (margin.empty()) {
            spanned = prev.size() == 1 && vd[prev[0]].factor;
          } else {
            spanned = std::all_of(margin.begin(), margin.end(), [&](int u) {
              return std::find(prev.begin(), prev.end(), u) != prev.end();
            });
          }
        }
        if (margin.empty() && f.intercept) spanned = true;
        for (size_t l = spanned ? 1 : 0; l < levels[v].size(); ++l) pc.push_back(static_cast<int>(l));
      }
      pieces.push_back(pc);
    }
    std::vector<size_t> pos(term.size(), 0);
    while (true) {
      Vec col(n, 1.0);
      std::string name;
      for (size_t q = 0; q < term.size(); ++q) {
        int v = term[q];
        int lev = pieces[q][pos[q]];
        if (q) name += ":";
        name += ts.vars[v];
        if (lev < 0) {
          for (size_t i = 0; i < n; ++i) col[i] *= vd[v].num[rows[i]];
        } else {
          name += levels[v][static_cast<size_t>(lev)];
          for (size_t i = 0; i < n; ++i) col[i] *= code[v][i] < 0 ? kNaN : (code[v][i] == lev ? 1.0 : 0.0);
        }
      }
      columns.push_back(std::move(col));
      names.push_back(name);
      asg.push_back(static_cast<int>(t + 1));
      size_t q = 0;  // odometer, first variable fastest
      while (q < term.size() && ++pos[q] == pieces[q].size()) pos[q++] = 0;
      if (q == term.size()) break;
    }
  }
  Matrix X(n, columns.size());
  for (size_t j = 0; j < columns.size(); ++j)
    for (size_t i = 0; i < n; ++i) X(i, j) = columns[j][i];
  if (colnames) *colnames = names;
  if (assign) *assign = asg;
  return X;
}

// The model frame of a fit: complete rows (R's na.omit over the response, the variables and the extra vectors),
// factor levels learned from those rows (R's drop.unused.levels), and the model matrix on those rows.
struct ModelFrame {
  Formula f;                       // with xlevels set
  Matrix X;
  Vec y;
  std::vector<size_t> rows;
  std::vector<std::string> names;
  std::vector<int> assign;
};

ModelFrame model_frame(const Formula& f, const Data& d, const std::vector<const Vec*>& extra, const char* fn) {
  ModelFrame mf;
  mf.f = f;
  Vec yall = f.response_values(d);
  size_t N = yall.size();
  TermStructure ts = term_structure(f);
  std::vector<VarData> vd;
  for (const auto& v : ts.vars) vd.push_back(eval_variable(v, d));
  std::vector<VarData> vd_extra;
  for (const auto& v : f.variables)  // variables outside the terms still take part in NA removal
    if (std::find(ts.vars.begin(), ts.vars.end(), v) == ts.vars.end() && v != f.response) {
      VarData extra_var = eval_variable(v, d);
      extra_var.factor = true;  // only the missing flags matter
      if (extra_var.miss.empty()) {
        extra_var.miss.resize(extra_var.num.size());
        for (size_t i = 0; i < extra_var.num.size(); ++i) extra_var.miss[i] = std::isnan(extra_var.num[i]);
      }
      vd_extra.push_back(std::move(extra_var));
    }
  std::vector<bool> mask(N);
  for (size_t i = 0; i < N; ++i) {
    bool ok = !std::isnan(yall[i]);
    for (const VarData& v : vd)
      if (ok) ok = v.factor ? !v.miss[i] : !std::isnan(v.num[i]);
    for (const VarData& v : vd_extra)
      if (ok) ok = !v.miss[i];
    for (const Vec* e : extra)
      if (ok && e && !e->empty()) ok = !std::isnan((*e)[i]);
    mask[i] = ok;
    if (ok) mf.rows.push_back(i);
  }
  if (mf.rows.empty()) throw Error(std::string(fn) + ": 0 (non-NA) cases");
  std::vector<std::vector<std::string>> levels(ts.vars.size());
  for (size_t k = 0; k < ts.vars.size(); ++k) {
    if (!vd[k].factor) continue;
    auto it = f.xlevels.find(ts.vars[k]);
    levels[k] = it != f.xlevels.end() ? it->second : learn_levels(vd[k], mask);
    mf.f.xlevels[ts.vars[k]] = levels[k];
  }
  mf.X = build_design(f, ts, vd, levels, mf.rows, &mf.names, &mf.assign);
  mf.y = take(yall, mf.rows);
  return mf;
}

}  // namespace

// ===================================================================================================================
// Formula
// ===================================================================================================================

Formula Formula::parse(const std::string& text, const Data* data) { return FormulaParser(text, data).run(); }

Vec Formula::response_values(const Data& d) const {
  if (response.empty()) throw Error("formula: the formula has no response");
  std::string inner;
  if (d.is_text(response) || factor_call(response, &inner))
    throw Error("formula: the response '" + response + "' is a factor (text); the response must be numeric");
  Vec y = eval_var(response, d);
  if (y.size() != d.nrow()) throw Error("formula: response length differs from the data");
  return y;
}

Matrix Formula::design(const Data& d, std::vector<std::string>* colnames, std::vector<int>* assign) const {
  size_t n = d.nrow();
  TermStructure ts = term_structure(*this);
  std::vector<VarData> vd;
  for (const auto& v : ts.vars) vd.push_back(eval_variable(v, d));
  std::vector<bool> all(n, true);
  std::vector<std::vector<std::string>> levels(ts.vars.size());
  for (size_t k = 0; k < ts.vars.size(); ++k) {
    if (!vd[k].factor) continue;
    auto it = xlevels.find(ts.vars[k]);
    levels[k] = it != xlevels.end() ? it->second : learn_levels(vd[k], all);
  }
  std::vector<size_t> rows(n);
  for (size_t i = 0; i < n; ++i) rows[i] = i;
  return build_design(*this, ts, vd, levels, rows, colnames, assign);
}

// ===================================================================================================================
// lm
// ===================================================================================================================

Vec LmFit::predict(const Data& newdata) const { return linear_predict(formula, coef, newdata); }

Vec LmFit::predict_se(const Data& newdata) const {
  std::vector<std::string> cn;
  Matrix X = formula.design(newdata, &cn);
  if (cn.size() != coef.names.size()) throw Error("predict: design does not match the fitted model");
  Vec out(X.rows);
  for (size_t i = 0; i < X.rows; ++i) {
    double q = 0;
    for (size_t a = 0; a < X.cols; ++a) {
      if (coef.aliased[a]) continue;
      for (size_t b = 0; b < X.cols; ++b)
        if (!coef.aliased[b]) q += X(i, a) * cov_unscaled(a, b) * X(i, b);
    }
    out[i] = std::sqrt(q) * sigma;
  }
  return out;
}

LmFit lm(const Formula& f, const Data& d, const Vec* weights) {
  size_t N = d.nrow();
  if (weights) {
    if (weights->size() != N) throw Error("lm: 'weights' must have the same length as the data");
    for (double w : *weights)
      if (w < 0 || std::isinf(w)) throw Error("lm: missing or negative weights not allowed");
  }
  ModelFrame mf = model_frame(f, d, {weights}, "lm");
  const std::vector<size_t>& rows = mf.rows;
  const std::vector<std::string>& names = mf.names;
  const Matrix& X = mf.X;
  const Vec& y = mf.y;
  Vec w = weights ? take(*weights, rows) : Vec(rows.size(), 1.0);
  check_finite_design("lm", X, y);
  size_t n = y.size(), p = X.cols;
  if (n == 0) throw Error("lm: 0 (non-NA) cases");

  // Weighted problem: rows scaled by sqrt(w); zero-weight rows are excluded from the fit.
  std::vector<size_t> pos;
  for (size_t i = 0; i < n; ++i)
    if (w[i] > 0) pos.push_back(i);
  Matrix Xw(pos.size(), p);
  Vec yw(pos.size());
  for (size_t k = 0; k < pos.size(); ++k) {
    double s = std::sqrt(w[pos[k]]);
    for (size_t j = 0; j < p; ++j) Xw(k, j) = X(pos[k], j) * s;
    yw[k] = y[pos[k]] * s;
  }
  QR qr(Xw, 1e-7);
  Vec effects = qr.qty(yw);
  Vec b = qr.solve_r(effects);
  Vec est(p, kNaN);
  std::vector<bool> aliased(p, true);
  for (size_t k = 0; k < qr.rank; ++k) {
    est[qr.pivot[k]] = b[k];
    aliased[qr.pivot[k]] = false;
  }

  LmFit fit;
  fit.formula = mf.f;
  fit.assign = mf.assign;
  fit.effects = effects;
  for (size_t k = 0; k < p; ++k) fit.qr_pivot.push_back(static_cast<int>(qr.pivot[k]));
  fit.n = static_cast<int>(pos.size());
  fit.rank = static_cast<int>(qr.rank);
  fit.df_residual = fit.n - fit.rank;
  fit.fitted.assign(n, 0.0);
  fit.residuals.assign(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    double s = 0;
    for (size_t j = 0; j < p; ++j)
      if (!aliased[j]) s += X(i, j) * est[j];
    fit.fitted[i] = s;
    fit.residuals[i] = y[i] - s;
  }
  fit.cov_unscaled = unpivot_cov(qr);

  double rss = 0, sw = 0, swf = 0, slogw = 0;
  for (size_t k : pos) {
    rss += w[k] * fit.residuals[k] * fit.residuals[k];
    sw += w[k];
    swf += w[k] * fit.fitted[k];
    slogw += std::log(w[k]);
  }
  fit.rss = rss;
  double rdf = fit.df_residual;
  double resvar = rdf > 0 ? rss / rdf : kNaN;
  fit.sigma = std::sqrt(resvar);
  fit.coef = make_table(names, est, aliased, fit.cov_unscaled, resvar, rdf, true);

  // ?summary.lm: R^2 = 1 - RSS / sum(w (y - y*)^2), y* the (weighted) mean if there is an intercept, else 0.
  int df_int = f.intercept ? 1 : 0;
  double mss = 0;
  double m = f.intercept ? swf / sw : 0.0;
  for (size_t k : pos) mss += w[k] * (fit.fitted[k] - m) * (fit.fitted[k] - m);
  if (fit.rank > 0) {
    fit.r_squared = mss / (mss + rss);
    fit.adj_r_squared = 1 - (1 - fit.r_squared) * ((fit.n - df_int) / rdf);
  } else {
    fit.r_squared = fit.adj_r_squared = 0;
  }
  if (fit.rank != df_int) {
    fit.f_numdf = fit.rank - df_int;
    fit.f_dendf = rdf;
    fit.fstatistic = (mss / fit.f_numdf) / resvar;
    fit.f_pvalue = (rdf > 0 && std::isfinite(fit.fstatistic)) ? Dist::f(fit.f_numdf, rdf).p(fit.fstatistic, false)
                                                              : kNaN;
  }
  // logLik.lm: 0.5 * (sum(log w) - N (log(2 pi) + 1 - log N + log(sum(w r^2)))), df = rank + 1.
  double Nn = fit.n;
  fit.loglik = 0.5 * (slogw - Nn * (std::log(2 * kPi) + 1 - std::log(Nn) + std::log(rss)));
  double dfll = fit.rank + 1;
  fit.aic = -2 * fit.loglik + 2 * dfll;
  fit.bic = -2 * fit.loglik + std::log(Nn) * dfll;
  return fit;
}

// ===================================================================================================================
// glm
// ===================================================================================================================

namespace {

enum class Link { Identity, Log, Inverse, Logit, Probit, Cloglog, Cauchit, Sqrt };
enum class Fam { Gaussian, Binomial, Poisson, Gamma };

const double kEpsD = kEps;

double qnorm01(double p) { return Dist::norm(0, 1).q(p); }
double pnorm01(double x) { return Dist::norm(0, 1).p(x); }

struct Family {
  Fam fam;
  Link link;

  double linkfun(double mu) const {
    switch (link) {
      case Link::Identity: return mu;
      case Link::Log: return std::log(mu);
      case Link::Inverse: return 1 / mu;
      case Link::Logit: return std::log(mu / (1 - mu));
      case Link::Probit: return qnorm01(mu);
      case Link::Cloglog: return std::log(-std::log1p(-mu));
      case Link::Cauchit: return std::tan(kPi * (mu - 0.5));
      case Link::Sqrt: return std::sqrt(mu);
    }
    return kNaN;
  }
  // Inverse links keep binomial means strictly inside (0, 1) and log means positive (numerical safeguards).
  double linkinv(double eta) const {
    switch (link) {
      case Link::Identity: return eta;
      case Link::Log: return std::max(std::exp(eta), kEpsD);
      case Link::Inverse: return 1 / eta;
      case Link::Logit: {
        double e = std::max(std::min(eta, 30.0), -30.0);
        return 1 / (1 + std::exp(-e));
      }
      case Link::Probit: {
        double th = -qnorm01(kEpsD);
        return pnorm01(std::max(std::min(eta, th), -th));
      }
      case Link::Cloglog: return std::max(std::min(-std::expm1(-std::exp(eta)), 1 - kEpsD), kEpsD);
      case Link::Cauchit: {
        double th = std::tan(kPi * (0.5 - kEpsD));  // -qcauchy(eps)
        double e = std::max(std::min(eta, th), -th);
        return 0.5 + std::atan(e) / kPi;
      }
      case Link::Sqrt: return eta * eta;
    }
    return kNaN;
  }
  double mu_eta(double eta) const {
    switch (link) {
      case Link::Identity: return 1;
      case Link::Log: return std::max(std::exp(eta), kEpsD);
      case Link::Inverse: return -1 / (eta * eta);
      case Link::Logit: {
        double e = std::exp(-std::fabs(eta));
        return std::max(e / ((1 + e) * (1 + e)), kEpsD);
      }
      case Link::Probit: return std::max(std::exp(-0.5 * eta * eta) / std::sqrt(2 * kPi), kEpsD);
      case Link::Cloglog: {
        double e = std::min(eta, 700.0);
        return std::max(std::exp(e - std::exp(e)), kEpsD);
      }
      case Link::Cauchit: return std::max(1 / (kPi * (1 + eta * eta)), kEpsD);
      case Link::Sqrt: return 2 * eta;
    }
    return kNaN;
  }
  bool valideta(double eta) const {
    if (link == Link::Inverse) return std::isfinite(eta) && eta != 0;
    if (link == Link::Sqrt) return std::isfinite(eta) && eta > 0;
    return true;
  }
  bool validmu(double mu) const {
    switch (fam) {
      case Fam::Gaussian: return true;
      case Fam::Binomial: return std::isfinite(mu) && mu > 0 && mu < 1;
      case Fam::Poisson: return std::isfinite(mu) && mu > 0;
      case Fam::Gamma: return std::isfinite(mu) && mu > 0;
    }
    return true;
  }
  double variance(double mu) const {
    switch (fam) {
      case Fam::Gaussian: return 1;
      case Fam::Binomial: return mu * (1 - mu);
      case Fam::Poisson: return mu;
      case Fam::Gamma: return mu * mu;
    }
    return kNaN;
  }
  static double ylogy(double y, double mu) { return y > 0 ? y * std::log(y / mu) : 0.0; }
  double dev_resid(double y, double mu, double wt) const {
    switch (fam) {
      case Fam::Gaussian: return wt * (y - mu) * (y - mu);
      case Fam::Binomial: return 2 * wt * (ylogy(y, mu) + ylogy(1 - y, 1 - mu));
      case Fam::Poisson: return 2 * wt * (y > 0 ? y * std::log(y / mu) - (y - mu) : mu);
      case Fam::Gamma: return -2 * wt * (std::log(y == 0 ? 1.0 : y / mu) - (y - mu) / mu);
    }
    return kNaN;
  }
  // -2 log-likelihood (plus 2 for an estimated dispersion), as R's family()$aic (see ?family).
  double aic(const Vec& y, const Vec& mu, const Vec& wt, double dev) const {
    size_t n = y.size();
    double s = 0;
    switch (fam) {
      case Fam::Gaussian: {
        double nobs = static_cast<double>(n), slw = 0;
        for (double w : wt) slw += std::log(w);
        return nobs * (std::log(2 * kPi * dev / nobs) + 1) + 2 - slw;
      }
      case Fam::Binomial:
        for (size_t i = 0; i < n; ++i) {
          if (!(wt[i] > 0)) continue;
          double m = wt[i];
          s += Dist::binom(std::round(m), mu[i]).d(std::round(m * y[i]), true);
        }
        return -2 * s;
      case Fam::Poisson:
        for (size_t i = 0; i < n; ++i)
          if (wt[i] != 0) s += Dist::pois(mu[i]).d(y[i], true) * wt[i];
        return -2 * s;
      case Fam::Gamma: {
        double sw = 0;
        for (double w : wt) sw += w;
        double disp = dev / sw;
        for (size_t i = 0; i < n; ++i)
          if (wt[i] != 0) s += Dist::gamma(1 / disp, 1 / (mu[i] * disp)).d(y[i], true) * wt[i];
        return -2 * s + 2;
      }
    }
    return kNaN;
  }
  bool dispersion_fixed() const { return fam == Fam::Binomial || fam == Fam::Poisson; }
};

Family make_family(const std::string& fam, const std::string& link_in) {
  struct Spec {
    const char* name;
    Fam fam;
    std::vector<std::string> links;
  };
  static const std::vector<Spec> specs = {
      {"gaussian", Fam::Gaussian, {"identity", "log", "inverse"}},
      {"binomial", Fam::Binomial, {"logit", "probit", "cloglog", "cauchit", "log"}},
      {"poisson", Fam::Poisson, {"log", "identity", "sqrt"}},
      {"Gamma", Fam::Gamma, {"inverse", "identity", "log"}},
  };
  for (const auto& s : specs) {
    if (fam != s.name) continue;
    std::string link = link_in == "default" || link_in.empty() ? s.links[0] : link_in;
    if (std::find(s.links.begin(), s.links.end(), link) == s.links.end()) {
      std::string avail;
      for (size_t i = 0; i < s.links.size(); ++i) avail += (i ? ", '" : "'") + s.links[i] + "'";
      throw Error("glm: link \"" + link + "\" not available for " + fam + " family; available links are " + avail);
    }
    static const std::vector<std::pair<std::string, Link>> lk = {
        {"identity", Link::Identity}, {"log", Link::Log},         {"inverse", Link::Inverse},
        {"logit", Link::Logit},       {"probit", Link::Probit},   {"cloglog", Link::Cloglog},
        {"cauchit", Link::Cauchit},   {"sqrt", Link::Sqrt}};
    for (const auto& kv : lk)
      if (kv.first == link) return Family{s.fam, kv.second};
  }
  throw Error("glm: unknown family '" + fam + "' (expected gaussian, binomial, poisson or Gamma)");
}

struct IrlsResult {
  Vec coef;                    // original order, NaN aliased
  std::vector<bool> aliased;
  size_t rank = 0;
  Vec eta, mu, w_work;         // w_work = squared IRLS weights of the final iteration
  double dev = kNaN;
  int iter = 0;
  bool conv = false;
  Matrix cov;
};

IrlsResult irls(const Family& fam, const Matrix& X, const Vec& y, const Vec& wt, const Vec& off, const Vec& mustart,
                int maxit, double epsilon) {
  size_t n = y.size(), p = X.cols;
  IrlsResult r;
  r.eta.resize(n);
  r.mu.resize(n);
  for (size_t i = 0; i < n; ++i) {
    r.eta[i] = fam.linkfun(mustart[i]);
    r.mu[i] = fam.linkinv(r.eta[i]);
  }
  for (size_t i = 0; i < n; ++i)
    if (!fam.valideta(r.eta[i]) || !fam.validmu(r.mu[i]))
      throw Error("glm: cannot find valid starting values: please specify some");
  auto deviance = [&](const Vec& mu) {
    double s = 0;
    for (size_t i = 0; i < n; ++i) s += fam.dev_resid(y[i], mu[i], wt[i]);
    return s;
  };
  double devold = deviance(r.mu);
  Vec coefold;  // empty until the first successful iteration
  double qtol = std::min(1e-7, epsilon / 1000);
  Vec start(p, 0.0);
  std::unique_ptr<QR> qr;
  for (int iter = 1; iter <= maxit; ++iter) {
    r.iter = iter;
    std::vector<size_t> good;
    Vec z, wsq;
    for (size_t i = 0; i < n; ++i) {
      if (!(wt[i] > 0)) continue;
      double me = fam.mu_eta(r.eta[i]);
      if (me == 0) continue;
      good.push_back(i);
      z.push_back(r.eta[i] - off[i] + (y[i] - r.mu[i]) / me);
      wsq.push_back(std::sqrt(wt[i] * me * me / fam.variance(r.mu[i])));
    }
    if (good.empty()) throw Error("glm: no observations informative at iteration " + std::to_string(iter));
    Matrix Xw(good.size(), p);
    Vec zw(good.size());
    for (size_t k = 0; k < good.size(); ++k) {
      for (size_t j = 0; j < p; ++j) Xw(k, j) = X(good[k], j) * wsq[k];
      zw[k] = z[k] * wsq[k];
    }
    qr.reset(new QR(Xw, qtol));
    Vec b = qr->solve_r(qr->qty(zw));
    std::fill(start.begin(), start.end(), 0.0);
    for (size_t k = 0; k < qr->rank; ++k) start[qr->pivot[k]] = b[k];
    for (double v : start)
      if (!std::isfinite(v)) throw Error("glm: non-finite coefficients at iteration " + std::to_string(iter));
    r.w_work.assign(n, 0.0);
    for (size_t k = 0; k < good.size(); ++k) r.w_work[good[k]] = wsq[k] * wsq[k];
    auto update = [&]() {
      for (size_t i = 0; i < n; ++i) {
        double s = off[i];
        for (size_t j = 0; j < p; ++j) s += X(i, j) * start[j];
        r.eta[i] = s;
        r.mu[i] = fam.linkinv(s);
      }
      return deviance(r.mu);
    };
    double dev = update();
    // Step halving when the deviance is not finite or the fit leaves the valid region (McCullagh & Nelder 2.5).
    auto invalid = [&]() {
      for (size_t i = 0; i < n; ++i)
        if (!fam.valideta(r.eta[i]) || !fam.validmu(r.mu[i])) return true;
      return false;
    };
    for (int ii = 0; !std::isfinite(dev) || invalid(); ++ii) {
      if (coefold.empty())
        throw Error("glm: no valid set of coefficients has been found: please supply starting values");
      if (ii >= maxit) throw Error("glm: inner loop; cannot correct step size");
      for (size_t j = 0; j < p; ++j) start[j] = (start[j] + coefold[j]) / 2;
      dev = update();
    }
    r.dev = dev;
    if (std::fabs(dev - devold) / (std::fabs(dev) + 0.1) < epsilon) {
      r.conv = true;
      break;
    }
    devold = dev;
    coefold = start;
  }
  r.rank = qr->rank;
  r.coef.assign(p, kNaN);
  r.aliased.assign(p, true);
  for (size_t k = 0; k < qr->rank; ++k) {
    r.coef[qr->pivot[k]] = start[qr->pivot[k]];
    r.aliased[qr->pivot[k]] = false;
  }
  r.cov = unpivot_cov(*qr);
  return r;
}

}  // namespace

namespace {
// Family initialisation and starting values (?family / ?glm); binomial y with zero weight is set to 0.
Vec glm_start(const Family& fam, Vec& y, const Vec& wt) {
  size_t n = y.size();
  Vec mustart(n);
  for (size_t i = 0; i < n; ++i) {
    switch (fam.fam) {
      case Fam::Gaussian: mustart[i] = y[i]; break;
      case Fam::Binomial:
        if (y[i] < 0 || y[i] > 1) throw Error("glm: y values must be 0 <= y <= 1");
        if (wt[i] == 0) y[i] = 0;
        mustart[i] = (wt[i] * y[i] + 0.5) / (wt[i] + 1);
        break;
      case Fam::Poisson:
        if (y[i] < 0) throw Error("glm: negative values not allowed for the 'Poisson' family");
        mustart[i] = y[i] + 0.1;
        break;
      case Fam::Gamma:
        if (y[i] <= 0) throw Error("glm: non-positive values not allowed for the 'Gamma' family");
        mustart[i] = y[i];
        break;
    }
  }
  return mustart;
}
}  // namespace

Vec GlmFit::predict(const Data& newdata, const std::string& type) const {
  Vec eta = linear_predict(formula, coef, newdata);
  if (type == "link") return eta;
  if (type != "response") throw Error("predict: type must be \"link\" or \"response\"");
  Family fam = make_family(family, link);
  for (double& e : eta) e = std::isnan(e) ? kNaN : fam.linkinv(e);
  return eta;
}

GlmFit glm(const Formula& f, const Data& d, const GlmOptions& opt) {
  Family fam = make_family(opt.family, opt.link);
  if (!(opt.epsilon > 0)) throw Error("glm: value of 'epsilon' must be > 0");
  if (opt.maxit < 1) throw Error("glm: maximum number of iterations must be > 0");
  size_t N = d.nrow();
  if (!opt.weights.empty()) {
    if (opt.weights.size() != N) throw Error("glm: 'weights' must have the same length as the data");
    for (double w : opt.weights)
      if (w < 0 || std::isinf(w)) throw Error("glm: negative weights not allowed");
  }
  if (!opt.offset.empty() && opt.offset.size() != N)
    throw Error("glm: number of offsets is " + std::to_string(opt.offset.size()) + ", should equal " +
                std::to_string(N) + " (number of observations)");
  ModelFrame mf = model_frame(f, d, {&opt.weights, &opt.offset}, "glm");
  const std::vector<size_t>& rows = mf.rows;
  const std::vector<std::string>& names = mf.names;
  const Matrix& X = mf.X;
  Vec y = mf.y;
  Vec wt = opt.weights.empty() ? Vec(rows.size(), 1.0) : take(opt.weights, rows);
  Vec off = opt.offset.empty() ? Vec(rows.size(), 0.0) : take(opt.offset, rows);
  check_finite_design("glm", X, y);
  size_t n = y.size();
  Vec mustart = glm_start(fam, y, wt);
  IrlsResult r = irls(fam, X, y, wt, off, mustart, opt.maxit, opt.epsilon);

  GlmFit g;
  g.formula = mf.f;
  g.assign = mf.assign;
  g.model_matrix = X;
  g.y = y;
  g.prior_weights = wt;
  g.offset = off;
  g.maxit = opt.maxit;
  g.epsilon = opt.epsilon;
  g.family = opt.family;
  g.link = opt.link == "default" || opt.link.empty()
               ? std::string(fam.fam == Fam::Gaussian ? "identity"
                             : fam.fam == Fam::Binomial ? "logit"
                             : fam.fam == Fam::Poisson  ? "log"
                                                        : "inverse")
               : opt.link;
  int n_ok = 0;
  for (double w : wt)
    if (w != 0) ++n_ok;
  g.n = n_ok;
  g.rank = static_cast<int>(r.rank);
  g.df_residual = n_ok - g.rank;
  g.df_null = n_ok - (f.intercept ? 1 : 0);
  g.iter = r.iter;
  g.converged = r.conv;
  g.deviance = r.dev;
  g.fitted = r.mu;
  g.linear_predictors = r.eta;
  g.weights_working = r.w_work;
  g.cov_unscaled = r.cov;
  g.residuals_deviance.resize(n);
  g.residuals_pearson.resize(n);
  g.residuals_working.resize(n);
  for (size_t i = 0; i < n; ++i) {
    double mu = r.mu[i];
    double dr = fam.dev_resid(y[i], mu, wt[i]);
    g.residuals_deviance[i] = (y[i] > mu ? 1.0 : (y[i] < mu ? -1.0 : 0.0)) * std::sqrt(std::max(dr, 0.0));
    g.residuals_pearson[i] = (y[i] - mu) * std::sqrt(wt[i]) / std::sqrt(fam.variance(mu));
    g.residuals_working[i] = (y[i] - mu) / fam.mu_eta(r.eta[i]);
  }
  // Null deviance: intercept-only fit (closed form without offset), or the zero linear predictor without intercept.
  bool has_off = false;
  for (double o : off)
    if (o != 0) has_off = true;
  if (f.intercept && has_off) {
    Matrix one(n, 1, 1.0);
    g.null_deviance = irls(fam, one, y, wt, off, mustart, opt.maxit, opt.epsilon).dev;
  } else {
    double sw = 0, swy = 0;
    for (size_t i = 0; i < n; ++i) {
      sw += wt[i];
      swy += wt[i] * y[i];
    }
    double nd = 0;
    for (size_t i = 0; i < n; ++i) {
      double mu0 = f.intercept ? swy / sw : fam.linkinv(off[i]);
      nd += fam.dev_resid(y[i], mu0, wt[i]);
    }
    g.null_deviance = nd;
  }
  // AIC over the observations with non-zero weight.
  {
    Vec ya, ma, wa;
    for (size_t i = 0; i < n; ++i)
      if (wt[i] != 0) {
        ya.push_back(y[i]);
        ma.push_back(r.mu[i]);
        wa.push_back(wt[i]);
      }
    g.aic = fam.aic(ya, ma, wa, g.deviance) + 2.0 * g.rank;
  }
  // Dispersion: 1 for binomial/poisson, otherwise Pearson chi-square / residual df.
  if (fam.dispersion_fixed()) {
    g.dispersion = 1;
  } else if (g.df_residual > 0) {
    double s = 0;
    for (size_t i = 0; i < n; ++i)
      if (wt[i] > 0) s += r.w_work[i] * g.residuals_working[i] * g.residuals_working[i];
    g.dispersion = s / g.df_residual;
  } else {
    g.dispersion = kNaN;
  }
  g.coef = make_table(names, r.coef, r.aliased, g.cov_unscaled, g.dispersion, g.df_residual,
                      !fam.dispersion_fixed());
  return g;
}

// ===================================================================================================================
// Analysis of variance and deviance tables
//
// Sequential (type I) sums of squares from the orthogonal effects Q'y of the QR decomposition: the squared effects of
// the columns of each term, in model order (Chambers & Hastie 1992, "Statistical Models in S", section 5.2; Seber &
// Lee, "Linear Regression Analysis", 2nd ed., section 4.5); analysis of deviance by refitting the nested sequence of models
// (McCullagh & Nelder 1989, section 2.3.2 and 4.4); layout and column names as documented in ?anova.lm and ?anova.glm.
// ===================================================================================================================

const Vec& AnovaTable::col(const std::string& name) const {
  for (size_t j = 0; j < col_names.size(); ++j)
    if (col_names[j] == name) return cols[j];
  throw Error("anova: no column '" + name + "'");
}

double AnovaTable::at(const std::string& row, const std::string& column) const {
  const Vec& c = col(column);
  for (size_t i = 0; i < row_names.size(); ++i)
    if (row_names[i] == row) return c[i];
  throw Error("anova: no row '" + row + "'");
}

namespace {

std::string formula_text(const Formula& f) {
  std::string s = f.response + " ~ ";
  if (f.term_labels.empty()) return s + (f.intercept ? "1" : "0");
  for (size_t i = 0; i < f.term_labels.size(); ++i) s += (i ? " + " : "") + f.term_labels[i];
  if (!f.intercept) s += " - 1";
  return s;
}

double pf_upper(double F, double df1, double df2) {
  if (std::isnan(F) || !(df1 > 0) || !(df2 > 0)) return kNaN;
  if (std::isinf(df2)) return Dist::chisq(df1).p(F * df1, false);
  return Dist::f(df1, df2).p(F, false);
}

}  // namespace

AnovaTable anova(const LmFit& fit) {
  size_t nt = fit.formula.term_labels.size();
  Vec ss(nt + 1, 0.0), df(nt + 1, 0.0);
  for (int k = 0; k < fit.rank; ++k) {
    int a = fit.assign[static_cast<size_t>(fit.qr_pivot[static_cast<size_t>(k)])];
    ss[static_cast<size_t>(a)] += fit.effects[static_cast<size_t>(k)] * fit.effects[static_cast<size_t>(k)];
    df[static_cast<size_t>(a)] += 1;
  }
  AnovaTable t;
  t.heading = "Analysis of Variance Table\n\nResponse: " + fit.formula.response + "\n";
  t.col_names = {"Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"};
  t.cols.assign(5, Vec());
  double rdf = fit.df_residual;
  double resvar = rdf > 0 ? fit.rss / rdf : kNaN;
  for (size_t a = 1; a <= nt; ++a) {
    if (df[a] == 0) continue;  // terms with no estimable columns are dropped, as R does
    double ms = ss[a] / df[a];
    double F = ms / resvar;
    t.row_names.push_back(fit.formula.term_labels[a - 1]);
    t.cols[0].push_back(df[a]);
    t.cols[1].push_back(ss[a]);
    t.cols[2].push_back(ms);
    t.cols[3].push_back(F);
    t.cols[4].push_back(pf_upper(F, df[a], rdf));
  }
  t.row_names.push_back("Residuals");
  t.cols[0].push_back(rdf);
  t.cols[1].push_back(fit.rss);
  t.cols[2].push_back(resvar);
  t.cols[3].push_back(kNaN);
  t.cols[4].push_back(kNaN);
  return t;
}

AnovaTable anova(const GlmFit& g, const std::string& test_in) {
  std::string test = test_in == "LRT" ? "Chisq" : test_in;
  if (test != "none" && !test.empty() && test != "Chisq" && test != "F")
    throw Error("anova: test must be \"none\", \"Chisq\", \"LRT\" or \"F\" (\"Rao\" and \"Cp\" are not supported)");
  Family fam = make_family(g.family, g.link);
  size_t nt = g.formula.term_labels.size();
  int n_ok = 0;
  for (double w : g.prior_weights)
    if (w != 0) ++n_ok;
  Vec resid_df(nt + 1), resid_dev(nt + 1);
  resid_df[0] = g.df_null;
  resid_dev[0] = g.null_deviance;
  for (size_t i = 1; i < nt; ++i) {
    std::vector<size_t> keep;
    for (size_t j = 0; j < g.assign.size(); ++j)
      if (g.assign[j] <= static_cast<int>(i)) keep.push_back(j);
    Matrix Xs(g.model_matrix.rows, keep.size());
    for (size_t r = 0; r < Xs.rows; ++r)
      for (size_t c = 0; c < keep.size(); ++c) Xs(r, c) = g.model_matrix(r, keep[c]);
    Vec y = g.y;
    Vec mustart = glm_start(fam, y, g.prior_weights);
    IrlsResult r = irls(fam, Xs, y, g.prior_weights, g.offset, mustart, g.maxit, g.epsilon);
    resid_df[i] = n_ok - static_cast<double>(r.rank);
    resid_dev[i] = r.dev;
  }
  if (nt > 0) {
    resid_df[nt] = g.df_residual;
    resid_dev[nt] = g.deviance;
  }
  AnovaTable t;
  t.heading = "Analysis of Deviance Table\n\nModel: " + g.family + ", link: " + g.link + "\n\nResponse: " +
              g.formula.response + "\n\nTerms added sequentially (first to last)\n\n";
  t.col_names = {"Df", "Deviance", "Resid. Df", "Resid. Dev"};
  t.cols.assign(4, Vec());
  t.row_names.push_back("NULL");
  for (size_t i = 0; i <= nt; ++i) {
    if (i > 0) t.row_names.push_back(g.formula.term_labels[i - 1]);
    t.cols[0].push_back(i ? resid_df[i - 1] - resid_df[i] : kNaN);
    t.cols[1].push_back(i ? resid_dev[i - 1] - resid_dev[i] : kNaN);
    t.cols[2].push_back(resid_df[i]);
    t.cols[3].push_back(resid_dev[i]);
  }
  if (test == "Chisq" || test == "F") {
    double disp = g.dispersion;
    double df_disp = fam.dispersion_fixed() ? kInf : static_cast<double>(g.df_residual);
    Vec fv, pv;
    for (size_t i = 0; i <= nt; ++i) {
      double df = t.cols[0][i], dev = t.cols[1][i];
      if (i == 0 || !(df > 0)) {
        fv.push_back(kNaN);
        pv.push_back(kNaN);
        continue;
      }
      if (test == "Chisq") {
        pv.push_back(Dist::chisq(df).p(dev / disp, false));
      } else {
        double F = dev / df / disp;
        fv.push_back(F);
        pv.push_back(pf_upper(F, df, df_disp));
      }
    }
    if (test == "F") {
      t.col_names.push_back("F");
      t.cols.push_back(fv);
    }
    t.col_names.push_back(test == "F" ? "Pr(>F)" : "Pr(>Chi)");
    t.cols.push_back(pv);
  }
  return t;
}

AnovaTable anova(const LmFit& m1, const LmFit& m2) {
  if (m1.n != m2.n) throw Error("anova: models were not all fitted to the same size of dataset");
  const LmFit& big = m2.df_residual <= m1.df_residual ? m2 : m1;
  double scale = big.rss / big.df_residual;
  AnovaTable t;
  t.heading = "Analysis of Variance Table\n\nModel 1: " + formula_text(m1.formula) + "\nModel 2: " +
              formula_text(m2.formula) + "\n";
  t.row_names = {"1", "2"};
  t.col_names = {"Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)"};
  double df = static_cast<double>(m1.df_residual) - m2.df_residual;
  double ss = m1.rss - m2.rss;
  double F = df != 0 ? ss / df / scale : kNaN;
  double p = df != 0 ? pf_upper(F, std::fabs(df), big.df_residual) : kNaN;
  t.cols = {{static_cast<double>(m1.df_residual), static_cast<double>(m2.df_residual)},
            {m1.rss, m2.rss},
            {kNaN, df},
            {kNaN, ss},
            {kNaN, F},
            {kNaN, p}};
  return t;
}

// ===================================================================================================================
// nls
// ===================================================================================================================

namespace {

void eval_model(const NlsModel& model, const Vec& par, Vec& out, size_t n) {
  out.assign(n, 0.0);
  model(par, out);
  if (out.size() != n) throw Error("nls: the model returned " + std::to_string(out.size()) + " values, expected " +
                                   std::to_string(n));
}

double rss_of(const Vec& y, const Vec& f) {
  double s = 0;
  for (size_t i = 0; i < y.size(); ++i) s += (y[i] - f[i]) * (y[i] - f[i]);
  return s;
}

bool all_finite(const Vec& v) {
  for (double x : v)
    if (!std::isfinite(x)) return false;
  return true;
}

// Forward-difference Jacobian (n x p) as R's numericDeriv: h = sqrt(eps) * |theta|, or sqrt(eps) when theta = 0;
// the difference is divided by the nominal h. `dir` gives -1 for parameters that must step downwards (upper bound).
Matrix jacobian(const NlsModel& model, const Vec& par, const Vec& f0, const std::vector<int>& dir) {
  size_t n = f0.size(), p = par.size();
  Matrix J(n, p);
  Vec pp = par, fv;
  const double se = std::sqrt(kEps);
  for (size_t j = 0; j < p; ++j) {
    double h = par[j] == 0 ? se : se * std::fabs(par[j]);
    if (!dir.empty() && dir[j] < 0) h = -h;
    pp[j] = par[j] + h;
    eval_model(model, pp, fv, n);
    pp[j] = par[j];
    for (size_t i = 0; i < n; ++i) J(i, j) = (fv[i] - f0[i]) / h;
  }
  for (double v : J.a)
    if (!std::isfinite(v)) throw Error("nls: Missing value or an infinity produced when evaluating the model");
  return J;
}

std::string fmt_g(double x) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%g", x);
  return buf;
}

void finish_nls(NlsFit& fit, const std::vector<std::string>& names, const Vec& par, const Matrix& J, const Vec& f,
                const Vec& y) {
  size_t n = y.size(), p = par.size();
  fit.n = static_cast<int>(n);
  fit.df_residual = static_cast<int>(n) - static_cast<int>(p);
  fit.fitted = f;
  fit.residuals.resize(n);
  for (size_t i = 0; i < n; ++i) fit.residuals[i] = y[i] - f[i];
  fit.rss = rss_of(y, f);
  double resvar = fit.df_residual > 0 ? fit.rss / fit.df_residual : kNaN;
  fit.sigma = std::sqrt(resvar);
  QR qr(J, 1e-7);
  if (qr.rank < p) throw Error("nls: singular gradient matrix at parameter estimates");
  fit.cov_unscaled = unpivot_cov(qr);
  fit.coef = make_table(names, par, std::vector<bool>(p, false), fit.cov_unscaled, resvar, fit.df_residual, true);
}

// Relative offset criterion (Bates & Watts 1988, 2.2.3): sqrt(|Q1'r|^2 / |Q2'r|^2) — the projection of the residual
// onto the tangent plane compared with its orthogonal part.
double relative_offset(const QR& qr, const Vec& resid) {
  Vec q = qr.qty(resid);
  double a = 0, b = 0;
  for (size_t i = 0; i < q.size(); ++i) (i < qr.rank ? a : b) += q[i] * q[i];
  return std::sqrt(a / b);
}

NlsFit nls_gauss_newton(const NlsModel& model, const Vec& y, const std::vector<std::string>& names, const Vec& start,
                        const NlsOptions& opt) {
  size_t n = y.size(), p = start.size();
  Vec par = start, f, fnew, resid(n);
  eval_model(model, par, f, n);
  if (!all_finite(f)) throw Error("nls: Missing value or an infinity produced when evaluating the model");
  double rss = rss_of(y, f);
  NlsFit fit;
  fit.message = "converged";
  double fac = 1.0;
  int iter = 0;
  Matrix J;
  for (;;) {
    // R's order (established black-box): the iteration limit is tested before the criterion is recomputed, so the
    // reported criterion is the one from before the last step.
    if (iter >= opt.maxiter) {
      fit.message = "number of iterations exceeded maximum of " + std::to_string(opt.maxiter);
      break;
    }
    for (size_t i = 0; i < n; ++i) resid[i] = y[i] - f[i];
    J = jacobian(model, par, f, {});
    QR qr(J, 1e-7);
    if (qr.rank < p)
      throw Error(iter == 0 ? "nls: singular gradient matrix at initial parameter estimates" : "nls: singular gradient");
    fit.conv_criterion = relative_offset(qr, resid);
    if (fit.conv_criterion <= opt.tol) {  // converged
      fit.converged = true;
      break;
    }
    Vec qtr = qr.qty(resid);
    Vec bp = qr.solve_r(qtr);
    Vec incr(p, 0.0);
    for (size_t k = 0; k < qr.rank; ++k) incr[qr.pivot[k]] = bp[k];
    bool accepted = false;
    Vec newpar(p);
    while (fac >= opt.min_factor) {
      for (size_t j = 0; j < p; ++j) newpar[j] = par[j] + fac * incr[j];
      eval_model(model, newpar, fnew, n);
      double newrss = rss_of(y, fnew);
      if (std::isfinite(newrss) && newrss <= rss) {
        par = newpar;
        f = fnew;
        rss = newrss;
        fac = std::min(2 * fac, 1.0);
        accepted = true;
        ++iter;
        break;
      }
      fac /= 2;
    }
    if (!accepted) {
      // As R with warnOnly = TRUE, the parameters are left at the last trial point.
      par = newpar;
      f = fnew;
      fit.message = "step factor " + fmt_g(fac) + " reduced below 'minFactor' of " + fmt_g(opt.min_factor);
      break;
    }
  }
  fit.iterations = iter;
  if (!fit.converged) J = jacobian(model, par, f, {});
  finish_nls(fit, names, par, J, f, y);
  return fit;
}

// Bounded least squares: projected Levenberg-Marquardt (Nocedal & Wright 10.3) with clipping to the box.
NlsFit nls_port(const NlsModel& model, const Vec& y, const std::vector<std::string>& names, const Vec& start,
                const NlsOptions& opt) {
  size_t n = y.size(), p = start.size();
  Vec lo = opt.lower.empty() ? Vec(p, -kInf) : opt.lower;
  Vec hi = opt.upper.empty() ? Vec(p, kInf) : opt.upper;
  if (lo.size() == 1 && p > 1) lo.assign(p, lo[0]);
  if (hi.size() == 1 && p > 1) hi.assign(p, hi[0]);
  if (lo.size() != p || hi.size() != p) throw Error("nls: 'lower' and 'upper' must have one value per parameter");
  for (size_t j = 0; j < p; ++j)
    if (!(lo[j] <= hi[j])) throw Error("nls: 'lower' must not exceed 'upper'");
  for (size_t j = 0; j < p; ++j)
    if (start[j] < lo[j] || start[j] > hi[j]) throw Error("nls: Convergence failure: initial par violates constraints");
  Vec par = start;
  Vec f, fnew;
  eval_model(model, par, f, n);
  if (!all_finite(f)) throw Error("nls: Missing value or an infinity produced when evaluating the model");
  double rss = rss_of(y, f);
  NlsFit fit;
  fit.message = "converged";
  double lambda = 1e-3;
  int iter = 0;
  const int maxit = std::max(opt.maxiter, 1) * 4;
  auto direction = [&](const Vec& par_) {
    std::vector<int> dir(p, 1);
    for (size_t j = 0; j < p; ++j) {
      double h = par_[j] == 0 ? std::sqrt(kEps) : std::sqrt(kEps) * std::fabs(par_[j]);
      if (par_[j] + h > hi[j]) dir[j] = -1;
    }
    return dir;
  };
  Matrix J = jacobian(model, par, f, direction(par));
  bool converged = false;
  int small_steps = 0;
  while (iter < maxit) {
    ++iter;
    Vec r(n);
    for (size_t i = 0; i < n; ++i) r[i] = y[i] - f[i];
    // gradient of rss/2 is -J'r
    Vec g(p, 0.0);
    for (size_t j = 0; j < p; ++j)
      for (size_t i = 0; i < n; ++i) g[j] -= J(i, j) * r[i];
    // free set: not at a bound where the negative gradient points outward
    std::vector<size_t> freev;
    for (size_t j = 0; j < p; ++j) {
      bool at_lo = par[j] <= lo[j] && g[j] > 0;
      bool at_hi = par[j] >= hi[j] && g[j] < 0;
      if (!at_lo && !at_hi) freev.push_back(j);
    }
    if (freev.empty()) {
      converged = true;
      break;
    }
    // Normal equations with Marquardt scaling on the free set: (A + lambda diag(A)) d = J'r
    size_t m = freev.size();
    Matrix A(m, m, 0.0);
    Vec rhs(m, 0.0);
    for (size_t a = 0; a < m; ++a) {
      for (size_t i = 0; i < n; ++i) rhs[a] += J(i, freev[a]) * r[i];
      for (size_t b = 0; b < m; ++b)
        for (size_t i = 0; i < n; ++i) A(a, b) += J(i, freev[a]) * J(i, freev[b]);
    }
    bool accepted = false;
    for (int tries = 0; tries < 60 && !accepted; ++tries) {
      // Solve via QR of the augmented system [J_free; sqrt(lambda D)] d = [r; 0] (numerically stable).
      Matrix Aug(n + m, m, 0.0);
      Vec rr(n + m, 0.0);
      for (size_t i = 0; i < n; ++i) {
        for (size_t a = 0; a < m; ++a) Aug(i, a) = J(i, freev[a]);
        rr[i] = r[i];
      }
      for (size_t a = 0; a < m; ++a) Aug(n + a, a) = std::sqrt(lambda * std::max(A(a, a), 1e-300));
      QR qa(Aug, 1e-12);
      Vec ds = qa.solve_r(qa.qty(rr));
      Vec d(m, 0.0);
      for (size_t k = 0; k < qa.rank; ++k) d[qa.pivot[k]] = ds[k];
      Vec newpar = par;
      for (size_t a = 0; a < m; ++a) {
        size_t j = freev[a];
        newpar[j] = std::min(std::max(par[j] + d[a], lo[j]), hi[j]);
      }
      eval_model(model, newpar, fnew, n);
      double newrss = rss_of(y, fnew);
      if (std::isfinite(newrss) && newrss <= rss) {
        double step = 0;
        for (size_t j = 0; j < p; ++j)
          step = std::max(step, std::fabs(newpar[j] - par[j]) / std::max(std::fabs(par[j]), 1e-10));
        double rel = (rss - newrss) / std::max(rss, 1e-300);
        par = newpar;
        f = fnew;
        rss = newrss;
        lambda = std::max(lambda / 10, 1e-12);
        accepted = true;
        if (rel < 1e-15 && step < 1e-10) ++small_steps;
        else small_steps = 0;
      } else {
        lambda *= 10;
      }
    }
    if (!accepted || small_steps >= 2) {
      converged = true;
      break;
    }
    J = jacobian(model, par, f, direction(par));
  }
  fit.iterations = iter;
  fit.converged = converged;
  if (!converged) fit.message = "iteration limit reached without convergence";
  J = jacobian(model, par, f, direction(par));
  {
    QR qr(J, 1e-7);
    Vec r(n);
    for (size_t i = 0; i < n; ++i) r[i] = y[i] - f[i];
    if (qr.rank == p) fit.conv_criterion = relative_offset(qr, r);
  }
  finish_nls(fit, names, par, J, f, y);
  return fit;
}

}  // namespace

NlsFit nls(const NlsModel& model, const Vec& y, const std::vector<std::string>& names, const Vec& start,
           const NlsOptions& opt) {
  if (!model) throw Error("nls: no model function");
  if (start.empty()) throw Error("nls: no starting values specified");
  if (names.size() != start.size()) throw Error("nls: 'names' and 'start' must have the same length");
  for (double s : start)
    if (!std::isfinite(s)) throw Error("nls: missing or non-finite starting values");
  for (double v : y)
    if (!std::isfinite(v)) throw Error("nls: missing or non-finite response values");
  if (y.size() < start.size()) throw Error("nls: fewer observations than parameters");
  if (y.size() == start.size()) throw Error("nls: number of observations equals the number of parameters");
  if (!(opt.tol > 0)) throw Error("nls: 'tol' must be positive");
  if (!(opt.min_factor > 0)) throw Error("nls: 'minFactor' must be positive");
  if (opt.algorithm == "default") return nls_gauss_newton(model, y, names, start, opt);
  if (opt.algorithm == "port") return nls_port(model, y, names, start, opt);
  throw Error("nls: unknown algorithm '" + opt.algorithm + "' (expected \"default\" or \"port\")");
}

}  // namespace rcli
