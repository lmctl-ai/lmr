// Hypothesis tests (a Test builder that runs into an Htest result) and the discrete Fourier transform (CVec).
#include <algorithm>
#include <cstdio>
#include <optional>

#include "host.hpp"
#include "rcli/dist.hpp"
#include "rcli/fft.hpp"
#include "rcli/htest.hpp"

namespace rcli::host {

namespace {

struct TestSpec {
  std::string kind;
  std::optional<Vec> x, y, counts, probs, successes, trials, p_vec, groups, blocks;
  std::optional<Dist> dist;
  std::optional<Matrix> table;
  TestOptions o;
  bool exact_set = false;
};

const char* kKinds = "t, wilcox, var, cor, shapiro, ks, chisq, fisher, binom, prop, kruskal, bartlett, fligner, oneway, "
                     "friedman, quade, mcnemar, ansari";

void need(const std::optional<Vec>& v, const std::string& kind, const char* what) {
  if (!v) throw dslvm::Error("test(\"" + kind + "\"): call ." + what + "(...) before run()");
}

// Blocked designs: .x(y).groups(g).blocks(b), or .table(values, nblocks, byrow) with blocks as rows and groups as
// columns (R's friedman.test(matrix)).
void block_design(const TestSpec& s, const std::string& k, Vec& y, Vec& g, Vec& bl) {
  if (s.table) {
    const Matrix& m = *s.table;
    for (size_t i = 0; i < m.rows; ++i)
      for (size_t j = 0; j < m.cols; ++j) {
        y.push_back(m(i, j));
        g.push_back(static_cast<double>(j));
        bl.push_back(static_cast<double>(i));
      }
    return;
  }
  need(s.x, k, "x");
  need(s.groups, k, "groups");
  need(s.blocks, k, "blocks (or .table(values, nblocks, byrow))");
  y = *s.x;
  g = *s.groups;
  bl = *s.blocks;
}

HtestResult run_test(const TestSpec& s) {
  const std::string& k = s.kind;
  if (k == "t") { need(s.x, k, "x"); return s.y ? t_test(*s.x, *s.y, s.o) : t_test(*s.x, s.o); }
  if (k == "wilcox") { need(s.x, k, "x"); return s.y ? wilcox_test(*s.x, *s.y, s.o) : wilcox_test(*s.x, s.o); }
  if (k == "var") { need(s.x, k, "x"); need(s.y, k, "y"); return var_test(*s.x, *s.y, s.o); }
  if (k == "cor") { need(s.x, k, "x"); need(s.y, k, "y"); return cor_test(*s.x, *s.y, s.o); }
  if (k == "shapiro") { need(s.x, k, "x"); return shapiro_test(*s.x); }
  if (k == "ks") {
    need(s.x, k, "x");
    if (s.dist) return ks_test(*s.x, *s.dist, s.o);
    need(s.y, k, "y (or .dist(...))");
    return ks_test(*s.x, *s.y, s.o);
  }
  if (k == "chisq") {
    if (s.table) return chisq_test(*s.table, s.o);
    need(s.counts, k, "counts (or .table(...))");
    Vec p = s.probs ? *s.probs : Vec(s.counts->size(), 1.0 / static_cast<double>(s.counts->size()));
    return chisq_test_gof(*s.counts, p);
  }
  if (k == "fisher") {
    if (!s.table) throw dslvm::Error("test(\"fisher\"): call .table(values, nrow, byrow) before run()");
    return fisher_test(*s.table, s.o);
  }
  if (k == "binom") {
    need(s.successes, k, "successes");
    need(s.trials, k, "trials");
    if (s.successes->size() != 1 || s.trials->size() != 1) throw dslvm::Error("test(\"binom\"): successes and trials take one number each");
    return binom_test((*s.successes)[0], (*s.trials)[0], s.o);
  }
  if (k == "prop") {
    need(s.successes, k, "successes");
    need(s.trials, k, "trials");
    return prop_test(*s.successes, *s.trials, s.p_vec ? *s.p_vec : Vec{}, s.o);
  }
  if (k == "kruskal" || k == "bartlett" || k == "fligner" || k == "oneway") {
    need(s.x, k, "x");
    need(s.groups, k, "groups");
    if (k == "kruskal") return kruskal_test(*s.x, *s.groups);
    if (k == "bartlett") return bartlett_test(*s.x, *s.groups);
    if (k == "fligner") return fligner_test(*s.x, *s.groups);
    return oneway_test(*s.x, *s.groups, s.o.var_equal);
  }
  if (k == "friedman" || k == "quade") {
    Vec y, g, bl;
    block_design(s, k, y, g, bl);
    return k == "friedman" ? friedman_test(y, g, bl) : quade_test(y, g, bl);
  }
  if (k == "mcnemar") {
    if (!s.table) throw dslvm::Error("test(\"mcnemar\"): call .table(values, nrow, byrow) before run()");
    return mcnemar_test(*s.table, s.o.correct);
  }
  if (k == "ansari") { need(s.x, k, "x"); need(s.y, k, "y"); return ansari_test(*s.x, *s.y, s.o); }
  throw dslvm::Error("test: unknown kind '" + k + "' (kinds: " + kKinds + ")");
}

Vec level_codes(const Strs& labels) {
  Strs levels(labels);
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  Vec out;
  for (const auto& l : labels) out.push_back(static_cast<double>(std::lower_bound(levels.begin(), levels.end(), l) - levels.begin()));
  return out;
}

std::string num(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

std::string summary_text(const HtestResult& r) {
  std::string out = r.method + "\n";
  std::string line;
  if (!r.statistic_name.empty()) line += r.statistic_name + " = " + num(r.statistic);
  for (const auto& p : r.parameters) line += (line.empty() ? "" : ", ") + p.first + " = " + num(p.second);
  line += (line.empty() ? "" : ", ") + std::string("p-value = ") + num(r.p_value);
  out += line + "\n";
  if (!r.null_name.empty()) {
    std::string rel = r.alternative == "less" ? "less than" : r.alternative == "greater" ? "greater than" : "not equal to";
    out += "alternative hypothesis: true " + r.null_name + " is " + rel + " " + num(r.null_value) + "\n";
  } else {
    out += "alternative hypothesis: " + r.alternative + "\n";
  }
  if (r.has_conf_int) out += num(100 * r.conf_level) + " percent confidence interval: " + num(r.conf_low) + " " + num(r.conf_high) + "\n";
  for (const auto& e : r.estimates) out += e.first + ": " + num(e.second) + "\n";
  return out;
}

double param_named(const HtestResult& r, const std::string& name) {
  for (const auto& p : r.parameters) if (p.first == name) return p.second;
  std::string have;
  for (const auto& p : r.parameters) have += (have.empty() ? "" : ", ") + p.first;
  throw dslvm::Error("Htest: no parameter '" + name + "' (" + (have.empty() ? std::string("this test has none") : "parameters: " + have) + ")");
}

}  // namespace

void bind_stats(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  // ---- Test builder ----
  b.fn("test", {"string"}, "Test", [&h](A& a) {
    static const char* kinds[] = {"t", "wilcox", "var", "cor", "shapiro", "ks", "chisq", "fisher", "binom", "prop",
                                  "kruskal", "bartlett", "fligner", "oneway", "friedman", "quade", "mcnemar", "ansari"};
    bool known = false;
    for (const char* k : kinds) known = known || a[0].s == k;
    if (!known) throw dslvm::Error("test: unknown kind '" + a[0].s + "' (kinds: " + kKinds + ")");
    TestSpec s;
    s.kind = a[0].s;
    return h.hold("Test", std::move(s));
  });
  auto setter = [&b](const char* name, const char* type, std::function<void(TestSpec&, const Value&)> f) {
    b.m("Test", name, {type}, "Test", [f](void* s, A& a) { f(self<TestSpec>(s), a[0]); return Value::object("Test", s); });
  };
  setter("x", "Vec", [](TestSpec& s, const Value& v) { s.x = obj<Vec>(v); });
  setter("y", "Vec", [](TestSpec& s, const Value& v) { s.y = obj<Vec>(v); });
  setter("dist", "Dist", [](TestSpec& s, const Value& v) { s.dist = obj<Dist>(v); });
  setter("alt", "string", [](TestSpec& s, const Value& v) {
    if (v.s != "two.sided" && v.s != "less" && v.s != "greater")
      throw dslvm::Error("Test.alt: alternative must be \"two.sided\", \"less\" or \"greater\"");
    s.o.alternative = v.s;
  });
  setter("mu", "float", [](TestSpec& s, const Value& v) { s.o.mu = v.f; });
  setter("paired", "bool", [](TestSpec& s, const Value& v) { s.o.paired = v.b; });
  setter("var_equal", "bool", [](TestSpec& s, const Value& v) { s.o.var_equal = v.b; });
  setter("conf", "float", [](TestSpec& s, const Value& v) {
    if (!(v.f > 0 && v.f < 1)) throw dslvm::Error("Test.conf: confidence level must be between 0 and 1");
    s.o.conf_level = v.f;
  });
  setter("exact", "bool", [](TestSpec& s, const Value& v) { s.o.exact = v.b ? 1 : 0; });
  setter("correct", "bool", [](TestSpec& s, const Value& v) { s.o.correct = v.b; });
  setter("conf_int", "bool", [](TestSpec& s, const Value& v) { s.o.conf_int = v.b; });
  setter("ratio", "float", [](TestSpec& s, const Value& v) { s.o.ratio = v.f; });
  setter("p", "float", [](TestSpec& s, const Value& v) { s.o.p = v.f; s.p_vec = Vec{v.f}; });
  setter("p", "Vec", [](TestSpec& s, const Value& v) { s.p_vec = obj<Vec>(v); });
  setter("method", "string", [](TestSpec& s, const Value& v) { s.o.method = v.s; });
  setter("odds_ratio", "float", [](TestSpec& s, const Value& v) { s.o.or_ = v.f; });
  setter("counts", "Vec", [](TestSpec& s, const Value& v) { s.counts = obj<Vec>(v); });
  setter("probs", "Vec", [](TestSpec& s, const Value& v) { s.probs = obj<Vec>(v); });
  setter("successes", "float", [](TestSpec& s, const Value& v) { s.successes = Vec{v.f}; });
  setter("successes", "Vec", [](TestSpec& s, const Value& v) { s.successes = obj<Vec>(v); });
  setter("trials", "float", [](TestSpec& s, const Value& v) { s.trials = Vec{v.f}; });
  setter("trials", "Vec", [](TestSpec& s, const Value& v) { s.trials = obj<Vec>(v); });
  setter("groups", "Vec", [](TestSpec& s, const Value& v) { s.groups = obj<Vec>(v); });
  // Text labels become codes in sorted (byte) order, as R orders the levels of a character vector in the C locale.
  setter("groups", "StrList", [](TestSpec& s, const Value& v) { s.groups = level_codes(obj<Strs>(v)); });
  setter("blocks", "Vec", [](TestSpec& s, const Value& v) { s.blocks = obj<Vec>(v); });
  setter("blocks", "StrList", [](TestSpec& s, const Value& v) { s.blocks = level_codes(obj<Strs>(v)); });
  b.m("Test", "table", {"Vec", "int", "bool"}, "Test", [](void* s, A& a) {
    const Vec& v = obj<Vec>(a[0]);
    int64_t nr = a[1].i;
    if (nr < 1 || v.empty() || v.size() % static_cast<size_t>(nr) != 0)
      throw dslvm::Error("Test.table: " + std::to_string(v.size()) + " values do not fill " + std::to_string(nr) + " rows");
    size_t rows = static_cast<size_t>(nr), cols = v.size() / rows;
    Matrix m(rows, cols);
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < cols; ++j) m(i, j) = a[2].b ? v[i * cols + j] : v[j * rows + i];  // R fills by column by default
    self<TestSpec>(s).table = m;
    return Value::object("Test", s);
  });
  b.m("Test", "run", {}, "Htest", [&h](void* s, A&) { return h.hold("Htest", run_test(self<TestSpec>(s))); });

  // ---- Htest result ----
  b.m("Htest", "method", {}, "string", [](void* s, A&) { return Value::string(self<HtestResult>(s).method); });
  b.m("Htest", "alternative", {}, "string", [](void* s, A&) { return Value::string(self<HtestResult>(s).alternative); });
  b.m("Htest", "statistic", {}, "float", [](void* s, A&) { return Value::real(self<HtestResult>(s).statistic); });
  b.m("Htest", "statistic_name", {}, "string", [](void* s, A&) { return Value::string(self<HtestResult>(s).statistic_name); });
  b.m("Htest", "p_value", {}, "float", [](void* s, A&) { return Value::real(self<HtestResult>(s).p_value); });
  b.m("Htest", "n_params", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<HtestResult>(s).parameters.size())); });
  b.m("Htest", "param_name", {"int"}, "string", [](void* s, A& a) {
    auto& r = self<HtestResult>(s);
    return Value::string(r.parameters[static_cast<size_t>(as_index(a[0].i, r.parameters.size(), "Htest.param_name"))].first);
  });
  b.m("Htest", "param", {"int"}, "float", [](void* s, A& a) {
    auto& r = self<HtestResult>(s);
    return Value::real(r.parameters[static_cast<size_t>(as_index(a[0].i, r.parameters.size(), "Htest.param"))].second);
  });
  b.m("Htest", "param", {"string"}, "float", [](void* s, A& a) { return Value::real(param_named(self<HtestResult>(s), a[0].s)); });
  b.m("Htest", "df", {}, "float", [](void* s, A&) { return Value::real(param_named(self<HtestResult>(s), "df")); });
  b.m("Htest", "has_ci", {}, "bool", [](void* s, A&) { return Value::boolean(self<HtestResult>(s).has_conf_int); });
  auto ci = [](const HtestResult& r) {
    if (!r.has_conf_int) throw dslvm::Error("Htest: this test has no confidence interval (" + r.method + ")");
  };
  b.m("Htest", "conf_low", {}, "float", [ci](void* s, A&) { ci(self<HtestResult>(s)); return Value::real(self<HtestResult>(s).conf_low); });
  b.m("Htest", "conf_high", {}, "float", [ci](void* s, A&) { ci(self<HtestResult>(s)); return Value::real(self<HtestResult>(s).conf_high); });
  b.m("Htest", "conf_level", {}, "float", [ci](void* s, A&) { ci(self<HtestResult>(s)); return Value::real(self<HtestResult>(s).conf_level); });
  b.m("Htest", "n_estimates", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<HtestResult>(s).estimates.size())); });
  b.m("Htest", "estimate_name", {"int"}, "string", [](void* s, A& a) {
    auto& r = self<HtestResult>(s);
    return Value::string(r.estimates[static_cast<size_t>(as_index(a[0].i, r.estimates.size(), "Htest.estimate_name"))].first);
  });
  b.m("Htest", "estimate", {"int"}, "float", [](void* s, A& a) {
    auto& r = self<HtestResult>(s);
    return Value::real(r.estimates[static_cast<size_t>(as_index(a[0].i, r.estimates.size(), "Htest.estimate"))].second);
  });
  b.m("Htest", "estimate", {}, "float", [](void* s, A&) {
    auto& r = self<HtestResult>(s);
    if (r.estimates.empty()) throw dslvm::Error("Htest: this test has no estimate (" + r.method + ")");
    return Value::real(r.estimates[0].second);
  });
  b.m("Htest", "null_value", {}, "float", [](void* s, A&) { return Value::real(self<HtestResult>(s).null_value); });
  b.m("Htest", "expected", {}, "Vec", [vec](void* s, A&) { return vec(self<HtestResult>(s).expected); });
  b.m("Htest", "summary", {}, "string", [](void* s, A&) { return Value::string(summary_text(self<HtestResult>(s))); });

  // ---- FFT ----
  auto cvec = [&h](CVec v) { return h.hold("CVec", std::move(v)); };
  b.fn("fft", {"Vec"}, "CVec", [cvec](A& a) { return cvec(fft_real(obj<Vec>(a[0]), false)); });
  b.fn("cvec", {"Vec", "Vec"}, "CVec", [cvec](A& a) {
    const Vec& re = obj<Vec>(a[0]);
    const Vec& im = obj<Vec>(a[1]);
    if (re.size() != im.size()) throw dslvm::Error("cvec: real and imaginary parts differ in length");
    CVec z(re.size());
    for (size_t i = 0; i < z.size(); ++i) z[i] = {re[i], im[i]};
    return cvec(z);
  });
  b.m("CVec", "fft", {"bool"}, "CVec", [cvec](void* s, A& a) { return cvec(fft(self<CVec>(s), a[0].b)); });
  b.m("CVec", "len", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<CVec>(s).size())); });
  auto part = [&b](const char* name, double (*f)(const std::complex<double>&)) {
    b.m("CVec", name, {"int"}, "float", [f, name](void* s, A& a) {
      auto& z = self<CVec>(s);
      return Value::real(f(z[static_cast<size_t>(as_index(a[0].i, z.size(), name))]));
    });
  };
  part("re", [](const std::complex<double>& z) { return z.real(); });
  part("im", [](const std::complex<double>& z) { return z.imag(); });
  part("mod", [](const std::complex<double>& z) { return std::abs(z); });
  part("arg", [](const std::complex<double>& z) { return std::arg(z); });
  auto whole = [&b, vec](const char* name, double (*f)(const std::complex<double>&)) {
    b.m("CVec", name, {}, "Vec", [f, vec](void* s, A&) {
      Vec out;
      for (const auto& z : self<CVec>(s)) out.push_back(f(z));
      return vec(out);
    });
  };
  whole("re_vec", [](const std::complex<double>& z) { return z.real(); });
  whole("im_vec", [](const std::complex<double>& z) { return z.imag(); });
  whole("mod_vec", [](const std::complex<double>& z) { return std::abs(z); });
  b.m("CVec", "scale", {"float"}, "CVec", [cvec](void* s, A& a) {
    CVec z = self<CVec>(s);
    for (auto& c : z) c *= a[0].f;
    return cvec(z);
  });
}

}  // namespace rcli::host
