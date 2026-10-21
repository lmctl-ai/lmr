// lm, glm (formula over a Frame) and nls (reverse communication: the script computes the model's fitted values).
#include <algorithm>
#include <cstdio>
#include <memory>
#include <set>

#include "host.hpp"
#include "rcli/dist.hpp"
#include "rcli/models.hpp"

namespace rcli::host {

namespace {

// Columns named in the formula, then complete cases only (R's default na.action = na.omit).
Data model_data(const Frame& f, const std::string& formula, Vec* weights, std::vector<size_t>* kept_rows = nullptr) {
  std::set<std::string> used;
  bool dot = false;
  size_t i = 0;
  auto ident_start = [](char c) { return std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '.'; };
  auto ident_char = [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '.'; };
  while (i < formula.size()) {
    if (ident_start(formula[i])) {
      size_t j = i;
      while (j < formula.size() && ident_char(formula[j])) ++j;
      std::string w = formula.substr(i, j - i);
      if (w == ".") dot = true;
      else if (f.find(w) >= 0) used.insert(w);
      i = j;
    } else {
      ++i;
    }
  }
  Data d;
  std::vector<const Column*> cols;
  for (size_t j = 0; j < f.names.size(); ++j) {
    const Column& c = f.cols[j];
    if (!(dot || used.count(f.names[j]))) continue;
    if (c.text) {
      if (used.count(f.names[j])) throw dslvm::Error("model: column '" + f.names[j] + "' is text; models take numeric columns");
      continue;  // "." skips text columns
    }
    d.names.push_back(f.names[j]);
    cols.push_back(&c);
  }
  if (weights && weights->size() != f.nrow())
    throw dslvm::Error("model: weights have " + std::to_string(weights->size()) + " values; the frame has " + std::to_string(f.nrow()) + " rows");
  std::vector<size_t> rows;
  for (size_t r = 0; r < f.nrow(); ++r) {
    bool ok = !(weights && std::isnan((*weights)[r]));
    for (const Column* c : cols) ok = ok && !std::isnan(c->num[r]);
    if (ok) rows.push_back(r);
  }
  for (const Column* c : cols) {
    Vec v;
    v.reserve(rows.size());
    for (size_t r : rows) v.push_back(c->num[r]);
    d.cols.push_back(std::move(v));
  }
  if (weights) {
    Vec w;
    for (size_t r : rows) w.push_back((*weights)[r]);
    *weights = std::move(w);
  }
  if (kept_rows) *kept_rows = rows;
  return d;
}

Data frame_data(const Frame& f) {
  Data d;
  for (size_t j = 0; j < f.names.size(); ++j)
    if (!f.cols[j].text) { d.names.push_back(f.names[j]); d.cols.push_back(f.cols[j].num); }
  return d;
}

size_t coef_index(const CoefTable& t, const Value& v, const std::string& what) {
  if (v.kind == Value::INT) return static_cast<size_t>(as_index(v.i, t.names.size(), what.c_str()));
  int k = t.index(v.s);
  if (k < 0) {
    std::string have;
    for (const auto& n : t.names) have += (have.empty() ? "" : ", ") + n;
    throw dslvm::Error(what + ": no coefficient '" + v.s + "' (coefficients: " + have + ")");
  }
  return static_cast<size_t>(k);
}

std::string fmt_cell(double v, int width) {
  char buf[64];
  if (std::isnan(v)) std::snprintf(buf, sizeof buf, "%*s", width, "NA");
  else if (v != 0 && (std::fabs(v) < 1e-4 || std::fabs(v) >= 1e6)) std::snprintf(buf, sizeof buf, "%*.3e", width, v);
  else std::snprintf(buf, sizeof buf, "%*.5g", width, v);
  return buf;
}

std::string coef_text(const CoefTable& t, const std::string& stat_name, const std::string& p_name) {
  size_t w = 12;
  for (const auto& n : t.names) w = std::max(w, n.size() + 1);
  std::string out = std::string(w, ' ');
  char head[200];
  std::snprintf(head, sizeof head, "%12s %12s %10s %12s\n", "Estimate", "Std. Error", stat_name.c_str(), p_name.c_str());
  out += head;
  for (size_t i = 0; i < t.names.size(); ++i) {
    std::string name = t.names[i];
    name.resize(w, ' ');
    out += name + fmt_cell(t.estimate[i], 12) + " " + fmt_cell(t.std_error[i], 12) + " " + fmt_cell(t.statistic[i], 10) + " " +
           fmt_cell(t.p_value[i], 12) + (t.aliased.size() > i && t.aliased[i] ? "  (aliased)" : "") + "\n";
  }
  return out;
}

// Registers the coefficient-table methods on a model type. `table` returns the CoefTable; `df` the residual degrees
// of freedom used for t intervals (negative for normal-theory intervals).
template <class T>
void bind_coefs(Binder& b, const std::string& type, std::function<const CoefTable&(T&)> table, std::function<double(T&)> df,
                const std::string& stat_name, const std::string& p_name) {
  using A = Args;
  Host& h = b.host;
  b.m(type, "n_coef", {}, "int", [table](void* s, A&) { return Value::integer(static_cast<int64_t>(table(self<T>(s)).names.size())); });
  b.m(type, "name", {"int"}, "string", [table, type](void* s, A& a) {
    const CoefTable& t = table(self<T>(s));
    return Value::string(t.names[static_cast<size_t>(as_index(a[0].i, t.names.size(), (type + ".name").c_str()))]);
  });
  b.m(type, "names", {}, "StrList", [table, &h](void* s, A&) { return h.hold("StrList", table(self<T>(s)).names); });
  b.m(type, "has_coef", {"string"}, "bool", [table](void* s, A& a) { return Value::boolean(table(self<T>(s)).index(a[0].s) >= 0); });
  auto field = [&b, table, type](const char* name, Vec CoefTable::*member) {
    for (const char* key : {"int", "string"}) {
      b.m(type, name, {key}, "float", [table, type, name, member](void* s, A& a) {
        const CoefTable& t = table(self<T>(s));
        return Value::real((t.*member)[coef_index(t, a[0], type + "." + name)]);
      });
    }
  };
  field("estimate", &CoefTable::estimate);
  field("se", &CoefTable::std_error);
  field("statistic", &CoefTable::statistic);
  field("p_value", &CoefTable::p_value);
  for (const char* key : {"int", "string"}) {
    b.m(type, "aliased", {key}, "bool", [table, type](void* s, A& a) {
      const CoefTable& t = table(self<T>(s));
      size_t i = coef_index(t, a[0], type + ".aliased");
      return Value::boolean(t.aliased.size() > i && t.aliased[i]);
    });
    for (int side = 0; side < 2; ++side) {
      b.m(type, side ? "conf_high" : "conf_low", {key, "float"}, "float", [table, df, type, side](void* s, A& a) {
        const CoefTable& t = table(self<T>(s));
        size_t i = coef_index(t, a[0], type + ".conf");
        double level = a[1].f;
        if (!(level > 0 && level < 1)) throw dslvm::Error(type + ".conf: level must be between 0 and 1");
        double d = df(self<T>(s));
        double q = d > 0 ? Dist::t(d).q(0.5 + level / 2) : Dist::norm(0, 1).q(0.5 + level / 2);
        return Value::real(t.estimate[i] + (side ? q : -q) * t.std_error[i]);
      });
    }
  }
  b.m(type, "estimates", {}, "Vec", [table, &h](void* s, A&) { return h.hold("Vec", table(self<T>(s)).estimate); });
  b.m(type, "coef_table", {}, "string", [table, stat_name, p_name](void* s, A&) {
    return Value::string(coef_text(table(self<T>(s)), stat_name, p_name));
  });
}

struct NlsJob {
  Vec y, start;
  std::vector<std::string> names;
  NlsOptions opt;
  std::vector<bool> has_lower, has_upper;
  Vec par, pred;
  std::vector<bool> told;
  size_t told_count = 0;
  bool pending = false, done = false;
  int evaluations = 0;
  NlsFit fit;
  Driver driver;  // last: destroyed first

  void model(const Vec& p, Vec& out) {
    par = p;
    std::fill(told.begin(), told.end(), false);
    told_count = 0;
    driver.yield();
    out = pred;
  }
  int index(const std::string& name) const {
    for (size_t i = 0; i < names.size(); ++i) if (names[i] == name) return static_cast<int>(i);
    return -1;
  }
};

}  // namespace

void bind_models(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  // ---- lm ----
  auto fit_lm = [&h](const Frame& f, const std::string& formula, const Vec* w) {
    Vec weights = w ? *w : Vec{};
    Data d = model_data(f, formula, w ? &weights : nullptr);
    Formula form = Formula::parse(formula, &d);
    return h.hold("Lm", lm(form, d, w ? &weights : nullptr));
  };
  b.fn("lm", {"Frame", "string"}, "Lm", [fit_lm](A& a) { return fit_lm(obj<Frame>(a[0]), a[1].s, nullptr); });
  b.fn("lm", {"Frame", "string", "Vec"}, "Lm", [fit_lm](A& a) { return fit_lm(obj<Frame>(a[0]), a[1].s, &obj<Vec>(a[2])); });
  bind_coefs<LmFit>(b, "Lm", [](LmFit& m) -> const CoefTable& { return m.coef; },
                    [](LmFit& m) { return static_cast<double>(m.df_residual); }, "t value", "Pr(>|t|)");
  auto lm_num = [&b](const char* name, std::function<double(const LmFit&)> f) {
    b.m("Lm", name, {}, "float", [f](void* s, A&) { return Value::real(f(self<LmFit>(s))); });
  };
  lm_num("sigma", [](const LmFit& m) { return m.sigma; });
  lm_num("r_squared", [](const LmFit& m) { return m.r_squared; });
  lm_num("adj_r_squared", [](const LmFit& m) { return m.adj_r_squared; });
  lm_num("f_statistic", [](const LmFit& m) { return m.fstatistic; });
  lm_num("f_df1", [](const LmFit& m) { return m.f_numdf; });
  lm_num("f_df2", [](const LmFit& m) { return m.f_dendf; });
  lm_num("f_p_value", [](const LmFit& m) { return m.f_pvalue; });
  lm_num("rss", [](const LmFit& m) { return m.rss; });
  lm_num("loglik", [](const LmFit& m) { return m.loglik; });
  lm_num("aic", [](const LmFit& m) { return m.aic; });
  lm_num("bic", [](const LmFit& m) { return m.bic; });
  b.m("Lm", "n", {}, "int", [](void* s, A&) { return Value::integer(self<LmFit>(s).n); });
  b.m("Lm", "rank", {}, "int", [](void* s, A&) { return Value::integer(self<LmFit>(s).rank); });
  b.m("Lm", "df_residual", {}, "int", [](void* s, A&) { return Value::integer(self<LmFit>(s).df_residual); });
  b.m("Lm", "residuals", {}, "Vec", [vec](void* s, A&) { return vec(self<LmFit>(s).residuals); });
  b.m("Lm", "fitted", {}, "Vec", [vec](void* s, A&) { return vec(self<LmFit>(s).fitted); });
  b.m("Lm", "predict", {"Frame"}, "Vec", [vec](void* s, A& a) { return vec(self<LmFit>(s).predict(frame_data(obj<Frame>(a[0])))); });
  b.m("Lm", "predict_se", {"Frame"}, "Vec", [vec](void* s, A& a) { return vec(self<LmFit>(s).predict_se(frame_data(obj<Frame>(a[0])))); });

  // ---- glm ----
  auto fit_glm = [&h](const Frame& f, const std::string& formula, const std::string& family, const std::string& link, const Vec* w) {
    GlmOptions o;
    o.family = family;
    o.link = link;
    Vec weights = w ? *w : Vec{};
    Data d = model_data(f, formula, w ? &weights : nullptr);
    if (w) o.weights = weights;
    Formula form = Formula::parse(formula, &d);
    return h.hold("Glm", glm(form, d, o));
  };
  b.fn("glm", {"Frame", "string", "string"}, "Glm", [fit_glm](A& a) { return fit_glm(obj<Frame>(a[0]), a[1].s, a[2].s, "default", nullptr); });
  b.fn("glm", {"Frame", "string", "string", "string"}, "Glm", [fit_glm](A& a) {
    return fit_glm(obj<Frame>(a[0]), a[1].s, a[2].s, a[3].s, nullptr);
  });
  b.fn("glm", {"Frame", "string", "string", "string", "Vec"}, "Glm", [fit_glm](A& a) {
    return fit_glm(obj<Frame>(a[0]), a[1].s, a[2].s, a[3].s, &obj<Vec>(a[4]));
  });
  // Prior weights (an empty Vec means none) and an offset added to the linear predictor, e.g. log(exposure).
  b.fn("glm", {"Frame", "string", "string", "string", "Vec", "Vec"}, "Glm", [&h](A& a) {
    const Frame& f = obj<Frame>(a[0]);
    const Vec& w = obj<Vec>(a[4]);
    const Vec& off = obj<Vec>(a[5]);
    if (off.size() != f.nrow())
      throw dslvm::Error("glm: offset has " + std::to_string(off.size()) + " values; the frame has " + std::to_string(f.nrow()) + " rows");
    // Carry the offset through the complete-case filter as a frame column so rows stay aligned.
    Frame g = f;
    Column oc;
    oc.num = off;
    g.set(".offset", std::move(oc));
    Vec weights = w.empty() ? Vec(f.nrow(), 1.0) : w;
    std::vector<size_t> rows;
    Data d = model_data(g, a[1].s + " + .offset", &weights, &rows);
    GlmOptions o;
    o.family = a[2].s;
    o.link = a[3].s;
    if (!w.empty()) o.weights = weights;
    o.offset = d.col(".offset");
    d.cols.erase(d.cols.begin() + (std::find(d.names.begin(), d.names.end(), ".offset") - d.names.begin()));
    d.names.erase(std::find(d.names.begin(), d.names.end(), ".offset"));
    return h.hold("Glm", glm(Formula::parse(a[1].s, &d), d, o));
  });
  bind_coefs<GlmFit>(b, "Glm", [](GlmFit& m) -> const CoefTable& { return m.coef; },
                     [](GlmFit& m) { return (m.family == "gaussian" || m.family == "Gamma") ? static_cast<double>(m.df_residual) : -1.0; },
                     "z or t", "p-value");
  auto glm_num = [&b](const char* name, std::function<double(const GlmFit&)> f) {
    b.m("Glm", name, {}, "float", [f](void* s, A&) { return Value::real(f(self<GlmFit>(s))); });
  };
  glm_num("deviance", [](const GlmFit& m) { return m.deviance; });
  glm_num("null_deviance", [](const GlmFit& m) { return m.null_deviance; });
  glm_num("aic", [](const GlmFit& m) { return m.aic; });
  glm_num("dispersion", [](const GlmFit& m) { return m.dispersion; });
  b.m("Glm", "family", {}, "string", [](void* s, A&) { return Value::string(self<GlmFit>(s).family); });
  b.m("Glm", "link", {}, "string", [](void* s, A&) { return Value::string(self<GlmFit>(s).link); });
  b.m("Glm", "n", {}, "int", [](void* s, A&) { return Value::integer(self<GlmFit>(s).n); });
  b.m("Glm", "df_residual", {}, "int", [](void* s, A&) { return Value::integer(self<GlmFit>(s).df_residual); });
  b.m("Glm", "df_null", {}, "int", [](void* s, A&) { return Value::integer(self<GlmFit>(s).df_null); });
  b.m("Glm", "iter", {}, "int", [](void* s, A&) { return Value::integer(self<GlmFit>(s).iter); });
  b.m("Glm", "converged", {}, "bool", [](void* s, A&) { return Value::boolean(self<GlmFit>(s).converged); });
  b.m("Glm", "fitted", {}, "Vec", [vec](void* s, A&) { return vec(self<GlmFit>(s).fitted); });
  b.m("Glm", "linear_predictors", {}, "Vec", [vec](void* s, A&) { return vec(self<GlmFit>(s).linear_predictors); });
  b.m("Glm", "residuals", {"string"}, "Vec", [vec](void* s, A& a) {
    const GlmFit& m = self<GlmFit>(s);
    if (a[0].s == "deviance") return vec(m.residuals_deviance);
    if (a[0].s == "pearson") return vec(m.residuals_pearson);
    if (a[0].s == "working") return vec(m.residuals_working);
    throw dslvm::Error("Glm.residuals: type must be deviance, pearson or working");
  });
  b.m("Glm", "predict", {"Frame", "string"}, "Vec", [vec](void* s, A& a) {
    if (a[1].s != "link" && a[1].s != "response") throw dslvm::Error("Glm.predict: type must be \"link\" or \"response\"");
    return vec(self<GlmFit>(s).predict(frame_data(obj<Frame>(a[0])), a[1].s));
  });

  // ---- nls (reverse communication) ----
  b.fn("nls", {"Vec"}, "Nls", [&h](A& a) {
    const Vec& y = obj<Vec>(a[0]);
    if (y.empty()) throw dslvm::Error("nls: the response is empty");
    for (double v : y) if (!std::isfinite(v)) throw dslvm::Error("nls: the response has non-finite values");
    auto [value, j] = h.emplace<NlsJob>("Nls");
    j->y = y;
    j->pred.assign(y.size(), kNaN);
    j->told.assign(y.size(), false);
    j->driver.set_algorithm([j = j] {
      if (j->names.empty()) throw rcli::Error("nls: declare parameters with .param(name, start) before the first next()");
      bool any_bound = false;
      for (size_t i = 0; i < j->names.size(); ++i) any_bound = any_bound || j->has_lower[i] || j->has_upper[i];
      if (any_bound && j->opt.algorithm != "port") throw rcli::Error("nls: bounds need .algorithm(\"port\")");
      if (!any_bound) { j->opt.lower.clear(); j->opt.upper.clear(); }
      j->fit = nls([j](const Vec& p, Vec& out) { j->model(p, out); }, j->y, j->names, j->start, j->opt);
    });
    return value;
  });
  auto configurable = [](NlsJob& j) {
    if (j.driver.started()) throw dslvm::Error("Nls: options and parameters must be set before the first next()");
  };
  b.m("Nls", "param", {"string", "float"}, "Nls", [configurable](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    configurable(j);
    if (j.index(a[0].s) >= 0) throw dslvm::Error("Nls.param: '" + a[0].s + "' is declared twice");
    if (!std::isfinite(a[1].f)) throw dslvm::Error("Nls.param: start value must be finite");
    j.names.push_back(a[0].s);
    j.start.push_back(a[1].f);
    j.has_lower.push_back(false);
    j.has_upper.push_back(false);
    j.opt.lower.push_back(-kInf);
    j.opt.upper.push_back(kInf);
    return Value::object("Nls", s);
  });
  for (int up = 0; up < 2; ++up) {
    b.m("Nls", up ? "upper" : "lower", {"string", "float"}, "Nls", [configurable, up](void* s, A& a) {
      auto& j = self<NlsJob>(s);
      configurable(j);
      int k = j.index(a[0].s);
      if (k < 0) throw dslvm::Error("Nls bound: declare '" + a[0].s + "' with .param first");
      (up ? j.has_upper : j.has_lower)[static_cast<size_t>(k)] = true;
      (up ? j.opt.upper : j.opt.lower)[static_cast<size_t>(k)] = a[1].f;
      return Value::object("Nls", s);
    });
  }
  b.m("Nls", "algorithm", {"string"}, "Nls", [configurable](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    configurable(j);
    if (a[0].s != "default" && a[0].s != "port") throw dslvm::Error("Nls.algorithm: must be \"default\" or \"port\"");
    j.opt.algorithm = a[0].s;
    return Value::object("Nls", s);
  });
  b.m("Nls", "maxiter", {"int"}, "Nls", [configurable](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    configurable(j);
    if (a[0].i < 1 || a[0].i > 10000) throw dslvm::Error("Nls.maxiter: must be between 1 and 10000");
    j.opt.maxiter = static_cast<int>(a[0].i);
    return Value::object("Nls", s);
  });
  b.m("Nls", "tol", {"float"}, "Nls", [configurable](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    configurable(j);
    j.opt.tol = a[0].f;
    return Value::object("Nls", s);
  });
  b.m("Nls", "next", {}, "bool", [](void* s, A&) {
    auto& j = self<NlsJob>(s);
    if (j.done) return Value::boolean(false);
    if (j.pending && j.told_count != j.y.size())
      throw dslvm::Error("Nls.next: tell(i, value) was called for " + std::to_string(j.told_count) + " of " + std::to_string(j.y.size()) +
                         " observations");
    bool more = j.driver.resume();
    j.pending = more;
    if (more) ++j.evaluations; else j.done = true;
    return Value::boolean(more);
  });
  auto current = [](NlsJob& j) -> const Vec& {
    if (!j.pending) throw dslvm::Error("Nls: no current parameter vector (call next() and check that it returned true)");
    return j.par;
  };
  b.m("Nls", "par", {"int"}, "float", [current](void* s, A& a) {
    const Vec& p = current(self<NlsJob>(s));
    return Value::real(p[static_cast<size_t>(as_index(a[0].i, p.size(), "Nls.par"))]);
  });
  b.m("Nls", "par", {"string"}, "float", [current](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    const Vec& p = current(j);
    int k = j.index(a[0].s);
    if (k < 0) throw dslvm::Error("Nls.par: no parameter '" + a[0].s + "'");
    return Value::real(p[static_cast<size_t>(k)]);
  });
  b.m("Nls", "n", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<NlsJob>(s).y.size())); });
  b.m("Nls", "tell", {"int", "float"}, "void", [current](void* s, A& a) {
    auto& j = self<NlsJob>(s);
    current(j);
    size_t i = static_cast<size_t>(as_index(a[0].i, j.y.size(), "Nls.tell"));
    if (j.told[i]) throw dslvm::Error("Nls.tell: observation " + std::to_string(i) + " already has a value");
    j.told[i] = true;
    ++j.told_count;
    j.pred[i] = a[1].f;
    return Value::none();
  });
  auto fit = [](void* s) -> NlsFit& {
    auto& j = self<NlsJob>(s);
    if (!j.done) throw dslvm::Error("Nls: results are available after next() returns false");
    return j.fit;
  };
  bind_coefs<NlsJob>(b, "Nls", [fit](NlsJob& j) -> const CoefTable& { return fit(&j).coef; },
                     [fit](NlsJob& j) { return static_cast<double>(fit(&j).df_residual); }, "t value", "Pr(>|t|)");
  b.m("Nls", "sigma", {}, "float", [fit](void* s, A&) { return Value::real(fit(s).sigma); });
  b.m("Nls", "rss", {}, "float", [fit](void* s, A&) { return Value::real(fit(s).rss); });
  b.m("Nls", "df_residual", {}, "int", [fit](void* s, A&) { return Value::integer(fit(s).df_residual); });
  b.m("Nls", "iterations", {}, "int", [fit](void* s, A&) { return Value::integer(fit(s).iterations); });
  b.m("Nls", "converged", {}, "bool", [fit](void* s, A&) { return Value::boolean(fit(s).converged); });
  b.m("Nls", "message", {}, "string", [fit](void* s, A&) { return Value::string(fit(s).message); });
  b.m("Nls", "conv_criterion", {}, "float", [fit](void* s, A&) { return Value::real(fit(s).conv_criterion); });
  b.m("Nls", "fitted", {}, "Vec", [fit, vec](void* s, A&) { return vec(fit(s).fitted); });
  b.m("Nls", "residuals", {}, "Vec", [fit, vec](void* s, A&) { return vec(fit(s).residuals); });
  b.m("Nls", "evaluations", {}, "int", [](void* s, A&) { return Value::integer(self<NlsJob>(s).evaluations); });
}

}  // namespace rcli::host
