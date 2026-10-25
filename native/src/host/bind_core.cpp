// Out (the context object), Vec, Frame, math and formatting functions, keep/kept across chunks.
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <numeric>

#include "host.hpp"
#include "rcli/dd.hpp"
#include "rcli/summary.hpp"

namespace rcli::host {

namespace {

Vec& V(const Value& v) { return obj<Vec>(v); }
Strs& S(const Value& v) { return strs_of(v); }

void same_length(const Vec& a, const Vec& b, const char* what) {
  if (a.size() != b.size())
    throw dslvm::Error(std::string(what) + ": lengths differ (" + std::to_string(a.size()) + " and " + std::to_string(b.size()) + ")");
}

template <class F>
Vec zip(const Vec& a, const Vec& b, const char* what, F f) {
  same_length(a, b, what);
  Vec out(a.size());
  for (size_t i = 0; i < a.size(); ++i) out[i] = f(a[i], b[i]);
  return out;
}

template <class F>
Vec map(const Vec& a, F f) {
  Vec out(a.size());
  for (size_t i = 0; i < a.size(); ++i) out[i] = f(a[i]);
  return out;
}

double round_digits(double x, int64_t digits) {
  if (!std::isfinite(x)) return x;
  if (digits < -15 || digits > 15) throw dslvm::Error("round: digits must be between -15 and 15");
  if (digits == 0) return std::nearbyint(x);
  double s = std::pow(10.0, static_cast<double>(digits));
  double y = x * s;
  if (!std::isfinite(y)) return x;
  return std::nearbyint(y) / s;
}

Vec parse_numbers(const std::string& text) {
  Vec out;
  size_t i = 0;
  while (i < text.size()) {
    while (i < text.size() && (text[i] == ' ' || text[i] == ',' || text[i] == '\t' || text[i] == '\n' || text[i] == ';')) ++i;
    if (i >= text.size()) break;
    size_t j = i;
    while (j < text.size() && !(text[j] == ' ' || text[j] == ',' || text[j] == '\t' || text[j] == '\n' || text[j] == ';')) ++j;
    std::string tok = text.substr(i, j - i);
    if (tok == "NA" || tok == "NaN") out.push_back(kNaN);
    else if (tok == "Inf" || tok == "+Inf") out.push_back(kInf);
    else if (tok == "-Inf") out.push_back(-kInf);
    else {
      double v;
      const char* b = tok.data();
      if (*b == '+') ++b;
      if (!parse_number(b, tok.data() + tok.size(), v)) throw dslvm::Error("parse_vec: not a number: '" + tok + "'");
      out.push_back(v);
    }
    i = j;
  }
  return out;
}

std::string vec_text(const Vec& v, const std::string& spec, const std::string& sep) {
  std::string out;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) out += sep;
    out += std::isnan(v[i]) ? "NA" : format_number(spec, v[i]);
  }
  return out;
}

Strs split_names(const std::string& list) {
  Strs out;
  std::string cur;
  for (char c : list) {
    if (c == ',') { out.push_back(cur); cur.clear(); }
    else if (c != ' ' && c != '\t') cur += c;
  }
  out.push_back(cur);
  for (const auto& n : out) if (n.empty()) throw dslvm::Error("select: empty column name in '" + list + "'");
  return out;
}

}  // namespace

void bind_core(Binder& b) {
  Host& h = b.host;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };
  auto strs = [&b](Strs v) { return str_list(b.reg, std::move(v)); };
  auto frame = [&h](Frame f) { return h.hold("Frame", std::move(f)); };
  using A = Args;

  // ---- Out: the single context object ----
  b.reg.context("out", "Out");
  b.m("Out", "print", {"string"}, "void", [](void*, A& a) {
    std::cerr << a[0].s << "\n";
    return Value::none();
  });
  b.m("Out", "check", {"bool", "string"}, "bool", [&h](void*, A& a) {
    ++h.out.checks;
    if (!a[0].b) {
      ++h.out.failures;
      if (h.out.first_failure.empty()) h.out.first_failure = a[1].s;
      std::cerr << "CHECK FAILED: " << a[1].s << "\n";
    }
    return Value::boolean(a[0].b);
  });
  b.m("Out", "check_close", {"float", "float", "float", "string"}, "bool", [&h](void*, A& a) {
    ++h.out.checks;
    double actual = a[0].f, expected = a[1].f, tol = a[2].f;
    bool ok = std::isfinite(actual) && std::fabs(actual - expected) <= tol;
    if (!ok) {
      ++h.out.failures;
      char buf[200];
      std::snprintf(buf, sizeof buf, ": actual %.10g, expected %.10g, tolerance %.3g", actual, expected, tol);
      std::string msg = a[3].s + buf;
      if (h.out.first_failure.empty()) h.out.first_failure = msg;
      std::cerr << "CHECK FAILED: " << msg << "\n";
    }
    return Value::boolean(ok);
  });
  b.m("Out", "failures", {}, "int", [&h](void*, A&) { return Value::integer(h.out.failures); });
  b.m("Out", "fail", {"string"}, "void", [](void*, A& a) -> Value { throw dslvm::Error("APP_FAIL", a[0].s); });
  b.m("Out", "set", {"string", "float"}, "void", [&h](void*, A& a) { h.out.set(a[0].s, json_number(a[1].f)); return Value::none(); });
  b.m("Out", "set", {"string", "int"}, "void", [&h](void*, A& a) { h.out.set(a[0].s, std::to_string(a[1].i)); return Value::none(); });
  b.m("Out", "set", {"string", "string"}, "void", [&h](void*, A& a) { h.out.set(a[0].s, dslvm::json_escape(a[1].s)); return Value::none(); });
  b.m("Out", "set", {"string", "bool"}, "void", [&h](void*, A& a) { h.out.set(a[0].s, a[1].b ? "true" : "false"); return Value::none(); });
  b.m("Out", "set", {"string", "Vec"}, "void", [&h](void*, A& a) { h.out.set(a[0].s, json_vec(V(a[1]))); return Value::none(); });
  b.m("Out", "set", {"string", "StrList"}, "void", [&h](void*, A& a) {
    std::string j = "[";
    for (size_t i = 0; i < S(a[1]).size(); ++i) j += (i ? ", " : "") + dslvm::json_escape(S(a[1])[i]);
    h.out.set(a[0].s, j + "]");
    return Value::none();
  });

  // ---- Vec construction ----
  b.fn("vec", {}, "Vec", [vec](A&) { return vec({}); });
  b.fn("zeros", {"int"}, "Vec", [vec](A& a) {
    if (a[0].i < 0 || a[0].i > 100000000) throw dslvm::Error("zeros: n must be between 0 and 1e8");
    return vec(Vec(static_cast<size_t>(a[0].i), 0.0));
  });
  b.fn("rep", {"float", "int"}, "Vec", [vec](A& a) {
    if (a[1].i < 0 || a[1].i > 100000000) throw dslvm::Error("rep: n must be between 0 and 1e8");
    return vec(Vec(static_cast<size_t>(a[1].i), a[0].f));
  });
  b.fn("seq", {"float", "float", "float"}, "Vec", [vec](A& a) {
    double from = a[0].f, to = a[1].f, by = a[2].f;
    if (!std::isfinite(from) || !std::isfinite(to) || !std::isfinite(by) || by == 0 || (to - from) / by < 0)
      throw dslvm::Error("seq: wrong sign or zero 'by' for from=" + std::to_string(from) + " to=" + std::to_string(to));
    double n = std::floor((to - from) / by + 1e-10);
    if (n > 1e8) throw dslvm::Error("seq: too many values");
    Vec v(static_cast<size_t>(n) + 1);
    for (size_t i = 0; i < v.size(); ++i) v[i] = from + static_cast<double>(i) * by;
    return vec(std::move(v));
  });
  b.fn("seq_len", {"int"}, "Vec", [vec](A& a) {
    if (a[0].i < 0 || a[0].i > 100000000) throw dslvm::Error("seq_len: n must be between 0 and 1e8");
    Vec v(static_cast<size_t>(a[0].i));
    for (size_t i = 0; i < v.size(); ++i) v[i] = static_cast<double>(i + 1);
    return vec(std::move(v));
  });
  b.fn("linspace", {"float", "float", "int"}, "Vec", [vec](A& a) {
    if (a[2].i < 2 || a[2].i > 100000000) throw dslvm::Error("linspace: n must be between 2 and 1e8");
    size_t n = static_cast<size_t>(a[2].i);
    Vec v(n);
    for (size_t i = 0; i < n; ++i) v[i] = a[0].f + (a[1].f - a[0].f) * static_cast<double>(i) / static_cast<double>(n - 1);
    v[n - 1] = a[1].f;
    return vec(std::move(v));
  });
  b.fn("parse_vec", {"string"}, "Vec", [vec](A& a) { return vec(parse_numbers(a[0].s)); });

  // ---- Vec access and mutation ----
  b.m("Vec", "len", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<Vec>(s).size())); });
  b.m("Vec", "at", {"int"}, "float", [](void* s, A& a) {
    auto& v = self<Vec>(s);
    return Value::real(v[static_cast<size_t>(as_index(a[0].i, v.size(), "Vec.at"))]);
  });
  b.m("Vec", "set", {"int", "float"}, "Vec", [](void* s, A& a) {
    auto& v = self<Vec>(s);
    v[static_cast<size_t>(as_index(a[0].i, v.size(), "Vec.set"))] = a[1].f;
    return Value::object("Vec", s);
  });
  b.m("Vec", "push", {"float"}, "Vec", [](void* s, A& a) { self<Vec>(s).push_back(a[0].f); return Value::object("Vec", s); });
  b.m("Vec", "append", {"Vec"}, "Vec", [](void* s, A& a) {
    Vec add = V(a[0]);  // copy first: appending a vector to itself
    self<Vec>(s).insert(self<Vec>(s).end(), add.begin(), add.end());
    return Value::object("Vec", s);
  });
  b.m("Vec", "copy", {}, "Vec", [vec](void* s, A&) { return vec(self<Vec>(s)); });
  b.m("Vec", "slice", {"int", "int"}, "Vec", [vec](void* s, A& a) {
    auto& v = self<Vec>(s);
    int64_t n = static_cast<int64_t>(v.size());
    if (a[0].i < 0 || a[1].i < a[0].i || a[1].i > n)
      throw dslvm::Error("Vec.slice: need 0 <= start <= end <= len (" + std::to_string(n) + ")");
    return vec(Vec(v.begin() + a[0].i, v.begin() + a[1].i));
  });
  b.m("Vec", "head", {"int"}, "Vec", [vec](void* s, A& a) {
    auto& v = self<Vec>(s);
    size_t k = static_cast<size_t>(std::clamp<int64_t>(a[0].i, 0, static_cast<int64_t>(v.size())));
    return vec(Vec(v.begin(), v.begin() + static_cast<long>(k)));
  });
  b.m("Vec", "tail", {"int"}, "Vec", [vec](void* s, A& a) {
    auto& v = self<Vec>(s);
    size_t k = static_cast<size_t>(std::clamp<int64_t>(a[0].i, 0, static_cast<int64_t>(v.size())));
    return vec(Vec(v.end() - static_cast<long>(k), v.end()));
  });
  b.m("Vec", "rev", {}, "Vec", [vec](void* s, A&) { Vec v = self<Vec>(s); std::reverse(v.begin(), v.end()); return vec(v); });
  b.m("Vec", "sort", {}, "Vec", [vec](void* s, A&) {
    Vec v = self<Vec>(s);
    for (double x : v) if (std::isnan(x)) throw dslvm::Error("Vec.sort: NaN values (use drop_nan first)");
    std::sort(v.begin(), v.end());
    return vec(v);
  });
  b.m("Vec", "order", {}, "Vec", [vec](void* s, A&) {
    const Vec& v = self<Vec>(s);
    for (double x : v) if (std::isnan(x)) throw dslvm::Error("Vec.order: NaN values (use drop_nan first)");
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });
    Vec out(v.size());
    for (size_t i = 0; i < idx.size(); ++i) out[i] = static_cast<double>(idx[i]);
    return vec(out);
  });
  b.m("Vec", "take", {"Vec"}, "Vec", [vec](void* s, A& a) {
    const Vec& v = self<Vec>(s);
    Vec out;
    for (double d : V(a[0])) {
      if (!(d >= 0) || d != std::floor(d) || d >= static_cast<double>(v.size()))
        throw dslvm::Error("Vec.take: index " + std::to_string(d) + " is not an integer in 0.." + std::to_string(v.size() - 1));
      out.push_back(v[static_cast<size_t>(d)]);
    }
    return vec(out);
  });
  b.m("Vec", "where", {"Vec"}, "Vec", [vec](void* s, A& a) {
    const Vec& v = self<Vec>(s);
    same_length(v, V(a[0]), "Vec.where");
    Vec out;
    for (size_t i = 0; i < v.size(); ++i) if (V(a[0])[i] != 0 && !std::isnan(V(a[0])[i])) out.push_back(v[i]);
    return vec(out);
  });
  b.m("Vec", "rank", {}, "Vec", [vec](void* s, A&) { return vec(summary::rank(self<Vec>(s))); });

  // ---- Vec summaries ----
  auto red = [&b](const char* name, double (*f)(const Vec&)) {
    b.m("Vec", name, {}, "float", [f](void* s, A&) { return Value::real(f(self<Vec>(s))); });
  };
  red("sum", summary::sum);
  red("mean", summary::mean);
  red("var", summary::var);
  red("sd", summary::sd);
  red("median", summary::median);
  b.m("Vec", "mad", {}, "float", [](void* s, A&) { return Value::real(summary::mad(self<Vec>(s))); });
  b.m("Vec", "min", {}, "float", [](void* s, A&) {
    const Vec& v = self<Vec>(s);
    if (v.empty()) throw dslvm::Error("Vec.min: empty vector");
    double m = v[0];
    for (double x : v) { if (std::isnan(x)) return Value::real(kNaN); m = std::min(m, x); }
    return Value::real(m);
  });
  b.m("Vec", "max", {}, "float", [](void* s, A&) {
    const Vec& v = self<Vec>(s);
    if (v.empty()) throw dslvm::Error("Vec.max: empty vector");
    double m = v[0];
    for (double x : v) { if (std::isnan(x)) return Value::real(kNaN); m = std::max(m, x); }
    return Value::real(m);
  });
  b.m("Vec", "prod", {}, "float", [](void* s, A&) {
    double p = 1;
    for (double x : self<Vec>(s)) p *= x;
    return Value::real(p);
  });
  b.m("Vec", "which_min", {}, "int", [](void* s, A&) {
    const Vec& v = self<Vec>(s);
    int64_t k = -1;
    for (size_t i = 0; i < v.size(); ++i) if (!std::isnan(v[i]) && (k < 0 || v[i] < v[static_cast<size_t>(k)])) k = static_cast<int64_t>(i);
    if (k < 0) throw dslvm::Error("Vec.which_min: no non-NaN values");
    return Value::integer(k);
  });
  b.m("Vec", "which_max", {}, "int", [](void* s, A&) {
    const Vec& v = self<Vec>(s);
    int64_t k = -1;
    for (size_t i = 0; i < v.size(); ++i) if (!std::isnan(v[i]) && (k < 0 || v[i] > v[static_cast<size_t>(k)])) k = static_cast<int64_t>(i);
    if (k < 0) throw dslvm::Error("Vec.which_max: no non-NaN values");
    return Value::integer(k);
  });
  b.m("Vec", "count_nan", {}, "int", [](void* s, A&) {
    int64_t k = 0;
    for (double x : self<Vec>(s)) k += std::isnan(x);
    return Value::integer(k);
  });
  b.m("Vec", "drop_nan", {}, "Vec", [vec](void* s, A&) {
    Vec out;
    for (double x : self<Vec>(s)) if (!std::isnan(x)) out.push_back(x);
    return vec(out);
  });
  b.m("Vec", "quantile", {"float"}, "float", [](void* s, A& a) { return Value::real(summary::quantile(self<Vec>(s), a[0].f, 7)); });
  b.m("Vec", "quantile", {"float", "int"}, "float", [](void* s, A& a) {
    return Value::real(summary::quantile(self<Vec>(s), a[0].f, static_cast<int>(a[1].i)));
  });
  b.m("Vec", "quantiles", {"Vec", "int"}, "Vec", [vec](void* s, A& a) {
    return vec(summary::quantiles(self<Vec>(s), V(a[0]), static_cast<int>(a[1].i)));
  });
  b.m("Vec", "cor", {"Vec"}, "float", [](void* s, A& a) { return Value::real(summary::cor(self<Vec>(s), V(a[0]), "pearson")); });
  b.m("Vec", "cor", {"Vec", "string"}, "float", [](void* s, A& a) { return Value::real(summary::cor(self<Vec>(s), V(a[0]), a[1].s)); });
  b.m("Vec", "cov", {"Vec"}, "float", [](void* s, A& a) { return Value::real(summary::cov(self<Vec>(s), V(a[0]))); });
  b.m("Vec", "dot", {"Vec"}, "float", [](void* s, A& a) {
    same_length(self<Vec>(s), V(a[0]), "Vec.dot");
    DD acc;  // double-double: the same bits on every platform
    for (size_t i = 0; i < V(a[0]).size(); ++i) acc += DD(self<Vec>(s)[i]) * DD(V(a[0])[i]);
    return Value::real(static_cast<double>(acc));
  });

  // ---- Vec arithmetic (new vectors; the receiver is unchanged) ----
  b.m("Vec", "add", {"Vec"}, "Vec", [vec](void* s, A& a) { return vec(zip(self<Vec>(s), V(a[0]), "Vec.add", [](double x, double y) { return x + y; })); });
  b.m("Vec", "add", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return x + c; })); });
  b.m("Vec", "sub", {"Vec"}, "Vec", [vec](void* s, A& a) { return vec(zip(self<Vec>(s), V(a[0]), "Vec.sub", [](double x, double y) { return x - y; })); });
  b.m("Vec", "sub", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return x - c; })); });
  b.m("Vec", "rsub", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return c - x; })); });
  b.m("Vec", "mul", {"Vec"}, "Vec", [vec](void* s, A& a) { return vec(zip(self<Vec>(s), V(a[0]), "Vec.mul", [](double x, double y) { return x * y; })); });
  b.m("Vec", "mul", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return x * c; })); });
  b.m("Vec", "div", {"Vec"}, "Vec", [vec](void* s, A& a) { return vec(zip(self<Vec>(s), V(a[0]), "Vec.div", [](double x, double y) { return x / y; })); });
  b.m("Vec", "div", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return x / c; })); });
  b.m("Vec", "rdiv", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return c / x; })); });
  b.m("Vec", "neg", {}, "Vec", [vec](void* s, A&) { return vec(map(self<Vec>(s), [](double x) { return -x; })); });
  b.m("Vec", "pow", {"float"}, "Vec", [vec](void* s, A& a) { double c = a[0].f; return vec(map(self<Vec>(s), [c](double x) { return std::pow(x, c); })); });
  b.m("Vec", "exp", {}, "Vec", [vec](void* s, A&) { return vec(map(self<Vec>(s), [](double x) { return std::exp(x); })); });
  b.m("Vec", "log", {}, "Vec", [vec](void* s, A&) { return vec(map(self<Vec>(s), [](double x) { return std::log(x); })); });
  b.m("Vec", "sqrt", {}, "Vec", [vec](void* s, A&) { return vec(map(self<Vec>(s), [](double x) { return std::sqrt(x); })); });
  b.m("Vec", "abs", {}, "Vec", [vec](void* s, A&) { return vec(map(self<Vec>(s), [](double x) { return std::fabs(x); })); });
  b.m("Vec", "round", {"int"}, "Vec", [vec](void* s, A& a) { int64_t d = a[0].i; return vec(map(self<Vec>(s), [d](double x) { return round_digits(x, d); })); });
  b.m("Vec", "cumsum", {}, "Vec", [vec](void* s, A&) {
    Vec v = self<Vec>(s);
    for (size_t i = 1; i < v.size(); ++i) v[i] += v[i - 1];
    return vec(v);
  });
  b.m("Vec", "diff", {}, "Vec", [vec](void* s, A&) {
    const Vec& v = self<Vec>(s);
    Vec out;
    for (size_t i = 1; i < v.size(); ++i) out.push_back(v[i] - v[i - 1]);
    return vec(out);
  });
  b.m("Vec", "diff", {"int"}, "Vec", [vec](void* s, A& a) {
    const Vec& v = self<Vec>(s);
    if (a[0].i < 1) throw dslvm::Error("Vec.diff: lag must be at least 1");
    size_t lag = static_cast<size_t>(a[0].i);
    Vec out;
    for (size_t i = lag; i < v.size(); ++i) out.push_back(v[i] - v[i - lag]);
    return vec(out);
  });
  b.m("Vec", "standardize", {}, "Vec", [vec](void* s, A&) {
    const Vec& v = self<Vec>(s);
    double m = summary::mean(v), sd = summary::sd(v);
    return vec(map(v, [m, sd](double x) { return (x - m) / sd; }));
  });
  b.m("Vec", "format", {"string", "string"}, "string", [](void* s, A& a) { return Value::string(vec_text(self<Vec>(s), a[0].s, a[1].s)); });
  b.m("Vec", "str", {}, "string", [](void* s, A&) {
    // Shortest text that reads back to the same double, so exact integers such as 2598960 print in full.
    std::string out = "[";
    const Vec& v = self<Vec>(s);
    for (size_t i = 0; i < v.size(); ++i) out += (i ? ", " : "") + (std::isnan(v[i]) ? std::string("NA") : json_number(v[i]));
    return Value::string(out + "]");
  });
  b.op("Vec", "+", "add");
  b.op("Vec", "-", "sub");
  b.op("Vec", "*", "mul");
  b.op("Vec", "neg", "neg");
  b.op("Vec", "r+", "add");
  b.op("Vec", "r-", "rsub");
  b.op("Vec", "r*", "mul");

  // ---- Frame ----
  b.fn("read_csv", {"string"}, "Frame", [frame](A& a) { return frame(read_csv(a[0].s)); });
  b.fn("frame", {}, "Frame", [frame](A&) { return frame(Frame{}); });
  b.m("Frame", "nrow", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<Frame>(s).nrow())); });
  b.m("Frame", "ncol", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<Frame>(s).cols.size())); });
  b.m("Frame", "names", {}, "StrList", [strs](void* s, A&) { return strs(self<Frame>(s).names); });
  b.m("Frame", "has", {"string"}, "bool", [](void* s, A& a) { return Value::boolean(self<Frame>(s).find(a[0].s) >= 0); });
  b.m("Frame", "is_text", {"string"}, "bool", [](void* s, A& a) { return Value::boolean(self<Frame>(s).col(a[0].s).text); });
  b.m("Frame", "col", {"string"}, "Vec", [vec](void* s, A& a) {
    const Column& c = self<Frame>(s).col(a[0].s);
    if (c.text) throw dslvm::Error("Frame.col: column '" + a[0].s + "' is text (use Frame.text)");
    return vec(c.num);
  });
  b.m("Frame", "text", {"string"}, "StrList", [strs](void* s, A& a) {
    const Column& c = self<Frame>(s).col(a[0].s);
    if (c.text) return strs(c.str);
    Strs out;
    for (double x : c.num) out.push_back(std::isnan(x) ? "NA" : json_number(x));
    return strs(out);
  });
  b.m("Frame", "set", {"string", "Vec"}, "Frame", [](void* s, A& a) {
    Column c;
    c.num = V(a[1]);
    self<Frame>(s).set(a[0].s, std::move(c));
    return Value::object("Frame", s);
  });
  b.m("Frame", "set", {"string", "StrList"}, "Frame", [](void* s, A& a) {
    Column c;
    c.text = true;
    c.str = S(a[1]);
    self<Frame>(s).set(a[0].s, std::move(c));
    return Value::object("Frame", s);
  });
  b.m("Frame", "select", {"string"}, "Frame", [frame](void* s, A& a) {
    const Frame& f = self<Frame>(s);
    Frame out;
    for (const auto& n : split_names(a[0].s)) out.set(n, f.col(n));
    return frame(out);
  });
  b.m("Frame", "drop", {"string"}, "Frame", [frame](void* s, A& a) {  // one column, by its exact name
    const Frame& f = self<Frame>(s);
    if (f.find(a[0].s) < 0) f.col(a[0].s);  // throws the usual "no column" error
    Frame out;
    for (size_t j = 0; j < f.cols.size(); ++j)
      if (f.names[j] != a[0].s) out.set(f.names[j], f.cols[j]);
    return frame(out);
  });
  b.m("Frame", "where", {"Vec"}, "Frame", [frame](void* s, A& a) {
    const Frame& f = self<Frame>(s);
    const Vec& mask = V(a[0]);
    if (mask.size() != f.nrow()) throw dslvm::Error("Frame.where: mask length " + std::to_string(mask.size()) + " differs from nrow " + std::to_string(f.nrow()));
    Frame out;
    for (size_t j = 0; j < f.cols.size(); ++j) {
      Column c;
      c.text = f.cols[j].text;
      for (size_t r = 0; r < mask.size(); ++r) {
        if (mask[r] == 0 || std::isnan(mask[r])) continue;
        if (c.text) c.str.push_back(f.cols[j].str[r]); else c.num.push_back(f.cols[j].num[r]);
      }
      out.names.push_back(f.names[j]);
      out.cols.push_back(std::move(c));
    }
    return frame(out);
  });
  b.m("Frame", "head", {"int"}, "Frame", [frame](void* s, A& a) {
    const Frame& f = self<Frame>(s);
    size_t k = static_cast<size_t>(std::clamp<int64_t>(a[0].i, 0, static_cast<int64_t>(f.nrow())));
    Frame out = f;
    for (auto& c : out.cols) { if (c.text) c.str.resize(k); else c.num.resize(k); }
    return frame(out);
  });
  b.m("Frame", "drop_nan", {}, "Frame", [frame](void* s, A&) {
    const Frame& f = self<Frame>(s);
    Frame out;
    std::vector<bool> keep(f.nrow(), true);
    for (const auto& c : f.cols) if (!c.text) for (size_t r = 0; r < c.num.size(); ++r) if (std::isnan(c.num[r])) keep[r] = false;
    for (size_t j = 0; j < f.cols.size(); ++j) {
      Column c;
      c.text = f.cols[j].text;
      for (size_t r = 0; r < keep.size(); ++r) {
        if (!keep[r]) continue;
        if (c.text) c.str.push_back(f.cols[j].str[r]); else c.num.push_back(f.cols[j].num[r]);
      }
      out.names.push_back(f.names[j]);
      out.cols.push_back(std::move(c));
    }
    return frame(out);
  });
  b.m("Frame", "write_csv", {"string"}, "void", [](void* s, A& a) { write_csv(self<Frame>(s), a[0].s); return Value::none(); });

  // ---- keep / kept: objects that outlive a chunk (lmr --dsl a.dslb b.dslb ...) ----
  auto keep = [&h](const std::string& name, const char* type, void* p) {
    h.kept[name] = h.owner(p);
    h.kept_type[name] = type;
  };
  auto kept = [&h](const std::string& name, const char* type) {
    auto it = h.kept.find(name);
    if (it == h.kept.end()) throw dslvm::Error("nothing kept under '" + name + "'");
    if (h.kept_type[name] != type) throw dslvm::Error("'" + name + "' holds a " + h.kept_type[name] + ", not a " + type);
    return Value::object(type, it->second.get());
  };
  b.fn("keep", {"string", "Vec"}, "void", [keep](A& a) { keep(a[0].s, "Vec", a[1].obj); return Value::none(); });
  b.fn("keep", {"string", "Frame"}, "void", [keep](A& a) { keep(a[0].s, "Frame", a[1].obj); return Value::none(); });
  b.fn("kept_vec", {"string"}, "Vec", [kept](A& a) { return kept(a[0].s, "Vec"); });
  b.fn("kept_frame", {"string"}, "Frame", [kept](A& a) { return kept(a[0].s, "Frame"); });
  b.fn("is_kept", {"string"}, "bool", [&h](A& a) { return Value::boolean(h.kept.count(a[0].s) > 0); });

  // ---- math on floats (the DSL itself has sin, cos, sqrt, abs, floor, pow, min, max, clamp) ----
  auto f1 = [&b](const char* name, double (*f)(double)) {
    b.fn(name, {"float"}, "float", [f](A& a) { return Value::real(f(a[0].f)); });
  };
  f1("exp", [](double x) { return std::exp(x); });
  f1("log", [](double x) { return std::log(x); });
  f1("log10", [](double x) { return std::log10(x); });
  f1("log2", [](double x) { return std::log2(x); });
  f1("log1p", [](double x) { return std::log1p(x); });
  f1("expm1", [](double x) { return std::expm1(x); });
  f1("tan", [](double x) { return std::tan(x); });
  f1("asin", [](double x) { return std::asin(x); });
  f1("acos", [](double x) { return std::acos(x); });
  f1("atan", [](double x) { return std::atan(x); });
  f1("sinh", [](double x) { return std::sinh(x); });
  f1("cosh", [](double x) { return std::cosh(x); });
  f1("tanh", [](double x) { return std::tanh(x); });
  f1("ceil", [](double x) { return std::ceil(x); });
  f1("trunc", [](double x) { return std::trunc(x); });
  f1("round", [](double x) { return std::nearbyint(x); });
  b.fn("round", {"float", "int"}, "float", [](A& a) { return Value::real(round_digits(a[0].f, a[1].i)); });
  b.fn("atan2", {"float", "float"}, "float", [](A& a) { return Value::real(std::atan2(a[0].f, a[1].f)); });
  b.fn("is_nan", {"float"}, "bool", [](A& a) { return Value::boolean(std::isnan(a[0].f)); });
  b.fn("is_finite", {"float"}, "bool", [](A& a) { return Value::boolean(std::isfinite(a[0].f)); });
  b.fn("nan", {}, "float", [](A&) { return Value::real(kNaN); });
  b.fn("inf", {}, "float", [](A&) { return Value::real(kInf); });
  b.fn("pi", {}, "float", [](A&) { return Value::real(3.141592653589793238462643383279502884); });

  // ---- formatting ----
  b.fn("fmt", {"string", "float"}, "string", [](A& a) { return Value::string(format_number(a[0].s, a[1].f)); });
  b.fn("fmt", {"string", "int"}, "string", [](A& a) { return Value::string(format_number(a[0].s, static_cast<double>(a[1].i))); });
  b.fn("fmt", {"string", "string"}, "string", [](A& a) { return Value::string(format_string(a[0].s, a[1].s)); });
  b.fn("lpad", {"string", "int"}, "string", [](A& a) {
    std::string s = a[0].s;
    if (a[1].i > 1000) throw dslvm::Error("lpad: width too large");
    while (static_cast<int64_t>(s.size()) < a[1].i) s.insert(s.begin(), ' ');
    return Value::string(s);
  });
  b.fn("rpad", {"string", "int"}, "string", [](A& a) {
    std::string s = a[0].s;
    if (a[1].i > 1000) throw dslvm::Error("rpad: width too large");
    while (static_cast<int64_t>(s.size()) < a[1].i) s.push_back(' ');
    return Value::string(s);
  });
}

}  // namespace rcli::host
