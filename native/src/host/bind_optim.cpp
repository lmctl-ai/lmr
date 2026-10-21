// Reverse-communication solvers: the script computes the objective itself.
//
//   let o = minimize(parse_vec("-1.2 1"), "BFGS");
//   while (o.next()) {
//     let a = o.x(0); let b = o.x(1);
//     o.tell(100.0 * pow(b - a * a, 2.0) + pow(1.0 - a, 2.0));
//   }
//   out.print(fmt("%.6f", o.par(0)));
//
// next() runs the native algorithm until it needs f at a new point (true) or has finished (false). Each next() that
// returned true must be answered by exactly one tell(value) before the following next().
#include <memory>

#include "host.hpp"
#include "rcli/optim.hpp"
#include "rcli/rng.hpp"

namespace rcli::host {

namespace {

// Shared protocol state. `driver` is declared last so it is destroyed first: the worker thread is cancelled and
// joined before the fields it touches go away.
template <class Point, class Result>
struct Solver {
  std::string what;
  Point x{};
  double told = kNaN;
  bool pending = false, answered = false, done = false;
  int evaluations = 0;
  Result result{};
  Driver driver;

  // The callback the native algorithm calls (on the worker thread).
  double evaluate(const Point& at) {
    x = at;
    answered = false;
    driver.yield();
    return told;
  }
  bool next() {
    if (done) return false;
    if (pending && !answered) throw dslvm::Error(what + ".next: tell(value) was not called for the current point");
    bool more = driver.resume();
    pending = more;
    if (more) ++evaluations; else done = true;
    return more;
  }
  void tell(double v) {
    if (!pending) throw dslvm::Error(what + ".tell: no point is waiting (call next() first)");
    if (answered) throw dslvm::Error(what + ".tell: the current point already has a value");
    told = v;
    answered = true;
  }
  const Point& point() const {
    if (!pending) throw dslvm::Error(what + ": no current point (call next() and check that it returned true)");
    return x;
  }
  const Result& res() const {
    if (!done) throw dslvm::Error(what + ": results are available after next() returns false");
    return result;
  }
  void configurable() const {
    if (driver.started()) throw dslvm::Error(what + ": options must be set before the first next()");
  }
};

struct Minimizer : Solver<Vec, OptimResult> {
  Vec start;
  OptimOptions opt;
  std::unique_ptr<Rng> rng;
};
struct Optimizer1 : Solver<double, OptimizeResult> {
  double lower = 0, upper = 0, tol = kEps25;
  bool maximum = false;
};
struct RootFinder : Solver<double, RootResult> {
  double lower = 0, upper = 0;
  RootOptions opt;
};
struct Integral : Solver<double, IntegrateResult> {
  double lower = 0, upper = 0;
  IntegrateOptions opt;
  bool abs_tol_set = false;  // R: abs.tol defaults to rel.tol
};

}  // namespace

void bind_optim(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  // ---- minimize: optim() over a vector ----
  b.fn("minimize", {"Vec", "string"}, "Minimizer", [&h](A& a) {
    static const char* methods[] = {"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"};
    bool known = false;
    for (const char* m : methods) known = known || a[1].s == m;
    if (!known) throw dslvm::Error("minimize: method must be Nelder-Mead, BFGS, CG, L-BFGS-B or SANN");
    if (obj<Vec>(a[0]).empty()) throw dslvm::Error("minimize: start must have at least one value");
    auto [value, m] = h.emplace<Minimizer>("Minimizer");
    m->what = "Minimizer";
    m->start = obj<Vec>(a[0]);
    m->opt.method = a[1].s;
    m->driver.set_algorithm([m = m] {
      if (m->opt.method == "SANN" && !m->rng) throw rcli::Error("minimize: SANN needs .seed(n) before the first next()");
      m->opt.rng = m->rng.get();
      m->result = optim([m](const Vec& p) { return m->evaluate(p); }, m->start, m->opt);
    });
    return value;
  });
  auto mset = [&b](const char* name, const char* type, std::function<void(Minimizer&, const Value&)> f) {
    b.m("Minimizer", name, {type}, "Minimizer", [f](void* s, A& a) {
      auto& m = self<Minimizer>(s);
      m.configurable();
      f(m, a[0]);
      return Value::object("Minimizer", s);
    });
  };
  mset("lower", "Vec", [](Minimizer& m, const Value& v) { m.opt.lower = obj<Vec>(v); });
  mset("upper", "Vec", [](Minimizer& m, const Value& v) { m.opt.upper = obj<Vec>(v); });
  mset("maxit", "int", [](Minimizer& m, const Value& v) {
    if (v.i < 1 || v.i > 100000) throw dslvm::Error("Minimizer.maxit: must be between 1 and 100000");
    m.opt.maxit = static_cast<int>(v.i);
  });
  mset("reltol", "float", [](Minimizer& m, const Value& v) { m.opt.reltol = v.f; });
  mset("ndeps", "float", [](Minimizer& m, const Value& v) { m.opt.ndeps = v.f; });
  mset("hessian", "bool", [](Minimizer& m, const Value& v) { m.opt.hessian = v.b; });
  mset("seed", "int", [](Minimizer& m, const Value& v) { m.rng = std::make_unique<Rng>(static_cast<uint64_t>(v.i)); });
  b.m("Minimizer", "next", {}, "bool", [](void* s, A&) { return Value::boolean(self<Minimizer>(s).next()); });
  b.m("Minimizer", "x", {"int"}, "float", [](void* s, A& a) {
    const Vec& p = self<Minimizer>(s).point();
    return Value::real(p[static_cast<size_t>(as_index(a[0].i, p.size(), "Minimizer.x"))]);
  });
  b.m("Minimizer", "point", {}, "Vec", [vec](void* s, A&) { return vec(self<Minimizer>(s).point()); });
  b.m("Minimizer", "tell", {"float"}, "void", [](void* s, A& a) { self<Minimizer>(s).tell(a[0].f); return Value::none(); });
  b.m("Minimizer", "par", {"int"}, "float", [](void* s, A& a) {
    const Vec& p = self<Minimizer>(s).res().par;
    return Value::real(p[static_cast<size_t>(as_index(a[0].i, p.size(), "Minimizer.par"))]);
  });
  b.m("Minimizer", "pars", {}, "Vec", [vec](void* s, A&) { return vec(self<Minimizer>(s).res().par); });
  b.m("Minimizer", "value", {}, "float", [](void* s, A&) { return Value::real(self<Minimizer>(s).res().value); });
  b.m("Minimizer", "fn_count", {}, "int", [](void* s, A&) { return Value::integer(self<Minimizer>(s).res().fn_count); });
  b.m("Minimizer", "gr_count", {}, "int", [](void* s, A&) { return Value::integer(self<Minimizer>(s).res().gr_count); });
  b.m("Minimizer", "evaluations", {}, "int", [](void* s, A&) { return Value::integer(self<Minimizer>(s).evaluations); });
  b.m("Minimizer", "convergence", {}, "int", [](void* s, A&) { return Value::integer(self<Minimizer>(s).res().convergence); });
  b.m("Minimizer", "message", {}, "string", [](void* s, A&) { return Value::string(self<Minimizer>(s).res().message); });
  b.m("Minimizer", "hessian", {"int", "int"}, "float", [](void* s, A& a) {
    const auto& r = self<Minimizer>(s).res();
    if (r.hessian.empty()) throw dslvm::Error("Minimizer.hessian: call .hessian(true) before the first next()");
    size_t n = r.par.size();
    return Value::real(r.hessian[static_cast<size_t>(as_index(a[0].i, n, "Minimizer.hessian")) * n +
                                 static_cast<size_t>(as_index(a[1].i, n, "Minimizer.hessian"))]);
  });

  // ---- optimize: one-dimensional minimum or maximum on an interval ----
  auto make_opt1 = [&h](double lower, double upper, bool maximum) {
    if (!(std::isfinite(lower) && std::isfinite(upper) && lower < upper))
      throw dslvm::Error("optimize: need finite lower < upper");
    auto [value, o] = h.emplace<Optimizer1>("Optimizer1");
    o->what = "Optimizer1";
    o->lower = lower;
    o->upper = upper;
    o->maximum = maximum;
    o->driver.set_algorithm([o = o] {
      o->result = optimize([o](double x) { return o->evaluate(x); }, o->lower, o->upper, o->maximum, o->tol);
    });
    return value;
  };
  b.fn("optimize", {"float", "float"}, "Optimizer1", [make_opt1](A& a) { return make_opt1(a[0].f, a[1].f, false); });
  b.fn("optimize", {"float", "float", "bool"}, "Optimizer1", [make_opt1](A& a) { return make_opt1(a[0].f, a[1].f, a[2].b); });
  b.m("Optimizer1", "tol", {"float"}, "Optimizer1", [](void* s, A& a) {
    auto& o = self<Optimizer1>(s);
    o.configurable();
    if (!(a[0].f > 0)) throw dslvm::Error("Optimizer1.tol: must be positive");
    o.tol = a[0].f;
    return Value::object("Optimizer1", s);
  });
  b.m("Optimizer1", "next", {}, "bool", [](void* s, A&) { return Value::boolean(self<Optimizer1>(s).next()); });
  b.m("Optimizer1", "x", {}, "float", [](void* s, A&) { return Value::real(self<Optimizer1>(s).point()); });
  b.m("Optimizer1", "tell", {"float"}, "void", [](void* s, A& a) { self<Optimizer1>(s).tell(a[0].f); return Value::none(); });
  b.m("Optimizer1", "argmin", {}, "float", [](void* s, A&) { return Value::real(self<Optimizer1>(s).res().minimum); });
  b.m("Optimizer1", "objective", {}, "float", [](void* s, A&) { return Value::real(self<Optimizer1>(s).res().objective); });
  b.m("Optimizer1", "evaluations", {}, "int", [](void* s, A&) { return Value::integer(self<Optimizer1>(s).evaluations); });

  // ---- uniroot ----
  b.fn("uniroot", {"float", "float"}, "RootFinder", [&h](A& a) {
    if (!(std::isfinite(a[0].f) && std::isfinite(a[1].f) && a[0].f < a[1].f))
      throw dslvm::Error("uniroot: need finite lower < upper");
    auto [value, r] = h.emplace<RootFinder>("RootFinder");
    r->what = "RootFinder";
    r->lower = a[0].f;
    r->upper = a[1].f;
    r->driver.set_algorithm([r = r] {
      r->result = uniroot([r](double x) { return r->evaluate(x); }, r->lower, r->upper, r->opt);
    });
    return value;
  });
  auto rset = [&b](const char* name, const char* type, std::function<void(RootFinder&, const Value&)> f) {
    b.m("RootFinder", name, {type}, "RootFinder", [f](void* s, A& a) {
      auto& r = self<RootFinder>(s);
      r.configurable();
      f(r, a[0]);
      return Value::object("RootFinder", s);
    });
  };
  rset("tol", "float", [](RootFinder& r, const Value& v) { r.opt.tol = v.f; });
  rset("maxiter", "int", [](RootFinder& r, const Value& v) {
    if (v.i < 1 || v.i > 100000) throw dslvm::Error("RootFinder.maxiter: must be between 1 and 100000");
    r.opt.maxiter = static_cast<int>(v.i);
  });
  rset("extend", "string", [](RootFinder& r, const Value& v) {
    if (v.s != "no" && v.s != "yes" && v.s != "upX" && v.s != "downX")
      throw dslvm::Error("RootFinder.extend: must be \"no\", \"yes\", \"upX\" or \"downX\"");
    r.opt.extend_int = v.s;
  });
  b.m("RootFinder", "next", {}, "bool", [](void* s, A&) { return Value::boolean(self<RootFinder>(s).next()); });
  b.m("RootFinder", "x", {}, "float", [](void* s, A&) { return Value::real(self<RootFinder>(s).point()); });
  b.m("RootFinder", "tell", {"float"}, "void", [](void* s, A& a) { self<RootFinder>(s).tell(a[0].f); return Value::none(); });
  b.m("RootFinder", "root", {}, "float", [](void* s, A&) { return Value::real(self<RootFinder>(s).res().root); });
  b.m("RootFinder", "f_root", {}, "float", [](void* s, A&) { return Value::real(self<RootFinder>(s).res().f_root); });
  b.m("RootFinder", "iter", {}, "int", [](void* s, A&) { return Value::integer(self<RootFinder>(s).res().iter); });
  b.m("RootFinder", "init_it", {}, "int", [](void* s, A&) { return Value::integer(self<RootFinder>(s).res().init_it); });
  b.m("RootFinder", "estim_prec", {}, "float", [](void* s, A&) { return Value::real(self<RootFinder>(s).res().estim_prec); });
  b.m("RootFinder", "converged", {}, "bool", [](void* s, A&) { return Value::boolean(self<RootFinder>(s).res().converged); });
  b.m("RootFinder", "evaluations", {}, "int", [](void* s, A&) { return Value::integer(self<RootFinder>(s).evaluations); });

  // ---- integrate ----
  b.fn("integrate", {"float", "float"}, "Integral", [&h](A& a) {
    if (std::isnan(a[0].f) || std::isnan(a[1].f)) throw dslvm::Error("integrate: bounds must not be NaN");
    auto [value, q] = h.emplace<Integral>("Integral");
    q->what = "Integral";
    q->lower = a[0].f;
    q->upper = a[1].f;
    q->driver.set_algorithm([q = q] {
      // The DSL has no try/catch, so a non-finite integrand value is reported (ier 6, message), not thrown.
      if (!q->abs_tol_set) q->opt.abs_tol = q->opt.rel_tol;
      try {
        q->result = integrate([q](double x) { return q->evaluate(x); }, q->lower, q->upper, q->opt);
      } catch (const rcli::Error& e) {
        q->result = IntegrateResult{};
        q->result.ier = 6;
        q->result.message = e.what();
        q->result.evaluations = q->evaluations;
      }
    });
    return value;
  });
  auto qset = [&b](const char* name, const char* type, std::function<void(Integral&, const Value&)> f) {
    b.m("Integral", name, {type}, "Integral", [f](void* s, A& a) {
      auto& q = self<Integral>(s);
      q.configurable();
      f(q, a[0]);
      return Value::object("Integral", s);
    });
  };
  qset("subdivisions", "int", [](Integral& q, const Value& v) {
    if (v.i < 1 || v.i > 100000) throw dslvm::Error("Integral.subdivisions: must be between 1 and 100000");
    q.opt.subdivisions = static_cast<int>(v.i);
  });
  qset("rel_tol", "float", [](Integral& q, const Value& v) { q.opt.rel_tol = v.f; });
  qset("abs_tol", "float", [](Integral& q, const Value& v) {
    q.opt.abs_tol = v.f;
    q.abs_tol_set = true;
  });
  b.m("Integral", "next", {}, "bool", [](void* s, A&) { return Value::boolean(self<Integral>(s).next()); });
  b.m("Integral", "x", {}, "float", [](void* s, A&) { return Value::real(self<Integral>(s).point()); });
  b.m("Integral", "tell", {"float"}, "void", [](void* s, A& a) { self<Integral>(s).tell(a[0].f); return Value::none(); });
  b.m("Integral", "value", {}, "float", [](void* s, A&) { return Value::real(self<Integral>(s).res().value); });
  b.m("Integral", "abs_error", {}, "float", [](void* s, A&) { return Value::real(self<Integral>(s).res().abs_error); });
  b.m("Integral", "subdivisions_used", {}, "int", [](void* s, A&) { return Value::integer(self<Integral>(s).res().subdivisions); });
  b.m("Integral", "evaluations", {}, "int", [](void* s, A&) { return Value::integer(self<Integral>(s).evaluations); });
  b.m("Integral", "ier", {}, "int", [](void* s, A&) { return Value::integer(self<Integral>(s).res().ier); });
  b.m("Integral", "ok", {}, "bool", [](void* s, A&) { return Value::boolean(self<Integral>(s).res().ier == 0); });
  b.m("Integral", "message", {}, "string", [](void* s, A&) { return Value::string(self<Integral>(s).res().message); });
}

}  // namespace rcli::host
