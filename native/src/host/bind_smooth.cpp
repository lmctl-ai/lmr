// Interpolating splines, piecewise approximation, smoothing splines and loess.
#include "host.hpp"
#include "rcli/interp.hpp"
#include "rcli/smooth.hpp"

namespace rcli::host {

void bind_smooth(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };
  auto eval = [vec](const Vec& xs, const std::function<double(double)>& f) {
    Vec out(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) out[i] = f(xs[i]);
    return vec(std::move(out));
  };

  // ---- spline_fun: R splinefun ----
  b.fn("spline_fun", {"Vec", "Vec"}, "Spline", [&h](A& a) { return h.hold("Spline", CubicSpline(obj<Vec>(a[0]), obj<Vec>(a[1]), "fmm")); });
  b.fn("spline_fun", {"Vec", "Vec", "string"}, "Spline", [&h](A& a) {
    return h.hold("Spline", CubicSpline(obj<Vec>(a[0]), obj<Vec>(a[1]), a[2].s));
  });
  b.m("Spline", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<CubicSpline>(s)(a[0].f)); });
  b.m("Spline", "deriv", {"float", "int"}, "float", [](void* s, A& a) {
    if (a[1].i < 1 || a[1].i > 3) throw dslvm::Error("Spline.deriv: order must be 1, 2 or 3");
    return Value::real(self<CubicSpline>(s).deriv(a[0].f, static_cast<int>(a[1].i)));
  });
  b.m("Spline", "eval", {"Vec"}, "Vec", [eval](void* s, A& a) {
    auto& sp = self<CubicSpline>(s);
    return eval(obj<Vec>(a[0]), [&sp](double x) { return sp(x); });
  });
  b.m("Spline", "method", {}, "string", [](void* s, A&) { return Value::string(self<CubicSpline>(s).method()); });

  // ---- approx_fun: R approxfun ----
  b.fn("approx_fun", {"Vec", "Vec"}, "ApproxFun", [&h](A& a) { return h.hold("ApproxFun", Approx(obj<Vec>(a[0]), obj<Vec>(a[1]))); });
  b.fn("approx_fun", {"Vec", "Vec", "string", "int"}, "ApproxFun", [&h](A& a) {
    if (a[3].i != 1 && a[3].i != 2) throw dslvm::Error("approx_fun: rule must be 1 or 2");
    int rule = static_cast<int>(a[3].i);
    return h.hold("ApproxFun", Approx(obj<Vec>(a[0]), obj<Vec>(a[1]), a[2].s, rule, rule));
  });
  b.fn("approx_fun", {"Vec", "Vec", "string", "int", "float"}, "ApproxFun", [&h](A& a) {
    if (a[3].i != 1 && a[3].i != 2) throw dslvm::Error("approx_fun: rule must be 1 or 2");
    int rule = static_cast<int>(a[3].i);
    return h.hold("ApproxFun", Approx(obj<Vec>(a[0]), obj<Vec>(a[1]), a[2].s, rule, rule, a[4].f));
  });
  b.m("ApproxFun", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Approx>(s)(a[0].f)); });
  b.m("ApproxFun", "eval", {"Vec"}, "Vec", [eval](void* s, A& a) {
    auto& f = self<Approx>(s);
    return eval(obj<Vec>(a[0]), [&f](double x) { return f(x); });
  });

  // ---- smooth_spline ----
  auto ss = [&h](const Vec& x, const Vec& y, SmoothSplineOptions o) { return h.hold("SmoothSpline", smooth_spline(x, y, o)); };
  b.fn("smooth_spline", {"Vec", "Vec"}, "SmoothSpline", [ss](A& a) { return ss(obj<Vec>(a[0]), obj<Vec>(a[1]), {}); });
  b.fn("smooth_spline", {"Vec", "Vec", "bool"}, "SmoothSpline", [ss](A& a) {
    SmoothSplineOptions o;
    o.cv = a[2].b;
    return ss(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.fn("smooth_spline_df", {"Vec", "Vec", "float"}, "SmoothSpline", [ss](A& a) {
    SmoothSplineOptions o;
    o.df = a[2].f;
    return ss(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.fn("smooth_spline_lambda", {"Vec", "Vec", "float"}, "SmoothSpline", [ss](A& a) {
    SmoothSplineOptions o;
    o.lambda = a[2].f;
    return ss(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.m("SmoothSpline", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<SmoothSplineFit>(s).predict(a[0].f)); });
  b.m("SmoothSpline", "deriv", {"float", "int"}, "float", [](void* s, A& a) {
    if (a[1].i < 1 || a[1].i > 3) throw dslvm::Error("SmoothSpline.deriv: order must be 1, 2 or 3");
    return Value::real(self<SmoothSplineFit>(s).predict(a[0].f, static_cast<int>(a[1].i)));
  });
  b.m("SmoothSpline", "eval", {"Vec"}, "Vec", [eval](void* s, A& a) {
    auto& f = self<SmoothSplineFit>(s);
    return eval(obj<Vec>(a[0]), [&f](double x) { return f.predict(x); });
  });
  b.m("SmoothSpline", "x", {}, "Vec", [vec](void* s, A&) { return vec(self<SmoothSplineFit>(s).x); });
  b.m("SmoothSpline", "y", {}, "Vec", [vec](void* s, A&) { return vec(self<SmoothSplineFit>(s).y); });
  b.m("SmoothSpline", "leverages", {}, "Vec", [vec](void* s, A&) { return vec(self<SmoothSplineFit>(s).lev); });
  b.m("SmoothSpline", "df", {}, "float", [](void* s, A&) { return Value::real(self<SmoothSplineFit>(s).df); });
  b.m("SmoothSpline", "lambda", {}, "float", [](void* s, A&) { return Value::real(self<SmoothSplineFit>(s).lambda); });
  b.m("SmoothSpline", "cv_crit", {}, "float", [](void* s, A&) { return Value::real(self<SmoothSplineFit>(s).cv_crit); });
  b.m("SmoothSpline", "pen_crit", {}, "float", [](void* s, A&) { return Value::real(self<SmoothSplineFit>(s).pen_crit); });

  // ---- loess ----
  auto lo = [&h](const Vec& x, const Vec& y, LoessOptions o) { return h.hold("Loess", loess(x, y, o)); };
  b.fn("loess", {"Vec", "Vec"}, "Loess", [lo](A& a) { return lo(obj<Vec>(a[0]), obj<Vec>(a[1]), {}); });
  b.fn("loess", {"Vec", "Vec", "float"}, "Loess", [lo](A& a) {
    LoessOptions o;
    o.span = a[2].f;
    return lo(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.fn("loess", {"Vec", "Vec", "float", "int"}, "Loess", [lo](A& a) {
    LoessOptions o;
    o.span = a[2].f;
    o.degree = static_cast<int>(a[3].i);
    return lo(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.fn("loess", {"Vec", "Vec", "float", "int", "string"}, "Loess", [lo](A& a) {
    LoessOptions o;
    o.span = a[2].f;
    o.degree = static_cast<int>(a[3].i);
    o.family = a[4].s;
    return lo(obj<Vec>(a[0]), obj<Vec>(a[1]), o);
  });
  b.m("Loess", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<LoessFit>(s).predict(a[0].f)); });
  b.m("Loess", "eval", {"Vec"}, "Vec", [eval](void* s, A& a) {
    auto& f = self<LoessFit>(s);
    return eval(obj<Vec>(a[0]), [&f](double x) { return f.predict(x); });
  });
  b.m("Loess", "fitted", {}, "Vec", [vec](void* s, A&) { return vec(self<LoessFit>(s).fitted); });
  b.m("Loess", "residuals", {}, "Vec", [vec](void* s, A&) { return vec(self<LoessFit>(s).residuals); });
  b.m("Loess", "robustness_weights", {}, "Vec", [vec](void* s, A&) { return vec(self<LoessFit>(s).robustness_weights); });
  b.m("Loess", "enp", {}, "float", [](void* s, A&) { return Value::real(self<LoessFit>(s).enp); });
  b.m("Loess", "s", {}, "float", [](void* s, A&) { return Value::real(self<LoessFit>(s).s); });
  b.m("Loess", "trace_hat", {}, "float", [](void* s, A&) { return Value::real(self<LoessFit>(s).trace_hat); });
  b.m("Loess", "one_delta", {}, "float", [](void* s, A&) { return Value::real(self<LoessFit>(s).one_delta); });
  b.m("Loess", "two_delta", {}, "float", [](void* s, A&) { return Value::real(self<LoessFit>(s).two_delta); });
}

}  // namespace rcli::host
