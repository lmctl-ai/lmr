// Kernel density estimation, bandwidth rules, the empirical CDF, histograms and pretty breakpoints.
#include "host.hpp"
#include "rcli/kde.hpp"

namespace rcli::host {

void bind_kde(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  b.fn("bw", {"Vec", "string"}, "float", [](A& a) { return Value::real(kde::bandwidth(obj<Vec>(a[0]), a[1].s)); });

  auto dens = [&h](const Vec& x, const std::string& rule, double bw, const std::string& kernel, double adjust) {
    kde::DensityOptions o;
    o.bw_rule = rule;
    o.bw = bw;
    o.kernel = kernel;
    o.adjust = adjust;
    return h.hold("Density", kde::density(x, o));
  };
  b.fn("density", {"Vec"}, "Density", [dens](A& a) { return dens(obj<Vec>(a[0]), "nrd0", kNaN, "gaussian", 1.0); });
  b.fn("density", {"Vec", "string"}, "Density", [dens](A& a) { return dens(obj<Vec>(a[0]), a[1].s, kNaN, "gaussian", 1.0); });
  b.fn("density", {"Vec", "string", "string"}, "Density", [dens](A& a) { return dens(obj<Vec>(a[0]), a[1].s, kNaN, a[2].s, 1.0); });
  b.fn("density", {"Vec", "float", "string"}, "Density", [dens](A& a) {
    if (!(a[1].f > 0)) throw dslvm::Error("density: bw must be positive");
    return dens(obj<Vec>(a[0]), "nrd0", a[1].f, a[2].s, 1.0);
  });
  b.fn("density", {"Vec", "string", "string", "float"}, "Density", [dens](A& a) {
    if (!(a[3].f > 0)) throw dslvm::Error("density: adjust must be positive");
    return dens(obj<Vec>(a[0]), a[1].s, kNaN, a[2].s, a[3].f);
  });
  b.m("Density", "x", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Density>(s).x); });
  b.m("Density", "y", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Density>(s).y); });
  b.m("Density", "bw", {}, "float", [](void* s, A&) { return Value::real(self<kde::Density>(s).bw); });
  b.m("Density", "kernel", {}, "string", [](void* s, A&) { return Value::string(self<kde::Density>(s).kernel); });
  b.m("Density", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<kde::Density>(s).at(a[0].f)); });
  b.m("Density", "eval", {"Vec"}, "Vec", [vec](void* s, A& a) {
    Vec out;
    for (double x : obj<Vec>(a[0])) out.push_back(self<kde::Density>(s).at(x));
    return vec(std::move(out));
  });

  b.fn("ecdf", {"Vec"}, "Ecdf", [&h](A& a) { return h.hold("Ecdf", kde::ecdf(obj<Vec>(a[0]))); });
  b.m("Ecdf", "at", {"float"}, "float", [](void* s, A& a) { return Value::real(self<kde::Ecdf>(s)(a[0].f)); });
  b.m("Ecdf", "eval", {"Vec"}, "Vec", [vec](void* s, A& a) { return vec(self<kde::Ecdf>(s)(obj<Vec>(a[0]))); });
  b.m("Ecdf", "knots", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Ecdf>(s).knots); });
  b.m("Ecdf", "quantile", {"float"}, "float", [](void* s, A& a) {
    if (!(a[0].f >= 0 && a[0].f <= 1)) throw dslvm::Error("Ecdf.quantile: p must be in [0, 1]");
    return Value::real(self<kde::Ecdf>(s).quantile(a[0].f));
  });

  b.fn("hist", {"Vec"}, "Hist", [&h](A& a) { return h.hold("Hist", kde::hist(obj<Vec>(a[0]))); });
  b.fn("hist", {"Vec", "string"}, "Hist", [&h](A& a) { return h.hold("Hist", kde::hist(obj<Vec>(a[0]), a[1].s)); });
  b.fn("hist", {"Vec", "Vec"}, "Hist", [&h](A& a) { return h.hold("Hist", kde::hist(obj<Vec>(a[0]), obj<Vec>(a[1]))); });
  b.m("Hist", "breaks", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Histogram>(s).breaks); });
  b.m("Hist", "counts", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Histogram>(s).counts); });
  b.m("Hist", "density", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Histogram>(s).density); });
  b.m("Hist", "mids", {}, "Vec", [vec](void* s, A&) { return vec(self<kde::Histogram>(s).mids); });

  b.fn("pretty", {"float", "float", "int"}, "Vec", [vec](A& a) {
    if (a[2].i < 1 || a[2].i > 10000) throw dslvm::Error("pretty: n must be between 1 and 10000");
    return vec(kde::pretty(a[0].f, a[1].f, static_cast<int>(a[2].i)));
  });
}

}  // namespace rcli::host
