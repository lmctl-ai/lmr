// Dist (d/p/q for 17 families), Rng (seeded draws and sampling) and the special functions.
#include "host.hpp"
#include "rcli/dist.hpp"
#include "rcli/rng.hpp"
#include "rcli/special.hpp"

namespace rcli::host {

void bind_dist(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto dist = [&h](Dist d) { return h.hold("Dist", std::move(d)); };
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  // ---- constructors, R parameterisations ----
  b.fn("norm", {}, "Dist", [dist](A&) { return dist(Dist::norm(0, 1)); });
  b.fn("norm", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::norm(a[0].f, a[1].f)); });
  b.fn("lnorm", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::lnorm(a[0].f, a[1].f)); });
  b.fn("t_dist", {"float"}, "Dist", [dist](A& a) { return dist(Dist::t(a[0].f)); });
  b.fn("chisq", {"float"}, "Dist", [dist](A& a) { return dist(Dist::chisq(a[0].f)); });
  b.fn("f_dist", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::f(a[0].f, a[1].f)); });
  // Noncentral laws (the last argument is the noncentrality) and the studentized range (p and q only).
  b.fn("t_dist", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::t(a[0].f, a[1].f)); });
  b.fn("chisq", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::chisq(a[0].f, a[1].f)); });
  b.fn("f_dist", {"float", "float", "float"}, "Dist", [dist](A& a) { return dist(Dist::f(a[0].f, a[1].f, a[2].f)); });
  b.fn("tukey", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::tukey(a[0].f, a[1].f)); });
  b.fn("tukey", {"float", "float", "float"}, "Dist", [dist](A& a) { return dist(Dist::tukey(a[0].f, a[1].f, a[2].f)); });
  b.fn("gamma_dist", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::gamma(a[0].f, a[1].f)); });
  b.fn("beta_dist", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::beta(a[0].f, a[1].f)); });
  b.fn("exp_dist", {"float"}, "Dist", [dist](A& a) { return dist(Dist::exp(a[0].f)); });
  b.fn("unif", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::unif(a[0].f, a[1].f)); });
  b.fn("weibull", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::weibull(a[0].f, a[1].f)); });
  b.fn("cauchy", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::cauchy(a[0].f, a[1].f)); });
  b.fn("logis", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::logis(a[0].f, a[1].f)); });
  b.fn("binom", {"int", "float"}, "Dist", [dist](A& a) { return dist(Dist::binom(static_cast<double>(a[0].i), a[1].f)); });
  b.fn("pois", {"float"}, "Dist", [dist](A& a) { return dist(Dist::pois(a[0].f)); });
  b.fn("geom", {"float"}, "Dist", [dist](A& a) { return dist(Dist::geom(a[0].f)); });
  b.fn("nbinom", {"float", "float"}, "Dist", [dist](A& a) { return dist(Dist::nbinom(a[0].f, a[1].f)); });
  b.fn("hyper", {"int", "int", "int"}, "Dist", [dist](A& a) {
    return dist(Dist::hyper(static_cast<double>(a[0].i), static_cast<double>(a[1].i), static_cast<double>(a[2].i)));
  });
  b.fn("dist", {"string", "float"}, "Dist", [dist](A& a) { return dist(Dist::by_name(a[0].s, a[1].f)); });
  b.fn("dist", {"string", "float", "float"}, "Dist", [dist](A& a) { return dist(Dist::by_name(a[0].s, a[1].f, a[2].f)); });
  b.fn("dist", {"string", "float", "float", "float"}, "Dist", [dist](A& a) {
    return dist(Dist::by_name(a[0].s, a[1].f, a[2].f, a[3].f));
  });

  // ---- evaluation ----
  b.m("Dist", "name", {}, "string", [](void* s, A&) { return Value::string(self<Dist>(s).name()); });
  b.m("Dist", "d", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).d(a[0].f)); });
  b.m("Dist", "log_d", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).d(a[0].f, true)); });
  b.m("Dist", "p", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).p(a[0].f)); });
  b.m("Dist", "p_upper", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).p(a[0].f, false)); });
  b.m("Dist", "log_p", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).p(a[0].f, true, true)); });
  b.m("Dist", "log_p_upper", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).p(a[0].f, false, true)); });
  b.m("Dist", "q", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).q(a[0].f)); });
  b.m("Dist", "q_upper", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).q(a[0].f, false)); });
  b.m("Dist", "q_log", {"float"}, "float", [](void* s, A& a) { return Value::real(self<Dist>(s).q(a[0].f, true, true)); });
  b.m("Dist", "mean", {}, "float", [](void* s, A&) { return Value::real(self<Dist>(s).mean()); });
  b.m("Dist", "var", {}, "float", [](void* s, A&) { return Value::real(self<Dist>(s).var()); });
  b.m("Dist", "sd", {}, "float", [](void* s, A&) { return Value::real(std::sqrt(self<Dist>(s).var())); });
  b.m("Dist", "is_discrete", {}, "bool", [](void* s, A&) { return Value::boolean(self<Dist>(s).discrete()); });
  auto vmap = [&b, vec](const char* name, double (*f)(const Dist&, double)) {
    b.m("Dist", name, {"Vec"}, "Vec", [vec, f](void* s, A& a) {
      const Dist& d = self<Dist>(s);
      Vec out = obj<Vec>(a[0]);
      for (double& x : out) x = f(d, x);
      return vec(std::move(out));
    });
  };
  vmap("d_vec", [](const Dist& d, double x) { return d.d(x); });
  vmap("p_vec", [](const Dist& d, double x) { return d.p(x); });
  vmap("q_vec", [](const Dist& d, double x) { return d.q(x); });

  // ---- Rng ----
  b.fn("rng", {"int"}, "Rng", [&h](A& a) { return h.hold("Rng", Rng(static_cast<uint64_t>(a[0].i))); });
  b.m("Rng", "unif", {}, "float", [](void* s, A&) { return Value::real(self<Rng>(s).unif()); });
  b.m("Rng", "norm", {}, "float", [](void* s, A&) { return Value::real(self<Rng>(s).norm()); });
  b.m("Rng", "below", {"int"}, "int", [](void* s, A& a) {
    if (a[0].i < 1) throw dslvm::Error("Rng.below: n must be at least 1");
    return Value::integer(static_cast<int64_t>(self<Rng>(s).below(static_cast<uint64_t>(a[0].i))));
  });
  b.m("Rng", "draw", {"Dist"}, "float", [](void* s, A& a) { return Value::real(self<Rng>(s).draw(obj<Dist>(a[0]))); });
  b.m("Rng", "draw", {"Dist", "int"}, "Vec", [vec](void* s, A& a) {
    if (a[1].i < 0 || a[1].i > 100000000) throw dslvm::Error("Rng.draw: n must be between 0 and 1e8");
    return vec(self<Rng>(s).draw(obj<Dist>(a[0]), static_cast<size_t>(a[1].i)));
  });
  auto sample = [vec](Rng& r, const Vec& x, int64_t size, bool replace, const Vec* prob) {
    if (size < 0) throw dslvm::Error("Rng.sample: size must be nonnegative");
    if (prob && prob->size() != x.size()) throw dslvm::Error("Rng.sample: prob must have the length of x");
    auto idx = r.sample(x.size(), static_cast<size_t>(size), replace, prob);
    Vec out(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) out[i] = x[idx[i]];
    return vec(out);
  };
  b.m("Rng", "sample", {"Vec", "int", "bool"}, "Vec", [sample](void* s, A& a) {
    return sample(self<Rng>(s), obj<Vec>(a[0]), a[1].i, a[2].b, nullptr);
  });
  b.m("Rng", "sample", {"Vec", "int", "bool", "Vec"}, "Vec", [sample](void* s, A& a) {
    return sample(self<Rng>(s), obj<Vec>(a[0]), a[1].i, a[2].b, &obj<Vec>(a[3]));
  });
  b.m("Rng", "permutation", {"int"}, "Vec", [vec](void* s, A& a) {
    if (a[0].i < 0) throw dslvm::Error("Rng.permutation: n must be nonnegative");
    size_t n = static_cast<size_t>(a[0].i);
    auto idx = self<Rng>(s).sample(n, n, false, nullptr);
    Vec out(n);
    for (size_t i = 0; i < n; ++i) out[i] = static_cast<double>(idx[i]);
    return vec(out);
  });

  // ---- special functions (R names, snake case for the Bessel family) ----
  auto s1 = [&b](const char* name, double (*f)(double)) {
    b.fn(name, {"float"}, "float", [f](A& a) { return Value::real(f(a[0].f)); });
  };
  auto s2 = [&b](const char* name, double (*f)(double, double)) {
    b.fn(name, {"float", "float"}, "float", [f](A& a) { return Value::real(f(a[0].f, a[1].f)); });
  };
  s1("gamma", special::gamma);
  s1("lgamma", special::lgamma);
  s1("digamma", special::digamma);
  s1("trigamma", special::trigamma);
  s1("factorial", special::factorial);
  s1("lfactorial", special::lfactorial);
  s1("erf", special::erf);
  s1("erfc", special::erfc);
  s2("beta", special::beta);
  s2("lbeta", special::lbeta);
  s2("choose", special::choose);
  s2("lchoose", special::lchoose);
  s2("bessel_j", special::bessel_j);
  s2("bessel_y", special::bessel_y);
  b.fn("bessel_i", {"float", "float"}, "float", [](A& a) { return Value::real(special::bessel_i(a[0].f, a[1].f)); });
  b.fn("bessel_k", {"float", "float"}, "float", [](A& a) { return Value::real(special::bessel_k(a[0].f, a[1].f)); });
  b.fn("bessel_i", {"float", "float", "bool"}, "float", [](A& a) { return Value::real(special::bessel_i(a[0].f, a[1].f, a[2].b)); });
  b.fn("bessel_k", {"float", "float", "bool"}, "float", [](A& a) { return Value::real(special::bessel_k(a[0].f, a[1].f, a[2].b)); });
  b.fn("gamma_p", {"float", "float"}, "float", [](A& a) { return Value::real(special::gamma_p(a[0].f, a[1].f)); });
  b.fn("beta_inc", {"float", "float", "float"}, "float", [](A& a) { return Value::real(special::beta_inc(a[0].f, a[1].f, a[2].f)); });
}

}  // namespace rcli::host
