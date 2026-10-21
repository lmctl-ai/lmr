// acf/pacf (Vec methods), ARIMA (spec builder -> fit -> forecast), arima_sim and STL.
#include "host.hpp"
#include "rcli/dist.hpp"
#include "rcli/models.hpp"
#include "rcli/rng.hpp"
#include "rcli/ts.hpp"

namespace rcli::host {

namespace {

struct ArimaModel {
  ArimaFit fit;
  CoefTable table;  // estimate, se, z = estimate / se, two-sided normal p-value
};

const Vec& as_ts(const Value& v) { return obj<Vec>(v); }

double level_q(double level, const char* what) {
  if (!(level > 0 && level < 1)) throw dslvm::Error(std::string(what) + ": level must be between 0 and 1");
  return Dist::norm(0, 1).q(0.5 + level / 2);
}

}  // namespace

void bind_ts(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  b.m("Vec", "acf", {}, "Vec", [vec](void* s, A&) { return vec(acf(self<Vec>(s))); });
  b.m("Vec", "acf", {"int"}, "Vec", [vec](void* s, A& a) { return vec(acf(self<Vec>(s), static_cast<int>(a[0].i))); });
  b.m("Vec", "pacf", {}, "Vec", [vec](void* s, A&) { return vec(pacf(self<Vec>(s))); });
  b.m("Vec", "pacf", {"int"}, "Vec", [vec](void* s, A& a) { return vec(pacf(self<Vec>(s), static_cast<int>(a[0].i))); });

  // ---- ARIMA ----
  b.fn("arima", {"int", "int", "int"}, "ArimaSpec", [&h](A& a) {
    if (a[0].i < 0 || a[1].i < 0 || a[2].i < 0 || a[0].i > 50 || a[1].i > 3 || a[2].i > 50)
      throw dslvm::Error("arima: need 0 <= p <= 50, 0 <= d <= 3, 0 <= q <= 50");
    ArimaSpec s;
    s.p = static_cast<int>(a[0].i);
    s.d = static_cast<int>(a[1].i);
    s.q = static_cast<int>(a[2].i);
    return h.hold("ArimaSpec", s);
  });
  b.m("ArimaSpec", "seasonal", {"int", "int", "int", "int"}, "ArimaSpec", [](void* s, A& a) {
    auto& sp = self<ArimaSpec>(s);
    if (a[0].i < 0 || a[1].i < 0 || a[2].i < 0 || a[0].i > 10 || a[1].i > 2 || a[2].i > 10 || a[3].i < 2)
      throw dslvm::Error("ArimaSpec.seasonal: need 0 <= P <= 10, 0 <= D <= 2, 0 <= Q <= 10 and period >= 2");
    sp.P = static_cast<int>(a[0].i);
    sp.D = static_cast<int>(a[1].i);
    sp.Q = static_cast<int>(a[2].i);
    sp.period = static_cast<int>(a[3].i);
    return Value::object("ArimaSpec", s);
  });
  b.m("ArimaSpec", "method", {"string"}, "ArimaSpec", [](void* s, A& a) {
    if (a[0].s != "CSS-ML" && a[0].s != "ML" && a[0].s != "CSS") throw dslvm::Error("ArimaSpec.method: must be CSS-ML, ML or CSS");
    self<ArimaSpec>(s).method = a[0].s;
    return Value::object("ArimaSpec", s);
  });
  b.m("ArimaSpec", "include_mean", {"bool"}, "ArimaSpec", [](void* s, A& a) {
    self<ArimaSpec>(s).include_mean = a[0].b;
    return Value::object("ArimaSpec", s);
  });
  b.m("ArimaSpec", "fit", {"Vec"}, "Arima", [&h](void* s, A& a) {
    ArimaModel m;
    m.fit = arima(as_ts(a[0]), self<ArimaSpec>(s));
    m.table.names = m.fit.coef_names;
    m.table.estimate = m.fit.coef;
    m.table.std_error = m.fit.se;
    for (size_t i = 0; i < m.fit.coef.size(); ++i) {
      double z = m.fit.coef[i] / m.fit.se[i];
      m.table.statistic.push_back(z);
      m.table.p_value.push_back(std::isfinite(z) ? 2 * Dist::norm(0, 1).p(-std::fabs(z)) : kNaN);
      m.table.aliased.push_back(false);
    }
    return h.hold("Arima", std::move(m));
  });
  b.m("Arima", "n_coef", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<ArimaModel>(s).table.names.size())); });
  b.m("Arima", "name", {"int"}, "string", [](void* s, A& a) {
    auto& t = self<ArimaModel>(s).table;
    return Value::string(t.names[static_cast<size_t>(as_index(a[0].i, t.names.size(), "Arima.name"))]);
  });
  for (const char* key : {"int", "string"}) {
    auto idx = [](const CoefTable& t, const Value& v) -> size_t {
      if (v.kind == Value::INT) return static_cast<size_t>(as_index(v.i, t.names.size(), "Arima"));
      int k = t.index(v.s);
      if (k < 0) {
        std::string have;
        for (const auto& n : t.names) have += (have.empty() ? "" : ", ") + n;
        throw dslvm::Error("Arima: no coefficient '" + v.s + "' (coefficients: " + have + ")");
      }
      return static_cast<size_t>(k);
    };
    b.m("Arima", "estimate", {key}, "float", [idx](void* s, A& a) { auto& t = self<ArimaModel>(s).table; return Value::real(t.estimate[idx(t, a[0])]); });
    b.m("Arima", "se", {key}, "float", [idx](void* s, A& a) { auto& t = self<ArimaModel>(s).table; return Value::real(t.std_error[idx(t, a[0])]); });
    b.m("Arima", "statistic", {key}, "float", [idx](void* s, A& a) { auto& t = self<ArimaModel>(s).table; return Value::real(t.statistic[idx(t, a[0])]); });
    b.m("Arima", "p_value", {key}, "float", [idx](void* s, A& a) { auto& t = self<ArimaModel>(s).table; return Value::real(t.p_value[idx(t, a[0])]); });
  }
  b.m("Arima", "has_coef", {"string"}, "bool", [](void* s, A& a) { return Value::boolean(self<ArimaModel>(s).table.index(a[0].s) >= 0); });
  b.m("Arima", "sigma2", {}, "float", [](void* s, A&) { return Value::real(self<ArimaModel>(s).fit.sigma2); });
  b.m("Arima", "loglik", {}, "float", [](void* s, A&) { return Value::real(self<ArimaModel>(s).fit.loglik); });
  b.m("Arima", "aic", {}, "float", [](void* s, A&) { return Value::real(self<ArimaModel>(s).fit.aic); });
  b.m("Arima", "nobs", {}, "int", [](void* s, A&) { return Value::integer(self<ArimaModel>(s).fit.nobs); });
  b.m("Arima", "convergence", {}, "int", [](void* s, A&) { return Value::integer(self<ArimaModel>(s).fit.convergence); });
  b.m("Arima", "residuals", {}, "Vec", [vec](void* s, A&) { return vec(self<ArimaModel>(s).fit.residuals); });
  b.m("Arima", "forecast", {"int"}, "Forecast", [&h](void* s, A& a) {
    if (a[0].i < 1 || a[0].i > 100000) throw dslvm::Error("Arima.forecast: horizon must be between 1 and 100000");
    return h.hold("Forecast", self<ArimaModel>(s).fit.forecast(static_cast<int>(a[0].i)));
  });
  b.m("Forecast", "len", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<ArimaForecast>(s).mean.size())); });
  b.m("Forecast", "mean", {}, "Vec", [vec](void* s, A&) { return vec(self<ArimaForecast>(s).mean); });
  b.m("Forecast", "se", {}, "Vec", [vec](void* s, A&) { return vec(self<ArimaForecast>(s).se); });
  b.m("Forecast", "mean_at", {"int"}, "float", [](void* s, A& a) {
    auto& f = self<ArimaForecast>(s);
    return Value::real(f.mean[static_cast<size_t>(as_index(a[0].i, f.mean.size(), "Forecast.mean_at"))]);
  });
  b.m("Forecast", "se_at", {"int"}, "float", [](void* s, A& a) {
    auto& f = self<ArimaForecast>(s);
    return Value::real(f.se[static_cast<size_t>(as_index(a[0].i, f.se.size(), "Forecast.se_at"))]);
  });
  for (int up = 0; up < 2; ++up) {
    b.m("Forecast", up ? "upper" : "lower", {"int", "float"}, "float", [up](void* s, A& a) {
      auto& f = self<ArimaForecast>(s);
      size_t i = static_cast<size_t>(as_index(a[0].i, f.mean.size(), "Forecast interval"));
      double q = level_q(a[1].f, "Forecast interval");
      return Value::real(f.mean[i] + (up ? q : -q) * f.se[i]);
    });
  }

  // ---- arima_sim ----
  b.fn("arima_sim", {"Rng", "int", "Vec", "Vec", "float"}, "Vec", [vec](A& a) {
    if (a[1].i < 1 || a[1].i > 10000000) throw dslvm::Error("arima_sim: n must be between 1 and 1e7");
    return vec(arima_sim(obj<Rng>(a[0]), static_cast<int>(a[1].i), obj<Vec>(a[2]), obj<Vec>(a[3]), a[4].f));
  });
  b.fn("arima_sim", {"Rng", "int", "Vec", "Vec", "float", "int"}, "Vec", [vec](A& a) {
    if (a[1].i < 1 || a[1].i > 10000000) throw dslvm::Error("arima_sim: n must be between 1 and 1e7");
    if (a[5].i < 0 || a[5].i > 3) throw dslvm::Error("arima_sim: d must be between 0 and 3");
    return vec(arima_sim(obj<Rng>(a[0]), static_cast<int>(a[1].i), obj<Vec>(a[2]), obj<Vec>(a[3]), a[4].f, static_cast<int>(a[5].i)));
  });

  // ---- STL ----
  auto run_stl = [&h](const Vec& x, int64_t period, StlOptions o) {
    if (period < 2 || period > 100000) throw dslvm::Error("stl: period must be at least 2");
    return h.hold("Stl", stl(x, static_cast<int>(period), o));
  };
  b.fn("stl", {"Vec", "int", "int"}, "Stl", [run_stl](A& a) {
    StlOptions o;
    o.s_window = static_cast<int>(a[2].i);
    return run_stl(obj<Vec>(a[0]), a[1].i, o);
  });
  b.fn("stl", {"Vec", "int", "int", "bool"}, "Stl", [run_stl](A& a) {
    StlOptions o;
    o.s_window = static_cast<int>(a[2].i);
    o.robust = a[3].b;
    return run_stl(obj<Vec>(a[0]), a[1].i, o);
  });
  b.fn("stl_periodic", {"Vec", "int"}, "Stl", [run_stl](A& a) {
    StlOptions o;
    o.periodic = true;
    return run_stl(obj<Vec>(a[0]), a[1].i, o);
  });
  b.fn("stl_periodic", {"Vec", "int", "bool"}, "Stl", [run_stl](A& a) {
    StlOptions o;
    o.periodic = true;
    o.robust = a[2].b;
    return run_stl(obj<Vec>(a[0]), a[1].i, o);
  });
  b.m("Stl", "seasonal", {}, "Vec", [vec](void* s, A&) { return vec(self<StlFit>(s).seasonal); });
  b.m("Stl", "trend", {}, "Vec", [vec](void* s, A&) { return vec(self<StlFit>(s).trend); });
  b.m("Stl", "remainder", {}, "Vec", [vec](void* s, A&) { return vec(self<StlFit>(s).remainder); });
  b.m("Stl", "weights", {}, "Vec", [vec](void* s, A&) { return vec(self<StlFit>(s).weights); });
}

}  // namespace rcli::host
