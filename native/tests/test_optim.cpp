// Tests for optimisation, root finding and integration. Expected values come from R 4.3.3 (Rscript, digits = 17);
// each expectation names the R expression it reproduces.
#include <cmath>
#include <memory>
#include <string>

#include "check.hpp"
#include "rcli/optim.hpp"
#include "rcli/rng.hpp"

using namespace rcli;

namespace {

const double kPi = 3.14159265358979323846;

double dnorm(double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * kPi); }

// fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2   (the ?optim Rosenbrock example)
double fr(const Vec& x) { return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0]); }
// grr <- function(x) c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]), 200 * (x[2] - x[1]^2))
Vec grr(const Vec& x) {
  return {-400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]), 200.0 * (x[1] - x[0] * x[0])};
}
// q <- function(x) 2*x[1]^2 + x[1]*x[2] + 3*x[2]^2 - x[1] + 2*x[2]   (Hessian [[4,1],[1,6]])
double quad(const Vec& x) { return 2 * x[0] * x[0] + x[0] * x[1] + 3 * x[1] * x[1] - x[0] + 2 * x[1]; }

// Error estimates at the roundoff level depend on whether the compiler contracts a*b+c into FMA instructions
// (GCC does at -O2 on aarch64, as did the R build used for the expectations), so they get a looser tolerance.
void check_integral(const std::string& what, const IntegrateResult& r, double value, double abs_error, int subdiv,
                    int ier, double err_tol = 0.0) {
  if (err_tol == 0.0) err_tol = abs_error < 1e-10 ? 1e-3 : 1e-6;
  check::rel(r.value, value, 1e-13, what + " value");
  check::rel(r.abs_error, abs_error, err_tol, what + " abs.error");
  check::ok(r.subdivisions == subdiv, what + " subdivisions " + std::to_string(r.subdivisions));
  check::ok(r.ier == ier, what + " ier " + std::to_string(r.ier));
}

void test_integrate() {
  // integrate(dnorm, -Inf, Inf): 0.99999999979415932, abs.error 9.3611240867000345e-05, 3 subdivisions
  auto r = integrate(dnorm, -kInf, kInf);
  check_integral("integrate dnorm (-Inf,Inf)", r, 0.99999999979415932, 9.3611240867000345e-05, 3, 0);
  check::ok(r.message == "OK", "integrate message OK");
  check::ok(r.evaluations == 150, "integrate dnorm (-Inf,Inf) evaluations (QAGI: 2 x 15 x 5)");
  // integrate(dnorm, 0, Inf)
  check_integral("integrate dnorm (0,Inf)", integrate(dnorm, 0, kInf), 0.49999999989707966, 4.6805620433500172e-05,
                 3, 0);
  // integrate(dnorm, -Inf, 0)
  check_integral("integrate dnorm (-Inf,0)", integrate(dnorm, -kInf, 0), 0.49999999989707966,
                 4.6805620433500172e-05, 3, 0);
  // integrate(dnorm, 0, -Inf): R maps a finite lower limit to (lower, Inf)
  check_integral("integrate dnorm (0,-Inf)", integrate(dnorm, 0, -kInf), 0.49999999989707966,
                 4.6805620433500172e-05, 3, 0);
  // integrate(function(x) 1/((x+1)*sqrt(x)), 0, Inf)  (= pi)
  check_integral("integrate 1/((x+1)sqrt(x))", integrate([](double x) { return 1 / ((x + 1) * std::sqrt(x)); }, 0, kInf),
                 3.1415926530279124, 2.6514999261539884e-05, 10, 0);
  // integrate(function(x) exp(-x^2/2), -1.96, 1.96)
  r = integrate([](double x) { return std::exp(-x * x / 2); }, -1.96, 1.96);
  check_integral("integrate exp(-x^2/2)", r, 2.3813074130614198, 2.6254671753903046e-11, 1, 0);
  check::ok(r.evaluations == 21, "integrate one QK21 panel = 21 evaluations");
  // integrate(function(x) x^-0.5, 0, 1)
  r = integrate([](double x) { return std::pow(x, -0.5); }, 0, 1);
  check_integral("integrate x^-0.5", r, 1.9999999999999993, 3.7747582837255322e-15, 6, 0);
  check::ok(r.evaluations == 42 * 6 - 21, "integrate x^-0.5 evaluations (QUADPACK neval)");
  // integrate(function(x) x^-0.9, 0, 1)
  check_integral("integrate x^-0.9", integrate([](double x) { return std::pow(x, -0.9); }, 0, 1), 10.000000000000002,
                 6.5014660322049167e-13, 6, 0, 0.1);  // 6.80e-13 without FMA contraction
  // integrate(function(x) log(x), 0, 1)
  check_integral("integrate log", integrate([](double x) { return std::log(x); }, 0, 1), -0.99999999999999989,
                 1.1102230246251563e-15, 6, 0);
  // integrate(function(x) exp(x), 0, 1, rel.tol = 1e-10)
  IntegrateOptions o;
  o.rel_tol = 1e-10;
  check_integral("integrate exp rel.tol", integrate([](double x) { return std::exp(x); }, 0, 1, o),
                 1.7182818284590453, 1.9076760487502457e-14, 1, 0);
  // integrate(function(x) 1/(1+x^2), -Inf, 3)
  check_integral("integrate cauchy (-Inf,3)", integrate([](double x) { return 1 / (1 + x * x); }, -kInf, 3),
                 2.819842099197106, 4.9157283347289835e-06, 5, 0);
  // integrate(dnorm, 1, 0)  (reversed limits)
  check_integral("integrate dnorm (1,0)", integrate(dnorm, 1, 0), -0.34134474606854293, 3.7896879642012371e-15, 1, 0);
  // integrate(function(x) cos(x)*exp(-x), 0, Inf)
  check_integral("integrate cos exp", integrate([](double x) { return std::cos(x) * std::exp(-x); }, 0, kInf),
                 0.50000000000705014, 4.3219218420461128e-07, 5, 0);
  // integrate(function(x) 1/x, 0, 1, stop.on.error = FALSE): "maximum number of subdivisions reached"
  r = integrate([](double x) { return 1 / x; }, 0, 1);
  check_integral("integrate 1/x", r, 76.334199703385366, 9.3505603731408655, 100, 1);
  check::ok(r.message == "maximum number of subdivisions reached", "integrate 1/x message");
  // integrate(function(x) sin(x)/x, 0, Inf, stop.on.error = FALSE): "the integral is probably divergent"
  r = integrate([](double x) { return std::sin(x) / x; }, 0, kInf);
  check_integral("integrate sin(x)/x", r, 1.861015277758469, 0.31605488593933129, 100, 5);
  check::ok(r.message == "the integral is probably divergent", "integrate divergent message");
  // integrate(function(x) sin(1/x), 0.001, 1, subdivisions = 5, stop.on.error = FALSE)
  o = IntegrateOptions();
  o.subdivisions = 5;
  check_integral("integrate sin(1/x) limit 5", integrate([](double x) { return std::sin(1 / x); }, 0.001, 1, o),
                 0.49713240188691843, 0.034023898151482046, 5, 1);
  // integrate(function(x) abs(sin(10*x)), 0, 10, subdivisions = 10, stop.on.error = FALSE)
  o.subdivisions = 10;
  check_integral("integrate |sin(10x)| limit 10",
                 integrate([](double x) { return std::fabs(std::sin(10 * x)); }, 0, 10, o), 6.3613080396164596,
                 2.6172287000863621, 10, 1);
  // integrate(function(x) abs(sin(10*x)), 0, 10, stop.on.error = FALSE)
  check_integral("integrate |sin(10x)|", integrate([](double x) { return std::fabs(std::sin(10 * x)); }, 0, 10),
                 6.3865681434123758, 0.0077874832381743531, 100, 1);
  // errors: non-finite function value; invalid parameters; NA limits
  check::throws([] { integrate([](double x) { return x > 0.5 ? kNaN : x; }, 0, 1); }, "integrate non-finite value");
  check::throws([] { integrate([](double x) { return x > 0.5 ? kInf : x; }, 0, 1); }, "integrate infinite value");
  o = IntegrateOptions();
  o.subdivisions = 0;
  check::throws([&] { integrate(dnorm, 0, 1, o); }, "integrate subdivisions < 1");
  o = IntegrateOptions();
  o.abs_tol = 0;
  o.rel_tol = 1e-20;
  check::throws([&] { integrate(dnorm, 0, 1, o); }, "integrate invalid tolerances");
  check::throws([] { integrate(dnorm, kNaN, 1); }, "integrate NaN limit");
}

void test_optimize() {
  // optimize(sin, c(0, 2*pi)): $minimum 4.7123910812000886, $objective -0.99999999999779332 (10 calls of f)
  auto r = optimize([](double x) { return std::sin(x); }, 0, 2 * kPi);
  check::close(r.minimum, 4.7123910812000886, 1e-12, "optimize sin minimum");
  check::close(r.objective, -0.99999999999779332, 1e-14, "optimize sin objective");
  check::ok(r.evaluations == 10, "optimize sin evaluations");
  // optimize(sin, c(0, 2*pi), maximum = TRUE): $maximum 1.5707942259794982
  r = optimize([](double x) { return std::sin(x); }, 0, 2 * kPi, true);
  check::close(r.minimum, 1.5707942259794982, 1e-12, "optimize sin maximum");
  check::close(r.objective, 0.99999999999779332, 1e-14, "optimize sin maximum objective");
  // g <- function(x) (x^2-2)^2 + x ; optimize(g, c(-3, 3)): -1.4729957803549456 (13 calls)
  auto g = [](double x) { return (x * x - 2) * (x * x - 2) + x; };
  r = optimize(g, -3, 3);
  check::close(r.minimum, -1.4729957803549456, 1e-12, "optimize quartic minimum");
  check::close(r.objective, -1.4441920665810002, 1e-14, "optimize quartic objective");
  check::ok(r.evaluations == 13, "optimize quartic evaluations");
  // optimize(g, c(-3, 3), tol = 1e-10): -1.4729976022914812 (14 calls)
  r = optimize(g, -3, 3, false, 1e-10);
  check::close(r.minimum, -1.4729976022914812, 1e-12, "optimize quartic tol=1e-10");
  check::ok(r.evaluations == 14, "optimize quartic tol=1e-10 evaluations");
  // optimize(g, c(0, 3)): the other local minimum 1.3470075772653312 (12 calls)
  r = optimize(g, 0, 3);
  check::close(r.minimum, 1.3470075772653312, 1e-12, "optimize quartic local minimum");
  check::ok(r.evaluations == 12, "optimize quartic (0,3) evaluations");
  // optimize(function(x) abs(x - 0.3), c(0, 1)): 0.29999278971096544 (15 calls)
  r = optimize([](double x) { return std::fabs(x - 0.3); }, 0, 1);
  check::close(r.minimum, 0.29999278971096544, 1e-12, "optimize |x-0.3|");
  check::ok(r.evaluations == 15, "optimize |x-0.3| evaluations");
  // optimize(function(x) x, c(1, 2)): boundary minimum 1.0000661069613521 (21 calls)
  r = optimize([](double x) { return x; }, 1, 2);
  check::close(r.minimum, 1.0000661069613521, 1e-12, "optimize monotone");
  check::ok(r.evaluations == 21, "optimize monotone evaluations");
  // non-finite values are replaced by the largest double, as R does (with a warning):
  // optimize(function(x) if (x > 0.5) Inf else (x-0.3)^2, c(0,1)): minimum 0.29999999999999999, objective 0
  r = optimize([](double x) { return x > 0.5 ? kInf : (x - 0.3) * (x - 0.3); }, 0, 1);
  check::close(r.minimum, 0.29999999999999999, 1e-12, "optimize with Inf values");
  // optimize(function(x) NaN, c(0,1)): minimum 0.9999338930386481, objective NaN
  r = optimize([](double) { return kNaN; }, 0, 1);
  check::close(r.minimum, 0.9999338930386481, 1e-12, "optimize all-NaN minimum");
  check::ok(std::isnan(r.objective), "optimize all-NaN objective NaN");
  check::throws([] { optimize([](double x) { return x; }, 2, 1); }, "optimize lower >= upper");
  check::throws([] { optimize([](double x) { return x; }, 0, kInf); }, "optimize infinite bound");
}

void test_uniroot() {
  auto cubic = [](double x) { return x * x * x - x - 1; };
  const double cubic_root = 1.3247179572447460;  // polyroot(c(-1, -1, 0, 1)), real root
  // uniroot(function(x) x^3 - x - 1, c(1, 2)): root 1.3247178251287683, iter 7, estim.prec 6.1035156250666134e-05
  auto r = uniroot(cubic, 1, 2);
  check::close(r.root, cubic_root, 1.3e-4, "uniroot cubic root within tol");
  check::close(r.root, 1.3247178251287683, 1e-4, "uniroot cubic root vs R");
  check::close(r.f_root, cubic(r.root), 1e-15, "uniroot cubic f.root");
  check::ok(r.iter >= 5 && r.iter <= 8, "uniroot cubic iter (R 7) " + std::to_string(r.iter));
  check::ok(r.init_it == -1, "uniroot cubic init.it NA");
  check::ok(r.estim_prec > 0 && r.estim_prec <= 1.3e-4, "uniroot cubic estim.prec");
  check::ok(r.converged, "uniroot cubic converged");
  check::ok(r.lower <= r.root && r.root <= r.upper, "uniroot cubic bracket");
  // uniroot(function(x) cos(x) - x, c(-pi, pi)): 0.73908430786659995, iter 7
  auto cosx = [](double x) { return std::cos(x) - x; };
  const double dottie = 0.73908513321516064;
  r = uniroot(cosx, -kPi, kPi);
  check::close(r.root, 0.73908430786659995, 1e-15, "uniroot cos(x)-x root (same path as R)");
  check::ok(r.iter == 7, "uniroot cos(x)-x iter 7");
  check::rel(r.estim_prec, 6.1035156250333067e-05, 1e-12, "uniroot cos(x)-x estim.prec");
  // uniroot(function(x) cos(x) - x, c(-pi, pi), tol = 1e-9): 0.73908513317531843, iter 8, estim.prec 5.0000037e-10
  RootOptions o;
  o.tol = 1e-9;
  r = uniroot(cosx, -kPi, kPi, o);
  check::close(r.root, dottie, 1e-9, "uniroot cos(x)-x tol=1e-9");
  check::close(r.root, 0.73908513317531843, 1e-15, "uniroot cos(x)-x tol=1e-9 vs R");
  check::ok(r.iter == 8, "uniroot cos(x)-x tol=1e-9 iter 8");
  check::rel(r.estim_prec, 5.0000037443709289e-10, 1e-12, "uniroot cos(x)-x tol=1e-9 estim.prec");
  // uniroot(function(x) x - 1/3, c(0, 1), tol = 1e-4): exact root on the first step, estim.prec 2/3, iter 1
  o.tol = 1e-4;
  r = uniroot([](double x) { return x - 1.0 / 3; }, 0, 1, o);
  check::close(r.root, 1.0 / 3, 1e-16, "uniroot x-1/3");
  check::ok(r.iter == 1 && r.f_root == 0, "uniroot x-1/3 iter 1, f.root 0");
  check::close(r.estim_prec, 0.66666666666666674, 1e-15, "uniroot x-1/3 estim.prec 2/3");
  // uniroot(function(x) x*(x^2-1) + .5, c(-2, 2), tol = 1e-10): -1.1914878839357335, iter 8
  o.tol = 1e-10;
  r = uniroot([](double x) { return x * (x * x - 1) + .5; }, -2, 2, o);
  check::close(r.root, -1.1914878839357335, 1e-14, "uniroot x(x^2-1)+.5");
  check::ok(r.iter == 8, "uniroot x(x^2-1)+.5 iter 8");
  // endpoint roots: uniroot(function(x) x - 1, c(1, 2)): root 1, iter 0, estim.prec 0
  r = uniroot([](double x) { return x - 1; }, 1, 2);
  check::ok(r.root == 1 && r.iter == 0 && r.estim_prec == 0, "uniroot root at lower end");
  r = uniroot([](double x) { return x - 2; }, 1, 2);
  check::ok(r.root == 2 && r.iter == 0 && r.estim_prec == 0, "uniroot root at upper end");
  // uniroot(..., maxiter = 3): warning "_NOT_ converged in 3 iterations", iter 3
  RootOptions m3;
  m3.maxiter = 3;
  r = uniroot(cubic, 1, 2, m3);
  check::ok(!r.converged && r.iter == 3, "uniroot maxiter 3 not converged");
  // extendInt = "yes": uniroot(f1, c(0,10), extendInt="yes"): root ~11, init.it 4 (R iter 12)
  RootOptions yes;
  yes.extend_int = "yes";
  r = uniroot([](double x) { return (121 - x * x) / (x * x + 1); }, 0, 10, yes);
  check::close(r.root, 11, 1e-5, "uniroot extendInt yes root");
  check::ok(r.init_it == 4, "uniroot extendInt yes init.it 4");
  check::ok(r.iter > r.init_it, "uniroot extendInt iter includes init.it");
  // uniroot(function(x) x^3 - x - 1, c(1.5, 2), extendInt = "yes"): init.it 4
  r = uniroot(cubic, 1.5, 2, yes);
  check::close(r.root, cubic_root, 1.3e-4, "uniroot cubic extendInt yes");
  check::ok(r.init_it == 4, "uniroot cubic extendInt yes init.it 4");
  // uniroot(function(x) exp(-x)*(x - 12), c(0, 2), extendInt = "upX"): root ~12, init.it 9
  RootOptions upx;
  upx.extend_int = "upX";
  r = uniroot([](double x) { return std::exp(-x) * (x - 12); }, 0, 2, upx);
  check::close(r.root, 12, 1e-5, "uniroot upX root");
  check::ok(r.init_it == 9, "uniroot upX init.it 9");
  // uniroot(function(x) x^3 - x - 1, c(-5, -3), extendInt = "upX"): init.it 8
  r = uniroot(cubic, -5, -3, upx);
  check::close(r.root, cubic_root, 1.3e-4, "uniroot cubic upX");
  check::ok(r.init_it == 8, "uniroot cubic upX init.it 8");
  // uniroot(function(x) -x + 5, c(0, 2), extendInt = "downX"): root 5 exactly, iter 9, init.it 8, estim.prec 5
  RootOptions downx;
  downx.extend_int = "downX";
  r = uniroot([](double x) { return -x + 5; }, 0, 2, downx);
  check::close(r.root, 5, 1e-15, "uniroot downX root");
  check::ok(r.init_it == 8 && r.iter == 9, "uniroot downX init.it 8, iter 9");
  check::close(r.estim_prec, 5, 1e-15, "uniroot downX estim.prec");
  // uniroot(function(x) 3 - x, c(5, 10), extendInt = "yes"): root 3, init.it 6
  r = uniroot([](double x) { return 3 - x; }, 5, 10, yes);
  check::ok(r.root == 3 && r.init_it == 6, "uniroot yes extends the lower end");
  // errors: try(uniroot(function(x) x^2+1, c(1, 2))) "f() values at end points not of opposite sign"
  check::throws([] { uniroot([](double x) { return x * x + 1; }, 1, 2); }, "uniroot same sign");
  // uniroot(function(x) 1, 0:1, extendInt = "yes"): "no sign change found in 1000 iterations"
  check::throws([&] { uniroot([](double) { return 1.0; }, 0, 1, yes); }, "uniroot no sign change");
  // uniroot(function(x) x - 5, c(0, 1), extendInt = "downX"): "no sign change found in 1000 iterations"
  check::throws([&] { uniroot([](double x) { return x - 5; }, 0, 1, downx); }, "uniroot downX wrong direction");
  check::throws([] { uniroot([](double x) { return x; }, 1, 1); }, "uniroot lower == upper");
  RootOptions bad;
  bad.extend_int = "sideways";
  check::throws([&] { uniroot([](double x) { return x; }, -1, 1, bad); }, "uniroot bad extendInt");
}

void check_optim(const std::string& what, const OptimResult& r, const Vec& par, double value, double par_tol,
                 double value_tol) {
  check::ok(r.par.size() == par.size(), what + " par size");
  for (size_t i = 0; i < par.size() && i < r.par.size(); ++i)
    check::close(r.par[i], par[i], par_tol, what + " par[" + std::to_string(i) + "]");
  check::close(r.value, value, value_tol, what + " value");
}

void test_optim() {
  OptimOptions nm;
  // optim(c(-1.2,1), fr): par 1.0002601387256695 1.000505999303765, value 8.8252410967227472e-08, counts 195 NA
  auto r = optim(fr, {-1.2, 1}, nm);
  check_optim("optim NM rosenbrock", r, {1.0002601387256695, 1.000505999303765}, 8.8252410967227472e-08, 1e-12,
              1e-14);
  check::ok(r.fn_count == 195 && r.gr_count == -1, "optim NM rosenbrock counts 195 NA");
  check::ok(r.convergence == 0, "optim NM convergence 0");
  // optim(c(-1.2,1), fr, control = list(maxit = 50)): value 3.1037095813242419, counts 51, convergence 1
  OptimOptions nm50;
  nm50.maxit = 50;
  r = optim(fr, {-1.2, 1}, nm50);
  check_optim("optim NM maxit 50", r, {-0.71460362434386515, 0.47018069267272222}, 3.1037095813242419, 1e-12, 1e-14);
  check::ok(r.fn_count == 51 && r.convergence == 1, "optim NM maxit 50 counts 51, convergence 1");
  // optim(1, function(x) (x-2)^2): 1.9000000000000008, counts 10 (R warns that 1-D Nelder-Mead is unreliable)
  r = optim([](const Vec& x) { return (x[0] - 2) * (x[0] - 2); }, {1}, nm);
  check_optim("optim NM 1-D", r, {1.9000000000000008}, 0.0099999999999998406, 1e-14, 1e-14);
  check::ok(r.fn_count == 10, "optim NM 1-D counts 10");
  // optim(c(1,2,3), function(x) sum((x-c(3,-1,0.5))^2) + x[1]*x[2]): value -7.3333331013270868, counts 122
  auto f3 = [](const Vec& x) {
    const double c[3] = {3, -1, 0.5};
    double s = 0;
    for (int i = 0; i < 3; ++i) s += (x[i] - c[i]) * (x[i] - c[i]);
    return s + x[0] * x[1];
  };
  r = optim(f3, {1, 2, 3}, nm);
  check_optim("optim NM 3-D", r, {4.666313381459533, -3.3328330048434491, 0.50018337356677944},
              -7.3333331013270868, 1e-12, 1e-14);
  check::ok(r.fn_count == 122, "optim NM 3-D counts 122");

  OptimOptions bfgs;
  bfgs.method = "BFGS";
  // optim(c(-1.2,1), fr, grr, method = "BFGS"): 0.99999999690491426 0.99999999379741933, 9.59e-18, counts 110 43
  r = optim(fr, {-1.2, 1}, bfgs, grr);
  check_optim("optim BFGS gr", r, {0.99999999690491426, 0.99999999379741933}, 9.5949545448136196e-18, 1e-12, 1e-14);
  check::ok(r.fn_count == 110 && r.gr_count == 43, "optim BFGS gr counts 110 43");
  // optim(c(-1.2,1), fr, method = "BFGS"): 0.99980443323134416 0.99960838062338098, counts 120 38
  r = optim(fr, {-1.2, 1}, bfgs);
  check_optim("optim BFGS numgrad", r, {0.99980443323134416, 0.99960838062338098}, 3.8273827581259985e-08, 1e-10,
              1e-14);
  check::ok(std::abs(r.fn_count - 120) <= 2 && r.gr_count == 38, "optim BFGS numgrad counts (R 120 38)");
  // optim(c(-1.2,1), fr, grr, method = "BFGS", control = list(maxit = 5)): value 2.8306330601569445, counts 16 5
  OptimOptions bfgs5 = bfgs;
  bfgs5.maxit = 5;
  r = optim(fr, {-1.2, 1}, bfgs5, grr);
  check_optim("optim BFGS maxit 5", r, {-0.59874761361855311, 0.41090472814195395}, 2.8306330601569445, 1e-12, 1e-14);
  check::ok(r.fn_count == 16 && r.gr_count == 5 && r.convergence == 1, "optim BFGS maxit 5 counts 16 5, conv 1");
  // optim(c(1,2,3), f3, method = "BFGS"): value -7.3333333159884209, counts 10 8
  r = optim(f3, {1, 2, 3}, bfgs);
  check_optim("optim BFGS 3-D", r, {4.6666942473343198, -3.3334230388788084, 0.49989506538323181},
              -7.3333333159884209, 1e-10, 1e-14);
  check::ok(r.fn_count == 10 && r.gr_count == 8, "optim BFGS 3-D counts 10 8");

  OptimOptions cg;
  cg.method = "CG";
  // optim(c(-1.2,1), fr, grr, method = "CG"): value 3.1065787964981539, counts 402 101, convergence 1
  // (R's $par there, -0.7648373 0.5927588, is the base point of the last search, not the point with $value;
  // we report the point whose value is returned.)
  r = optim(fr, {-1.2, 1}, cg, grr);
  check::close(r.value, 3.1065787964981539, 1e-12, "optim CG FR value");
  check::close(fr(r.par), r.value, 1e-15, "optim CG FR par consistent with value");
  check::ok(r.fn_count == 402 && r.gr_count == 101 && r.convergence == 1, "optim CG FR counts 402 101, conv 1");
  // optim(c(-1.2,1), fr, method = "CG"): value 3.1064747173049079, counts 402 101
  r = optim(fr, {-1.2, 1}, cg);
  check::close(r.value, 3.1064747173049079, 1e-10, "optim CG numgrad value");
  check::ok(r.fn_count == 402 && r.gr_count == 101, "optim CG numgrad counts 402 101");
  // optim(c(-1.2,1), fr, grr, method = "CG", control = list(type = 3)): value 0.010044747270463623, counts 355 101
  OptimOptions cg3 = cg;
  cg3.cg_type = 3;
  r = optim(fr, {-1.2, 1}, cg3, grr);
  check::close(r.value, 0.010044747270463623, 1e-9, "optim CG Beale-Sorenson value");
  check::ok(r.fn_count == 355 && r.gr_count == 101, "optim CG Beale-Sorenson counts 355 101");
  // optim(c(-1.2,1), fr, grr, method = "CG", control = list(maxit = 5000)): value 2.2214464856014689e-07,
  // counts 19957 5001
  OptimOptions cg5000 = cg;
  cg5000.maxit = 5000;
  r = optim(fr, {-1.2, 1}, cg5000, grr);
  check::close(r.value, 2.2214464856014689e-07, 1e-6, "optim CG maxit 5000 value");
  check::ok(r.fn_count == 19957 && r.gr_count == 5001, "optim CG maxit 5000 counts 19957 5001");
  // optim(c(-1.2,1), fr, grr, method = "CG", control = list(type = 2)): counts 385 101, convergence 1
  OptimOptions cg2 = cg;
  cg2.cg_type = 2;
  r = optim(fr, {-1.2, 1}, cg2, grr);
  check::ok(r.gr_count == 101 && r.convergence == 1 && std::abs(r.fn_count - 385) <= 6,
            "optim CG Polak-Ribiere counts (R 385 101) " + std::to_string(r.fn_count));

  OptimOptions lb;
  lb.method = "L-BFGS-B";
  // optim(c(-1.2,1), fr, grr, method = "L-BFGS-B"): 0.99999971674932298 0.99999947177995274, counts 47 47
  r = optim(fr, {-1.2, 1}, lb, grr);
  check_optim("optim L-BFGS-B gr", r, {1, 1}, 0, 1e-5, 1e-10);
  check::ok(r.convergence == 0 && r.message == "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",
            "optim L-BFGS-B message");
  check::ok(r.fn_count == r.gr_count && r.fn_count > 30 && r.fn_count < 70, "optim L-BFGS-B counts (R 47 47)");
  // optim(c(-1.2,1), fr, method = "L-BFGS-B"): 0.99980003783382454 0.99960011579791486, value 3.998e-08, 49 49
  r = optim(fr, {-1.2, 1}, lb);
  check_optim("optim L-BFGS-B numgrad", r, {0.99980003783382454, 0.99960011579791486}, 3.9984867903696438e-08, 3e-4,
              1e-7);
  // bounded: optim(c(-1.2,1), fr, grr, method = "L-BFGS-B", lower = c(-2,-2), upper = c(0.5, 2)):
  // par 0.5 0.25, value 0.25, counts 31 31, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"
  OptimOptions lbb = lb;
  lbb.lower = {-2, -2};
  lbb.upper = {0.5, 2};
  r = optim(fr, {-1.2, 1}, lbb, grr);
  check_optim("optim L-BFGS-B bounded", r, {0.5, 0.25}, 0.25, 1e-10, 1e-12);
  check::ok(r.par[0] == 0.5, "optim L-BFGS-B active bound exactly at 0.5");
  check::ok(r.message == "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL", "optim L-BFGS-B bounded message");
  r = optim(fr, {-1.2, 1}, lbb);
  check_optim("optim L-BFGS-B bounded numgrad", r, {0.5, 0.25}, 0.25, 1e-10, 1e-12);
  // ?optim's 25-dimensional box example: flb <- function(x) { p <- length(x); sum(c(1, rep(4, p-1)) *
  // (x - c(1, x[-p])^2)^2) }; optim(rep(3, 25), flb, NULL, method = "L-BFGS-B", lower = rep(2, 25),
  // upper = rep(4, 25)): par 2 (x23) 2.1090930615370751 4, value 368.10591287434056
  auto flb = [](const Vec& x) {
    double s = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      double prev = i == 0 ? 1.0 : x[i - 1] * x[i - 1];
      s += (i == 0 ? 1.0 : 4.0) * (x[i] - prev) * (x[i] - prev);
    }
    return s;
  };
  OptimOptions lb25 = lb;
  lb25.lower = {2};
  lb25.upper = {4};
  r = optim(flb, Vec(25, 3.0), lb25);
  Vec flb_par(25, 2.0);
  flb_par[23] = 2.1090930615370751;
  flb_par[24] = 4;
  check_optim("optim L-BFGS-B flb 25", r, flb_par, 368.10591287434056, 1e-6, 1e-12);
  // optim(c(-1.2,1), fr, grr, method = "L-BFGS-B", control = list(maxit = 5)): convergence 1, message "NEW_X"
  OptimOptions lb5 = lb;
  lb5.maxit = 5;
  r = optim(fr, {-1.2, 1}, lb5, grr);
  check::ok(r.convergence == 1 && r.message == "NEW_X", "optim L-BFGS-B maxit 5");
  // bounds with Nelder-Mead switch to L-BFGS-B, as in R
  OptimOptions nmb;
  nmb.lower = {-2, -2};
  nmb.upper = {0.5, 2};
  r = optim(fr, {-1.2, 1}, nmb, grr);
  check::ok(r.par[0] == 0.5 && r.message.find("CONVERGENCE") == 0, "optim bounds switch to L-BFGS-B");

  // quadratic with Hessian [[4,1],[1,6]]: minimum at (8/23, -9/23), value -13/23
  const Vec qpar = {8.0 / 23, -9.0 / 23};
  const double qval = -13.0 / 23;
  const char* methods[] = {"Nelder-Mead", "BFGS", "CG", "L-BFGS-B"};
  for (const char* m : methods) {
    OptimOptions o;
    o.method = m;
    o.hessian = true;
    r = optim(quad, {1, 1}, o);
    std::string w = std::string("optim quad ") + m;
    // NM stops on function values: par is accurate to ~1e-4 only (R: 0.34782209680132437 -0.39139372014616802)
    check_optim(w, r, qpar, qval, std::string(m) == "Nelder-Mead" ? 3e-4 : 1e-6, 1e-7);
    check::ok(r.hessian.size() == 4, w + " hessian size");
    check::close(r.hessian[0], 4, 1e-8, w + " hessian[1,1]");
    check::close(r.hessian[1], 1, 1e-8, w + " hessian[1,2]");
    check::close(r.hessian[2], 1, 1e-8, w + " hessian[2,1]");
    check::close(r.hessian[3], 6, 1e-8, w + " hessian[2,2]");
  }
  {
    // optim(c(1,1), q, method = "Nelder-Mead") exact path: counts 69, value -0.56521736695366087
    r = optim(quad, {1, 1}, nm);
    check::ok(r.fn_count == 69, "optim quad NM counts 69");
    check::close(r.value, -0.56521736695366087, 1e-14, "optim quad NM value");
    // method "BFGS": counts 10 6; "CG": counts 31 13; "L-BFGS-B": counts 6 6
    r = optim(quad, {1, 1}, bfgs);
    check::ok(r.fn_count == 10 && r.gr_count == 6, "optim quad BFGS counts 10 6");
    r = optim(quad, {1, 1}, cg);
    check::ok(r.fn_count == 31 && r.gr_count == 13, "optim quad CG counts 31 13");
    check::close(r.value, -0.56521739130432502, 1e-14, "optim quad CG value");
    r = optim(quad, {1, 1}, lb);
    check::ok(r.fn_count == 6 && r.gr_count == 6, "optim quad L-BFGS-B counts 6 6");
  }
  // optim(c(-1.2,1), fr, hessian = TRUE)$hessian: 802.42281442667115 -400.10405549025830 199.99999999998931
  OptimOptions nmh;
  nmh.hessian = true;
  r = optim(fr, {-1.2, 1}, nmh);
  check::close(r.hessian[0], 802.42281442667115, 1e-10, "optim NM hessian[1,1]");
  check::close(r.hessian[1], -400.10405549025830, 1e-10, "optim NM hessian[1,2]");
  check::close(r.hessian[3], 199.99999999998931, 1e-10, "optim NM hessian[2,2]");
  // optim(c(-1.2,1), fr, grr, method = "BFGS", hessian = TRUE)$hessian: 802.00039505280961 -399.99999876196171 200
  OptimOptions bh = bfgs;
  bh.hessian = true;
  r = optim(fr, {-1.2, 1}, bh, grr);
  check::close(r.hessian[0], 802.00039505280961, 1e-10, "optim BFGS gr hessian[1,1]");
  check::close(r.hessian[1], -399.99999876196171, 1e-10, "optim BFGS gr hessian[1,2]");
  check::close(r.hessian[3], 200.00000000000017, 1e-10, "optim BFGS gr hessian[2,2]");
  // optimHess(c(1,1), fr): 802.0008 -400 200
  Vec H = numeric_hessian(fr, {1, 1});
  check::close(H[0], 802.0008, 1e-9, "optimHess [1,1]");
  check::close(H[1], -400, 1e-9, "optimHess [1,2]");
  check::close(H[2], -400, 1e-9, "optimHess [2,1]");
  check::close(H[3], 200, 1e-9, "optimHess [2,2]");

  // optim(1, function(x) (x-2)^2, method = "Brent", lower = 0, upper = 5): par 2, value 0, counts NA NA
  OptimOptions br;
  br.method = "Brent";
  br.lower = {0};
  br.upper = {5};
  r = optim([](const Vec& x) { return (x[0] - 2) * (x[0] - 2); }, {1}, br);
  check::close(r.par[0], 2, 1e-8, "optim Brent par");
  check::ok(r.fn_count == -1 && r.gr_count == -1, "optim Brent counts NA NA");

  // invalid input
  check::throws([] { optim(fr, {}, OptimOptions()); }, "optim empty par");
  OptimOptions unk;
  unk.method = "Newton";
  check::throws([&] { optim(fr, {1, 1}, unk); }, "optim unknown method");
  check::throws([&] { optim([](const Vec&) { return kNaN; }, {1, 1}, nm); }, "optim NaN at initial parameters");
  check::throws([&] { optim([](const Vec& x) { return x[0] < -0.5 ? kInf : x[0]; }, {0.0}, lb); },
                "optim L-BFGS-B non-finite fn");
  OptimOptions badb = lb;
  badb.lower = {1, 1};
  badb.upper = {0, 0};
  check::throws([&] { optim(fr, {0.5, 0.5}, badb); }, "optim lower > upper");
  check::throws([&] { optim(fr, {1, 1}, bfgs, [](const Vec&) { return Vec{1.0}; }); },
                "optim wrong gradient length");
  OptimOptions sann_norng;
  sann_norng.method = "SANN";
  check::throws([&] { optim(fr, {1, 1}, sann_norng); }, "optim SANN without Rng");
}

void test_sann() {
  // optim(c(-1.2,1), fr, method = "SANN"): counts 10000 NA; the value is random (R typically ~1e-3 .. 1e-2).
  Rng rng(20231);
  OptimOptions o;
  o.method = "SANN";
  o.rng = &rng;
  auto r = optim(fr, {-1.2, 1}, o);
  check::ok(r.fn_count == 10000 && r.gr_count == -1 && r.convergence == 0, "optim SANN counts 10000 NA");
  check::ok(r.value < 0.5 && r.value == fr(r.par), "optim SANN improves on fr(-1.2, 1) = 24.2");
  Rng rng2(20231);
  o.rng = &rng2;
  auto r2 = optim(fr, {-1.2, 1}, o);
  check::ok(r2.par == r.par && r2.value == r.value, "optim SANN deterministic for a seed");
  o.maxit = 100;
  r = optim(fr, {-1.2, 1}, o);
  check::ok(r.fn_count == 100, "optim SANN maxit 100 evaluations");
}

// A callback that throws on its 5th evaluation must unwind out of every entry point (run under ASan/valgrind
// to check that nothing leaks).
struct Abandoned {
  int code;
};
struct Thrower {
  std::shared_ptr<int> calls = std::make_shared<int>(0);
  void tick() const {
    if (++*calls == 5) throw Abandoned{5};
  }
};

template <class F>
void expect_abandoned(F&& f, const std::string& what) {
  bool got = false;
  try {
    f();
  } catch (const Abandoned& a) {
    got = a.code == 5;
  }
  check::ok(got, what + " propagates a callback exception");
}

void test_exceptions() {
  {
    Thrower t;
    expect_abandoned([&] { integrate([t](double x) { t.tick(); return x; }, 0, 1); }, "integrate");
  }
  {
    Thrower t;
    expect_abandoned([&] { integrate([t](double x) { t.tick(); return dnorm(x); }, -kInf, kInf); }, "integrate QAGI");
  }
  {
    Thrower t;
    expect_abandoned([&] { optimize([t](double x) { t.tick(); return std::sin(x); }, 0, 6); }, "optimize");
  }
  {
    Thrower t;
    expect_abandoned([&] { uniroot([t](double x) { t.tick(); return x * x * x - x - 1; }, 1, 2); }, "uniroot");
  }
  {
    Thrower t;
    RootOptions y;
    y.extend_int = "yes";
    expect_abandoned([&] { uniroot([t](double x) { t.tick(); return x - 50; }, 0, 1, y); }, "uniroot extendInt");
  }
  Rng rng(1);
  const char* methods[] = {"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"};
  for (const char* m : methods) {
    for (int analytic = 0; analytic < 2; ++analytic) {
      Thrower t;
      OptimOptions o;
      o.method = m;
      o.rng = &rng;
      std::function<Vec(const Vec&)> g;
      if (analytic) g = grr;
      expect_abandoned([&] { optim([t](const Vec& x) { t.tick(); return fr(x); }, {-1.2, 1}, o, g); },
                       std::string("optim ") + m + (analytic ? " (gr)" : " (numgrad)"));
    }
    // a gradient that throws
    Thrower t;
    OptimOptions o;
    o.method = m;
    o.rng = &rng;
    if (std::string(m) != "Nelder-Mead" && std::string(m) != "SANN")
      expect_abandoned([&] { optim(fr, {-1.2, 1}, o, [t](const Vec& x) { t.tick(); return grr(x); }); },
                       std::string("optim ") + m + " gradient");
  }
  {
    Thrower t;
    OptimOptions o;
    o.method = "Brent";
    o.lower = {0};
    o.upper = {5};
    expect_abandoned([&] { optim([t](const Vec& x) { t.tick(); return x[0] * x[0]; }, {1}, o); }, "optim Brent");
  }
  {
    // the Hessian evaluation after convergence
    int calls = 0;
    OptimOptions o;
    o.hessian = true;
    auto f = [&calls](const Vec& x) {
      ++calls;
      if (calls > 195) throw Abandoned{5};
      return fr(x);
    };
    expect_abandoned([&] { optim(f, {-1.2, 1}, o); }, "optim hessian");
  }
  {
    Thrower t;
    expect_abandoned([&] { numeric_hessian([t](const Vec& x) { t.tick(); return fr(x); }, {1, 1}); },
                     "numeric_hessian");
  }
}

}  // namespace

int main() {
  test_integrate();
  test_optimize();
  test_uniroot();
  test_optim();
  test_sann();
  test_exceptions();
  return check::report("optim");
}
