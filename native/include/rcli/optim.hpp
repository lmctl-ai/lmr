// Optimisation, root finding and quadrature over caller-supplied functions (callbacks).
// The DSL host turns a script loop into these callbacks (reverse communication), so every algorithm here must be
// exception-safe: a callback may throw at any evaluation and the algorithm must unwind without leaking.
// Algorithms come from the published literature (Nelder & Mead 1965 with Nash's variant, Nash's BFGS and CG
// (Compact Numerical Methods, 1990), Byrd, Lu, Nocedal & Zhu 1995 for L-BFGS-B, Brent 1973 fmin and zeroin,
// QUADPACK QAGS/QAGI (public domain) for integration).
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

class Rng;

struct IntegrateOptions {
  int subdivisions = 100;
  double rel_tol = kEps25;
  double abs_tol = kEps25;
};
struct IntegrateResult {
  double value = kNaN;
  double abs_error = kNaN;
  int subdivisions = 0;
  int evaluations = 0;
  int ier = 0;            // 0 ok; 1 max subdivisions; 2 roundoff; 3 bad integrand behaviour; 4 no convergence; 5 divergent
  std::string message;    // "OK" or R's wording for ier ("maximum number of subdivisions reached", ...)
};
// lower/upper may be -Inf/+Inf (QAGI transformation). A non-finite f value is an error (throws), as in R.
IntegrateResult integrate(const Fn1& f, double lower, double upper, const IntegrateOptions& opt = {});

struct OptimizeResult {
  double minimum = kNaN;  // the argument (R's $minimum or $maximum)
  double objective = kNaN;
  int evaluations = 0;
};
OptimizeResult optimize(const Fn1& f, double lower, double upper, bool maximum = false, double tol = kEps25);

struct RootOptions {
  double tol = kEps25;
  int maxiter = 1000;
  std::string extend_int = "no";  // "no", "yes", "upX", "downX" (R uniroot extendInt)
};
struct RootResult {
  double root = kNaN, f_root = kNaN;
  int iter = 0;            // zeroin iterations (R $iter, includes extension steps)
  int init_it = -1;        // extension steps, -1 if none (R $init.it is NA then)
  double estim_prec = kNaN;
  double lower = kNaN, upper = kNaN;  // final bracket
  bool converged = true;   // false when maxiter was reached (R warns "_NOT_ converged in <maxiter> iterations")
};
RootResult uniroot(const Fn1& f, double lower, double upper, const RootOptions& opt = {});

struct OptimOptions {
  std::string method = "Nelder-Mead";  // "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"
  Vec lower, upper;          // L-BFGS-B only; empty = unbounded
  int maxit = -1;            // -1: R defaults (500 Nelder-Mead, 100 BFGS/CG/L-BFGS-B, 10000 SANN)
  double reltol = 1.490116119384765625e-8;  // sqrt(eps)
  double abstol = -kInf;
  double ndeps = 1e-3;       // finite-difference step for numerical gradients (central differences)
  double alpha = 1.0, beta = 0.5, gamma = 2.0;  // Nelder-Mead reflection, contraction, expansion
  int cg_type = 1;           // 1 Fletcher-Reeves, 2 Polak-Ribiere, 3 Beale-Sorenson
  int lmm = 5;
  double factr = 1e7, pgtol = 0.0;
  double temp = 10.0;
  int tmax = 10;
  bool hessian = false;      // numerically differentiated Hessian at the optimum (like optimHess)
  Rng* rng = nullptr;        // required for SANN
};
struct OptimResult {
  Vec par;
  double value = kNaN;
  int fn_count = 0, gr_count = 0;  // -1 where R reports NA (gradient for Nelder-Mead/SANN, both for Brent)
  int convergence = 0;       // 0 ok, 1 maxit reached, 10 Nelder-Mead degeneracy, 51/52 L-BFGS-B warning/error
  std::string message;
  Vec hessian;               // row-major n x n when requested
};
// gr may be empty: then gradients use central differences with step ndeps (counting fn evaluations, not gr_count
// increments beyond R's convention: gr_count counts gradient evaluations either way).
OptimResult optim(const FnN& f, const Vec& start, const OptimOptions& opt = {},
                  const std::function<Vec(const Vec&)>& gr = {});
// Numerical Hessian by central differences of the (numerical) gradient, R optimHess style.
Vec numeric_hessian(const FnN& f, const Vec& x, double ndeps = 1e-3);

}  // namespace rcli
