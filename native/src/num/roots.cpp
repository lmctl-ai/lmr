// One-dimensional minimisation (optimize) and root finding (uniroot).
// Both implement the procedures published in R. P. Brent, "Algorithms for Minimization without Derivatives",
// Prentice-Hall 1973: chapter 5 (localmin: golden section with successive parabolic interpolation) and chapter 4
// (zero: bisection, secant and inverse quadratic interpolation), with the tolerance conventions of the
// Forsythe-Malcolm-Moler versions (fmin: tol1 = sqrt(eps)|x| + tol/3; zeroin: tol1 = 2 eps |b| + tol/2).
// uniroot's interval extension follows R's documented extendInt semantics.
#include <cfloat>
#include <cmath>
#include <string>

#include "rcli/optim.hpp"

namespace rcli {

OptimizeResult optimize(const Fn1& f, double lower, double upper, bool maximum, double tol) {
  if (!f) throw Error("optimize: 'f' is not a function");
  if (!std::isfinite(lower) || !std::isfinite(upper)) throw Error("optimize: invalid 'xmin' or 'xmax' value");
  if (lower >= upper) throw Error("optimize: 'xmin' not less than 'xmax'");
  if (!(tol > 0.0)) throw Error("optimize: invalid 'tol' value");
  OptimizeResult out;
  auto raw = [&](double x) {
    ++out.evaluations;
    return f(x);
  };
  auto fn = [&](double x) {
    double v = raw(x);
    if (maximum) v = -v;
    return std::isfinite(v) ? v : DBL_MAX;  // as R: "NA/Inf replaced by maximum positive value"
  };
  const double c = 0.5 * (3.0 - std::sqrt(5.0));  // squared inverse of the golden ratio
  const double eps = std::sqrt(DBL_EPSILON);
  double a = lower, b = upper;
  double v = a + c * (b - a), w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = fn(x), fv = fx, fw = fx;
  const double tol3 = tol / 3.0;
  for (;;) {
    double xm = 0.5 * (a + b);
    double tol1 = eps * std::fabs(x) + tol3;
    double t2 = 2.0 * tol1;
    if (std::fabs(x - xm) <= t2 - 0.5 * (b - a)) break;  // stopping criterion
    bool golden = true;
    if (std::fabs(e) > tol1) {
      // fit a parabola through x, v, w
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      else q = -q;
      r = e;
      e = d;
      if (std::fabs(p) < std::fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
        // parabolic interpolation step
        d = p / q;
        double u = x + d;
        if (u - a < t2 || b - u < t2) d = (x < xm) ? tol1 : -tol1;  // f must not be evaluated too close to a, b
        golden = false;
      }
    }
    if (golden) {
      e = (x < xm) ? b - x : a - x;
      d = c * e;
    }
    // f must not be evaluated too close to x
    double u = x + (std::fabs(d) >= tol1 ? d : (d > 0.0 ? tol1 : -tol1));
    double fu = fn(u);
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  out.minimum = x;
  // R reports $objective as a fresh evaluation f(minimum); mirror it (and its evaluation count).
  out.objective = raw(x);
  return out;
}

namespace {

struct ZeroResult {
  double root, froot, estim_prec, c;
  int iter;  // function evaluations used; -1 if not converged
};

// Brent's zero: a, b bracket a root (fa, fb of opposite sign or zero); at most maxit evaluations of f.
template <class F>
ZeroResult zeroin(F& f, double a, double b, double fa, double fb, double tol, int maxit) {
  double c = a, fc = fa;
  double d = b - a, e = d;
  int used = 0;
  // the endpoints themselves may be roots
  if (fa == 0.0) return {a, fa, 0.0, a, 0};
  if (fb == 0.0) return {b, fb, 0.0, b, 0};
  for (;;) {
    if (std::fabs(fc) < std::fabs(fb)) {  // make b the best approximation
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    double tol1 = 2.0 * DBL_EPSILON * std::fabs(b) + 0.5 * tol;
    double m = 0.5 * (c - b);
    if (std::fabs(m) <= tol1 || fb == 0.0) return {b, fb, std::fabs(c - b), c, used};
    if (used >= maxit) return {b, fb, std::fabs(c - b), c, -1};
    if (std::fabs(e) < tol1 || std::fabs(fa) <= std::fabs(fb)) {
      d = e = m;  // bisection
    } else {
      double s = fb / fa, p, q;
      if (a == c) {  // linear interpolation
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {  // inverse quadratic interpolation
        double qq = fa / fc, r = fb / fc;
        p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
        q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      else p = -p;
      s = e;
      e = d;
      if (2.0 * p < 3.0 * m * q - std::fabs(tol1 * q) && p < std::fabs(0.5 * s * q)) {
        d = p / q;
      } else {
        d = e = m;
      }
    }
    a = b;
    fa = fb;
    b += (std::fabs(d) > tol1) ? d : (m > 0.0 ? tol1 : -tol1);
    fb = f(b);
    ++used;
    if ((fb > 0.0 && fc > 0.0) || (fb <= 0.0 && fc <= 0.0)) {
      c = a; fc = fa;
      d = e = b - a;
    }
  }
}

}  // namespace

RootResult uniroot(const Fn1& f, double lower, double upper, const RootOptions& opt) {
  if (!f) throw Error("uniroot: 'f' is not a function");
  if (!std::isfinite(lower) || !std::isfinite(upper)) throw Error("uniroot: lower and upper must be finite");
  if (!(lower < upper)) throw Error("uniroot: lower < upper  is not fulfilled");
  if (!(opt.tol > 0.0)) throw Error("uniroot: invalid 'tol' value");
  if (opt.maxiter < 1) throw Error("uniroot: invalid 'maxiter' value");
  int mode;  // 0 no, 1 yes, 2 upX, 3 downX
  if (opt.extend_int == "no") mode = 0;
  else if (opt.extend_int == "yes") mode = 1;
  else if (opt.extend_int == "upX") mode = 2;
  else if (opt.extend_int == "downX") mode = 3;
  else throw Error("uniroot: 'extendInt' must be one of \"no\", \"yes\", \"upX\", \"downX\"");

  auto fn = [&](double x) {
    double v = f(x);
    if (std::isnan(v)) throw Error("uniroot: f() value is NA");
    return v;
  };
  double flo = fn(lower), fup = fn(upper);
  RootResult out;
  int init_it = -1;
  if (mode == 0) {
    if (flo * fup > 0.0) throw Error("uniroot: f() values at end points not of opposite sign");
  } else {
    // Widen geometrically: steps start at 1% of max(1e-4, |end|) and double at each iteration.
    double dlo = 0.01 * std::fmax(1e-4, std::fabs(lower));
    double dup = 0.01 * std::fmax(1e-4, std::fabs(upper));
    double S = (mode == 3) ? -1.0 : 1.0;
    int it = 0;
    auto needs = [&](bool& lo, bool& up) {
      if (mode == 1) {
        lo = up = flo * fup > 0.0;
      } else {
        lo = S * flo > 0.0;
        up = S * fup < 0.0;
      }
      return lo || up;
    };
    bool lo = false, up = false;
    while (needs(lo, up)) {
      if (++it > opt.maxiter)
        throw Error("uniroot: no sign change found in " + std::to_string(opt.maxiter) + " iterations");
      if (lo) { lower -= dlo; flo = fn(lower); dlo *= 2.0; }
      if (up) { upper += dup; fup = fn(upper); dup *= 2.0; }
      if (!std::isfinite(lower) || !std::isfinite(upper))
        throw Error("uniroot: no sign change found: the search interval became infinite");
    }
    if (it > 0) init_it = it;
  }
  int remaining = opt.maxiter - (init_it > 0 ? init_it : 0);
  int calls = 0;
  auto counted = [&](double x) { ++calls; return fn(x); };
  ZeroResult z = zeroin(counted, lower, upper, flo, fup, opt.tol, remaining);
  out.root = z.root;
  out.f_root = fn(z.root);  // R re-evaluates f at the root for $f.root
  out.init_it = init_it;
  out.iter = (z.iter < 0 ? calls : z.iter) + (init_it > 0 ? init_it : 0);
  out.estim_prec = z.estim_prec;
  out.converged = z.iter >= 0;
  out.lower = std::fmin(z.root, z.c);
  out.upper = std::fmax(z.root, z.c);
  return out;
}

}  // namespace rcli
