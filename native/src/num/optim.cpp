// General-purpose optimisation (R's optim methods) and numerical Hessians.
// Sources (published algorithms; no R code consulted):
//   Nelder-Mead : J. C. Nash, "Compact Numerical Methods for Computers", 2nd ed. 1990, Algorithm 19
//                 (Nelder & Mead 1965 with Nash's choices: contraction after any reflection that does not improve on
//                 the best point, shrink only when the contraction fails and the reflection was worst).
//   BFGS        : Nash 1990, Algorithm 21 (variable metric method with backtracking acceptable-point search,
//                 inverse-Hessian BFGS update and periodic resets).
//   CG          : Nash 1990, Algorithm 22 (Fletcher-Reeves / Polak-Ribiere / Beale-Sorenson conjugate gradients
//                 with backtracking and a quadratic interpolation step).
//   L-BFGS-B    : R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited memory algorithm for bound constrained
//                 optimization", SIAM J. Sci. Comput. 16 (1995) 1190-1208: compact limited-memory BFGS matrix,
//                 generalised Cauchy point along the projected gradient path, subspace minimisation (direct primal
//                 method with backtracking into the box) and a line search enforcing the strong Wolfe conditions
//                 (Nocedal & Wright, "Numerical Optimization", 2006, Algorithms 3.5 and 3.6).
//   SANN        : C. J. P. Belisle, "Convergence theorems for a class of simulated annealing algorithms on R^d",
//                 J. Appl. Prob. 29 (1992) 885-895, with the schedule and kernel documented in R's ?optim.
//   Brent       : the one-dimensional optimize() of this library.
// Termination tests and convergence codes follow R's documentation (?optim).
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "rcli/optim.hpp"
#include "rcli/rng.hpp"

namespace rcli {
namespace {

constexpr double kBig = 1.0e35;  // Nash's stand-in for an uncomputable function value

using GradFn = std::function<Vec(const Vec&)>;

double dot(const Vec& a, const Vec& b) {
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Objective and gradient with R's counting conventions and bound-aware central differences.
struct Problem {
  const FnN& f;
  const GradFn& gr;
  size_t n;
  double ndeps;
  const Vec* lower = nullptr;  // finite-difference steps are kept inside these (L-BFGS-B)
  const Vec* upper = nullptr;
  int fn_count = 0, gr_count = 0;

  double value(const Vec& x) {
    ++fn_count;
    return f(x);
  }
  Vec grad(const Vec& x) {
    ++gr_count;
    return raw_grad(x);
  }
  Vec raw_grad(const Vec& x) const {
    if (gr) {
      Vec g = gr(x);
      if (g.size() != n)
        throw Error("optim: gradient in optim evaluated to length " + std::to_string(g.size()) + " not " +
                    std::to_string(n));
      return g;
    }
    Vec g(n), xx = x;
    for (size_t i = 0; i < n; ++i) {
      double hi = x[i] + ndeps, lo = x[i] - ndeps;
      if (upper && hi > (*upper)[i]) hi = (*upper)[i];
      if (lower && lo < (*lower)[i]) lo = (*lower)[i];
      xx[i] = hi;
      double f1 = f(xx);
      xx[i] = lo;
      double f2 = f(xx);
      xx[i] = x[i];
      double d = (f1 - f2) / (hi - lo);
      if (!std::isfinite(d)) throw Error("optim: non-finite finite-difference value [" + std::to_string(i + 1) + "]");
      g[i] = d;
    }
    return g;
  }
};

// ---------------------------------------------------------------------------------------------------------------
// Nelder-Mead (Nash Algorithm 19)
void nelder_mead(Problem& P, Vec& b, const OptimOptions& o, int maxit, OptimResult& out) {
  const size_t n = P.n;
  const double alpha = o.alpha, beta = o.beta, gamma = o.gamma;
  auto fval = [&](const Vec& x) {
    double v = P.value(x);
    return std::isfinite(v) ? v : kBig;
  };
  double f0 = P.value(b);
  if (!std::isfinite(f0)) throw Error("optim: function cannot be evaluated at initial parameters");
  if (maxit <= 0) {
    out.par = b;
    out.value = f0;
    out.fn_count = 0;
    out.convergence = 0;
    return;
  }
  const double convtol = o.reltol * (std::fabs(f0) + o.reltol);
  const size_t n1 = n + 1;
  // simplex vertices (columns) with their function values
  std::vector<Vec> p(n1, b);
  Vec fv(n1, 0.0);
  fv[0] = f0;
  double step = 0.0;
  for (size_t i = 0; i < n; ++i) step = std::max(step, 0.1 * std::fabs(b[i]));
  if (step == 0.0) step = 0.1;
  double size = 0.0;
  for (size_t j = 1; j < n1; ++j) {
    double trystep = step;
    while (p[j][j - 1] == b[j - 1]) {  // make sure the step changes the coordinate
      p[j][j - 1] = b[j - 1] + trystep;
      trystep *= 10.0;
    }
    size += trystep;
  }
  double oldsize = size;
  bool calcvl = true;
  size_t L = 0;
  int fail = 0;
  Vec c(n), r(n), t(n);
  for (;;) {
    if (calcvl) {
      for (size_t j = 0; j < n1; ++j)
        if (j != L) fv[j] = fval(p[j]);
      calcvl = false;
    }
    double VL = fv[L], VH = VL;
    size_t H = L;
    for (size_t j = 0; j < n1; ++j) {
      if (j == L) continue;
      if (fv[j] < VL) { L = j; VL = fv[j]; }
      if (fv[j] > VH) { H = j; VH = fv[j]; }
    }
    if (VH <= VL + convtol || VL <= o.abstol) break;
    // centroid of all vertices except the highest
    for (size_t i = 0; i < n; ++i) {
      double s = -p[H][i];
      for (size_t j = 0; j < n1; ++j) s += p[j][i];
      c[i] = s / static_cast<double>(n);
    }
    for (size_t i = 0; i < n; ++i) r[i] = (1.0 + alpha) * c[i] - alpha * p[H][i];
    double VR = fval(r);
    if (VR < VL) {
      // try an expansion
      for (size_t i = 0; i < n; ++i) t[i] = gamma * r[i] + (1.0 - gamma) * c[i];
      double fe = fval(t);
      if (fe < VR) { p[H] = t; fv[H] = fe; }
      else { p[H] = r; fv[H] = VR; }
    } else {
      if (VR < VH) { p[H] = r; fv[H] = VR; }  // reflection kept, then contract towards the centroid
      for (size_t i = 0; i < n; ++i) t[i] = (1.0 - beta) * p[H][i] + beta * c[i];
      double fc = fval(t);
      if (fc < fv[H]) {
        p[H] = t;
        fv[H] = fc;
      } else if (VR >= VH) {
        // shrink towards the lowest vertex
        calcvl = true;
        size = 0.0;
        for (size_t j = 0; j < n1; ++j) {
          if (j == L) continue;
          for (size_t i = 0; i < n; ++i) {
            p[j][i] = beta * (p[j][i] - p[L][i]) + p[L][i];
            size += std::fabs(p[j][i] - p[L][i]);
          }
        }
        if (size < oldsize) {
          oldsize = size;
        } else {
          fail = 10;  // the simplex has degenerated
          break;
        }
      }
    }
    if (P.fn_count > maxit) {  // tested after the step, as R does; L is then the best vertex before the step
      fail = 1;
      break;
    }
  }
  out.par = p[L];
  out.value = fv[L];
  out.convergence = fail;
}

// ---------------------------------------------------------------------------------------------------------------
// BFGS variable metric method (Nash Algorithm 21)
void bfgs(Problem& P, Vec& b, const OptimOptions& o, int maxit, OptimResult& out) {
  const size_t n = P.n;
  const double stepredn = 0.2, acctol = 0.0001, reltest = 10.0;
  double f = P.value(b);
  if (!std::isfinite(f)) throw Error("optim: initial value in 'vmmin' is not finite");
  if (maxit <= 0) {
    out.par = b;
    out.value = f;
    out.fn_count = 0;
    out.gr_count = 0;
    return;
  }
  double fmin = f;
  Vec xbest = b;
  Vec g = P.grad(b);
  int iter = 1;
  int gradcount = 1;
  int ilast = gradcount;
  std::vector<Vec> B(n, Vec(n, 0.0));
  Vec t(n), X(n), c(n);
  size_t count = 0;
  int fail = 1;
  do {
    if (ilast == gradcount) {
      for (size_t i = 0; i < n; ++i) {
        std::fill(B[i].begin(), B[i].end(), 0.0);
        B[i][i] = 1.0;
      }
    }
    X = b;
    c = g;
    double gradproj = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (size_t j = 0; j < n; ++j) s -= B[i][j] * g[j];
      t[i] = s;
      gradproj += s * g[i];
    }
    if (gradproj < 0.0) {  // downhill: backtrack to an acceptable point
      double steplength = 1.0;
      bool accpoint = false;
      do {
        count = 0;
        for (size_t i = 0; i < n; ++i) {
          b[i] = X[i] + steplength * t[i];
          if (reltest + X[i] == reltest + b[i]) ++count;
        }
        if (count < n) {
          f = P.value(b);
          accpoint = std::isfinite(f) && f <= fmin + gradproj * steplength * acctol;
          if (!accpoint) steplength *= stepredn;
        }
      } while (!(count == n || accpoint));
      if (accpoint) {
        xbest = b;
        // relative convergence test on the function value (R's reltol / abstol)
        bool enough = (f > o.abstol) && std::fabs(f - fmin) > o.reltol * (std::fabs(fmin) + o.reltol);
        if (!enough) {
          count = n;
          fmin = f;
        }
      }
      if (count < n) {  // update the inverse Hessian approximation
        fmin = f;
        g = P.grad(b);
        ++gradcount;
        ++iter;
        double D1 = 0.0;
        for (size_t i = 0; i < n; ++i) {
          t[i] *= steplength;
          c[i] = g[i] - c[i];
          D1 += t[i] * c[i];
        }
        if (D1 > 0.0) {
          double D2 = 0.0;
          for (size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (size_t j = 0; j < n; ++j) s += B[i][j] * c[j];
            X[i] = s;
            D2 += s * c[i];
          }
          D2 = 1.0 + D2 / D1;
          for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) B[i][j] += (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
        } else {
          ilast = gradcount;  // the update would lose positive definiteness: reset
        }
      } else {
        b = xbest;  // the failed search leaves b at its last trial point
        if (ilast < gradcount) {  // no progress: retry from steepest descent
          count = 0;
          ilast = gradcount;
        }
      }
    } else {  // uphill search direction
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
    }
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * static_cast<int>(n)) ilast = gradcount;  // periodic restart
  } while (count != n || ilast != gradcount);
  if (iter < maxit) fail = 0;
  out.par = xbest;
  out.value = fmin;
  out.convergence = fail;
}

// ---------------------------------------------------------------------------------------------------------------
// Conjugate gradients (Nash Algorithm 22)
void conjugate_gradients(Problem& P, Vec& b, const OptimOptions& o, int maxit, OptimResult& out) {
  const size_t n = P.n;
  const double setstep = 1.7, acctol = 0.0001, stepredn = 0.2, reltest = 10.0;
  const int type = o.cg_type;
  if (type < 1 || type > 3) throw Error("optim: unknown 'type' in \"CG\" method");
  const double tol = o.reltol * static_cast<double>(n) * std::sqrt(o.reltol);
  double f = P.value(b);
  if (!std::isfinite(f)) throw Error("optim: function cannot be evaluated at initial parameters");
  if (maxit <= 0) {
    out.par = b;
    out.value = f;
    out.fn_count = 0;
    out.gr_count = 0;
    return;
  }
  double fmin = f;
  Vec xbest = b;
  Vec c(n, 0.0), t(n), X(n), g(n);
  int gradcount = 0;
  int fail = 0;
  size_t count = 0;
  int cycle = 0;
  double G1 = 0.0;
  do {
    std::fill(t.begin(), t.end(), 0.0);
    std::fill(c.begin(), c.end(), 0.0);
    cycle = 0;
    double oldstep = 1.0;
    count = 0;
    do {
      ++cycle;
      ++gradcount;
      if (gradcount > maxit) {
        fail = 1;
        break;
      }
      g = P.grad(b);
      G1 = 0.0;
      double G2 = 0.0;
      for (size_t i = 0; i < n; ++i) {
        X[i] = b[i];
        switch (type) {
          case 1:  // Fletcher-Reeves
            G1 += g[i] * g[i];
            G2 += c[i] * c[i];
            break;
          case 2:  // Polak-Ribiere
            G1 += g[i] * (g[i] - c[i]);
            G2 += c[i] * c[i];
            break;
          default:  // Beale-Sorenson
            G1 += g[i] * (g[i] - c[i]);
            G2 += t[i] * (g[i] - c[i]);
            break;
        }
        c[i] = g[i];
      }
      double steplength = oldstep;
      if (G1 > tol) {
        double G3 = (G2 > 0.0) ? G1 / G2 : 1.0;
        double gradproj = 0.0;
        for (size_t i = 0; i < n; ++i) {
          t[i] = t[i] * G3 - g[i];
          gradproj += t[i] * g[i];
        }
        bool accpoint = false;
        do {
          count = 0;
          for (size_t i = 0; i < n; ++i) {
            b[i] = X[i] + steplength * t[i];
            if (reltest + X[i] == reltest + b[i]) ++count;
          }
          if (count < n) {
            f = P.value(b);
            accpoint = std::isfinite(f) && f <= fmin + gradproj * steplength * acctol;
            if (!accpoint) steplength *= stepredn;
            else fmin = f;  // as observed in R: the accepted value becomes the reference before interpolating
          }
        } while (!(count == n || accpoint));
        if (count < n) {
          // "quadratic interpolation" step (with the reference value already updated this proposes half the
          // accepted step, which is what R's CG evaluates)
          double newstep = 2.0 * (f - fmin - gradproj * steplength);
          fmin = f;
          xbest = b;
          if (newstep > 0.0) {
            newstep = -(gradproj * steplength * steplength / newstep);
            for (size_t i = 0; i < n; ++i) b[i] = X[i] + newstep * t[i];
            f = P.value(b);
            if (std::isfinite(f) && f < fmin) {
              fmin = f;
              xbest = b;
            } else {
              b = xbest;
            }
          }
        } else {
          b = xbest;
        }
      }
      oldstep = setstep * steplength;
      if (oldstep > 1.0) oldstep = 1.0;
    } while (count != n && G1 > tol && cycle != static_cast<int>(n));
    if (fail) break;
  } while (!(cycle == 1 && (count == n || G1 <= tol)));
  out.par = xbest;
  out.value = fmin;
  out.convergence = fail;
}

// ---------------------------------------------------------------------------------------------------------------
// Small dense linear algebra for L-BFGS-B (row-major square matrices).

// Solve A x = rhs by Gaussian elimination with partial pivoting; false if A is (numerically) singular.
bool solve_dense(std::vector<double> A, Vec& rhs, size_t m) {
  for (size_t k = 0; k < m; ++k) {
    size_t piv = k;
    double best = std::fabs(A[k * m + k]);
    for (size_t i = k + 1; i < m; ++i)
      if (std::fabs(A[i * m + k]) > best) { best = std::fabs(A[i * m + k]); piv = i; }
    if (!(best > 0.0) || !std::isfinite(best)) return false;
    if (piv != k) {
      for (size_t j = 0; j < m; ++j) std::swap(A[k * m + j], A[piv * m + j]);
      std::swap(rhs[k], rhs[piv]);
    }
    for (size_t i = k + 1; i < m; ++i) {
      double l = A[i * m + k] / A[k * m + k];
      if (l == 0.0) continue;
      for (size_t j = k; j < m; ++j) A[i * m + j] -= l * A[k * m + j];
      rhs[i] -= l * rhs[k];
    }
  }
  for (size_t k = m; k-- > 0;) {
    double s = rhs[k];
    for (size_t j = k + 1; j < m; ++j) s -= A[k * m + j] * rhs[j];
    rhs[k] = s / A[k * m + k];
  }
  for (double v : rhs)
    if (!std::isfinite(v)) return false;
  return true;
}

// Compact limited-memory BFGS matrix B = theta I - W M W' with W = [Y, theta S] (Byrd, Nocedal & Schnabel 1994).
struct LbfgsMatrix {
  size_t n = 0, m = 0;
  std::vector<Vec> S, Y;  // oldest first
  double theta = 1.0;
  std::vector<double> Minv;  // 2k x 2k middle matrix inverse [[-D, L'], [L, theta S'S]]
  std::vector<double> M;     // its inverse

  size_t k() const { return S.size(); }
  void reset() {
    S.clear();
    Y.clear();
    theta = 1.0;
    Minv.clear();
    M.clear();
  }
  // Returns false (and resets) if the middle matrix is singular.
  bool rebuild() {
    size_t kk = k(), m2 = 2 * kk;
    Minv.assign(m2 * m2, 0.0);
    for (size_t i = 0; i < kk; ++i) {
      Minv[i * m2 + i] = -dot(S[i], Y[i]);
      for (size_t j = 0; j < kk; ++j) {
        if (i > j) {
          double lij = dot(S[i], Y[j]);
          Minv[(kk + i) * m2 + j] = lij;  // L
          Minv[j * m2 + kk + i] = lij;    // L'
        }
        Minv[(kk + i) * m2 + kk + j] = theta * dot(S[i], S[j]);
      }
    }
    M.assign(m2 * m2, 0.0);
    for (size_t col = 0; col < m2; ++col) {
      Vec e(m2, 0.0);
      e[col] = 1.0;
      if (!solve_dense(Minv, e, m2)) {
        reset();
        return false;
      }
      for (size_t row = 0; row < m2; ++row) M[row * m2 + col] = e[row];
    }
    return true;
  }
  void add(const Vec& s, const Vec& y) {
    if (S.size() == m) {
      S.erase(S.begin());
      Y.erase(Y.begin());
    }
    S.push_back(s);
    Y.push_back(y);
    theta = dot(y, y) / dot(s, y);
    rebuild();
  }
  // W' v (length 2k)
  Vec Wt(const Vec& v) const {
    size_t kk = k();
    Vec w(2 * kk);
    for (size_t i = 0; i < kk; ++i) {
      w[i] = dot(Y[i], v);
      w[kk + i] = theta * dot(S[i], v);
    }
    return w;
  }
  // row i of W
  double W(size_t row, size_t col) const {
    size_t kk = k();
    return col < kk ? Y[col][row] : theta * S[col - kk][row];
  }
  Vec Mtimes(const Vec& v) const {
    size_t m2 = v.size();
    Vec r(m2, 0.0);
    for (size_t i = 0; i < m2; ++i)
      for (size_t j = 0; j < m2; ++j) r[i] += M[i * m2 + j] * v[j];
    return r;
  }
  Vec times(const Vec& v) const {
    Vec r(v.size());
    for (size_t i = 0; i < v.size(); ++i) r[i] = theta * v[i];
    size_t kk = k();
    if (kk == 0) return r;
    Vec u = Mtimes(Wt(v));
    for (size_t i = 0; i < v.size(); ++i) {
      double s = 0.0;
      for (size_t j = 0; j < 2 * kk; ++j) s += W(i, j) * u[j];
      r[i] -= s;
    }
    return r;
  }
};

double proj_grad_norm(const Vec& x, const Vec& g, const Vec& l, const Vec& u) {
  double nrm = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    double gi = g[i];
    if (gi < 0.0) gi = std::max(x[i] - u[i], gi);
    else gi = std::min(x[i] - l[i], gi);
    nrm = std::max(nrm, std::fabs(gi));
  }
  return nrm;
}

// Generalised Cauchy point: first local minimiser of the quadratic model along the projected steepest-descent
// path x(t) = P(x - t g) (Byrd et al. 1995, section 4). Returns the Cauchy point.
Vec cauchy_point(const Vec& x, const Vec& g, const Vec& l, const Vec& u, const LbfgsMatrix& B) {
  const size_t n = x.size();
  Vec tb(n), d(n, 0.0);
  std::vector<size_t> order;
  for (size_t i = 0; i < n; ++i) {
    double t = kInf;
    if (g[i] < 0.0) t = std::isfinite(u[i]) ? (x[i] - u[i]) / g[i] : kInf;
    else if (g[i] > 0.0) t = std::isfinite(l[i]) ? (x[i] - l[i]) / g[i] : kInf;
    tb[i] = t;
    if (t > 0.0) d[i] = -g[i];
    if (t > 0.0 && std::isfinite(t)) order.push_back(i);
  }
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return tb[a] < tb[b]; });
  Vec xc = x, z(n, 0.0);
  double tcur = 0.0;
  size_t next = 0;
  for (;;) {
    bool any = false;
    for (double v : d)
      if (v != 0.0) { any = true; break; }
    if (!any) break;
    Vec Bz = B.times(z), Bd = B.times(d);
    double fp = dot(g, d) + dot(d, Bz);
    double fpp = dot(d, Bd);
    if (fp >= 0.0) break;  // the model increases along this segment: stop at its start
    double dtmin = fpp > 0.0 ? -fp / fpp : kInf;
    double tnext = next < order.size() ? tb[order[next]] : kInf;
    if (dtmin < tnext - tcur) {
      for (size_t i = 0; i < n; ++i) z[i] += dtmin * d[i];
      break;
    }
    if (!std::isfinite(tnext)) break;  // unbounded model along the path (cannot happen with B positive definite)
    double dt = tnext - tcur;
    for (size_t i = 0; i < n; ++i) z[i] += dt * d[i];
    tcur = tnext;
    // the variables whose breakpoint is reached become fixed at their bounds
    while (next < order.size() && tb[order[next]] <= tcur) {
      size_t i = order[next++];
      z[i] = (d[i] > 0.0 ? u[i] : l[i]) - x[i];
      d[i] = 0.0;
    }
  }
  for (size_t i = 0; i < n; ++i) xc[i] = std::min(std::max(x[i] + z[i], l[i]), u[i]);
  for (size_t k = 0; k < next; ++k) {
    size_t i = order[k];
    xc[i] = (-g[i] > 0.0) ? u[i] : l[i];
  }
  return xc;
}

// Subspace minimisation over the variables that are free at the Cauchy point, by the direct primal method:
// the unconstrained minimiser of the model in the subspace, truncated to stay in the box (section 5.1).
Vec subspace_min(const Vec& x, const Vec& g, const Vec& xc, const Vec& l, const Vec& u, const LbfgsMatrix& B) {
  const size_t n = x.size();
  std::vector<size_t> F;
  for (size_t i = 0; i < n; ++i)
    if (xc[i] > l[i] && xc[i] < u[i]) F.push_back(i);
  if (F.empty()) return xc;
  Vec diff(n);
  for (size_t i = 0; i < n; ++i) diff[i] = xc[i] - x[i];
  Vec r = B.times(diff);
  for (size_t i = 0; i < n; ++i) r[i] += g[i];
  const size_t nf = F.size();
  const double th = B.theta;
  Vec du(nf);
  for (size_t a = 0; a < nf; ++a) du[a] = -r[F[a]] / th;
  size_t kk = B.k();
  if (kk > 0) {
    size_t m2 = 2 * kk;
    // N = M^{-1} - (1/theta) W_F' W_F ;  v = N^{-1} W_F' r_F ; du = -(1/theta) r_F - (1/theta^2) W_F v
    std::vector<double> N = B.Minv;
    Vec v(m2, 0.0);
    for (size_t p = 0; p < m2; ++p) {
      for (size_t q = 0; q < m2; ++q) {
        double s = 0.0;
        for (size_t a = 0; a < nf; ++a) s += B.W(F[a], p) * B.W(F[a], q);
        N[p * m2 + q] -= s / th;
      }
      double s = 0.0;
      for (size_t a = 0; a < nf; ++a) s += B.W(F[a], p) * r[F[a]];
      v[p] = s;
    }
    if (!solve_dense(N, v, m2)) return xc;
    for (size_t a = 0; a < nf; ++a) {
      double s = 0.0;
      for (size_t p = 0; p < m2; ++p) s += B.W(F[a], p) * v[p];
      du[a] -= s / (th * th);
    }
  }
  // backtrack towards the Cauchy point to stay feasible
  double alpha = 1.0;
  for (size_t a = 0; a < nf; ++a) {
    size_t i = F[a];
    if (du[a] > 0.0 && std::isfinite(u[i])) alpha = std::min(alpha, (u[i] - xc[i]) / du[a]);
    else if (du[a] < 0.0 && std::isfinite(l[i])) alpha = std::min(alpha, (l[i] - xc[i]) / du[a]);
  }
  alpha = std::max(alpha, 0.0);
  Vec xb = xc;
  for (size_t a = 0; a < nf; ++a) {
    size_t i = F[a];
    xb[i] = std::min(std::max(xc[i] + alpha * du[a], l[i]), u[i]);
  }
  return xb;
}

struct LinePoint {
  double a, f, dphi;
  Vec x, g;
};

// Line search for the strong Wolfe conditions (Nocedal & Wright 2006, Algorithms 3.5 and 3.6) with cubic
// interpolation, step bounded by stpmax. Returns false if no acceptable step was found within maxls evaluations.
template <class Eval>
bool wolfe_search(Eval& eval, const Vec& x, double f0, const Vec& d, double dphi0, double a1, double stpmax,
                  LinePoint& out) {
  const double c1 = 1e-3, c2 = 0.9;
  const int maxls = 20;
  int nls = 0;
  auto at = [&](double a) {
    LinePoint p;
    p.a = a;
    p.x.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) p.x[i] = x[i] + a * d[i];
    eval(p.x, p.f, p.g);
    p.dphi = dot(p.g, d);
    ++nls;
    return p;
  };
  auto sufficient = [&](const LinePoint& p) { return p.f <= f0 + c1 * p.a * dphi0; };
  auto curvature = [&](const LinePoint& p) { return std::fabs(p.dphi) <= -c2 * dphi0; };
  // minimiser of the cubic interpolating (a, fa, da) and (b, fb, db), safeguarded into the inner 80% of [a, b]
  auto interpolate = [](const LinePoint& p, const LinePoint& q) {
    double lo = std::min(p.a, q.a), hi = std::max(p.a, q.a), w = hi - lo;
    double d1 = p.dphi + q.dphi - 3.0 * (p.f - q.f) / (p.a - q.a);
    double disc = d1 * d1 - p.dphi * q.dphi;
    double t = 0.5 * (lo + hi);
    if (disc >= 0.0) {
      double d2 = std::copysign(std::sqrt(disc), q.a - p.a);
      double denom = q.dphi - p.dphi + 2.0 * d2;
      if (denom != 0.0) {
        double c = q.a - (q.a - p.a) * (q.dphi + d2 - d1) / denom;
        if (std::isfinite(c)) t = c;
      }
    }
    return std::min(std::max(t, lo + 0.1 * w), hi - 0.1 * w);
  };
  LinePoint prev{0.0, f0, dphi0, x, {}};
  bool prev_is_origin = true;
  double a = std::min(a1, stpmax);
  LinePoint lo, hi;
  bool zoom = false;
  for (;;) {
    LinePoint cur = at(a);
    if (!sufficient(cur) || (!prev_is_origin && cur.f >= prev.f)) {
      lo = prev; hi = cur; zoom = true;
      break;
    }
    if (curvature(cur)) { out = std::move(cur); return true; }
    if (cur.dphi >= 0.0) {
      lo = cur; hi = prev; zoom = true;
      break;
    }
    if (a >= stpmax || nls >= maxls) {  // cannot extrapolate further: accept a point with sufficient decrease
      out = std::move(cur);
      return true;
    }
    double anext = std::min(stpmax, a + 4.0 * (a - prev.a));
    prev = std::move(cur);
    prev_is_origin = false;
    a = anext;
  }
  if (!zoom) return false;
  // lo always satisfies sufficient decrease and has the lower value; hi brackets a Wolfe point with it.
  while (nls < maxls) {
    double aj = interpolate(lo, hi);
    if (!(std::fabs(hi.a - lo.a) > 1e-16 * std::max(1.0, std::fabs(lo.a)))) break;
    LinePoint cur = at(aj);
    if (!sufficient(cur) || cur.f >= lo.f) {
      hi = std::move(cur);
    } else {
      if (curvature(cur)) { out = std::move(cur); return true; }
      if (cur.dphi * (hi.a - lo.a) >= 0.0) hi = lo;
      lo = std::move(cur);
    }
  }
  if (lo.a > 0.0 && lo.f < f0) {  // best effort: a strict decrease satisfying the Armijo condition
    out = std::move(lo);
    return true;
  }
  return false;
}

void lbfgsb(Problem& P, Vec& x, const Vec& l, const Vec& u, const OptimOptions& o, int maxit, OptimResult& out) {
  const size_t n = P.n;
  if (o.lmm < 1) throw Error("optim: 'lmm' must be positive");
  if (!(o.factr >= 0.0) || !(o.pgtol >= 0.0)) throw Error("optim: invalid 'factr' or 'pgtol'");
  for (size_t i = 0; i < n; ++i) x[i] = std::min(std::max(x[i], l[i]), u[i]);
  bool constrained = false;
  for (size_t i = 0; i < n; ++i)
    if (std::isfinite(l[i]) || std::isfinite(u[i])) constrained = true;
  auto eval = [&](const Vec& xx, double& fv, Vec& gv) {
    fv = P.value(xx);
    if (!std::isfinite(fv)) throw Error("optim: L-BFGS-B needs finite values of 'fn'");
    gv = P.grad(xx);
  };
  double f;
  Vec g;
  eval(x, f, g);
  LbfgsMatrix B;
  B.n = n;
  B.m = static_cast<size_t>(o.lmm);
  const double tolf = o.factr * DBL_EPSILON;
  int iter = 0;
  out.convergence = 0;
  if (proj_grad_norm(x, g, l, u) <= o.pgtol) {
    out.message = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
  } else {
    for (;;) {
      if (iter >= maxit) {
        out.convergence = 1;
        out.message = "NEW_X";
        break;
      }
      Vec xc = cauchy_point(x, g, l, u, B);
      Vec xb = subspace_min(x, g, xc, l, u, B);
      Vec d(n);
      for (size_t i = 0; i < n; ++i) d[i] = xb[i] - x[i];
      double dphi0 = dot(g, d);
      LinePoint np;
      bool ok = false;
      if (dphi0 < 0.0) {
        double dnorm = std::sqrt(dot(d, d));
        double stpmax = 1e10;
        if (constrained) {
          if (iter == 0) {  // largest feasible step along d
            for (size_t i = 0; i < n; ++i) {
              if (d[i] > 0.0 && std::isfinite(u[i])) stpmax = std::min(stpmax, (u[i] - x[i]) / d[i]);
              else if (d[i] < 0.0 && std::isfinite(l[i])) stpmax = std::min(stpmax, (l[i] - x[i]) / d[i]);
            }
          } else {
            stpmax = 1.0;
          }
        }
        double a1 = (iter == 0 && B.k() == 0) ? std::min(1.0 / dnorm, stpmax) : std::min(1.0, stpmax);
        auto clamped_eval = [&](Vec& xx, double& fv, Vec& gv) {
          for (size_t i = 0; i < n; ++i) xx[i] = std::min(std::max(xx[i], l[i]), u[i]);
          eval(xx, fv, gv);
        };
        ok = wolfe_search(clamped_eval, x, f, d, dphi0, a1, stpmax, np);
      }
      if (!ok) {
        if (B.k() > 0) {  // discard the limited-memory matrix and restart from steepest descent
          B.reset();
          continue;
        }
        out.convergence = 52;
        out.message = "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH";
        break;
      }
      ++iter;
      Vec s(n), y(n);
      for (size_t i = 0; i < n; ++i) {
        s[i] = np.x[i] - x[i];
        y[i] = np.g[i] - g[i];
      }
      double fold = f;
      x = std::move(np.x);
      f = np.f;
      g = std::move(np.g);
      if (proj_grad_norm(x, g, l, u) <= o.pgtol) {
        out.message = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
        break;
      }
      if (fold - f <= tolf * std::max({std::fabs(fold), std::fabs(f), 1.0})) {
        out.message = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH";
        break;
      }
      double sy = dot(s, y);
      if (sy > DBL_EPSILON * dot(y, y)) B.add(s, y);  // keep the matrix positive definite
    }
  }
  out.par = x;
  out.value = f;
}

// ---------------------------------------------------------------------------------------------------------------
// Simulated annealing (Belisle 1992; schedule and kernel as documented in ?optim)
void sann(Problem& P, Vec& x, const OptimOptions& o, int maxit, OptimResult& out) {
  if (!o.rng) throw Error("optim: method \"SANN\" needs a random number generator");
  if (o.tmax < 1) throw Error("optim: 'tmax' is not a positive integer");
  const size_t n = P.n;
  const double E1 = 1.7182818284590452354;  // exp(1) - 1
  double y = P.value(x);
  if (!std::isfinite(y)) throw Error("optim: function cannot be evaluated at initial parameters");
  out.par = x;
  out.value = y;
  if (maxit <= 0) {
    out.fn_count = 0;
    return;
  }
  Vec best = x, cand(n);
  double ybest = y;
  int its = 1;
  while (its < maxit) {
    // temperature temp / log(((its - 1) %/% tmax) * tmax + e)
    double t = o.temp / std::log(static_cast<double>((its - 1) / o.tmax * o.tmax) + E1 + 1.0);
    for (int k = 1; k <= o.tmax && its < maxit; ++k) {
      for (size_t i = 0; i < n; ++i) cand[i] = x[i] + t * o.rng->norm();  // Gaussian kernel, scale = temperature
      double ytry = P.value(cand);
      if (!std::isfinite(ytry)) ytry = kBig;
      double dy = ytry - y;
      if (dy <= 0.0 || o.rng->unif() < std::exp(-dy / t)) {  // Metropolis acceptance
        x = cand;
        y = ytry;
        if (y <= ybest) {
          best = x;
          ybest = y;
        }
      }
      ++its;
    }
  }
  out.par = best;
  out.value = ybest;
}

}  // namespace

Vec numeric_hessian(const FnN& f, const Vec& x, double ndeps) {
  if (!f) throw Error("optimHess: 'fn' is not a function");
  if (!(ndeps > 0.0)) throw Error("optimHess: 'ndeps' must be positive");
  GradFn none;
  Problem P{f, none, x.size(), ndeps};
  const size_t n = x.size();
  Vec H(n * n);
  Vec xx = x;
  for (size_t i = 0; i < n; ++i) {
    xx[i] = x[i] + ndeps;
    Vec g1 = P.raw_grad(xx);
    xx[i] = x[i] - ndeps;
    Vec g2 = P.raw_grad(xx);
    xx[i] = x[i];
    for (size_t j = 0; j < n; ++j) H[i * n + j] = (g1[j] - g2[j]) / (2.0 * ndeps);
  }
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < i; ++j) {
      double s = 0.5 * (H[i * n + j] + H[j * n + i]);
      H[i * n + j] = H[j * n + i] = s;
    }
  return H;
}

OptimResult optim(const FnN& f, const Vec& start, const OptimOptions& opt, const GradFn& gr) {
  if (!f) throw Error("optim: 'fn' is not a function");
  const size_t n = start.size();
  if (n == 0) throw Error("optim: 'par' must have positive length");
  for (double v : start)
    if (!std::isfinite(v)) throw Error("optim: non-finite value supplied by 'par'");
  if (!(opt.ndeps > 0.0)) throw Error("optim: 'ndeps' must be positive");
  if (std::isnan(opt.reltol) || std::isnan(opt.abstol)) throw Error("optim: invalid 'reltol' or 'abstol'");
  auto expand = [&](const Vec& v, double fill, const char* what) {
    if (v.empty()) return Vec(n, fill);
    if (v.size() == 1) return Vec(n, v[0]);
    if (v.size() != n) throw Error(std::string("optim: '") + what + "' must have length 1 or length(par)");
    return v;
  };
  Vec lower = expand(opt.lower, -kInf, "lower"), upper = expand(opt.upper, kInf, "upper");
  for (size_t i = 0; i < n; ++i)
    if (std::isnan(lower[i]) || std::isnan(upper[i]) || lower[i] > upper[i])
      throw Error("optim: invalid bounds ('lower' must not exceed 'upper')");
  bool bounded = false;
  for (size_t i = 0; i < n; ++i)
    if (lower[i] > -kInf || upper[i] < kInf) bounded = true;

  std::string method = opt.method;
  // As in R: bounds are only honoured by L-BFGS-B (and Brent); other methods switch to L-BFGS-B.
  if (bounded && method != "L-BFGS-B" && method != "Brent") method = "L-BFGS-B";
  OptimOptions o = opt;
  OptimResult out;
  Vec x = start;

  if (method == "Brent") {
    if (n != 1) throw Error("optim: method = \"Brent\" is only available for one-dimensional optimization");
    if (!std::isfinite(lower[0]) || !std::isfinite(upper[0]))
      throw Error("optim: 'lower' and 'upper' must be finite values");
    OptimizeResult r = optimize([&](double v) { return f(Vec{v}); }, lower[0], upper[0], false, opt.reltol);
    out.par = Vec{r.minimum};
    out.value = r.objective;
    out.fn_count = -1;
    out.gr_count = -1;
    out.convergence = 0;
  } else {
    Problem P{f, gr, n, opt.ndeps};
    int maxit = opt.maxit;
    if (method == "Nelder-Mead") {
      if (maxit < 0) maxit = 500;
      nelder_mead(P, x, o, maxit, out);
      if (maxit > 0) out.fn_count = P.fn_count;
      out.gr_count = -1;
    } else if (method == "BFGS") {
      if (maxit < 0) maxit = 100;
      bfgs(P, x, o, maxit, out);
      if (maxit > 0) {
        out.fn_count = P.fn_count;
        out.gr_count = P.gr_count;
      }
    } else if (method == "CG") {
      if (maxit < 0) maxit = 100;
      conjugate_gradients(P, x, o, maxit, out);
      if (maxit > 0) {
        out.fn_count = P.fn_count;
        out.gr_count = P.gr_count + (out.convergence == 1 ? 1 : 0);  // R counts the refused gradient call
      }
    } else if (method == "L-BFGS-B") {
      if (maxit < 0) maxit = 100;
      P.lower = &lower;
      P.upper = &upper;
      lbfgsb(P, x, lower, upper, o, maxit, out);
      out.fn_count = P.fn_count;
      out.gr_count = P.gr_count;
    } else if (method == "SANN") {
      if (maxit < 0) maxit = 10000;
      sann(P, x, o, maxit, out);
      if (maxit > 0) out.fn_count = P.fn_count;
      out.gr_count = -1;
    } else {
      throw Error("optim: unknown 'method' \"" + opt.method + "\"");
    }
  }
  if (opt.hessian) {
    Problem P{f, gr, n, opt.ndeps};
    if (gr) {
      Vec H(n * n);
      Vec xx = out.par;
      for (size_t i = 0; i < n; ++i) {
        xx[i] = out.par[i] + opt.ndeps;
        Vec g1 = P.raw_grad(xx);
        xx[i] = out.par[i] - opt.ndeps;
        Vec g2 = P.raw_grad(xx);
        xx[i] = out.par[i];
        for (size_t j = 0; j < n; ++j) H[i * n + j] = (g1[j] - g2[j]) / (2.0 * opt.ndeps);
      }
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j) H[i * n + j] = H[j * n + i] = 0.5 * (H[i * n + j] + H[j * n + i]);
      out.hessian = std::move(H);
    } else {
      out.hessian = numeric_hessian(f, out.par, opt.ndeps);
    }
  }
  return out;
}

}  // namespace rcli
