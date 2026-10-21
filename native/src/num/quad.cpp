// Numerical integration: QUADPACK QAGS (finite range) and QAGI (infinite range).
// Implemented from the published algorithms in R. Piessens, E. de Doncker-Kapenga, C. W. Ueberhuber and
// D. K. Kahaner, "QUADPACK: A Subroutine Package for Automatic Integration", Springer 1983 (public domain):
//   - 21-point Gauss-Kronrod rule (QK21) for finite subintervals,
//   - 15-point Gauss-Kronrod rule on the transformation x = bound + s(1-t)/t of (0,1] (QK15I) for infinite ranges,
//   - globally adaptive bisection with an error-ordered interval list (QPSRT) and
//   - Wynn's epsilon algorithm (QELG) for extrapolation.
// The error flags follow the book; ier values above 2 are shifted down by one as in R, whose messages are reused.
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <vector>

#include "rcli/optim.hpp"

namespace rcli {
namespace {

constexpr double kEpmach = DBL_EPSILON;
constexpr double kUflow = DBL_MIN;
constexpr double kOflow = DBL_MAX;

// Counts evaluations and rejects non-finite values (R: "non-finite function value").
struct Integrand {
  const Fn1& f;
  int evals = 0;
  double operator()(double x) {
    ++evals;
    double v = f(x);
    if (!std::isfinite(v)) throw Error("integrate: non-finite function value");
    return v;
  }
};

struct RuleResult {
  double result, abserr, resabs, resasc;
};

// Common error scaling of the Gauss-Kronrod rules (QUADPACK book, section 3.4).
double scale_error(double abserr, double resabs, double resasc) {
  if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(200.0 * abserr / resasc, 1.5));
  if (resabs > kUflow / (50.0 * kEpmach)) abserr = std::max(kEpmach * 50.0 * resabs, abserr);
  return abserr;
}

// QK21: 21-point Kronrod rule with the embedded 10-point Gauss rule on [a, b].
RuleResult qk21(Integrand& f, double a, double b) {
  static constexpr double xgk[11] = {
      0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
      0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
      0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
      0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
      0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
      0.000000000000000000000000000000000};
  static constexpr double wgk[11] = {
      0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
      0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
      0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
      0.123491976262065851077208063032617, 0.134709217311473325928054001771707,
      0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
      0.149445554002916905664936468389821};
  static constexpr double wg[5] = {
      0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
      0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
      0.295524224714752870173892994651338};
  double centr = 0.5 * (a + b), hlgth = 0.5 * (b - a), dhlgth = std::fabs(hlgth);
  double fv1[10], fv2[10];
  double resg = 0.0;
  double fc = f(centr);
  double resk = wgk[10] * fc;
  double resabs = std::fabs(resk);
  for (int j = 0; j < 5; ++j) {  // Gauss abscissae (odd Kronrod indices)
    int jtw = 2 * j + 1;
    double absc = hlgth * xgk[jtw];
    double f1 = f(centr - absc), f2 = f(centr + absc);
    fv1[jtw] = f1; fv2[jtw] = f2;
    double fsum = f1 + f2;
    resg += wg[j] * fsum;
    resk += wgk[jtw] * fsum;
    resabs += wgk[jtw] * (std::fabs(f1) + std::fabs(f2));
  }
  for (int j = 0; j < 5; ++j) {  // Kronrod-only abscissae
    int jtwm1 = 2 * j;
    double absc = hlgth * xgk[jtwm1];
    double f1 = f(centr - absc), f2 = f(centr + absc);
    fv1[jtwm1] = f1; fv2[jtwm1] = f2;
    double fsum = f1 + f2;
    resk += wgk[jtwm1] * fsum;
    resabs += wgk[jtwm1] * (std::fabs(f1) + std::fabs(f2));
  }
  double reskh = resk * 0.5;
  double resasc = wgk[10] * std::fabs(fc - reskh);
  for (int j = 0; j < 10; ++j) resasc += wgk[j] * (std::fabs(fv1[j] - reskh) + std::fabs(fv2[j] - reskh));
  RuleResult r;
  r.result = resk * hlgth;
  resabs *= dhlgth;
  resasc *= dhlgth;
  r.resabs = resabs;
  r.resasc = resasc;
  r.abserr = scale_error(std::fabs((resk - resg) * hlgth), resabs, resasc);
  return r;
}

// QK15I: 15-point Kronrod rule (7-point Gauss) for the transformed integrand on [a, b] within (0, 1].
// inf = 1: (bound, +Inf); inf = -1: (-Inf, bound); inf = 2: (-Inf, +Inf) with bound = 0.
RuleResult qk15i(Integrand& f, double boun, int inf, double a, double b) {
  static constexpr double xgk[8] = {
      0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
      0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
      0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
      0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
  static constexpr double wgk[8] = {
      0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
      0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
      0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
      0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
  static constexpr double wg[8] = {
      0.0, 0.129484966168869693270611432679082, 0.0, 0.279705391489276667901467771423780,
      0.0, 0.381830050505118944950369775488975, 0.0, 0.417959183673469387755102040816327};
  double dinf = std::min(1, inf);
  double centr = 0.5 * (a + b), hlgth = 0.5 * (b - a);
  double tabsc1 = boun + dinf * (1.0 - centr) / centr;
  double fval1 = f(tabsc1);
  if (inf == 2) fval1 += f(-tabsc1);
  double fc = (fval1 / centr) / centr;
  double resg = wg[7] * fc, resk = wgk[7] * fc;
  double resabs = std::fabs(resk);
  double fv1[7], fv2[7];
  for (int j = 0; j < 7; ++j) {
    double absc = hlgth * xgk[j];
    double absc1 = centr - absc, absc2 = centr + absc;
    double t1 = boun + dinf * (1.0 - absc1) / absc1;
    double t2 = boun + dinf * (1.0 - absc2) / absc2;
    double f1 = f(t1), f2 = f(t2);
    if (inf == 2) { f1 += f(-t1); f2 += f(-t2); }
    f1 = (f1 / absc1) / absc1;
    f2 = (f2 / absc2) / absc2;
    fv1[j] = f1; fv2[j] = f2;
    double fsum = f1 + f2;
    resg += wg[j] * fsum;
    resk += wgk[j] * fsum;
    resabs += wgk[j] * (std::fabs(f1) + std::fabs(f2));
  }
  double reskh = resk * 0.5;
  double resasc = wgk[7] * std::fabs(fc - reskh);
  for (int j = 0; j < 7; ++j) resasc += wgk[j] * (std::fabs(fv1[j] - reskh) + std::fabs(fv2[j] - reskh));
  RuleResult r;
  r.result = resk * hlgth;
  resasc *= hlgth;
  resabs *= hlgth;
  r.resabs = resabs;
  r.resasc = resasc;
  r.abserr = scale_error(std::fabs((resk - resg) * hlgth), resabs, resasc);
  return r;
}

// QPSRT: maintain the descending ordering of error estimates (1-based arrays as in the book).
void qpsrt(int limit, int last, int& maxerr, double& ermax, const std::vector<double>& elist, std::vector<int>& iord,
           int& nrmax) {
  if (last <= 2) {
    iord[1] = 1;
    iord[2] = 2;
  } else {
    double errmax = elist[maxerr];
    if (nrmax != 1) {
      int ido = nrmax - 1;
      for (int i = 1; i <= ido; ++i) {
        int isucc = iord[nrmax - 1];
        if (errmax <= elist[isucc]) break;
        iord[nrmax] = isucc;
        --nrmax;
      }
    }
    int jupbn = last;
    if (last > limit / 2 + 2) jupbn = limit + 3 - last;
    double errmin = elist[last];
    int jbnd = jupbn - 1;
    int ibeg = nrmax + 1;
    int i = ibeg;
    bool inserted = false;
    for (; i <= jbnd; ++i) {
      int isucc = iord[i];
      if (errmax >= elist[isucc]) { inserted = true; break; }
      iord[i - 1] = isucc;
    }
    if (!inserted) {
      iord[jbnd] = maxerr;
      iord[jupbn] = last;
    } else {
      iord[i - 1] = maxerr;
      int k = jbnd;
      bool placed = false;
      for (int j = i; j <= jbnd; ++j) {
        int isucc = iord[k];
        if (errmin < elist[isucc]) { iord[k + 1] = last; placed = true; break; }
        iord[k + 1] = isucc;
        --k;
      }
      if (!placed) iord[i] = last;
    }
  }
  maxerr = iord[nrmax];
  ermax = elist[maxerr];
}

// QELG: Wynn's epsilon algorithm. epstab is 1-based with room for 52 + 2 entries; n is updated in place.
struct EpsilonTable {
  std::array<double, 56> epstab{};
  std::array<double, 4> res3la{};
  int nres = 0;
};

void qelg(int& n, EpsilonTable& t, double& result, double& abserr) {
  auto& e = t.epstab;
  ++t.nres;
  abserr = kOflow;
  result = e[n];
  if (n < 3) {
    abserr = std::max(abserr, 5.0 * kEpmach * std::fabs(result));
    return;
  }
  const int limexp = 50;
  e[n + 2] = e[n];
  int newelm = (n - 1) / 2;
  e[n] = kOflow;
  int num = n;
  int k1 = n;
  for (int i = 1; i <= newelm; ++i) {
    int k2 = k1 - 1, k3 = k1 - 2;
    double res = e[k1 + 2];
    double e0 = e[k3], e1 = e[k2], e2 = res;
    double e1abs = std::fabs(e1);
    double delta2 = e2 - e1, err2 = std::fabs(delta2);
    double tol2 = std::max(std::fabs(e2), e1abs) * kEpmach;
    double delta3 = e1 - e0, err3 = std::fabs(delta3);
    double tol3 = std::max(e1abs, std::fabs(e0)) * kEpmach;
    if (err2 <= tol2 && err3 <= tol3) {
      // e0, e1 and e2 are equal to within machine accuracy: convergence is assumed.
      result = res;
      abserr = err2 + err3;
      abserr = std::max(abserr, 5.0 * kEpmach * std::fabs(result));
      return;
    }
    double e3 = e[k1];
    e[k1] = e1;
    double delta1 = e1 - e3, err1 = std::fabs(delta1);
    double tol1 = std::max(e1abs, std::fabs(e3)) * kEpmach;
    double ss = 0.0;
    bool irregular = err1 <= tol1 || err2 <= tol2 || err3 <= tol3;
    if (!irregular) {
      ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
      irregular = std::fabs(ss * e1) <= 1e-4;
    }
    if (irregular) {  // omit part of the table
      n = i + i - 1;
      break;
    }
    res = e1 + 1.0 / ss;
    e[k1] = res;
    k1 -= 2;
    double error = err2 + std::fabs(res - e2) + err3;
    if (error <= abserr) {
      abserr = error;
      result = res;
    }
  }
  // shift the table
  if (n == limexp) n = 2 * (limexp / 2) - 1;
  int ib = (num % 2 == 0) ? 2 : 1;
  int ie = newelm + 1;
  for (int i = 1; i <= ie; ++i) {
    int ib2 = ib + 2;
    e[ib] = e[ib2];
    ib = ib2;
  }
  if (num != n) {
    int indx = num - n + 1;
    for (int i = 1; i <= n; ++i) {
      e[i] = e[indx];
      ++indx;
    }
  }
  if (t.nres < 4) {
    t.res3la[t.nres] = result;
    abserr = kOflow;
  } else {
    auto& r = t.res3la;
    abserr = std::fabs(result - r[3]) + std::fabs(result - r[2]) + std::fabs(result - r[1]);
    r[1] = r[2];
    r[2] = r[3];
    r[3] = result;
  }
  abserr = std::max(abserr, 5.0 * kEpmach * std::fabs(result));
}

// The adaptive driver shared by QAGS and QAGI (the book's two drivers differ only in the rule and the initial
// interval). rule(a, b) integrates over a subinterval.
template <class Rule>
IntegrateResult adapt(Rule&& rule, double a, double b, double epsabs, double epsrel, int limit) {
  IntegrateResult out;
  std::vector<double> alist(limit + 2), blist(limit + 2), rlist(limit + 2), elist(limit + 2);
  std::vector<int> iord(limit + 2);
  EpsilonTable tab;
  int ier = 0, last = 0;
  double result = 0.0, abserr = 0.0;
  alist[1] = a;
  blist[1] = b;

  int ierro = 0;
  RuleResult r0 = rule(a, b);
  result = r0.result;
  abserr = r0.abserr;
  double defabs = r0.resabs, resabs = r0.resasc;
  double dres = std::fabs(result);
  double errbnd = std::max(epsabs, epsrel * dres);
  last = 1;
  rlist[1] = result;
  elist[1] = abserr;
  iord[1] = 1;
  if (abserr <= 100.0 * kEpmach * defabs && abserr > errbnd) ier = 2;
  if (limit == 1) ier = 1;
  bool done = ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0;

  if (!done) {
    tab.epstab[1] = result;
    double errmax = abserr;
    int maxerr = 1;
    double area = result, errsum = abserr;
    abserr = kOflow;
    int nrmax = 1, numrl2 = 2, ktmin = 0;
    bool extrap = false, noext = false;
    int iroff1 = 0, iroff2 = 0, iroff3 = 0;
    int ksgn = (dres >= (1.0 - 50.0 * kEpmach) * defabs) ? 1 : -1;
    double small = 0.0, erlarg = 0.0, ertest = 0.0, correc = 0.0;
    bool sum_up = false;  // label 115 of the book: result = sum of the subinterval contributions

    for (last = 2; last <= limit; ++last) {
      // bisect the subinterval with the largest error estimate
      double a1 = alist[maxerr], b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
      double a2 = b1, b2 = blist[maxerr];
      double erlast = errmax;
      RuleResult r1 = rule(a1, b1), r2 = rule(a2, b2);
      double area1 = r1.result, error1 = r1.abserr, defab1 = r1.resasc;
      double area2 = r2.result, error2 = r2.abserr, defab2 = r2.resasc;
      double area12 = area1 + area2, erro12 = error1 + error2;
      errsum = errsum + erro12 - errmax;
      area = area + area12 - rlist[maxerr];
      if (defab1 != error1 && defab2 != error2) {
        if (std::fabs(rlist[maxerr] - area12) <= 1e-5 * std::fabs(area12) && erro12 >= 0.99 * errmax) {
          if (extrap) ++iroff2;
          else ++iroff1;
        }
        if (last > 10 && erro12 > errmax) ++iroff3;
      }
      rlist[maxerr] = area1;
      rlist[last] = area2;
      errbnd = std::max(epsabs, epsrel * std::fabs(area));
      if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;  // roundoff error
      if (iroff2 >= 5) ierro = 3;
      if (last == limit) ier = 1;  // maximum number of subdivisions
      // bad integrand behaviour at a point of the integration range
      if (std::max(std::fabs(a1), std::fabs(b2)) <= (1.0 + 100.0 * kEpmach) * (std::fabs(a2) + 1000.0 * kUflow))
        ier = 4;
      // append the newly created intervals to the list
      if (error2 > error1) {
        alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
      } else {
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
      }
      qpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
      if (errsum <= errbnd) { sum_up = true; break; }
      if (ier != 0) break;
      if (last == 2) {
        small = std::fabs(b - a) * 0.375;
        erlarg = errsum;
        ertest = errbnd;
        tab.epstab[2] = area;
        continue;
      }
      if (noext) continue;
      erlarg -= erlast;
      if (std::fabs(b1 - a1) > small) erlarg += erro12;
      if (!extrap) {
        // is the interval to be bisected next the smallest one?
        if (std::fabs(blist[maxerr] - alist[maxerr]) > small) continue;
        extrap = true;
        nrmax = 2;
      }
      bool do_extrapolate = true;
      if (ierro != 3 && erlarg > ertest) {
        // The smallest interval has the largest error: before bisecting, decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        int id = nrmax, jupbnd = last;
        if (last > 2 + limit / 2) jupbnd = limit + 3 - last;
        for (int k = id; k <= jupbnd; ++k) {
          maxerr = iord[nrmax];
          errmax = elist[maxerr];
          if (std::fabs(blist[maxerr] - alist[maxerr]) > small) { do_extrapolate = false; break; }
          ++nrmax;
        }
      }
      if (!do_extrapolate) continue;
      // perform extrapolation
      ++numrl2;
      tab.epstab[numrl2] = area;
      double reseps, abseps;
      qelg(numrl2, tab, reseps, abseps);
      ++ktmin;
      if (ktmin > 5 && abserr < 1e-3 * errsum) ier = 5;
      if (abseps < abserr) {
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::fabs(reseps));
        if (abserr <= ertest) break;
      }
      // prepare bisection of the smallest interval
      if (numrl2 == 1) noext = true;
      if (ier == 5) break;
      maxerr = iord[1];
      errmax = elist[maxerr];
      nrmax = 1;
      extrap = false;
      small *= 0.5;
      erlarg = errsum;
    }
    if (last > limit) last = limit;

    // set the final result and error estimate
    if (!sum_up) {
      bool divergence_test = false;
      if (abserr == kOflow) {
        sum_up = true;
      } else if (ier + ierro != 0) {
        if (ierro == 3) abserr += correc;
        if (ier == 0) ier = 3;
        if (result != 0.0 && area != 0.0) {
          if (abserr / std::fabs(result) > errsum / std::fabs(area)) sum_up = true;
          else divergence_test = true;
        } else if (abserr > errsum) {
          sum_up = true;
        } else if (area != 0.0) {
          divergence_test = true;
        }
      } else {
        divergence_test = true;
      }
      if (divergence_test) {
        if (!(ksgn == -1 && std::max(std::fabs(result), std::fabs(area)) <= defabs * 0.01)) {
          if (0.01 > result / area || result / area > 100.0 || errsum > std::fabs(area)) ier = 6;
        }
      }
    }
    if (sum_up) {
      result = 0.0;
      for (int k = 1; k <= last; ++k) result += rlist[k];
      abserr = errsum;
    }
    if (ier > 2) --ier;
  }
  out.value = result;
  out.abs_error = abserr;
  out.subdivisions = last;
  out.ier = ier;
  return out;
}

const char* ier_message(int ier) {
  switch (ier) {
    case 0: return "OK";
    case 1: return "maximum number of subdivisions reached";
    case 2: return "roundoff error was detected";
    case 3: return "extremely bad integrand behaviour";
    case 4: return "roundoff error is detected in the extrapolation table";
    case 5: return "the integral is probably divergent";
    default: return "the input is invalid";
  }
}

}  // namespace

IntegrateResult integrate(const Fn1& f, double lower, double upper, const IntegrateOptions& opt) {
  if (!f) throw Error("integrate: 'f' is not a function");
  if (std::isnan(lower) || std::isnan(upper)) throw Error("integrate: a limit is NA or NaN");
  if (opt.subdivisions < 1 || std::isnan(opt.rel_tol) || std::isnan(opt.abs_tol) ||
      (opt.abs_tol <= 0.0 && opt.rel_tol < std::max(50.0 * kEpmach, 0.5e-28)))
    throw Error("integrate: invalid parameter values");
  Integrand fn{f};
  IntegrateResult r;
  if (std::isfinite(lower) && std::isfinite(upper)) {
    r = adapt([&](double a, double b) { return qk21(fn, a, b); }, lower, upper, opt.abs_tol, opt.rel_tol,
              opt.subdivisions);
  } else {
    // R's mapping: a finite lower limit means (lower, +Inf), else a finite upper limit means (-Inf, upper),
    // else the whole line (so, as in R, integrate(f, 0, -Inf) integrates over (0, Inf)).
    int inf;
    double bound;
    if (std::isfinite(lower)) { inf = 1; bound = lower; }
    else if (std::isfinite(upper)) { inf = -1; bound = upper; }
    else { inf = 2; bound = 0.0; }
    r = adapt([&](double a, double b) { return qk15i(fn, bound, inf, a, b); }, 0.0, 1.0, opt.abs_tol, opt.rel_tol,
              opt.subdivisions);
  }
  r.evaluations = fn.evals;
  r.message = ier_message(r.ier);
  return r;
}

}  // namespace rcli
