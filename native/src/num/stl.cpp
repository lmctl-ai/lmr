// STL: Seasonal-Trend decomposition based on Loess.
// Implemented from R. B. Cleveland, W. S. Cleveland, J. E. McRae & I. Terpenning (1990), "STL: A Seasonal-Trend
// Decomposition Procedure Based on Loess", Journal of Official Statistics 6(1), 3-73:
//   inner loop (sec. 2.2): detrending, cycle-subseries smoothing (extended one period at each end), low-pass filter
//   of the smoothed cycle-subseries (moving averages of length period, period, 3 followed by loess), detrending of the
//   smoothed cycle-subseries, deseasonalising, trend smoothing;
//   outer loop (sec. 2.3): bisquare robustness weights on |R| / (6 median |R|);
//   loess (sec. 2.1): q nearest neighbours with tricube weights, local constant or local linear fit; when q exceeds
//   the number of points the neighbourhood distance is enlarged (sec. 2.1, last paragraph);
//   the "jump" speed-ups (sec. 2.4 / ?stl): loess evaluated every jump-th point and linearly interpolated.
// The R-level conventions (defaults, the periodic case and the returned weights) follow ?stl. Tiny-distance /
// tiny-residual weights are exactly 1 and near-boundary ones exactly 0 (cut-offs at 0.001 and 0.999 of the scale),
// as R's returned robustness weights show.
// Known deviation: for robust fits of even-length series R occasionally (at some outer iterations) uses a robustness
// scale that is not 6 * median|R| (it averages the lower middle order statistic with a different one); this
// implementation always uses the median, so such fits can differ from R by ~1e-4..1e-2. Odd-length series and
// non-robust fits agree with R to ~1e-13.
#include <algorithm>
#include <cmath>
#include "rcli/ts.hpp"

namespace rcli {

namespace {

// Loess estimate at abscissa xs (1-based, may lie outside 1..n) from y[nleft..nright] (1-based, inclusive).
// q is the nominal neighbourhood size. Returns false when all weights vanish.
bool loess_at(const Vec& y, int n, int q, int degree, double xs, int nleft, int nright, Vec& w, bool use_rw,
              const Vec& rw, double& ys) {
  double h = std::max(xs - nleft, nright - xs);
  if (q > n) h += static_cast<double>((q - n) / 2);
  const double h9 = 0.999 * h, h1 = 0.001 * h;
  double a = 0.0;
  for (int j = nleft; j <= nright; ++j) {
    double r = std::fabs(j - xs);
    double wj = 0.0;
    if (r <= h9) {
      if (r <= h1) wj = 1.0;
      else {
        double u = r / h;
        u = 1.0 - u * u * u;
        wj = u * u * u;
      }
      if (use_rw) wj *= rw[static_cast<size_t>(j - 1)];
    }
    w[static_cast<size_t>(j - 1)] = wj;
    a += wj;
  }
  if (a <= 0.0) return false;
  for (int j = nleft; j <= nright; ++j) w[static_cast<size_t>(j - 1)] /= a;
  if (h > 0.0 && degree > 0) {
    // Local linear fit written as a weighted combination: with normalised weights w, centre c = sum w_j j and
    // spread b = sum w_j (j - c)^2, the fitted value at xs is sum w_j (1 + (xs - c)(j - c) / b) y_j.
    double c = 0.0;
    for (int j = nleft; j <= nright; ++j) c += w[static_cast<size_t>(j - 1)] * j;
    double b = xs - c, s = 0.0;
    for (int j = nleft; j <= nright; ++j) s += w[static_cast<size_t>(j - 1)] * (j - c) * (j - c);
    double range = static_cast<double>(n) - 1.0;
    if (std::sqrt(s) > 0.001 * range) {
      b /= s;
      for (int j = nleft; j <= nright; ++j) w[static_cast<size_t>(j - 1)] *= b * (j - c) + 1.0;
    }
  }
  double v = 0.0;
  for (int j = nleft; j <= nright; ++j) v += w[static_cast<size_t>(j - 1)] * y[static_cast<size_t>(j - 1)];
  ys = v;
  return true;
}

// Loess smooth of y (length n) evaluated at 1..n, every `jump`-th point with linear interpolation in between.
Vec loess_smooth(const Vec& y, int q, int degree, int jump, bool use_rw, const Vec& rw) {
  const int n = static_cast<int>(y.size());
  Vec ys(y.size());
  if (n < 2) { ys = y; return ys; }
  Vec w(y.size());
  int nj = std::min(jump, n - 1);
  int nleft = 1, nright = n;
  auto est = [&](int i) {
    double v;
    if (loess_at(y, n, q, degree, i, nleft, nright, w, use_rw, rw, v)) ys[static_cast<size_t>(i - 1)] = v;
    else ys[static_cast<size_t>(i - 1)] = y[static_cast<size_t>(i - 1)];
  };
  if (q >= n) {
    nleft = 1; nright = n;
    for (int i = 1; i <= n; i += nj) est(i);
  } else if (nj == 1) {
    int nsh = (q + 1) / 2;
    nleft = 1; nright = q;
    for (int i = 1; i <= n; ++i) {
      if (i > nsh && nright != n) { ++nleft; ++nright; }
      est(i);
    }
  } else {
    int nsh = (q + 1) / 2;
    for (int i = 1; i <= n; i += nj) {
      if (i < nsh) { nleft = 1; nright = q; }
      else if (i >= n - nsh + 1) { nleft = n - q + 1; nright = n; }
      else { nleft = i - nsh + 1; nright = q + i - nsh; }
      est(i);
    }
  }
  if (nj != 1) {
    int last = 1;
    for (int i = 1; i <= n - nj; i += nj) {
      double delta = (ys[static_cast<size_t>(i + nj - 1)] - ys[static_cast<size_t>(i - 1)]) / nj;
      for (int j = i + 1; j < i + nj; ++j) ys[static_cast<size_t>(j - 1)] = ys[static_cast<size_t>(i - 1)] + delta * (j - i);
      last = i + nj;
    }
    if (last != n) {
      // evaluate the last point directly, then interpolate the gap from the last evaluated point
      if (q >= n) { nleft = 1; nright = n; }
      else { nleft = n - q + 1; nright = n; }
      est(n);
      if (last != n - 1) {
        double delta = (ys[static_cast<size_t>(n - 1)] - ys[static_cast<size_t>(last - 1)]) / (n - last);
        for (int j = last + 1; j < n; ++j) ys[static_cast<size_t>(j - 1)] = ys[static_cast<size_t>(last - 1)] + delta * (j - last);
      }
    }
  }
  return ys;
}

// Moving average of length len: output length y.size() - len + 1.
Vec moving_average(const Vec& y, int len) {
  const size_t n = y.size(), L = static_cast<size_t>(len);
  Vec out(n - L + 1);
  double s = 0.0;
  for (size_t i = 0; i < L; ++i) s += y[i];
  out[0] = s / len;
  for (size_t j = 1; j + L <= n; ++j) {
    s = s - y[j - 1] + y[j + L - 1];
    out[j] = s / len;
  }
  return out;
}

struct StlParams {
  int np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump;
};

// Cycle-subseries smoothing; output has length n + 2 * np (each subseries extended one value at each end).
Vec subseries_smooth(const Vec& y, const StlParams& p, bool use_rw, const Vec& rw) {
  const int n = static_cast<int>(y.size()), np = p.np;
  Vec c(static_cast<size_t>(n + 2 * np), 0.0);
  for (int j = 0; j < np; ++j) {
    int k = (n - j - 1) / np + 1;  // number of values in this subseries
    Vec sub(static_cast<size_t>(k)), srw(static_cast<size_t>(k), 1.0);
    for (int i = 0; i < k; ++i) {
      sub[static_cast<size_t>(i)] = y[static_cast<size_t>(i * np + j)];
      if (use_rw) srw[static_cast<size_t>(i)] = rw[static_cast<size_t>(i * np + j)];
    }
    Vec sm = loess_smooth(sub, p.ns, p.isdeg, p.nsjump, use_rw, srw);
    Vec w(static_cast<size_t>(k));
    double v;
    // extension at position 0
    int nright = std::min(p.ns, k);
    double first = sub[0];
    if (loess_at(sub, k, p.ns, p.isdeg, 0.0, 1, nright, w, use_rw, srw, v)) first = v;
    // extension at position k + 1
    int nleft = std::max(1, k - p.ns + 1);
    double lastv = sub[static_cast<size_t>(k - 1)];
    if (loess_at(sub, k, p.ns, p.isdeg, k + 1.0, nleft, k, w, use_rw, srw, v)) lastv = v;
    c[static_cast<size_t>(j)] = first;
    for (int i = 0; i < k; ++i) c[static_cast<size_t>((i + 1) * np + j)] = sm[static_cast<size_t>(i)];
    c[static_cast<size_t>((k + 1) * np + j)] = lastv;
  }
  return c;
}

void inner_loop(const Vec& y, const StlParams& p, int inner, bool use_rw, const Vec& rw, Vec& season, Vec& trend) {
  const size_t n = y.size();
  const Vec ones(n, 1.0);
  for (int it = 0; it < inner; ++it) {
    Vec detr(n);
    for (size_t i = 0; i < n; ++i) detr[i] = y[i] - trend[i];
    Vec c = subseries_smooth(detr, p, use_rw, rw);
    Vec l = moving_average(moving_average(moving_average(c, p.np), p.np), 3);
    l = loess_smooth(l, p.nl, p.ildeg, p.nljump, false, ones);
    for (size_t i = 0; i < n; ++i) season[i] = c[static_cast<size_t>(p.np) + i] - l[i];
    Vec des(n);
    for (size_t i = 0; i < n; ++i) des[i] = y[i] - season[i];
    trend = loess_smooth(des, p.nt, p.itdeg, p.ntjump, use_rw, rw);
  }
}

double median_of(Vec v) {
  const size_t n = v.size();
  size_t m = n / 2;
  std::nth_element(v.begin(), v.begin() + static_cast<long>(m), v.end());
  double hi = v[m];
  if (n % 2 == 1) return hi;
  double lo = *std::max_element(v.begin(), v.begin() + static_cast<long>(m));
  return 0.5 * (lo + hi);
}

// Bisquare robustness weights B(|R| / h), h = 6 median |R| (Cleveland et al. 1990, sec. 2.3).
void robustness_weights(const Vec& y, const Vec& fit, Vec& rw) {
  const size_t n = y.size();
  Vec r(n);
  for (size_t i = 0; i < n; ++i) r[i] = std::fabs(y[i] - fit[i]);
  double cmad = 6.0 * median_of(r);
  double c9 = 0.999 * cmad, c1 = 0.001 * cmad;
  for (size_t i = 0; i < n; ++i) {
    if (r[i] <= c1) rw[i] = 1.0;
    else if (r[i] <= c9) {
      double u = r[i] / cmad;
      u = 1.0 - u * u;
      rw[i] = u * u;
    } else rw[i] = 0.0;
  }
}

int next_odd(double x) {
  long v = std::lround(x);
  if (v % 2 == 0) ++v;
  return static_cast<int>(v);
}

}  // namespace

StlFit stl(const Vec& x, int period, const StlOptions& opt) {
  const int n = static_cast<int>(x.size());
  for (double v : x)
    if (!std::isfinite(v)) throw Error("stl: 'x' contains missing or non-finite values");
  if (period < 2 || n <= 2 * period) throw Error("stl: series is not periodic or has less than two periods");
  bool periodic = opt.periodic;
  if (!periodic && opt.s_window <= 0) throw Error("stl: 's.window' must be given (or periodic = true)");
  if (opt.s_degree < 0 || opt.s_degree > 1) throw Error("stl: 's.degree' must be 0 or 1");
  if (opt.t_degree < 0 || opt.t_degree > 1) throw Error("stl: 't.degree' must be 0 or 1");
  StlParams p{};
  p.np = period;
  int s_degree = opt.s_degree;
  double s_window;
  if (periodic) {
    s_window = 10.0 * n + 1.0;
    s_degree = 0;
  } else {
    s_window = opt.s_window;
  }
  int t_window = opt.t_window > 0 ? opt.t_window : next_odd(std::ceil(1.5 * period / (1.0 - 1.5 / s_window)));
  int l_window = opt.l_window > 0 ? opt.l_window : next_odd(period);
  p.ns = static_cast<int>(s_window);
  p.nt = t_window;
  p.nl = l_window;
  p.isdeg = s_degree;
  p.itdeg = opt.t_degree;
  p.ildeg = opt.t_degree;  // l.degree defaults to t.degree
  p.nsjump = static_cast<int>(std::ceil(s_window / 10.0));
  p.ntjump = static_cast<int>(std::ceil(t_window / 10.0));
  p.nljump = static_cast<int>(std::ceil(l_window / 10.0));
  // windows must be odd and at least 3; jumps at least 1
  auto fix_window = [](int& w) {
    if (w < 3) w = 3;
    if (w % 2 == 0) ++w;
  };
  fix_window(p.ns);
  fix_window(p.nt);
  fix_window(p.nl);
  p.nsjump = std::max(1, p.nsjump);
  p.ntjump = std::max(1, p.ntjump);
  p.nljump = std::max(1, p.nljump);
  const int inner = opt.robust ? 1 : 2;
  const int outer = opt.robust ? 15 : 0;

  const size_t N = x.size();
  Vec season(N, 0.0), trend(N, 0.0), rw(N, 1.0);
  bool use_rw = false;
  for (int k = 0;; ++k) {
    inner_loop(x, p, inner, use_rw, rw, season, trend);
    if (k + 1 > outer) break;
    Vec fit(N);
    for (size_t i = 0; i < N; ++i) fit[i] = season[i] + trend[i];
    robustness_weights(x, fit, rw);
    use_rw = true;
  }
  if (outer <= 0) std::fill(rw.begin(), rw.end(), 1.0);
  if (periodic) {
    // ?stl: in the periodic case the seasonal component is the mean of each cycle position
    for (int j = 0; j < period; ++j) {
      double s = 0.0;
      int cnt = 0;
      for (int i = j; i < n; i += period) { s += season[static_cast<size_t>(i)]; ++cnt; }
      s /= cnt;
      for (int i = j; i < n; i += period) season[static_cast<size_t>(i)] = s;
    }
  }
  StlFit out;
  out.seasonal = season;
  out.trend = trend;
  out.remainder.resize(N);
  for (size_t i = 0; i < N; ++i) out.remainder[i] = x[i] - season[i] - trend[i];
  out.weights = rw;
  return out;
}

}  // namespace rcli
