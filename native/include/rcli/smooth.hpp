// Scatterplot smoothers: a penalised cubic smoothing spline and local regression (loess).
// These are clean-room implementations from the literature (Reinsch 1967 / Green & Silverman 1994 for the spline,
// Cleveland 1979 / Cleveland, Grosse & Shyu 1992 for loess). They match R's behaviour within documented tolerances,
// not bit for bit: R's smooth.spline places knots and scales spar in implementation-specific ways.
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

struct SmoothSplineOptions {
  double df = kNaN;    // target equivalent degrees of freedom (trace of the smoother matrix), 1 < df <= nx
  double lambda = kNaN;  // explicit penalty (takes precedence over df)
  bool cv = false;     // when df and lambda are both absent: false = GCV (R default), true = ordinary leave-one-out CV
};

struct SmoothSplineFit {
  Vec x;               // unique sorted x
  Vec y;               // fitted values at x
  Vec w;               // weights per unique x (counts of ties)
  Vec lev;             // leverages at x
  double df = kNaN;    // trace of the smoother matrix
  double lambda = kNaN;
  double cv_crit = kNaN;  // GCV or CV score at the chosen lambda
  double pen_crit = kNaN; // weighted residual sum of squares
  double predict(double x0, int deriv = 0) const;  // natural cubic spline through the fit; linear beyond the ends
  Vec coef_b, coef_c, coef_d;  // piecewise cubic coefficients on intervals of x (for predict)
};

SmoothSplineFit smooth_spline(const Vec& x, const Vec& y, const SmoothSplineOptions& opt = {});

struct LoessOptions {
  double span = 0.75;
  int degree = 2;                    // 0, 1 or 2
  std::string family = "gaussian";  // or "symmetric" (4 bisquare robustness iterations, as R)
};

struct LoessFit {
  Vec x, y;            // data as given
  Vec fitted, residuals;
  double enp = kNaN;   // equivalent number of parameters, tr(L'L) as R's $enp
  double trace_hat = kNaN;  // tr(L)
  double one_delta = kNaN, two_delta = kNaN;  // trace(I-L)'(I-L) and trace((I-L)'(I-L))^2
  double s = kNaN;     // residual standard error
  Vec robustness_weights;
  LoessOptions opt;
  double predict(double x0) const;  // direct local fit at x0 (NaN outside the data range, like R)
};

// Uses exact local fitting at every point (R's surface = "direct", statistics = "exact"); R's default "interpolate"
// surface is a kd-tree vertex interpolation that can differ by several percent of the noise level on small noisy data.
// Exact statistics cost O(n^3): about 1 s at n = 1000.
LoessFit loess(const Vec& x, const Vec& y, const LoessOptions& opt = {});

}  // namespace rcli
