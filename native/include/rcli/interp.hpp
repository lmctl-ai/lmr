// Interpolation with R spline()/splinefun() and approx()/approxfun() semantics.
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

// Cubic interpolating spline. method: "fmm" (Forsythe, Malcolm & Moler end conditions; R's default), "natural",
// "periodic" (the last y is replaced by the first, as R does), "monoH.FC" (Fritsch-Carlson monotone Hermite).
// x need not be sorted; duplicated x are collapsed by the mean of their y (R ties=mean). At least 2 distinct x.
// Outside the data range: the end cubic is extrapolated (fmm), linear (natural, monoH.FC), or wrapped (periodic).
class CubicSpline {
 public:
  CubicSpline(const Vec& x, const Vec& y, const std::string& method = "fmm");
  double operator()(double x) const;
  double deriv(double x, int order = 1) const;  // order 1..3
  const std::string& method() const { return method_; }

 private:
  std::string method_;
  Vec x_, y_, b_, c_, d_;
};

// Piecewise linear or constant interpolation. rule 1: NaN outside [min x, max x]; rule 2: the nearest end value.
// For "constant", f in [0, 1] mixes the left (0) and right (1) values, as R approx(f=).
class Approx {
 public:
  Approx(const Vec& x, const Vec& y, const std::string& method = "linear", int rule_left = 1, int rule_right = 1,
         double f = 0.0);
  double operator()(double x) const;

 private:
  Vec x_, y_;
  bool constant_;
  int rule_left_, rule_right_;
  double f_;
};

}  // namespace rcli
