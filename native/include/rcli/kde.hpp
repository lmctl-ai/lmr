// Kernel density estimation, bandwidth selectors, the empirical CDF and histogram helpers (clean-room).
//
// Conventions follow R's documented behaviour (?density, ?bw.nrd, ?ecdf, ?hist, ?nclass.Sturges, ?pretty); the
// algorithms come from the literature cited in kde.cpp. NaN or infinite data are errors (callers drop NA first).
//
// Agreement with R 4.3.3 (measured by the tests in tests/test_kde.cpp):
//   - bw_nrd0, bw_nrd: to the last bits.
//   - bw_ucv, bw_bcv, bw_sj ("ste" and "dpi") with the default nb = 1000: R bins the data and minimises / solves with
//     coarse tolerances (tol = 0.1 * lower); we reproduce the binning, search intervals and solver iterates (found by
//     black-box testing), so the results agree to 1e-13 relative. nb = 0 gives exact pairwise sums (not R's numbers:
//     differences of ~1e-3..1e-2 relative at R's tolerances).
//   - density: evaluated exactly at the grid points (R uses linear binning + FFT + interpolation), so values differ from
//     R by R's binning error (up to ~1e-3 * max(y) for the smooth kernels at n = 512); the grid (x) is R's exactly.
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli::kde {

// ---- bandwidth selectors (Gaussian-kernel scale: the bandwidth is the kernel standard deviation, as in R) ----------

// Silverman (1986, eq. 3.31): 0.9 * min(sd, IQR/1.34) * n^-1/5. If that spread is 0, it falls back to sd, then |x[0]|,
// then 1 (R's documented "positive result guaranteed"). Needs n >= 2.
double bw_nrd0(const Vec& x);
// Scott (1992): 1.06 * min(sd, IQR/1.34) * n^-1/5 (0 when the quartiles coincide, like R). Needs n >= 2.
double bw_nrd(const Vec& x);

struct BwOptions {
  int nb = 1000;          // bins for the pairwise distances (R's nb); 0 = exact pairwise distances (not R's numbers)
  double lower = kNaN;    // search interval; default 0.1 * hmax .. hmax with hmax = 1.144 * s * n^-1/5, where
                          // s = sd for ucv/bcv and min(sd, IQR/1.349) for SJ (as R)
  double upper = kNaN;
  double tol = kNaN;      // optimiser / root tolerance; default 0.1 * lower (R's default)
};
// Unbiased (least-squares) cross-validation (Rudemo 1982; Bowman 1984; Scott & Terrell 1987), Gaussian kernel,
// minimised over [lower, upper] by Brent's fmin like R's optimize().
double bw_ucv(const Vec& x, const BwOptions& opt = {});
// Biased cross-validation (Scott & Terrell 1987), Gaussian kernel.
double bw_bcv(const Vec& x, const BwOptions& opt = {});
// Sheather & Jones (1991): method "ste" (solve-the-equation) or "dpi" (direct plug-in). Throws "sample is too sparse to
// find TD" when the pilot estimate of int f'''^2 is not positive (e.g. IQR = 0).
double bw_sj(const Vec& x, const std::string& method = "ste", const BwOptions& opt = {});
// R's density(bw = <string>) names, case-insensitive: "nrd0", "nrd", "ucv", "bcv", "sj" (= "sj-ste"), "sj-ste",
// "sj-dpi".
double bandwidth(const Vec& x, const std::string& rule);

// ---- kernel density estimate --------------------------------------------------------------------------------------

// Kernels scaled so that bw is the kernel's standard deviation (R's convention): "gaussian", "epanechnikov",
// "rectangular", "triangular", "biweight", "cosine", "optcosine" (a unique prefix is accepted, as in R).
std::string kernel_name(const std::string& kernel);  // canonical name; throws for an unknown kernel
double kernel_value(const std::string& kernel, double u, double bw);  // K_bw(u)
double kernel_rk(const std::string& kernel);  // R(K) = int K^2 at bw = 1 (R's density(give.Rkern = TRUE))

struct DensityOptions {
  std::string bw_rule = "nrd0";  // used when bw is NaN
  double bw = kNaN;              // explicit bandwidth (kernel sd)
  double adjust = 1.0;           // the bandwidth used is adjust * bw
  std::string kernel = "gaussian";
  Vec weights;                   // empty: 1/n each; else non-negative, same length as x (not renormalised, like R)
  int n = 512;                   // grid points
  double from = kNaN, to = kNaN; // default min(x) - cut * bw, max(x) + cut * bw
  double cut = 3.0;
};
struct Density {
  Vec x, y;        // grid and density values (exact kernel sums)
  double bw = kNaN;
  size_t n = 0;    // sample size
  std::string kernel;
  Vec data, w;     // sorted observations and their weights, for at()
  double at(double x0) const;  // exact density at x0
};
Density density(const Vec& x, const DensityOptions& opt = {});
// Exact kernel estimate sum_i w_i K_bw(x0 - x_i) (w_i = 1/n when weights is empty).
double kde_at(const Vec& x, double x0, double bw, const std::string& kernel = "gaussian", const Vec& weights = {});

// ---- empirical CDF ------------------------------------------------------------------------------------------------

struct Ecdf {
  Vec sorted;   // the sample, sorted
  Vec knots;    // unique values (R's knots())
  Vec cumprob;  // Fn at each knot
  double operator()(double t) const;  // #{x_i <= t} / n
  Vec operator()(const Vec& t) const;
  // R's quantile(Fn, probs, type) equals quantile(x, probs, type); type 7 is R's default.
  double quantile(double p, int type = 7) const;
};
Ecdf ecdf(const Vec& x);

// ---- histograms ---------------------------------------------------------------------------------------------------

// Counts in the cells defined by breaks (strictly increasing): right = true gives (a, b] cells, the first closed when
// include_lowest; right = false gives [a, b), the last closed when include_lowest. Like R's hist(), the breaks are
// moved outward by fuzz * s so values within rounding of a break count as on it, with s = median(diff(breaks)) for
// >= 5 cells, min(diff(breaks)) for 3-4 cells and diff(range(x)) for 1-2 cells (black-box). A value outside all cells
// is an error, as in hist().
Vec hist_counts(const Vec& x, const Vec& breaks, bool right = true, bool include_lowest = true, double fuzz = 1e-7);

int nclass_sturges(const Vec& x);  // ceiling(log2(n) + 1)
int nclass_scott(const Vec& x);    // Scott (1979): h = 3.5 sd n^-1/3; 1 when h = 0
int nclass_fd(const Vec& x, int digits = 5);  // Freedman & Diaconis (1981) with R's documented fallbacks

// R's pretty(c(lo, hi), n, min.n, shrink.sml, high.u.bias, u5.bias): about n + 1 equally spaced round values covering
// [lo, hi]. The unit choice and rounding were characterised by black-box testing (see kde.cpp).
Vec pretty(double lo, double hi, int n = 5, int min_n = -1 /* n / 3 */, double shrink_sml = 0.75,
           double high_u_bias = 1.5, double u5_bias = kNaN /* 0.5 + 1.5 * high_u_bias */);

struct Histogram {
  Vec breaks, counts, density, mids;
  bool equidist = true;
};
// R's hist(x, breaks = rule, plot = FALSE) for rule "sturges", "scott", "fd" / "freedman-diaconis" (case-insensitive):
// breaks = pretty(range(x), nclass, min.n = 1).
Vec hist_breaks(const Vec& x, const std::string& rule = "sturges");
Histogram hist(const Vec& x, const std::string& rule = "sturges");
Histogram hist(const Vec& x, const Vec& breaks, bool right = true, bool include_lowest = true);

}  // namespace rcli::kde
