// Tests for rcli::summary against R 4.3.3 (values from Rscript, printed with %.17g).
// Data: set.seed(42); x <- round(rnorm(23)*10, 3); y <- round(x*0.5 + rnorm(23)*5, 2); y[5] <- y[6]; x[7] <- x[8]
#include "rcli/summary.hpp"

#include "check.hpp"

using rcli::Vec;
namespace s = rcli::summary;

int main() {
  const Vec x = {13.710000000000001, -5.6470000000000002, 3.6309999999999998, 6.3289999999999997, 4.0430000000000001,
                 -1.0609999999999999, -0.94699999999999995, -0.94699999999999995, 20.184000000000001, -0.627,
                 13.048999999999999, 22.866, -13.888999999999999, -2.7879999999999998, -1.333, 6.3600000000000003,
                 -2.843, -26.565000000000001, -24.405000000000001, 13.201000000000001, -3.0659999999999998,
                 -17.812999999999999, -1.7190000000000001};
  const Vec y = {12.93, 6.6500000000000004, -0.34000000000000002, 1.8799999999999999, 1.77, 1.77, 4.3600000000000003,
                 1.8, 13.619999999999999, 4.8600000000000003, 3.48, 13.960000000000001, -15.529999999999999,
                 -5.3200000000000003, -4.9199999999999999, -8.8900000000000006, -1.24, -12.25, -14.01,
                 10.390000000000001, -5.1699999999999999, -15.75, 1.3};
  const double tol = 1e-15;

  // quantile(x, c(0, .1, .25, .5, .33, .9, 1), type = t)
  const Vec p = {0, 0.1, 0.25, 0.5, 0.33, 0.9, 1};
  const double qx[9][7] = {
      {-26.565000000000001, -17.812999999999999, -3.0659999999999998, -0.94699999999999995, -2.7879999999999998, 13.710000000000001, 22.866},
      {-26.565000000000001, -17.812999999999999, -3.0659999999999998, -0.94699999999999995, -2.7879999999999998, 13.710000000000001, 22.866},
      {-26.565000000000001, -24.405000000000001, -3.0659999999999998, -0.94699999999999995, -2.7879999999999998, 13.710000000000001, 22.866},
      {-26.565000000000001, -22.427399999999999, -3.7112500000000002, -1.004, -2.8105500000000001, 13.5573, 22.866},
      {-26.565000000000001, -19.131399999999999, -3.0102500000000001, -0.94699999999999995, -2.6917900000000001, 15.004799999999996, 22.866},
      {-26.565000000000001, -21.768199999999997, -3.0659999999999998, -0.94699999999999995, -2.7923999999999998, 17.594400000000007, 22.866},
      {-26.565000000000001, -17.028199999999998, -2.9544999999999999, -0.94699999999999995, -2.5100599999999984, 13.608200000000002, 22.866},
      {-26.565000000000001, -20.010333333333328, -3.028833333333333, -0.94699999999999995, -2.7523666666666649, 15.868000000000016, 22.866},
      {-26.565000000000001, -19.790599999999998, -3.0241875, -0.94699999999999995, -2.7372225000000006, 15.652200000000004, 22.866}};
  for (int t = 1; t <= 9; ++t) {
    const Vec q = s::quantiles(x, p, t);
    for (size_t i = 0; i < p.size(); ++i)
      check::close(q[i], qx[t - 1][i], tol, "quantile(x, " + std::to_string(p[i]) + ", type=" + std::to_string(t) + ")");
  }
  // quantile(1:10, c(0.1, 0.5, 0.25, 0.75, 0.05), type = t): exercises the 4*eps fuzz of the discontinuous types
  Vec one10;
  for (int i = 1; i <= 10; ++i) one10.push_back(i);
  const Vec p2 = {0.1, 0.5, 0.25, 0.75, 0.05};
  const double q10[9][5] = {{1, 5, 3, 8, 1},
                            {1.5, 5.5, 3, 8, 1},
                            {1, 5, 2, 8, 1},
                            {1, 5, 2.5, 7.5, 1},
                            {1.5, 5.5, 3, 8, 1},
                            {1.1000000000000001, 5.5, 2.75, 8.25, 1},
                            {1.8999999999999999, 5.5, 3.25, 7.75, 1.45},
                            {1.3666666666666665, 5.4999999999999991, 2.9166666666666665, 8.0833333333333321, 1},
                            {1.4000000000000001, 5.5, 2.9375, 8.0625, 1}};
  for (int t = 1; t <= 9; ++t) {
    const Vec q = s::quantiles(one10, p2, t);
    for (size_t i = 0; i < p2.size(); ++i)
      check::close(q[i], q10[t - 1][i], tol, "quantile(1:10, " + std::to_string(p2[i]) + ", type=" + std::to_string(t) + ")");
  }
  check::close(s::quantile(x, 0.5), -0.94699999999999995, tol, "quantile(x, 0.5)  # default type 7");
  check::ok(std::isnan(s::quantile(x, rcli::kNaN)), "quantile(x, NA) is NA");
  check::close(s::quantile(x, 1 + 1e-15, 7), 22.866, tol, "quantile(x, 1+1e-15) accepted (fuzz 2e-14)");
  check::throws([&] { s::quantile(x, 1.1); }, "quantile(x, 1.1) errors");
  check::throws([&] { s::quantile(x, 0.5, 10); }, "quantile(type=10) errors");

  check::close(s::mean(x), -0.012043478260869513, tol, "mean(x)");
  check::close(s::var(x), 156.79484749802373, tol, "var(x)");
  check::close(s::sd(x), 12.521774934010901, tol, "sd(x)");
  check::close(s::median(x), -0.94699999999999995, tol, "median(x)");
  check::close(s::median(y), 1.77, tol, "median(y)");
  check::close(s::median({0.1, 0.7, 0.3, 0.2}), 0.25, tol, "median(c(0.1,0.7,0.3,0.2))");
  check::close(s::mad(x), 7.398174, tol, "mad(x)");
  check::close(s::cov(x, y), 94.044802628458498, tol, "cov(x, y)");
  check::close(s::cor(x, y), 0.8387982647505724, tol, "cor(x, y)");
  check::close(s::cor(x, y, "spearman"), 0.75729115175481954, tol, "cor(x, y, method='spearman')");
  check::close(s::cor(x, y, "kendall"), 0.62301587301587302, tol, "cor(x, y, method='kendall')");
  check::close(s::cor({1, 1, 2, 2, 3}, {2, 1, 2, 3, 3}, "kendall"), 0.75, tol, "cor(c(1,1,2,2,3), c(2,1,2,3,3), method='kendall')");
  check::close(s::cor({1, 1, 2, 2, 3}, {2, 1, 2, 3, 3}, "spearman"), 0.80555555555555558, tol,
               "cor(c(1,1,2,2,3), c(2,1,2,3,3), method='spearman')");
  check::close(s::mean({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}), 0.40000000000000002, tol, "mean(c(.1,.2,...,.7))");
  check::close(s::mean({1e20, 1, -1e20, 3}), 1, tol, "mean(c(1e20,1,-1e20,3))");
  check::close(s::sum({1e20, 1, -1e20, 3}), 4, tol, "sum(c(1e20,1,-1e20,3)) exact by compensation (R gives 0)");

  // rank(y)
  const Vec ry = {21, 19, 10, 15, 12.5, 12.5, 17, 14, 22, 18, 16, 23, 2, 6, 8, 5, 9, 4, 3, 20, 7, 1, 11};
  const Vec r = s::rank(y);
  for (size_t i = 0; i < ry.size(); ++i) check::close(r[i], ry[i], 0, "rank(y)[" + std::to_string(i + 1) + "]");

  // Degenerate inputs: R returns NA / NaN
  check::ok(std::isnan(s::var({1.0})), "var(1) is NA");
  check::ok(std::isnan(s::mean({})), "mean(numeric(0)) is NaN");
  check::ok(std::isnan(s::cor({1, 1, 1}, {1, 2, 3})), "cor with zero sd is NA");
  check::throws([] { s::mean({1, rcli::kNaN}); }, "mean with NaN errors");
  check::throws([] { s::cor({1, 2}, {1, 2, 3}); }, "cor length mismatch errors");
  check::throws([] { s::cor({1, 2}, {1, 2}, "foo"); }, "cor unknown method errors");
  return check::report("test_summary");
}
