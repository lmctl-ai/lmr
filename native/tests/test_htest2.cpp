// Tests for the k-sample, block, symmetry and scale tests (htest2.cpp) against R 4.3.3 used as a black-box oracle.
// The block between BEGIN/END GENERATED was produced by evaluating each quoted R expression with Rscript and printing
// every component with 17 significant digits. Data: datasets::airquality, InsectSprays, PlantGrowth, warpbreaks;
// RoundingTimes (?friedman.test), dataFreq = the store x brand matrix y of ?quade.test (Conover 1999, p. 375f),
// the Performance table of ?mcnemar.test, ramsay / jung.parekh of ?ansari.test (Hollander & Wolfe 1973, p. 86f),
// mucx / mucg the Hollander & Wolfe (1973, p. 116) mucociliary data of ?kruskal.test.
// u1 <- {set.seed(7); round(rnorm(40), 3)}; u2 <- round(rnorm(45, 0, 2), 3) (continuing the stream);
// big1 <- {set.seed(11); rnorm(60)}; big2 <- rnorm(55, 0, 1.5).
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "check.hpp"
#include "rcli/htest.hpp"

using namespace rcli;

namespace {

struct O {
  TestOptions o;
  O& alt(const std::string& a) { o.alternative = a; return *this; }
  O& exact(int v) { o.exact = v; return *this; }
  O& ci() { o.conf_int = true; return *this; }
  operator TestOptions() const { return o; }  // NOLINT: implicit by design
};

Matrix M(size_t r, size_t c, const Vec& v) {
  Matrix m(r, c);
  m.a = v;
  return m;
}
Vec head(const Vec& a, size_t n) { return Vec(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(n)); }

using Named = std::vector<std::pair<std::string, double>>;
struct Exp {
  std::string method, alternative, stat_name;
  double stat;
  Named params;
  double p;
  std::string null_name;
  double null_value;
};

std::map<std::string, double> worst;
double reldev(double a, double e) {
  if (std::isnan(e)) return std::isnan(a) ? 0.0 : kInf;
  if (std::isinf(e)) return a == e ? 0.0 : kInf;
  if (std::fabs(e) < 1e-300) return std::fabs(a - e);
  return std::fabs(a - e) / std::fabs(e);
}
void track(const std::string& key, double a, double e) {
  double& w = worst[key];
  w = std::max(w, reldev(a, e));
}

constexpr double kTolStat = 1e-12, kTolP = 1e-12;

void expect(const std::string& fam, const std::string& what, const HtestResult& r, const Exp& e) {
  check::ok(r.method == e.method, what + ": method '" + r.method + "' vs '" + e.method + "'");
  if (!e.alternative.empty())
    check::ok(r.alternative == e.alternative, what + ": alternative '" + r.alternative + "'");
  check::ok(r.statistic_name == e.stat_name, what + ": statistic name '" + r.statistic_name + "'");
  check::close(r.statistic, e.stat, kTolStat, what + ": statistic", 1e-300);
  track(fam + " statistic", r.statistic, e.stat);
  check::ok(r.parameters.size() == e.params.size(), what + ": parameter count");
  for (size_t i = 0; i < std::min(r.parameters.size(), e.params.size()); ++i) {
    check::ok(r.parameters[i].first == e.params[i].first, what + ": parameter name " + r.parameters[i].first);
    check::close(r.parameters[i].second, e.params[i].second, kTolStat, what + ": parameter " + e.params[i].first);
    track(fam + " parameter", r.parameters[i].second, e.params[i].second);
  }
  check::close(r.p_value, e.p, kTolP, what + ": p-value", 1e-300);
  track(fam + " p-value", r.p_value, e.p);
  check::ok(!r.has_conf_int && r.estimates.empty(), what + ": no interval or estimates");
  check::ok(r.null_name == e.null_name, what + ": null name '" + r.null_name + "'");
  if (!e.null_name.empty()) check::close(r.null_value, e.null_value, 0, what + ": null value");
}

// BEGIN GENERATED
const Vec oz = {41, 36, 12, 18, kNaN, 28, 23, 19, 8, kNaN, 7, 16, 11, 14, 18, 14, 34, 6, 30, 11, 1, 11, 4, 32, kNaN, kNaN, kNaN, 23, 45, 115, 37, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, 29, kNaN, 71, 39, kNaN, kNaN, 23, kNaN, kNaN, 21, 37, 20, 12, 13, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, 135, 49, 32, kNaN, 64, 40, 77, 97, 97, 85, kNaN, 10, 27, kNaN, 7, 48, 35, 61, 79, 63, 16, kNaN, kNaN, 80, 108, 20, 52, 82, 50, 64, 59, 39, 9, 16, 78, 35, 66, 122, 89, 110, kNaN, kNaN, 44, 28, 65, kNaN, 22, 59, 23, 31, 44, 21, 9, kNaN, 45, 168, 73, kNaN, 76, 118, 84, 85, 96, 78, 73, 91, 47, 32, 20, 23, 21, 24, 44, 21, 28, 9, 13, 46, 18, 13, 24, 16, 13, 23, 36, 7, 14, 30, kNaN, 14, 18, 20};
const Vec ozm = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9};
const Vec ins = {10, 7, 20, 14, 14, 12, 10, 23, 17, 20, 14, 13, 11, 17, 21, 11, 16, 14, 17, 17, 19, 21, 7, 13, 0, 1, 7, 2, 3, 1, 2, 1, 3, 0, 1, 4, 3, 5, 12, 6, 4, 3, 5, 5, 5, 5, 2, 4, 3, 5, 3, 5, 3, 6, 1, 1, 3, 2, 6, 4, 11, 9, 15, 22, 15, 16, 13, 10, 26, 26, 24, 13};
const Vec insg = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};
const Vec mucx = {2.8999999999999999, 3, 2.5, 2.6000000000000001, 3.2000000000000002, 3.7999999999999998, 2.7000000000000002, 4, 2.3999999999999999, 2.7999999999999998, 3.3999999999999999, 3.7000000000000002, 2.2000000000000002, 2};
const Vec mucg = {1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3};
const Vec rt = {5.4000000000000004, 5.5, 5.5499999999999998, 5.8499999999999996, 5.7000000000000002, 5.75, 5.2000000000000002, 5.5999999999999996, 5.5, 5.5499999999999998, 5.5, 5.4000000000000004, 5.9000000000000004, 5.8499999999999996, 5.7000000000000002, 5.4500000000000002, 5.5499999999999998, 5.5999999999999996, 5.4000000000000004, 5.4000000000000004, 5.3499999999999996, 5.4500000000000002, 5.5, 5.3499999999999996, 5.25, 5.1500000000000004, 5, 5.8499999999999996, 5.7999999999999998, 5.7000000000000002, 5.25, 5.2000000000000002, 5.0999999999999996, 5.6500000000000004, 5.5499999999999998, 5.4500000000000002, 5.5999999999999996, 5.3499999999999996, 5.4500000000000002, 5.0499999999999998, 5, 4.9500000000000002, 5.5, 5.5, 5.4000000000000004, 5.4500000000000002, 5.5499999999999998, 5.5, 5.5499999999999998, 5.5499999999999998, 5.3499999999999996, 5.4500000000000002, 5.5, 5.5499999999999998, 5.5, 5.4500000000000002, 5.25, 5.6500000000000004, 5.5999999999999996, 5.4000000000000004, 5.7000000000000002, 5.6500000000000004, 5.5499999999999998, 6.2999999999999998, 6.2999999999999998, 6.25};
const Vec rtg = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};
const Vec rtb = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22};
const Vec qy = {5, 4, 7, 10, 12, 1, 3, 1, 0, 2, 16, 12, 22, 22, 35, 5, 4, 3, 5, 4, 10, 9, 7, 13, 10, 19, 18, 28, 37, 58, 10, 7, 6, 8, 7};
const Vec qyg = {1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5};
const Vec qyb = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7};
const Vec wbx = {44.555555555555557, 28.222222222222221, 24, 28.777777777777779, 24.555555555555557, 18.777777777777779};
const Vec wbw = {1, 2, 1, 2, 1, 2};
const Vec wbt = {1, 1, 2, 2, 3, 3};
const Vec ramsay = {111, 107, 100, 99, 102, 106, 109, 108, 104, 99, 101, 96, 97, 102, 107, 113, 116, 113, 110, 98};
const Vec jung = {107, 108, 106, 98, 105, 103, 110, 105, 104, 100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99};
const Vec pg = {4.1699999999999999, 5.5800000000000001, 5.1799999999999997, 6.1100000000000003, 4.5, 4.6100000000000003, 5.1699999999999999, 4.5300000000000002, 5.3300000000000001, 5.1399999999999997, 4.8099999999999996, 4.1699999999999999, 4.4100000000000001, 3.5899999999999999, 5.8700000000000001, 3.8300000000000001, 6.0300000000000002, 4.8899999999999997, 4.3200000000000003, 4.6900000000000004, 6.3099999999999996, 5.1200000000000001, 5.54, 5.5, 5.3700000000000001, 5.29, 4.9199999999999999, 6.1500000000000004, 5.7999999999999998, 5.2599999999999998};
const Vec pgg = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
const Vec ax = {1.3, 2.7000000000000002, 0.40000000000000002, 5.0999999999999996, 3.2999999999999998};
const Vec ay = {2.2000000000000002, 2.8999999999999999, 3, 2.5, 1.8999999999999999, 3.1000000000000001};
const Vec tx = {1, 2, 2, 3, 5, 7, 7, 8};
const Vec ty = {2, 3, 3, 4, 4, 4, 5, 6};
const Vec u1 = {2.2869999999999999, -1.1970000000000001, -0.69399999999999995, -0.41199999999999998, -0.97099999999999997, -0.94699999999999995, 0.748, -0.11700000000000001, 0.153, 2.1899999999999999, 0.35699999999999998, 2.7170000000000001, 2.2810000000000001, 0.32400000000000001, 1.8959999999999999, 0.46800000000000003, -0.89400000000000002, -0.307, -0.0050000000000000001, 0.98799999999999999, 0.83999999999999997, 0.70499999999999996, 1.306, -1.3879999999999999, 1.2729999999999999, 0.184, 0.752, 0.59199999999999997, -0.98299999999999998, -0.27600000000000002, -0.871, 0.71899999999999997, 0.111, -0.078, -0.41999999999999998, -0.56200000000000006, 0.998, -1.105, -0.14199999999999999, 0.315};
const Vec u2 = {2.4369999999999998, -1.399, -0.57099999999999995, -2.6230000000000002, -0.78200000000000003, -0.80300000000000005, 2.7010000000000001, 1.1819999999999999, 0.20100000000000001, 1.8620000000000001, -0.52500000000000002, -0.014999999999999999, 0.73399999999999999, 3.4140000000000001, 1.4470000000000001, 0.96199999999999997, -3.1360000000000001, 0.63700000000000001, 0.33200000000000002, -1.8, 0.153, 0.318, 1.087, 1.4099999999999999, 0.63800000000000001, 2.218, 1.538, 2.3069999999999999, 2.5209999999999999, 1.401, 0.86499999999999999, -1.845, -1.2310000000000001, -1.7330000000000001, -3.2789999999999999, -2.6520000000000001, -1.778, -1.115, -0.125, 4.8449999999999998, 0.68500000000000005, 0.0080000000000000002, 0.058000000000000003, -0.78700000000000003, -1.585};
const Vec big1 = {-0.59103110258436853, 0.026594369016167005, -1.516553097081865, -1.3626533492958086, 1.1784891560316197, -0.93415131967336484, 1.3236056462120183, 0.6249177899447137, -0.045722955706094734, -1.00412057567228, -0.82843323666395796, -0.34835172517298091, -1.5382933970870489, -0.25556524546495329, -1.1499450327075542, 0.012326967679467856, -0.22296954074762165, 0.88777164791130303, -0.59215527987674588, -0.65571811858088613, -0.6825176222150191, -0.015858192728274354, -0.44260478535504522, 0.35255749944881615, 0.073170582291403163, 0.007158800380645901, -0.18760011060225232, -0.76570064552538997, -0.22105682074495833, -0.98358858747019573, -1.1042840421509679, -0.93815021417836408, 0.67862424418257394, -1.5774978652948257, -0.86993845822564153, 0.48467704535710093, -0.18605269864439061, 1.545554700183664, -0.61138007012819984, -0.34775648733634534, -1.6365163055394207, 0.020381440100037087, 0.89174268016784064, -0.87274968325225077, 0.89005082983947781, -0.34387434612934781, -2.1867813683280022, 0.8800581814920434, 0.7238565561464253, 0.21985268026050242, 0.78987057400075966, -0.22999390007756476, -0.81850248262772096, 0.49973415728563192, 0.1591923480860179, 0.54262642817766871, -0.15664505418271715, 0.43879332064263904, 1.4878705952295634, 0.060165097795432131};
const Vec big2 = {-1.273519301898475, 3.5095395858109519, -0.18180442831108407, -2.92531105658645, 0.80806723546281323, 2.5402722252420777, -1.1864523307534225, -1.6128909001315503, -0.91181268355782319, 1.131602490883274, 0.68021409377282316, -0.18515052066939075, -1.1446452334256405, 0.3424051320819157, 1.6791928747056528, 0.23485977544149633, -1.033158201123568, 0.67942440622451983, -1.6013199750749589, 0.60234772311677931, -0.097160536390880944, 0.47324443515071923, -0.90852233199684029, -1.361377037372935, 3.3924134636731691, -0.90484006229896052, -1.9467931434109416, 0.75967680414573313, -1.2800013935829968, -2.2590476840227343, 1.8034933490887606, -1.5417980785870513, 1.4074049366775274, -0.81473198762043531, 0.76964268759068699, -0.52888632377902822, 1.9897999299076052, -1.7102277710627307, 2.119664906938989, -0.90326817578302676, -0.66250292829168356, 0.36545815947264737, 0.89517172732811523, -0.18006107096857388, -3.1046152143601295, 0.87843837801041436, 0.73408713044207841, -1.5098944412861153, 1.9044799362105933, 1.6858735818636279, 1.2637256447030683, 2.4309144782528267, 0.67149876362346783, -3.4483581545456676, -0.11874548367208004};

void generated_cases() {
  // R: kruskal.test(Ozone ~ Month, data = airquality)
  expect("kruskal", "kruskal.test(Ozone ~ Month, data = airquality)", kruskal_test(oz, ozm),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 29.266576306116939, {{"df", 4}}, 6.9007141185467839e-06, "", kNaN});
  // R: kruskal.test(count ~ spray, data = InsectSprays)
  expect("kruskal", "kruskal.test(count ~ spray, data = InsectSprays)", kruskal_test(ins, insg),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 54.691344622371446, {{"df", 5}}, 1.5108444394185103e-10, "", kNaN});
  // R: kruskal.test(mucx, mucg)
  expect("kruskal", "kruskal.test(mucx, mucg)", kruskal_test(mucx, mucg),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 0.77142857142857224, {{"df", 2}}, 0.67996477357889362, "", kNaN});
  // R: kruskal.test(c(1, 2, 2, 3, 3, 3, 4), c(5, 5, -1, -1, 2, 2, 2))
  expect("kruskal", "kruskal.test(c(1, 2, 2, 3, 3, 3, 4), c(5, 5, -1, -1, 2, 2, 2))", kruskal_test(Vec{1, 2, 2, 3, 3, 3, 4}, Vec{5, 5, -1, -1, 2, 2, 2}),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 4.37254901960784, {{"df", 2}}, 0.11233447099424955, "", kNaN});
  // R: kruskal.test(c(1.5, 2, 9), c(0, 1, 1))
  expect("kruskal", "kruskal.test(c(1.5, 2, 9), c(0, 1, 1))", kruskal_test(Vec{1.5, 2, 9}, Vec{0, 1, 1}),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 1.5, {{"df", 1}}, 0.22067136191984676, "", kNaN});
  // R: kruskal.test(weight ~ group, data = PlantGrowth)
  expect("kruskal", "kruskal.test(weight ~ group, data = PlantGrowth)", kruskal_test(pg, pgg),
         Exp{"Kruskal-Wallis rank sum test", "", "Kruskal-Wallis chi-squared", 7.9882287494437154, {{"df", 2}}, 0.01842375573147197, "", kNaN});
  // R: bartlett.test(count ~ spray, data = InsectSprays)
  expect("bartlett", "bartlett.test(count ~ spray, data = InsectSprays)", bartlett_test(ins, insg),
         Exp{"Bartlett test of homogeneity of variances", "", "Bartlett's K-squared", 25.959825320368687, {{"df", 5}}, 9.0851223329453131e-05, "", kNaN});
  // R: bartlett.test(Ozone ~ Month, data = airquality)
  expect("bartlett", "bartlett.test(Ozone ~ Month, data = airquality)", bartlett_test(oz, ozm),
         Exp{"Bartlett test of homogeneity of variances", "", "Bartlett's K-squared", 13.450927363290358, {{"df", 4}}, 0.0092702786812983162, "", kNaN});
  // R: bartlett.test(weight ~ group, data = PlantGrowth)
  expect("bartlett", "bartlett.test(weight ~ group, data = PlantGrowth)", bartlett_test(pg, pgg),
         Exp{"Bartlett test of homogeneity of variances", "", "Bartlett's K-squared", 2.8785737872360935, {{"df", 2}}, 0.23709677363455817, "", kNaN});
  // R: bartlett.test(c(1, 2, 4, 8, 3, 3.5), c(2, 2, 2, 1, 1, 1))
  expect("bartlett", "bartlett.test(c(1, 2, 4, 8, 3, 3.5), c(2, 2, 2, 1, 1, 1))", bartlett_test(Vec{1, 2, 4, 8, 3, 3.5}, Vec{2, 2, 2, 1, 1, 1}),
         Exp{"Bartlett test of homogeneity of variances", "", "Bartlett's K-squared", 0.52622177345778243, {{"df", 1}}, 0.46819965007668479, "", kNaN});
  // R: bartlett.test(mucx, mucg)
  expect("bartlett", "bartlett.test(mucx, mucg)", bartlett_test(mucx, mucg),
         Exp{"Bartlett test of homogeneity of variances", "", "Bartlett's K-squared", 3.2794144046012046, {{"df", 2}}, 0.19403684751681766, "", kNaN});
  // R: fligner.test(InsectSprays$count, InsectSprays$spray)
  expect("fligner", "fligner.test(InsectSprays$count, InsectSprays$spray)", fligner_test(ins, insg),
         Exp{"Fligner-Killeen test of homogeneity of variances", "", "Fligner-Killeen:med chi-squared", 14.482781038458608, {{"df", 5}}, 0.012816779189709193, "", kNaN});
  // R: fligner.test(Ozone ~ Month, data = airquality)
  expect("fligner", "fligner.test(Ozone ~ Month, data = airquality)", fligner_test(oz, ozm),
         Exp{"Fligner-Killeen test of homogeneity of variances", "", "Fligner-Killeen:med chi-squared", 19.34070748250301, {{"df", 4}}, 0.00067359297535243064, "", kNaN});
  // R: fligner.test(weight ~ group, data = PlantGrowth)
  expect("fligner", "fligner.test(weight ~ group, data = PlantGrowth)", fligner_test(pg, pgg),
         Exp{"Fligner-Killeen test of homogeneity of variances", "", "Fligner-Killeen:med chi-squared", 2.3499371311581991, {{"df", 2}}, 0.30882868738638652, "", kNaN});
  // R: fligner.test(mucx, mucg)
  expect("fligner", "fligner.test(mucx, mucg)", fligner_test(mucx, mucg),
         Exp{"Fligner-Killeen test of homogeneity of variances", "", "Fligner-Killeen:med chi-squared", 5.2091832088217354, {{"df", 2}}, 0.073933325076689232, "", kNaN});
  // R: fligner.test(c(1, 2, 3, 4, 5, 7), c(1, 1, 1, 1, 1, 2))
  expect("fligner", "fligner.test(c(1, 2, 3, 4, 5, 7), c(1, 1, 1, 1, 1, 2))", fligner_test(Vec{1, 2, 3, 4, 5, 7}, Vec{1, 1, 1, 1, 1, 2}),
         Exp{"Fligner-Killeen test of homogeneity of variances", "", "Fligner-Killeen:med chi-squared", 1.3222262418600783, {{"df", 1}}, 0.25019290041039738, "", kNaN});
  // R: oneway.test(count ~ spray, data = InsectSprays)
  expect("oneway", "oneway.test(count ~ spray, data = InsectSprays)", oneway_test(ins, insg, false),
         Exp{"One-way analysis of means (not assuming equal variances)", "", "F", 36.065443893577253, {{"num df", 5}, {"denom df", 30.04256050876738}}, 7.9993794556733478e-12, "", kNaN});
  // R: oneway.test(count ~ spray, data = InsectSprays, var.equal = TRUE)
  expect("oneway", "oneway.test(count ~ spray, data = InsectSprays, var.equal = TRUE)", oneway_test(ins, insg, true),
         Exp{"One-way analysis of means", "", "F", 34.702282055491708, {{"num df", 5}, {"denom df", 66}}, 3.1825837261451404e-17, "", kNaN});
  // R: oneway.test(Ozone ~ Month, data = airquality)
  expect("oneway", "oneway.test(Ozone ~ Month, data = airquality)", oneway_test(oz, ozm, false),
         Exp{"One-way analysis of means (not assuming equal variances)", "", "F", 8.026676183748096, {{"num df", 4}, {"denom df", 42.668201053387996}}, 6.4390842025272748e-05, "", kNaN});
  // R: oneway.test(Ozone ~ Month, data = airquality, var.equal = TRUE)
  expect("oneway", "oneway.test(Ozone ~ Month, data = airquality, var.equal = TRUE)", oneway_test(oz, ozm, true),
         Exp{"One-way analysis of means", "", "F", 8.5356065886138452, {{"num df", 4}, {"denom df", 111}}, 4.8270645341147376e-06, "", kNaN});
  // R: oneway.test(weight ~ group, data = PlantGrowth)
  expect("oneway", "oneway.test(weight ~ group, data = PlantGrowth)", oneway_test(pg, pgg, false),
         Exp{"One-way analysis of means (not assuming equal variances)", "", "F", 5.1809724081131892, {{"num df", 2}, {"denom df", 17.128418616644129}}, 0.017392821490170001, "", kNaN});
  // R: oneway.test(x ~ g, data.frame(x = mucx, g = factor(mucg)))
  expect("oneway", "oneway.test(x ~ g, data.frame(x = mucx, g = factor(mucg)))", oneway_test(mucx, mucg, false),
         Exp{"One-way analysis of means (not assuming equal variances)", "", "F", 0.39190547838578599, {{"num df", 2}, {"denom df", 5.5187015751634521}}, 0.69318547683735066, "", kNaN});
  // R: oneway.test(x ~ g, data.frame(x = c(1, 2, 3, 4, 5, 6, 7), g = factor(c(1, 1, 2, 2, 3, 3, 3))))
  expect("oneway", "oneway.test(x ~ g, data.frame(x = c(1, 2, 3, 4, 5, 6, 7), g = factor(c(1, 1, 2, 2, 3, 3, 3))))", oneway_test(Vec{1, 2, 3, 4, 5, 6, 7}, Vec{1, 1, 2, 2, 3, 3, 3}, false),
         Exp{"One-way analysis of means (not assuming equal variances)", "", "F", 13.687296416938112, {{"num df", 2}, {"denom df", 2.4820512820512821}}, 0.045646404288774985, "", kNaN});
  // R: oneway.test(x ~ g, data.frame(x = c(1, 2, 3, 5), g = factor(c(1, 1, 2, 2))), var.equal = TRUE)
  expect("oneway", "oneway.test(x ~ g, data.frame(x = c(1, 2, 3, 5), g = factor(c(1, 1, 2, 2))), var.equal = TRUE)", oneway_test(Vec{1, 2, 3, 5}, Vec{1, 1, 2, 2}, true),
         Exp{"One-way analysis of means", "", "F", 5, {{"num df", 1}, {"denom df", 2}}, 0.1548457452714834, "", kNaN});
  // R: friedman.test(RoundingTimes)
  expect("friedman", "friedman.test(RoundingTimes)", friedman_test(rt, rtg, rtb),
         Exp{"Friedman rank sum test", "", "Friedman chi-squared", 11.142857142857142, {{"df", 2}}, 0.003805040775511363, "", kNaN});
  // R: friedman.test(wb$x, wb$w, wb$t)
  expect("friedman", "friedman.test(wb$x, wb$w, wb$t)", friedman_test(wbx, wbw, wbt),
         Exp{"Friedman rank sum test", "", "Friedman chi-squared", 0.33333333333333331, {{"df", 1}}, 0.56370286165077299, "", kNaN});
  // R: friedman.test(dataFreq)
  expect("friedman", "friedman.test(dataFreq)", friedman_test(qy, qyg, qyb),
         Exp{"Friedman rank sum test", "", "Friedman chi-squared", 8.3283582089552244, {{"df", 4}}, 0.080264090581825515, "", kNaN});
  // R: friedman.test(c(1, NA, 3, 4, 5, 6, 2, 2, 1), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))
  expect("friedman", "friedman.test(c(1, NA, 3, 4, 5, 6, 2, 2, 1), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))", friedman_test(Vec{1, kNaN, 3, 4, 5, 6, 2, 2, 1}, Vec{1, 2, 3, 1, 2, 3, 1, 2, 3}, Vec{1, 1, 1, 2, 2, 2, 3, 3, 3}),
         Exp{"Friedman rank sum test", "", "Friedman chi-squared", 0.2857142857142857, {{"df", 2}}, 0.86687789975018159, "", kNaN});
  // R: friedman.test(c(3, 1, 2, 9, 7, 8), c(30, 10, 20, 30, 10, 20), c(2, 2, 2, 1, 1, 1))
  expect("friedman", "friedman.test(c(3, 1, 2, 9, 7, 8), c(30, 10, 20, 30, 10, 20), c(2, 2, 2, 1, 1, 1))", friedman_test(Vec{3, 1, 2, 9, 7, 8}, Vec{30, 10, 20, 30, 10, 20}, Vec{2, 2, 2, 1, 1, 1}),
         Exp{"Friedman rank sum test", "", "Friedman chi-squared", 4, {{"df", 2}}, 0.1353352832366127, "", kNaN});
  // R: quade.test(dataFreq)
  expect("quade", "quade.test(dataFreq)", quade_test(qy, qyg, qyb),
         Exp{"Quade test", "", "Quade F", 3.8292515841753727, {{"num df", 4}, {"denom df", 24}}, 0.015189020073274591, "", kNaN});
  // R: quade.test(RoundingTimes)
  expect("quade", "quade.test(RoundingTimes)", quade_test(rt, rtg, rtb),
         Exp{"Quade test", "", "Quade F", 7.0422651149297435, {{"num df", 2}, {"denom df", 42}}, 0.002304253635080359, "", kNaN});
  // R: quade.test(wb$x, wb$w, wb$t)
  expect("quade", "quade.test(wb$x, wb$w, wb$t)", quade_test(wbx, wbw, wbt),
         Exp{"Quade test", "", "Quade F", 1.2307692307692306, {{"num df", 1}, {"denom df", 2}}, 0.38278660015163246, "", kNaN});
  // R: quade.test(c(1, 2, 3, 4, 5, 7, 1, 1, 1), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))
  expect("quade", "quade.test(c(1, 2, 3, 4, 5, 7, 1, 1, 1), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))", quade_test(Vec{1, 2, 3, 4, 5, 7, 1, 1, 1}, Vec{1, 2, 3, 1, 2, 3, 1, 2, 3}, Vec{1, 1, 1, 2, 2, 2, 3, 3, 3}),
         Exp{"Quade test", "", "Quade F", 3.5714285714285721, {{"num df", 2}, {"denom df", 4}}, 0.12886259040105194, "", kNaN});
  // R: quade.test(c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2))
  expect("quade", "quade.test(c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2))", quade_test(Vec{1, 2, 3, 4, 5, 6}, Vec{1, 2, 3, 1, 2, 3}, Vec{1, 1, 1, 2, 2, 2}),
         Exp{"Quade test", "", "Quade F", kNaN, {{"num df", kNaN}, {"denom df", kNaN}}, 0.16666666666666666, "", kNaN});
  // R: quade.test(c(1, 1, 1, 4, 4, 4, 2, 2, 2), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))
  expect("quade", "quade.test(c(1, 1, 1, 4, 4, 4, 2, 2, 2), c(1, 2, 3, 1, 2, 3, 1, 2, 3), c(1, 1, 1, 2, 2, 2, 3, 3, 3))", quade_test(Vec{1, 1, 1, 4, 4, 4, 2, 2, 2}, Vec{1, 2, 3, 1, 2, 3, 1, 2, 3}, Vec{1, 1, 1, 2, 2, 2, 3, 3, 3}),
         Exp{"Quade test", "", "Quade F", kNaN, {{"num df", kNaN}, {"denom df", kNaN}}, 0.027777777777777776, "", kNaN});
  // R: mcnemar.test(matrix(c(794, 86, 150, 570), nrow = 2))
  expect("mcnemar", "mcnemar.test(matrix(c(794, 86, 150, 570), nrow = 2))", mcnemar_test(M(2, 2, {794, 150, 86, 570}), true),
         Exp{"McNemar's Chi-squared test with continuity correction", "", "McNemar's chi-squared", 16.817796610169491, {{"df", 1}}, 4.1145622813459449e-05, "", kNaN});
  // R: mcnemar.test(matrix(c(794, 86, 150, 570), nrow = 2), correct = FALSE)
  expect("mcnemar", "mcnemar.test(matrix(c(794, 86, 150, 570), nrow = 2), correct = FALSE)", mcnemar_test(M(2, 2, {794, 150, 86, 570}), false),
         Exp{"McNemar's Chi-squared test", "", "McNemar's chi-squared", 17.35593220338983, {{"df", 1}}, 3.0992934410452189e-05, "", kNaN});
  // R: mcnemar.test(matrix(c(10, 5, 5, 10), 2))
  expect("mcnemar", "mcnemar.test(matrix(c(10, 5, 5, 10), 2))", mcnemar_test(M(2, 2, {10, 5, 5, 10}), true),
         Exp{"McNemar's Chi-squared test", "", "McNemar's chi-squared", 0, {{"df", 1}}, 1, "", kNaN});
  // R: mcnemar.test(matrix(c(10, 0, 0, 10), 2))
  expect("mcnemar", "mcnemar.test(matrix(c(10, 0, 0, 10), 2))", mcnemar_test(M(2, 2, {10, 0, 0, 10}), true),
         Exp{"McNemar's Chi-squared test", "", "McNemar's chi-squared", kNaN, {{"df", 1}}, kNaN, "", kNaN});
  // R: mcnemar.test(matrix(c(10, 5, 3, 10, 1, 2, 4, 4, 4), 3))
  expect("mcnemar", "mcnemar.test(matrix(c(10, 5, 3, 10, 1, 2, 4, 4, 4), 3))", mcnemar_test(M(3, 3, {10, 10, 4, 5, 1, 4, 3, 2, 4}), true),
         Exp{"McNemar's Chi-squared test", "", "McNemar's chi-squared", 2.4761904761904763, {{"df", 3}}, 0.4796093688043172, "", kNaN});
  // R: mcnemar.test(matrix(c(1, 2.5, 3, 4), 2))
  expect("mcnemar", "mcnemar.test(matrix(c(1, 2.5, 3, 4), 2))", mcnemar_test(M(2, 2, {1, 3, 2.5, 4}), true),
         Exp{"McNemar's Chi-squared test with continuity correction", "", "McNemar's chi-squared", 0.045454545454545456, {{"df", 1}}, 0.83117040954176236, "", kNaN});
  // R: mcnemar.test(matrix(c(12, 3, 1, 7, 20, 2, 9, 5, 30, 4, 6, 8, 2, 1, 3, 15), 4))
  expect("mcnemar", "mcnemar.test(matrix(c(12, 3, 1, 7, 20, 2, 9, 5, 30, 4, 6, 8, 2, 1, 3, 15), 4))", mcnemar_test(M(4, 4, {12, 20, 30, 2, 3, 2, 4, 1, 1, 9, 6, 3, 7, 5, 8, 15}), true),
         Exp{"McNemar's Chi-squared test", "", "McNemar's chi-squared", 49.334498289617507, {{"df", 6}}, 6.3905344805301044e-09, "", kNaN});
  // R: ansari.test(ramsay, jung.parekh)
  expect("ansari", "ansari.test(ramsay, jung.parekh)", ansari_test(ramsay, jung),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 185.5, {}, 0.18145819972867069, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "two.sided")
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"two.sided\")", ansari_test(ax, ay, O().alt("two.sided")),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 12, {}, 0.16883116883116883, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "two.sided", exact = FALSE)
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"two.sided\", exact = FALSE)", ansari_test(ax, ay, O().alt("two.sided").exact(0)),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 12, {}, 0.11549160949968168, "ratio of scales", 1});
  // R: ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = "two.sided")
  expect("ansari", "ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = \"two.sided\")", ansari_test(tx, ty, O().alt("two.sided")),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 24.5, {}, 0.021622223939692077, "ratio of scales", 1});
  // R: ansari.test(u1, u2, alternative = "two.sided")
  expect("ansari", "ansari.test(u1, u2, alternative = \"two.sided\")", ansari_test(u1, u2, O().alt("two.sided")),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 1014.5, {}, 0.011023485353862661, "ratio of scales", 1});
  // R: ansari.test(big1, big2, alternative = "two.sided")
  expect("ansari", "ansari.test(big1, big2, alternative = \"two.sided\")", ansari_test(big1, big2, O().alt("two.sided")),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 2107, {}, 8.1571370481370309e-05, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "less")
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"less\")", ansari_test(ax, ay, O().alt("less")),
         Exp{"Ansari-Bradley test", "less", "AB", 12, {}, 0.96103896103896103, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "less", exact = FALSE)
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"less\", exact = FALSE)", ansari_test(ax, ay, O().alt("less").exact(0)),
         Exp{"Ansari-Bradley test", "less", "AB", 12, {}, 0.94225419525015919, "ratio of scales", 1});
  // R: ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = "less")
  expect("ansari", "ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = \"less\")", ansari_test(tx, ty, O().alt("less")),
         Exp{"Ansari-Bradley test", "less", "AB", 24.5, {}, 0.98918888803015392, "ratio of scales", 1});
  // R: ansari.test(u1, u2, alternative = "less")
  expect("ansari", "ansari.test(u1, u2, alternative = \"less\")", ansari_test(u1, u2, O().alt("less")),
         Exp{"Ansari-Bradley test", "less", "AB", 1014.5, {}, 0.0055117426769313305, "ratio of scales", 1});
  // R: ansari.test(big1, big2, alternative = "less")
  expect("ansari", "ansari.test(big1, big2, alternative = \"less\")", ansari_test(big1, big2, O().alt("less")),
         Exp{"Ansari-Bradley test", "less", "AB", 2107, {}, 4.0785685240685154e-05, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "greater")
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"greater\")", ansari_test(ax, ay, O().alt("greater")),
         Exp{"Ansari-Bradley test", "greater", "AB", 12, {}, 0.084415584415584416, "ratio of scales", 1});
  // R: ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = "greater", exact = FALSE)
  expect("ansari", "ansari.test(c(1.3, 2.7, 0.4, 5.1, 3.3), c(2.2, 2.9, 3.0, 2.5, 1.9, 3.1), alternative = \"greater\", exact = FALSE)", ansari_test(ax, ay, O().alt("greater").exact(0)),
         Exp{"Ansari-Bradley test", "greater", "AB", 12, {}, 0.05774580474984084, "ratio of scales", 1});
  // R: ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = "greater")
  expect("ansari", "ansari.test(c(1, 2, 2, 3, 5, 7, 7, 8), c(2, 3, 3, 4, 4, 4, 5, 6), alternative = \"greater\")", ansari_test(tx, ty, O().alt("greater")),
         Exp{"Ansari-Bradley test", "greater", "AB", 24.5, {}, 0.010811111969846039, "ratio of scales", 1});
  // R: ansari.test(u1, u2, alternative = "greater")
  expect("ansari", "ansari.test(u1, u2, alternative = \"greater\")", ansari_test(u1, u2, O().alt("greater")),
         Exp{"Ansari-Bradley test", "greater", "AB", 1014.5, {}, 0.99448825732306867, "ratio of scales", 1});
  // R: ansari.test(big1, big2, alternative = "greater")
  expect("ansari", "ansari.test(big1, big2, alternative = \"greater\")", ansari_test(big1, big2, O().alt("greater")),
         Exp{"Ansari-Bradley test", "greater", "AB", 2107, {}, 0.99995921431475931, "ratio of scales", 1});
  // R: ansari.test(u2, u1)
  expect("ansari", "ansari.test(u2, u1)", ansari_test(u2, u1),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 834.5, {}, 0.01102348535386272, "ratio of scales", 1});
  // R: ansari.test(u1[1:7], u2[1:6])
  expect("ansari", "ansari.test(u1[1:7], u2[1:6])", ansari_test(head(u1, 7), head(u2, 6)),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 27, {}, 0.97552447552447563, "ratio of scales", 1});
  // R: ansari.test(c(1, 2, 3), c(1.5, 2.5))
  expect("ansari", "ansari.test(c(1, 2, 3), c(1.5, 2.5))", ansari_test(Vec{1, 2, 3}, Vec{1.5, 2.5}),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 5, {}, 1, "ratio of scales", 1});
  // R: ansari.test(c(1, Inf, 3), c(1.5, 2.5, NaN))
  expect("ansari", "ansari.test(c(1, Inf, 3), c(1.5, 2.5, NaN))", ansari_test(Vec{1, kInf, 3}, Vec{1.5, 2.5, kNaN}),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 4, {}, 0.40000000000000002, "ratio of scales", 1});
  // R: ansari.test(c(4, 1, 3, 2), c(2.5, 0))
  expect("ansari", "ansari.test(c(4, 1, 3, 2), c(2.5, 0))", ansari_test(Vec{4, 1, 3, 2}, Vec{2.5, 0}),
         Exp{"Ansari-Bradley test", "two.sided", "AB", 8, {}, 1, "ratio of scales", 1});
}
// END GENERATED

void error_cases() {
  // R: kruskal.test(c(1, 2, 3), c(1, 1, 1)) -> "all observations are in the same group"
  check::throws([] { kruskal_test(Vec{1, 2, 3}, Vec{1, 1, 1}); }, "kruskal.test: one group");
  check::throws([] { kruskal_test(Vec{1, 2, 3}, Vec{1, 2}); }, "kruskal.test: lengths differ");
  check::throws([] { kruskal_test(Vec{}, Vec{}); }, "kruskal.test: empty");
  // R: kruskal.test(c(1, 2, 3, 4), c(1, 1, NA, NA)) -> one group after dropping NA
  check::throws([] { kruskal_test(Vec{1, 2, 3, 4}, Vec{1, 1, kNaN, kNaN}); }, "kruskal.test: one group after NaN");
  check::throws([] { kruskal_test(Vec{1, 2, kNaN, kNaN}, Vec{1, 1, 2, 2}); }, "kruskal.test: one group after NaN x");
  // R: bartlett.test(c(1, 2, 3, 4, 5), c(1, 1, 2, 2, 3)) -> "there must be at least 2 observations in each group"
  check::throws([] { bartlett_test(Vec{1, 2, 3, 4, 5}, Vec{1, 1, 2, 2, 3}); }, "bartlett.test: singleton group");
  check::throws([] { bartlett_test(Vec{1, 2, 3}, Vec{1, 1, 1}); }, "bartlett.test: one group");
  check::throws([] { fligner_test(Vec{1, 2, 3}, Vec{1, 1, 1}); }, "fligner.test: one group");
  check::throws([] { fligner_test(Vec{kNaN, kNaN}, Vec{1, 2}); }, "fligner.test: all NaN");
  // R: oneway.test(x ~ g, ...) with a group of one -> "not enough observations"; one group -> "not enough groups"
  check::throws([] { oneway_test(Vec{1, 2, 3, 4, 5}, Vec{1, 1, 2, 2, 3}, false); }, "oneway.test: singleton (Welch)");
  check::throws([] { oneway_test(Vec{1, 2, 3, 4, 5}, Vec{1, 1, 2, 2, 3}, true); }, "oneway.test: singleton (F)");
  check::throws([] { oneway_test(Vec{1, 2, 3}, Vec{1, 1, 1}, true); }, "oneway.test: one group");
  // R: friedman.test(c(1, 2, 3, 4), c(1, 2, 1, 1), c(1, 1, 2, 2)) -> "not an unreplicated complete block design"
  check::throws([] { friedman_test(Vec{1, 2, 3, 4}, Vec{1, 2, 1, 1}, Vec{1, 1, 2, 2}); }, "friedman.test: replicated");
  check::throws([] { friedman_test(Vec{1, 2, 3}, Vec{1, 2, 1}, Vec{1, 1, 2}); }, "friedman.test: incomplete");
  check::throws([] { friedman_test(Vec{1, 2, 3}, Vec{1, 2}, Vec{1, 1, 1}); }, "friedman.test: lengths differ");
  check::throws([] { friedman_test(Vec{1, 2, 3, 4}, Vec{1, 2, 1, kNaN}, Vec{1, 1, 2, 2}); }, "friedman.test: NaN group");
  check::throws([] { friedman_test(Vec{1, 2, 3}, Vec{1, 1, 1}, Vec{1, 2, 3}); }, "friedman.test: one group");
  check::throws([] { friedman_test(Vec{1, 2, 3}, Vec{1, 2, 3}, Vec{1, 1, 1}); }, "friedman.test: one block");
  // R errors here too (one block left after removing the block with NA)
  check::throws([] { friedman_test(Vec{1, 2, kNaN, 4}, Vec{1, 2, 1, 2}, Vec{1, 1, 2, 2}); }, "friedman.test: one complete block");
  check::throws([] { quade_test(Vec{1, 2, 3, 4}, Vec{1, 2, 1, 1}, Vec{1, 1, 2, 2}); }, "quade.test: replicated");
  check::throws([] { quade_test(Vec{1, 2, kNaN, 4}, Vec{1, 2, 1, 2}, Vec{1, 1, 2, 2}); }, "quade.test: one complete block");
  // R: mcnemar.test(matrix(1:6, 2)) -> "'x' must be square with at least two rows and columns"
  check::throws([] { mcnemar_test(M(2, 3, {1, 2, 3, 4, 5, 6}), true); }, "mcnemar.test: not square");
  check::throws([] { mcnemar_test(M(1, 1, {1}), true); }, "mcnemar.test: 1 x 1");
  // R: mcnemar.test(matrix(c(1, -2, 3, 4), 2)) -> "all entries of 'x' must be nonnegative and finite"
  check::throws([] { mcnemar_test(M(2, 2, {1, 3, -2, 4}), true); }, "mcnemar.test: negative");
  check::throws([] { mcnemar_test(M(2, 2, {1, 3, kNaN, 4}), true); }, "mcnemar.test: NaN");
  // R: ansari.test(numeric(0), 1:2) -> "not enough 'x' observations"
  check::throws([] { ansari_test(Vec{}, Vec{1, 2}); }, "ansari.test: empty x");
  check::throws([] { ansari_test(Vec{1, 2}, Vec{kNaN, kNaN}); }, "ansari.test: all y NaN");
  check::throws([] { ansari_test(Vec{1, 2}, Vec{3, 4}, O().alt("sideways")); }, "ansari.test: bad alternative");
  check::throws([] { ansari_test(Vec{1, 2}, Vec{3, 4}, O().ci()); }, "ansari.test: conf.int unsupported");
}

void corner_cases() {
  {  // R: kruskal.test(c(1, 1, 1, 1), c(1, 1, 2, 2)) -> statistic NaN, df 1, p NaN (all tied)
    const HtestResult r = kruskal_test(Vec{1, 1, 1, 1}, Vec{1, 1, 2, 2});
    check::ok(std::isnan(r.statistic) && std::isnan(r.p_value) && r.parameters[0].second == 1, "kruskal all ties: NaN");
  }
  {  // R: kruskal.test(c(1, 2, 3, NaN), c(1, 2, 3, 4)) -> 2, df 2, p 0.36787944117144233 (group 4 dropped)
    const HtestResult r = kruskal_test(Vec{1, 2, 3, kNaN}, Vec{1, 2, 3, 4});
    check::close(r.statistic, 2.0, 1e-15, "kruskal NaN drop: statistic");
    check::close(r.parameters[0].second, 2.0, 0, "kruskal NaN drop: df");
    check::close(r.p_value, 0.36787944117144233, 1e-14, "kruskal NaN drop: p");
  }
  {  // R: kruskal.test(c(1, Inf, 3, 4), c(1, 1, 2, 2)) -> 0, p 1 (Inf is kept and ranked)
    const HtestResult r = kruskal_test(Vec{1, kInf, 3, 4}, Vec{1, 1, 2, 2});
    check::close(r.statistic, 0.0, 1e-15, "kruskal Inf kept: statistic");
    check::close(r.p_value, 1.0, 1e-15, "kruskal Inf kept: p");
  }
  {  // R: bartlett.test(c(1, 1, 3, 4), c(1, 1, 2, 2)) -> Inf, p 0
    const HtestResult r = bartlett_test(Vec{1, 1, 3, 4}, Vec{1, 1, 2, 2});
    check::ok(r.statistic == kInf && r.p_value == 0.0, "bartlett zero variance: Inf, p 0");
  }
  {  // R: bartlett.test(c(1, NaN, 3, 4, 5, 7), c(1, 1, 1, 2, 2, 2)) -> 0.0056017348696526525, p 0.94033824550260336
    const HtestResult r = bartlett_test(Vec{1, kNaN, 3, 4, 5, 7}, Vec{1, 1, 1, 2, 2, 2});
    check::close(r.statistic, 0.0056017348696526525, 1e-12, "bartlett NaN drop: statistic");
    check::close(r.p_value, 0.94033824550260336, 1e-12, "bartlett NaN drop: p");
  }
  {  // R: fligner.test(c(1, 2, 3, 4, 5, 6), c(1, 1, 2, 2, 3, 3)) -> NaN (all scores equal)
    const HtestResult r = fligner_test(Vec{1, 2, 3, 4, 5, 6}, Vec{1, 1, 2, 2, 3, 3});
    check::ok(std::isnan(r.statistic) && std::isnan(r.p_value), "fligner constant scores: NaN");
  }
  {  // R: oneway.test(x ~ g, data.frame(x = c(1, 1, 3, 4, 5, 6), g = factor(c(1, 1, 2, 2, 3, 3)))) -> NaN, df (2, NaN)
    const HtestResult r = oneway_test(Vec{1, 1, 3, 4, 5, 6}, Vec{1, 1, 2, 2, 3, 3}, false);
    check::ok(std::isnan(r.statistic) && std::isnan(r.p_value) && r.parameters[0].second == 2 &&
                  std::isnan(r.parameters[1].second), "oneway Welch zero variance: NaN");
  }
  {  // R: oneway.test(..., var.equal = TRUE) on c(1, 1, 3, 4, 5, 6) -> 30.500000000000004, df (2, 3), p 0.010148735200598885
    const HtestResult r = oneway_test(Vec{1, 1, 3, 4, 5, 6}, Vec{1, 1, 2, 2, 3, 3}, true);
    check::close(r.statistic, 30.500000000000004, 1e-13, "oneway F zero variance group: F");
    check::close(r.p_value, 0.010148735200598885, 1e-12, "oneway F zero variance group: p");
  }
  {  // R: friedman.test(c(1, 1, 3, 3), c(1, 2, 1, 2), c(1, 1, 2, 2)) -> NaN (every block tied)
    const HtestResult r = friedman_test(Vec{1, 1, 3, 3}, Vec{1, 2, 1, 2}, Vec{1, 1, 2, 2});
    check::ok(std::isnan(r.statistic) && std::isnan(r.p_value), "friedman all ties: NaN");
  }
  {  // groups given in any order and with arbitrary labels: same as the sorted design
    const HtestResult a = friedman_test(Vec{9, 7, 8, 3, 1, 2}, Vec{-1, 5, 2.5, -1, 5, 2.5}, Vec{7, 7, 7, 0, 0, 0});
    const HtestResult b = friedman_test(Vec{3, 1, 2, 9, 7, 8}, Vec{1, 2, 3, 1, 2, 3}, Vec{1, 1, 1, 2, 2, 2});
    check::close(a.statistic, b.statistic, 0, "friedman label order invariance");
  }
  {  // R: mcnemar.test(matrix(c(10, 5, 4, 10), 2)) -> corrected, statistic 0, p 1
    const HtestResult r = mcnemar_test(M(2, 2, {10, 4, 5, 10}), true);
    check::ok(r.method == "McNemar's Chi-squared test with continuity correction", "mcnemar |b - c| = 1: corrected");
    check::close(r.statistic, 0.0, 0, "mcnemar |b - c| = 1: statistic");
    check::close(r.p_value, 1.0, 0, "mcnemar |b - c| = 1: p");
  }
  {  // ansari.test(big1, big2, exact = TRUE): AB = 2107. R gives p 6.2266638376051731e-05; the exact rational value
     // 2 P(AB >= 2107) (big-integer enumeration of the m-subset score sums, C(115, 60) total) is
     // 6.226663837674315e-05, so R is off by 1.1e-11 relative here and we check against the exact value.
    const HtestResult r = ansari_test(big1, big2, O().exact(1));
    check::close(r.statistic, 2107, 0, "ansari exact large: AB");
    check::rel(r.p_value, 6.226663837674315e-05, 1e-14, "ansari exact large: p (exact rational)");
    check::rel(r.p_value, 6.2266638376051731e-05, 2e-11, "ansari exact large: p (R)");
  }
}

}  // namespace

int main() {
  error_cases();
  corner_cases();
  generated_cases();
  if (std::getenv("HTEST_VERBOSE"))
    for (const auto& kv : worst) std::fprintf(stderr, "  max rel dev %-22s %.3g\n", kv.first.c_str(), kv.second);
  return check::report("htest2");
}
