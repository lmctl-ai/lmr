// Evaluates fuzz cases read from stdin (one per line) with the native core; prints one result per line.
// Line formats: "dist NAME p1 p2 p3 FN X LOWER LOG"   FN in d,p,q
//               "quant TYPE P x1 x2 ..."  "cor METHOD n x... y..."  "test KIND ALT n x... [m y...]"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include "rcli/dist.hpp"
#include "rcli/htest.hpp"
#include "rcli/summary.hpp"
using namespace rcli;
// Tokens are parsed with strtod: istream >> double rejects "nan" and "inf".
static double num(std::istringstream& in) { std::string t; in >> t; return std::strtod(t.c_str(), nullptr); }
static Vec read_vec(std::istringstream& in) { size_t n; in >> n; Vec v(n); for (auto& x : v) x = num(in); return v; }
int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string kind;
    in >> kind;
    try {
      if (kind == "dist") {
        std::string name, fn; int lower, lg;
        in >> name;
        double a = num(in), b = num(in), c = num(in);
        in >> fn;
        double x = num(in);
        in >> lower >> lg;
        Dist d = Dist::by_name(name, a, b, c);
        double r = fn == "d" ? d.d(x, lg) : fn == "p" ? d.p(x, lower, lg) : d.q(x, lower, lg);
        std::printf("%.17g\n", r);
      } else if (kind == "quant") {
        int type; in >> type; double p = num(in);
        std::printf("%.17g\n", summary::quantile(read_vec(in), p, type));
      } else if (kind == "cor") {
        std::string m; in >> m; Vec x = read_vec(in), y = read_vec(in);
        std::printf("%.17g\n", summary::cor(x, y, m));
      } else if (kind == "test") {
        std::string k, alt; in >> k >> alt;
        TestOptions o; o.alternative = alt;
        Vec x = read_vec(in);
        HtestResult r;
        if (k == "t1") r = t_test(x, o);
        else if (k == "welch") r = t_test(x, read_vec(in), o);
        else if (k == "w1") r = wilcox_test(x, o);
        else if (k == "w2") r = wilcox_test(x, read_vec(in), o);
        else if (k == "kendall" || k == "spearman" || k == "pearson") { o.method = k; r = cor_test(x, read_vec(in), o); }
        std::printf("%.17g %.17g\n", r.statistic, r.p_value);
      }
    } catch (const std::exception& e) {
      std::printf("ERROR\n");
    }
    std::fflush(stdout);
  }
}
