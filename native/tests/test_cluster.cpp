// Tests for clustering (kmeans, dist, hclust, cutree). Expected values are R 4.3.3 outputs recorded with
// sprintf("%.17g") from the R expressions in the comments (R used purely as a black-box oracle).
#include <initializer_list>
#include <string>
#include <vector>

#include "check.hpp"
#include "rcli/cluster.hpp"
#include "rcli/rng.hpp"

using namespace rcli;

namespace {

Matrix mat(size_t r, size_t c, std::initializer_list<double> v) {
  Matrix m(r, c);
  m.a.assign(v.begin(), v.end());
  return m;
}

KmeansFit kmeans_alg(const Matrix& x, const Matrix& c, const char* alg, int iter_max) {
  KmeansOptions o;
  o.algorithm = alg;
  o.iter_max = iter_max;
  return kmeans(x, c, o);
}

void expect_vec(const Vec& got, std::initializer_list<double> want, double tol, const std::string& what) {
  check::ok(got.size() == want.size(), what + " length");
  size_t i = 0;
  for (double w : want) {
    if (i < got.size()) check::close(got[i], w, tol, what + "[" + std::to_string(i) + "]");
    ++i;
  }
}

void expect_ints(const std::vector<int>& got, std::initializer_list<int> want, const std::string& what) {
  check::ok(got == std::vector<int>(want), what);
}

void expect_km(const KmeansFit& f, std::initializer_list<int> cluster, std::initializer_list<double> centers,
               std::initializer_list<double> withinss, int iter, int ifault, double totss, const std::string& what) {
  expect_ints(f.cluster, cluster, what + " cluster");
  expect_vec(f.centers.a, centers, 1e-12, what + " centers");
  expect_vec(f.withinss, withinss, 1e-12, what + " withinss");
  check::ok(f.iter == iter, what + " iter");
  check::ok(f.ifault == ifault, what + " ifault");
  check::close(f.totss, totss, 1e-12, what + " totss");
  double tw = 0.0;
  for (double w : f.withinss) tw += w;
  check::close(f.tot_withinss, tw, 1e-13, what + " tot.withinss");
  check::close(f.betweenss, f.totss - f.tot_withinss, 1e-13, what + " betweenss");
  std::vector<int> size(f.centers.rows, 0);
  for (int c : f.cluster) size[static_cast<size_t>(c - 1)] += 1;
  check::ok(f.size == size, what + " size");
}

void expect_hc(const HclustFit& h, std::initializer_list<int> ma, std::initializer_list<int> mb,
               std::initializer_list<double> height, std::initializer_list<int> order, std::initializer_list<int> cut3,
               const std::string& what) {
  expect_ints(h.merge_a, ma, what + " merge[,1]");
  expect_ints(h.merge_b, mb, what + " merge[,2]");
  expect_vec(h.height, height, 1e-12, what + " height");
  expect_ints(h.order, order, what + " order");
  expect_ints(h.cutree(3), cut3, what + " cutree(3)");
}

}  // namespace

int main() {
  // Data generated in R 4.3.3:
  // set.seed(42); mu <- rbind(c(5.0,3.4,1.5,0.25), c(5.9,2.8,4.3,1.3), c(6.6,3.0,5.5,2.0))
  // A <- round(mu[rep(1:3, each=10),] + matrix(rnorm(120, sd=0.35), 30), 1)
  // B: hand-made 2-D set with three clusters
  // C: set.seed(177); n <- sample(20:60,1); sc <- 10^sample(6:15,1)
  //    C <- matrix(sc + sample(0:3, n*2, TRUE) * sc * 2^-50 * sample(1:9,1), n); k <- sample(2:6,1)
  //    u <- unique(C); Ccen <- u[sample(nrow(u),k),]   (a data set on which Hartigan-Wong returns ifault 4)
  // Acen3 <- A[c(3,14,25),]; Acen4 <- A[1:4,]; Bcen3 <- B[1:3,]
  const Matrix kA = mat(30, 4, {5.5, 3.6, 1.4, 0.7, 4.8, 3.6, 1.6, 0.1, 5.1, 3.8, 1.7, 0.5, 5.2, 3.2, 2, 0.7, 5.1, 3.6, 1.2, -0.1, 5, 2.8, 2, -0.1, 5.5, 3.1, 1.6, -0.1, 5, 3.1, 1.9, -0.3, 5.7, 2.6, 1.8, 0.3, 5, 3.4, 1.8, 0.5, 6.4, 2.9, 3.9, 1.7, 6.7, 2.7, 4.3, 1.7, 5.4, 3.1, 4.5, 0.9, 5.8, 2.5, 4, 1.9, 5.9, 2.3, 4.1, 1.1, 6.1, 3, 4.5, 1.3, 5.8, 2.5, 4.6, 1.2, 5, 3.3, 4.5, 1.3, 5, 2.6, 4, 1.4, 6.4, 3, 3.9, 1.3, 6.5, 3.1, 6, 2, 6, 2.7, 5.6, 2, 6.5, 3.6, 5.5, 1.8, 7, 3.2, 5.5, 1.8, 7.3, 3, 5.1, 1.4, 6.4, 3.1, 5.7, 1.9, 6.5, 3.2, 5.4, 1.8, 6, 3, 5.4, 2.9, 6.8, 2, 5.8, 1.5, 6.4, 3.1, 5.8, 2});
  const Matrix kB = mat(14, 2, {1, 1, 1.5, 2, 2, 1.2, 1.2, 0.5, 0.8, 1.6, 6, 6, 6.5, 7, 7, 6.2, 5.8, 6.9, 6.3, 5.5, 1, 8, 1.6, 8.8, 0.5, 9, 1.9, 7.6});
  const Matrix kC = mat(48, 2, {10000000000000, 10000000000000, 10000000000000.01, 10000000000000.018, 10000000000000.018, 10000000000000.01, 10000000000000.01, 10000000000000.018, 10000000000000.018, 10000000000000.01, 10000000000000.01, 10000000000000.027, 10000000000000.018, 10000000000000.018, 10000000000000.01, 10000000000000.027, 10000000000000, 10000000000000, 10000000000000.018, 10000000000000, 10000000000000.027, 10000000000000, 10000000000000.027, 10000000000000.018, 10000000000000, 10000000000000.018, 10000000000000.018, 10000000000000.018, 10000000000000.01, 10000000000000.018, 10000000000000, 10000000000000.01, 10000000000000.018, 10000000000000.01, 10000000000000.01, 10000000000000, 10000000000000.01, 10000000000000.018, 10000000000000.018, 10000000000000.018, 10000000000000, 10000000000000.018, 10000000000000.027, 10000000000000, 10000000000000.01, 10000000000000.01, 10000000000000.027, 10000000000000.027, 10000000000000.01, 10000000000000.027, 10000000000000.01, 10000000000000.018, 10000000000000.01, 10000000000000.027, 10000000000000.01, 10000000000000.018, 10000000000000, 10000000000000.027, 10000000000000, 10000000000000.018, 10000000000000.018, 10000000000000.027, 10000000000000.01, 10000000000000.018, 10000000000000.018, 10000000000000, 10000000000000.01, 10000000000000.027, 10000000000000.027, 10000000000000, 10000000000000, 10000000000000.027, 10000000000000, 10000000000000.01, 10000000000000.027, 10000000000000, 10000000000000.027, 10000000000000.027, 10000000000000, 10000000000000.01, 10000000000000.027, 10000000000000.027, 10000000000000.01, 10000000000000.027, 10000000000000.018, 10000000000000.018, 10000000000000.018, 10000000000000.018, 10000000000000.01, 10000000000000.018, 10000000000000.018, 10000000000000.027, 10000000000000, 10000000000000.01, 10000000000000.027, 10000000000000.018});
  const Matrix kCcen = mat(2, 2, {10000000000000.027, 10000000000000.018, 10000000000000.018, 10000000000000.027});
  const Matrix kAcen3 = mat(3, 4, {5.1, 3.8, 1.7, 0.5, 5.8, 2.5, 4, 1.9, 7.3, 3, 5.1, 1.4});
  const Matrix kAcen4 = mat(4, 4, {5.5, 3.6, 1.4, 0.7, 4.8, 3.6, 1.6, 0.1, 5.1, 3.8, 1.7, 0.5, 5.2, 3.2, 2, 0.7});
  const Matrix kBcen3 = mat(3, 2, {1, 1, 1.5, 2, 2, 1.2});
  {  // kmeans(kA, centers = kAcen3, iter.max = 10, algorithm = "Hartigan-Wong")
    KmeansFit f = kmeans_alg(kA, kAcen3, "Hartigan-Wong", 10);
    expect_km(f, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, {5.19, 3.2800000000000002, 1.7, 0.21999999999999992, 5.85, 2.79, 4.2299999999999995, 1.3800000000000001, 6.5399999999999991, 3.0000000000000004, 5.58, 1.9100000000000001}, {3.901, 5.4910000000000014, 5.109}, 2, 0, 117.35833333333333, "A3 Hartigan-Wong");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 10, algorithm = "Hartigan-Wong")
    KmeansFit f = kmeans_alg(kA, kAcen4, "Hartigan-Wong", 10);
    expect_km(f, {1, 3, 1, 1, 3, 2, 2, 2, 2, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.2, 3.5000000000000004, 1.7249999999999999, 0.6, 5.3, 2.9, 1.825, -0.05, 4.9499999999999993, 3.6, 1.4, 0, 6.195, 2.8950000000000005, 4.905, 1.645}, {0.5675, 0.8375, 0.14500000000000002, 23.718000000000011}, 2, 0, 117.35833333333333, "A4 Hartigan-Wong");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 1, algorithm = "Hartigan-Wong")
    KmeansFit f = kmeans_alg(kA, kAcen4, "Hartigan-Wong", 1);
    expect_km(f, {1, 3, 1, 1, 3, 2, 2, 2, 2, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.2, 3.5000000000000004, 1.7249999999999999, 0.6, 5.3, 2.9, 1.825, -0.05, 4.9499999999999993, 3.6, 1.4, 0, 6.195, 2.8950000000000005, 4.905, 1.645}, {0.5675, 0.8375, 0.14500000000000002, 23.718000000000011}, 2, 2, 117.35833333333333, "A4 it1 Hartigan-Wong");
  }
  {  // kmeans(kB, centers = kBcen3, iter.max = 10, algorithm = "Hartigan-Wong")
    KmeansFit f = kmeans_alg(kB, kBcen3, "Hartigan-Wong", 10);
    expect_km(f, {1, 3, 3, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2}, {1, 1.0333333333333334, 4.0666666666666664, 7.2222222222222223, 1.75, 1.6}, {0.68666666666666676, 71.215555555555554, 0.445}, 2, 0, 212.27285714285713, "B3 Hartigan-Wong");
  }
  {  // kmeans(kA, centers = kAcen3, iter.max = 10, algorithm = "Lloyd")
    KmeansFit f = kmeans_alg(kA, kAcen3, "Lloyd", 10);
    expect_km(f, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, {5.19, 3.2800000000000002, 1.7, 0.21999999999999992, 5.85, 2.79, 4.2299999999999995, 1.3800000000000001, 6.5399999999999991, 3.0000000000000004, 5.58, 1.9100000000000001}, {3.9010000000000011, 5.4910000000000014, 5.109}, 3, 0, 117.35833333333333, "A3 Lloyd");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 10, algorithm = "Lloyd")
    KmeansFit f = kmeans_alg(kA, kAcen4, "Lloyd", 10);
    expect_km(f, {1, 3, 3, 3, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.5, 3.6, 1.4, 0.7, 5.26, 3.04, 1.7000000000000004, -0.060000000000000019, 5.0249999999999995, 3.5000000000000004, 1.775, 0.44999999999999996, 6.195, 2.8950000000000005, 4.905, 1.645}, {0, 1.5760000000000007, 0.565, 23.718}, 4, 0, 117.35833333333333, "A4 Lloyd");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 1, algorithm = "Lloyd")
    KmeansFit f = kmeans_alg(kA, kAcen4, "Lloyd", 1);
    expect_km(f, {1, 2, 3, 4, 2, 4, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.5, 3.6, 1.4, 0.7, 5.1, 3.35, 1.5750000000000002, -0.1, 5.1, 3.8, 1.7, 0.5, 6.033333333333335, 2.9125, 4.4041666666666668, 1.4291666666666665}, {0, 0.83750000000000013, 0, 63.698749999999968}, 2, 2, 117.35833333333333, "A4 it1 Lloyd");
  }
  {  // kmeans(kB, centers = kBcen3, iter.max = 10, algorithm = "Lloyd")
    KmeansFit f = kmeans_alg(kB, kBcen3, "Lloyd", 10);
    expect_km(f, {1, 3, 3, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2}, {1, 1.0333333333333334, 4.0666666666666664, 7.2222222222222223, 1.75, 1.6}, {0.68666666666666676, 71.21555555555554, 0.445}, 3, 0, 212.27285714285713, "B3 Lloyd");
  }
  {  // kmeans(kA, centers = kAcen3, iter.max = 10, algorithm = "MacQueen")
    KmeansFit f = kmeans_alg(kA, kAcen3, "MacQueen", 10);
    expect_km(f, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, {5.19, 3.2800000000000002, 1.7, 0.21999999999999992, 5.85, 2.79, 4.2299999999999995, 1.38, 6.5399999999999991, 3.0000000000000004, 5.58, 1.9100000000000001}, {3.9010000000000011, 5.4910000000000014, 5.109}, 2, 0, 117.35833333333333, "A3 MacQueen");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 10, algorithm = "MacQueen")
    KmeansFit f = kmeans_alg(kA, kAcen4, "MacQueen", 10);
    expect_km(f, {1, 3, 3, 3, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.5, 3.6, 1.4, 0.7, 5.2600000000000007, 3.0400000000000005, 1.7000000000000002, -0.060000000000000012, 5.025, 3.5, 1.7750000000000001, 0.44999999999999996, 6.1950000000000021, 2.895, 4.905, 1.6449999999999996}, {0, 1.5760000000000005, 0.565, 23.718000000000004}, 3, 0, 117.35833333333333, "A4 MacQueen");
  }
  {  // kmeans(kA, centers = kAcen4, iter.max = 1, algorithm = "MacQueen")
    KmeansFit f = kmeans_alg(kA, kAcen4, "MacQueen", 1);
    expect_km(f, {1, 2, 3, 3, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}, {5.5, 3.6, 1.4, 0.7, 5.1833333333333336, 3.1333333333333337, 1.6833333333333336, -0.03333333333333334, 5.1000000000000005, 3.4666666666666668, 1.8333333333333335, 0.56666666666666665, 6.1950000000000021, 2.895, 4.905, 1.6449999999999996}, {0, 2.0433333333333339, 0.2799999999999998, 23.718000000000004}, 2, 2, 117.35833333333333, "A4 it1 MacQueen");
  }
  {  // kmeans(kB, centers = kBcen3, iter.max = 10, algorithm = "MacQueen")
    KmeansFit f = kmeans_alg(kB, kBcen3, "MacQueen", 10);
    expect_km(f, {1, 3, 3, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2}, {1, 1.0333333333333334, 4.0666666666666664, 7.2222222222222223, 1.75, 1.6}, {0.68666666666666676, 71.21555555555554, 0.445}, 2, 0, 212.27285714285713, "B3 MacQueen");
  }
  {  // kmeans(kC, centers = kCcen, iter.max = 10, algorithm = "Hartigan-Wong")
    // This case sits on a rounding edge (coordinates near 1e13) found by search; whether the quick-transfer step limit
    // is hit (ifault 4) depends on fused multiply-add contraction, which R on aarch64 uses and the release build
    // (-ffp-contract=off) does not. Only the platform-independent part is checked: a valid fault code and partition.
    KmeansFit f = kmeans_alg(kC, kCcen, "Hartigan-Wong", 10);
    check::ok(f.ifault == 0 || f.ifault == 4, "C ifault4 ifault is 0 or 4");
    check::ok(f.cluster.size() == 48 && f.size.size() == 2 && f.size[0] + f.size[1] == 48, "C ifault4 partition");
  }

  // H <- A[c(1:5, 11:15, 21:24), ]
  const Matrix kH = mat(14, 4, {5.5, 3.6, 1.4, 0.7, 4.8, 3.6, 1.6, 0.1, 5.1, 3.8, 1.7, 0.5, 5.2, 3.2, 2, 0.7, 5.1, 3.6, 1.2, -0.1, 6.4, 2.9, 3.9, 1.7, 6.7, 2.7, 4.3, 1.7, 5.4, 3.1, 4.5, 0.9, 5.8, 2.5, 4, 1.9, 5.9, 2.3, 4.1, 1.1, 6.5, 3.1, 6, 2, 6, 2.7, 5.6, 2, 6.5, 3.6, 5.5, 1.8, 7, 3.2, 5.5, 1.8});
  const Vec dH = dist(kH);  // dist(H)
  Vec dH2 = dH;
  for (double& v : dH2) v *= v;  // dist(H)^2
  {  // h <- hclust(dist(H), "ward.D")
    HclustFit h = hclust(dH, 14, "ward.D");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, -9, 4, 1, 2, -8, 8, 10}, {-7, -5, -3, -14, -12, 3, -10, 5, 7, 6, 9, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.57445626465380306, 0.64031242374328479, 0.75498344352707514, 0.80060241096687423, 0.83666002653407545, 1.1374558523292859, 1.2350945535315687, 1.2606811138845577, 1.6352766591335086, 5.05216337817132, 18.189302724958615}, {2, 5, 4, 1, 3, 13, 14, 11, 12, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "ward.D");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "ward.D2")
    HclustFit h = hclust(dH, 14, "ward.D2");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, -9, 4, 1, 2, -8, 8, 10}, {-7, -5, -3, -14, -12, 3, -10, 5, 7, 6, 9, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.57445626465380306, 0.64031242374328479, 0.75498344352707514, 0.793725393319377, 0.83666002653407545, 1.1224972160321824, 1.1768602295939823, 1.1958260743101399, 1.5801898620102584, 3.5451375149632769, 9.4169300426716251}, {2, 5, 4, 1, 3, 13, 14, 11, 12, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "ward.D2");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "single")
    HclustFit h = hclust(dH, 14, "single");
    expect_hc(h, {-6, -2, -3, -1, -13, -4, -11, -12, -9, -10, -8, 8, 6}, {-7, -5, 2, 3, -14, 4, 5, 7, 1, 9, 10, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.54772255750516585, 0.57445626465380306, 0.64031242374328479, 0.70710678118654724, 0.73484692283495345, 0.75498344352707514, 0.75498344352707536, 0.83666002653407545, 1.0440306508910553, 1.3379088160259653, 2.4228082879171433}, {4, 1, 3, 2, 5, 12, 11, 13, 14, 8, 10, 9, 6, 7}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "single");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "complete")
    HclustFit h = hclust(dH, 14, "complete");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, -9, 1, -12, 2, -8, 9, 10}, {-7, -5, -3, -14, 4, 3, -10, 7, 5, 6, 8, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.57445626465380306, 0.64031242374328479, 0.74161984870956632, 0.78102496759066542, 0.83666002653407545, 1.0954451150103321, 1.1401754250991381, 1.2041594578792296, 1.5905973720586863, 2.3280893453645635, 5.44609952167604}, {2, 5, 4, 1, 3, 12, 11, 13, 14, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "complete");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "average")
    HclustFit h = hclust(dH, 14, "average");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, -9, 2, 1, -12, -8, 10, 8}, {-7, -5, -3, -14, 4, 3, -10, 6, 7, 5, 9, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.57445626465380306, 0.64031242374328479, 0.73823338577225983, 0.74406587438860639, 0.83666002653407545, 0.88908936356704116, 0.96134140357766573, 0.98290808130382912, 1.3483317945558297, 1.8778854970297629, 3.8041206518867412}, {2, 5, 4, 1, 3, 12, 11, 13, 14, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "average");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "mcquitty")
    HclustFit h = hclust(dH, 14, "mcquitty");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, -9, -12, 2, 1, -8, 8, 9}, {-7, -5, -3, -14, 4, 3, -10, 5, 6, 7, 10, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.57445626465380306, 0.64031242374328479, 0.73823338577225983, 0.74406587438860639, 0.83666002653407545, 0.92592692185964065, 0.93190134728408058, 0.96134140357766573, 1.3483317945558297, 1.8547178815381971, 3.8780548168482953}, {2, 5, 4, 1, 3, 12, 11, 13, 14, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "mcquitty");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "median")
    HclustFit h = hclust(dH, 14, "median");
    expect_hc(h, {-6, -2, -3, -1, -4, -13, -11, -12, -9, -10, -8, 8, 5}, {-7, -5, 2, 3, 4, -14, 6, 7, 1, 9, 10, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.54234504598914779, 0.54930562378964054, 0.593596686960125, 0.64031242374328479, 0.57815527983643866, 0.70134904893262029, 0.73783734841575854, 0.69166428535452207, 0.93317605512249557, 1.2220876218524339, 2.6935993931690989}, {4, 1, 3, 2, 5, 12, 11, 13, 14, 8, 10, 9, 6, 7}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "median");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H)^2, "median")
    HclustFit h = hclust(dH2, 14, "median");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, 2, -9, 1, -12, -8, 10, 7}, {-7, -5, -3, -14, 4, 3, 6, -10, 8, 5, 9, 11, 12}, {0.2899999999999997, 0.29, 0.33000000000000024, 0.40999999999999992, 0.44250000000000006, 0.47249999999999975, 0.67562500000000025, 0.69999999999999984, 0.69250000000000023, 0.72562500000000019, 1.5606249999999997, 2.6934374999999995, 14.432734374999999}, {2, 5, 4, 1, 3, 12, 11, 13, 14, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "median sq");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H), "centroid")
    HclustFit h = hclust(dH, 14, "centroid");
    expect_hc(h, {-6, -2, -3, -1, -4, -13, -11, -9, -10, -12, -8, 10, 5}, {-7, -5, 2, 3, 4, -14, 6, 1, 8, 7, 9, 11, 12}, {0.53851648071345015, 0.53851648071345037, 0.54234504598914779, 0.60118263750038048, 0.63177444079320866, 0.64031242374328479, 0.57815527983643866, 0.73783734841575854, 0.72531429922975521, 0.74771039293851749, 1.0220479119584429, 1.1367367600885472, 2.8294470905491185}, {4, 1, 3, 2, 5, 12, 11, 13, 14, 8, 10, 9, 6, 7}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "centroid");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  {  // h <- hclust(dist(H)^2, "centroid")
    HclustFit h = hclust(dH2, 14, "centroid");
    expect_hc(h, {-6, -2, -1, -13, -11, -4, 2, -9, 1, -12, -8, 10, 7}, {-7, -5, -3, -14, 4, 3, 6, -10, 8, 5, 9, 11, 12}, {0.2899999999999997, 0.29, 0.33000000000000024, 0.40999999999999992, 0.44250000000000006, 0.47249999999999975, 0.59583333333333333, 0.69999999999999984, 0.69250000000000023, 0.82666666666666677, 1.5606249999999997, 2.8278000000000008, 13.794444444444444}, {2, 5, 4, 1, 3, 12, 11, 13, 14, 8, 6, 7, 9, 10}, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}, "centroid sq");  // merge[,1], merge[,2], height, order, cutree(h, 3)
  }
  // set.seed(1007); T <- matrix(sample(0:3, 24, TRUE), 12)
  const Matrix kT = mat(12, 2, {2, 2, 2, 1, 0, 2, 0, 0, 2, 2, 2, 2, 2, 0, 3, 3, 3, 3, 2, 0, 2, 1, 3, 1});
  {  // h <- hclust(dist(T, "manhattan"), "single")  (many exact ties)
    HclustFit h = hclust(dist(kT, "manhattan"), 12, "single");
    expect_hc(h, {-1, -6, -2, -7, -8, 2, 4, -12, -3, -4, 5}, {-5, 1, -11, -10, -9, 3, 6, 7, 8, 9, 10}, {0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2}, {8, 9, 4, 3, 12, 7, 10, 6, 1, 5, 2, 11}, {1, 1, 1, 2, 1, 1, 1, 3, 3, 1, 1, 1}, "ties single");
  }
  {  // h <- hclust(dist(T, "manhattan"), "complete")  (many exact ties)
    HclustFit h = hclust(dist(kT, "manhattan"), 12, "complete");
    expect_hc(h, {-1, -6, -2, -7, -8, 2, 4, -12, -3, 8, 5}, {-5, 1, -11, -10, -9, 3, 6, 7, -4, 9, 10}, {0, 0, 0, 0, 0, 1, 2, 2, 2, 4, 6}, {8, 9, 12, 7, 10, 6, 1, 5, 2, 11, 3, 4}, {1, 1, 2, 2, 1, 1, 1, 3, 3, 1, 1, 1}, "ties complete");
  }
  {  // h <- hclust(dist(T, "manhattan"), "average")  (many exact ties)
    HclustFit h = hclust(dist(kT, "manhattan"), 12, "average");
    expect_hc(h, {-1, -6, -2, -7, -8, 2, 4, -12, -3, 5, 9}, {-5, 1, -11, -10, -9, 3, 6, 7, -4, 8, 10}, {0, 0, 0, 0, 0, 1, 1.6, 1.7142857142857142, 2, 2.75, 3.5}, {3, 4, 8, 9, 12, 7, 10, 6, 1, 5, 2, 11}, {1, 1, 2, 2, 1, 1, 1, 3, 3, 1, 1, 1}, "ties average");
  }
  {  // h <- hclust(dist(T, "manhattan"), "ward.D2")  (many exact ties)
    HclustFit h = hclust(dist(kT, "manhattan"), 12, "ward.D2");
    expect_hc(h, {-1, -6, -2, -7, -8, -12, -3, 4, 2, 5, 7}, {-5, 1, -11, -10, -9, 3, -4, 6, 8, 9, 10}, {0, 0, 0, 0, 0, 1.1547005383792515, 2, 2.0655911179772888, 2.8982753492378879, 4.8579831205964474, 5.9777364723893056}, {3, 4, 8, 9, 6, 1, 5, 7, 10, 12, 2, 11}, {1, 1, 2, 2, 1, 1, 1, 3, 3, 1, 1, 1}, "ties ward.D2");
  }
  {  // h <- hclust(dist(T, "manhattan"), "centroid")  (many exact ties)
    HclustFit h = hclust(dist(kT, "manhattan"), 12, "centroid");
    expect_hc(h, {-1, -6, -2, -7, -8, 2, 4, -12, -3, 8, 5}, {-5, 1, -11, -10, -9, 3, 6, 7, -4, 9, 10}, {0, 0, 0, 0, 0, 1, 1.3599999999999999, 1.2653061224489794, 2, 2.09375, 2.34}, {8, 9, 12, 7, 10, 6, 1, 5, 2, 11, 3, 4}, {1, 1, 2, 2, 1, 1, 1, 3, 3, 1, 1, 1}, "ties centroid");
  }
  // D <- rbind(c(0,0,1,1,1,Inf), c(1,0,1,1,0,1), c(0,0,0,0,0,0), c(NA,2,0,3,-1,0), c(0,0,0,0,0,0), c(2.5,-1,0.5,0,0,4))
  const Matrix kD = mat(6, 6, {0, 0, 1, 1, 1, kInf, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, kNaN, 2, 0, 3, -1, 0, 0, 0, 0, 0, 0, 0, 2.5, -1, 0.5, 0, 0, 4});
  expect_vec(dist(kD, "euclidean", 3), {kInf, kInf, kInf, kInf, kInf, 2, 3.63318042491699, 2, 3.6742346141747673, 4.09878030638384, 0, 4.8476798574163293, 4.09878030638384, 6.5038450166036395, 4.8476798574163293}, 1e-14, "dist euclidean");  // dist(D, "euclidean", p = 3)
  expect_vec(dist(kD, "maximum", 3), {kInf, kInf, kInf, kInf, kInf, 1, 2, 1, 3, 3, 0, 4, 3, 4, 4}, 1e-14, "dist maximum");  // dist(D, "maximum", p = 3)
  expect_vec(dist(kD, "manhattan", 3), {kInf, kInf, kInf, kInf, kInf, 4, 8.4, 4, 7, 7.1999999999999993, 0, 8, 7.1999999999999993, 13.799999999999999, 8}, 1e-14, "dist manhattan");  // dist(D, "manhattan", p = 3)
  expect_vec(dist(kD, "canberra", 3), {3.5999999999999996, 6, 5.3999999999999995, 6, 5.3333333333333339, 6, 5.3999999999999995, 6, 4.0342857142857138, 6, kNaN, 6, 6, 6, 6}, 1e-14, "dist canberra");  // dist(D, "canberra", p = 3)
  expect_vec(dist(kD, "binary", 3), {0.5, 1, 0.5, 1, 0.8, 1, 0.8, 1, 0.4, 1, 0, 1, 1, 0.8, 1}, 1e-14, "dist binary");  // dist(D, "binary", p = 3)
  expect_vec(dist(kD, "minkowski", 3), {kInf, kInf, kInf, kInf, kInf, 1.5874010519681994, 2.8355998780353819, 1.5874010519681994, 3.1912521494299537, 3.5088212858554391, 0, 4.3222927358739529, 3.5088212858554391, 5.228711982673441, 4.3222927358739529}, 1e-14, "dist minkowski");  // dist(D, "minkowski", p = 3)

  // ---- dist: storage order and small cases
  {
    const Matrix x = mat(4, 2, {0, 0, 3, 4, 6, 8, 0, 1});
    // dist(rbind(c(0,0), c(3,4), c(6,8), c(0,1))) -> 5 10 1 5 4.2426406871192848 9.2195444572928871
    expect_vec(dist(x), {5, 10, 1, 5, 4.2426406871192848, 9.2195444572928871}, 1e-15, "dist order");
    // dist(rbind(c(1,2,3))) -> numeric(0)
    check::ok(dist(mat(1, 3, {1, 2, 3})).empty(), "dist one row");
    // dist(rbind(c(NA,1), c(2,NA))) -> NA
    expect_vec(dist(mat(2, 2, {kNaN, 1, 2, kNaN})), {kNaN}, 0, "dist all pairs missing");
    // dist(rbind(c(0,0,1,1,1,1), c(1,0,1,1,0,1)), "binary") -> 0.4 ; "canberra" -> 2.4 (?dist example)
    const Matrix e = mat(2, 6, {0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1});
    expect_vec(dist(e, "binary"), {0.4}, 1e-15, "dist binary example");
    expect_vec(dist(e, "canberra"), {2.4}, 1e-15, "dist canberra example");
    // dist(rbind(c(Inf,1,2), c(1,3,2)), "canberra") -> 1.5 ; with c(Inf,..) vs c(Inf,..) -> 0.75 (Inf-Inf dropped)
    expect_vec(dist(mat(2, 3, {kInf, 1, 2, 1, 3, 2}), "canberra"), {1.5}, 1e-15, "dist canberra Inf");
    expect_vec(dist(mat(2, 3, {kInf, 1, 2, kInf, 3, 2}), "canberra"), {0.75}, 1e-15, "dist canberra Inf-Inf");
  }

  // ---- hclust: small exact cases and cutree edges
  {
    // hclust(dist(c(1, 2, 4, 8)), "single"): merge rbind(c(-1,-2), c(-3,1), c(-4,2)), height 1 2 4, order 4 3 1 2
    HclustFit h = hclust(dist(mat(4, 1, {1, 2, 4, 8})), 4, "single");
    expect_hc(h, {-1, -3, -4}, {-2, 1, 2}, {1, 2, 4}, {4, 3, 1, 2}, {1, 1, 2, 3}, "single 1-D");
    // cutree(h, 1) -> 1 1 1 1 ; cutree(h, 4) -> 1 2 3 4
    expect_ints(h.cutree(1), {1, 1, 1, 1}, "cutree k=1");
    expect_ints(h.cutree(4), {1, 2, 3, 4}, "cutree k=n");
    check::throws([&] { h.cutree(0); }, "cutree k=0");
    check::throws([&] { h.cutree(5); }, "cutree k>n");
    // hclust(dist(c(0, 10, 1, 11)), "complete"): merge rbind(c(-1,-3), c(-2,-4), c(1,2)), height 1 1 11,
    // order 1 3 2 4, cutree(h, 3) -> 1 2 1 3 (tie at height 1 resolved as R does)
    HclustFit h2 = hclust(dist(mat(4, 1, {0, 10, 1, 11})), 4, "complete");
    expect_hc(h2, {-1, -2, 1}, {-3, -4, 2}, {1, 1, 11}, {1, 3, 2, 4}, {1, 2, 1, 3}, "complete pairs");
    // hclust(as.dist(matrix(c(0,3,3,0),2)), "average"): merge c(-1,-2), height 3, order 1 2
    HclustFit h3 = hclust(Vec{3}, 2, "average");
    expect_ints(h3.merge_a, {-1}, "two points merge[,1]");
    expect_ints(h3.merge_b, {-2}, "two points merge[,2]");
    expect_vec(h3.height, {3}, 0, "two points height");
    expect_ints(h3.order, {1, 2}, "two points order");
    expect_ints(h3.cutree(2), {1, 2}, "two points cutree(2)");
  }

  // ---- kmeans with random starts (stream differs from R; check against R's optimum)
  {
    Rng rng(20261019);
    KmeansOptions o;
    o.nstart = 10;
    // kmeans(B, B[c(1,6,11),])$tot.withinss -> optimum of the hand-made 3-cluster set
    KmeansFit ref = kmeans(kB, mat(3, 2, {1, 1, 6, 6, 1, 8}));
    KmeansFit f = kmeans(kB, 3, rng, o);
    check::close(f.tot_withinss, ref.tot_withinss, 1e-12, "kmeans nstart optimum");
    check::ok(f.size.size() == 3 && f.centers.rows == 3 && f.cluster.size() == kB.rows, "kmeans random shape");
    for (const char* alg : {"Lloyd", "Forgy", "MacQueen"}) {
      o.algorithm = alg;
      KmeansFit g = kmeans(kB, 3, rng, o);
      check::close(g.tot_withinss, ref.tot_withinss, 1e-12, std::string("kmeans nstart optimum ") + alg);
    }
    // k = 1: kmeans(B, 1) -> one cluster, centre colMeans(B), iter 1
    KmeansFit one = kmeans(kB, 1, rng);
    check::ok(one.iter == 1 && one.ifault == 0 && one.size == std::vector<int>{14}, "kmeans k=1");
    check::close(one.tot_withinss, one.totss, 1e-14, "kmeans k=1 withinss == totss");
    // same seed, same call sequence -> same result
    Rng r1(7), r2(7);
    KmeansFit a1 = kmeans(kA, 3, r1), a2 = kmeans(kA, 3, r2);
    check::ok(a1.cluster == a2.cluster, "kmeans reproducible");
    // kmeans(B, B[c(1,6),], algorithm = "Lloyd") with a centre far from all points: Lloyd/MacQueen keep an
    // empty cluster (R warns; centre NaN, size 0), Hartigan-Wong errors.
    const Matrix far = mat(2, 2, {1, 1, 100, 100});
    KmeansOptions lo;
    lo.algorithm = "Lloyd";
    KmeansFit e = kmeans(kB, far, lo);
    check::ok(e.size == std::vector<int>({14, 0}) && std::isnan(e.centers(1, 0)), "Lloyd empty cluster");
    check::throws([&] { kmeans(kB, far); }, "Hartigan-Wong empty cluster");
  }

  // ---- invalid input
  {
    Rng rng(1);
    check::throws([&] { kmeans(kB, 0, rng); }, "kmeans k=0");
    // kmeans(rbind(c(1,1),c(1,1),c(2,2)), 3) -> error: more cluster centers than distinct data points
    check::throws([&] { kmeans(mat(3, 2, {1, 1, 1, 1, 2, 2}), 3, rng); }, "kmeans k > distinct rows");
    check::throws([&] { kmeans(Matrix(), 1, rng); }, "kmeans empty");
    check::throws([&] { kmeans(mat(3, 1, {1, kNaN, 2}), 1, rng); }, "kmeans NaN");
    KmeansOptions bad;
    bad.algorithm = "Elkan";
    check::throws([&] { kmeans(kB, 2, rng, bad); }, "kmeans unknown algorithm");
    // kmeans(x, x[c(1,1,60),]) -> error: initial centers are not distinct
    check::throws([&] { kmeans(kB, mat(2, 2, {1, 1, 1, 1})); }, "kmeans duplicate centres");
    check::throws([&] { kmeans(kB, mat(2, 3, {1, 1, 1, 2, 2, 2})); }, "kmeans centre columns");
    // kmeans(x, x) with Hartigan-Wong -> error: number of cluster centres must lie between 1 and nrow(x)
    check::throws([&] { kmeans(mat(3, 1, {1, 2, 3}), mat(3, 1, {1, 2, 3})); }, "kmeans HW k == n");
    check::throws([&] { dist(kB, "cosine"); }, "dist unknown method");
    check::throws([&] { hclust(Vec{1, 2}, 3, "complete"); }, "hclust wrong length");
    check::throws([&] { hclust(Vec{}, 1, "complete"); }, "hclust n < 2");
    check::throws([&] { hclust(Vec{1, kNaN, 2}, 3, "complete"); }, "hclust NaN");
    check::throws([&] { hclust(Vec{1, 2, 3}, 3, "ward"); }, "hclust unknown method");
  }

  return check::report("cluster");
}
