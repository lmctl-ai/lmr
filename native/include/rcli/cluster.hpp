// Clustering: k-means (Hartigan & Wong 1979 AS 136, Lloyd, MacQueen) and agglomerative hierarchical clustering
// (Lance-Williams updates). Data are row-major matrices: rows are observations, columns are variables.
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli {

class Rng;

struct KmeansOptions {
  int nstart = 1;
  int iter_max = 10;
  std::string algorithm = "Hartigan-Wong";  // or "Lloyd", "MacQueen" ("Forgy" is an alias of Lloyd)
};
struct KmeansFit {
  std::vector<int> cluster;   // 1-based cluster per row (R convention)
  Matrix centers;             // k x p
  Vec withinss;               // per cluster
  std::vector<int> size;
  double totss = kNaN, tot_withinss = kNaN, betweenss = kNaN;
  int iter = 0;
  int ifault = 0;             // 0 ok, 2 did not converge in iter_max, 4 Hartigan-Wong ran out of quick-transfer steps
};
// Initial centres: k distinct rows chosen with rng (different stream from R, so labels may be permuted vs R).
KmeansFit kmeans(const Matrix& x, int k, Rng& rng, const KmeansOptions& opt = {});
// Explicit initial centres (k x p): deterministic, for testing against R with the same centres.
KmeansFit kmeans(const Matrix& x, const Matrix& centers, const KmeansOptions& opt = {});

// Distances between rows, R dist() order: d[(i, j)] for i < j stored column-wise like R's lower triangle,
// index = n*i - i*(i+1)/2 + j - i - 1 (0-based i < j). method: "euclidean", "maximum", "manhattan", "canberra",
// "binary", "minkowski" (with p).
Vec dist(const Matrix& x, const std::string& method = "euclidean", double p = 2.0);

struct HclustFit {
  std::vector<int> merge_a, merge_b;  // R $merge: negative = singleton observation (1-based), positive = earlier step
  Vec height;                         // n - 1 merge heights
  std::vector<int> order;             // 1-based leaf order for a dendrogram without crossings (R $order)
  std::string method;
  size_t n = 0;
  std::vector<int> cutree(int k) const;  // 1-based group per observation, groups numbered by first appearance
};
// method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
// Ties between equal distances are broken as R does (the pair found first in its scan order) where documented;
// otherwise merge order with exact ties may differ, heights never do.
HclustFit hclust(const Vec& d, size_t n, const std::string& method = "complete");

}  // namespace rcli
