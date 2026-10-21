// k-means and hierarchical clustering over the numeric columns of a Frame (rows = observations).
#include "host.hpp"
#include "rcli/cluster.hpp"
#include "rcli/rng.hpp"

namespace rcli::host {

namespace {

struct KmeansModel {
  KmeansFit fit;
  Strs columns;
};

Matrix frame_matrix(const Frame& f, const char* what, Strs* names = nullptr) {
  for (size_t j = 0; j < f.cols.size(); ++j)
    if (f.cols[j].text)
      throw dslvm::Error(std::string(what) + ": column '" + f.names[j] + "' is text; select the numeric columns first (Frame.select)");
  if (f.cols.empty() || f.nrow() == 0) throw dslvm::Error(std::string(what) + ": the frame is empty");
  Matrix m(f.nrow(), f.cols.size());
  for (size_t j = 0; j < f.cols.size(); ++j)
    for (size_t i = 0; i < f.nrow(); ++i) {
      double v = f.cols[j].num[i];
      if (std::isnan(v)) throw dslvm::Error(std::string(what) + ": column '" + f.names[j] + "' has missing values (Frame.drop_nan)");
      m(i, j) = v;
    }
  if (names) *names = f.names;
  return m;
}

size_t label_index(int64_t label, size_t k, const char* what) {
  if (label < 1 || static_cast<uint64_t>(label) > k)
    throw dslvm::Error(std::string(what) + ": cluster label " + std::to_string(label) + " is outside 1.." + std::to_string(k));
  return static_cast<size_t>(label - 1);
}

}  // namespace

void bind_cluster(Binder& b) {
  Host& h = b.host;
  using A = Args;
  auto vec = [&h](Vec v) { return h.hold("Vec", std::move(v)); };

  auto run = [&h](const Frame& f, int64_t k, Rng& rng, int64_t nstart, int64_t iter_max, const std::string& algorithm) {
    if (nstart < 1 || nstart > 10000) throw dslvm::Error("kmeans: nstart must be between 1 and 10000");
    if (iter_max < 1 || iter_max > 100000) throw dslvm::Error("kmeans: iter_max must be between 1 and 100000");
    if (k < 1) throw dslvm::Error("kmeans: k must be at least 1");
    KmeansModel m;
    Matrix x = frame_matrix(f, "kmeans", &m.columns);
    KmeansOptions o;
    o.nstart = static_cast<int>(nstart);
    o.iter_max = static_cast<int>(iter_max);
    o.algorithm = algorithm;
    m.fit = kmeans(x, static_cast<int>(k), rng, o);
    return h.hold("Kmeans", std::move(m));
  };
  b.fn("kmeans", {"Frame", "int", "Rng"}, "Kmeans", [run](A& a) { return run(obj<Frame>(a[0]), a[1].i, obj<Rng>(a[2]), 1, 10, "Hartigan-Wong"); });
  b.fn("kmeans", {"Frame", "int", "Rng", "int"}, "Kmeans", [run](A& a) {
    return run(obj<Frame>(a[0]), a[1].i, obj<Rng>(a[2]), a[3].i, 10, "Hartigan-Wong");
  });
  b.fn("kmeans", {"Frame", "int", "Rng", "int", "int", "string"}, "Kmeans", [run](A& a) {
    return run(obj<Frame>(a[0]), a[1].i, obj<Rng>(a[2]), a[3].i, a[4].i, a[5].s);
  });
  b.m("Kmeans", "k", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<KmeansModel>(s).fit.size.size())); });
  b.m("Kmeans", "n", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<KmeansModel>(s).fit.cluster.size())); });
  b.m("Kmeans", "cluster", {"int"}, "int", [](void* s, A& a) {
    auto& c = self<KmeansModel>(s).fit.cluster;
    return Value::integer(c[static_cast<size_t>(as_index(a[0].i, c.size(), "Kmeans.cluster"))]);
  });
  b.m("Kmeans", "clusters", {}, "Vec", [vec](void* s, A&) {
    Vec out;
    for (int c : self<KmeansModel>(s).fit.cluster) out.push_back(c);
    return vec(out);
  });
  b.m("Kmeans", "size", {"int"}, "int", [](void* s, A& a) {
    auto& f = self<KmeansModel>(s).fit;
    return Value::integer(f.size[label_index(a[0].i, f.size.size(), "Kmeans.size")]);
  });
  b.m("Kmeans", "withinss", {"int"}, "float", [](void* s, A& a) {
    auto& f = self<KmeansModel>(s).fit;
    return Value::real(f.withinss[label_index(a[0].i, f.withinss.size(), "Kmeans.withinss")]);
  });
  b.m("Kmeans", "center", {"int", "int"}, "float", [](void* s, A& a) {
    auto& f = self<KmeansModel>(s).fit;
    size_t i = label_index(a[0].i, f.centers.rows, "Kmeans.center");
    return Value::real(f.centers(i, static_cast<size_t>(as_index(a[1].i, f.centers.cols, "Kmeans.center column"))));
  });
  b.m("Kmeans", "center", {"int", "string"}, "float", [](void* s, A& a) {
    auto& m = self<KmeansModel>(s);
    size_t i = label_index(a[0].i, m.fit.centers.rows, "Kmeans.center");
    for (size_t j = 0; j < m.columns.size(); ++j) if (m.columns[j] == a[1].s) return Value::real(m.fit.centers(i, j));
    throw dslvm::Error("Kmeans.center: no column '" + a[1].s + "'");
  });
  b.m("Kmeans", "tot_withinss", {}, "float", [](void* s, A&) { return Value::real(self<KmeansModel>(s).fit.tot_withinss); });
  b.m("Kmeans", "betweenss", {}, "float", [](void* s, A&) { return Value::real(self<KmeansModel>(s).fit.betweenss); });
  b.m("Kmeans", "totss", {}, "float", [](void* s, A&) { return Value::real(self<KmeansModel>(s).fit.totss); });
  b.m("Kmeans", "iter", {}, "int", [](void* s, A&) { return Value::integer(self<KmeansModel>(s).fit.iter); });
  b.m("Kmeans", "ifault", {}, "int", [](void* s, A&) { return Value::integer(self<KmeansModel>(s).fit.ifault); });

  auto hc = [&h](const Frame& f, const std::string& method, const std::string& dist_method, double p) {
    Matrix x = frame_matrix(f, "hclust");
    if (x.rows < 2) throw dslvm::Error("hclust: need at least 2 rows");
    Vec d = dist(x, dist_method, p);
    if (method == "ward.D2" || method == "centroid" || method == "median") {
      // These methods are applied to whatever dissimilarities are given, as in R; ward.D2 squares them internally.
    }
    return h.hold("Hclust", hclust(d, x.rows, method));
  };
  b.fn("hclust", {"Frame", "string"}, "Hclust", [hc](A& a) { return hc(obj<Frame>(a[0]), a[1].s, "euclidean", 2.0); });
  b.fn("hclust", {"Frame", "string", "string"}, "Hclust", [hc](A& a) { return hc(obj<Frame>(a[0]), a[1].s, a[2].s, 2.0); });
  b.fn("hclust", {"Frame", "string", "string", "float"}, "Hclust", [hc](A& a) { return hc(obj<Frame>(a[0]), a[1].s, a[2].s, a[3].f); });
  b.m("Hclust", "n", {}, "int", [](void* s, A&) { return Value::integer(static_cast<int64_t>(self<HclustFit>(s).n)); });
  b.m("Hclust", "height", {"int"}, "float", [](void* s, A& a) {
    auto& f = self<HclustFit>(s);
    return Value::real(f.height[static_cast<size_t>(as_index(a[0].i, f.height.size(), "Hclust.height"))]);
  });
  b.m("Hclust", "heights", {}, "Vec", [vec](void* s, A&) { return vec(self<HclustFit>(s).height); });
  b.m("Hclust", "merge_a", {"int"}, "int", [](void* s, A& a) {
    auto& f = self<HclustFit>(s);
    return Value::integer(f.merge_a[static_cast<size_t>(as_index(a[0].i, f.merge_a.size(), "Hclust.merge_a"))]);
  });
  b.m("Hclust", "merge_b", {"int"}, "int", [](void* s, A& a) {
    auto& f = self<HclustFit>(s);
    return Value::integer(f.merge_b[static_cast<size_t>(as_index(a[0].i, f.merge_b.size(), "Hclust.merge_b"))]);
  });
  b.m("Hclust", "order", {}, "Vec", [vec](void* s, A&) {
    Vec out;
    for (int o : self<HclustFit>(s).order) out.push_back(o);
    return vec(out);
  });
  b.m("Hclust", "cut", {"int"}, "Vec", [vec](void* s, A& a) {
    auto& f = self<HclustFit>(s);
    if (a[0].i < 1 || static_cast<size_t>(a[0].i) > f.n) throw dslvm::Error("Hclust.cut: k must be between 1 and n");
    Vec out;
    for (int g : f.cutree(static_cast<int>(a[0].i))) out.push_back(g);
    return vec(out);
  });
}

}  // namespace rcli::host
