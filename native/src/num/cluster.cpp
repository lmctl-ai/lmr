// Clustering: k-means and agglomerative hierarchical clustering.
//
// Clean-room implementation from published descriptions:
//  - Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-Means Clustering Algorithm.
//    Applied Statistics 28(1), 100-108. (Optimal-transfer and quick-transfer stages with live sets.)
//  - Lloyd, S. P. (1957, publ. 1982) / Forgy, E. W. (1965): alternate assignment and centroid steps.
//  - MacQueen, J. (1967): online centroid updates after each reassignment.
//  - Lance, G. N. and Williams, W. T. (1967). A general theory of classificatory sorting strategies.
//    Computer Journal 9, 373-380. (Dissimilarity update formulae.)
//  - Murtagh, F. (1985). Multidimensional Clustering Algorithms (nearest-neighbour list scheme for
//    stored-dissimilarity agglomeration).
// Behavioural details (iteration counts, ifault codes, merge-row ordering, tie handling) were determined by
// black-box comparison with R 4.3.3 outputs and from R's help pages (?kmeans, ?hclust, ?dist, ?cutree).
#include "rcli/cluster.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

#include "rcli/rng.hpp"

namespace rcli {

namespace {

void check_finite(const Matrix& x, const char* fn) {
  if (x.rows == 0 || x.cols == 0) throw Error(std::string(fn) + ": empty matrix");
  if (x.a.size() != x.rows * x.cols) throw Error(std::string(fn) + ": matrix storage size does not match dimensions");
  for (double v : x.a)
    if (!std::isfinite(v)) throw Error(std::string(fn) + ": NA/NaN/Inf in data");
}

inline double sqdist(const Matrix& x, size_t i, const Matrix& c, size_t l) {
  double s = 0.0;
  const double* xi = &x.a[i * x.cols];
  const double* cl = &c.a[l * c.cols];
  for (size_t j = 0; j < x.cols; ++j) {
    double d = xi[j] - cl[j];
    s += d * d;
  }
  return s;
}

// Row indices of the first occurrence of each distinct row, in row order.
std::vector<size_t> distinct_rows(const Matrix& x) {
  std::map<std::vector<double>, size_t> seen;
  std::vector<size_t> out;
  for (size_t i = 0; i < x.rows; ++i) {
    std::vector<double> r(x.a.begin() + i * x.cols, x.a.begin() + (i + 1) * x.cols);
    for (double& v : r)
      if (v == 0.0) v = 0.0;  // -0 and +0 are the same row
    if (seen.emplace(std::move(r), i).second) out.push_back(i);
  }
  return out;
}

// Final summary: withinss, size, totss etc. from x, cluster (0-based) and centres.
void finish(const Matrix& x, const std::vector<int>& cl0, KmeansFit& f) {
  const size_t m = x.rows, p = x.cols, k = f.centers.rows;
  f.withinss.assign(k, 0.0);
  f.size.assign(k, 0);
  f.cluster.resize(m);
  for (size_t i = 0; i < m; ++i) {
    size_t l = static_cast<size_t>(cl0[i]);
    f.cluster[i] = cl0[i] + 1;
    f.size[l] += 1;
    f.withinss[l] += sqdist(x, i, f.centers, l);
  }
  // Total sum of squares about the column means.
  // Mean with one refinement pass and a compensated (Neumaier) sum, so that data with a large common offset
  // still give an accurate total (R accumulates these sums in extended precision).
  double tss = 0.0, comp = 0.0;
  auto acc = [&](double v) {
    double t = tss + v;
    comp += std::fabs(tss) >= std::fabs(v) ? (tss - t) + v : (v - t) + tss;
    tss = t;
  };
  for (size_t j = 0; j < p; ++j) {
    double mu = 0.0;
    for (size_t i = 0; i < m; ++i) mu += x(i, j);
    mu /= static_cast<double>(m);
    double corr = 0.0;
    for (size_t i = 0; i < m; ++i) corr += x(i, j) - mu;
    mu += corr / static_cast<double>(m);
    for (size_t i = 0; i < m; ++i) {
      double d = x(i, j) - mu;
      acc(d * d);
    }
  }
  f.totss = tss + comp;
  f.tot_withinss = 0.0;
  for (double w : f.withinss) f.tot_withinss += w;
  f.betweenss = f.totss - f.tot_withinss;
}

// ---------------------------------------------------------------------------------------------------------------
// Hartigan & Wong (1979), AS 136.
// State follows the paper: ic1/ic2 (closest and second-closest cluster), nc (sizes), an1 = n/(n-1),
// an2 = n/(n+1), itran (cluster updated in the last quick-transfer stage), ncp (step of last update),
// live (live-set boundary), d (an1-weighted distance of a point to its own cluster).
// Indices i are 1-based inside the live-set/step comparisons, as in the published algorithm.
struct HW {
  const Matrix& a;
  Matrix& c;
  size_t m, k, p;
  std::vector<int> ic1, ic2, nc, itran, live;
  std::vector<long long> ncp;
  Vec an1, an2, d;
  long long indx = 0;
  static constexpr double kBig = 1.0e30;

  HW(const Matrix& x, Matrix& cen) : a(x), c(cen), m(x.rows), k(cen.rows), p(x.cols) {}

  void transfer(size_t i, int l1, int l2) {
    double al1 = nc[l1], alw = al1 - 1.0, al2 = nc[l2], alt = al2 + 1.0;
    for (size_t j = 0; j < p; ++j) {
      c(l1, j) = (c(l1, j) * al1 - a(i, j)) / alw;
      c(l2, j) = (c(l2, j) * al2 + a(i, j)) / alt;
    }
    nc[l1] -= 1;
    nc[l2] += 1;
    an2[l1] = alw / al1;
    an1[l1] = kBig;
    if (alw > 1.0) an1[l1] = alw / (alw - 1.0);
    an1[l2] = alt / al2;
    an2[l2] = alt / (alt + 1.0);
    ic1[i] = l2;
    ic2[i] = l1;
  }

  // Optimal-transfer stage.
  void optra() {
    const long long M = static_cast<long long>(m);
    for (size_t l = 0; l < k; ++l)
      if (itran[l] == 1) live[l] = static_cast<int>(M + 1);
    for (size_t i = 0; i < m; ++i) {
      const long long ii = static_cast<long long>(i) + 1;
      ++indx;
      int l1 = ic1[i];
      if (nc[l1] != 1) {
        // Re-compute D(I) if cluster l1 was updated in this stage.
        if (ncp[l1] != 0) d[i] = sqdist(a, i, c, l1) * an1[l1];
        int l2 = ic2[i];
        const int ll = l2;
        double r2 = sqdist(a, i, c, l2) * an2[l2];
        for (size_t l = 0; l < k; ++l) {
          const int li = static_cast<int>(l);
          // If l1 is not in the live set, only clusters in the live set are candidates.
          if ((ii >= live[l1] && ii >= live[l]) || li == l1 || li == ll) continue;
          const double rr = r2 / an2[l];
          double dc = 0.0;
          bool beaten = false;
          for (size_t j = 0; j < p; ++j) {
            double dd = a(i, j) - c(l, j);
            dc += dd * dd;
            if (dc >= rr) { beaten = true; break; }
          }
          if (beaten) continue;
          r2 = dc * an2[l];
          l2 = li;
        }
        if (r2 >= d[i]) {
          ic2[i] = l2;  // no transfer: l2 is the new second-closest cluster
        } else {
          indx = 0;
          live[l1] = static_cast<int>(M + ii);
          live[l2] = static_cast<int>(M + ii);
          ncp[l1] = ii;
          ncp[l2] = ii;
          transfer(i, l1, l2);
        }
      }
      if (indx == M) return;
    }
    for (size_t l = 0; l < k; ++l) {
      itran[l] = 0;
      live[l] -= static_cast<int>(M);
    }
  }

  // Quick-transfer stage. Returns false if the step limit was exceeded.
  bool qtran(long long max_steps) {
    const long long M = static_cast<long long>(m);
    long long icoun = 0, istep = 0;
    for (;;) {
      for (size_t i = 0; i < m; ++i) {
        ++icoun;
        ++istep;
        if (istep >= max_steps) return false;
        const int l1 = ic1[i], l2 = ic2[i];
        if (nc[l1] != 1) {
          // Distance to own cluster needs refreshing if l1 was updated within the last m steps.
          if (istep <= ncp[l1]) d[i] = sqdist(a, i, c, l1) * an1[l1];
          if (istep < ncp[l1] || istep < ncp[l2]) {
            const double r2 = d[i] / an2[l2];
            double dd = 0.0;
            bool beaten = false;
            for (size_t j = 0; j < p; ++j) {
              double de = a(i, j) - c(l2, j);
              dd += de * de;
              if (dd >= r2) { beaten = true; break; }
            }
            if (!beaten) {
              icoun = 0;
              indx = 0;
              itran[l1] = 1;
              itran[l2] = 1;
              ncp[l1] = istep + M;
              ncp[l2] = istep + M;
              transfer(i, l1, l2);
            }
          }
        }
        if (icoun == M) return true;
      }
    }
  }

  void run(int iter_max, KmeansFit& f) {
    ic1.assign(m, 0);
    ic2.assign(m, 0);
    nc.assign(k, 0);
    itran.assign(k, 1);
    live.assign(k, 0);
    ncp.assign(k, -1);
    an1.assign(k, 0.0);
    an2.assign(k, 0.0);
    d.assign(m, 0.0);
    // Step 1: closest and second-closest centres for each point.
    for (size_t i = 0; i < m; ++i) {
      int c1 = 0, c2 = 1;
      double dt1 = sqdist(a, i, c, 0), dt2 = sqdist(a, i, c, 1);
      if (dt1 > dt2) { std::swap(c1, c2); std::swap(dt1, dt2); }
      for (size_t l = 2; l < k; ++l) {
        double db = sqdist(a, i, c, l);
        if (db >= dt2) continue;
        if (db >= dt1) { dt2 = db; c2 = static_cast<int>(l); }
        else { dt2 = dt1; c2 = c1; dt1 = db; c1 = static_cast<int>(l); }
      }
      ic1[i] = c1;
      ic2[i] = c2;
    }
    // Step 2: centres as cluster means.
    std::fill(c.a.begin(), c.a.end(), 0.0);
    for (size_t i = 0; i < m; ++i) {
      nc[ic1[i]] += 1;
      for (size_t j = 0; j < p; ++j) c(ic1[i], j) += a(i, j);
    }
    for (size_t l = 0; l < k; ++l)
      if (nc[l] == 0) throw Error("kmeans: empty cluster: try a better set of initial centers");
    for (size_t l = 0; l < k; ++l) {
      double aa = nc[l];
      for (size_t j = 0; j < p; ++j) c(l, j) /= aa;
      an2[l] = aa / (aa + 1.0);
      an1[l] = kBig;
      if (aa > 1.0) an1[l] = aa / (aa - 1.0);
    }
    indx = 0;
    const long long max_qtr = 50LL * static_cast<long long>(m);
    f.ifault = 2;
    f.iter = iter_max + 1;
    for (int ij = 1; ij <= iter_max; ++ij) {
      optra();
      if (indx == static_cast<long long>(m)) { f.ifault = 0; f.iter = ij; break; }
      if (!qtran(max_qtr)) { f.ifault = 4; f.iter = ij; break; }
      if (k == 2) { f.ifault = 0; f.iter = ij; break; }
      std::fill(ncp.begin(), ncp.end(), 0);
    }
    // Final centres and within-cluster sums of squares, recomputed from the memberships.
    std::fill(c.a.begin(), c.a.end(), 0.0);
    std::vector<int> cnt(k, 0);
    for (size_t i = 0; i < m; ++i) {
      cnt[ic1[i]] += 1;
      for (size_t j = 0; j < p; ++j) c(ic1[i], j) += a(i, j);
    }
    for (size_t l = 0; l < k; ++l)
      for (size_t j = 0; j < p; ++j) c(l, j) /= static_cast<double>(cnt[l]);
    finish(a, ic1, f);
  }
};

// Nearest centre (first minimum, strict comparison); NaN centres are never chosen.
inline int nearest(const Matrix& x, size_t i, const Matrix& c) {
  int best = 0;
  double bd = kInf;
  for (size_t l = 0; l < c.rows; ++l) {
    double dd = sqdist(x, i, c, l);
    if (dd < bd) { bd = dd; best = static_cast<int>(l); }
  }
  return best;
}

void lloyd(const Matrix& x, Matrix& c, int iter_max, KmeansFit& f) {
  const size_t m = x.rows, k = c.rows, p = x.cols;
  std::vector<int> cl(m, -1), cnt(k);
  int it = 0;
  for (it = 0; it < iter_max; ++it) {
    bool updated = false;
    for (size_t i = 0; i < m; ++i) {
      int b = nearest(x, i, c);
      if (cl[i] != b) { updated = true; cl[i] = b; }
    }
    if (!updated) break;
    std::fill(c.a.begin(), c.a.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), 0);
    for (size_t i = 0; i < m; ++i) {
      cnt[cl[i]] += 1;
      for (size_t j = 0; j < p; ++j) c(cl[i], j) += x(i, j);
    }
    for (size_t l = 0; l < k; ++l)
      for (size_t j = 0; j < p; ++j) c(l, j) /= static_cast<double>(cnt[l]);  // empty cluster -> NaN centre
  }
  f.iter = it + 1;
  f.ifault = f.iter > iter_max ? 2 : 0;
  finish(x, cl, f);
}

void macqueen(const Matrix& x, Matrix& c, int iter_max, KmeansFit& f) {
  const size_t m = x.rows, k = c.rows, p = x.cols;
  std::vector<int> cl(m), cnt(k, 0);
  for (size_t i = 0; i < m; ++i) cl[i] = nearest(x, i, c);
  std::fill(c.a.begin(), c.a.end(), 0.0);
  for (size_t i = 0; i < m; ++i) {
    cnt[cl[i]] += 1;
    for (size_t j = 0; j < p; ++j) c(cl[i], j) += x(i, j);
  }
  for (size_t l = 0; l < k; ++l)
    for (size_t j = 0; j < p; ++j) c(l, j) /= static_cast<double>(cnt[l]);
  int it = 0;
  for (it = 0; it < iter_max; ++it) {
    bool updated = false;
    for (size_t i = 0; i < m; ++i) {
      int b = nearest(x, i, c);
      int o = cl[i];
      if (b == o) continue;
      updated = true;
      cl[i] = b;
      cnt[o] -= 1;
      cnt[b] += 1;
      for (size_t j = 0; j < p; ++j) {
        if (cnt[o] > 0) c(o, j) += (c(o, j) - x(i, j)) / cnt[o];
        else c(o, j) = kNaN;
        c(b, j) += (x(i, j) - c(b, j)) / cnt[b];
      }
    }
    if (!updated) break;
  }
  f.iter = it + 1;
  f.ifault = f.iter > iter_max ? 2 : 0;
  finish(x, cl, f);
}

std::string kmeans_algorithm(const std::string& a) {
  if (a == "Hartigan-Wong" || a == "Lloyd" || a == "MacQueen") return a;
  if (a == "Forgy") return "Lloyd";
  throw Error("kmeans: unknown algorithm '" + a + "'");
}

KmeansFit kmeans_one(const Matrix& x, const Matrix& centers, int iter_max, const std::string& alg) {
  KmeansFit f;
  f.centers = centers;
  const size_t k = centers.rows;
  if (k == 1) {
    // k = 1: the centre is the column mean.
    for (size_t j = 0; j < x.cols; ++j) {
      double s = 0.0;
      for (size_t i = 0; i < x.rows; ++i) s += x(i, j);
      f.centers(0, j) = s / static_cast<double>(x.rows);
    }
    f.iter = 1;
    f.ifault = 0;
    finish(x, std::vector<int>(x.rows, 0), f);
    return f;
  }
  if (alg == "Hartigan-Wong") {
    if (k >= x.rows) throw Error("kmeans: number of cluster centres must lie between 1 and nrow(x)");
    HW hw(x, f.centers);
    hw.run(iter_max, f);
  } else if (alg == "Lloyd") {
    lloyd(x, f.centers, iter_max, f);
  } else {
    macqueen(x, f.centers, iter_max, f);
  }
  return f;
}

void kmeans_check_opts(const KmeansOptions& opt) {
  if (opt.iter_max < 1) throw Error("kmeans: 'iter_max' must be positive");
  if (opt.nstart < 1) throw Error("kmeans: 'nstart' must be positive");
}

}  // namespace

KmeansFit kmeans(const Matrix& x, const Matrix& centers, const KmeansOptions& opt) {
  check_finite(x, "kmeans");
  kmeans_check_opts(opt);
  const std::string alg = kmeans_algorithm(opt.algorithm);
  if (centers.rows == 0) throw Error("kmeans: no initial centres");
  if (centers.cols != x.cols) throw Error("kmeans: must have same number of columns in 'x' and 'centers'");
  if (centers.a.size() != centers.rows * centers.cols) throw Error("kmeans: centres storage size does not match");
  for (double v : centers.a)
    if (!std::isfinite(v)) throw Error("kmeans: NA/NaN/Inf in initial centres");
  if (centers.rows > x.rows) throw Error("kmeans: more cluster centers than data points");
  if (distinct_rows(centers).size() != centers.rows) throw Error("kmeans: initial centers are not distinct");
  return kmeans_one(x, centers, opt.iter_max, alg);
}

KmeansFit kmeans(const Matrix& x, int k, Rng& rng, const KmeansOptions& opt) {
  check_finite(x, "kmeans");
  kmeans_check_opts(opt);
  const std::string alg = kmeans_algorithm(opt.algorithm);
  if (k < 1) throw Error("kmeans: number of clusters k must be at least 1");
  const std::vector<size_t> uniq = distinct_rows(x);
  if (static_cast<size_t>(k) > uniq.size()) throw Error("kmeans: more cluster centers than distinct data points");
  const size_t kk = static_cast<size_t>(k);
  auto draw = [&]() {
    std::vector<size_t> idx = rng.sample(uniq.size(), kk, false);
    Matrix c(kk, x.cols);
    for (size_t l = 0; l < kk; ++l)
      for (size_t j = 0; j < x.cols; ++j) c(l, j) = x(uniq[idx[l]], j);
    return c;
  };
  KmeansFit best = kmeans_one(x, draw(), opt.iter_max, alg);
  for (int s = 1; s < opt.nstart; ++s) {
    KmeansFit f = kmeans_one(x, draw(), opt.iter_max, alg);
    if (f.tot_withinss < best.tot_withinss) best = std::move(f);
  }
  return best;
}

// ---------------------------------------------------------------------------------------------------------------
// dist(): distances between rows. Missing (NaN) terms, and terms whose contribution is NaN (e.g. Inf - Inf), are
// dropped; for euclidean, manhattan, canberra and minkowski the sum is scaled up by (columns / columns used);
// if no term is usable the distance is NaN (as ?dist documents).

Vec dist(const Matrix& x, const std::string& method, double p) {
  enum { EUC, MAXI, MAN, CAN, BIN, MINK } meth;
  if (method == "euclidean") meth = EUC;
  else if (method == "maximum") meth = MAXI;
  else if (method == "manhattan") meth = MAN;
  else if (method == "canberra") meth = CAN;
  else if (method == "binary") meth = BIN;
  else if (method == "minkowski") meth = MINK;
  else throw Error("dist: unknown method '" + method + "'");
  if (meth == MINK && !(p > 0.0 && std::isfinite(p))) throw Error("dist: 'p' must be a positive finite number");
  if (x.a.size() != x.rows * x.cols) throw Error("dist: matrix storage size does not match dimensions");
  const size_t n = x.rows, nc = x.cols;
  Vec out(n < 2 ? 0 : n * (n - 1) / 2);
  size_t idx = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      double s = 0.0;
      size_t cnt = 0;
      double res;
      if (meth == BIN) {
        size_t total = 0, diff = 0;
        for (size_t c = 0; c < nc; ++c) {
          double u = x(i, c), v = x(j, c);
          if (!std::isfinite(u) || !std::isfinite(v)) continue;  // non-finite pairs excluded
          ++cnt;
          bool bu = u != 0.0, bv = v != 0.0;
          if (bu || bv) {
            ++total;
            if (!(bu && bv)) ++diff;
          }
        }
        if (cnt == 0) res = kNaN;
        else if (total == 0) res = 0.0;
        else res = static_cast<double>(diff) / static_cast<double>(total);
      } else {
        for (size_t c = 0; c < nc; ++c) {
          double u = x(i, c), v = x(j, c);
          if (std::isnan(u) || std::isnan(v)) continue;
          double dev;
          switch (meth) {
            case EUC: dev = (u - v) * (u - v); break;
            case MAXI: dev = std::fabs(u - v); break;
            case MAN: dev = std::fabs(u - v); break;
            case CAN: {
              double num = std::fabs(u - v), den = std::fabs(u) + std::fabs(v);
              if (num == 0.0 && den == 0.0) continue;  // 0/0 term omitted, treated as missing
              // An infinite numerator over an infinite denominator counts as 1 (R treats Inf like a huge
              // finite value here); Inf - Inf gives NaN and the term is dropped below.
              dev = (std::isinf(num) && std::isinf(den)) ? 1.0 : num / den;
              break;
            }
            default: dev = std::pow(std::fabs(u - v), p); break;
          }
          if (std::isnan(dev)) continue;
          ++cnt;
          if (meth == MAXI) { if (dev > s || cnt == 1) s = dev; }
          else s += dev;
        }
        if (cnt == 0) res = kNaN;
        else if (meth == MAXI) res = s;
        else {
          if (cnt != nc) s /= (static_cast<double>(cnt) / static_cast<double>(nc));
          if (meth == EUC) res = std::sqrt(s);
          else if (meth == MINK) res = std::pow(s, 1.0 / p);
          else res = s;
        }
      }
      out[idx++] = res;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// hclust(): stored-dissimilarity agglomeration with nearest-neighbour lists (Murtagh 1985) and the
// Lance-Williams updates as listed in ?hclust. Each cluster lives in the slot of its smallest-index member
// slot i2 = min of the pair; slot j2 is retired. The nearest-neighbour list of slot i scans only slots j > i.

HclustFit hclust(const Vec& d, size_t n, const std::string& method) {
  enum { WARD, WARD2, SINGLE, COMPLETE, AVERAGE, MCQUITTY, MEDIAN, CENTROID } meth;
  if (method == "ward.D") meth = WARD;
  else if (method == "ward.D2") meth = WARD2;
  else if (method == "single") meth = SINGLE;
  else if (method == "complete") meth = COMPLETE;
  else if (method == "average") meth = AVERAGE;
  else if (method == "mcquitty") meth = MCQUITTY;
  else if (method == "median") meth = MEDIAN;
  else if (method == "centroid") meth = CENTROID;
  else throw Error("hclust: unknown method '" + method + "'");
  if (n < 2) throw Error("hclust: must have n >= 2 objects to cluster");
  if (d.size() != n * (n - 1) / 2) throw Error("hclust: dissimilarities of improper length");
  for (double v : d)
    if (std::isnan(v)) throw Error("hclust: NA/NaN in dissimilarities");

  // Full symmetric working matrix (n x n), diss(i, j) for i != j.
  std::vector<double> D(n * n, 0.0);
  {
    size_t idx = 0;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) {
        double v = d[idx++];
        if (meth == WARD2) v *= v;
        D[i * n + j] = D[j * n + i] = v;
      }
  }
  auto at = [&](size_t i, size_t j) -> double& { return D[i * n + j]; };

  std::vector<char> active(n, 1);
  std::vector<double> membr(n, 1.0);
  std::vector<size_t> nn(n, 0);
  std::vector<double> disnn(n, kInf);
  std::vector<int> laststep(n, 0);  // 0: singleton; else 1-based step that formed the cluster in this slot

  auto recompute = [&](size_t i) {
    double best = kInf;
    size_t bj = i;
    for (size_t j = i + 1; j < n; ++j) {
      if (!active[j]) continue;
      if (at(i, j) < best) { best = at(i, j); bj = j; }
    }
    nn[i] = bj;
    disnn[i] = best;
  };
  for (size_t i = 0; i + 1 < n; ++i) recompute(i);

  HclustFit f;
  f.method = method;
  f.n = n;
  f.merge_a.resize(n - 1);
  f.merge_b.resize(n - 1);
  f.height.resize(n - 1);

  for (size_t step = 1; step < n; ++step) {
    // Least dissimilarity among the nearest-neighbour lists (first minimum in slot order).
    double dmin = kInf;
    size_t im = n, jm = n;
    for (size_t i = 0; i + 1 < n; ++i) {
      if (!active[i] || nn[i] == i) continue;
      if (disnn[i] < dmin) { dmin = disnn[i]; im = i; jm = nn[i]; }
    }
    if (im == n) {
      // Only possible with infinite dissimilarities: take the first remaining active pair.
      for (size_t i = 0; i < n && im == n; ++i)
        if (active[i])
          for (size_t j = i + 1; j < n; ++j)
            if (active[j]) { im = i; jm = j; dmin = at(i, j); break; }
    }
    const size_t i2 = std::min(im, jm), j2 = std::max(im, jm);

    // Merge row in R's convention.
    int ea = laststep[i2] ? laststep[i2] : -static_cast<int>(i2 + 1);
    int eb = laststep[j2] ? laststep[j2] : -static_cast<int>(j2 + 1);
    if ((ea > 0 && eb < 0) || (ea > 0 && eb > 0 && ea > eb)) std::swap(ea, eb);
    f.merge_a[step - 1] = ea;
    f.merge_b[step - 1] = eb;
    f.height[step - 1] = meth == WARD2 ? std::sqrt(dmin) : dmin;

    // Lance-Williams update of dissimilarities from the new cluster (slot i2) to every other active slot.
    const double ni = membr[i2], nj = membr[j2], dij = at(i2, j2);
    for (size_t k = 0; k < n; ++k) {
      if (!active[k] || k == i2 || k == j2) continue;
      const double nk = membr[k], dki = at(k, i2), dkj = at(k, j2);
      double v;
      switch (meth) {
        case WARD:
        case WARD2: v = ((ni + nk) * dki + (nj + nk) * dkj - nk * dij) / (ni + nj + nk); break;
        case SINGLE: v = std::min(dki, dkj); break;
        case COMPLETE: v = std::max(dki, dkj); break;
        case AVERAGE: v = (ni * dki + nj * dkj) / (ni + nj); break;
        case MCQUITTY: v = 0.5 * dki + 0.5 * dkj; break;
        case MEDIAN: v = 0.5 * dki + 0.5 * dkj - 0.25 * dij; break;
        default: v = (ni * dki + nj * dkj - ni * nj * dij / (ni + nj)) / (ni + nj); break;
      }
      at(k, i2) = at(i2, k) = v;
    }
    membr[i2] = ni + nj;
    active[j2] = 0;
    laststep[i2] = static_cast<int>(step);

    // Update nearest-neighbour lists.
    for (size_t i = 0; i + 1 < n; ++i) {
      if (!active[i]) continue;
      if (i == i2 || nn[i] == i2 || nn[i] == j2) recompute(i);
      else if (i < i2 && at(i, i2) < disnn[i]) { disnn[i] = at(i, i2); nn[i] = i2; }
    }
  }

  // Leaf order: left-to-right traversal of the merge tree from the last merge.
  f.order.clear();
  f.order.reserve(n);
  // Stack entries: positive = merge step, negative = leaf (observation, 1-based).
  std::vector<int> stack{static_cast<int>(n - 1)};
  while (!stack.empty()) {
    const int s = stack.back();
    stack.pop_back();
    if (s < 0) {
      f.order.push_back(-s);
    } else {
      stack.push_back(f.merge_b[s - 1]);  // right part visited after the left part
      stack.push_back(f.merge_a[s - 1]);
    }
  }
  return f;
}

std::vector<int> HclustFit::cutree(int k) const {
  if (n < 1) throw Error("cutree: empty tree");
  if (k < 1 || static_cast<size_t>(k) > n) throw Error("cutree: k must be between 1 and the number of observations");
  // Union the first n - k merges.
  std::vector<size_t> parent(n);
  std::iota(parent.begin(), parent.end(), size_t{0});
  auto find = [&](size_t v) {
    while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
    return v;
  };
  std::vector<size_t> rep(n > 0 ? n - 1 : 0);  // an observation in the cluster formed at each step
  const size_t steps = n - static_cast<size_t>(k);
  for (size_t s = 0; s < steps; ++s) {
    auto obs = [&](int e) { return e < 0 ? static_cast<size_t>(-e - 1) : rep[static_cast<size_t>(e - 1)]; };
    size_t u = find(obs(merge_a[s])), v = find(obs(merge_b[s]));
    parent[v] = u;
    rep[s] = u;
  }
  std::vector<int> group(n, 0), label(n, 0);
  int next = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t r = find(i);
    if (label[r] == 0) label[r] = ++next;
    group[i] = label[r];
  }
  return group;
}

}  // namespace rcli
