// Tests for the noncentral t, chi-square and F distributions and the studentized range (ncdist.cpp).
//
// Block 1: expectations from R 4.3.3 used as a black-box oracle (the R expression is in each line's comment).
// Block 2: cases where R's own value is less accurate (R documents that its noncentral code "will not be highly
// accurate, especially in the tails", that pt/pchisq compute the lower tail only so the upper tail cancels, and that
// ptukey uses a fixed 16-point Legendre rule). There the expectation comes from an independent reference computed in
// R by black-box numerical integration or summation of central quantities (helpers below, our own code), and R's
// value is quoted for comparison:
//   lint <- function(lf) { # log of the integral of exp(lf(s)) over s > 0, in u = log(s), split at the grid maximum
//     lu <- function(u) sapply(u, function(uu) lf(exp(uu)) + uu)
//     g <- seq(-150, 25, length.out = 35001); v <- lu(g); i <- which.max(v); M <- v[i]
//     ok <- which(v > M - 60); lo <- g[max(1, min(ok) - 1)]; hi <- g[min(length(g), max(ok) + 1)]
//     f <- function(u) exp(lu(u) - M)
//     M + log(integrate(f, lo, g[i], rel.tol = 1e-12, subdivisions = 2000, stop.on.error = FALSE)$value +
//             integrate(f, g[i], hi, rel.tol = 1e-12, subdivisions = 2000, stop.on.error = FALSE)$value) }
//   ptref <- function(t, n, d, lower = TRUE)   # log P(T <= t) or log P(T > t): T = (Z + d) / sqrt(S / n)
//     lint(function(s) pnorm(t * sqrt(s / n) - d, lower.tail = lower, log.p = TRUE) + dchisq(s, n, log = TRUE))
//   dtref <- function(t, n, d) lint(function(s) dnorm(t * sqrt(s / n) - d, log = TRUE) + 0.5 * log(s / n) + dchisq(s, n, log = TRUE))
//   pcref <- function(x, k, l, lower = TRUE) { j <- 0:3000; sum(dpois(j, l / 2) * pchisq(x, k + 2 * j, lower.tail = lower)) }
//   dcref <- function(x, k, l) { j <- 0:3000; sum(dpois(j, l / 2) * dchisq(x, k + 2 * j)) }
//   pfref <- function(x, m, n, l, lower = TRUE) { j <- 0:3000; y <- m * x / (m * x + n)
//     sum(dpois(j, l / 2) * pbeta(y, m / 2 + j, n / 2, lower.tail = lower)) }
//   dfref <- function(x, m, n, l) { j <- 0:3000; y <- m * x / (m * x + n); sum(dpois(j, l / 2) * dbeta(y, m / 2 + j, n / 2)) * m * n / (m * x + n)^2 }
//   Wref <- function(w, k) integrate(function(z) k * dnorm(z) * (pnorm(z) - pnorm(z - w))^(k - 1), -Inf, Inf,
//                                    rel.tol = 1e-13, subdivisions = 1000)$value
//   ptkref <- function(q, k, df, r = 1) { if (is.infinite(df)) return(Wref(q, k)^r)
//     f <- function(x) sapply(x, function(xx) dchisq(xx, df) * Wref(q * sqrt(xx / df), k)^r)
//     integrate(f, 0, df, rel.tol = 1e-12, subdivisions = 1000)$value + integrate(f, df, Inf, rel.tol = 1e-12, subdivisions = 1000)$value }
// Block 3: round trips, identities, moments, edge cases and the API.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "check.hpp"
#include "rcli/dist.hpp"

using rcli::Dist;
using rcli::DistKind;
using rcli::kInf;
using rcli::kNaN;

namespace {

void r_oracle_cases() {
  check::rel(Dist::t(1, 0.5).p(-2, true, false), 0.07457768951886379, 1e-09, "pt(-2,1,ncp=0.5)");  // R: pt(-2,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).p(-2, false, false), 0.9254223104811362, 1e-09, "pt(-2,1,ncp=0.5,lower.tail=FALSE)");  // R: pt(-2,1,ncp=0.5,lower.tail=FALSE)
  check::rel(Dist::t(1, 0.5).p(1, true, false), 0.5927477364399605, 1e-09, "pt(1,1,ncp=0.5)");  // R: pt(1,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).p(1, false, false), 0.4072522635600395, 1e-09, "pt(1,1,ncp=0.5,lower.tail=FALSE)");  // R: pt(1,1,ncp=0.5,lower.tail=FALSE)
  check::rel(Dist::t(1, 0.5).p(5, true, false), 0.8906260124443279, 1e-09, "pt(5,1,ncp=0.5)");  // R: pt(5,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).p(5, false, false), 0.10937398755567207, 1e-09, "pt(5,1,ncp=0.5,lower.tail=FALSE)");  // R: pt(5,1,ncp=0.5,lower.tail=FALSE)
  check::rel(Dist::t(1, 0.5).p(-2, true, true), -2.5959138846139536, 1e-09, "pt(-2,1,ncp=0.5,log.p=TRUE)");  // R: pt(-2,1,ncp=0.5,log.p=TRUE)
  check::rel(Dist::t(1, 0.5).p(5, false, true), -2.212982190981057, 1e-09, "pt(5,1,ncp=0.5,lower.tail=FALSE,log.p=TRUE)");  // R: pt(5,1,ncp=0.5,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(1, 0.5).d(-2, false), 0.033396244149602516, 1e-09, "dt(-2,1,ncp=0.5)");  // R: dt(-2,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).d(1, false), 0.22501149476431062, 1e-09, "dt(1,1,ncp=0.5)");  // R: dt(1,1,ncp=0.5)
  check::rel(Dist::t(3, 1).p(-1, true, false), 0.03597657076124339, 1e-09, "pt(-1,3,ncp=1)");  // R: pt(-1,3,ncp=1)
  check::rel(Dist::t(3, 1).p(-1, false, false), 0.9640234292387566, 1e-09, "pt(-1,3,ncp=1,lower.tail=FALSE)");  // R: pt(-1,3,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(3, 1).p(0.5, true, false), 0.29811049337095613, 1e-09, "pt(0.5,3,ncp=1)");  // R: pt(0.5,3,ncp=1)
  check::rel(Dist::t(3, 1).p(0.5, false, false), 0.7018895066290438, 1e-09, "pt(0.5,3,ncp=1,lower.tail=FALSE)");  // R: pt(0.5,3,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(3, 1).p(2, true, false), 0.7436298684076462, 1e-09, "pt(2,3,ncp=1)");  // R: pt(2,3,ncp=1)
  check::rel(Dist::t(3, 1).p(2, false, false), 0.2563701315923538, 1e-09, "pt(2,3,ncp=1,lower.tail=FALSE)");  // R: pt(2,3,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(3, 1).p(6, true, false), 0.9775253578225593, 1e-09, "pt(6,3,ncp=1)");  // R: pt(6,3,ncp=1)
  check::rel(Dist::t(3, 1).p(6, false, false), 0.022474642177440662, 1e-09, "pt(6,3,ncp=1,lower.tail=FALSE)");  // R: pt(6,3,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(3, 1).p(-1, true, true), -3.324887364583893, 1e-09, "pt(-1,3,ncp=1,log.p=TRUE)");  // R: pt(-1,3,ncp=1,log.p=TRUE)
  check::rel(Dist::t(3, 1).p(6, false, true), -3.795367619665957, 1e-09, "pt(6,3,ncp=1,lower.tail=FALSE,log.p=TRUE)");  // R: pt(6,3,ncp=1,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(3, 1).d(-1, false), 0.051784087936924135, 1e-09, "dt(-1,3,ncp=1)");  // R: dt(-1,3,ncp=1)
  check::rel(Dist::t(3, 1).d(0.5, false), 0.3258914646538911, 1e-09, "dt(0.5,3,ncp=1)");  // R: dt(0.5,3,ncp=1)
  check::rel(Dist::t(3, -2).p(-3, true, false), 0.3199728259450725, 1e-09, "pt(-3,3,ncp=-2)");  // R: pt(-3,3,ncp=-2)
  check::rel(Dist::t(3, -2).p(-3, false, false), 0.6800271740549275, 1e-09, "pt(-3,3,ncp=-2,lower.tail=FALSE)");  // R: pt(-3,3,ncp=-2,lower.tail=FALSE)
  check::rel(Dist::t(3, -2).p(1, true, false), 0.9969947524045499, 1e-09, "pt(1,3,ncp=-2)");  // R: pt(1,3,ncp=-2)
  check::rel(Dist::t(3, -2).p(1, false, false), 0.003005247595450067, 1e-09, "pt(1,3,ncp=-2,lower.tail=FALSE)");  // R: pt(1,3,ncp=-2,lower.tail=FALSE)
  check::rel(Dist::t(3, -2).p(4, true, false), 0.9998923719471835, 1e-09, "pt(4,3,ncp=-2)");  // R: pt(4,3,ncp=-2)
  check::rel(Dist::t(3, -2).p(-3, true, true), -1.139519205715829, 1e-09, "pt(-3,3,ncp=-2,log.p=TRUE)");  // R: pt(-3,3,ncp=-2,log.p=TRUE)
  check::rel(Dist::t(3, -2).d(-3, false), 0.181695020545246, 1e-09, "dt(-3,3,ncp=-2)");  // R: dt(-3,3,ncp=-2)
  check::rel(Dist::t(3, -2).d(1, false), 0.005278237808599153, 1e-09, "dt(1,3,ncp=-2)");  // R: dt(1,3,ncp=-2)
  check::rel(Dist::t(10, 2).p(-1, true, false), 0.0018164334105896796, 1e-09, "pt(-1,10,ncp=2)");  // R: pt(-1,10,ncp=2)
  check::rel(Dist::t(10, 2).p(-1, false, false), 0.9981835665894103, 1e-09, "pt(-1,10,ncp=2,lower.tail=FALSE)");  // R: pt(-1,10,ncp=2,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).p(1, true, false), 0.158515059537495, 1e-09, "pt(1,10,ncp=2)");  // R: pt(1,10,ncp=2)
  check::rel(Dist::t(10, 2).p(1, false, false), 0.841484940462505, 1e-09, "pt(1,10,ncp=2,lower.tail=FALSE)");  // R: pt(1,10,ncp=2,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).p(2.5, true, false), 0.6480249588698735, 1e-09, "pt(2.5,10,ncp=2)");  // R: pt(2.5,10,ncp=2)
  check::rel(Dist::t(10, 2).p(2.5, false, false), 0.35197504113012645, 1e-09, "pt(2.5,10,ncp=2,lower.tail=FALSE)");  // R: pt(2.5,10,ncp=2,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).p(5, true, false), 0.9761872541463537, 1e-09, "pt(5,10,ncp=2)");  // R: pt(5,10,ncp=2)
  check::rel(Dist::t(10, 2).p(5, false, false), 0.023812745853646344, 1e-09, "pt(5,10,ncp=2,lower.tail=FALSE)");  // R: pt(5,10,ncp=2,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).p(-1, true, true), -6.3108803650593455, 1e-09, "pt(-1,10,ncp=2,log.p=TRUE)");  // R: pt(-1,10,ncp=2,log.p=TRUE)
  check::rel(Dist::t(10, 2).p(5, false, true), -3.737534301585939, 1e-09, "pt(5,10,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: pt(5,10,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(10, 2).d(1, false), 0.2413718676179922, 1e-09, "dt(1,10,ncp=2)");  // R: dt(1,10,ncp=2)
  check::rel(Dist::t(10, -1.5).p(-4, true, false), 0.03378569625204486, 1e-09, "pt(-4,10,ncp=-1.5)");  // R: pt(-4,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).p(-4, false, false), 0.9662143037479551, 1e-09, "pt(-4,10,ncp=-1.5,lower.tail=FALSE)");  // R: pt(-4,10,ncp=-1.5,lower.tail=FALSE)
  check::rel(Dist::t(10, -1.5).p(-1, true, false), 0.6958906710772366, 1e-09, "pt(-1,10,ncp=-1.5)");  // R: pt(-1,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).p(-1, false, false), 0.30410932892276343, 1e-09, "pt(-1,10,ncp=-1.5,lower.tail=FALSE)");  // R: pt(-1,10,ncp=-1.5,lower.tail=FALSE)
  check::rel(Dist::t(10, -1.5).p(1, true, false), 0.9922209046455744, 1e-09, "pt(1,10,ncp=-1.5)");  // R: pt(1,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).p(1, false, false), 0.007779095354425558, 1e-09, "pt(1,10,ncp=-1.5,lower.tail=FALSE)");  // R: pt(1,10,ncp=-1.5,lower.tail=FALSE)
  check::rel(Dist::t(10, -1.5).p(-4, true, true), -3.387717753871719, 1e-09, "pt(-4,10,ncp=-1.5,log.p=TRUE)");  // R: pt(-4,10,ncp=-1.5,log.p=TRUE)
  check::rel(Dist::t(10, -1.5).p(1, false, true), -4.856315225904159, 1e-09, "pt(1,10,ncp=-1.5,lower.tail=FALSE,log.p=TRUE)");  // R: pt(1,10,ncp=-1.5,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(10, -1.5).d(-4, false), 0.0428399056476561, 1e-09, "dt(-4,10,ncp=-1.5)");  // R: dt(-4,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).d(-1, false), 0.3410749319241869, 1e-09, "dt(-1,10,ncp=-1.5)");  // R: dt(-1,10,ncp=-1.5)
  check::rel(Dist::t(50, 3).p(1, true, false), 0.02301939422974286, 1e-09, "pt(1,50,ncp=3)");  // R: pt(1,50,ncp=3)
  check::rel(Dist::t(50, 3).p(1, false, false), 0.9769806057702571, 1e-09, "pt(1,50,ncp=3,lower.tail=FALSE)");  // R: pt(1,50,ncp=3,lower.tail=FALSE)
  check::rel(Dist::t(50, 3).p(3, true, false), 0.4941233065764572, 1e-09, "pt(3,50,ncp=3)");  // R: pt(3,50,ncp=3)
  check::rel(Dist::t(50, 3).p(3, false, false), 0.5058766934235428, 1e-09, "pt(3,50,ncp=3,lower.tail=FALSE)");  // R: pt(3,50,ncp=3,lower.tail=FALSE)
  check::rel(Dist::t(50, 3).p(5, true, false), 0.9616909178206544, 1e-09, "pt(5,50,ncp=3)");  // R: pt(5,50,ncp=3)
  check::rel(Dist::t(50, 3).p(5, false, false), 0.03830908217934559, 1e-09, "pt(5,50,ncp=3,lower.tail=FALSE)");  // R: pt(5,50,ncp=3,lower.tail=FALSE)
  check::rel(Dist::t(50, 3).p(1, true, true), -3.7714181909893427, 1e-09, "pt(1,50,ncp=3,log.p=TRUE)");  // R: pt(1,50,ncp=3,log.p=TRUE)
  check::rel(Dist::t(50, 3).p(5, false, true), -3.2620682782870434, 1e-09, "pt(5,50,ncp=3,lower.tail=FALSE,log.p=TRUE)");  // R: pt(5,50,ncp=3,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(50, 3).d(1, false), 0.05506503983125509, 1e-09, "dt(1,50,ncp=3)");  // R: dt(1,50,ncp=3)
  check::rel(Dist::t(50, 3).d(3, false), 0.3802560938148181, 1e-09, "dt(3,50,ncp=3)");  // R: dt(3,50,ncp=3)
  check::rel(Dist::t(1000, 1).p(0.5, true, false), 0.30850453647433435, 1e-09, "pt(0.5,1000,ncp=1)");  // R: pt(0.5,1000,ncp=1)
  check::rel(Dist::t(1000, 1).p(0.5, false, false), 0.6914954635256656, 1e-09, "pt(0.5,1000,ncp=1,lower.tail=FALSE)");  // R: pt(0.5,1000,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(1000, 1).p(2, true, false), 0.8409820769557533, 1e-09, "pt(2,1000,ncp=1)");  // R: pt(2,1000,ncp=1)
  check::rel(Dist::t(1000, 1).p(2, false, false), 0.15901792304424667, 1e-09, "pt(2,1000,ncp=1,lower.tail=FALSE)");  // R: pt(2,1000,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(1000, 1).p(0.5, true, true), -1.176018730803695, 1e-09, "pt(0.5,1000,ncp=1,log.p=TRUE)");  // R: pt(0.5,1000,ncp=1,log.p=TRUE)
  check::rel(Dist::t(1000, 1).p(2, false, true), -1.838738359566006, 1e-09, "pt(2,1000,ncp=1,lower.tail=FALSE,log.p=TRUE)");  // R: pt(2,1000,ncp=1,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(1000, 1).d(2, false), 0.24178943015817908, 1e-09, "dt(2,1000,ncp=1)");  // R: dt(2,1000,ncp=1)
  check::rel(Dist::t(10, 37).p(30, true, false), 0.1269306380074432, 1e-09, "pt(30,10,ncp=37)");  // R: pt(30,10,ncp=37)
  check::rel(Dist::t(10, 37).p(30, false, false), 0.8730693619925568, 1e-09, "pt(30,10,ncp=37,lower.tail=FALSE)");  // R: pt(30,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::t(10, 37).p(37, true, false), 0.44112559676822116, 1e-09, "pt(37,10,ncp=37)");  // R: pt(37,10,ncp=37)
  check::rel(Dist::t(10, 37).p(37, false, false), 0.5588744032317788, 1e-09, "pt(37,10,ncp=37,lower.tail=FALSE)");  // R: pt(37,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::t(10, 37).p(50, true, false), 0.8561387907735566, 1e-09, "pt(50,10,ncp=37)");  // R: pt(50,10,ncp=37)
  check::rel(Dist::t(10, 37).p(50, false, false), 0.1438612092264434, 1e-09, "pt(50,10,ncp=37,lower.tail=FALSE)");  // R: pt(50,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::t(10, 37).p(30, true, true), -2.064114499143544, 1e-09, "pt(30,10,ncp=37,log.p=TRUE)");  // R: pt(30,10,ncp=37,log.p=TRUE)
  check::rel(Dist::t(10, 37).p(50, false, true), -1.938906269000307, 1e-09, "pt(50,10,ncp=37,lower.tail=FALSE,log.p=TRUE)");  // R: pt(50,10,ncp=37,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(10, 37).d(30, false), 0.03520098849219861, 1e-09, "dt(30,10,ncp=37)");  // R: dt(30,10,ncp=37)
  check::rel(Dist::t(10, 37).d(37, false), 0.04708077081032277, 1e-09, "dt(37,10,ncp=37)");  // R: dt(37,10,ncp=37)
  check::rel(Dist::t(0.5, 4).p(2, true, false), 0.08440094200967573, 1e-09, "pt(2,0.5,ncp=4)");  // R: pt(2,0.5,ncp=4)
  check::rel(Dist::t(0.5, 4).p(2, false, false), 0.9155990579903243, 1e-09, "pt(2,0.5,ncp=4,lower.tail=FALSE)");  // R: pt(2,0.5,ncp=4,lower.tail=FALSE)
  check::rel(Dist::t(0.5, 4).p(30, true, false), 0.7178096722204167, 1e-09, "pt(30,0.5,ncp=4)");  // R: pt(30,0.5,ncp=4)
  check::rel(Dist::t(0.5, 4).p(30, false, false), 0.2821903277795833, 1e-09, "pt(30,0.5,ncp=4,lower.tail=FALSE)");  // R: pt(30,0.5,ncp=4,lower.tail=FALSE)
  check::rel(Dist::t(0.5, 4).p(2, true, true), -2.47217671619043, 1e-09, "pt(2,0.5,ncp=4,log.p=TRUE)");  // R: pt(2,0.5,ncp=4,log.p=TRUE)
  check::rel(Dist::t(0.5, 4).p(30, false, true), -1.2651735144967153, 1e-09, "pt(30,0.5,ncp=4,lower.tail=FALSE,log.p=TRUE)");  // R: pt(30,0.5,ncp=4,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(0.5, 4).d(2, false), 0.09983033554732408, 1e-09, "dt(2,0.5,ncp=4)");  // R: dt(2,0.5,ncp=4)
  check::rel(Dist::t(0.5, 4).d(30, false), 0.004684384584377413, 1e-09, "dt(30,0.5,ncp=4)");  // R: dt(30,0.5,ncp=4)
  check::rel(Dist::t(25, 10).p(8, true, false), 0.08501918644033023, 1e-09, "pt(8,25,ncp=10)");  // R: pt(8,25,ncp=10)
  check::rel(Dist::t(25, 10).p(8, false, false), 0.9149808135596698, 1e-09, "pt(8,25,ncp=10,lower.tail=FALSE)");  // R: pt(8,25,ncp=10,lower.tail=FALSE)
  check::rel(Dist::t(25, 10).p(12, true, false), 0.8305657760648368, 1e-09, "pt(12,25,ncp=10)");  // R: pt(12,25,ncp=10)
  check::rel(Dist::t(25, 10).p(12, false, false), 0.16943422393516316, 1e-09, "pt(12,25,ncp=10,lower.tail=FALSE)");  // R: pt(12,25,ncp=10,lower.tail=FALSE)
  check::rel(Dist::t(25, 10).p(8, true, true), -2.4648783251359707, 1e-09, "pt(8,25,ncp=10,log.p=TRUE)");  // R: pt(8,25,ncp=10,log.p=TRUE)
  check::rel(Dist::t(25, 10).p(12, false, true), -1.7752904868557264, 1e-09, "pt(12,25,ncp=10,lower.tail=FALSE,log.p=TRUE)");  // R: pt(12,25,ncp=10,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(25, 10).d(8, false), 0.1144847141667174, 1e-09, "dt(8,25,ncp=10)");  // R: dt(8,25,ncp=10)
  check::rel(Dist::t(25, 10).d(12, false), 0.11587186832456316, 1e-09, "dt(12,25,ncp=10)");  // R: dt(12,25,ncp=10)
  check::rel(Dist::t(7, 1.5).p(0, true, false), 0.06680720126885807, 1e-09, "pt(0,7,ncp=1.5)");  // R: pt(0,7,ncp=1.5)
  check::rel(Dist::t(7, 1.5).p(0, false, false), 0.9331927987311419, 1e-09, "pt(0,7,ncp=1.5,lower.tail=FALSE)");  // R: pt(0,7,ncp=1.5,lower.tail=FALSE)
  check::rel(Dist::t(7, 1.5).d(0, false), 0.12498842442456637, 1e-09, "dt(0,7,ncp=1.5)");  // R: dt(0,7,ncp=1.5)
  check::rel(Dist::t(7, 1.5).d(0, false), 0.12498842442456637, 1e-09, "dt(0,7,ncp=1.5)");  // R: dt(0,7,ncp=1.5)
  check::rel(Dist::t(kInf, 0.7).p(1.2, true, false), 0.691462461274013, 1e-09, "pt(1.2,Inf,ncp=0.7)");  // R: pt(1.2,Inf,ncp=0.7)
  check::rel(Dist::t(kInf, 0.7).p(1.2, false, true), -1.1759117615936185, 1e-09, "pt(1.2,Inf,ncp=0.7,lower.tail=FALSE,log.p=TRUE)");  // R: pt(1.2,Inf,ncp=0.7,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::t(kInf, 0.7).d(1.2, false), 0.3520653267642995, 1e-09, "dt(1.2,Inf,ncp=0.7)");  // R: dt(1.2,Inf,ncp=0.7)
  check::rel(Dist::t(kInf, 0.7).q(0.2, true, false), -0.14162123357291423, 1e-08, "qt(0.2,Inf,ncp=0.7)");  // R: qt(0.2,Inf,ncp=0.7)
  check::rel(Dist::t(3, 1).q(0.05, true, false), -0.7755096654293823, 1e-08, "qt(0.05,3,ncp=1)");  // R: qt(0.05,3,ncp=1)
  check::rel(Dist::t(3, 1).q(0.5, true, false), 1.091152546580517, 1e-08, "qt(0.5,3,ncp=1)");  // R: qt(0.5,3,ncp=1)
  check::rel(Dist::t(3, 1).q(0.9, true, false), 3.2973384698807138, 1e-08, "qt(0.9,3,ncp=1)");  // R: qt(0.9,3,ncp=1)
  check::rel(Dist::t(3, 1).q(0.2, false, false), 2.3273700504384607, 1e-08, "qt(0.2,3,ncp=1,lower.tail=FALSE)");  // R: qt(0.2,3,ncp=1,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).q(0.05, true, false), 0.35881009732044333, 1e-08, "qt(0.05,10,ncp=2)");  // R: qt(0.05,10,ncp=2)
  check::rel(Dist::t(10, 2).q(0.5, true, false), 2.053691151119523, 1e-08, "qt(0.5,10,ncp=2)");  // R: qt(0.5,10,ncp=2)
  check::rel(Dist::t(10, 2).q(0.9, true, false), 3.746633570263782, 1e-08, "qt(0.9,10,ncp=2)");  // R: qt(0.9,10,ncp=2)
  check::rel(Dist::t(10, 2).q(0.2, false, false), 3.0985180550852363, 1e-08, "qt(0.2,10,ncp=2,lower.tail=FALSE)");  // R: qt(0.2,10,ncp=2,lower.tail=FALSE)
  check::rel(Dist::t(10, -1.5).q(0.05, true, false), -3.689789189411755, 1e-08, "qt(0.05,10,ncp=-1.5)");  // R: qt(0.05,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).q(0.5, true, false), -1.5393384473886869, 1e-08, "qt(0.5,10,ncp=-1.5)");  // R: qt(0.5,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).q(0.9, true, false), -0.22238510472424622, 1e-08, "qt(0.9,10,ncp=-1.5)");  // R: qt(0.9,10,ncp=-1.5)
  check::rel(Dist::t(10, -1.5).q(0.2, false, false), -0.6657921758567511, 1e-08, "qt(0.2,10,ncp=-1.5,lower.tail=FALSE)");  // R: qt(0.2,10,ncp=-1.5,lower.tail=FALSE)
  check::rel(Dist::t(50, 3).q(0.05, true, false), 1.3470148589301658, 1e-08, "qt(0.05,50,ncp=3)");  // R: qt(0.05,50,ncp=3)
  check::rel(Dist::t(50, 3).q(0.5, true, false), 3.0154599326513534, 1e-08, "qt(0.5,50,ncp=3)");  // R: qt(0.5,50,ncp=3)
  check::rel(Dist::t(50, 3).q(0.9, true, false), 4.421952935126072, 1e-08, "qt(0.9,50,ncp=3)");  // R: qt(0.9,50,ncp=3)
  check::rel(Dist::t(50, 3).q(0.2, false, false), 3.923606311353282, 1e-08, "qt(0.2,50,ncp=3,lower.tail=FALSE)");  // R: qt(0.2,50,ncp=3,lower.tail=FALSE)
  check::rel(Dist::t(25, 10).q(0.05, true, false), 7.6377261033596255, 1e-08, "qt(0.05,25,ncp=10)");  // R: qt(0.05,25,ncp=10)
  check::rel(Dist::t(25, 10).q(0.5, true, false), 10.123862571849855, 1e-08, "qt(0.5,25,ncp=10)");  // R: qt(0.5,25,ncp=10)
  check::rel(Dist::t(25, 10).q(0.9, true, false), 12.736320921254531, 1e-08, "qt(0.9,25,ncp=10)");  // R: qt(0.9,25,ncp=10)
  check::rel(Dist::t(25, 10).q(0.2, false, false), 11.753148107720364, 1e-08, "qt(0.2,25,ncp=10,lower.tail=FALSE)");  // R: qt(0.2,25,ncp=10,lower.tail=FALSE)
  check::rel(Dist::t(1, 0.5).q(0.05, true, false), -3.0786368317735366, 1e-08, "qt(0.05,1,ncp=0.5)");  // R: qt(0.05,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).q(0.5, true, false), 0.6369699194979432, 1e-08, "qt(0.5,1,ncp=0.5)");  // R: qt(0.5,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).q(0.9, true, false), 5.485009833233306, 1e-08, "qt(0.9,1,ncp=0.5)");  // R: qt(0.9,1,ncp=0.5)
  check::rel(Dist::t(1, 0.5).q(0.2, false, false), 2.6160225389996725, 1e-08, "qt(0.2,1,ncp=0.5,lower.tail=FALSE)");  // R: qt(0.2,1,ncp=0.5,lower.tail=FALSE)
  check::rel(Dist::t(10, 2).q(-3, true, true), 0.3567441349903646, 1e-08, "qt(-3,10,ncp=2,log.p=TRUE)");  // R: qt(-3,10,ncp=2,log.p=TRUE)
  check::rel(Dist::t(10, 2).q(-0.1, false, true), 0.6921584153844549, 1e-08, "qt(-0.1,10,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: qt(-0.1,10,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(0.5, 0.1).p(0.01, true, false), 0.2788428323090405, 1e-09, "pchisq(0.01,0.5,ncp=0.1)");  // R: pchisq(0.01,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).p(0.01, false, false), 0.7211571676909595, 1e-09, "pchisq(0.01,0.5,ncp=0.1,lower.tail=FALSE)");  // R: pchisq(0.01,0.5,ncp=0.1,lower.tail=FALSE)
  check::rel(Dist::chisq(0.5, 0.1).p(1, true, false), 0.8187702252927079, 1e-09, "pchisq(1,0.5,ncp=0.1)");  // R: pchisq(1,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).p(1, false, false), 0.18122977470729199, 1e-09, "pchisq(1,0.5,ncp=0.1,lower.tail=FALSE)");  // R: pchisq(1,0.5,ncp=0.1,lower.tail=FALSE)
  check::rel(Dist::chisq(0.5, 0.1).p(10, true, false), 0.9989100575661441, 1e-09, "pchisq(10,0.5,ncp=0.1)");  // R: pchisq(10,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).p(10, false, false), 0.00108994243385588, 1e-09, "pchisq(10,0.5,ncp=0.1,lower.tail=FALSE)");  // R: pchisq(10,0.5,ncp=0.1,lower.tail=FALSE)
  check::rel(Dist::chisq(0.5, 0.1).p(0.01, true, true), -1.2771069809451467, 1e-09, "pchisq(0.01,0.5,ncp=0.1,log.p=TRUE)");  // R: pchisq(0.01,0.5,ncp=0.1,log.p=TRUE)
  check::rel(Dist::chisq(0.5, 0.1).p(10, false, true), -6.821630397111996, 1e-09, "pchisq(10,0.5,ncp=0.1,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(10,0.5,ncp=0.1,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(0.5, 0.1).d(0.01, false), 6.948792614065083, 1e-09, "dchisq(0.01,0.5,ncp=0.1)");  // R: dchisq(0.01,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).d(1, false), 0.1473290008970214, 1e-09, "dchisq(1,0.5,ncp=0.1)");  // R: dchisq(1,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).d(10, false), 0.0005561266944461309, 1e-09, "dchisq(10,0.5,ncp=0.1)");  // R: dchisq(10,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).q(0.01, true, false), 1.6488199333785078e-08, 1e-09, "qchisq(0.01,0.5,ncp=0.1)");  // R: qchisq(0.01,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).q(0.5, true, false), 0.10660191351855172, 1e-09, "qchisq(0.5,0.5,ncp=0.1)");  // R: qchisq(0.5,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).q(0.95, true, false), 2.9008937533671437, 1e-09, "qchisq(0.95,0.5,ncp=0.1)");  // R: qchisq(0.95,0.5,ncp=0.1)
  check::rel(Dist::chisq(0.5, 0.1).q(0.0001, false, false), 14.740929949424684, 1e-09, "qchisq(0.0001,0.5,ncp=0.1,lower.tail=FALSE)");  // R: qchisq(0.0001,0.5,ncp=0.1,lower.tail=FALSE)
  check::rel(Dist::chisq(1, 1).p(0.5, true, false), 0.3409007826677747, 1e-09, "pchisq(0.5,1,ncp=1)");  // R: pchisq(0.5,1,ncp=1)
  check::rel(Dist::chisq(1, 1).p(0.5, false, false), 0.6590992173322253, 1e-09, "pchisq(0.5,1,ncp=1,lower.tail=FALSE)");  // R: pchisq(0.5,1,ncp=1,lower.tail=FALSE)
  check::rel(Dist::chisq(1, 1).p(3, true, false), 0.7647841496310316, 1e-09, "pchisq(3,1,ncp=1)");  // R: pchisq(3,1,ncp=1)
  check::rel(Dist::chisq(1, 1).p(3, false, false), 0.23521585036896836, 1e-09, "pchisq(3,1,ncp=1,lower.tail=FALSE)");  // R: pchisq(3,1,ncp=1,lower.tail=FALSE)
  check::rel(Dist::chisq(1, 1).p(20, true, false), 0.9997418103789227, 1e-09, "pchisq(20,1,ncp=1)");  // R: pchisq(20,1,ncp=1)
  check::rel(Dist::chisq(1, 1).p(20, false, false), 0.0002581896210773188, 1e-09, "pchisq(20,1,ncp=1,lower.tail=FALSE)");  // R: pchisq(20,1,ncp=1,lower.tail=FALSE)
  check::rel(Dist::chisq(1, 1).p(0.5, true, true), -1.0761638039542192, 1e-09, "pchisq(0.5,1,ncp=1,log.p=TRUE)");  // R: pchisq(0.5,1,ncp=1,log.p=TRUE)
  check::rel(Dist::chisq(1, 1).p(20, false, true), -8.26181627758151, 1e-09, "pchisq(20,1,ncp=1,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(20,1,ncp=1,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(1, 1).d(0.5, false), 0.33595313069891597, 1e-09, "dchisq(0.5,1,ncp=1)");  // R: dchisq(0.5,1,ncp=1)
  check::rel(Dist::chisq(1, 1).d(3, false), 0.09085233082365742, 1e-09, "dchisq(3,1,ncp=1)");  // R: dchisq(3,1,ncp=1)
  check::rel(Dist::chisq(1, 1).d(20, false), 0.00010753594205112994, 1e-09, "dchisq(20,1,ncp=1)");  // R: dchisq(20,1,ncp=1)
  check::rel(Dist::chisq(1, 1).q(0.01, true, false), 0.0004269867137283649, 1e-09, "qchisq(0.01,1,ncp=1)");  // R: qchisq(0.01,1,ncp=1)
  check::rel(Dist::chisq(1, 1).q(0.5, true, false), 1.1036433113367627, 1e-09, "qchisq(0.5,1,ncp=1)");  // R: qchisq(0.5,1,ncp=1)
  check::rel(Dist::chisq(1, 1).q(0.95, true, false), 7.002086262339543, 1e-09, "qchisq(0.95,1,ncp=1)");  // R: qchisq(0.95,1,ncp=1)
  check::rel(Dist::chisq(1, 1).q(0.0001, false, false), 22.269244310730897, 1e-09, "qchisq(0.0001,1,ncp=1,lower.tail=FALSE)");  // R: qchisq(0.0001,1,ncp=1,lower.tail=FALSE)
  check::rel(Dist::chisq(3, 5).p(1, true, false), 0.025450001822339443, 1e-09, "pchisq(1,3,ncp=5)");  // R: pchisq(1,3,ncp=5)
  check::rel(Dist::chisq(3, 5).p(1, false, false), 0.9745499981776606, 1e-09, "pchisq(1,3,ncp=5,lower.tail=FALSE)");  // R: pchisq(1,3,ncp=5,lower.tail=FALSE)
  check::rel(Dist::chisq(3, 5).p(8, true, false), 0.5734922507093005, 1e-09, "pchisq(8,3,ncp=5)");  // R: pchisq(8,3,ncp=5)
  check::rel(Dist::chisq(3, 5).p(8, false, false), 0.42650774929069946, 1e-09, "pchisq(8,3,ncp=5,lower.tail=FALSE)");  // R: pchisq(8,3,ncp=5,lower.tail=FALSE)
  check::rel(Dist::chisq(3, 5).p(30, true, false), 0.9984709363618772, 1e-09, "pchisq(30,3,ncp=5)");  // R: pchisq(30,3,ncp=5)
  check::rel(Dist::chisq(3, 5).p(30, false, false), 0.0015290636381227606, 1e-09, "pchisq(30,3,ncp=5,lower.tail=FALSE)");  // R: pchisq(30,3,ncp=5,lower.tail=FALSE)
  check::rel(Dist::chisq(3, 5).p(1, true, true), -3.6710394643809146, 1e-09, "pchisq(1,3,ncp=5,log.p=TRUE)");  // R: pchisq(1,3,ncp=5,log.p=TRUE)
  check::rel(Dist::chisq(3, 5).p(30, false, true), -6.48309973215394, 1e-09, "pchisq(30,3,ncp=5,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(30,3,ncp=5,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(3, 5).d(1, false), 0.04108035197915606, 1e-09, "dchisq(1,3,ncp=5)");  // R: dchisq(1,3,ncp=5)
  check::rel(Dist::chisq(3, 5).d(8, false), 0.07485124312304355, 1e-09, "dchisq(8,3,ncp=5)");  // R: dchisq(8,3,ncp=5)
  check::rel(Dist::chisq(3, 5).d(30, false), 0.0004669184428665265, 1e-09, "dchisq(30,3,ncp=5)");  // R: dchisq(30,3,ncp=5)
  check::rel(Dist::chisq(3, 5).q(0.01, true, false), 0.5574215949255873, 1e-09, "qchisq(0.01,3,ncp=5)");  // R: qchisq(0.01,3,ncp=5)
  check::rel(Dist::chisq(3, 5).q(0.5, true, false), 7.061953885191706, 1e-09, "qchisq(0.5,3,ncp=5)");  // R: qchisq(0.5,3,ncp=5)
  check::rel(Dist::chisq(3, 5).q(0.95, true, false), 17.69461637828771, 1e-09, "qchisq(0.95,3,ncp=5)");  // R: qchisq(0.95,3,ncp=5)
  check::rel(Dist::chisq(3, 5).q(0.0001, false, false), 38.61736992444639, 1e-09, "qchisq(0.0001,3,ncp=5,lower.tail=FALSE)");  // R: qchisq(0.0001,3,ncp=5,lower.tail=FALSE)
  check::rel(Dist::chisq(10, 37).p(20, true, false), 0.0066392829688228, 1e-09, "pchisq(20,10,ncp=37)");  // R: pchisq(20,10,ncp=37)
  check::rel(Dist::chisq(10, 37).p(20, false, false), 0.9933607170311772, 1e-09, "pchisq(20,10,ncp=37,lower.tail=FALSE)");  // R: pchisq(20,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::chisq(10, 37).p(47, true, false), 0.5296418293078341, 1e-09, "pchisq(47,10,ncp=37)");  // R: pchisq(47,10,ncp=37)
  check::rel(Dist::chisq(10, 37).p(47, false, false), 0.4703581706921659, 1e-09, "pchisq(47,10,ncp=37,lower.tail=FALSE)");  // R: pchisq(47,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::chisq(10, 37).p(90, true, false), 0.9975674496228606, 1e-09, "pchisq(90,10,ncp=37)");  // R: pchisq(90,10,ncp=37)
  check::rel(Dist::chisq(10, 37).p(90, false, false), 0.0024325503771393643, 1e-09, "pchisq(90,10,ncp=37,lower.tail=FALSE)");  // R: pchisq(90,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::chisq(10, 37).p(20, true, true), -5.014751307947842, 1e-09, "pchisq(20,10,ncp=37,log.p=TRUE)");  // R: pchisq(20,10,ncp=37,log.p=TRUE)
  check::rel(Dist::chisq(10, 37).p(90, false, true), -6.018815034092103, 1e-09, "pchisq(90,10,ncp=37,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(90,10,ncp=37,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(10, 37).d(20, false), 0.0022753416751545353, 1e-09, "dchisq(20,10,ncp=37)");  // R: dchisq(20,10,ncp=37)
  check::rel(Dist::chisq(10, 37).d(47, false), 0.030544084340896584, 1e-09, "dchisq(47,10,ncp=37)");  // R: dchisq(47,10,ncp=37)
  check::rel(Dist::chisq(10, 37).d(90, false), 0.00041242920369497345, 1e-09, "dchisq(90,10,ncp=37)");  // R: dchisq(90,10,ncp=37)
  check::rel(Dist::chisq(10, 37).q(0.01, true, false), 21.245833890302094, 1e-09, "qchisq(0.01,10,ncp=37)");  // R: qchisq(0.01,10,ncp=37)
  check::rel(Dist::chisq(10, 37).q(0.5, true, false), 46.03666119185563, 1e-09, "qchisq(0.5,10,ncp=37)");  // R: qchisq(0.5,10,ncp=37)
  check::rel(Dist::chisq(10, 37).q(0.95, true, false), 69.84351712558431, 1e-09, "qchisq(0.95,10,ncp=37)");  // R: qchisq(0.95,10,ncp=37)
  check::rel(Dist::chisq(10, 37).q(0.0001, false, false), 107.37960238130573, 1e-09, "qchisq(0.0001,10,ncp=37,lower.tail=FALSE)");  // R: qchisq(0.0001,10,ncp=37,lower.tail=FALSE)
  check::rel(Dist::chisq(100, 200).p(250, true, false), 0.0513734231090317, 1e-09, "pchisq(250,100,ncp=200)");  // R: pchisq(250,100,ncp=200)
  check::rel(Dist::chisq(100, 200).p(250, false, false), 0.9486265768909683, 1e-09, "pchisq(250,100,ncp=200,lower.tail=FALSE)");  // R: pchisq(250,100,ncp=200,lower.tail=FALSE)
  check::rel(Dist::chisq(100, 200).p(300, true, false), 0.511779823794299, 1e-09, "pchisq(300,100,ncp=200)");  // R: pchisq(300,100,ncp=200)
  check::rel(Dist::chisq(100, 200).p(300, false, false), 0.488220176205701, 1e-09, "pchisq(300,100,ncp=200,lower.tail=FALSE)");  // R: pchisq(300,100,ncp=200,lower.tail=FALSE)
  check::rel(Dist::chisq(100, 200).p(380, true, false), 0.9916244413097184, 1e-09, "pchisq(380,100,ncp=200)");  // R: pchisq(380,100,ncp=200)
  check::rel(Dist::chisq(100, 200).p(380, false, false), 0.008375558690281615, 1e-09, "pchisq(380,100,ncp=200,lower.tail=FALSE)");  // R: pchisq(380,100,ncp=200,lower.tail=FALSE)
  check::rel(Dist::chisq(100, 200).p(250, true, true), -2.9686343003782762, 1e-09, "pchisq(250,100,ncp=200,log.p=TRUE)");  // R: pchisq(250,100,ncp=200,log.p=TRUE)
  check::rel(Dist::chisq(100, 200).p(380, false, true), -4.782437494208615, 1e-09, "pchisq(380,100,ncp=200,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(380,100,ncp=200,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(100, 200).d(250, false), 0.0036902758810205675, 1e-09, "dchisq(250,100,ncp=200)");  // R: dchisq(250,100,ncp=200)
  check::rel(Dist::chisq(100, 200).d(300, false), 0.012601350317345012, 1e-09, "dchisq(300,100,ncp=200)");  // R: dchisq(300,100,ncp=200)
  check::rel(Dist::chisq(100, 200).d(380, false), 0.0006324065726016581, 1e-09, "dchisq(380,100,ncp=200)");  // R: dchisq(380,100,ncp=200)
  check::rel(Dist::chisq(100, 200).q(0.01, true, false), 230.60705430144574, 1e-09, "qchisq(0.01,100,ncp=200)");  // R: qchisq(0.01,100,ncp=200)
  check::rel(Dist::chisq(100, 200).q(0.5, true, false), 299.06627982828905, 1e-09, "qchisq(0.5,100,ncp=200)");  // R: qchisq(0.5,100,ncp=200)
  check::rel(Dist::chisq(100, 200).q(0.95, true, false), 353.56112982719253, 1e-09, "qchisq(0.95,100,ncp=200)");  // R: qchisq(0.95,100,ncp=200)
  check::rel(Dist::chisq(100, 200).q(0.0001, false, false), 429.54962189605897, 1e-09, "qchisq(0.0001,100,ncp=200,lower.tail=FALSE)");  // R: qchisq(0.0001,100,ncp=200,lower.tail=FALSE)
  check::rel(Dist::chisq(0, 2).p(0.1, true, false), 0.3860447697801761, 1e-09, "pchisq(0.1,0,ncp=2)");  // R: pchisq(0.1,0,ncp=2)
  check::rel(Dist::chisq(0, 2).p(0.1, false, false), 0.6139552302198239, 1e-09, "pchisq(0.1,0,ncp=2,lower.tail=FALSE)");  // R: pchisq(0.1,0,ncp=2,lower.tail=FALSE)
  check::rel(Dist::chisq(0, 2).p(2, true, false), 0.6542541612768356, 1e-09, "pchisq(2,0,ncp=2)");  // R: pchisq(2,0,ncp=2)
  check::rel(Dist::chisq(0, 2).p(2, false, false), 0.34574583872316444, 1e-09, "pchisq(2,0,ncp=2,lower.tail=FALSE)");  // R: pchisq(2,0,ncp=2,lower.tail=FALSE)
  check::rel(Dist::chisq(0, 2).p(10, true, false), 0.9766500547706445, 1e-09, "pchisq(10,0,ncp=2)");  // R: pchisq(10,0,ncp=2)
  check::rel(Dist::chisq(0, 2).p(10, false, false), 0.02334994522935553, 1e-09, "pchisq(10,0,ncp=2,lower.tail=FALSE)");  // R: pchisq(10,0,ncp=2,lower.tail=FALSE)
  check::rel(Dist::chisq(0, 2).p(0.1, true, true), -0.9518019323564506, 1e-09, "pchisq(0.1,0,ncp=2,log.p=TRUE)");  // R: pchisq(0.1,0,ncp=2,log.p=TRUE)
  check::rel(Dist::chisq(0, 2).p(10, false, true), -3.7571606405078604, 1e-09, "pchisq(10,0,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(10,0,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(0, 2).d(0.1, false), 0.17937970053137609, 1e-09, "dchisq(0.1,0,ncp=2)");  // R: dchisq(0.1,0,ncp=2)
  check::rel(Dist::chisq(0, 2).d(2, false), 0.10763464462446837, 1e-09, "dchisq(2,0,ncp=2)");  // R: dchisq(2,0,ncp=2)
  check::rel(Dist::chisq(0, 2).d(10, false), 0.008315333956012472, 1e-09, "dchisq(10,0,ncp=2)");  // R: dchisq(10,0,ncp=2)
  check::rel(Dist::chisq(0, 2).q(0.5, true, false), 0.7934451320402568, 1e-09, "qchisq(0.5,0,ncp=2)");  // R: qchisq(0.5,0,ncp=2)
  check::rel(Dist::chisq(0, 2).q(0.95, true, false), 7.836073876991101, 1e-09, "qchisq(0.95,0,ncp=2)");  // R: qchisq(0.95,0,ncp=2)
  check::rel(Dist::chisq(0, 2).q(0.0001, false, false), 24.51279047135543, 1e-09, "qchisq(0.0001,0,ncp=2,lower.tail=FALSE)");  // R: qchisq(0.0001,0,ncp=2,lower.tail=FALSE)
  check::rel(Dist::chisq(2, 3).p(0.2, true, false), 0.02283893536602521, 1e-09, "pchisq(0.2,2,ncp=3)");  // R: pchisq(0.2,2,ncp=3)
  check::rel(Dist::chisq(2, 3).p(0.2, false, false), 0.9771610646339748, 1e-09, "pchisq(0.2,2,ncp=3,lower.tail=FALSE)");  // R: pchisq(0.2,2,ncp=3,lower.tail=FALSE)
  check::rel(Dist::chisq(2, 3).p(5, true, false), 0.5940608030781969, 1e-09, "pchisq(5,2,ncp=3)");  // R: pchisq(5,2,ncp=3)
  check::rel(Dist::chisq(2, 3).p(5, false, false), 0.40593919692180314, 1e-09, "pchisq(5,2,ncp=3,lower.tail=FALSE)");  // R: pchisq(5,2,ncp=3,lower.tail=FALSE)
  check::rel(Dist::chisq(2, 3).p(0.2, true, true), -3.779288507976652, 1e-09, "pchisq(0.2,2,ncp=3,log.p=TRUE)");  // R: pchisq(0.2,2,ncp=3,log.p=TRUE)
  check::rel(Dist::chisq(2, 3).p(5, false, true), -0.9015518918720911, 1e-09, "pchisq(5,2,ncp=3,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(5,2,ncp=3,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(2, 3).d(0.2, false), 0.1166678849616928, 1e-09, "dchisq(0.2,2,ncp=3)");  // R: dchisq(0.2,2,ncp=3)
  check::rel(Dist::chisq(2, 3).d(5, false), 0.09277869381135642, 1e-09, "dchisq(5,2,ncp=3)");  // R: dchisq(5,2,ncp=3)
  check::rel(Dist::chisq(2, 3).d(5, false), 0.09277869381135642, 1e-09, "dchisq(5,2,ncp=3)");  // R: dchisq(5,2,ncp=3)
  check::rel(Dist::chisq(2, 3).q(0.01, true, false), 0.08867605289078992, 1e-09, "qchisq(0.01,2,ncp=3)");  // R: qchisq(0.01,2,ncp=3)
  check::rel(Dist::chisq(2, 3).q(0.5, true, false), 4.059801668031781, 1e-09, "qchisq(0.5,2,ncp=3)");  // R: qchisq(0.5,2,ncp=3)
  check::rel(Dist::chisq(2, 3).q(0.95, true, false), 12.806398856206382, 1e-09, "qchisq(0.95,2,ncp=3)");  // R: qchisq(0.95,2,ncp=3)
  check::rel(Dist::chisq(2, 3).q(0.0001, false, false), 31.415284526968637, 1e-09, "qchisq(0.0001,2,ncp=3,lower.tail=FALSE)");  // R: qchisq(0.0001,2,ncp=3,lower.tail=FALSE)
  check::rel(Dist::chisq(0, 2).p(0, true, false), 0.36787944117144233, 1e-09, "pchisq(0,0,ncp=2)");  // R: pchisq(0,0,ncp=2)
  check::rel(Dist::chisq(0, 2).p(0, false, true), -0.45867514538708193, 1e-09, "pchisq(0,0,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: pchisq(0,0,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::chisq(2, 3).d(0, false), 0.11156508007421491, 1e-09, "dchisq(0,2,ncp=3)");  // R: dchisq(0,2,ncp=3)
  check::rel(Dist::chisq(3, 10000.0).q(-0.001, true, true), 10629.657233054342, 1e-09, "qchisq(-0.001,3,ncp=10000.0,log.p=TRUE)");  // R: qchisq(-0.001,3,ncp=10000.0,log.p=TRUE)
  check::rel(Dist::chisq(3, 10000.0).p(10000.0, true, false), 0.49601057719790786, 1e-09, "pchisq(10000.0,3,ncp=10000.0)");  // R: pchisq(10000.0,3,ncp=10000.0)
  check::rel(Dist::chisq(3, 10000.0).p(10000.0, false, false), 0.5039894228020921, 1e-09, "pchisq(10000.0,3,ncp=10000.0,lower.tail=FALSE)");  // R: pchisq(10000.0,3,ncp=10000.0,lower.tail=FALSE)
  check::rel(Dist::chisq(3, 10000.0).q(0.01, true, false), 9542.118907944698, 1e-09, "qchisq(0.01,3,ncp=10000.0)");  // R: qchisq(0.01,3,ncp=10000.0)
  check::rel(Dist::f(1, 5, 0.5).p(0.1, false, false), 0.8148838928763276, 1e-09, "pf(0.1,1,5,ncp=0.5,lower.tail=FALSE)");  // R: pf(0.1,1,5,ncp=0.5,lower.tail=FALSE)
  check::rel(Dist::f(1, 5, 0.5).p(10, true, false), 0.9518680437606696, 1e-09, "pf(10,1,5,ncp=0.5)");  // R: pf(10,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).p(0.1, true, true), -1.6867720446922354, 1e-09, "pf(0.1,1,5,ncp=0.5,log.p=TRUE)");  // R: pf(0.1,1,5,ncp=0.5,log.p=TRUE)
  check::rel(Dist::f(1, 5, 0.5).d(0.1, false), 0.9070478696927038, 1e-09, "df(0.1,1,5,ncp=0.5)");  // R: df(0.1,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).d(2, false), 0.11209586474498372, 1e-09, "df(2,1,5,ncp=0.5)");  // R: df(2,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).d(10, false), 0.007786494801945317, 1e-09, "df(10,1,5,ncp=0.5)");  // R: df(10,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).q(0.05, true, false), 0.007161130613414133, 1e-08, "qf(0.05,1,5,ncp=0.5)");  // R: qf(0.05,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).q(0.5, true, false), 0.8457714506430282, 1e-08, "qf(0.5,1,5,ncp=0.5)");  // R: qf(0.5,1,5,ncp=0.5)
  check::rel(Dist::f(1, 5, 0.5).q(0.9, true, false), 6.091760782448988, 1e-08, "qf(0.9,1,5,ncp=0.5)");  // R: qf(0.9,1,5,ncp=0.5)
  check::rel(Dist::f(3, 10, 4).p(0.5, false, false), 0.9231732441911956, 1e-09, "pf(0.5,3,10,ncp=4,lower.tail=FALSE)");  // R: pf(0.5,3,10,ncp=4,lower.tail=FALSE)
  check::rel(Dist::f(3, 10, 4).d(0.5, false), 0.2299913002153038, 1e-09, "df(0.5,3,10,ncp=4)");  // R: df(0.5,3,10,ncp=4)
  check::rel(Dist::f(3, 10, 4).d(2, false), 0.22713302262485235, 1e-09, "df(2,3,10,ncp=4)");  // R: df(2,3,10,ncp=4)
  check::rel(Dist::f(3, 10, 4).d(8, false), 0.016577779547626136, 1e-09, "df(8,3,10,ncp=4)");  // R: df(8,3,10,ncp=4)
  check::rel(Dist::f(3, 10, 4).q(0.05, true, false), 0.37627831793751626, 1e-08, "qf(0.05,3,10,ncp=4)");  // R: qf(0.05,3,10,ncp=4)
  check::rel(Dist::f(3, 10, 4).q(0.5, true, false), 2.152168771403392, 1e-08, "qf(0.5,3,10,ncp=4)");  // R: qf(0.5,3,10,ncp=4)
  check::rel(Dist::f(3, 10, 4).q(0.9, true, false), 6.016929090813807, 1e-08, "qf(0.9,3,10,ncp=4)");  // R: qf(0.9,3,10,ncp=4)
  check::rel(Dist::f(10, 30, 37).d(2, false), 0.04816253533468992, 1e-09, "df(2,10,30,ncp=37)");  // R: df(2,10,30,ncp=37)
  check::rel(Dist::f(10, 30, 37).d(4, false), 0.23621397687347864, 1e-09, "df(4,10,30,ncp=37)");  // R: df(4,10,30,ncp=37)
  check::rel(Dist::f(10, 30, 37).d(10, false), 0.014368680304645045, 1e-09, "df(10,10,30,ncp=37)");  // R: df(10,10,30,ncp=37)
  check::rel(Dist::f(10, 30, 37).q(0.05, true, false), 2.4513032849730685, 1e-08, "qf(0.05,10,30,ncp=37)");  // R: qf(0.05,10,30,ncp=37)
  check::rel(Dist::f(10, 30, 37).q(0.5, true, false), 4.696438417913255, 1e-08, "qf(0.5,10,30,ncp=37)");  // R: qf(0.5,10,30,ncp=37)
  check::rel(Dist::f(10, 30, 37).q(0.9, true, false), 7.635814355872371, 1e-08, "qf(0.9,10,30,ncp=37)");  // R: qf(0.9,10,30,ncp=37)
  check::rel(Dist::f(0.5, 2, 1).p(0.05, true, false), 0.20342875183793616, 1e-09, "pf(0.05,0.5,2,ncp=1)");  // R: pf(0.05,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).p(0.05, false, false), 0.7965712481620638, 1e-09, "pf(0.05,0.5,2,ncp=1,lower.tail=FALSE)");  // R: pf(0.05,0.5,2,ncp=1,lower.tail=FALSE)
  check::rel(Dist::f(0.5, 2, 1).p(0.05, true, true), -1.5924394492221352, 1e-09, "pf(0.05,0.5,2,ncp=1,log.p=TRUE)");  // R: pf(0.05,0.5,2,ncp=1,log.p=TRUE)
  check::rel(Dist::f(0.5, 2, 1).d(0.05, false), 1.0293910320722524, 1e-09, "df(0.05,0.5,2,ncp=1)");  // R: df(0.05,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).d(1, false), 0.12551560896495112, 1e-09, "df(1,0.5,2,ncp=1)");  // R: df(1,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).d(20, false), 0.0048836100248081825, 1e-09, "df(20,0.5,2,ncp=1)");  // R: df(20,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).q(0.05, true, false), 0.00018471787289098512, 1e-08, "qf(0.05,0.5,2,ncp=1)");  // R: qf(0.05,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).q(0.5, true, false), 1.4751742148885694, 1e-08, "qf(0.5,0.5,2,ncp=1)");  // R: qf(0.5,0.5,2,ncp=1)
  check::rel(Dist::f(0.5, 2, 1).q(0.9, true, false), 25.190235813980838, 5e-08, "qf(0.9,0.5,2,ncp=1)");  // R: qf(0.9,0.5,2,ncp=1) (R's qf with ncp solves to ~1e-9 in p: pfref(R's q) - 0.9 = 9e-10)
  check::rel(Dist::f(4, 0.7, 2).d(1, false), 0.17211539231677028, 1e-09, "df(1,4,0.7,ncp=2)");  // R: df(1,4,0.7,ncp=2)
  check::rel(Dist::f(4, 0.7, 2).d(50, false), 0.0014942504905974513, 1e-09, "df(50,4,0.7,ncp=2)");  // R: df(50,4,0.7,ncp=2)
  check::rel(Dist::f(4, 0.7, 2).d(50, false), 0.0014942504905974513, 1e-09, "df(50,4,0.7,ncp=2)");  // R: df(50,4,0.7,ncp=2)
  check::rel(Dist::f(4, 0.7, 2).q(0.05, true, false), 0.19403313516662182, 1e-08, "qf(0.05,4,0.7,ncp=2)");  // R: qf(0.05,4,0.7,ncp=2)
  check::rel(Dist::f(4, 0.7, 2).q(0.5, true, false), 4.091708242292133, 1e-08, "qf(0.5,4,0.7,ncp=2)");  // R: qf(0.5,4,0.7,ncp=2)
  check::rel(Dist::f(4, 0.7, 2).q(0.9, true, false), 451.35636651780777, 5e-08, "qf(0.9,4,0.7,ncp=2)");  // R: qf(0.9,4,0.7,ncp=2) (R's qf with ncp solves to ~1e-9 in p: pfref(R's q) - 0.9 = 9e-10)
  check::rel(Dist::f(5, kInf, 2).p(0.5, true, false), 0.11282703768006931, 1e-09, "pf(0.5,5,Inf,ncp=2)");  // R: pf(0.5,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).p(0.5, false, false), 0.8871729623199307, 1e-09, "pf(0.5,5,Inf,ncp=2,lower.tail=FALSE)");  // R: pf(0.5,5,Inf,ncp=2,lower.tail=FALSE)
  check::rel(Dist::f(5, kInf, 2).p(1, true, false), 0.3750814647441277, 1e-09, "pf(1,5,Inf,ncp=2)");  // R: pf(1,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).p(1, false, false), 0.6249185352558723, 1e-09, "pf(1,5,Inf,ncp=2,lower.tail=FALSE)");  // R: pf(1,5,Inf,ncp=2,lower.tail=FALSE)
  check::rel(Dist::f(5, kInf, 2).p(3, true, false), 0.9485676384659193, 1e-09, "pf(3,5,Inf,ncp=2)");  // R: pf(3,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).p(3, false, false), 0.05143236153408065, 1e-09, "pf(3,5,Inf,ncp=2,lower.tail=FALSE)");  // R: pf(3,5,Inf,ncp=2,lower.tail=FALSE)
  check::rel(Dist::f(5, kInf, 2).p(0.5, true, true), -2.181899272895319, 1e-09, "pf(0.5,5,Inf,ncp=2,log.p=TRUE)");  // R: pf(0.5,5,Inf,ncp=2,log.p=TRUE)
  check::rel(Dist::f(5, kInf, 2).p(3, false, true), -2.9674877028067352, 1e-09, "pf(3,5,Inf,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: pf(3,5,Inf,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::f(5, kInf, 2).d(0.5, false), 0.44268440258768593, 1e-09, "df(0.5,5,Inf,ncp=2)");  // R: df(0.5,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).d(1, false), 0.5458092883475683, 1e-09, "df(1,5,Inf,ncp=2)");  // R: df(1,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).d(3, false), 0.07683352999137283, 1e-09, "df(3,5,Inf,ncp=2)");  // R: df(3,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).q(0.05, true, false), 0.33648255592236415, 1e-08, "qf(0.05,5,Inf,ncp=2)");  // R: qf(0.05,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).q(0.5, true, false), 1.2375710085883065, 1e-08, "qf(0.5,5,Inf,ncp=2)");  // R: qf(0.5,5,Inf,ncp=2)
  check::rel(Dist::f(5, kInf, 2).q(0.9, true, false), 2.5424652531559246, 1e-08, "qf(0.9,5,Inf,ncp=2)");  // R: qf(0.9,5,Inf,ncp=2)
  check::rel(Dist::f(20, 100, 10).d(1, false), 0.5953264055956622, 1e-09, "df(1,20,100,ncp=10)");  // R: df(1,20,100,ncp=10)
  check::rel(Dist::f(20, 100, 10).d(1.5, false), 0.7985594272819241, 1e-09, "df(1.5,20,100,ncp=10)");  // R: df(1.5,20,100,ncp=10)
  check::rel(Dist::f(20, 100, 10).d(3, false), 0.03145773630831058, 1e-09, "df(3,20,100,ncp=10)");  // R: df(3,20,100,ncp=10)
  check::rel(Dist::f(20, 100, 10).q(0.05, true, false), 0.8125927489662308, 1e-08, "qf(0.05,20,100,ncp=10)");  // R: qf(0.05,20,100,ncp=10)
  check::rel(Dist::f(20, 100, 10).q(0.5, true, false), 1.4675926172371812, 1e-08, "qf(0.5,20,100,ncp=10)");  // R: qf(0.5,20,100,ncp=10)
  check::rel(Dist::f(20, 100, 10).q(0.9, true, false), 2.207665871957964, 1e-08, "qf(0.9,20,100,ncp=10)");  // R: qf(0.9,20,100,ncp=10)
  check::rel(Dist::f(2, 5, 1).d(0, false), 0.6065306597126335, 1e-09, "df(0,2,5,ncp=1)");  // R: df(0,2,5,ncp=1)
  check::rel(Dist::f(3, 5, 1).d(0, false), 0.0, 1e-09, "df(0,3,5,ncp=1)");  // R: df(0,3,5,ncp=1)
  check::rel(Dist::f(3, 4, 2).q(-2, false, true), 5.5996341710375805, 1e-08, "qf(-2,3,4,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // R: qf(-2,3,4,ncp=2,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::tukey(3, kInf, 1).p(1, true, false), 0.24071264122928304, 1e-08, "ptukey(1,3,Inf,1)");  // R: ptukey(1,3,Inf,1)
  check::rel(Dist::tukey(3, kInf, 1).p(3.3, true, false), 0.9486890935810137, 1e-08, "ptukey(3.3,3,Inf,1)");  // R: ptukey(3.3,3,Inf,1)
  check::rel(Dist::tukey(3, kInf, 1).p(5, true, false), 0.9988160087690132, 1e-08, "ptukey(5,3,Inf,1)");  // R: ptukey(5,3,Inf,1)
  check::rel(Dist::tukey(3, kInf, 1).p(5, false, false), 0.0011839912309867806, 1e-08, "ptukey(5,3,Inf,1,lower.tail=FALSE)");  // R: ptukey(5,3,Inf,1,lower.tail=FALSE)
  check::rel(Dist::tukey(5, kInf, 1).p(2, true, false), 0.38155051944800633, 1e-08, "ptukey(2,5,Inf,1)");  // R: ptukey(2,5,Inf,1)
  check::rel(Dist::tukey(5, kInf, 1).p(4, true, false), 0.9623039317799011, 1e-08, "ptukey(4,5,Inf,1)");  // R: ptukey(4,5,Inf,1)
  check::rel(Dist::tukey(5, kInf, 1).p(4, false, false), 0.037696068220098855, 1e-08, "ptukey(4,5,Inf,1,lower.tail=FALSE)");  // R: ptukey(4,5,Inf,1,lower.tail=FALSE)
  check::rel(Dist::tukey(10, kInf, 2).p(4.5, true, false), 0.9076261866993814, 1e-08, "ptukey(4.5,10,Inf,2)");  // R: ptukey(4.5,10,Inf,2)
  check::rel(Dist::tukey(10, kInf, 2).p(4.5, false, false), 0.09237381330061856, 1e-08, "ptukey(4.5,10,Inf,2,lower.tail=FALSE)");  // R: ptukey(4.5,10,Inf,2,lower.tail=FALSE)
  check::rel(Dist::tukey(3, 10, 1).p(1, true, false), 0.23510918942120207, 1e-08, "ptukey(1,3,10,1)");  // R: ptukey(1,3,10,1)
  check::rel(Dist::tukey(3, 10, 1).p(3, true, false), 0.8650165848100901, 1e-08, "ptukey(3,3,10,1)");  // R: ptukey(3,3,10,1)
  check::rel(Dist::tukey(3, 10, 1).p(5, true, false), 0.9863843190832841, 1e-08, "ptukey(5,3,10,1)");  // R: ptukey(5,3,10,1)
  check::rel(Dist::tukey(3, 10, 1).p(5, false, false), 0.013615680916715944, 1e-08, "ptukey(5,3,10,1,lower.tail=FALSE)");  // R: ptukey(5,3,10,1,lower.tail=FALSE)
  check::rel(Dist::tukey(4, 20, 1).p(3.9, true, false), 0.9456231021833177, 1e-08, "ptukey(3.9,4,20,1)");  // R: ptukey(3.9,4,20,1)
  check::rel(Dist::tukey(4, 20, 1).p(3.9, false, false), 0.05437689781668231, 1e-08, "ptukey(3.9,4,20,1,lower.tail=FALSE)");  // R: ptukey(3.9,4,20,1,lower.tail=FALSE)
  check::rel(Dist::tukey(5, 60, 1).p(2, true, false), 0.37892713546238044, 1e-08, "ptukey(2,5,60,1)");  // R: ptukey(2,5,60,1)
  check::rel(Dist::tukey(5, 60, 1).p(4.4, true, false), 0.9770192192284746, 1e-08, "ptukey(4.4,5,60,1)");  // R: ptukey(4.4,5,60,1)
  check::rel(Dist::tukey(10, 120, 1).p(4.5, true, false), 0.9437760987476767, 1e-08, "ptukey(4.5,10,120,1)");  // R: ptukey(4.5,10,120,1)
  check::rel(Dist::tukey(10, 120, 1).p(6, true, false), 0.9982658799569858, 1e-08, "ptukey(6,10,120,1)");  // R: ptukey(6,10,120,1)
  check::rel(Dist::tukey(10, 120, 1).p(6, false, false), 0.001734120043014209, 1e-08, "ptukey(6,10,120,1,lower.tail=FALSE)");  // R: ptukey(6,10,120,1,lower.tail=FALSE)
  check::rel(Dist::tukey(20, 10000.0, 1).p(5, true, false), 0.9484756317707885, 1e-08, "ptukey(5,20,10000.0,1)");  // R: ptukey(5,20,10000.0,1)
  check::rel(Dist::tukey(3, 30, 3).p(3.5, true, false), 0.8645107916154758, 1e-08, "ptukey(3.5,3,30,3)");  // R: ptukey(3.5,3,30,3)
  check::rel(Dist::tukey(3, 30, 3).p(3.5, false, false), 0.13548920838452416, 1e-08, "ptukey(3.5,3,30,3,lower.tail=FALSE)");  // R: ptukey(3.5,3,30,3,lower.tail=FALSE)
  check::rel(Dist::tukey(6, 15, 2).p(4, true, false), 0.8126721241659863, 1e-08, "ptukey(4,6,15,2)");  // R: ptukey(4,6,15,2)
  check::rel(Dist::tukey(6, 15, 2).p(6, true, false), 0.985566394125999, 1e-08, "ptukey(6,6,15,2)");  // R: ptukey(6,6,15,2)
  check::rel(Dist::tukey(6, 15, 2).p(6, false, false), 0.014433605874000999, 1e-08, "ptukey(6,6,15,2,lower.tail=FALSE)");  // R: ptukey(6,6,15,2,lower.tail=FALSE)
  check::rel(Dist::tukey(3, 10, 1).p(3, true, true), -0.14500659904320534, 1e-08, "ptukey(3,3,10,1,log.p=TRUE)");  // R: ptukey(3,3,10,1,log.p=TRUE)
  check::rel(Dist::tukey(3, 10, 1).p(5, false, true), -4.296533141861965, 1e-08, "ptukey(5,3,10,1,lower.tail=FALSE,log.p=TRUE)");  // R: ptukey(5,3,10,1,lower.tail=FALSE,log.p=TRUE)
  check::rel(Dist::tukey(3, kInf, 1).q(0.95, true, false), 3.3144931575530063, 0.0001, "qtukey(0.95,3,Inf,1)");  // R: qtukey(0.95,3,Inf,1)
  check::rel(Dist::tukey(3, 10, 1).q(0.95, true, false), 3.8767767491915737, 0.0001, "qtukey(0.95,3,10,1)");  // R: qtukey(0.95,3,10,1)
  check::rel(Dist::tukey(5, 60, 1).q(0.5, true, false), 2.2701236938497984, 0.0001, "qtukey(0.5,5,60,1)");  // R: qtukey(0.5,5,60,1)
  check::rel(Dist::tukey(10, 120, 1).q(0.99, true, false), 5.299213945252483, 0.0001, "qtukey(0.99,10,120,1)");  // R: qtukey(0.99,10,120,1)
  check::rel(Dist::tukey(4, 20, 2).q(0.9, true, false), 3.915294940258141, 0.0001, "qtukey(0.9,4,20,2)");  // R: qtukey(0.9,4,20,2)
}

void reference_cases() {
  check::rel(Dist::t(3, -2).p(4, false, false), 0.00010762805262109493, 1e-10, "pt(4,3,ncp=-2,lower.tail=FALSE)");  // ref: exp(ptref(4,3,-2,FALSE))  (R: pt(4,3,ncp=-2,lower.tail=FALSE) = 0.00010762805281649968)
  check::rel(Dist::t(3, -2).p(4, false, true), -9.136829232124121, 1e-10, "pt(4,3,ncp=-2,lower.tail=FALSE,log.p=TRUE)");  // ref: ptref(4,3,-2,FALSE)  (R: pt(4,3,ncp=-2,lower.tail=FALSE,log.p=TRUE) = -9.136829230308566)
  check::rel(Dist::t(10, 2).d(-1, false), 0.004809024120553366, 1e-10, "dt(-1,10,ncp=2)");  // ref: exp(dtref(-1,10,2))  (R: dt(-1,10,ncp=2) = 0.0048090241220455346)
  check::rel(Dist::t(1000, 1).d(0.5, false), 0.35198282067492026, 1e-10, "dt(0.5,1000,ncp=1)");  // ref: exp(dtref(0.5,1000,1))  (R: dt(0.5,1000,ncp=1) = 0.3519828205709172)
  check::rel(Dist::t(10, 10).p(-1, true, false), 2.6351020387524e-27, 1e-10, "pt(-1,10,ncp=10)");  // ref: exp(ptref(-1,10,10,TRUE))  (R: pt(-1,10,ncp=10) = 0.0)
  check::rel(Dist::t(50, 37).p(8, true, false), 6.318051647727921e-101, 1e-10, "pt(8,50,ncp=37)");  // ref: exp(ptref(8,50,37,TRUE))  (R: pt(8,50,ncp=37) = 1.2213222328608778e-166)
  check::rel(Dist::t(50, 37).p(60, false, false), 2.9513085244114688e-05, 1e-10, "pt(60,50,ncp=37,lower.tail=FALSE)");  // ref: exp(ptref(60,50,37,FALSE))  (R: pt(60,50,ncp=37,lower.tail=FALSE) = 2.9513086169941083e-05)
  check::rel(Dist::t(10, 1).p(40, false, false), 2.2663166945001668e-11, 1e-10, "pt(40,10,ncp=1,lower.tail=FALSE)");  // ref: exp(ptref(40,10,1,FALSE))  (R: pt(40,10,ncp=1,lower.tail=FALSE) = 2.2855717318748248e-11)
  check::rel(Dist::t(50, 10).p(-5, true, true), -97.04423689496527, 1e-10, "pt(-5,50,ncp=10,log.p=TRUE)");  // ref: ptref(-5,50,10,TRUE)  (R: pt(-5,50,ncp=10,log.p=TRUE) = -31.056627960660034)
  check::rel(Dist::t(50, 1).d(40, false), 6.39550252992921e-37, 1e-10, "dt(40,50,ncp=1)");  // ref: exp(dtref(40,50,1))  (R: dt(40,50,ncp=1) = 0.0)
  check::rel(Dist::t(50, 37).d(3, false), 3.013022130548561e-226, 1e-10, "dt(3,50,ncp=37)");  // ref: exp(dtref(3,50,37))  (R: dt(3,50,ncp=37) = 2.1883685724362705e-254)
  check::close(Dist::chisq(0, 2).q(0.01, true, false), 0, 0, "qchisq(0.01,0,ncp=2)");  // exact: P(X = 0) = exp(-1) > 0.01, so the quantile is 0 (R: qchisq(0.01,0,ncp=2) gives 1.5192533514048113e-308)
  check::rel(Dist::chisq(3, 5).p(300, false, false), 7.969988036933611e-51, 1e-10, "pchisq(300,3,ncp=5,lower.tail=FALSE)");  // ref: pcref(300,3,5,FALSE)  (R: pchisq(300,3,ncp=5,lower.tail=FALSE) = 7.373250990247633e-51)
  check::rel(Dist::chisq(1, 1).p(300, true, true), -3.527300724008134e-60, 1e-10, "pchisq(300,1,ncp=1,log.p=TRUE)");  // ref: log1p(-pcref(300,1,1,FALSE))
  check::rel(Dist::f(1, 5, 0.5).p(0.1, true, false), 0.18511610715020363, 1e-10, "pf(0.1,1,5,ncp=0.5)");  // ref: pfref(0.1,1,5,0.5,TRUE)  (R: pf(0.1,1,5,ncp=0.5) = 0.18511610712367244)
  check::rel(Dist::f(1, 5, 0.5).p(2, true, false), 0.6908589470965399, 1e-10, "pf(2,1,5,ncp=0.5)");  // ref: pfref(2,1,5,0.5,TRUE)  (R: pf(2,1,5,ncp=0.5) = 0.6908589462737124)
  check::rel(Dist::f(1, 5, 0.5).p(2, false, false), 0.3091410529034601, 1e-10, "pf(2,1,5,ncp=0.5,lower.tail=FALSE)");  // ref: pfref(2,1,5,0.5,FALSE)  (R: pf(2,1,5,ncp=0.5,lower.tail=FALSE) = 0.30914105372628764)
  check::rel(Dist::f(1, 5, 0.5).p(10, false, false), 0.04813195618258438, 1e-10, "pf(10,1,5,ncp=0.5,lower.tail=FALSE)");  // ref: pfref(10,1,5,0.5,FALSE)  (R: pf(10,1,5,ncp=0.5,lower.tail=FALSE) = 0.04813195623933039)
  check::rel(Dist::f(1, 5, 0.5).p(10, false, true), -3.0338089527659724, 1e-10, "pf(10,1,5,ncp=0.5,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(10,1,5,0.5,FALSE))  (R: pf(10,1,5,ncp=0.5,lower.tail=FALSE,log.p=TRUE) = -3.033808951587005)
  check::rel(Dist::f(3, 10, 4).p(0.5, true, false), 0.0768267558807746, 1e-10, "pf(0.5,3,10,ncp=4)");  // ref: pfref(0.5,3,10,4,TRUE)  (R: pf(0.5,3,10,ncp=4) = 0.07682675580880441)
  check::rel(Dist::f(3, 10, 4).p(2, true, false), 0.4663642160480802, 1e-10, "pf(2,3,10,ncp=4)");  // ref: pfref(2,3,10,4,TRUE)  (R: pf(2,3,10,ncp=4) = 0.4663642159696639)
  check::rel(Dist::f(3, 10, 4).p(2, false, false), 0.5336357839519199, 1e-10, "pf(2,3,10,ncp=4,lower.tail=FALSE)");  // ref: pfref(2,3,10,4,FALSE)  (R: pf(2,3,10,ncp=4,lower.tail=FALSE) = 0.5336357840303362)
  check::rel(Dist::f(3, 10, 4).p(8, true, false), 0.9508059738873698, 1e-10, "pf(8,3,10,ncp=4)");  // ref: pfref(8,3,10,4,TRUE)  (R: pf(8,3,10,ncp=4) = 0.9508059729996408)
  check::rel(Dist::f(3, 10, 4).p(8, false, false), 0.049194026112630386, 1e-10, "pf(8,3,10,ncp=4,lower.tail=FALSE)");  // ref: pfref(8,3,10,4,FALSE)  (R: pf(8,3,10,ncp=4,lower.tail=FALSE) = 0.04919402700035924)
  check::rel(Dist::f(3, 10, 4).p(0.5, true, true), -2.5662023156358256, 1e-10, "pf(0.5,3,10,ncp=4,log.p=TRUE)");  // ref: log(pfref(0.5,3,10,4,TRUE))  (R: pf(0.5,3,10,ncp=4,log.p=TRUE) = -2.566202316572611)
  check::rel(Dist::f(3, 10, 4).p(8, false, true), -3.0119830833309296, 1e-10, "pf(8,3,10,ncp=4,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(8,3,10,4,FALSE))  (R: pf(8,3,10,ncp=4,lower.tail=FALSE,log.p=TRUE) = -3.011983065285469)
  check::rel(Dist::f(10, 30, 37).p(2, true, false), 0.016771898084784964, 1e-10, "pf(2,10,30,ncp=37)");  // ref: pfref(2,10,30,37,TRUE)  (R: pf(2,10,30,ncp=37) = 0.016771897696511714)
  check::rel(Dist::f(10, 30, 37).p(2, false, false), 0.9832281019152151, 1e-10, "pf(2,10,30,ncp=37,lower.tail=FALSE)");  // ref: pfref(2,10,30,37,FALSE)  (R: pf(2,10,30,ncp=37,lower.tail=FALSE) = 0.9832281023034883)
  check::rel(Dist::f(10, 30, 37).p(4, true, false), 0.3381391273745535, 1e-10, "pf(4,10,30,ncp=37)");  // ref: pfref(4,10,30,37,TRUE)  (R: pf(4,10,30,ncp=37) = 0.33813912690564923)
  check::rel(Dist::f(10, 30, 37).p(4, false, false), 0.6618608726254467, 1e-10, "pf(4,10,30,ncp=37,lower.tail=FALSE)");  // ref: pfref(4,10,30,37,FALSE)  (R: pf(4,10,30,ncp=37,lower.tail=FALSE) = 0.6618608730943507)
  check::rel(Dist::f(10, 30, 37).p(10, true, false), 0.9766745635669052, 1e-10, "pf(10,10,30,ncp=37)");  // ref: pfref(10,10,30,37,TRUE)  (R: pf(10,10,30,ncp=37) = 0.9766745632015471)
  check::rel(Dist::f(10, 30, 37).p(10, false, false), 0.0233254364330949, 1e-10, "pf(10,10,30,ncp=37,lower.tail=FALSE)");  // ref: pfref(10,10,30,37,FALSE)  (R: pf(10,10,30,ncp=37,lower.tail=FALSE) = 0.023325436798452847)
  check::rel(Dist::f(10, 30, 37).p(2, true, true), -4.088050526201446, 1e-10, "pf(2,10,30,ncp=37,log.p=TRUE)");  // ref: log(pfref(2,10,30,37,TRUE))  (R: pf(2,10,30,ncp=37,log.p=TRUE) = -4.0880505493516734)
  check::rel(Dist::f(10, 30, 37).p(10, false, true), -3.7582108214657963, 1e-10, "pf(10,10,30,ncp=37,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(10,10,30,37,FALSE))  (R: pf(10,10,30,ncp=37,lower.tail=FALSE,log.p=TRUE) = -3.758210805802298)
  check::rel(Dist::f(0.5, 2, 1).p(1, true, false), 0.4482700320176827, 1e-10, "pf(1,0.5,2,ncp=1)");  // ref: pfref(1,0.5,2,1,TRUE)  (R: pf(1,0.5,2,ncp=1) = 0.44827003144618377)
  check::rel(Dist::f(0.5, 2, 1).p(1, false, false), 0.5517299679823174, 1e-10, "pf(1,0.5,2,ncp=1,lower.tail=FALSE)");  // ref: pfref(1,0.5,2,1,FALSE)  (R: pf(1,0.5,2,ncp=1,lower.tail=FALSE) = 0.5517299685538163)
  check::rel(Dist::f(0.5, 2, 1).p(20, true, false), 0.8790498044654728, 1e-10, "pf(20,0.5,2,ncp=1)");  // ref: pfref(20,0.5,2,1,TRUE)  (R: pf(20,0.5,2,ncp=1) = 0.8790498038348014)
  check::rel(Dist::f(0.5, 2, 1).p(20, false, false), 0.1209501955345273, 1e-10, "pf(20,0.5,2,ncp=1,lower.tail=FALSE)");  // ref: pfref(20,0.5,2,1,FALSE)  (R: pf(20,0.5,2,ncp=1,lower.tail=FALSE) = 0.12095019616519861)
  check::rel(Dist::f(0.5, 2, 1).p(20, false, true), -2.112376425271539, 1e-10, "pf(20,0.5,2,ncp=1,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(20,0.5,2,1,FALSE))  (R: pf(20,0.5,2,ncp=1,lower.tail=FALSE,log.p=TRUE) = -2.1123764200572333)
  check::rel(Dist::f(4, 0.7, 2).p(1, true, false), 0.25846531380770255, 1e-10, "pf(1,4,0.7,ncp=2)");  // ref: pfref(1,4,0.7,2,TRUE)  (R: pf(1,4,0.7,ncp=2) = 0.2584653135533275)
  check::rel(Dist::f(4, 0.7, 2).p(1, false, false), 0.7415346861922975, 1e-10, "pf(1,4,0.7,ncp=2,lower.tail=FALSE)");  // ref: pfref(1,4,0.7,2,FALSE)  (R: pf(1,4,0.7,ncp=2,lower.tail=FALSE) = 0.7415346864466725)
  check::rel(Dist::f(4, 0.7, 2).p(50, true, false), 0.784608132780013, 1e-10, "pf(50,4,0.7,ncp=2)");  // ref: pfref(50,4,0.7,2,TRUE)  (R: pf(50,4,0.7,ncp=2) = 0.7846081322670034)
  check::rel(Dist::f(4, 0.7, 2).p(50, false, false), 0.21539186721998704, 1e-10, "pf(50,4,0.7,ncp=2,lower.tail=FALSE)");  // ref: pfref(50,4,0.7,2,FALSE)  (R: pf(50,4,0.7,ncp=2,lower.tail=FALSE) = 0.21539186773299657)
  check::rel(Dist::f(4, 0.7, 2).p(1, true, true), -1.3529937765990177, 1e-10, "pf(1,4,0.7,ncp=2,log.p=TRUE)");  // ref: log(pfref(1,4,0.7,2,TRUE))  (R: pf(1,4,0.7,ncp=2,log.p=TRUE) = -1.3529937775831924)
  check::rel(Dist::f(4, 0.7, 2).p(50, false, true), -1.5352962716115204, 1e-10, "pf(50,4,0.7,ncp=2,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(50,4,0.7,2,FALSE))  (R: pf(50,4,0.7,ncp=2,lower.tail=FALSE,log.p=TRUE) = -1.5352962692297707)
  check::rel(Dist::f(20, 100, 10).p(1, true, false), 0.1366539835409977, 1e-10, "pf(1,20,100,ncp=10)");  // ref: pfref(1,20,100,10,TRUE)  (R: pf(1,20,100,ncp=10) = 0.13665398342723067)
  check::rel(Dist::f(20, 100, 10).p(1, false, false), 0.8633460164590024, 1e-10, "pf(1,20,100,ncp=10,lower.tail=FALSE)");  // ref: pfref(1,20,100,10,FALSE)  (R: pf(1,20,100,ncp=10,lower.tail=FALSE) = 0.8633460165727693)
  check::rel(Dist::f(20, 100, 10).p(1.5, true, false), 0.5261319329527165, 1e-10, "pf(1.5,20,100,ncp=10)");  // ref: pfref(1.5,20,100,10,TRUE)  (R: pf(1.5,20,100,ncp=10) = 0.5261319323127243)
  check::rel(Dist::f(20, 100, 10).p(1.5, false, false), 0.47386806704728357, 1e-10, "pf(1.5,20,100,ncp=10,lower.tail=FALSE)");  // ref: pfref(1.5,20,100,10,FALSE)  (R: pf(1.5,20,100,ncp=10,lower.tail=FALSE) = 0.47386806768727574)
  check::rel(Dist::f(20, 100, 10).p(3, true, false), 0.9897909829990946, 1e-10, "pf(3,20,100,ncp=10)");  // ref: pfref(3,20,100,10,TRUE)  (R: pf(3,20,100,ncp=10) = 0.9897909827712881)
  check::rel(Dist::f(20, 100, 10).p(3, false, false), 0.010209017000905454, 1e-10, "pf(3,20,100,ncp=10,lower.tail=FALSE)");  // ref: pfref(3,20,100,10,FALSE)  (R: pf(3,20,100,ncp=10,lower.tail=FALSE) = 0.010209017228711962)
  check::rel(Dist::f(20, 100, 10).p(1, true, true), -1.990303215616945, 1e-10, "pf(1,20,100,ncp=10,log.p=TRUE)");  // ref: log(pfref(1,20,100,10,TRUE))  (R: pf(1,20,100,ncp=10,log.p=TRUE) = -1.9903032164494636)
  check::rel(Dist::f(20, 100, 10).p(3, false, true), -4.584483929510578, 1e-10, "pf(3,20,100,ncp=10,lower.tail=FALSE,log.p=TRUE)");  // ref: log(pfref(3,20,100,10,FALSE))  (R: pf(3,20,100,ncp=10,lower.tail=FALSE,log.p=TRUE) = -4.584483907196334)
  check::rel(Dist::f(3, 1000000.0, 1).p(2, true, false), 0.7876243783284491, 1e-10, "pf(2,3,1000000.0,ncp=1)");  // ref: pfref(2,3,1000000.0,1,TRUE)  (R: pf(2,3,1000000.0,ncp=1) = 0.7876243782048199)
  check::rel(Dist::tukey(5, 60, 1).p(4.4, false, false), 0.02298078074364751, 1e-10, "ptukey(4.4,5,60,1,lower.tail=FALSE)");  // ref: (1-ptkref(4.4,5,60,1))  (R: ptukey(4.4,5,60,1,lower.tail=FALSE) = 0.022980780771525433)
  check::rel(Dist::tukey(20, 10000.0, 1).p(5, false, false), 0.05152436860337073, 1e-10, "ptukey(5,20,10000.0,1,lower.tail=FALSE)");  // ref: (1-ptkref(5,20,10000.0,1))  (R: ptukey(5,20,10000.0,1,lower.tail=FALSE) = 0.05152436822921147)
  check::rel(Dist::tukey(50, 10, 3).p(10, true, false), 0.9865434560158026, 1e-10, "ptukey(10,50,10,3)");  // ref: ptkref(10,50,10,3)  (R: ptukey(10,50,10,3) = 0.9865501813627222)
  check::rel(Dist::tukey(50, 10, 3).p(10, false, false), 0.013456543984197378, 1e-10, "ptukey(10,50,10,3,lower.tail=FALSE)");  // ref: (1-ptkref(10,50,10,3))  (R: ptukey(10,50,10,3,lower.tail=FALSE) = 0.013449818637277788)
  check::rel(Dist::tukey(2, 2, 1).p(1, true, false), 0.44721359549995787, 1e-10, "ptukey(1,2,2,1)");  // ref: ptkref(1,2,2,1)  (R: ptukey(1,2,2,1) = 0.4472265936580529)
  check::rel(Dist::tukey(2, 2, 1).p(10, true, false), 0.9805806756909199, 1e-10, "ptukey(10,2,2,1)");  // ref: ptkref(10,2,2,1)  (R: ptukey(10,2,2,1) = 0.9807140305343182)
  check::rel(Dist::tukey(2, 2, 1).p(10, false, false), 0.01941932430908011, 1e-10, "ptukey(10,2,2,1,lower.tail=FALSE)");  // ref: (1-ptkref(10,2,2,1))  (R: ptukey(10,2,2,1,lower.tail=FALSE) = 0.01928596946568184)
}

std::string str(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.17g", v);
  return b;
}

// q then p returns the probability (on the scale and in the tail it was given).
void round_trip(const Dist& D, double p, bool lower, bool lg, double tol, const std::string& what) {
  const double x = D.q(p, lower, lg);
  check::rel(D.p(x, lower, lg), p, tol, "round trip " + what + " p=" + str(p) + " x=" + str(x));
}

void round_trips() {
  struct Case { Dist d; const char* name; };
  const std::vector<Case> cs = {
      {Dist::t(0.5, 4), "t(0.5, 4)"}, {Dist::t(3, -2), "t(3, -2)"}, {Dist::t(10, 2), "t(10, 2)"},
      {Dist::t(50, 37), "t(50, 37)"}, {Dist::t(1e6, 1.5), "t(1e6, 1.5)"},
      {Dist::chisq(0.5, 0.1), "chisq(0.5, 0.1)"}, {Dist::chisq(3, 5), "chisq(3, 5)"}, {Dist::chisq(10, 37), "chisq(10, 37)"},
      {Dist::chisq(100, 200), "chisq(100, 200)"},
      {Dist::f(1, 5, 0.5), "f(1, 5, 0.5)"}, {Dist::f(3, 10, 4), "f(3, 10, 4)"}, {Dist::f(10, 30, 37), "f(10, 30, 37)"},
      {Dist::f(4, 0.7, 2), "f(4, 0.7, 2)"}, {Dist::f(5, kInf, 2), "f(5, Inf, 2)"}};
  for (const auto& c : cs) {
    for (double p : {1e-12, 1e-6, 0.05, 0.5, 0.95}) {
      round_trip(c.d, p, true, false, 1e-11, c.name);
      round_trip(c.d, p, false, false, 1e-11, c.name);
    }
    round_trip(c.d, -30.0, true, true, 1e-11, std::string(c.name) + " log");
    round_trip(c.d, -1e-3, false, true, 1e-10, std::string(c.name) + " log upper");
  }
  // Studentized range (slower: a double integral per evaluation).
  round_trip(Dist::tukey(3, 10), 0.95, true, false, 1e-10, "tukey(3, 10)");
  round_trip(Dist::tukey(5, kInf), 1e-3, false, false, 1e-9, "tukey(5, Inf)");
  round_trip(Dist::tukey(4, 20, 2), 0.1, true, false, 1e-10, "tukey(4, 20, 2)");
}

void identities() {
  // Both tails are computed separately; they must add to 1.
  for (double x : {-3.0, -0.5, 0.7, 2.0, 9.0}) {
    const Dist d = Dist::t(4, 1.3);
    check::close(d.p(x, true, false) + d.p(x, false, false), 1, 1e-14, "pt lower + upper at " + str(x));
  }
  for (double x : {0.3, 4.0, 12.0}) {
    const Dist c = Dist::chisq(3, 5), f = Dist::f(3, 10, 4);
    check::close(c.p(x, true, false) + c.p(x, false, false), 1, 1e-14, "pchisq lower + upper at " + str(x));
    check::close(f.p(x, true, false) + f.p(x, false, false), 1, 1e-14, "pf lower + upper at " + str(x));
  }
  // Two means: the range is |Z1 - Z2|, so Q = sqrt(2) |t_df| and ptukey(q, 2, df) = 2 pt(q / sqrt(2), df) - 1.
  for (double df : {2.0, 5.0, 40.0}) {
    for (double q : {0.5, 2.0, 7.0}) {
      const double want = 2 * Dist::t(df).p(q / std::sqrt(2.0)) - 1;
      check::rel(Dist::tukey(2, df).p(q), want, 1e-11, "ptukey(" + str(q) + ", 2, " + str(df) + ") = 2 pt - 1");
    }
  }
  // Noncentral F with ncp is chisq(df1, ncp) / df1 when df2 = Inf; t(Inf, ncp) is N(ncp, 1).
  check::rel(Dist::f(5, kInf, 2).p(1.3), Dist::chisq(5, 2).p(6.5), 1e-14, "pf(1.3, 5, Inf, 2) = pchisq(6.5, 5, 2)");
  check::rel(Dist::t(kInf, 0.4).q(0.8), 0.4 + Dist::norm(0, 1).q(0.8), 1e-14, "qt(0.8, Inf, 0.4)");
  // F with df1 = 1 and the square of a noncentral t: P(F <= x) = P(-sqrt(x) <= T <= sqrt(x)) with ncp^2.
  {
    const Dist t = Dist::t(7, 1.5), f = Dist::f(1, 7, 2.25);
    const double x = 3.0, r = std::sqrt(x);
    check::rel(f.p(x), t.p(r) - t.p(-r), 1e-12, "pf(3, 1, 7, 2.25) = pt(sqrt 3) - pt(-sqrt 3) with ncp 1.5");
  }
  // Tiny ncp approaches the central distribution.
  check::rel(Dist::chisq(3, 1e-12).p(2.0), Dist::chisq(3).p(2.0), 1e-11, "pchisq ncp 1e-12");
  check::rel(Dist::t(6, 1e-12).p(1.1), Dist::t(6).p(1.1), 1e-11, "pt ncp 1e-12");
}

void moments() {
  // mean = ncp sqrt(df/2) Gamma((df-1)/2) / Gamma(df/2), var = df (1 + ncp^2) / (df - 2) - mean^2.
  check::rel(Dist::t(10, 2).mean(), 2.1674446158770424, 1e-11, "nct mean");  // R: integrate(function(x) x*dt(x,10,ncp=2),-Inf,Inf,rel.tol=1e-12)
  check::rel(Dist::t(10, 2).var(), 1.5521838370894185, 1e-10, "nct var");    // R: integrate of x^2 dt(x,10,ncp=2) - mean^2
  check::rel(Dist::t(3, -1.5).mean(), -2.0729648968292107, 1e-11, "nct mean df 3");  // R: integrate(function(x) x*dt(x,3,ncp=-1.5),-Inf,Inf,rel.tol=1e-12)
  check::ok(std::isinf(Dist::t(2, 1).var()) && std::isfinite(Dist::t(2, 1).mean()), "nct df 2: finite mean, Inf var");
  check::ok(std::isnan(Dist::t(1, 1).mean()) && std::isnan(Dist::t(1, 1).var()), "nct df 1: NaN moments");
  check::close(Dist::t(kInf, 1.5).mean(), 1.5, 0, "nct df Inf mean");
  check::close(Dist::t(kInf, 1.5).var(), 1, 0, "nct df Inf var");
  check::close(Dist::chisq(3, 5).mean(), 8, 0, "ncchisq mean");  // R: integrate(function(x) x*dchisq(x,3,ncp=5),0,Inf) = 8
  check::close(Dist::chisq(3, 5).var(), 26, 0, "ncchisq var");   // 2 (df + 2 ncp)
  check::rel(Dist::f(3, 10, 4).mean(), 2.9166666666666665, 1e-14, "ncf mean");  // R: integrate(function(x) x*df(x,3,10,ncp=4),0,Inf,rel.tol=1e-12)
  check::rel(Dist::f(3, 10, 4).var(), 7.9282407407407405, 1e-12, "ncf var");    // R: integrate of x^2 df(x,3,10,ncp=4) - mean^2
  check::rel(Dist::f(5, kInf, 2).mean(), 1.4, 1e-15, "ncf df2 Inf mean");
  check::rel(Dist::f(5, kInf, 2).var(), 0.72, 1e-15, "ncf df2 Inf var");
  check::ok(std::isinf(Dist::f(3, 2, 1).mean()) && std::isnan(Dist::f(3, 2, 1).var()), "ncf df2 2: Inf mean, NaN var");
  check::ok(std::isinf(Dist::f(3, 4, 1).var()), "ncf df2 4: Inf var");
  check::ok(std::isnan(Dist::tukey(3, 10).mean()) && std::isnan(Dist::tukey(3, 10).var()), "tukey moments NaN");
}

void edges_and_api() {
  // ncp = 0 gives the central kinds.
  check::ok(Dist::t(5, 0).kind() == DistKind::T, "t(5, 0) is central");
  check::ok(Dist::chisq(5, 0).kind() == DistKind::Chisq, "chisq(5, 0) is central");
  check::ok(Dist::f(5, 3, 0).kind() == DistKind::F, "f(5, 3, 0) is central");
  check::ok(Dist::t(5, 1).kind() == DistKind::NcT && Dist::t(5, 1).name() == "t", "NcT kind and name");
  check::ok(Dist::chisq(5, 1).kind() == DistKind::NcChisq && Dist::chisq(5, 1).name() == "chisq", "NcChisq kind and name");
  check::ok(Dist::f(5, 3, 1).kind() == DistKind::NcF && Dist::f(5, 3, 1).name() == "f", "NcF kind and name");
  check::ok(Dist::tukey(3, 10).kind() == DistKind::Tukey && Dist::tukey(3, 10).name() == "tukey", "Tukey kind and name");
  check::ok(!Dist::t(5, 1).discrete() && !Dist::tukey(3, 10).discrete(), "continuous");
  check::close(Dist::tukey(3, 10).param(2), 1, 0, "tukey nranges default 1");
  // by_name
  check::ok(Dist::by_name("t", 5).kind() == DistKind::T, "by_name t central");
  check::ok(Dist::by_name("t", 5, 1).kind() == DistKind::NcT, "by_name t ncp");
  check::ok(Dist::by_name("chisq", 5, 1).kind() == DistKind::NcChisq, "by_name chisq ncp");
  check::ok(Dist::by_name("f", 5, 3).kind() == DistKind::F, "by_name f central");
  check::ok(Dist::by_name("f", 5, 3, 1).kind() == DistKind::NcF, "by_name f ncp");
  check::ok(Dist::by_name("tukey", 4, 12).kind() == DistKind::Tukey && Dist::by_name("tukey", 4, 12).param(2) == 1, "by_name tukey");
  check::close(Dist::by_name("tukey", 4, 12, 2).param(2), 2, 0, "by_name tukey nranges");
  check::rel(Dist::by_name("t", 5, 1).p(1.7), Dist::t(5, 1).p(1.7), 0, "by_name t same values");
  check::throws([] { Dist::by_name("tukey", 4); }, "by_name tukey needs df");
  // Parameter errors (R returns NaN with a warning; the constructors throw).
  check::throws([] { Dist::t(0, 1); }, "t df 0");
  check::throws([] { Dist::t(5, kNaN); }, "t ncp NaN");
  check::throws([] { Dist::t(5, kInf); }, "t ncp Inf");
  check::throws([] { Dist::chisq(3, -1); }, "chisq ncp < 0");
  check::throws([] { Dist::chisq(3, kInf); }, "chisq ncp Inf");
  check::throws([] { Dist::chisq(-1, 1); }, "chisq df < 0");
  check::throws([] { Dist::f(3, 4, -1); }, "f ncp < 0");
  check::throws([] { Dist::f(kInf, 4, 1); }, "f df1 Inf with ncp");  // R: pf(2, Inf, 5, ncp = 1) is NaN
  check::throws([] { Dist::tukey(1, 10); }, "tukey nmeans 1");   // R: ptukey(3, 1, 10) is NaN
  check::throws([] { Dist::tukey(3, 1.5); }, "tukey df 1.5");     // R: ptukey(3, 2, 1.5) is NaN
  check::throws([] { Dist::tukey(3, 10, 0.5); }, "tukey nranges 0.5");
  check::throws([] { Dist::tukey(3, 10).d(1.0); }, "dtukey throws");
  {
    std::string msg;
    try { (void)Dist::tukey(3, 10).d(1.0); } catch (const rcli::Error& e) { msg = e.what(); }
    check::ok(msg == "dtukey is not available (R has no studentized range density)", "dtukey message: " + msg);
  }
  // Non-integer nmeans and nranges are allowed, as in R.
  check::rel(Dist::tukey(2.5, 10).p(3), 0.90326074214023211, 1e-8, "ptukey(3,2.5,10)");  // R: ptukey(3,2.5,10)
  check::rel(Dist::tukey(3, 10, 1.5).p(3), 0.81043303445291104, 1e-8, "ptukey(3,3,10,1.5)");  // R: ptukey(3,3,10,1.5)
  // NaN and out-of-range arguments.
  check::ok(std::isnan(Dist::t(5, 1).p(kNaN)) && std::isnan(Dist::chisq(3, 1).d(kNaN)) && std::isnan(Dist::f(3, 4, 1).q(kNaN)), "NaN in, NaN out");
  check::ok(std::isnan(Dist::t(5, 1).q(1.5)) && std::isnan(Dist::chisq(3, 1).q(-0.1)) && std::isnan(Dist::tukey(3, 10).q(0.1, true, true)), "q outside [0, 1] is NaN");
  // Support boundaries.
  check::close(Dist::t(5, 1).q(0), -kInf, 0, "qt(0, 5, ncp=1)");
  check::close(Dist::t(5, 1).q(1), kInf, 0, "qt(1, 5, ncp=1)");
  check::close(Dist::chisq(3, 1).q(0), 0, 0, "qchisq(0, 3, ncp=1)");
  check::close(Dist::chisq(3, 1).q(0, false), kInf, 0, "qchisq(0, 3, ncp=1, lower.tail=FALSE)");
  check::close(Dist::f(3, 4, 1).q(1), kInf, 0, "qf(1, 3, 4, ncp=1)");
  check::close(Dist::tukey(3, 10).q(0), 0, 0, "qtukey(0, 3, 10)");
  check::close(Dist::tukey(3, 10).q(1), kInf, 0, "qtukey(1, 3, 10)");
  check::close(Dist::tukey(3, 10).p(0), 0, 0, "ptukey(0, 3, 10)");
  check::close(Dist::tukey(3, 10).p(-1, false), 1, 0, "ptukey(-1, 3, 10, lower.tail=FALSE)");
  check::close(Dist::tukey(3, 10).p(kInf, false, true), -kInf, 0, "ptukey(Inf, 3, 10, lower.tail=FALSE, log.p=TRUE)");
  check::close(Dist::t(5, 1).p(kInf), 1, 0, "pt(Inf, 5, ncp=1)");
  check::close(Dist::t(5, 1).p(-kInf), 0, 0, "pt(-Inf, 5, ncp=1)");
  check::close(Dist::t(5, 1).d(kInf), 0, 0, "dt(Inf, 5, ncp=1)");
  check::close(Dist::chisq(3, 1).p(-1), 0, 0, "pchisq(-1, 3, ncp=1)");
  check::close(Dist::chisq(3, 1).d(-1), 0, 0, "dchisq(-1, 3, ncp=1)");
  check::close(Dist::f(3, 4, 1).d(-1, true), -kInf, 0, "df(-1, 3, 4, ncp=1, log=TRUE)");
  check::close(Dist::chisq(1, 3).d(0), kInf, 0, "dchisq(0, 1, ncp=3)");  // R: dchisq(0,1,ncp=3) = Inf
  check::close(Dist::f(1, 5, 1).d(0), kInf, 0, "df(0, 1, 5, ncp=1)");    // R: df(0,1,5,ncp=1) = Inf
  check::close(Dist::chisq(0, 2).d(0), kInf, 0, "dchisq(0, 0, ncp=2)");  // R: dchisq(0,0,ncp=2) = Inf
  // df = 0: point mass exp(-ncp/2) at 0.
  check::rel(Dist::chisq(0, 2).p(0), std::exp(-1.0), 1e-15, "pchisq(0, 0, ncp=2) = exp(-1)");
  check::close(Dist::chisq(0, 2).q(0.3), 0, 0, "qchisq(0.3, 0, ncp=2) inside the point mass");
  check::ok(Dist::chisq(0, 2).q(0.37) > 0, "qchisq(0.37, 0, ncp=2) beyond the point mass");
  // Far tails stay finite and positive on the log scale where R underflows to 0 / -Inf.
  check::rel(Dist::t(10, 37).p(-1, true, true), Dist::t(10, -37).p(1, false, true), 1e-14, "pt reflection log");
  check::ok(Dist::t(10, 37).p(-1, true, true) > -800 && Dist::t(10, 37).p(-1, true, true) < -600, "pt(-1, 10, ncp=37, log.p=TRUE) finite (R: -Inf)");
}

}  // namespace

int main() {
  r_oracle_cases();
  reference_cases();
  round_trips();
  identities();
  moments();
  edges_and_api();
  return check::report("test_ncdist");
}
