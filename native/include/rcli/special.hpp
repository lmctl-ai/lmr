// Special functions. Accuracy target: 1e-13 relative where R's own functions reach it, 1e-10 for Bessel functions.
#pragma once
#include "rcli/common.hpp"

namespace rcli::special {

double gamma(double x);                 // R gamma(); NaN at non-positive integers
double lgamma(double x);                // log|gamma(x)|, R lgamma()
double beta(double a, double b);
double lbeta(double a, double b);
double digamma(double x);
double trigamma(double x);
double choose(double n, double k);      // R choose(): k rounded to an integer, n real (generalised binomial)
double lchoose(double n, double k);     // log|choose(n, k)|
double factorial(double x);             // gamma(x + 1)
double lfactorial(double x);

double bessel_j(double x, double nu);   // R besselJ(x, nu), x >= 0
double bessel_y(double x, double nu);   // R besselY(x, nu), x >= 0
double bessel_i(double x, double nu, bool expon_scaled = false);  // R besselI
double bessel_k(double x, double nu, bool expon_scaled = false);  // R besselK

double erf(double x);
double erfc(double x);
double log1pmx(double x);               // log(1 + x) - x, accurate near 0

// Regularized incomplete gamma P(a, x) and Q(a, x) = 1 - P, with log variants.
double gamma_p(double a, double x);
double gamma_q(double a, double x);
// Regularized incomplete beta I_x(a, b) and its complement 1 - I_x(a, b).
double beta_inc(double a, double b, double x);
double beta_inc_c(double a, double b, double x);

// --- Backward-compatible additions (dist module kernels; usable by other modules) ---
double lgamma1p(double a);              // log(gamma(1 + a)), accurate for small |a|
double stirlerr(double n);              // lgamma(n + 1) - (n + 1/2) log(n) + n - log(sqrt(2 pi)), n > 0
double bd0(double x, double np);        // x log(x / np) + np - x, accurate when x is close to np
double dpois_raw(double x, double lambda, bool give_log);                  // lambda^x e^-lambda / gamma(x + 1), real x >= 0
double dbinom_raw(double x, double n, double p, double q, bool give_log);  // q = 1 - p supplied for accuracy
// Regularized incomplete gamma with tail and log selection: lower ? P(a, x) : Q(a, x), log if log_p.
double gamma_inc(double a, double x, bool lower, bool log_p);
// Regularized incomplete beta with y = 1 - x supplied separately: lower ? I_x(a, b) : 1 - I_x(a, b), log if log_p.
double beta_inc_tail(double a, double b, double x, double y, bool lower, bool log_p);
// Standard normal lower and upper tail probabilities of z, both accurate (log scale if log_p).
void pnorm_both(double z, double* lower, double* upper, bool log_p);

}  // namespace rcli::special
