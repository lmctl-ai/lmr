// Descriptive statistics. NaN inputs are errors unless a function says otherwise (callers drop NA first).
#pragma once
#include <string>
#include "rcli/common.hpp"

namespace rcli::summary {

double sum(const Vec& x);               // compensated (Neumaier) summation
double mean(const Vec& x);              // two-pass mean as in R (mean then mean of residuals correction)
double var(const Vec& x);               // n - 1 denominator; NaN for n < 2
double sd(const Vec& x);
double median(const Vec& x);
double mad(const Vec& x, double constant = 1.4826);
double cov(const Vec& x, const Vec& y);
// method "pearson", "spearman" or "kendall" (Kendall's tau-b with ties, as R cor()).
double cor(const Vec& x, const Vec& y, const std::string& method = "pearson");
Vec rank(const Vec& x);                 // average ranks for ties, 1-based like R rank()
// Hyndman & Fan sample quantile types 1..9 (R quantile(type=)); type 7 is R's default.
double quantile(const Vec& x, double prob, int type = 7);
Vec quantiles(const Vec& x, const Vec& probs, int type = 7);

}  // namespace rcli::summary
