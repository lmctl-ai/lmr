// Shared conventions for the native numeric core (no R, no third-party libraries, C++17 standard library only).
// - Missing values are NaN; functions that cannot accept NaN say so and throw.
// - Invalid arguments throw rcli::Error with a message that names the function and the problem.
// - Indices exposed to scripts are 0-based unless a field says "1-based (R convention)".
#pragma once
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace rcli {

struct Error : std::runtime_error {
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};

using Vec = std::vector<double>;
using Fn1 = std::function<double(double)>;
using FnN = std::function<double(const Vec&)>;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kEps = std::numeric_limits<double>::epsilon();  // 2.220446e-16, R's .Machine$double.eps
inline const double kEps25 = std::pow(kEps, 0.25);              // 1.220703e-4, R's default tol for optimize/uniroot/integrate

// Dense row-major matrix.
struct Matrix {
  size_t rows = 0, cols = 0;
  Vec a;
  Matrix() = default;
  Matrix(size_t r, size_t c, double fill = 0.0) : rows(r), cols(c), a(r * c, fill) {}
  double& operator()(size_t i, size_t j) { return a[i * cols + j]; }
  double operator()(size_t i, size_t j) const { return a[i * cols + j]; }
};

}  // namespace rcli
