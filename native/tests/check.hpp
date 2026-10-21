// Minimal unit-test helpers for the native core. Each test file is its own executable that returns nonzero on failure.
#pragma once
#include <cmath>
#include <cstdio>
#include <string>
#include "rcli/common.hpp"

namespace check {
inline int failures = 0;
inline int passes = 0;

inline void ok(bool cond, const std::string& what) {
  if (cond) { ++passes; return; }
  ++failures;
  std::fprintf(stderr, "FAIL %s\n", what.c_str());
}
// |actual - expected| <= tol * max(|expected|, floor); NaN == NaN and equal infinities pass.
inline void close(double actual, double expected, double tol, const std::string& what, double floor = 1.0) {
  bool good;
  if (std::isnan(expected)) good = std::isnan(actual);
  else if (std::isinf(expected)) good = actual == expected;
  else good = std::fabs(actual - expected) <= tol * std::fmax(std::fabs(expected), floor);
  if (good) { ++passes; return; }
  ++failures;
  std::fprintf(stderr, "FAIL %s: actual %.17g expected %.17g (tol %g)\n", what.c_str(), actual, expected, tol);
}
// Relative closeness without an absolute floor: for tail probabilities that must be right in relative terms.
inline void rel(double actual, double expected, double tol, const std::string& what) {
  close(actual, expected, tol, what, 0.0);
}
template <class F>
void throws(F&& f, const std::string& what) {
  try { f(); } catch (const rcli::Error&) { ++passes; return; } catch (...) {}
  ++failures;
  std::fprintf(stderr, "FAIL %s: expected rcli::Error\n", what.c_str());
}
inline int report(const char* name) {
  std::fprintf(stderr, "%s: %d passed, %d failed\n", name, passes, failures);
  return failures ? 1 : 0;
}
}  // namespace check
