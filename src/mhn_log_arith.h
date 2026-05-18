// Log-space arithmetic helpers used by the Psi series, CDF series, and
// any other code that accumulates products of terms whose logarithms span
// a wide dynamic range.

#ifndef MHN_LOG_ARITH_H
#define MHN_LOG_ARITH_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mhn {

// log(sum(exp(log_x))) using the max-subtraction trick: factor out the
// largest log_x to keep every exp() near unity, avoiding overflow when
// the terms span more than ~700 in log-magnitude.
inline double log_sum_exp(const std::vector<double>& log_x) {
  if (log_x.empty()) return -std::numeric_limits<double>::infinity();
  double m = -std::numeric_limits<double>::infinity();
  for (double v : log_x) if (v > m) m = v;
  if (!std::isfinite(m)) return -std::numeric_limits<double>::infinity();
  double s = 0.0;
  for (double v : log_x) s += std::exp(v - m);
  return m + std::log(s);
}

// log(exp(a) + exp(b)) with numerical stability.  Pairs the larger value
// with std::log1p(exp(smaller - larger)) so the exp() argument is always
// non-positive (no overflow) and small differences are preserved by
// log1p rather than collapsing to log(1.0) = 0.
inline double log_add_exp(double a, double b) {
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

}  // namespace mhn

#endif  // MHN_LOG_ARITH_H
