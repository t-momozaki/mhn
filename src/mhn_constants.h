// Package-wide numerical constant sqrt(eps), the default working tolerance
// and the width of rmhn's special-case windows.

#ifndef MHN_CONSTANTS_H
#define MHN_CONSTANTS_H

#include <cmath>
#include <limits>

namespace mhn {

// sqrt(2), log(2) and pi, the doubles nearest them as <cmath>'s M_SQRT2,
// M_LN2 and M_PI give them where it defines those (they are not ISO C++),
// so that a header read before R's needs neither.
constexpr double kSqrt2 = 1.41421356237309504880;
constexpr double kLn2 = 0.693147180559945309417;
constexpr double kPi = 3.14159265358979323846;

// Square root of machine epsilon, mirroring R's `.Machine$double.eps^0.5`.
// It is the default tolerance of the series and quadratures, and the width
// of the special-case windows of rmhn (is_sqrt_gamma and is_truncated_normal
// in mhn_special_cases.cpp): a scale-free tilt |gamma| / sqrt(beta) below it
// is drawn as gamma = 0, and a shape within it of 1 as alpha = 1, the whole
// law replaced by the closed form (sqrt-Gamma or truncated normal).  dmhn,
// pmhn, qmhn and the normalizing-constant dispatcher take the closed forms
// only at gamma = 0 and alpha = 1 exactly (is_sqrt_gamma_exact).
inline double mhn_eps() {
  static const double v = std::sqrt(std::numeric_limits<double>::epsilon());
  return v;
}

// x itself, read through a volatile so that the compiler cannot relate it to
// the values it was formed from.  GCC 14 and 15 (14.2, 14.3 and 15.2 on
// x86_64, 14.3 and 15.1 on aarch64) took 2 beta > DBL_MAX to mean
// beta = +Inf, and replaced beta by +Inf on that branch, which a finite beta
// above DBL_MAX / 2 is not.  On aarch64 they did so where the test singled
// out +Inf (x > DBL_MAX, x == +Inf) of a constant times one variable, and not
// for std::isfinite, std::isinf, !(x <= DBL_MAX) or a value formed from two
// variables.  A comparison of a computed sum or product with the largest
// double is therefore written on opaque_value(x).  Where the test is of a
// power of two times a value, it is made on the value against the largest
// double divided by that power instead, which is exact.
inline double opaque_value(double x) {
  volatile double v = x;
  return v;
}

}  // namespace mhn

#endif  // MHN_CONSTANTS_H
