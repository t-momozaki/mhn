// Package-wide numerical constants.  Currently a single threshold used
// by every special-case dispatcher; centralising it here keeps the
// "is alpha exactly 1?" test consistent across dmhn / pmhn / qmhn / rmhn.

#ifndef MHN_CONSTANTS_H
#define MHN_CONSTANTS_H

#include <cmath>
#include <limits>

namespace mhn {

// Square root of machine epsilon, mirroring R's `.Machine$double.eps^0.5`.
// Used as the tolerance for special-case dispatch: parameters within this
// distance of the canonical boundary (gamma == 0, alpha == 1, alpha == 2)
// are routed to the closed-form branch instead of the general series.
inline double mhn_eps() {
  static const double v = std::sqrt(std::numeric_limits<double>::epsilon());
  return v;
}

}  // namespace mhn

#endif  // MHN_CONSTANTS_H
