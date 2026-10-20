// Mode of the MHN distribution.  The cases follow Sun et al. (2023),
// Lemma 3b-d, except at alpha = 1, where the mode of the truncated normal of
// Lemma 6b is returned.  Lemma 3d says that no mode exists for gamma < 0 and
// alpha <= 1; at alpha = 1 the density is bounded and decreases from its
// supremum at the boundary x = 0, and that point, 0, is returned.

#include "mhn_check.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::export(.mhn_mode_cpp, rng = false)]]
double mhn_mode_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);

  if (alpha > 1.0) {
    // Sun et al. (2023) Lemma 3b: unique interior mode.
    return mhn::positive_root(beta, gamma, alpha - 1.0);
  }
  if (alpha == 1.0) {
    // Sun et al. (2023) Lemma 6b: alpha = 1 truncated-normal mode.
    //
    // Only alpha = 1 itself takes this branch.  The Lemma 3b root above tends
    // to this value as alpha -> 1+, and so does the Lemma 3c one below for
    // gamma > 0, so no tolerance band around 1 is needed, and one would be
    // wrong on both sides.  Near gamma = 0 the Lemma 3b mode grows like
    // sqrt((alpha - 1) / (2 beta)), so snapping alpha = 1 + 1e-8 to 1 returns 0
    // for a mode of 7.1e-5 at beta = 1; and for alpha just below 1 with
    // gamma < 0 (Lemma 3d) or gamma = 0 (where the derivative of log f is
    // negative) there is no mode at all.  2 beta overflows for beta above
    // half the largest double, where the mode is (gamma / 2) / beta instead;
    // the test is made on beta (see opaque_value in mhn_constants.h).
    return std::max(0.0, (beta <= 0.5 * mhn::kDoubleMax)
                             ? gamma / (2.0 * beta)
                             : 0.5 * gamma / beta);
  }
  // 0 < alpha < 1.  Here d/dx log f = -(2 beta x^2 - gamma x + 1 - alpha) / x.
  // A local maximum exists only for gamma > 0 and
  // alpha > 1 - gamma^2 / (8 beta) (Sun et al. 2023, Lemma 3c).  The lemma
  // prints >=, but at equality the derivative is
  // -(2 beta / x) (x - gamma / (4 beta))^2, negative on both sides of its
  // zero: the two stationary points have merged into an inflection point and
  // the density is strictly decreasing, so the strict inequality is used.
  // In every other case the density decreases from its unbounded value at the
  // origin (Lemma 3d, which lists gamma < 0 and gamma > 0; at gamma = 0 the
  // derivative is plainly negative) and there is no mode.
  const double m = mhn::kernel_interior_peak(alpha, beta, gamma);
  return (m > 0.0) ? m : NA_REAL;
}
