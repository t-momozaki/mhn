// Mode of the MHN distribution.  Cases follow Sun et al. (2023),
// Lemma 3b-d and Lemma 6b (alpha = 1, which is the truncated-normal
// boundary not directly covered by Lemma 3).  Mirrors R-side
// mhn/R/mhn_mode.R.

#include "mhn_check.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export(.mhn_mode_cpp)]]
double mhn_mode_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);

  const double eps = mhn::mhn_eps();

  if (alpha > 1.0 + eps) {
    // Sun et al. (2023) Lemma 3b: unique interior mode.
    return (gamma + std::sqrt(gamma * gamma + 8.0 * beta * (alpha - 1.0)))
           / (4.0 * beta);
  }
  if (std::fabs(alpha - 1.0) <= eps) {
    // Sun et al. (2023) Lemma 6b: alpha = 1 truncated-normal mode.
    return std::max(0.0, gamma / (2.0 * beta));
  }
  // alpha < 1 (Sun et al. 2023, Lemma 3c/3d).
  if (gamma > 0.0 && alpha >= 1.0 - gamma * gamma / (8.0 * beta) - eps) {
    const double disc = gamma * gamma + 8.0 * beta * (alpha - 1.0);
    if (disc >= 0.0) {
      return (gamma + std::sqrt(disc)) / (4.0 * beta);
    }
  }
  // Density is monotonically decreasing: no interior mode.
  return NA_REAL;
}
