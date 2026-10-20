// Parameter validation for MHN(alpha, beta, gamma).  The parameter space is
// that of Sun et al. (2023), Eq. (1): alpha > 0, beta > 0 and gamma
// real; the shape and the rate must also be finite.
//   * params_valid / invalid_result - element by element, for dmhn, pmhn,
//     qmhn and rmhn, which follow base R's distribution functions: an element
//     whose parameters lie outside the space is NaN (the R wrappers of
//     dmhn, pmhn and qmhn, and rmhn itself, warn of it), and one with an NA
//     or NaN parameter is NA or NaN, with no error for the call.
//   * check_params_scalar - throws on any invalid value.  Used by the moment
//     functions, mhn_mode and the internal sampler entry points and
//     diagnostic hooks in mhn_sun.cpp and mhn_rtdr.cpp, which take one
//     triple.

#ifndef MHN_CHECK_H
#define MHN_CHECK_H

#include <Rcpp.h>
#include <cmath>

namespace mhn {

// Whether (alpha, beta, gamma) lies in the parameter space.  ISNAN, R's own
// test, is made first and apart from the comparisons, so that an NA or NaN is
// refused by it rather than by how the comparisons below treat a NaN.
inline bool params_valid(double alpha, double beta, double gamma) {
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma)) return false;
  return alpha > 0.0 && beta > 0.0 && std::isfinite(alpha) &&
         std::isfinite(beta) && std::isfinite(gamma);
}

// The value of a d/p/q element whose parameters are not valid, as base R's
// distribution functions give it: NA if any input is NA, NaN otherwise.
// R_IsNA is tested on each, since NA_REAL is itself a NaN.
inline double invalid_result(double x, double alpha, double beta,
                             double gamma) {
  return (R_IsNA(x) || R_IsNA(alpha) || R_IsNA(beta) || R_IsNA(gamma))
             ? NA_REAL : R_NaN;
}

// Whether every element of alpha, beta and gamma is valid on its own, so that
// every recycled triple is.  The evaluators then skip the test per element.
bool all_params_valid(const Rcpp::NumericVector& alpha,
                      const Rcpp::NumericVector& beta,
                      const Rcpp::NumericVector& gamma);

// Throws Rcpp::exception unless (alpha, beta, gamma) is valid.  Its messages
// are those asserted in tests/testthat/test-input-validation.R.
void check_params_scalar(double alpha, double beta, double gamma);

}  // namespace mhn

#endif  // MHN_CHECK_H
