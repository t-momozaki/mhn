// Parameter validation for MHN(alpha, beta, gamma).  Three variants:
//   * scalar               - throws on any invalid value (alpha <= 0,
//                            beta <= 0, or non-finite gamma).
//   * vector               - element-wise; NA in any field is treated
//                            as invalid.  Used by dmhn / pmhn / qmhn.
//   * vector with NA pass  - element-wise but propagates NA into the
//                            caller's output.  Used by rmhn, where
//                            invalid rows yield NA samples rather than
//                            aborting the entire call.

#ifndef MHN_CHECK_H
#define MHN_CHECK_H

#include <Rcpp.h>

namespace mhn {

// C++-callable parameter validation.  Throws Rcpp::exception on failure
// with the same error messages as the R-side `.check_mhn_params` so
// existing regex-based expect_error() tests apply to both paths.
void check_params_scalar(double alpha, double beta, double gamma);

// Element-wise validation for the vectorized dmhn entry point.
// Each input must be non-empty; every element of alpha/beta must be > 0
// and every element of gamma must be finite.  Error messages match the
// scalar validator so existing error regex tests still apply.
void check_params_vector(const Rcpp::NumericVector& alpha,
                         const Rcpp::NumericVector& beta,
                         const Rcpp::NumericVector& gamma);

// Lenient variant for the rmhn entry point. NA elements are allowed and
// will propagate to NA_REAL in the caller's output. Non-NA elements of
// alpha/beta must be > 0; gamma is accepted in any finite or non-finite
// form (Inf/-Inf/NaN are passed through and the caller treats them as NA).
void check_params_vector_allow_na(const Rcpp::NumericVector& alpha,
                                  const Rcpp::NumericVector& beta,
                                  const Rcpp::NumericVector& gamma);

}  // namespace mhn

#endif  // MHN_CHECK_H
