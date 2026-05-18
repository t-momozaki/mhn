// Implementation of the parameter-validation helpers declared in
// mhn_check.h.  Error messages are kept identical across the scalar /
// vector / NA-aware paths so that the regex-based expect_error() tests
// in tests/testthat/test-input-validation.R match every branch.

#include "mhn_check.h"

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

// Replicates R's `is.numeric(x) && length(x) == 1L && !is.na(x)`.
// `expect_finite_for_gamma` toggles the "finite" requirement (gamma must be finite,
// while alpha/beta only need to be a non-NA scalar before the >0 check).
bool is_scalar_numeric(SEXP x, bool require_finite) {
  if (Rf_length(x) != 1) return false;
  switch (TYPEOF(x)) {
    case REALSXP: {
      double v = REAL(x)[0];
      if (R_IsNA(v) || R_IsNaN(v)) return false;
      if (require_finite && !R_finite(v)) return false;
      return true;
    }
    case INTSXP: {
      int v = INTEGER(x)[0];
      return v != NA_INTEGER;
    }
    case LGLSXP: {
      // R's is.numeric() returns FALSE for logical, mirror that.
      return false;
    }
    default:
      return false;
  }
}

double as_double_unchecked(SEXP x) {
  if (TYPEOF(x) == REALSXP) return REAL(x)[0];
  if (TYPEOF(x) == INTSXP) return static_cast<double>(INTEGER(x)[0]);
  return NA_REAL;
}

}  // namespace

namespace mhn {

void check_params_scalar(double alpha, double beta, double gamma) {
  if (!(alpha > 0.0)) Rcpp::stop("alpha must be positive");
  if (!(beta > 0.0))  Rcpp::stop("beta must be positive");
  if (!std::isfinite(gamma)) Rcpp::stop("gamma must be a finite numeric value");
}

void check_params_vector(const Rcpp::NumericVector& alpha,
                         const Rcpp::NumericVector& beta,
                         const Rcpp::NumericVector& gamma) {
  if (alpha.size() == 0) Rcpp::stop("alpha must be positive");
  if (beta.size()  == 0) Rcpp::stop("beta must be positive");
  if (gamma.size() == 0) Rcpp::stop("gamma must be a finite numeric value");
  for (R_xlen_t i = 0; i < alpha.size(); ++i) {
    const double a = alpha[i];
    if (Rcpp::NumericVector::is_na(a) || !(a > 0.0)) {
      Rcpp::stop("alpha must be positive");
    }
  }
  for (R_xlen_t i = 0; i < beta.size(); ++i) {
    const double b = beta[i];
    if (Rcpp::NumericVector::is_na(b) || !(b > 0.0)) {
      Rcpp::stop("beta must be positive");
    }
  }
  for (R_xlen_t i = 0; i < gamma.size(); ++i) {
    if (!std::isfinite(gamma[i])) {
      Rcpp::stop("gamma must be a finite numeric value");
    }
  }
}

void check_params_vector_allow_na(const Rcpp::NumericVector& alpha,
                                  const Rcpp::NumericVector& beta,
                                  const Rcpp::NumericVector& gamma) {
  if (alpha.size() == 0) Rcpp::stop("alpha must be positive");
  if (beta.size()  == 0) Rcpp::stop("beta must be positive");
  if (gamma.size() == 0) Rcpp::stop("gamma must be a finite numeric value");
  for (R_xlen_t i = 0; i < alpha.size(); ++i) {
    const double a = alpha[i];
    if (Rcpp::NumericVector::is_na(a)) continue;
    if (!(a > 0.0)) Rcpp::stop("alpha must be positive");
  }
  for (R_xlen_t i = 0; i < beta.size(); ++i) {
    const double b = beta[i];
    if (Rcpp::NumericVector::is_na(b)) continue;
    if (!(b > 0.0)) Rcpp::stop("beta must be positive");
  }
  // gamma: NA passes (caller emits NA_REAL); Inf/-Inf/NaN also pass
  // (caller treats non-finite gamma as NA, see mhn_rmhn.cpp design J).
  (void)gamma;
}

}  // namespace mhn

// [[Rcpp::export(.check_mhn_params)]]
SEXP check_mhn_params_R(SEXP alpha, SEXP beta, SEXP gamma) {
  if (!is_scalar_numeric(alpha, /*require_finite=*/false) ||
      as_double_unchecked(alpha) <= 0.0) {
    Rcpp::stop("alpha must be positive");
  }
  if (!is_scalar_numeric(beta, /*require_finite=*/false) ||
      as_double_unchecked(beta) <= 0.0) {
    Rcpp::stop("beta must be positive");
  }
  if (!is_scalar_numeric(gamma, /*require_finite=*/true)) {
    Rcpp::stop("gamma must be a finite numeric value");
  }
  return R_NilValue;
}

// [[Rcpp::export(.convert_to_gw)]]
Rcpp::List convert_to_gw_R(double alpha, double beta, double gamma) {
  return Rcpp::List::create(
    Rcpp::Named("lambda_gw") = alpha,
    Rcpp::Named("alpha_gw")  = beta,
    Rcpp::Named("beta_gw")   = -gamma
  );
}
