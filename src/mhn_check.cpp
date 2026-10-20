// Implementation of the parameter-validation helpers declared in
// mhn_check.h.

#include "mhn_check.h"

#include <Rcpp.h>
#include <cmath>

namespace mhn {

void check_params_scalar(double alpha, double beta, double gamma) {
  // alpha > 0 and beta > 0 are true of Inf, so finiteness is checked
  // explicitly: an infinite shape or rate lies outside the parameter space,
  // and passed on to the evaluators it would come out as NaN, as an infinite
  // log normalising constant or as a meaningless finite value rather than as
  // an error.
  if (!(alpha > 0.0) || !std::isfinite(alpha)) {
    Rcpp::stop("alpha must be positive and finite");
  }
  if (!(beta > 0.0) || !std::isfinite(beta)) {
    Rcpp::stop("beta must be positive and finite");
  }
  if (!std::isfinite(gamma)) Rcpp::stop("gamma must be a finite numeric value");
}

bool all_params_valid(const Rcpp::NumericVector& alpha,
                      const Rcpp::NumericVector& beta,
                      const Rcpp::NumericVector& gamma) {
  for (R_xlen_t i = 0; i < alpha.size(); ++i) {
    if (!params_valid(alpha[i], 1.0, 0.0)) return false;
  }
  for (R_xlen_t i = 0; i < beta.size(); ++i) {
    if (!params_valid(1.0, beta[i], 0.0)) return false;
  }
  for (R_xlen_t i = 0; i < gamma.size(); ++i) {
    if (!params_valid(1.0, 1.0, gamma[i])) return false;
  }
  return true;
}

}  // namespace mhn
