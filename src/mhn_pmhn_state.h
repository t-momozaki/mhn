#ifndef MHN_PMHN_STATE_H
#define MHN_PMHN_STATE_H

// Shared evaluation state for the MHN CDF, used by both pmhn (linear CDF
// evaluation per element) and qmhn (the same routine called repeatedly
// inside the root finder).
//
// Defined inline in a header so that both translation units can link
// without pulling in a dedicated .cpp file.

#include "mhn_cdf_integrate.h"
#include "mhn_cdf_series.h"
#include "mhn_constants.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <limits>

namespace mhn {

enum CdfKind {
  CDF_GENERAL = 0,
  CDF_SQRT_GAMMA = 1,
  CDF_TRUNCATED_NORMAL = 2
};

struct CdfState {
  CdfKind kind = CDF_GENERAL;
  double alpha = 1.0;
  double beta = 1.0;
  double gamma = 0.0;
  // General-case state:
  double log_psi = 0.0;
  double log_prefactor = 0.0;       // log(2) + (alpha/2) log(beta) - log Psi
  // Truncated-normal state (alpha == 1):
  double tn_mu = 0.0;
  double tn_sigma = 0.0;
  double tn_log_denom = 0.0;        // log Phi(mu/sigma)

  void recompute(double a, double b, double g) {
    alpha = a; beta = b; gamma = g;
    if (is_sqrt_gamma(g)) {
      kind = CDF_SQRT_GAMMA;
      return;
    }
    if (is_truncated_normal(a)) {
      kind = CDF_TRUNCATED_NORMAL;
      tn_mu = g / (2.0 * b);
      tn_sigma = 1.0 / std::sqrt(2.0 * b);
      tn_log_denom = R::pnorm(tn_mu / tn_sigma, 0.0, 1.0, 1, 1);
      return;
    }
    kind = CDF_GENERAL;
    log_psi = mhn_log_normalizing_const(a, b, g, -1.0);
    log_prefactor = std::log(2.0) + 0.5 * a * std::log(b) - log_psi;
  }

  // Returns F(qi) in [0, 1].  Boundary handling: qi <= 0 -> 0,
  // qi = +Inf -> 1, NA -> NA_REAL.
  double cdf_linear(double qi) const {
    if (Rcpp::NumericVector::is_na(qi)) return NA_REAL;
    if (qi <= 0.0) return 0.0;
    if (qi == R_PosInf) return 1.0;
    switch (kind) {
      case CDF_SQRT_GAMMA:
        return R::pgamma(qi * qi, alpha / 2.0, /*scale=*/1.0 / beta, 1, 0);
      case CDF_TRUNCATED_NORMAL: {
        const double z_q = (qi - tn_mu) / tn_sigma;
        const double log_upper = R::pnorm(z_q, 0.0, 1.0, 0, 1) - tn_log_denom;
        const double upper = std::exp(log_upper);
        return 1.0 - upper;
      }
      case CDF_GENERAL:
      default: {
        // Always try the Sun et al. (2023) Lemma 1b series first (the
        // paper's prescribed CDF formula, truncated at Lemma 10's K).
        // The series returns NaN to signal that double-precision
        // alternating-sign cancellation has crossed the safety margin
        // for gamma < 0; in that case fall back to the Boost.Math
        // integration of the unnormalised density.  See
        // src/mhn_cdf_series.cpp for the cancellation analysis.
        const double tol = mhn_eps();
        double log_F = log_cdf_series(alpha, beta, gamma, qi, log_psi, tol);
        if (!R_finite(log_F)) {
          const double log_I = log_cdf_integrate(alpha, beta, gamma, qi, tol);
          log_F = log_prefactor + log_I;
        }
        double F = std::exp(log_F);
        if (F < 0.0) F = 0.0;
        if (F > 1.0) F = 1.0;
        return F;
      }
    }
  }
};

}  // namespace mhn

#endif  // MHN_PMHN_STATE_H
