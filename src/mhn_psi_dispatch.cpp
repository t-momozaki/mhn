// Dispatcher for log Psi[alpha/2, gamma/sqrt(beta)].  Routes each call
// to the cheapest evaluation method allowed by the parameter region,
// using the closed forms and the series / integration paths from Sun
// et al. (2023) Supplementary Lemmas 9-11.

#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>

namespace mhn {

double resolve_psi_tol(double tol) {
  return (tol > 0.0) ? tol : mhn_eps();
}

double mhn_log_normalizing_const(double alpha, double beta, double gamma,
                                 double tol) {
  const double eps = mhn_eps();
  if (is_sqrt_gamma(gamma, beta)) {
    // Negligible tilt: Psi[alpha/2, 0] = Gamma(alpha/2).  The test is on the
    // scale-free Delta = gamma / sqrt(beta), not on gamma alone -- see
    // is_sqrt_gamma in mhn_special_cases.h.
    return R::lgammafn(alpha / 2.0);
  }
  if (std::fabs(alpha - 1.0) < eps) {
    // alpha = 1 (a = 1/2): closed form for all gamma
    // (Sun et al. 2023 Supplementary, Lemma 9c).
    return psi_alpha1(gamma, beta);
  }
  if (std::fabs(alpha - 2.0) < eps && gamma >= 0.0) {
    // alpha = 2 (a = 1), gamma >= 0: closed form via the ratio in the
    // proof of Sun et al. (2023) Lemma 9c.
    return psi_alpha2(gamma, beta);
  }
  if (gamma < 0.0) {
    // gamma < 0: numerical integration
    // (Sun et al. 2023 Supplementary, Lemma 11).
    return psi_integrate(alpha, beta, gamma, resolve_psi_tol(tol));
  }
  // gamma > 0: series expansion
  // (Sun et al. 2023 Supplementary, Lemma 10).  The series declines when its
  // truncation length would be unaffordable, in which case the same quadrature
  // used for gamma < 0 answers it -- Psi depends on the tilt only through
  // z = gamma / sqrt(beta), and the integrand is well behaved for either sign.
  const double series = psi_series(alpha, beta, gamma, resolve_psi_tol(tol));
  if (std::isnan(series)) {
    return psi_integrate(alpha, beta, gamma, resolve_psi_tol(tol));
  }
  return series;
}

}  // namespace mhn

// ---------------------------------------------------------------------------
// R-visible exports.  Each accepts a sentinel tol = -1.0 meaning "use the
// default, sqrt(.Machine$double.eps)".
// ---------------------------------------------------------------------------

// [[Rcpp::export(.mhn_log_normalizing_const, rng = false)]]
double mhn_log_normalizing_const_R(double alpha, double beta, double gamma,
                                   double tol = -1.0) {
  return mhn::mhn_log_normalizing_const(alpha, beta, gamma, tol);
}

// [[Rcpp::export(.psi_series, rng = false)]]
double psi_series_R(double alpha, double beta, double gamma,
                    double tol = -1.0) {
  return mhn::psi_series(alpha, beta, gamma, mhn::resolve_psi_tol(tol));
}

// [[Rcpp::export(.psi_integrate, rng = false)]]
double psi_integrate_R(double alpha, double beta, double gamma,
                       double tol = -1.0) {
  return mhn::psi_integrate(alpha, beta, gamma, mhn::resolve_psi_tol(tol));
}
