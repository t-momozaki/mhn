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
  // The same with gamma^2/(4 beta) removed, paired with the shifted integral.
  bool shifted = false;
  double log_prefactor_shifted = 0.0;
  // Truncated-normal state (alpha == 1):
  double tn_mu = 0.0;
  double tn_sigma = 0.0;
  double tn_log_denom = 0.0;        // log Phi(mu/sigma)

  void recompute(double a, double b, double g) {
    alpha = a; beta = b; gamma = g;
    if (is_sqrt_gamma(g, b)) {
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

    // For a large positive tilt, log Psi and the log of the integral over
    // [0, x] are each about gamma^2/(4 beta) while the log probability they
    // combine into is of order one.  Assembled directly the result carries
    // their rounding: at gamma = 1e7 both are 2.5e13, and pmhn at the mode came
    // out 1.1e-3 away from 0.5 -- and on another platform, where the quadrature
    // rounds differently, far enough away to be clamped to 1.  Past a threshold
    // well inside the safe range, take both with that term removed.
    const double z = g / std::sqrt(b);
    shifted = (z > 1.0e3);
    if (shifted) {
      log_prefactor_shifted = std::log(2.0) + 0.5 * a * std::log(b)
                              - psi_integrate_shifted(a, b, g, mhn_eps());
    }
  }

  // Returns F(qi) in [0, 1].  Boundary handling: qi <= 0 -> 0,
  // qi = +Inf -> 1, NA -> NA_REAL.
  double cdf_linear(double qi) const {
    // NA and NaN are distinct in the base R d/p/q contract, and
    // Rcpp::NumericVector::is_na does not separate them -- it is ISNAN.
    if (R_IsNA(qi)) return NA_REAL;
    if (R_IsNaN(qi)) return R_NaN;
    if (qi <= 0.0) return 0.0;
    if (qi == R_PosInf) return 1.0;
    switch (kind) {
      case CDF_SQRT_GAMMA:
        return R::pgamma(qi * qi, alpha / 2.0, /*scale=*/1.0 / beta, 1, 0);
      case CDF_TRUNCATED_NORMAL: {
        // Derived from the upper tail this returned exactly 0 once the lower
        // tail fell below the rounding of 1 -- pmhn(3, 1, 1, 18) gave 0 against
        // a true 1.3e-17.  The scalar path never had that problem because it
        // computes whichever tail is the smaller directly; take the smaller one
        // here too, by the same reasoning, so the two paths agree.
        const double z_q = (qi - tn_mu) / tn_sigma;
        if (z_q > 0.0) {
          return 1.0 - std::exp(R::pnorm(z_q, 0.0, 1.0, 0, 1) - tn_log_denom);
        }
        // z_q <= 0: the lower tail is the small one.  Phi(z_q) - Phi(-mu/sigma)
        // as an ordered difference of two lower-tail normals in log space.
        const double lo_q = R::pnorm(z_q, 0.0, 1.0, 1, 1);
        const double lo_0 = R::pnorm(-tn_mu / tn_sigma, 0.0, 1.0, 1, 1);
        if (!(lo_q > lo_0)) return 0.0;
        return std::exp(lo_q + std::log1p(-std::exp(lo_0 - lo_q)) - tn_log_denom);
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
        double log_F;
        if (shifted) {
          // Where the shift is needed the series must not be used at all, even
          // if its truncation length happens to be affordable.  It returns
          // log_sum - log Psi, and for a large positive tilt those are each
          // about gamma^2/(4 beta) -- 2.5e9 at gamma = 1e5 -- so the difference
          // carries their rounding just as the unshifted quadrature assembly
          // did.  Whether the series runs is decided by a truncation length
          // computed from lgamma, which differs between platforms by an ulp:
          // here it exceeds the ceiling and defers, on Windows it did not, and
          // pmhn returned 1 for every q at gamma = 1e5.
          log_F = log_prefactor_shifted
                  + log_cdf_integrate_shifted(alpha, beta, gamma, qi, tol);
        } else {
          log_F = log_cdf_series(alpha, beta, gamma, qi, log_psi, tol);
          if (!R_finite(log_F)) {
            log_F = log_prefactor + log_cdf_integrate(alpha, beta, gamma, qi, tol);
          }
        }
        double F = std::exp(log_F);
        if (F < 0.0) F = 0.0;
        if (F > 1.0) F = 1.0;
        return F;
      }
    }
  }

  // log P(X > qi), computed from the upper tail rather than derived as
  // log(1 - F).  The derived form loses the survival probability entirely once
  // F rounds to 1, which happens as soon as it falls below about 1e-16; base R
  // exposes lower.tail precisely so that its distribution functions do not have
  // that limitation, and this is what lets pmhn and qmhn match.
  double log_upper_tail(double qi) const {
    if (R_IsNA(qi)) return NA_REAL;
    if (R_IsNaN(qi)) return R_NaN;
    if (qi <= 0.0) return 0.0;                 // log 1
    if (qi == R_PosInf) return R_NegInf;
    switch (kind) {
      case CDF_SQRT_GAMMA:
        return R::pgamma(qi * qi, alpha / 2.0, /*scale=*/1.0 / beta,
                         /*lower_tail=*/0, /*log_p=*/1);
      case CDF_TRUNCATED_NORMAL: {
        const double z_q = (qi - tn_mu) / tn_sigma;
        return R::pnorm(z_q, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/1)
               - tn_log_denom;
      }
      case CDF_GENERAL:
      default: {
        const double tol = mhn_eps();
        double log_Q = shifted
            ? log_prefactor_shifted
                  + log_ccdf_integrate_shifted(alpha, beta, gamma, qi, tol)
            : log_prefactor + log_ccdf_integrate(alpha, beta, gamma, qi, tol);
        if (log_Q > 0.0) log_Q = 0.0;
        return log_Q;
      }
    }
  }
};

}  // namespace mhn

#endif  // MHN_PMHN_STATE_H
