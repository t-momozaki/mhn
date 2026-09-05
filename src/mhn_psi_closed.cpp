// Closed-form Fox-Wright Psi values from Sun et al. (2023)
// Supplementary, Lemma 9(c).
//
// Lemma 9(c) states  Psi[1/2, x] = 2 sqrt(pi) exp(x^2 / 4) (1 - Phi(-x / sqrt(2))).
// Rewriting in log-form and using the standard normal symmetry
// 1 - Phi(-t) = Phi(t):
//
//   alpha = 1:                     log Psi[1/2, z]
//                                    = log(2) + 0.5 * log(pi) + z^2 / 4
//                                      + log Phi(z / sqrt(2))
//
// The alpha = 2 form below is derived from the ratio identity in
// Lemma 9(c):
//   Psi[1, x] / Psi[1/2, x] = x/2 + exp(-x^2 / 4) /
//                                   (2 sqrt(pi) (1 - Phi(-x / sqrt(2)))).
// Multiplying by Psi[1/2, x] above and simplifying gives:
//
//   alpha = 2 with gamma >= 0:     log Psi[1, z]
//                                    = log(1 + sqrt(pi) z exp(z^2 / 4)
//                                                       Phi(z / sqrt(2)))
//
// Both are evaluated entirely in log space.  The alpha = 2 form is a log of
// (1 + T) where T itself overflows a double for z > 53.3, so T is formed as a
// logarithm and recombined with log(1 + exp(.)); see psi_alpha2 below.

#include "mhn_psi.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <cmath>

namespace mhn {

double psi_alpha1(double gamma, double beta) {
  const double z = gamma / std::sqrt(beta);
  if (z <= -4.0) {
    // For z well below zero the two large terms below cancel: z^2/4 grows
    // without bound while log Phi(z/sqrt(2)) falls at the same rate, and their
    // sum is only O(log|z|).  At z = -1e8 that costs every significant digit.
    // Writing Phi(z/sqrt(2)) = erfc(|z|/2) / 2 = exp(-z^2/4) erfcx(|z|/2) / 2
    // cancels the exponential analytically, leaving
    //     log Psi[1/2, z] = 0.5 log(pi) + log erfcx(|z|/2),
    // in which nothing large appears at all.
    return 0.5 * std::log(M_PI) + std::log(erfcx_cf(-z / 2.0));
  }
  const double log_pnorm = R::pnorm(z / M_SQRT2, 0.0, 1.0,
                                    /*lower_tail=*/1, /*log_p=*/1);
  return std::log(2.0) + 0.5 * std::log(M_PI) + z * z / 4.0 + log_pnorm;
}

double psi_alpha2(double gamma, double beta) {
  const double z = gamma / std::sqrt(beta);
  if (z <= 0.0) {
    // The bracketed term is z * (something positive), so for z <= 0 it does
    // not dominate and the direct form is both accurate and in range.
    const double pnorm_val = R::pnorm(z / M_SQRT2, 0.0, 1.0,
                                      /*lower_tail=*/1, /*log_p=*/0);
    return std::log1p(std::sqrt(M_PI) * z * std::exp(z * z / 4.0) * pnorm_val);
  }
  // z > 0.  Form log of the second term rather than the term itself: the
  // z^2 / 4 in the exponent overflows a double once z exceeds about 53.3,
  // which would return +Inf and drive the whole density to zero.
  const double log_pnorm = R::pnorm(z / M_SQRT2, 0.0, 1.0,
                                    /*lower_tail=*/1, /*log_p=*/1);
  const double log_term = 0.5 * std::log(M_PI) + std::log(z)
                          + z * z / 4.0 + log_pnorm;
  return log1p_exp(log_term);
}

}  // namespace mhn
