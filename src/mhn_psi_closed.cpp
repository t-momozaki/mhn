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

#include "mhn_psi.h"

#include <Rcpp.h>
#include <cmath>

namespace mhn {

double psi_alpha1(double gamma, double beta) {
  const double z = gamma / std::sqrt(beta);
  const double log_pnorm = R::pnorm(z / M_SQRT2, 0.0, 1.0,
                                    /*lower_tail=*/1, /*log_p=*/1);
  return std::log(2.0) + 0.5 * std::log(M_PI) + z * z / 4.0 + log_pnorm;
}

double psi_alpha2(double gamma, double beta) {
  const double z = gamma / std::sqrt(beta);
  const double pnorm_val = R::pnorm(z / M_SQRT2, 0.0, 1.0,
                                    /*lower_tail=*/1, /*log_p=*/0);
  return std::log(1.0 + std::sqrt(M_PI) * z * std::exp(z * z / 4.0) * pnorm_val);
}

}  // namespace mhn
