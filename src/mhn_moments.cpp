// Moment functions for the MHN distribution.  All formulas come from
// Sun et al. (2023) Lemma 2 (mean, variance) and the recurrence in
// Lemma 2b (higher raw moments).  Mirrors R-side mhn/R/mhn_moments.R.

#include "mhn_check.h"
#include "mhn_psi.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// Raw moments E(X^0..X^k_max) via the Sun et al. (2023) Lemma 2b
// recurrence.  Mirrors R-side .mhn_raw_moments in mhn/R/mhn_moments.R.
std::vector<double> mhn_raw_moments(double alpha, double beta, double gamma,
                                    int k_max) {
  const double log_psi_num = mhn::mhn_log_normalizing_const(alpha + 1.0,
                                                            beta, gamma, -1.0);
  const double log_psi_den = mhn::mhn_log_normalizing_const(alpha, beta,
                                                            gamma, -1.0);
  const double mu = std::exp(log_psi_num - 0.5 * std::log(beta) - log_psi_den);

  std::vector<double> moments(static_cast<size_t>(k_max + 1), 0.0);
  moments[0] = 1.0;
  if (k_max >= 1) moments[1] = mu;

  for (int k = 0; k <= k_max - 2; ++k) {
    moments[static_cast<size_t>(k + 2)] =
        (alpha + k) / (2.0 * beta) * moments[static_cast<size_t>(k)]
        + gamma / (2.0 * beta) * moments[static_cast<size_t>(k + 1)];
  }
  return moments;
}

}  // namespace

// [[Rcpp::export(.mhn_mean_cpp)]]
double mhn_mean_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const double log_psi_num = mhn::mhn_log_normalizing_const(alpha + 1.0,
                                                            beta, gamma, -1.0);
  const double log_psi_den = mhn::mhn_log_normalizing_const(alpha, beta,
                                                            gamma, -1.0);
  return std::exp(log_psi_num - 0.5 * std::log(beta) - log_psi_den);
}

// [[Rcpp::export(.mhn_var_cpp)]]
double mhn_var_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const double mu = mhn_mean_cpp(alpha, beta, gamma);
  const double v = alpha / (2.0 * beta) + mu * (gamma / (2.0 * beta) - mu);
  return std::max(0.0, v);
}

// [[Rcpp::export(.mhn_skewness_cpp)]]
double mhn_skewness_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const auto m = mhn_raw_moments(alpha, beta, gamma, 3);
  const double mu = m[1];
  const double sigma2 = std::max(0.0, m[2] - mu * mu);
  const double sigma = std::sqrt(sigma2);
  const double central3 = m[3] - 3.0 * mu * m[2] + 2.0 * mu * mu * mu;
  return central3 / (sigma * sigma * sigma);
}

// [[Rcpp::export(.mhn_kurtosis_cpp)]]
double mhn_kurtosis_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const auto m = mhn_raw_moments(alpha, beta, gamma, 4);
  const double mu = m[1];
  const double sigma2 = std::max(0.0, m[2] - mu * mu);
  const double central4 = m[4] - 4.0 * mu * m[3]
                          + 6.0 * mu * mu * m[2] - 3.0 * mu * mu * mu * mu;
  return central4 / (sigma2 * sigma2) - 3.0;
}
