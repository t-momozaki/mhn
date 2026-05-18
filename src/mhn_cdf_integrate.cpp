// Numerical-integration fallback for the MHN CDF.
//
// Computes log integral_0^x t^(alpha-1) exp(-beta t^2 + gamma t) dt.
// The integrand is scaled by an estimate of its peak value so that the
// quantity passed to Boost's quadrature routines stays near unity.
//
// Quadrature selection mirrors mhn_psi_integrate.cpp:
//   * alpha >= 1: Gauss-Kronrod 15-point (smooth integrand).
//   * alpha <  1: tanh-sinh (handles t^(alpha-1) singularity at t = 0).

#include "mhn_cdf_integrate.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// Peak location of the log-integrand on (0, x].  Solves the same
// first-order condition as Sun et al. (2023) Lemma 3:
// 2 beta t^2 - gamma t - (alpha - 1) = 0.  Falls back to a small
// positive value when the density is monotone decreasing on (0, x]
// (alpha < 1 with non-positive discriminant, or gamma <= 0).
double compute_log_peak(double alpha, double beta, double gamma, double x) {
  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double t_star;
  const double disc = gamma * gamma + 8.0 * beta * (alpha - 1.0);
  if (alpha > 1.0) {
    // Sun et al. (2023) Lemma 3b mode.
    t_star = (gamma + std::sqrt(disc)) / (4.0 * beta);
  } else if (gamma > 0.0 && disc > 0.0) {
    // Sun et al. (2023) Lemma 3c interior local maximum
    // (alpha < 1, gamma > 0, alpha >= 1 - gamma^2 / (8 beta)).
    t_star = (gamma + std::sqrt(disc)) / (4.0 * beta);
  } else {
    // Monotone-decreasing branch (Sun et al. 2023, Lemma 3d): peak is
    // at the left endpoint.
    t_star = sqrt_eps;
  }
  t_star = std::min(t_star, x);
  t_star = std::max(t_star, sqrt_eps);
  return (alpha - 1.0) * std::log(t_star) - beta * t_star * t_star
         + gamma * t_star;
}

}  // namespace

double log_cdf_integrate(double alpha, double beta, double gamma,
                         double x, double tol) {
  if (!(x > 0.0)) return -std::numeric_limits<double>::infinity();
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();

  const double log_peak = compute_log_peak(alpha, beta, gamma, x);

  auto scaled_kernel = [alpha, beta, gamma, log_peak](double t) -> double {
    if (t <= 0.0) return 0.0;
    const double log_g = (alpha - 1.0) * std::log(t)
                         - beta * t * t + gamma * t;
    const double delta = log_g - log_peak;
    if (delta < -700.0) return 0.0;  // exp underflow guard
    return std::exp(delta);
  };

  double val;
  if (alpha < 1.0) {
    boost::math::quadrature::tanh_sinh<double> ts;
    val = ts.integrate(scaled_kernel, 0.0, x, /*tol=*/tol_eff);
  } else {
    using boost::math::quadrature::gauss_kronrod;
    val = gauss_kronrod<double, 15>::integrate(
        scaled_kernel, 0.0, x, /*max_depth=*/15, /*tol=*/tol_eff);
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return -std::numeric_limits<double>::infinity();
  }
  return log_peak + std::log(val);
}

}  // namespace mhn
