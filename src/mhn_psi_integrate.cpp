// Numerical integration of the Fox-Wright Psi function for gamma < 0
// (Sun et al. 2023 Supplementary, Lemma 11).
//
// Using u = sqrt(beta) * x:
//   Psi[alpha/2, z] = 2 * integral_0^Inf u^(alpha-1) exp(-u^2 - |z|*u) du
//
// Lemma 11 supplies a finite upper limit for the integral that depends
// on the requested tolerance; we then evaluate the truncated integral
// with a peak-normalized integrand and a Boost.Math quadrature rule.

#include "mhn_psi.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// Mode of u^(alpha-1) * exp(-u^2 - |z|*u).  Solves 2u^2 + |z|*u - (alpha-1) = 0.
// Mirrors R-side .integrand_mode.
double integrand_mode(double alpha, double abs_z) {
  double m;
  if (alpha > 1.0) {
    m = (-abs_z + std::sqrt(abs_z * abs_z + 8.0 * (alpha - 1.0))) / 4.0;
  } else {
    m = alpha / (abs_z + 2.0);
  }
  return std::max(m, std::numeric_limits<double>::epsilon());
}

// Upper integration limit M_u in u-space.  Mirrors R-side .psi_upper_limit.
//   alpha >= 3: Lemma 11 reformulated in u-space (error-guaranteed).
//   alpha <  3: heuristic based on exp(-u^2) decay.
double psi_upper_limit(double alpha, double abs_z, double m_u, double tol) {
  const double mach_eps = std::numeric_limits<double>::epsilon();
  const double tail_reach = std::sqrt(-std::log(mach_eps));
  const double heuristic = std::max(m_u + tail_reach, abs_z / 2.0 + tail_reach);

  if (alpha >= 3.0) {
    const double a_L = alpha * (m_u + abs_z) / (2.0 * m_u + abs_z);
    const double b_L = m_u * m_u + abs_z * m_u;
    const double log_ga = R::lgammafn(a_L);

    const double v = std::exp(2.0 * log_ga - a_L * std::log(b_L))
                     - tol * std::exp(log_ga) * (2.0 * m_u + abs_z) /
                       (2.0 * std::pow(m_u, alpha) * (m_u + abs_z));

    const double p = v / std::exp(log_ga);
    if (p > 0.0 && p < 1.0) {
      const double M_u = R::qgamma(p, a_L, /*scale=*/1.0,
                                   /*lower_tail=*/1, /*log_p=*/0) / b_L;
      return std::max(M_u, heuristic);
    }
  }

  return heuristic;
}

}  // namespace

double psi_integrate(double alpha, double beta, double gamma, double tol) {
  const double abs_z = std::fabs(gamma / std::sqrt(beta));
  const double m_u = integrand_mode(alpha, abs_z);
  const double M_u = psi_upper_limit(alpha, abs_z, m_u, tol);

  // Integrand in log space, scaled by exp(-log_peak) to keep magnitudes near 1.
  const double log_peak = (alpha - 1.0) * std::log(m_u) - m_u * m_u - abs_z * m_u;

  auto scaled_fn = [alpha, abs_z, log_peak](double u) -> double {
    if (u <= 0.0) return 0.0;
    const double log_integrand = (alpha - 1.0) * std::log(u) - u * u - abs_z * u;
    return std::exp(log_integrand - log_peak);
  };

  double val;
  if (alpha < 1.0) {
    // alpha < 1: integrand has weak endpoint singularity u^(alpha-1) at u=0.
    // tanh_sinh handles such singularities natively; matches R's QUADPACK dqags
    // (Wynn epsilon-acceleration) precision on the existing test grid.
    boost::math::quadrature::tanh_sinh<double> ts;
    val = ts.integrate(scaled_fn, 0.0, M_u, /*tol=*/tol);
  } else {
    // alpha >= 1: smooth integrand, Gauss-Kronrod is fast and accurate.
    using boost::math::quadrature::gauss_kronrod;
    val = gauss_kronrod<double, 15>::integrate(
        scaled_fn, 0.0, M_u, /*max_depth=*/15, /*tol=*/tol);
  }

  return std::log(2.0) + log_peak + std::log(val);
}

}  // namespace mhn
