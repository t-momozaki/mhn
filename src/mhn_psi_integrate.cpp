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
#include "mhn_stable.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// Peak of u^(alpha-1) exp(-u^2 + z u), the solution of
// 2u^2 - z u - (alpha-1) = 0.  z carries its sign here: the same quadrature
// serves both tilts, and for z > 0 the peak moves out to about z/2.
double integrand_mode(double alpha, double z) {
  double m;
  if (alpha > 1.0) {
    m = mhn::positive_root(1.0, z, alpha - 1.0);
  } else if (z > 0.0) {
    // alpha <= 1: the kernel is unbounded at the origin, but a positive tilt
    // still produces an interior maximum once z^2 > 8(1 - alpha).  Use it when
    // it exists, since that is where the mass is.
    const double disc = z * z + 8.0 * (alpha - 1.0);
    m = (disc > 0.0) ? (z + std::sqrt(disc)) / 4.0
                     : std::numeric_limits<double>::epsilon();
  } else {
    m = alpha / (-z + 2.0);
  }
  return std::max(m, std::numeric_limits<double>::epsilon());
}

// Upper integration limit M_u in u-space.
//   alpha >= 3: Lemma 11 reformulated in u-space (error-guaranteed).
//   alpha <  3: heuristic based on exp(-u^2) decay.
double psi_upper_limit(double alpha, double z, double m_u, double tol) {
  const double abs_z = std::fabs(z);
  const double mach_eps = std::numeric_limits<double>::epsilon();
  const double tail_reach = std::sqrt(-std::log(mach_eps));

  // How far past the peak the integrand needs to travel before it is
  // negligible.  Beyond the peak the exponent falls by at least
  // s^2 + |z| s, so ask for a drop of L and solve s^2 + |z| s = L:
  //
  //     s = 2 L / ( |z| + sqrt(z^2 + 4 L) ),
  //
  // written in the form that only ever adds.  For |z| = 0 this is the
  // Gaussian reach sqrt(L); for large |z| it is L / |z|, which is the point:
  // the integrand exp(-u^2 - |z| u) decays *faster* as the tilt grows, so the
  // range must shrink like 1/|z|, not grow like |z|.  Extending it to
  // |z| / 2 -- the peak of exp(-u^2 + |z| u), the opposite sign -- left the
  // whole mass inside the first fraction of a percent of the interval, where
  // a fixed-order rule does not see it.
  //
  // Writing s for the distance past the peak, the log integrand falls by
  //     (alpha-1) [ log(1 + s/m) - s/m ] - s^2,
  // the z having cancelled through the stationarity condition at m.  For
  // z >= 0 the peak sits at or beyond sqrt((alpha-1)/2) and the s^2 term
  // carries the decay, so sqrt(L) is enough.  For z < 0 the peak is squeezed
  // towards the origin and the bracketed term dominates instead, which is what
  // the expression below captures.
  const double L = -std::log(tol) + 2.0 * std::max(0.0, alpha - 1.0) + 40.0;
  double past_peak;
  if (z >= 0.0) {
    past_peak = std::sqrt(L);
  } else {
    past_peak = 2.0 * L / (abs_z + std::sqrt(abs_z * abs_z + 4.0 * L));
    past_peak = std::max(past_peak, tail_reach / (1.0 + abs_z));
  }
  const double heuristic = m_u + past_peak;

  if (alpha >= 3.0 && z < 0.0) {
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

// Shared core.  Returns log Psi minus `shift`, where `shift` is either 0 or
// z^2/4 -- the two are computed by the same quadrature, differing only in how
// the peak exponent is formed.
double psi_integrate_core(double alpha, double beta, double gamma, double tol,
                          bool shifted) {
  // Psi[alpha/2, z] = 2 int_0^inf u^(alpha-1) exp(-u^2 + z u) du with
  // z = gamma / sqrt(beta), so the quadrature depends on the tilt only through
  // z and serves either sign of gamma.
  const double z = gamma / std::sqrt(beta);
  const double m_u = integrand_mode(alpha, z);
  const double M_u = psi_upper_limit(alpha, z, m_u, tol);

  // Integrand in log space, scaled by exp(-log_peak) to keep magnitudes near 1.
  const double log_peak = (alpha - 1.0) * std::log(m_u) - m_u * m_u + z * m_u;

  // The same exponent with z^2/4 removed analytically.  Completing the square,
  //     -m^2 + z m - z^2/4 = -(m - z/2)^2,
  // and the offset m - z/2 is itself a difference of two nearly equal
  // quantities for large z.  With m = (z + D)/4 and D = sqrt(z^2 + 8(alpha-1)),
  //     m - z/2 = (D - z)/4 = 2 (alpha-1) / (D + z),
  // which only ever adds for z > 0.  Forming the difference directly would put
  // back exactly the cancellation this is here to avoid.
  double log_peak_shifted = log_peak;
  double peak_offset = 0.0;   // m_u - z/2, formed without cancelling
  if (shifted) {
    const double disc = z * z + 8.0 * (alpha - 1.0);
    double offset;
    if (z > 0.0 && disc > 0.0) {
      const double D = std::sqrt(disc);
      offset = 2.0 * (alpha - 1.0) / (D + z);
      // Guard the case where integrand_mode did not take the root branch.
      if (std::fabs((m_u - z / 2.0) - offset) > 1e-6 * std::fabs(offset) + 1e-12) {
        offset = m_u - z / 2.0;
      }
    } else {
      offset = m_u - z / 2.0;
    }
    log_peak_shifted = (alpha - 1.0) * std::log(m_u) - offset * offset;
    peak_offset = offset;
  }

  // On the shifted path, integrate in the variable centred on the peak.
  //
  //   log g(m + s) - log g(m) = (alpha-1) log(1 + s/m) - s (2m - z) - s^2,
  //
  // with 2m - z = 2 (m - z/2) = 4 (alpha-1) / (D + z), small and stable.  The
  // uncentred form evaluates -u^2 + z u at u near z/2 and subtracts the same
  // pair at the peak: at z = 3e7 both are about 1e14, so the exponent carries
  // an absolute error of 0.02 and the integrand is 2 percent wrong.  Nothing
  // large appears in the centred form.  The peak sits at m = z/2 or beyond, and
  // the range is a few units wide, so the origin -- and the u^(alpha-1)
  // singularity there -- is far outside it and needs no separate treatment.
  if (shifted) {
    const double two_m_minus_z = 2.0 * peak_offset;
    auto centred = [alpha, m_u, two_m_minus_z](double sv) -> double {
      const double t = sv / m_u;
      if (t <= -1.0) return 0.0;
      const double e = (alpha - 1.0) * std::log1p(t) - sv * two_m_minus_z
                       - sv * sv;
      return (e < -700.0) ? 0.0 : std::exp(e);
    };
    const double L = -std::log(tol) + 2.0 * std::max(0.0, alpha - 1.0) + 40.0;
    const double reach = std::min(std::sqrt(L) + 2.0, m_u);
    using boost::math::quadrature::gauss_kronrod;
    const double half = tol / 10.0;
    const double val_c =
        gauss_kronrod<double, 31>::integrate(centred, -reach, 0.0, 15, half)
      + gauss_kronrod<double, 31>::integrate(centred, 0.0, reach, 15, half);
    if (!(val_c > 0.0) || !std::isfinite(val_c)) {
      return -std::numeric_limits<double>::infinity();
    }
    return std::log(2.0) + log_peak_shifted + std::log(val_c);
  }

  auto scaled_fn = [alpha, z, log_peak](double u) -> double {
    if (u <= 0.0) return 0.0;
    const double log_integrand = (alpha - 1.0) * std::log(u) - u * u + z * u;
    return std::exp(log_integrand - log_peak);
  };

  // Split at the peak.  The scaled integrand is 1 at m_u and falls away on
  // both sides, so each panel is monotone and the rule sees the mass whatever
  // the ratio M_u / m_u happens to be.  Integrating [0, M_u] in one go left a
  // fixed-order rule sampling only the flat tail once the tilt was large.
  // Each panel is asked for a tenth of the caller's tolerance, so that their
  // sum still meets it.
  const double panel_tol = tol / 10.0;

  // Left panel, [0, m_u].  The integrand behaves like u^(alpha-1) at the
  // origin, whose derivative is unbounded for every alpha < 2, so tanh_sinh
  // is the right rule here whatever alpha is -- Gauss-Kronrod converges only
  // slowly against an endpoint of that shape.
  //
  // Quadrature alone cannot reach the head of that range as alpha approaches
  // zero -- see mhn::kernel_head, which does it in closed form.  Left to
  // tanh_sinh this returned Psi = 1.3e3 where the true value is 2.0e6, and
  // every density, probability and moment for gamma < 0 was scaled by that.
  const mhn::KernelHead head = mhn::kernel_head(alpha, 1.0, z, m_u, log_peak);
  double val = head.value;

  boost::math::quadrature::tanh_sinh<double> ts;
  if (head.x0 < m_u) val += ts.integrate(scaled_fn, head.x0, m_u, panel_tol);

  // Right panel, [m_u, M_u].  Smooth and monotone decreasing.
  if (alpha < 1.0) {
    val += ts.integrate(scaled_fn, m_u, M_u, panel_tol);
  } else {
    using boost::math::quadrature::gauss_kronrod;
    val += gauss_kronrod<double, 15>::integrate(
        scaled_fn, m_u, M_u, /*max_depth=*/15, panel_tol);
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::log(2.0) + log_peak_shifted + std::log(val);
}

double psi_integrate(double alpha, double beta, double gamma, double tol) {
  return psi_integrate_core(alpha, beta, gamma, tol, /*shifted=*/false);
}

// log Psi[alpha/2, z] - z^2/4, for z = gamma / sqrt(beta) > 0.
//
// log Psi is about z^2/4 for a large positive tilt, and so is the kernel
// exponent -beta x^2 + gamma x near the mode, while the log density they
// combine into is of order one.  Assembling it as
// (-log Psi) + (-beta x^2 + gamma x) therefore subtracts two quantities of
// size gamma^2/(4 beta) to leave something small: the absolute error of the
// result is about eps z^2/4, which passes the working tolerance at |z| = 2e4
// and reaches 5 percent by |z| = 3e7.  Removing z^2/4 from both sides
// analytically -- here, and by completing the square in the kernel -- leaves
// nothing large in the sum.
double psi_integrate_shifted(double alpha, double beta, double gamma,
                             double tol) {
  return psi_integrate_core(alpha, beta, gamma, tol, /*shifted=*/true);
}

}  // namespace mhn
