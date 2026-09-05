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
#include "mhn_stable.h"

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
double compute_log_peak(double alpha, double beta, double gamma, double x,
                        double* anchor_out = nullptr) {
  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double t_star;
  const double disc = gamma * gamma + 8.0 * beta * (alpha - 1.0);
  if (alpha > 1.0) {
    // Sun et al. (2023) Lemma 3b mode.
    t_star = mhn::positive_root(beta, gamma, alpha - 1.0);
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
  if (anchor_out != nullptr) *anchor_out = t_star;
  return (alpha - 1.0) * std::log(t_star) - beta * t_star * t_star
         + gamma * t_star;
}

}  // namespace

// The integral as two pieces: the log kernel at the anchor, and the log of the
// peak-normalised quadrature.  They are kept apart because the first is about
// gamma^2/(4 beta) for a large positive tilt -- 2.5e13 at gamma = 1e7 -- and any
// caller that wants the value with that term removed must never form their sum.
struct CdfParts {
  double log_peak;
  double log_val;
  bool ok;
  double log_peak_shifted = 0.0;   // the same with gamma^2/(4 beta) removed
};

CdfParts log_cdf_integrate_parts(double alpha, double beta, double gamma,
                                 double x, double tol) {
  if (!(x > 0.0)) return CdfParts{0.0, 0.0, false};
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();

  double anchor_t = 0.0;
  const double log_peak = compute_log_peak(alpha, beta, gamma, x, &anchor_t);

  // The anchor the kernel is measured from.  Where it is the interior peak the
  // exponent below is written in the centred variable, because -beta t^2 and
  // gamma t are each about gamma^2/(4 beta) there -- 5e9 at gamma = 1e5 -- and
  // subtracting the same pair at the peak leaves the difference carrying their
  // rounding.  Using the stationarity condition 2 beta m^2 - gamma m = alpha-1,
  //
  //   log g(m+s) - log g(m) = (alpha-1) [log(1 + s/m) - s/m] - beta s^2,
  //
  // in which nothing large appears.  Outside that case the direct form is used,
  // since there is no peak to centre on.
  // The interior maximum, on the same branches compute_log_peak uses: for
  // alpha > 1 always, and for alpha <= 1 only when a positive tilt creates one
  // (Sun et al. 2023, Lemma 3c).  Omitting that second case left the shifted
  // exponent forming anchor - gamma/(2 beta) at a point where the two are
  // equal, and F came out 0.9126 where it should be 1.
  const double disc_c = gamma * gamma + 8.0 * beta * (alpha - 1.0);
  const double peak_c = (alpha > 1.0)
      ? mhn::positive_root(beta, gamma, alpha - 1.0)
      : ((gamma > 0.0 && disc_c > 0.0)
             ? (gamma + std::sqrt(disc_c)) / (4.0 * beta)
             : 0.0);
  const double centre = (peak_c > 0.0) ? std::min(peak_c, x) : 0.0;
  const bool centred = (centre > 0.0)
      && std::fabs(gamma) / std::sqrt(beta) > 1.0e2;
  // The linear coefficient of the expansion about `centre`.  At the peak the
  // stationarity condition 2 beta m - gamma = (alpha-1)/m gives it without
  // subtraction; elsewhere the subtraction is between quantities whose
  // difference is of order one, so it is harmless.
  const double lin = (centre == peak_c && peak_c > 0.0)
      ? (alpha - 1.0) / centre
      : 2.0 * beta * centre - gamma;

  auto scaled_kernel = [alpha, beta, gamma, log_peak, centre, centred, lin]
                       (double t) -> double {
    if (t <= 0.0) return 0.0;
    double delta;
    if (centred) {
      // log g(c + s) - log g(c) = (alpha-1) log(1 + s/c) - beta s^2 - s (2 beta c - gamma).
      // Expanding about the peak the last term cancels against the first's
      // linear part, which is the classical form; away from it the term is what
      // makes the expansion valid at all, and dropping it -- as an earlier
      // version did -- put F(mode - 1 sd) at 0.303 where the normal limit it
      // must approach gives 0.1587.
      const double sd = t - centre;
      const double u = sd / centre;
      if (u <= -1.0) return 0.0;
      delta = (alpha - 1.0) * std::log1p(u) - beta * sd * sd - sd * lin;
    } else {
      const double log_g = (alpha - 1.0) * std::log(t)
                           - beta * t * t + gamma * t;
      delta = log_g - log_peak;
    }
    if (delta < -700.0) return 0.0;  // exp underflow guard
    return std::exp(delta);
  };

  // Split the range at the peak of the kernel when it falls strictly inside.
  // Integrated in one piece, a fixed-order rule samples [0, x] uniformly and
  // simply misses a peak that occupies a vanishing fraction of it: at
  // gamma = 1e5 the mass sits within about one unit of x = 5e4, and the CDF
  // came back as 0.  Splitting there leaves each panel monotone.
  const double peak_x = (alpha > 1.0)
      ? mhn::positive_root(beta, gamma, alpha - 1.0)
      : ((gamma > 0.0 && gamma * gamma + 8.0 * beta * (alpha - 1.0) > 0.0)
            ? (gamma + std::sqrt(gamma * gamma + 8.0 * beta * (alpha - 1.0)))
                  / (4.0 * beta)
            : 0.0);
  // On [0, x] the kernel is largest at whichever of the peak and x comes
  // first, since it increases up to the peak.  Anchoring on that point matters
  // twice over: the range is split there so each panel is monotone, and the
  // range is trimmed below it, because a fixed-order rule spread over
  // [0, 15000] cannot see mass confined to the last few units.  Below the
  // anchor the log kernel falls by at least beta s^2, so anything further than
  // sqrt(L / beta) below contributes less than exp(-L) of the total.
  //
  // When there is no interior peak -- alpha <= 1 with a non-positive tilt --
  // the kernel decreases from an unbounded value at the origin, so the mass is
  // at the left end and the range must be kept whole for tanh_sinh.
  const bool has_peak = (peak_x > 0.0);
  const double anchor = has_peak ? std::min(peak_x, x) : x;
  const bool split = has_peak && (peak_x < x);

  // Trim both ends to the reach over which the kernel is not yet negligible.
  // Trimming only below the anchor was not enough: the panel above it still ran
  // out to x, so once x left the support the rule saw nothing there either and
  // F fell back to the value of its first panel alone.  That made the
  // distribution function non-monotone -- pmhn(q, 2.5, 1, -30) held at 1 until
  // q = 3600 and then dropped to 0.301 for every larger q.
  const double L = -std::log(tol_eff) + 2.0 * std::max(0.0, alpha - 1.0) + 40.0;
  // Anchored at the peak the slope is zero and this is the Gaussian reach;
  // anchored at x, when x falls short of the peak, it is the reach the slope
  // there actually gives.
  const double reach = mhn::kernel_decay_reach_at(
      beta, mhn::kernel_decay_slope(alpha, beta, gamma, anchor), L);
  double lo = 0.0;
  if (has_peak && alpha >= 1.0) lo = std::max(0.0, anchor - reach);
  const double hi = split ? std::min(x, anchor + reach) : x;

  double val;
  if (alpha < 1.0) {
    // lo is 0 whenever alpha < 1, and the head of that range is beyond the
    // reach of any quadrature -- see mhn::kernel_head.  Left to tanh_sinh alone
    // this returned F = 6.3e-4 at (1e-6, 1, -1000, q = 1e-3) where the true
    // value is 1.  It was invisible until the normalising constant was fixed,
    // because until then the same shortfall appeared in both.
    boost::math::quadrature::tanh_sinh<double> ts;
    const mhn::KernelHead head =
        mhn::kernel_head(alpha, beta, gamma, anchor, log_peak);
    val = head.value;
    if (head.x0 < anchor) val += ts.integrate(scaled_kernel, head.x0, anchor, tol_eff);
    if (hi > anchor) val += ts.integrate(scaled_kernel, anchor, hi, tol_eff);
  } else {
    using boost::math::quadrature::gauss_kronrod;
    val = gauss_kronrod<double, 15>::integrate(
        scaled_kernel, lo, anchor, /*max_depth=*/15, tol_eff);
    if (hi > anchor) {
      val += gauss_kronrod<double, 15>::integrate(
          scaled_kernel, anchor, hi, /*max_depth=*/15, tol_eff);
    }
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return CdfParts{0.0, 0.0, false};
  }
  // The shifted exponent is formed here, at the anchor the quadrature was
  // actually normalised on.  Recomputing that anchor outside got it wrong for
  // alpha < 1 with a positive tilt, where compute_log_peak uses the interior
  // maximum of Lemma 3c and a caller guessing "the peak if alpha > 1, else x"
  // does not.
  const double d_shift = (anchor_t == peak_c && peak_c > 0.0 && gamma > 0.0)
      ? 2.0 * (alpha - 1.0) / (std::sqrt(disc_c) + gamma)
      : anchor_t - gamma / (2.0 * beta);
  return CdfParts{log_peak, std::log(val), true,
                  (alpha - 1.0) * std::log(anchor_t) - beta * d_shift * d_shift};
}

// log g(a) with gamma^2/(4 beta) removed, by completing the square:
//   (alpha-1) log a - beta a^2 + gamma a - gamma^2/(4 beta)
//       = (alpha-1) log a - beta (a - gamma/(2 beta))^2.
// At the interior peak the offset is itself a difference of two nearly equal
// quantities, which the conjugate form 2(alpha-1)/(D + gamma) avoids.
double shifted_log_peak(double alpha, double beta, double gamma, double a,
                        double x) {
  const double peak_x = (alpha > 1.0)
      ? mhn::positive_root(beta, gamma, alpha - 1.0) : 0.0;
  double d;
  if (alpha > 1.0 && gamma > 0.0 && a == peak_x && peak_x <= x) {
    const double D = std::sqrt(gamma * gamma + 8.0 * beta * (alpha - 1.0));
    d = 2.0 * (alpha - 1.0) / (D + gamma);
  } else {
    d = a - gamma / (2.0 * beta);
  }
  return (alpha - 1.0) * std::log(a) - beta * d * d;
}

double log_cdf_integrate(double alpha, double beta, double gamma,
                         double x, double tol) {
  const CdfParts p = log_cdf_integrate_parts(alpha, beta, gamma, x, tol);
  if (!p.ok) return -std::numeric_limits<double>::infinity();
  return p.log_peak + p.log_val;
}

double log_cdf_integrate_shifted(double alpha, double beta, double gamma,
                                 double x, double tol) {
  const CdfParts p = log_cdf_integrate_parts(alpha, beta, gamma, x, tol);
  if (!p.ok) return -std::numeric_limits<double>::infinity();
  return p.log_peak_shifted + p.log_val;
}

CdfParts log_ccdf_integrate_parts(double alpha, double beta, double gamma,
                                  double x, double tol) {
  if (!(x > 0.0)) x = 0.0;
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();

  // Integrate on the log axis.  With t = exp(u) the Jacobian folds into the
  // power and
  //
  //     int_x^inf t^(alpha-1) e^{-beta t^2 + gamma t} dt
  //         = int_{log x}^{inf} exp(alpha u - beta e^{2u} + gamma e^u) du,
  //
  // an integrand that is bounded, smooth and free of the t^(alpha-1)
  // singularity.  In t-space this routine needed a different range rule for
  // every regime and got two of them wrong: anchored at the mode the slope is
  // zero, so the range collapsed to the Gaussian width and the panel past the
  // mode returned nothing; and for alpha < 1 with a tiny x the range spanned
  // sixty decades that no fixed-order rule can cover.  On the log axis one rule
  // serves every case.
  //
  // The stationary point solves alpha - 2 beta e^{2u} + gamma e^u = 0, that is
  // 2 beta v^2 - gamma v - alpha = 0 with v = e^u -- note alpha, not alpha - 1,
  // because the Jacobian shifted it.  The curvature there is
  // -(4 beta v^2 - gamma v) = -(2 beta v^2 + alpha), using the stationarity
  // condition to remove the subtraction.
  const double v_peak = mhn::positive_root(beta, gamma, alpha);
  const double u_peak = std::log(v_peak);
  const double curv = 2.0 * beta * v_peak * v_peak + alpha;
  const double L = -std::log(tol_eff) + 40.0;

  const double lo = (x > 0.0) ? std::log(x)
                              : -std::numeric_limits<double>::infinity();
  const bool bad = false;

  // The integrand is largest at whichever of the peak and the lower limit comes
  // last, since past the peak it only falls.  Everything -- the scaling, the
  // upper limit and the panel widths -- is measured from that point.  Scaling by
  // the peak when the range starts above it makes every value underflow: at
  // (2.5, 1, 1) with q = 30 the scaled integrand is exp(-861) throughout, so the
  // total came back as zero and the log upper tail as -Inf.
  const double u_anchor = std::max(u_peak, lo);
  const double v_anchor = std::exp(u_anchor);
  const double log_peak =
      alpha * u_anchor - beta * v_anchor * v_anchor + gamma * v_anchor;

  // Measured from the anchor rather than formed and subtracted.  With w = e^u,
  //
  //   L(u) - L(u_a) = alpha log(w/v) - (w - v) [ beta (w - v) + (2 beta v - gamma) ],
  //
  // and at the peak 2 beta v - gamma = alpha / v by stationarity, so no large
  // quantity appears.  Formed directly, beta w^2 and gamma w are each about
  // gamma^2/(4 beta) -- 2.5e9 at gamma = 1e5 -- and their difference carries the
  // rounding of both.
  const double lin_u = (u_anchor == u_peak)
      ? alpha / v_anchor
      : 2.0 * beta * v_anchor - gamma;
  auto kernel = [alpha, beta, u_anchor, v_anchor, lin_u](double u) -> double {
    const double w = std::exp(u);
    const double dw = w - v_anchor;
    // alpha (u - u_anchor), not alpha log(w / v_anchor): the two are equal, but
    // the second is +Inf once w overflows while the rest of the exponent is
    // -Inf, and their sum is NaN rather than the zero the guard below expects.
    const double lg = alpha * (u - u_anchor) - dw * (beta * dw + lin_u);
    return (lg < -700.0 || !std::isfinite(lg)) ? 0.0 : std::exp(lg);
  };

  // How far past the anchor the integrand stays visible.  Over a distance t the
  // drop is at least slope * t + curv t^2 / 2, so
  //
  //     t = 2 L / ( slope + sqrt( slope^2 + 2 L curv ) ),
  //
  // which is the curvature reach sqrt(2 L / curv) when the slope vanishes and
  // L / slope when the slope dominates, and never divides by either.  Writing
  // it as the larger of those two separately does not work: at the peak the
  // slope is analytically zero but leaves a rounding residue of order 1e-13,
  // and L / slope then puts the upper limit 1e15 away, which no quadrature
  // survives.  The curvature is taken at the peak, which understates it above
  // the peak and so only lengthens the range.
  const double slope =
      std::fabs(alpha - 2.0 * beta * v_anchor * v_anchor + gamma * v_anchor);
  const double reach =
      2.0 * L / (slope + std::sqrt(slope * slope + 2.0 * L * curv));
  const double hi = u_anchor + reach;
  if (lo >= hi) return CdfParts{0.0, 0.0, false};

  // Below U the omitted factor exp(-beta e^{2u} + gamma e^u) is within rounding
  // of one, so that stretch integrates in closed form as exp(alpha u)/alpha and
  // needs no abscissa.  This is what makes a small x harmless however far below
  // the anchor it lies.
  const double U = std::min(
      u_anchor - reach,
      std::log(4.0e-14 / (beta + std::fabs(gamma) + 1.0)));

  double val = 0.0;
  if (lo < U) {
    // (exp(alpha U) - exp(alpha lo)) / alpha, scaled by exp(-log_peak) and
    // written so that neither exponential is formed at full size.
    const double a_lo = alpha * lo - log_peak;
    const double a_hi = alpha * U - log_peak;
    if (a_hi > -700.0) {
      val += std::exp(a_hi) * (-std::expm1(a_lo - a_hi)) / alpha;
    }
  }

  using boost::math::quadrature::gauss_kronrod;
  const double q_lo = std::max(lo, U);
  if (hi > q_lo) {
    // Panels doubling in width away from the anchor, in both directions.
    //
    // One panel per side is not enough when the range is wide in units of the
    // anchor's own scale: at (0.1, 1e-4, 100) the range is 49 wide and the peak
    // 1.4e-4, and a fixed-order rule laid over the whole of it returns half the
    // mass.  Doubling costs a number of panels logarithmic in that ratio, and
    // each sees the part of the integrand that varies on its own scale.
    const double centre = std::min(std::max(u_anchor, q_lo), hi);
    double edge = centre;
    // Start narrow and double, in both directions.  Starting at the full reach
    // put the whole range in one panel whenever the reach was a generous
    // over-estimate, and a fixed-order rule then saw nothing: at
    // (1e-6, 1, -1, q = 0.01) the reach is 10390 while the mass lies within
    // five units of the anchor, and the upper tail came back as zero.  One unit
    // of u multiplies the argument by e, so the far field never varies on a
    // scale wider than that; nearer in, the reach itself is the smaller number.
    double width = std::min(reach, 1.0);
    while (edge < hi) {
      const double next = (width > 0.0) ? std::min(hi, edge + width) : hi;
      val += gauss_kronrod<double, 31>::integrate(kernel, edge, next, 15, tol_eff);
      edge = next;
      width = (width > 0.0) ? width * 2.0 : hi - edge;
      if (!(width > 0.0)) break;
    }
    edge = centre;
    width = std::min(reach, 1.0);
    while (edge > q_lo) {
      const double next = std::max(q_lo, edge - width);
      val += gauss_kronrod<double, 31>::integrate(kernel, next, edge, 15, tol_eff);
      edge = next;
      width *= 2.0;
    }
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return CdfParts{0.0, 0.0, false};
  }
  // The shifted exponent is built here, where v_anchor is still in hand, so
  // that no caller has to recover log(val) by subtracting two large numbers --
  // which would put back the cancellation this exists to remove.
  const double d_shift = (u_anchor == u_peak && gamma > 0.0)
      ? 2.0 * alpha / (std::sqrt(gamma * gamma + 8.0 * alpha * beta) + gamma)
      : v_anchor - gamma / (2.0 * beta);
  return CdfParts{log_peak, std::log(val), true,
                  alpha * u_anchor - beta * d_shift * d_shift};
}

double log_ccdf_integrate(double alpha, double beta, double gamma,
                          double x, double tol) {
  const CdfParts p = log_ccdf_integrate_parts(alpha, beta, gamma, x, tol);
  if (!p.ok) return -std::numeric_limits<double>::infinity();
  return p.log_peak + p.log_val;
}

double log_ccdf_integrate_shifted(double alpha, double beta, double gamma,
                                  double x, double tol) {
  const CdfParts p = log_ccdf_integrate_parts(alpha, beta, gamma, x, tol);
  if (!p.ok) return -std::numeric_limits<double>::infinity();
  return p.log_peak_shifted + p.log_val;
}

}  // namespace mhn
