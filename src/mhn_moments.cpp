// Moment functions for the MHN distribution.
//
// The mean is the ratio of Fox-Wright values of Sun et al. (2023), p. 1597,
// Lemma 2a, E X = Psi[(alpha+1)/2, z] / Psi[alpha/2, z], for z <= 12 wherever
// the error bound of that ratio is small and, for alpha >= 1, the ratio is
// not below the mode (Lemma 4a for alpha > 1; at alpha = 1 the mean of the
// truncated normal exceeds its mode), except under a strong negative tilt
// (z < -30), where it serves only below alpha = 6 and is capped at the upper
// bound of Lemma 4a (ratio_mean, mean_impl).  For z > 12 the ratio also
// serves where the quadrature finds the mean below a thousandth of the centre
// it is measured from, under a shape so small that most of the mass lies by
// the origin (kQuadMeanMinFrac).  The variance is
// the Lemma 2c form built on it wherever that form neither cancels nor
// magnifies the ratio's error past the working tolerance and, for
// alpha >= 1, does not exceed the bound of Lemma 4c (var_impl).  Everywhere
// else, and for the skewness and kurtosis always, the moments come from
// quadrature of the unnormalised kernel, centred on its peak or, where it has
// none, taken from the origin in log x (central_moments); that needs no
// Fox-Wright value and forms no difference of large terms.  The Lemma 2b
// recurrence for higher raw moments serves only where no quadrature can be
// set up.
//
// Everything is computed for beta = 1 and z = gamma / sqrt(beta), and then
// rescaled.  Lemma 2a gives E X^k = beta^(-k/2) Psi[(alpha+k)/2, z] /
// Psi[alpha/2, z], so the k-th central moment is beta^(-k/2) times a function
// of (alpha, z) alone, and skewness and kurtosis do not depend on beta at all.
// A quadrature in absolute units would carry that power of beta inside its
// integrals: the fourth-moment integral scales like beta^(-5/2) and leaves the
// double range once |log10 beta| passes about 124.

#include "mhn_check.h"
#include "mhn_psi.h"
#include "mhn_quadrature.h"
#include "mhn_special_cases.h"
#include "mhn_stable.h"
#include "mhn_vanishing.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <vector>

namespace {

const double kEps = std::numeric_limits<double>::epsilon();

// Every quadrature window below is cut where the kernel has fallen this many
// nats below its peak; exp(-45) = 2.9e-20.
constexpr double kDrop = 45.0;

// exp(-600).  Below this length scale an integrand taken in y = log x, which
// carries the factor x of dx = x dy, is rescaled (see moments_monotone).
constexpr double kTinyScale = 2.6503965530043108e-261;

// Largest z at which the mean is taken from the Lemma 2a ratio (see
// mean_impl).
constexpr double kRatioMaxZ = 12.0;

// The ratio is used as the mean only while its error bound (RatioMean) is at
// most this, a tenth of the working tolerance of 1e-8.
constexpr double kRatioMaxRelErr = 1e-9;

// Smallest z at which the variance is taken from the Lemma 2c form (see
// var_impl).
constexpr double kVarFormMinZ = -30.0;

// For z in [kVarFormMinZ, kRatioMaxZ] other than 0 the logarithms in the ratio
// grow like alpha log alpha, and from these values of alpha on the bound of
// ratio_mean rules the ratio out as the mean (kRatioMaxRelErr) and the
// Lemma 2c form as the variance; there the ratio is not evaluated at all (see
// mean_impl and var_impl).
constexpr double kMeanRatioMaxAlpha = 1e5;
constexpr double kVarFormMaxAlpha = 1000.0;

// For z < kVarFormMinZ the ratio serves as the mean only below this alpha,
// and there capped at the upper bound of Lemma 4a (see mean_impl).
constexpr double kFarNegRatioMaxAlpha = 6.0;

// For z > kRatioMaxZ the quadrature mean is kept only while it is at least
// this fraction of the centre it is measured from (see mean_impl).
constexpr double kQuadMeanMinFrac = 1e-3;

// Relative tolerance of every quadrature.  Gauss-Kronrod stops once |K - G|,
// the gap between its Kronrod and Gauss estimates, falls below it, and on the
// smooth panels used here that gap overstates the error of K by orders of
// magnitude, so sqrt(eps) returns results good to about 1e-15.
const double kGkTol = std::sqrt(kEps);

// E X for X ~ MHN(alpha, 1, z) by Lemma 2a, E X = Psi[(alpha+1)/2, z] /
// Psi[alpha/2, z], together with a bound on its relative error.
//
// The ratio is formed as exp(log Psi_num - log Psi_den).  Each logarithm is at
// best accurate to eps times its own size, and an absolute error in the
// exponent is a relative error in the ratio, so the ratio is only as good as
// eps (|log Psi_num| + |log Psi_den|), however modest the ratio itself.  Those
// logarithms grow like (alpha/2) log(alpha/2) along the shape axis and like
// z^2/4 under a positive tilt.  The evaluator of log Psi adds errors of its
// own.  For z > 0, whether it is a series or a quadrature, the factor 8
// covers them against high-precision references.  For z < 0 log Psi comes
// from a quadrature whose error, a few units in the last place of
// max(1, |log Psi|) (see psi_log_axis in mhn_psi_integrate.cpp), does not
// scale with the logarithms, and can exceed the bound where they are small:
// for alpha up to 1.2 and z from -30 to 0 the ratio was within 1.7e-14 of
// arbitrary-precision references at 1600 points, against a bound of 1.0e-14
// at the worst of them.  That is far inside the tolerance of the mean, and
// var_impl, which magnifies it, allows for it.
//
// One case is handled apart.  At z = 0, X^2 ~ Gamma(alpha/2, 1) (Sun et al.
// 2023, Lemma 6a), so E X = Gamma((alpha+1)/2) / Gamma(alpha/2), which is
// sqrt(pi) / B(alpha/2, 1/2).  For a large alpha R's lbeta forms the logarithm
// of that beta function without the cancellation (it expands both lgamma
// values by Stirling's series and cancels their leading terms analytically),
// whereas the difference of the two lgamma values loses eps times their size:
// at alpha = 1e8 the bound above is 3e-6, and 2e-14 with lbeta, whose result
// is good to a few units in the last place and is charged 8 eps (1 + |log|).
// Below alpha = 1 the general path is kept: its logarithms are modest there,
// and lbeta's Gamma(alpha/2) overflows for an alpha near the smallest doubles.
struct RatioMean {
  double mu;
  double rel_err;
};

// log_psi_alpha, when finite, is log Psi[alpha/2, z] as the caller already
// holds it, and is used in place of a second evaluation of the same value.
RatioMean ratio_mean(double alpha, double z,
                     double log_psi_alpha =
                         std::numeric_limits<double>::quiet_NaN()) {
  if (z == 0.0 && alpha >= 1.0) {
    const double lm = M_LN_SQRT_PI - R::lbeta(0.5 * alpha, 0.5);
    return RatioMean{std::exp(lm), 8.0 * kEps * (1.0 + std::fabs(lm))};
  }
  const RatioMean failed{std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()};
  try {
    const double lpn = mhn::mhn_log_normalizing_const(alpha + 1.0, 1.0, z, -1.0);
    const double lpd = std::isfinite(log_psi_alpha)
        ? log_psi_alpha
        : mhn::mhn_log_normalizing_const(alpha, 1.0, z, -1.0);
    // A logarithm that overflowed leaves no ratio at all, not a ratio of 0.
    if (!std::isfinite(lpn) || !std::isfinite(lpd)) return failed;
    const double rel_err = 8.0 * kEps * (std::fabs(lpn) + std::fabs(lpd));
    return RatioMean{std::exp(lpn - lpd), rel_err};
  } catch (const std::exception&) {
    // A Fox-Wright value that cannot be evaluated leaves the callers to the
    // quadrature, which does not need it.
    return failed;
  }
}

// Raw moments E X^0 .. E X^k_max of MHN(alpha, 1, z) by the Lemma 2b
// recurrence E X^(k+2) = ((alpha + k) E X^k + z E X^(k+1)) / 2.
std::vector<double> raw_moments(double alpha, double z, int k_max) {
  std::vector<double> moments(static_cast<size_t>(k_max + 1), 0.0);
  moments[0] = 1.0;
  if (k_max >= 1) moments[1] = ratio_mean(alpha, z).mu;
  for (int k = 0; k <= k_max - 2; ++k) {
    moments[static_cast<size_t>(k + 2)] =
        0.5 * (alpha + k) * moments[static_cast<size_t>(k)]
        + 0.5 * z * moments[static_cast<size_t>(k + 1)];
  }
  return moments;
}

// Central moments E[(X - mu)^k], k = 2, 3, 4, of X ~ MHN(alpha, 1, z), by
// quadrature of the kernel x^(alpha-1) exp(-x^2 + z x).
//
// Forming them from raw moments -- m3 - 3 mu m2 + 2 mu^3, and the fourth-order
// analogue -- cancels catastrophically as soon as sigma is small next to mu,
// which is exactly what a large tilt produces.  At (alpha, z) = (3, 1000) the
// mean is about 500 and the standard deviation about 0.7, so mu^4 is 6e10
// while the fourth central moment is under 1: every significant digit is lost,
// and the raw-moment form gives an excess kurtosis of -701 for a quantity that
// cannot fall below -2.
//
// Quadrature of (x - mu)^k against the kernel has no such subtraction.  The
// kernel is divided by its value at a point near its peak so that it never
// overflows, and that factor cancels between numerator and denominator, so Psi
// is not needed.  Each offset x - mu is divided by a length h of the order of
// the integration window before it is raised to a power, so every integral is
// of the order of the window's width rather than of its (k+1)-th power, and
// the moments are returned in units of h.  The scale of the law runs from
// alpha / |z| under a strong negative tilt to order one, and raised to the
// fifth power the former leaves the double range by |z| = 1e62.
//
// `order` is the highest moment needed: 1 for the mean alone, 2 for the
// variance, 4 for skewness and kurtosis.
struct CentralMoments {
  double mu;          // E X
  double h;           // unit of the central moments below
  double c2, c3, c4;  // E[(X - mu)^k] / h^k
  bool ok;
  double centre = 0.0;  // mu = centre + h du, the two parts kept for callers
  double du = 0.0;      // that rescale them (see mean_impl)
};

// d^k for the small k used here.
inline double ipow(double d, int k) {
  double r = 1.0;
  for (int i = 0; i < k; ++i) r *= d;
  return r;
}

// Moments from `integrate(k, shift)`, which returns the integral of
// (u - shift)^k against the kernel over the whole range, u being the offset
// from `centre` in units of h.
template <class Integrate>
CentralMoments moments_from(double centre, double h, int order,
                            const Integrate& integrate) {
  CentralMoments out{0.0, h, 0.0, 0.0, 0.0, false};
  const double n0 = integrate(0, 0.0);
  if (!(n0 > 0.0) || !std::isfinite(n0)) return out;
  // The mean as an offset from the centre, which is the accurate way to hold it
  // when the centre is the peak and the offset small.
  const double du = integrate(1, 0.0) / n0;
  out.mu = centre + h * du;
  out.centre = centre;
  out.du = du;
  // Where the mean lies far below a positive centre -- for alpha < 1 under a
  // positive tilt with alpha so small that the mass by the origin outweighs
  // the peak's by far, as at alpha = 1e-120, z = 10 -- it is the small
  // difference of centre and -h du, and holds their rounding, some 1e-15
  // centre.  Measured from it, the second central moment, about mu centre
  // for such a law, is swamped once mu falls below about 1e-22 centre.  Below
  // 1e-12 centre the callers fall back on the raw moments instead, which
  // serve such a law well, its standard deviation being many times its mean.
  out.ok = std::isfinite(out.mu) && out.mu > 0.0 &&
           !(out.mu < 1e-12 * centre);
  if (order < 2 || !out.ok) return out;
  out.c2 = integrate(2, du) / n0;
  out.ok = std::isfinite(out.c2) && out.c2 > 0.0;
  if (order < 4 || !out.ok) return out;
  out.c3 = integrate(3, du) / n0;
  out.c4 = integrate(4, du) / n0;
  out.ok = std::isfinite(out.c3) && std::isfinite(out.c4);
  return out;
}

// The part of the range next to the origin, [0, x0], in closed form.  The
// factor x^(alpha-1) puts mass all the way down to x = exp(-1/alpha), below
// the smallest double for small alpha, where no quadrature can place an
// abscissa.  mhn::kernel_head integrates x^(a-1) exp(-x^2 + z x) over [0, x0]
// in closed form, so taken at a = alpha + j it gives the integral of x^j times
// the kernel; a moment integrand ((x - c)/h)^k, expanded in powers of x, is a
// sum of those, sum_j C(k, j) (-c/h)^(k-j) int (x/h)^j kernel dx.  Each term
// after the first carries a further factor of about x0/c, so the sum does not
// cancel.
struct OriginHead {
  double x0;         // where the closed form stops
  double moment[5];  // int_0^x0 (x/h)^j kernel dx / exp(log_scale), j = 0..4
};

OriginHead origin_head(double alpha, double z, double x_limit, double h,
                       double log_scale) {
  OriginHead out{0.0, {0.0, 0.0, 0.0, 0.0, 0.0}};
  const double log_h = std::log(h);
  for (int j = 0; j <= 4; ++j) {
    const mhn::KernelHead head =
        mhn::kernel_head(alpha + j, 1.0, z, x_limit, log_scale + j * log_h);
    out.x0 = head.x0;  // the same for every j: it depends on z and x_limit only
    out.moment[j] = head.value;
  }
  return out;
}

// The head's share of the integral of ((x - c)/h)^k against the kernel, with
// a = -c/h.
double head_share(const OriginHead& head, int k, double a) {
  static const double binom[5][5] = {{1, 0, 0, 0, 0}, {1, 1, 0, 0, 0},
                                     {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0},
                                     {1, 4, 6, 4, 1}};
  double s = 0.0;
  for (int j = 0; j <= k; ++j) s += binom[k][j] * ipow(a, k - j) * head.moment[j];
  return s;
}

// The kernel decreases from the origin: no interior peak, so alpha <= 1 with
// z <= 0, or alpha < 1 with 0 < z <= sqrt(8 (1 - alpha)) (Sun et al. 2023,
// Lemma 3d, p. 1598, for z < 0 and for z > 0 with z^2 < 8 (1 - alpha); at
// z = 0 the derivative of the log kernel, -(2 x^2 + 1 - alpha) / x, is
// negative, and at z^2 = 8 (1 - alpha) Lemma 3c's two stationary points
// merge into an inflection point; see kernel_interior_peak).  (It also
// serves a peak below the smallest normal double, see central_moments.)  As
// the tilt grows more negative the law tends to Gamma(alpha, |z|), so its
// scale is alpha / |z| rather than anything Gaussian, and that is the
// reference length used here.
//
// The stretch [0, x0] next to the origin is taken in closed form
// (origin_head); without it the variance at (1e-6, 1, -1) is 1582 times too
// large.  The rest, [x0, hi0], is integrated in y = log x, in which
// x^(alpha-1) dx becomes x^alpha dy: smooth, and for small alpha nearly flat
// across the many decades over which the mass then spreads, where in x the
// integrand rises like 1/x towards x0.
CentralMoments moments_monotone(double alpha, double z, int order) {
  const double ref = alpha / (std::fabs(z) + std::sqrt(2.0 * alpha));
  // A scale alpha / |z| below the smallest normal double.  The quadrature
  // below would then run over abscissae that hold only some of their digits,
  // or none once they underflow.  Under a tilt strong enough that the law is
  // Gamma(alpha, rate |z|) to double precision -- its relative corrections,
  // at most 8 (alpha + 1) / z^2 (in the excess kurtosis), below eps / 4 --
  // the moments are those of the limit, returned in the unit h = 1 / |z|:
  // mean alpha h, central moments alpha h^2, 2 alpha h^3 and
  // 3 alpha (alpha + 2) h^4.  The mean itself may underflow; it is also
  // returned as centre 0 plus h du, du = alpha, for callers that rescale it.
  if (!(ref >= std::numeric_limits<double>::min()) && z < 0.0 &&
      z * z * std::numeric_limits<double>::epsilon() > 32.0 * (alpha + 1.0)) {
    const double h = 1.0 / -z;
    return CentralMoments{alpha * h, h, alpha, 2.0 * alpha,
                          3.0 * alpha * (alpha + 2.0), true, 0.0, alpha};
  }
  if (!(ref > 0.0) || !std::isfinite(ref)) {
    return CentralMoments{0.0, 1.0, 0.0, 0.0, 0.0, false};
  }
  // Where the factor exp(z x - x^2) has fallen by kDrop: the positive root of
  // x^2 - z x - kDrop = 0, for either sign of z.  positive_root solves
  // 2 beta x^2 - gamma x - c = 0, so it is taken with beta = 1/2; it forms
  // the root without overflow at every z, about kDrop / |z| for a strong
  // negative tilt, where the Gamma tail runs out to many times alpha / |z|.
  const double R = mhn::positive_root(0.5, z, kDrop);
  const double hi0 = std::max(R, 40.0 * ref);
  // The kernel at x = ref, by which every integral below is divided.  In
  // y = log x the integrand carries the factor x of dx = x dy, so relative to
  // this it is about x: ref itself where the mass is, and x0 = 4.5e-9 / |z| at
  // the edge of the head.  Once |z| passes about 1e295 that falls below the
  // exp(-700) at which the integrand is set to 0, and the lower part of the
  // range is cut off: at (alpha, z) = (1e-3, -1e300) the skewness would come
  // out 0.74 of its value.  Wherever ref is below exp(-600) the factor ref is
  // therefore divided out as well, which leaves the integrand about 1 where
  // the mass is.
  double peak0 = (alpha - 1.0) * std::log(ref) - ref * ref + z * ref;
  if (ref < kTinyScale) peak0 += std::log(ref);
  // The unit the offsets are measured in before they are raised to the
  // fourth power.  ref is the scale of the law, but the range runs out to
  // hi0, about kDrop / |z| under a strong negative tilt, and for a tiny alpha
  // that is far past ref: (hi0 / ref)^4 overflows once alpha is below about
  // 1e-75 (at (1e-100, -1e10) hi0 / ref is 4.5e101), and the skewness then
  // comes out as noise.  The unit is therefore at least 1e-70 hi0, which
  // keeps (hi0 / h)^4 below 1e280; the moments in that unit, of order
  // (ref / h)^k, stay far inside the range down to the smallest alpha.
  const double h = std::max(ref, 1e-70 * hi0);
  const OriginHead head = origin_head(alpha, z, hi0, h, peak0);
  if (!(head.x0 > 0.0)) return CentralMoments{0.0, 1.0, 0.0, 0.0, 0.0, false};
  const double y0 = std::log(head.x0);
  const double y1 = std::log(hi0);
  auto integrate = [&](int k, double shift) -> double {
    double v = head_share(head, k, -shift);
    if (head.x0 < hi0) {
      v += mhn::gk_unit<31>(
          [&](double y) {
            const double x = std::exp(y);
            const double e = alpha * y - x * x + z * x - peak0;
            return (e < -700.0) ? 0.0 : ipow(x / h - shift, k) * std::exp(e);
          },
          y0, y1, kGkTol);
    }
    return v;
  };
  return moments_from(0.0, h, order, integrate);
}

// The kernel has an interior peak m (Sun et al. 2023, Lemma 3b for alpha > 1,
// 3c for alpha < 1, p. 1598; Lemma 6b at alpha = 1, p. 1599).
//
// Work in s = x - m.  Using the stationarity condition at the peak,
// z - 2 m = -(alpha-1)/m, the log kernel relative to its peak is
//
//     L(m + s) - L(m) = (alpha-1) [ log(1 + s/m) - s/m ] - s^2,
//
// in which nothing large appears.  Evaluated as written, L(x) - L(m) subtracts
// x^2 and z x, both of order 1e13 at z = 1e7, to leave a quantity of order 1,
// and the exponent would carry an absolute error of 5e-3.
//
// The window reaches from the peak as far on each side as kernel_decay_reach
// gives, which follows the decay the kernel actually has: a negative tilt
// squeezes the peak towards the origin and the linear term takes over, so at
// (1.01, -100) the curvature alone would set a window some four hundred times
// too wide.  For alpha >= 1 the kernel has fallen by at least kDrop nats at
// both edges.  The window is split at the peak, and each side is then
// monotone for Gauss-Kronrod.
//
// For alpha < 1 the factor x^(alpha-1) is unbounded at the origin and the
// bracket above is positive: left of the peak the kernel falls to a local
// minimum at x_min = (1 - alpha) / (2 m) -- the other root of
// 2 x^2 - z x + 1 - alpha = 0, whose two roots multiply to (1 - alpha) / 2 --
// and then rises without bound.  The window stops at x_min, so that its left
// side is still monotone (and smooth, unless alpha is close to 1: see the
// last paragraph below), and the stretch [0, x_min] is integrated on its
// own, in log x rather than s: s = x - m holds x only to within eps m, while
// the mass by the origin sits at x down to exp(-1/alpha).  That stretch is
// the closed-form head of origin_head over [0, x0] and the rest, as in
// moments_monotone, in y = log x.
//
// kernel_decay_reach then gives the Gaussian reach sqrt(kDrop), and the
// positive bracket leaves an edge at s only kDrop - (1-alpha) (t - log(1+t))
// nats below the peak, t = s/m.  On the right that is at least 37.8: the peak
// lies at m >= sqrt((1 - alpha)/2), its value where Lemma 3c's condition
// holds with equality, so t <= sqrt(kDrop) / sqrt((1 - alpha)/2), and the
// product (1-alpha) (t - log(1+t)) is then largest as alpha -> 0, where
// t = 9.49 and the product is 7.14.  On the left the window stops at -reach
// only when m - reach exceeds x_min, so 1 + t > x_min / m, and the product
// is largest, 3.53, as alpha -> 0 with m - reach just above x_min: the edge
// lies at least 41.4 nats down.  Both are negligible.
//
// The origin can carry mass even when the window stops short of x_min.  For
// 0 < x <= x_min the kernel relative to its peak is at most
// (m/x)^(1-alpha) exp(-(m - x_min)^2), which integrates to at most
// m^(1-alpha) x_min^alpha / alpha exp(-(m - x_min)^2), while right of the
// peak the kernel is at least exp(-s^2) and holds at least sqrt(pi)/2.  (The
// kernel rises from x_min to the peak, so between x_min and the window it
// stays below its value at the window's edge and holds at most m exp(-41.4).)
// So [0, x_min] is needed only when the window reaches x_min, or when the
// bound above exceeds exp(-kDrop), which takes a small alpha: below 1e-22 at
// z = 20, below 1e-76 at z = 30.
//
// For 1 < alpha < 3 the factor x^(alpha-1) is bounded, but it has a corner at
// the origin: its first derivative is unbounded there for alpha < 2, its
// second for alpha < 3 (alpha = 2 excepted).  And when the peak is close to
// the origin, (1 + s/m)^(alpha-1) keeps changing on the right across many
// multiples of m, as far as the tail runs: reach / m grows like
// kDrop / (alpha - 1).  Neither is smooth on the scale of a window that
// reaches the origin, and Gauss-Kronrod in s bisects towards them level after
// level: at (1.001, -1000), m = 1e-6 and reach = 0.045.  In y = log x both are
// smooth.  So where the corner matters, the window keeps only [m/2, 2m] around
// the peak, and [0, m/2] (head and y, as above) and [2m, m + reach] (in y) are
// taken apart.  For alpha < 2 the corner matters wherever its weight, about
// exp((alpha-1) - m^2) of the peak, exceeds eps: at (1.5, 8), with a weight of
// 1e-7, the window alone leaves the skewness 1e-8 and the excess kurtosis 1e-7
// adrift, relative.  For 2 < alpha < 3 the milder corner leaves the window
// alone up to 3e-10 off in the skewness and 6.3e-10 in the excess kurtosis,
// absolute, against 300-bit references (alpha from 2.001 to 2.99, z from 1.5
// to 12): relative errors of up to 4.2e-9 in the skewness and 2.2e-10 in
// c4 / c2^2, though 1.7e-8 of an excess kurtosis as small as 0.035.  Taking
// the corner apart costs about twice the time at most such points, so it is
// done only where the corner also costs bisections: where it weighs as much
// as the peak, m^2 < alpha - 1, which a negative tilt brings about.
//
// For alpha just below 1 the trouble lies next to x_min, which then sits far
// below the peak.  Near it the factor (x/m)^(alpha-1) bends sharply: its
// log-derivative (alpha-1)/x is -2 m at x_min against (alpha-1)/m at the
// peak.  A window reaching down to x_min takes that bend in s, and against
// 300- and 400-bit references it leaves the skewness and excess kurtosis up
// to 6e-9 off -- 2e-6 relative at (0.999, 6.5), where the excess kurtosis is
// -8.3e-4.  The error falls with the weight of the stretch, about
// exp((alpha-1) - m^2) of the peak as above, but at a weight of 1e-10 it is
// still up to 6e-13, and 4.6e-8 of an excess kurtosis of 2.8e-6 near
// (0.999, 9.6); only below a weight of 1e-12 is it at most about 1e-14.  So
// where x_min < m/64 and the stretch weighs more than 1e-12 of the peak, the
// split moves up to m/2, and [x_min, m/2] joins [0, x_min] in y = log x,
// where the bend is smooth.  The split costs up to 2.5 times the time of the
// window alone.  Above m/64 the bend is gentle on the scale of the window,
// and moving the split there as well would cost up to 1.4 times the time.
//
// The factor bends right of the peak too, on the scale m, and just past the
// tilt z = sqrt(8 (1 - alpha)) at which the peak appears m is far smaller
// than the Gaussian reach sqrt(kDrop).  A panel [0, reach] in s then holds
// the bend of (1 + s/m)^(alpha-1) and can meet the tolerance of
// Gauss-Kronrod without resolving it: for 1 - alpha below 1e-3 and z below
// about 0.3 that leaves the skewness and excess kurtosis up to 1e-7 off.  So
// wherever the stretch by the origin is taken apart and the reach exceeds
// m, the window stops at 2 m on the right as well, and [2 m, m + reach] is
// taken in y = log x, as for 1 < alpha < 3 above, at up to about 1.5 times
// the time.  Against 320-bit references, for alpha from 1e-6 to 1 - 1e-8
// and z up to 40, the skewness and excess kurtosis are then within 1.4e-13
// where the split moves up to m/2 and within 1.2e-11 where it stays at
// x_min, relative to the larger of 1 and their size.
CentralMoments moments_about_peak(double alpha, double z, double m, int order) {
  const double reach = mhn::kernel_decay_reach(alpha, 1.0, m, kDrop);
  // `origin`: the stretch [0, x_split] by the origin is integrated on its own,
  // and the window stops there on the left.  `far_right`: so is
  // [2 m, m + reach], and the window stops at 2 m on the right.
  double x_split = 0.0;
  bool origin = false;
  bool far_right = false;
  if (alpha < 1.0) {
    x_split = (1.0 - alpha) / (2.0 * m);
    origin = x_split - m >= -reach ||
             (1.0 - alpha) * std::log(m) + alpha * std::log(x_split)
                     - std::log(alpha) - (m - x_split) * (m - x_split) > -kDrop;
    // The bend at x_min, for alpha just below 1 (see above).
    if (origin && x_split < m / 64.0 &&
        (alpha - 1.0) - m * m > std::log(1e-12)) {
      x_split = 0.5 * m;
    }
    // The bend right of the peak, where the window would reach past it.
    if (origin) far_right = reach > m;
  } else if (alpha > 1.0 && alpha < 3.0 && alpha != 2.0) {
    // Next to the origin the log kernel relative to its peak is
    // (alpha-1) (log(x/m) + 1) - m^2 + O(x): the corner x^(alpha-1) enters
    // with a weight of about exp(corner).
    const double corner = (alpha - 1.0) - m * m;
    origin = corner > (alpha < 2.0 ? std::log(kEps) : 0.0);
    if (origin) {
      x_split = 0.5 * m;
      far_right = reach > m;
    }
  }
  const double s_lo = std::max(x_split - m, -reach);
  const double s_hi = far_right ? m : reach;
  if (!(s_hi > s_lo)) return CentralMoments{0.0, 1.0, 0.0, 0.0, 0.0, false};
  const double h = reach;
  const double inv_h = 1.0 / h;

  auto w = [alpha, m](double s) -> double {
    const double t = s / m;
    if (t <= -1.0) return 0.0;
    const double e = (alpha - 1.0) * mhn::log1p_minus_u(t) - s * s;
    return (e < -700.0) ? 0.0 : std::exp(e);
  };
  auto window = [&](int k, double shift) -> double {
    auto f = [&](double s) { return ipow(s * inv_h - shift, k) * w(s); };
    return mhn::gk_unit<31>(f, s_lo, 0.0, kGkTol) + mhn::gk_unit<31>(f, 0.0, s_hi, kGkTol);
  };

  if (!origin) return moments_from(m, h, order, window);

  const double log_peak = (alpha - 1.0) * std::log(m) - m * m + z * m;
  // Relative to the peak, the integrand in y = log x below carries the factor
  // x of dx = x dy, and is about x there, as is the head relative to its
  // width.  Under a negative tilt beyond about 1e290 (alpha - 1) the peak m,
  // about (alpha - 1)/|z|, is so small that parts of both fall below the
  // exp(-700) at which they are set to 0, as in moments_monotone.  Wherever m
  // is below exp(-600) both are therefore taken divided by m, and multiplied
  // by it once integrated (jac); elsewhere jac is 1.
  const bool tiny = m < kTinyScale;
  const double log_jac = tiny ? std::log(m) : 0.0;
  const double jac = tiny ? std::exp(log_jac) : 1.0;
  const OriginHead head = origin_head(alpha, z, x_split, h, log_peak + log_jac);
  if (!(head.x0 > 0.0)) return CentralMoments{0.0, 1.0, 0.0, 0.0, 0.0, false};
  // The same log kernel in x, (alpha-1) [log(x/m) - (x - m)/m] - (x - m)^2,
  // taken in y = log x with the factor x of dx = x dy.
  const double log_m = std::log(m);
  auto in_log_x = [&](int k, double shift, double ya, double yb) -> double {
    return mhn::gk_unit<31>(
        [&](double y) {
          const double x = std::exp(y);
          const double d = x - m;
          const double e = alpha * y - (alpha - 1.0) * (log_m + d / m) - d * d
                           - log_jac;
          return (e < -700.0) ? 0.0 : ipow(d * inv_h - shift, k) * std::exp(e);
        },
        ya, yb, kGkTol);
  };
  const double y0 = std::log(head.x0);
  const double y1 = std::log(x_split);
  auto integrate = [&](int k, double shift) -> double {
    // u - shift = (x - c) / h with c = m + h shift.
    double v =
        jac * head_share(head, k, -(m * inv_h + shift)) + window(k, shift);
    if (head.x0 < x_split) v += jac * in_log_x(k, shift, y0, y1);
    if (far_right) {
      v += jac * in_log_x(k, shift, std::log(2.0 * m), std::log(m + reach));
    }
    return v;
  };
  return moments_from(m, h, order, integrate);
}

// A peak below the smallest normal double is taken from the origin as well.
// The expansion of moments_about_peak uses the stationarity condition at m,
// which holds for the exact peak; under a negative tilt near the largest
// double and a shape just above 1 the peak, about (alpha - 1)/|z|, is
// subnormal and holds only some of its digits -- 17 bits at
// (1 + 1e-10, -1e308), where it is 1e-318 -- and a peak that far off acts on
// the law like a tilt that far off (the mean would come out 1.3e-6 low).
// Taken in y = log x from the origin, as for a kernel without a peak, the
// integrand is the kernel itself, which needs no stationarity, and the
// narrow rise to such a peak is smooth in y.
CentralMoments central_moments(double alpha, double z, int order) {
  try {
    const double m = mhn::kernel_interior_peak(alpha, 1.0, z);
    if (m >= std::numeric_limits<double>::min() && std::isfinite(m)) {
      return moments_about_peak(alpha, z, m, order);
    }
    return moments_monotone(alpha, z, order);
  } catch (const std::exception&) {
    // Boost's quadrature reports an unusable range by throwing; the callers
    // then fall back on the closed forms.
    return CentralMoments{0.0, 1.0, 0.0, 0.0, 0.0, false};
  }
}

// The moments where z = gamma / sqrt(beta) itself overflows, that is where
// |gamma| exceeds the largest double times sqrt(beta), which needs beta < 1.
// Under such a negative tilt the law of X is Gamma(alpha, rate |gamma|) to
// within relative corrections of order (alpha + 1) / z^2 (see the Gamma-limit
// tests), below 1e-307 for every alpha: mean alpha / |gamma|, variance
// alpha / gamma^2, skewness 2 / sqrt(alpha) and excess kurtosis 6 / alpha.
// Under such a positive one it is normal about its mode m to within
// corrections of order 1 / (beta m^2), which is 0 in double precision:
// mean m, variance 1 / (2 beta + (alpha - 1) / m^2), and skewness and excess
// kurtosis 0.  m, from mhn::kernel_interior_peak, is +Inf where it lies beyond
// the largest double, and so is the mean.
struct LimitMoments {
  double mean, var, skew, kurt;
};

LimitMoments infinite_tilt_moments(double alpha, double beta, double gamma) {
  if (gamma < 0.0) {
    const double rate = -gamma;
    const double mean = alpha / rate;
    return LimitMoments{mean, mean / rate, 2.0 / std::sqrt(alpha), 6.0 / alpha};
  }
  const double m = mhn::kernel_interior_peak(alpha, beta, gamma);
  return LimitMoments{m, 1.0 / (2.0 * beta + ((alpha - 1.0) / m) / m), 0.0,
                      0.0};
}

// The Lemma 2a ratio of one triple, evaluated at most once for its mean and
// its variance.
struct RatioCache {
  double log_psi_alpha = std::numeric_limits<double>::quiet_NaN();
  bool have = false;
  RatioMean value{0.0, 0.0};
  const RatioMean& get(double alpha, double z) {
    if (!have) {
      value = ratio_mean(alpha, z, log_psi_alpha);
      have = true;
    }
    return value;
  }
};

double mean_impl(double alpha, double beta, double gamma, RatioCache& rc);
double var_impl(double alpha, double beta, double gamma, RatioCache& rc);

}  // namespace

// [[Rcpp::export(.mhn_mean_cpp, rng = false)]]
double mhn_mean_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  if (mhn::vanishing_candidate(alpha)) {
    mhn::VanishingShape v;
    if (v.set(alpha, beta, gamma)) return mhn::vanishing_moments(v).mean;
  }
  RatioCache rc;
  return mean_impl(alpha, beta, gamma, rc);
}

// [[Rcpp::export(.mhn_var_cpp, rng = false)]]
double mhn_var_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  if (mhn::vanishing_candidate(alpha)) {
    mhn::VanishingShape v;
    if (v.set(alpha, beta, gamma)) return mhn::vanishing_moments(v).var;
  }
  RatioCache rc;
  return var_impl(alpha, beta, gamma, rc);
}

// E X and Var X of one triple for the search bracket of qmhn
// (mhn_qmhn.cpp), the same doubles as mhn_mean_cpp and mhn_var_cpp.  The
// Lemma 2a ratio that both may need is evaluated once, and its denominator
// log Psi[alpha/2, z] is taken from log_psi_alpha where that is finite: the
// caller's distribution-function state already holds it for the same triple.
// It is the value the ratio would compute, from the same z = gamma /
// sqrt(beta), except where z underflows to 0 while gamma does not: the
// dispatcher, which takes Psi[alpha/2, 0] = Gamma(alpha/2) at a zero tilt
// only, then sees a tilt there and none here, and the value is evaluated
// afresh.
void mhn_mean_var_cpp(double alpha, double beta, double gamma,
                      double log_psi_alpha, double* mean, double* var) {
  mhn::check_params_scalar(alpha, beta, gamma);
  RatioCache rc;
  const double z = gamma / std::sqrt(beta);
  if (!(z == 0.0 && gamma != 0.0)) rc.log_psi_alpha = log_psi_alpha;
  *mean = mean_impl(alpha, beta, gamma, rc);
  *var = var_impl(alpha, beta, gamma, rc);
}

namespace {

double mean_impl(double alpha, double beta, double gamma, RatioCache& rc) {
  const double sqrt_beta = std::sqrt(beta);
  const double z = gamma / sqrt_beta;
  if (std::isinf(z)) return infinite_tilt_moments(alpha, beta, gamma).mean;

  // The ratio of Lemma 2a (ratio_mean) is cheap and accurate while both
  // logarithms are modest.  Under a positive tilt each is about z^2/4 and at
  // best good to eps times that, so the ratio carries a relative error of
  // about eps z^2/2, which the Lemma 2c form of the variance then multiplies
  // by about z^2/2 (see var_impl); along the shape axis the logarithms grow
  // like alpha log alpha.  The quadrature forms no such ratio: it holds the
  // mean as the peak plus a small offset.  It is accurate to rounding once the
  // origin, where x^(alpha-1) is not smooth, carries no weight: for alpha >= 1
  // the kernel there is at most about exp(-z^2/4) of its peak, below eps for
  // z > 12, and for alpha < 1 the origin is integrated on its own
  // (moments_about_peak).  Measured for alpha from 1e-6 to 3000, the
  // quadrature mean agrees with 400-bit references to 1e-15 relative from
  // z = 10 on, and costs a few microseconds.  So the ratio is used up to
  // z = 12 while its error bound is at most kRatioMaxRelErr, and the
  // quadrature otherwise.  For z in [-30, 12] other than 0 that bound passes
  // kRatioMaxRelErr by alpha = 6.3e4, so from alpha = 1e5 on the ratio is not
  // even evaluated there.  (At z = 0 it comes from lbeta, see ratio_mean, and
  // serves at every alpha where it is not below the mode; see below.)
  //
  // Two bounds of Sun et al. (2023, Lemma 4a) hold the mean closer than the
  // ratio's error can.  For alpha > 1 the mean is at least the mode, and at
  // alpha = 1, where the law is a normal truncated to x > 0, it exceeds the
  // mode, z / 2 for z > 0 and 0 otherwise; under a positive tilt the two draw
  // together faster than the error of the ratio falls, and just above
  // alpha = 1, from z of about 10.7 (11.1 at alpha = 1 itself) to 12, the
  // ratio lands below the mode by up to 3.3e-14 of the mean.  At z = 0 the
  // mean exceeds the mode
  // sqrt((alpha - 1)/2) by a relative 1/(4 alpha) or so, which the few-ulp
  // error of the lbeta value covers from alpha of about 1.6e14.  So a ratio
  // below the mode is refused for alpha >= 1, at z = 0 too, and the
  // quadrature, which holds the mean as the peak plus a small offset, is
  // taken instead.
  //
  // The upper bound (gamma + sqrt(gamma^2 + 8 alpha beta)) / (4 beta), which
  // the text after the lemma notes holds for every alpha > 0, lies above a
  // mean of about alpha/|z| (1 - 2 (alpha + 1)/z^2) under a strong negative
  // tilt by only a relative 2/z^2.  For z < -30 the error of the ratio, which
  // comes mostly from the quadrature of log Psi (see ratio_mean), grows with
  // alpha -- 1e-13 at alpha = 10, 1e-12 at 300 and 1e-10 at 3e6 against
  // 200-bit references -- and, once 2/z^2 falls below it, places the mean
  // above that bound, while the quadrature mean is within 3e-15 of the
  // references from alpha = 6 to 3e6 and z = -30 to -1e14.  Below alpha = 6
  // (kFarNegRatioMaxAlpha) the ratio is kept there: for |z| up to 1e14 its
  // error stays below 5e-14, where the quadrature's reaches 3.5e-13 at some
  // shapes between 3 and 4.5.  Beyond, the ratio's error grows with log |z|,
  // as both logarithms are about -alpha log |z|, to 7e-13 near the end of the
  // double range.  That error exceeds 2/z^2 once |z| passes about 1e7, and a
  // ratio above the bound is then replaced by the bound, written as
  // 2 alpha / (sqrt(z^2 + 8 alpha) - z) so that it does not cancel.  The mean
  // lies below the bound, so this can only bring the result closer to it:
  // where the cap acts the result is within 1.2e-14 of the mean, and within
  // 3.3e-16 once |z| passes 1e8, where the gap is no larger than the rounding
  // of the bound itself.
  //
  // The mean at beta = 1 is about alpha / |z| under a strong negative tilt,
  // and can fall below the smallest normal double, 2.2e-308, where the mean
  // itself, divided by sqrt(beta) < 1, need not.  There the ratio, whose
  // exponential has already lost its digits, is not used, and the quadrature
  // mean is rescaled before it is assembled, as centre / sqrt(beta) +
  // (h / sqrt(beta)) du.
  const bool try_ratio = z <= kRatioMaxZ &&
      (z == 0.0 || alpha < (z < kVarFormMinZ ? kFarNegRatioMaxAlpha
                                             : kMeanRatioMaxAlpha));
  RatioMean ratio{0.0, std::numeric_limits<double>::infinity()};
  if (try_ratio) {
    ratio = rc.get(alpha, z);
    if (ratio.rel_err <= kRatioMaxRelErr &&
        ratio.mu >= std::numeric_limits<double>::min() &&
        !(alpha >= 1.0 &&
          ratio.mu < mhn::kernel_interior_peak(alpha, 1.0, z))) {
      if (z >= kVarFormMinZ) return ratio.mu / sqrt_beta;
      // The upper bound of Lemma 4a at beta = 1, alpha / half_den with
      // half_den = (sqrt(z^2 + 8 alpha) - z) / 2, each half taken apart so
      // that the sum cannot overflow for |z| near the largest double.
      const double half_den =
          0.5 * std::hypot(z, std::sqrt(8.0 * alpha)) - 0.5 * z;
      return std::min(ratio.mu, alpha / half_den) / sqrt_beta;
    }
  }
  const CentralMoments m = central_moments(alpha, z, 1);
  if (m.ok) {
    // The quadrature holds the mean as the centre plus an offset.  Under a
    // positive tilt and a shape so small that most of the mass lies by the
    // origin, the mean is far below the centre, the two nearly cancel, and the
    // mean keeps the rounding of the centre: a relative error of about
    // eps centre / mean, 1.1e-4 at alpha = 1e-27, z = 12.1, where the mean is
    // 2.3e-12 of the centre.  The Lemma 2a ratio forms no such difference.
    // Its error there is that of its two logarithms, each about z^2/4, and z
    // stays below about 54 for a mean that far below the centre under a shape
    // a double can hold: against quadratures of the two integrals of the
    // lemma, at 572 points with z from 12.001 to 50 and the mean from 1e-13
    // to 0.6 of the centre, the ratio was within 1.5e-13.  It is taken below
    // kQuadMeanMinFrac of the centre, where the rounding of the quadrature
    // mean, 3e-13 at that fraction, becomes the larger of the two, and where
    // its own bound allows.
    if (z > kRatioMaxZ && m.mu < kQuadMeanMinFrac * m.centre) {
      const RatioMean& r = rc.get(alpha, z);
      if (r.rel_err <= kRatioMaxRelErr &&
          r.mu >= std::numeric_limits<double>::min()) {
        return r.mu / sqrt_beta;
      }
    }
    if (m.mu >= std::numeric_limits<double>::min()) return m.mu / sqrt_beta;
    return m.centre / sqrt_beta + (m.h / sqrt_beta) * m.du;
  }
  if (!try_ratio) ratio = rc.get(alpha, z);
  return ratio.mu / sqrt_beta;
}

double var_impl(double alpha, double beta, double gamma, RatioCache& rc) {
  const double z = gamma / std::sqrt(beta);
  if (std::isinf(z)) return infinite_tilt_moments(alpha, beta, gamma).var;

  // The Lemma 2c form, at beta = 1, is Var = t1 + t2 - t3 with t1 = alpha/2,
  // t2 = z mu / 2 and t3 = mu^2.  It is exact algebra, but it inherits two
  // errors.  Rounding its terms costs about eps max(t1, |t2|, t3) absolute,
  // which is what defeats it under a strong negative tilt: t1 and t2 are each
  // about alpha/2 in size and cancel, while the answer is about alpha/z^2, and
  // at z = -1e7 the form is 17 percent high.  (Under a strong positive tilt t2
  // and t3, each about z^2/4, cancel likewise to leave about 1/2, but there
  // the second error is the larger.)  And an error delta in the
  // relative accuracy of the mean moves it by |d Var / d mu| mu delta =
  // |t2 - 2 t3| delta, where |t2 - 2 t3| / Var is about z^2/2 under a strong
  // tilt of either sign and 4 alpha at z = 0: a mean good to 1e-12 gives a
  // variance good only to 4e-9 at alpha = 1000.  The form is used only where
  // the sum of the two, with delta the bound of ratio_mean, stays below 1e-8
  // of the result, and only where that bound is at most kRatioMaxRelErr: the
  // second term is linear in delta and is formed from the mean itself, so it
  // means something only while delta is small.  Otherwise the second central
  // moment comes from the quadrature, which forms no such difference.
  //
  // Where the form would be rejected anyway the ratio is not evaluated: for
  // z > 12, where the mean itself comes from the quadrature (see
  // mean_impl), and for z != 0 from alpha = 1000 on.  There the bound on
  // delta grows like alpha log alpha and its magnification like 4 alpha, and
  // at every z in [-30, 12] their product passes 1e-8 of the variance by
  // alpha = 600 and is three times that at alpha = 1000.  (At z = 0 the mean
  // comes from lbeta, see ratio_mean, and the form serves up to alpha of about
  // 2e5.)  Nor is the form used for z < -30.  For z < 0 log Psi comes from a
  // quadrature that, for alpha up to 1.2, leaves the ratio up to about 2e-14
  // adrift (see ratio_mean), which can exceed the bound charged for it.  For
  // such alpha |t2 - 2 t3| / Var is at most 457 for z >= -30, so that costs
  // below 1e-11; beyond, it grows like z^2/2.
  //
  // For alpha >= 1 the variance is at most 1/2 at beta = 1 (Sun et al. 2023,
  // Lemma 4c), and under a positive tilt it comes within far less than the
  // error allowed for above of that bound: just above alpha = 1, from z of
  // about 10.7 (11.1 at alpha = 1 itself) to 12, the form, whose t2 and t3
  // cancel there, exceeds 1/2 by up to 1.2e-12.  A form above 1/2 is
  // therefore refused for alpha >= 1, and the quadrature taken instead.  The
  // quadrature too can land above 1/2 where the variance lies within its
  // error of that bound, next to alpha = 1 for z of about 11 to 13.5 (by
  // 1.4e-14 at alpha = 1 + 2^-52, z = 11.4489), so for alpha >= 1 its result
  // is capped at 1/2, which is then the nearer value.
  const bool use_form = z >= kVarFormMinZ && z <= kRatioMaxZ &&
                        (z == 0.0 || alpha < kVarFormMaxAlpha);
  double v = std::numeric_limits<double>::quiet_NaN();
  if (use_form) {
    const RatioMean ratio = rc.get(alpha, z);
    const double mu = ratio.mu;
    const double t1 = 0.5 * alpha;
    const double t2 = 0.5 * z * mu;
    const double t3 = mu * mu;
    v = t1 + t2 - t3;
    const double largest = std::max(t1, std::max(std::fabs(t2), t3));
    const double err =
        kEps * largest + std::fabs(t2 - 2.0 * t3) * ratio.rel_err;
    if (v > 0.0 && ratio.rel_err <= kRatioMaxRelErr && err < 1e-8 * v &&
        !(alpha >= 1.0 && v > 0.5)) {
      return v / beta;
    }
  }

  const CentralMoments m = central_moments(alpha, z, 2);
  if (m.ok) {
    // The variance is h^2 c2.  Without an interior peak h is about
    // alpha / |z| and c2 about 1 / alpha, and near z = -1e154 h^2 falls below
    // the smallest normal double while alpha / z^2 is still representable:
    // at (1e-8, -1e154) h^2 rounds to 0 against a variance of 1e-316.  There
    // h is applied to c2 one factor at a time, so that only the final product
    // can lose digits to the subnormal range.  Where even that product falls
    // below the normal range, the variance of X, var1 / beta, can still lie
    // inside it when beta < 1 (at (0.5, 1e-300, -1e50), z = -1e200, it is
    // 5e-101 while var1 underflows to 0), so h is then taken in the units of
    // X, h / sqrt(beta), before it is applied.
    const double h2 = m.h * m.h;
    double var1 = h2 >= std::numeric_limits<double>::min()
                      ? h2 * m.c2
                      : m.h * (m.h * m.c2);
    if (alpha >= 1.0) var1 = std::min(var1, 0.5);
    if (var1 < std::numeric_limits<double>::min()) {
      const double hx = m.h / std::sqrt(beta);
      return hx * (hx * m.c2);
    }
    return var1 / beta;
  }
  // Neither is usable: return the form whatever its error, kept to [0, 1/2]
  // for alpha >= 1 and to [0, inf) otherwise, or NaN where the ratio could
  // not be evaluated either.
  if (!use_form) {
    const double mu = rc.get(alpha, z).mu;
    v = 0.5 * alpha + 0.5 * z * mu - mu * mu;
  }
  if (std::isnan(v)) return v;
  v = std::max(0.0, v);
  if (alpha >= 1.0) v = std::min(v, 0.5);
  return v / beta;
}

}  // namespace

// [[Rcpp::export(.mhn_skewness_cpp, rng = false)]]
double mhn_skewness_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  if (mhn::vanishing_candidate(alpha)) {
    mhn::VanishingShape v;
    if (v.set(alpha, beta, gamma)) return mhn::vanishing_moments(v).skewness;
  }
  const double z = gamma / std::sqrt(beta);
  if (std::isinf(z)) return infinite_tilt_moments(alpha, beta, gamma).skew;
  const CentralMoments m = central_moments(alpha, z, 4);
  if (m.ok) {
    // For a tiny alpha, c2 in the unit of moments_monotone can lie so far
    // below 1 that its power 3/2 leaves the normal range; the ratio is then
    // taken one factor at a time.
    const double c2_32 = std::pow(m.c2, 1.5);
    if (c2_32 >= std::numeric_limits<double>::min()) return m.c3 / c2_32;
    return (m.c3 / m.c2) / std::sqrt(m.c2);
  }

  // Fall back on the raw-moment expansion where the quadrature could not be
  // set up.  It is the exact algebra, and accurate whenever sigma is not small
  // against the mean.
  const auto r = raw_moments(alpha, z, 3);
  const double mu = r[1];
  const double sigma2 = std::max(0.0, r[2] - mu * mu);
  const double central3 = r[3] - 3.0 * mu * r[2] + 2.0 * mu * mu * mu;
  // For a tiny alpha sigma2 is about alpha E[X^2] at alpha = 1, and its
  // power 3/2 can leave the normal range while the skewness, about
  // 1 / sqrt(alpha), does not; the ratio is then taken a factor at a time.
  const double sigma3 = std::pow(sigma2, 1.5);
  if (sigma3 >= std::numeric_limits<double>::min()) return central3 / sigma3;
  return (central3 / sigma2) / std::sqrt(sigma2);
}

// [[Rcpp::export(.mhn_kurtosis_cpp, rng = false)]]
double mhn_kurtosis_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  if (mhn::vanishing_candidate(alpha)) {
    mhn::VanishingShape v;
    if (v.set(alpha, beta, gamma)) return mhn::vanishing_moments(v).kurtosis;
  }
  const double z = gamma / std::sqrt(beta);
  if (std::isinf(z)) return infinite_tilt_moments(alpha, beta, gamma).kurt;
  const CentralMoments m = central_moments(alpha, z, 4);
  if (m.ok) {
    // As for the skewness, c2^2 can leave the normal range for a tiny alpha.
    const double c2_2 = m.c2 * m.c2;
    if (c2_2 >= std::numeric_limits<double>::min()) return m.c4 / c2_2 - 3.0;
    return (m.c4 / m.c2) / m.c2 - 3.0;
  }

  const auto r = raw_moments(alpha, z, 4);
  const double mu = r[1];
  const double sigma2 = std::max(0.0, r[2] - mu * mu);
  const double central4 = r[4] - 4.0 * mu * r[3]
                          + 6.0 * mu * mu * r[2] - 3.0 * mu * mu * mu * mu;
  // As for the skewness, sigma2^2 can leave the normal range.
  const double sigma4 = sigma2 * sigma2;
  if (sigma4 >= std::numeric_limits<double>::min()) {
    return central4 / sigma4 - 3.0;
  }
  return (central4 / sigma2) / sigma2 - 3.0;
}
