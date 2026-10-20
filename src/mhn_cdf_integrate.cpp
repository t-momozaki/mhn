// Numerical integration of the MHN CDF.
//
// Computes log integral_0^x t^(alpha-1) exp(-beta t^2 + gamma t) dt and its
// complement over [x, Inf).  The integrand is scaled by its value at an
// anchor point, or for a tiny shape by the size of its part next to the
// origin where that is far larger, so that the quantity passed to Boost's
// quadrature routines stays within the range of a double.  pmhn and qmhn
// use these for a shape and a tilt z = gamma / sqrt(beta) up to 1e6, while
// their terms stay below about 1e7 in size (mhn_cdf_peak.h,
// use_peak_frame); beyond, they measure both tails from the peak of the
// law in exact offsets instead.  The closed forms
// below for a law narrower than the spacing of the doubles at its peak are
// then reached only by the diagnostic hook .pmhn_force, in the lower tail.
// The upper-tail form narrow_tail also serves ordinary pmhn and qmhn calls
// far out in the upper tail, where the kernel falls by L within half a
// spacing of the doubles in log x at the lower limit.
//
// Over [0, x]:
//   * up to the point x0 of mhn::kernel_head, next to the origin, the
//     leading terms of the kernel in closed form: the whole integral for
//     x <= x0 at any shape, and its head below the quadrature for alpha < 1.
//   * alpha >= 1: Gauss-Kronrod 15-point (smooth integrand), on panels split
//     at the kernel's interior peak (Sun et al. 2023, Lemma 3b; at
//     alpha = 1, Lemma 6b) and trimmed to where the kernel carries mass.
//   * alpha <  1: tanh-sinh from x0 up to x or up to the interior peak that
//     a positive tilt can create (Lemma 3c), whichever comes first, which
//     handles the t^(alpha-1) singularity at t = 0; past that peak, where
//     the kernel is smooth and log-concave, Gauss-Kronrod 15-point.
// Over [x, Inf): Gauss-Kronrod 31-point panels on the axis u = log t.
// Over either range: the Gamma limit, by pgamma, for a tilt below
// -1e150 sqrt(beta) wherever exp(-beta t^2) is 1 to within 2^-60 over the
// part that carries the integral, and closed forms for a kernel narrower than
// the spacing of the doubles.

#include "mhn_cdf_integrate.h"
#include "mhn_constants.h"
#include "mhn_fma.h"
#include "mhn_quadrature.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// The point t* at which the kernel is normalised before integration, and
// log g(t*).
//
// Where the kernel has an interior maximum m (Sun et al. 2023, Lemma 3b for
// alpha > 1, Lemma 6b for alpha = 1 with a positive tilt, Lemma 3c for
// alpha < 1 with a positive tilt; see mhn::kernel_interior_peak) it
// increases up to m, so its largest value on
// (0, x] is at min(m, x), and that is t*.  It must be exactly the point the
// centred kernel below is expanded about, however close to the origin it lies:
// the quadrature returns the integral divided by g(t*), and the caller
// multiplies back by g(t*), so normalising at one point and reporting another
// scales F by the ratio of g at the two.
//
// Without an interior maximum (Lemma 3d, gamma = 0, or the inflection point
// at alpha = 1 - gamma^2 / (8 beta); see mhn::kernel_interior_peak) the
// kernel falls from its value at the origin, unbounded for alpha < 1, and has
// no interior largest value to normalise by.  t* is then the peak of the
// density of log X, the positive root of 2 beta t^2 - gamma t - alpha = 0, or
// x if that is smaller.  It lies where the mass is and moves with the scale
// of the distribution, so g(t) / g(t*) stays within the range of a double
// over the whole quadrature range at every scale.
double compute_log_peak(double alpha, double beta, double gamma, double x,
                        double* anchor_out = nullptr) {
  const double peak = mhn::kernel_interior_peak(alpha, beta, gamma);
  const double t_star = std::min(
      x, (peak > 0.0) ? peak : mhn::positive_root(beta, gamma, alpha));
  if (anchor_out != nullptr) *anchor_out = t_star;
  return (alpha - 1.0) * std::log(t_star) - beta * t_star * t_star
         + gamma * t_star;
}

// log of int_0^inf exp(-2 a t - r^2 t^2) dt = (sqrt(pi) / (2 r)) erfcx(a / r),
// for a > 0 and r >= 0: the integral, from the end it starts at, of a kernel
// whose logarithm falls there with slope 2a and curvature 2 r^2.  Half the
// slope and the square root of half the curvature are passed, so that
// neither is squared.  erfcx(y) = exp(y^2) erfc(y) is taken from its
// continued fraction for y >= 2 and as exp(y^2) 2 Phi(-sqrt(2) y), on the log
// scale, below.  Where r is below 1e-16 a, y exceeds 1e16 and the integral is
// the pure exponential's 1 / (2a) to within 1/(2 y^2) < 5e-33 relative; the
// callers pass r = 0 for a curvature that is negative but as small as that
// next to the slope (see narrow_tail).
double log_quadratic_tail(double a, double r) {
  if (!(r >= 1e-16 * a)) return -M_LN2 - std::log(a);
  const double y = a / r;
  const double log_erfcx = (y >= 2.0)
      ? std::log(mhn::erfcx_cf(y))
      : y * y + M_LN2 + R::pnorm(-M_SQRT2 * y, 0.0, 1.0, 1, 1);
  return 0.5 * std::log(M_PI) - M_LN2 - std::log(r) + log_erfcx;
}

// log of the integral of the kernel g(t) = t^(alpha-1) exp(-beta t^2 + gamma t)
// from x away from its peak -- down to the origin for lower = true, out to
// infinity otherwise -- divided by g(x), where the kernel falls by far more
// than L within half a spacing of the doubles at x.  The stretch that
// carries the integral is then at most half a spacing long: t / x is below
// about 1.1e-16 below the peak, and 5.6e-14 past the log-axis peak, where the
// spacing is that of log x.  The cubic term of the log kernel,
// (alpha - 1) t^3 / (3 x^3), is then at most (alpha - 1) 6e-41 nats, below
// the rounding that log_kernel_resolvable admits wherever this form is used,
// and the integral is the quadratic tail of log_quadratic_tail
// with half the slope a = |L'(x)| / 2 and half the curvature
// r^2 = beta + (alpha - 1) / (2 x^2).  L'(x) = (alpha - 1)/x + gamma -
// 2 beta x is formed with gamma/2 - beta x by fma, in one rounding, since
// gamma and 2 beta x cancel near the peak.
//
// r^2 is positive for alpha >= 1, and r is formed as a hypot.  For alpha < 1
// it is negative below x = sqrt((1 - alpha) / (2 beta)).  There the
// curvature is too small to make the kernel narrow, so the callers reach
// this only where the slope alone carries the kernel down by L within half a
// spacing of the doubles -- in x for the lower tail, and in log x, whose
// spacing is at most 1.1e-13, for the upper one.  Then |L'(x)| x exceeds
// 1e15, the negative part, at most 1 / (2 x^2) in size, is below 2e-30 a^2,
// and over the stretch t < 60 / a that carries the integral it moves the
// exponent by less than 1e-26; r = 0, the pure exponential, is passed.
//
// a and r are rates per unit of x.  a overflows for an x below
// (alpha - 1) / DBL_MAX, and r for one below sqrt(alpha) / DBL_MAX, while the
// integral itself, about 1 / (2a), is still inside the range; there both are
// multiplied by x (with t = x tau), which leaves the dimensionless a x and
// r x, and log x is added back.  Returns NaN where the computed a is not
// positive, which the callers' branches rule out but for rounding; they then
// report the tail as empty, which it is to within 3e-10.
double narrow_tail(double alpha, double beta, double gamma, double x,
                   bool lower) {
  const double half_power = 0.5 * (alpha - 1.0);
  const double c_neg = (alpha < 1.0) ? std::sqrt(-half_power) : 0.0;
  auto root_for = [&](double b_x, double p_x) -> double {
    // sqrt(b_x^2 + p_x^2) for alpha >= 1, and sqrt(b_x^2 - p_x^2), or 0 where
    // that is not positive, for alpha < 1; b_x = sqrt(beta) times the unit,
    // p_x = sqrt(|alpha - 1| / 2) over it.
    if (alpha >= 1.0) return std::hypot(b_x, p_x);
    return (b_x > p_x) ? std::sqrt(b_x - p_x) * std::sqrt(b_x + p_x) : 0.0;
  };
  const double sb = std::sqrt(beta);
  const double sp = (alpha >= 1.0) ? std::sqrt(half_power) : c_neg;
  // Half the slope toward the far end: L'(x)/2 below the peak, -L'(x)/2
  // past it.
  double a = lower ? half_power / x + mhn::fma(-beta, x, 0.5 * gamma)
                   : mhn::fma(beta, x, -0.5 * gamma) - half_power / x;
  double r = root_for(sb, sp / x);
  double log_scale = 0.0;
  if (!std::isfinite(a) || !std::isfinite(r)) {
    a = lower ? half_power + x * mhn::fma(-beta, x, 0.5 * gamma)
              : x * mhn::fma(beta, x, -0.5 * gamma) - half_power;
    r = root_for(sb * x, sp);
    log_scale = std::log(x);
  }
  if (!(a > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  return log_scale + log_quadratic_tail(a, r);
}

// Whether the log probability at x, a tail beyond x of a law narrower than
// the spacing of the doubles, can be assembled from the log kernel.  The
// caller of the routines below forms it as the log kernel at the anchor they
// return plus log 2 + (alpha/2) log beta - log Psi, two quantities that
// nearly cancel.  Each is a sum of terms up to about alpha times the
// logarithms of x, of the peak m and of beta in size, and carries their
// rounding: 4 eps (|alpha - 1| + 1) times the sum of those logarithms bounds
// it.  That is below a nat up to a shape of about 1e12.  Past it, the
// probability itself is below exp(-D), D the fall of the log kernel from its
// peak to x, formed here relative to the peak as
// beta d^2 - (alpha - 1) (log(1 + d/m) - d/m), d = x - m, which nothing
// large enters.  Where the rounding stays below D - 19, the assembled log
// probability stays below -19, and the probability below 1e-8, even at the
// extreme of the rounding, and its logarithm is still right to within a
// relative rounding / D.  Otherwise the sum is noise that can put the
// probability anywhere in [0, 1] -- at alpha = 1e33 the rounding is some
// 1e19 nats -- and the tail is given as empty instead, which it is to within
// 3e-10.
bool log_kernel_resolvable(double alpha, double beta, double x, double m) {
  const double size = std::fabs(std::log(x))
                      + std::fabs(std::log(
                            std::max(m, std::numeric_limits<double>::min())))
                      + std::fabs(std::log(beta)) + 1.0;
  const double rounding = 4.0 * std::numeric_limits<double>::epsilon()
                          * (std::fabs(alpha - 1.0) + 1.0) * size;
  if (rounding <= 1.0) return true;
  if (!(m > 0.0)) return false;
  const double d = x - m;
  const double fall =
      beta * d * d - (alpha - 1.0) * mhn::log1p_minus_u(d / m);
  return rounding <= fall - 19.0;
}

// The integral as two pieces: the log of the scale the integrand is measured
// against (the kernel at the anchor, or the part next to the origin, see
// log_cdf_integrate_parts), and the log of the scaled quadrature.  They are
// kept apart because the first is about gamma^2/(4 beta) for a large positive
// tilt -- 2.5e13 at gamma = 1e7 -- and any caller that wants the value with
// that term removed must never form their sum.
struct CdfParts {
  double log_peak;
  double log_val;
  bool ok;
  double log_peak_shifted = 0.0;   // the same with gamma^2/(4 beta) removed
};

// The Gamma limit.  Under a tilt below z = gamma / sqrt(beta) = -1e150 the law
// sits near alpha / |gamma|, and wherever beta t^2 <= 2^-60 over the range
// that carries the integral, the factor exp(-beta t^2) is 1 to within that,
// so the integral of t^(alpha-1) exp(-beta t^2 + gamma t) is
// Gamma(alpha) |gamma|^(-alpha) times the regularised incomplete gamma
// function at |gamma| t, which R's pgamma gives in either tail and on the log
// scale.  This is needed, not only cheaper: past |z| of about 1e290 that
// scale approaches the bottom of the double range, where the abscissae of a
// quadrature in t are subnormal and hold only some of their digits (at
// (0.5, 1, -1e300) such a quadrature leaves the lower tail up to 19 percent
// low).  For the lower tail over [0, x] the condition is beta x^2 <= 2^-60.
// For the upper tail from x it is taken at
// t_c = (max(|gamma| x, 2 alpha) + 90) / |gamma|: past 2 alpha the
// Gamma(alpha, 1) log-density falls at least half as fast as its argument, so
// beyond |gamma| t_c lies at most 2 exp(-45) of the tail from x, and
// exp(-beta t^2) can only lower that.  (Returns false where the condition
// fails.)
bool gamma_limit_parts(double alpha, double beta, double gamma, double x,
                       bool lower, CdfParts* out) {
  if (!(gamma < -1e150 * std::sqrt(beta))) return false;
  const double rate = -gamma;
  const double y = x * rate;
  const double t_c = lower ? x : (std::max(y, 2.0 * alpha) + 90.0) / rate;
  if (!(beta * t_c * t_c <= 0x1p-60)) return false;
  const double log_scale = std::lgamma(alpha) - alpha * std::log(rate);
  const double log_val = R::pgamma(y, alpha, 1.0, lower ? 1 : 0, 1);
  *out = CdfParts{log_scale, log_val, std::isfinite(log_val), log_scale};
  return true;
}

CdfParts log_cdf_integrate_parts(double alpha, double beta, double gamma,
                                 double x, double tol) {
  if (!(x > 0.0)) return CdfParts{0.0, 0.0, false};
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();
  {
    CdfParts limit;
    if (gamma_limit_parts(alpha, beta, gamma, x, true, &limit)) return limit;
  }

  // Up to the point x0 of mhn::kernel_head, exp(-beta t^2 + gamma t) agrees
  // with 1 - beta t^2 + gamma t to within 1e-17, so for x <= x0 the whole
  // integral is the closed form there: x^alpha / alpha times a factor within
  // about 1e-8 of one.  Taken so it is exact however small x is, whereas a
  // quadrature range anchored at such an x is sized by a slope of about
  // (alpha - 1) / x, whose square overflows once x is below about 1e-154.
  // With gamma^2/(4 beta) removed, as the shifted form asks, the result is
  // formed directly: that form is only used for z > 1e3, where F(x) at such
  // an x is below exp(-z^2/4) and underflows whatever its rounding.
  {
    const double log_lead = alpha * std::log(x) - std::log(alpha);
    const mhn::KernelHead whole =
        mhn::kernel_head(alpha, beta, gamma, x, log_lead);
    if (whole.x0 >= x && whole.value > 0.0) {
      return CdfParts{log_lead, std::log(whole.value), true,
                      log_lead - gamma * (0.25 * gamma / beta)};
    }
  }

  double anchor_t = 0.0;
  const double log_peak = compute_log_peak(alpha, beta, gamma, x, &anchor_t);
  // The integrand is measured against exp(log_scale), the kernel at the
  // anchor, unless the part by the origin is larger by more than e^600.  For
  // alpha < 1 that part is int_0^x0 t^(alpha-1) dt = x0^alpha / alpha (times
  // a factor near 1, see mhn::kernel_head), about 1 / alpha for a tiny shape,
  // while the kernel at an anchor far from the origin can be small: at
  // alpha = 1e-300 with beta = 1e-300, an anchor x = 1e145 has
  // x^(alpha-1) = 1e-145, and the head measured against it was e^1024, past
  // the largest double, so that the integral was infinite and log F -Inf.
  // The head's own size is then the scale.  Both scales give the same
  // integral, each term divided by the same factor; outside that case the
  // kernel at the anchor stays the scale, and the values are as they were.
  // Under the head's scale the guard in scaled_kernel drops more terms, each
  // below e^-700 of the head, over a range that ends within
  // (max(z, 0) + sqrt(L)) / sqrt(beta) of the origin, at most e^387 for
  // z <= 1e6 and 1 / sqrt(beta) <= 2^537: together below e^-313 of the head,
  // which the integral holds.  The head
  // is cut at x here and at the range's anchor in kernel_head below.  The
  // density of log X rises up to its peak t*, so the kernel at the anchor
  // is at least x0^alpha / t*, and the switch needs t* / alpha > e^600, with
  // t* below e^386: alpha below e^-214 (1.2e-93), where the two cuts change
  // the head by a factor within 1e-89 of one.
  double log_scale = log_peak;
  if (alpha < 1.0) {
    const double log_head =
        alpha * std::log(std::min(mhn::kernel_head_x0(beta, gamma), x)) -
        std::log(alpha);
    if (log_head - log_peak > 600.0) log_scale = log_head;
  }
  const double scale_shift = log_scale - log_peak;

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
  // The interior maximum, on the same branches compute_log_peak uses (Sun et
  // al. 2023): for alpha > 1 always (Lemma 3b), at alpha = 1 for gamma > 0
  // (Lemma 6b), and for alpha < 1 only when a positive tilt creates one
  // (Lemma 3c).  The last two cases are needed as much as the first: for a
  // large tilt their peak lies within rounding of gamma/(2 beta), and the
  // shifted exponent formed at the end can take their offset from the
  // conjugate form 2(alpha-1)/(sqrt(D) + gamma), rather than from the
  // difference of the two, only if the peak is recognised here.
  // root_c = sqrt(D), D = gamma^2 + 8 beta (alpha - 1), formed without
  // squaring gamma where that would overflow.
  const double root_c = mhn::root_discriminant(beta, gamma, alpha - 1.0);
  const double peak_c = mhn::kernel_interior_peak(alpha, beta, gamma);
  const double centre = (peak_c > 0.0) ? std::min(peak_c, x) : 0.0;
  const bool centred = (centre > 0.0)
      && std::fabs(gamma) / std::sqrt(beta) > 1.0e2;
  // The linear coefficient of the expansion about `centre`.  At the peak the
  // stationarity condition 2 beta m - gamma = (alpha-1)/m gives it without
  // subtraction; elsewhere the subtraction is between quantities whose
  // difference is of order one, so it is harmless.
  double lin = (centre == peak_c && peak_c > 0.0)
      ? (alpha - 1.0) / centre
      : 2.0 * beta * centre - gamma;
  // 2 beta overflows for beta above about 9e307; beta c is then taken first.
  if (!std::isfinite(lin)) lin = 2.0 * (beta * centre) - gamma;

  auto scaled_kernel = [alpha, beta, gamma, log_scale, scale_shift, centre,
                        centred, lin](double t) -> double {
    if (t <= 0.0) return 0.0;
    double delta;
    if (centred) {
      // log g(c + s) - log g(c) = (alpha-1) log(1 + s/c) - beta s^2 - s (2 beta c - gamma).
      // Expanding about the peak the last term cancels against the first's
      // linear part, which is the classical form; away from it the term is what
      // makes the expansion valid at all, and without it F(mode - 1 sd) would
      // come out 0.303 where the normal limit it must approach gives 0.1587.
      const double sd = t - centre;
      const double u = sd / centre;
      if (u <= -1.0) return 0.0;
      delta = (alpha - 1.0) * std::log1p(u) - beta * sd * sd - sd * lin -
              scale_shift;
    } else {
      const double log_g = (alpha - 1.0) * std::log(t)
                           - beta * t * t + gamma * t;
      delta = log_g - log_scale;
    }
    if (delta < -700.0) return 0.0;  // exp underflow guard
    return std::exp(delta);
  };

  // The range is anchored on the point where the mass is, split there, and
  // trimmed on both sides to where the kernel has fallen by L from its value
  // at the anchor, so that each panel is monotone and a fixed-order rule is
  // never spread over a range much wider than the part of it that carries
  // mass (at gamma = 1e5 that part is about one unit wide and lies next to
  // x = 5e4).
  //
  // With an interior peak m the kernel increases up to m, so on [0, x] it is
  // largest at min(m, x), and that is the anchor.  Past a peak the reach is
  // that of mhn::kernel_decay_reach: the smaller of the Gaussian reach
  // sqrt(L / beta) and the distance over which the power-law factor alone
  // falls by L, about L m / (alpha - 1) for large L.  A strong negative tilt
  // squeezes the peak towards the origin, and the second is then thousands of
  // times the shorter: at (2.5, 1, -3e4) it is 2e-3 against 7.6.  Anchored at
  // x short of the peak, the reach is the one the slope there gives.
  //
  // Without an interior peak (see mhn::kernel_interior_peak) the kernel falls
  // from the origin, where the mass is.  Once beta t^2 - gamma t reaches L the
  // rest of the range is negligible: beyond that point exp(-beta t^2 +
  // gamma t) is below exp(-L) and t^(alpha-1) only falls, so the neglected
  // part is below alpha exp(1 - L) of the total.  The range stops there, or
  // at x if that is smaller; run out to an x far past the support, tanh_sinh
  // would spend its abscissae where the integrand has vanished.
  //
  // L is a drop from the value at the anchor, and both reaches count the
  // factor t^(alpha-1) themselves, so L needs no term in alpha; 40 nats past
  // -log(tol) keeps what is cut off far below the tolerance.  A margin that
  // grew with alpha would stretch the range with it: at alpha = 1e8 the mass
  // is about one unit wide next to x = 7071, and over a range reaching back
  // to the origin the 15-point rule would place no abscissa within it.
  const double L = -std::log(tol_eff) + 40.0;
  const bool has_peak = (peak_c > 0.0);
  const bool split = has_peak && (peak_c < x);
  const bool at_peak = has_peak && (peak_c <= x);
  const double anchor = has_peak
      ? std::min(peak_c, x)
      : std::min(x, mhn::positive_root(0.5 * beta, gamma, L));
  const double reach = at_peak
      ? mhn::kernel_decay_reach(alpha, beta, peak_c, L)
      : mhn::kernel_decay_reach_at(
            beta, mhn::kernel_decay_slope(alpha, beta, gamma, anchor), L);
  // The shifted exponent is formed at the anchor the result is normalised
  // on.  A caller recomputing that anchor would have to repeat
  // compute_log_peak's choice exactly -- for alpha < 1 with a positive tilt it
  // is the interior maximum of Lemma 3c, not x -- or F would be scaled by the
  // ratio of the kernel at the two points.
  const double d_shift = (anchor_t == peak_c && peak_c > 0.0 && gamma > 0.0)
      ? (alpha - 1.0) / (0.5 * root_c + 0.5 * gamma)
      : anchor_t - 0.5 * gamma / beta;
  const double log_peak_shifted =
      (alpha - 1.0) * std::log(anchor_t) - beta * d_shift * d_shift;

  // A kernel narrower than the spacing of the doubles at the anchor.  Past a
  // tilt of about 3e17 in scale units, or a shape of about 4e34, the kernel
  // falls by L within less than half a spacing of the doubles at its peak
  // (and sooner at an anchor x below it), so that anchor - reach rounds to
  // the anchor: every abscissa of a panel would round onto its ends, and the
  // panel has no width.  The integral is then taken in closed form.  Below
  // the peak it is the quadratic tail from x (narrow_tail).  At or past the
  // peak m it is the Laplace integral sqrt(2 pi / c) Phi((x - m) sqrt(c)),
  // with c = 2 beta + (alpha - 1) / m^2 the curvature at the maximum, its
  // root formed as a hypot for alpha >= 1: the normal factor is 1 once x lies
  // a spacing or more past m, and 1/2 at x = m.  m here is the double nearest
  // the peak, which the true peak can miss by up to half a spacing, many of
  // the law's standard deviations; so F(m) = 1/2 is a convention of this
  // form, not the law's value there, which mhn_cdf_peak.h resolves from the
  // exact offset of the peak.
  //
  // Below the peak, the kernel falls by at least L / 3 from the peak to x
  // (it falls by L within half a spacing of x, and x lies at least that far
  // from the peak), so F(x) is below Phi(-sqrt(2 L / 3)), 3e-10.  Where the
  // log probability cannot be assembled from the log kernel
  // (log_kernel_resolvable), that is the answer the caller is given, as an
  // empty integral: F = 0.
  if (has_peak && anchor - reach == anchor) {
    double log_val;
    if (at_peak) {
      const double root = (alpha >= 1.0)
          ? std::hypot(M_SQRT2 * std::sqrt(beta),
                       std::sqrt(alpha - 1.0) / peak_c)
          : std::sqrt(2.0 * beta + ((alpha - 1.0) / peak_c) / peak_c);
      log_val = 0.5 * std::log(2.0 * M_PI) - std::log(root)
          + ((x == peak_c) ? -M_LN2
                           : R::pnorm((x - peak_c) * root, 0.0, 1.0, 1, 1));
    } else {
      if (!log_kernel_resolvable(alpha, beta, anchor, peak_c)) {
        return CdfParts{0.0, 0.0, false};
      }
      log_val = narrow_tail(alpha, beta, gamma, anchor, /*lower=*/true);
      if (std::isnan(log_val)) return CdfParts{0.0, 0.0, false};
    }
    return CdfParts{log_peak, log_val, true, log_peak_shifted};
  }

  double lo = 0.0;
  if (has_peak && alpha >= 1.0) lo = std::max(0.0, anchor - reach);
  const double hi = split ? std::min(x, anchor + reach) : anchor;

  double val;
  if (alpha < 1.0) {
    // lo is 0 whenever alpha < 1, and the head of that range is beyond the
    // reach of any quadrature -- see mhn::kernel_head.  Left to tanh_sinh alone
    // it would be lost: F would come out 6.3e-4 at (1e-6, 1, -1000, q = 1e-3),
    // where the true value is 1.
    boost::math::quadrature::tanh_sinh<double> ts;
    const mhn::KernelHead head =
        mhn::kernel_head(alpha, beta, gamma, anchor, log_scale);
    val = head.value;
    // Below a Lemma 3c peak the kernel is not monotone: it falls from the
    // origin to a local minimum at the smaller root of
    // 2 beta t^2 - gamma t - (alpha - 1) = 0, t_min = 2 (1 - alpha) /
    // (gamma + sqrt(D)), and rises from there.  Each tanh_sinh range is kept
    // to one monotone piece.  The routine drops end abscissae whose
    // contribution is negligible, and if the integrand still increases
    // towards the end it dropped, it takes that for a singularity it has cut
    // off and raises an error, even though the part dropped is negligible.
    //
    // The rising piece is smooth, and for a large tilt its mass lies within a
    // few units of the anchor a while the piece itself spans [t_min, a], some
    // gamma / (2 beta) wide, too wide for tanh_sinh to find it reliably.  It
    // is trimmed like the alpha > 1 range.  With u = s / a <= 1/2,
    //   log g(a) - log g(a - s) = L'(a) s + beta s^2 - (1 - alpha)(-log(1 - u) - u)
    //                          >= L'(a) s + (beta - (1 - alpha) / a^2) s^2,
    // because -log(1 - u) - u <= u^2 there.  Where a^2 >= 2 (1 - alpha) / beta
    // the coefficient of s^2 is at least beta / 2, and since the kernel rises
    // on the piece, everything further below a than the distance at which
    // this bound reaches L is below exp(-L) g(a).
    double start = head.x0;
    if (has_peak) {
      const double t_min = (1.0 - alpha) / (0.5 * gamma + 0.5 * root_c);
      if (t_min > start && t_min < anchor) {
        val += ts.integrate(scaled_kernel, start, t_min, tol_eff);
        start = t_min;
      }
      if (t_min < anchor && beta * anchor * anchor >= 2.0 * (1.0 - alpha)) {
        double slope = (alpha - 1.0) / anchor + gamma - 2.0 * beta * anchor;
        if (!std::isfinite(slope)) {
          slope = (alpha - 1.0) / anchor + gamma - 2.0 * (beta * anchor);
        }
        const double reach_left = mhn::kernel_decay_reach_at(0.5 * beta, slope, L);
        if (reach_left <= 0.5 * anchor) {
          start = std::max(start, anchor - reach_left);
        }
      }
    }
    if (start < anchor) val += ts.integrate(scaled_kernel, start, anchor, tol_eff);
    // Past a Lemma 3c peak the kernel is smooth and log-concave: the second
    // derivative of its logarithm, (1 - alpha) / t^2 - 2 beta, is negative at
    // the peak and falls from there.  The Gauss-Kronrod rule of the
    // alpha >= 1 case then serves.  tanh-sinh does not: on a range many
    // orders of magnitude narrower than the kernel, such as [peak, x] for an
    // x past the peak by 1e-15 to 1e-11 of the peak's abscissa, its error
    // estimate need not settle, and it can refine to its limit, about 2 ms,
    // for a value that is the width of the range to double precision.
    if (hi > anchor) {
      val += mhn::gk_unit<15>(scaled_kernel, anchor, hi, tol_eff);
    }
  } else {
    // Both panels scale with 1 / sqrt(beta), and are mapped onto [-1, 1] so
    // that Boost's stopping rule sees the same integrand at every beta (see
    // mhn_quadrature.h).
    val = mhn::gk_unit<15>(scaled_kernel, lo, anchor, tol_eff);
    if (hi > anchor) {
      val += mhn::gk_unit<15>(scaled_kernel, anchor, hi, tol_eff);
    }
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return CdfParts{0.0, 0.0, false};
  }
  return CdfParts{log_scale, std::log(val), true,
                  log_peak_shifted + scale_shift};
}

}  // namespace

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

namespace {

CdfParts log_ccdf_integrate_parts(double alpha, double beta, double gamma,
                                  double x, double tol) {
  if (!(x > 0.0)) x = 0.0;
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();
  {
    CdfParts limit;
    if (gamma_limit_parts(alpha, beta, gamma, x, false, &limit)) return limit;
  }

  // Integrate on the log axis.  With t = exp(u) the Jacobian folds into the
  // power and
  //
  //     int_x^inf t^(alpha-1) e^{-beta t^2 + gamma t} dt
  //         = int_{log x}^{inf} exp(alpha u - beta e^{2u} + gamma e^u) du,
  //
  // an integrand that is bounded, smooth and free of the t^(alpha-1)
  // singularity.  In t-space the range would need a different rule for every
  // regime: anchored at the mode the slope is zero, so a range set by the
  // slope collapses to the Gaussian width and the panel past the mode sees
  // nothing; and for alpha < 1 with a tiny x the range spans sixty decades that
  // no fixed-order rule can cover.  On the log axis one rule serves every case.
  //
  // The stationary point solves alpha - 2 beta e^{2u} + gamma e^u = 0, that is
  // 2 beta v^2 - gamma v - alpha = 0 with v = e^u -- note alpha, not alpha - 1,
  // because the Jacobian shifted it.  The curvature there is
  // -(4 beta v^2 - gamma v) = -(2 beta v^2 + alpha), using the stationarity
  // condition to remove the subtraction.
  //
  // Under a positive tilt at a small beta the peak, about gamma / (2 beta),
  // can lie beyond the largest double, where the whole law does, and v_peak
  // is then +Inf; u_peak is still taken from its logarithm.
  const double v_peak = mhn::positive_root(beta, gamma, alpha);
  const double u_peak = (mhn::opaque_value(v_peak) <= mhn::kDoubleMax)
      ? std::log(v_peak)
      : mhn::log_positive_root(beta, gamma, alpha);
  const double curv = 2.0 * beta * v_peak * v_peak + alpha;
  const double L = -std::log(tol_eff) + 40.0;

  const double lo = (x > 0.0) ? std::log(x)
                              : -std::numeric_limits<double>::infinity();

  // The integrand is largest at whichever of the peak and the lower limit comes
  // last, since past the peak it only falls.  Everything -- the scaling, the
  // upper limit and the panel widths -- is measured from that point.  Scaling by
  // the peak when the range starts above it would make every value underflow:
  // at (2.5, 1, 1) with q = 30 the scaled integrand would be exp(-861)
  // throughout, the total zero and the log upper tail -Inf.
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
  //
  // The identity holds with 2 beta v - gamma for any v = e^{u_a}; alpha / v
  // equals it only at the true peak, and v_anchor is the peak passed through
  // its logarithm.  For a normal double that moves it by a relative 1e-13 at
  // most, but a subnormal peak is resolved only to its spacing, a relative
  // 1 / (2 k) for a peak k spacings above zero, and alpha / v then misses
  // 2 beta v - gamma by as much: the dip of the integrand moves along the
  // axis, and log Q was 3.5e-7 off at (1.41e-217, 7.65e208, -1.04e106),
  // q = 1.5e-323.  A subnormal peak lies under a strong negative tilt or a
  // small positive one, where 2 beta v - gamma = alpha / v cancels little:
  // it exceeds 1, alpha being a normal double, while a positive gamma is
  // below 2 beta v < 2 DBL_MAX DBL_MIN, about 8, so that at most 3 bits are
  // lost.  There the direct form is taken.
  double lin_u = (u_anchor == u_peak &&
                  v_anchor >= std::numeric_limits<double>::min())
      ? alpha / v_anchor
      : 2.0 * beta * v_anchor - gamma;
  // 2 beta overflows for beta above about 9e307; beta v is then taken first.
  if (!std::isfinite(lin_u)) lin_u = 2.0 * (beta * v_anchor) - gamma;
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
  //
  // Under a tilt beyond about 1e154 in scale units, beta v^2 and gamma v
  // overflow, and the slope formed above is Inf - Inf.  It is then formed as
  // alpha + v (gamma - 2 beta v), and at the peak, where the terms cancel
  // exactly, as 0.  The curvature and the slope's square can overflow as
  // well, and the reach is then written with everything halved,
  // L / (slope/2 + hypot(slope/2, sqrt(L/2) sqrt(curv))), with sqrt(curv)
  // itself a hypot of sqrt(2 beta) v and sqrt(alpha).  sqrt(curv), about
  // z / sqrt(2), can pass the largest double where sqrt(curv / 2) does not,
  // and the product is then formed as sqrt(L) sqrt(curv / 2).
  double slope =
      std::fabs(alpha - 2.0 * beta * v_anchor * v_anchor + gamma * v_anchor);
  if (!std::isfinite(slope)) {
    slope = (u_anchor == u_peak)
        ? 0.0
        : std::fabs(alpha + v_anchor * (gamma - 2.0 * (beta * v_anchor)));
  }
  const double reach_d2 = slope * slope + 2.0 * L * curv;
  double reach;
  if (mhn::opaque_value(reach_d2) <= mhn::kDoubleMax) {
    reach = 2.0 * L / (slope + std::sqrt(reach_d2));
  } else {
    const double root_curv =
        std::hypot(M_SQRT2 * std::sqrt(beta) * v_peak, std::sqrt(alpha));
    const double scaled_root = (mhn::opaque_value(root_curv) <= mhn::kDoubleMax)
        ? std::sqrt(0.5 * L) * root_curv
        : std::sqrt(L) * std::hypot(std::sqrt(beta) * v_peak,
                                    std::sqrt(0.5 * alpha));
    const double half_slope = 0.5 * slope;
    reach = L / (half_slope + std::hypot(half_slope, scaled_root));
  }
  // The shifted exponent is built here, where v_anchor is still in hand, so
  // that no caller has to recover log(val) by subtracting two large numbers --
  // which would put back the cancellation this exists to remove.
  const double d_shift = (u_anchor == u_peak && gamma > 0.0)
      ? alpha / (0.5 * mhn::root_discriminant(beta, gamma, alpha) + 0.5 * gamma)
      : v_anchor - 0.5 * gamma / beta;
  const double log_peak_shifted = alpha * u_anchor - beta * d_shift * d_shift;

  // A kernel narrower than the spacing of the doubles at the anchor, as in
  // the lower tail: u_anchor + reach rounds to u_anchor, and a panel would
  // have no width.  On the log axis this happens once the tilt passes about
  // 5e15 in scale units, a shape of about 6e30, or far out in the upper tail,
  // where the slope at the lower limit exceeds 2 L over the spacing of u.
  // Anchored at the peak (the range starts below it), the integral is the
  // Laplace integral sqrt(2 pi / curv) Phi((u_peak - lo) sqrt(curv)), with
  // sqrt(curv) formed as a hypot.  u_peak - lo is taken as
  // log1p((v_peak - x) / x), from the exact difference of the two doubles:
  // log x itself rounds onto u_peak for an x several of the law's standard
  // deviations away from the peak.  Anchored at x past the peak, it is the
  // quadratic tail of the kernel in x itself (narrow_tail), from x rather
  // than from exp(log x), which can differ from x by more than the kernel's
  // width; the kernel's value is then taken at x as well, its exponent as
  // (alpha - 1) log x + x (gamma - beta x), which stays finite where beta x^2
  // and gamma x overflow.  As below the peak in the lower tail, S(x) is then
  // below 3e-10, and where the log probability cannot be assembled from the
  // log kernel the caller is given an empty integral, S = 0.
  if (u_anchor + reach == u_anchor) {
    if (u_anchor == u_peak) {
      // sqrt(curv) can pass the largest double; it is then taken as
      // sqrt(2) times sqrt(curv / 2), on the log scale.  Where the peak
      // itself lies beyond the largest double, sqrt(beta) v_peak outweighs
      // sqrt(alpha / 2) by far, and the log of the root is formed from
      // u_peak; every finite x then lies many standard deviations below it.
      const double root_curv =
          std::hypot(M_SQRT2 * std::sqrt(beta) * v_peak, std::sqrt(alpha));
      const double depth = (x > 0.0)
          ? std::log1p((v_peak - x) / x)
          : std::numeric_limits<double>::infinity();
      double log_val;
      if (mhn::opaque_value(root_curv) <= mhn::kDoubleMax) {
        log_val = 0.5 * std::log(2.0 * M_PI) - std::log(root_curv)
            + R::pnorm(depth * root_curv, 0.0, 1.0, 1, 1);
      } else {
        const double half_root_curv =
            std::hypot(std::sqrt(beta) * v_peak, std::sqrt(0.5 * alpha));
        const double log_half_root_curv =
            (mhn::opaque_value(half_root_curv) <= mhn::kDoubleMax)
            ? std::log(half_root_curv)
            : 0.5 * std::log(beta) + u_peak;
        log_val = 0.5 * std::log(M_PI) - log_half_root_curv
            + R::pnorm(M_SQRT2 * (depth * half_root_curv), 0.0, 1.0, 1, 1);
      }
      return CdfParts{log_peak, log_val, true, log_peak_shifted};
    }
    if (!log_kernel_resolvable(alpha, beta, x,
                               mhn::kernel_interior_peak(alpha, beta, gamma))) {
      return CdfParts{0.0, 0.0, false};
    }
    const double log_val = narrow_tail(alpha, beta, gamma, x, /*lower=*/false);
    if (std::isnan(log_val)) return CdfParts{0.0, 0.0, false};
    const double log_x = std::log(x);
    const double d = x - 0.5 * gamma / beta;
    return CdfParts{(alpha - 1.0) * log_x + x * mhn::fma(-beta, x, gamma),
                    log_val, true, (alpha - 1.0) * log_x - beta * d * d};
  }

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
    // written so that neither exponential is formed at full size.  The
    // argument of expm1 is alpha (lo - U) itself.  As the difference of
    // alpha lo - log_peak and alpha U - log_peak it is lost to rounding once
    // alpha (U - lo) falls below the spacing of the doubles at log_peak, as
    // it does for a shape below about 1e-11, and this stretch, which then
    // holds most of the upper tail, would count for nothing.
    const double a_hi = alpha * U - log_peak;
    if (a_hi > -700.0) {
      val += std::exp(a_hi) * (-std::expm1(alpha * (lo - U))) / alpha;
    }
  }

  const double q_lo = std::max(lo, U);
  if (hi > q_lo) {
    // Panels doubling in width away from the anchor, in both directions.
    //
    // One panel per side is not enough when the range is wide in units of the
    // anchor's own scale: at (0.1, 1e-4, 100) the range is 49 wide and the peak
    // 1.4e-4, and a fixed-order rule laid over the whole of it returns half the
    // mass.  Doubling costs a number of panels logarithmic in that ratio, and
    // each sees the part of the integrand that varies on its own scale.
    // Every panel is mapped onto [-1, 1] (mhn_quadrature.h): near the peak of a
    // strongly tilted kernel they are far narrower than one unit, and on the
    // raw interval Boost's stopping rule could not be met there.
    const double centre = std::min(std::max(u_anchor, q_lo), hi);
    double edge = centre;
    // Start narrow and double, in both directions.  Starting at the full reach
    // would put the whole range in one panel whenever the reach is a generous
    // over-estimate, and a fixed-order rule would then see nothing: at
    // (1e-6, 1, -1, q = 0.01) the reach is about 4.6e3 while the mass lies within
    // five units of the anchor, and the upper tail would come out as zero.
    // One unit of u multiplies the argument by e, so the far field never
    // varies on a scale wider than that; nearer in, the reach itself is the
    // smaller number.  A reach that has underflowed to 0 gives a width of 0,
    // which could never carry either loop to its end; the rest of the range
    // is then taken as one panel, as the upward loop does.
    double width = std::min(reach, 1.0);
    while (edge < hi) {
      const double next = (width > 0.0) ? std::min(hi, edge + width) : hi;
      val += mhn::gk_unit<31>(kernel, edge, next, tol_eff);
      edge = next;
      width = (width > 0.0) ? width * 2.0 : hi - edge;
      if (!(width > 0.0)) break;
    }
    edge = centre;
    width = std::min(reach, 1.0);
    while (edge > q_lo) {
      const double next = (width > 0.0) ? std::max(q_lo, edge - width) : q_lo;
      val += mhn::gk_unit<31>(kernel, next, edge, tol_eff);
      edge = next;
      width *= 2.0;
    }
  }

  if (!(val > 0.0) || !std::isfinite(val)) {
    return CdfParts{0.0, 0.0, false};
  }
  return CdfParts{log_peak, std::log(val), true, log_peak_shifted};
}

}  // namespace

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
