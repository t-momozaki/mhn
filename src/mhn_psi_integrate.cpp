// Numerical integration of the Fox-Wright Psi function, for either sign of
// gamma.  With u = sqrt(beta) * x and z = gamma / sqrt(beta), the
// normalisation of the density (Sun et al. 2023, Lemma 1a) gives
//
//   Psi[alpha/2, z] = 2 * integral_0^Inf u^(alpha-1) exp(-u^2 + z*u) du.
//
// The dispatcher (mhn_psi_dispatch.cpp) uses psi_integrate for gamma < 0, for
// z > 20 and for alpha > 100, and psi_integrate_shifted (below) serves the
// density and the distribution function at a large positive tilt.
//
// The integral is taken in one of two variables.  Under a positive tilt
// psi_integrate takes it in u, on panels of at most sqrt(L0) = 7.6 units on
// either side of the peak, L0 = -log(tol) + 40 (psi_integrate_shifted
// integrates in the offset from the peak).  It is taken in v = log u
// instead (psi_log_axis below), where the origin and the power u^(alpha-1)
// drop out of the integrand: for z <= 0; for z > 0 at a shape below 6 while
// the peak lies within sqrt(L0) of the origin; and for z beyond 1e3 or a
// shape beyond 1e6.  In u the terms alpha |log u|, u^2 and z u of the
// exponent (about z^2/4 and z^2/2 at the peak) round by more than a panel's
// tolerance from z of a few thousand and from alpha near 1e7, and the
// adaptive rule then subdivides to its depth limit, at about a thousand
// times the cost, where the integral in log u forms no such terms.  In
// both, the integrand is divided by its value at the peak, split there, and
// cut off where it has fallen far enough that what is left out is negligible
// at the requested tolerance; the panels are integrated with Boost.Math
// quadrature rules.

#include "mhn_fma.h"
#include "mhn_psi.h"
#include "mhn_quadrature.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// Peak of u^(alpha-1) exp(-u^2 + z u) for z > 0, the larger solution of
// 2u^2 - z u - (alpha-1) = 0, which moves out to about z/2 as z grows.  For
// alpha < 1 the kernel is unbounded at the origin, but a positive tilt still
// produces an interior maximum once z^2 > 8(1 - alpha) (Sun et al. 2023,
// Lemma 3c), and that is where the mass is; at alpha = 1 it is bounded and
// every positive tilt puts the maximum at z/2 (Lemma 6b).
// kernel_interior_peak gives both roots without squaring z where z^2 would
// overflow.
double integrand_mode(double alpha, double z) {
  const double m = mhn::kernel_interior_peak(alpha, 1.0, z);
  return std::max(m, std::numeric_limits<double>::epsilon());
}

// alpha log u - u^2 + z u, for a positive normal double u.
//
// Where the peak value of the integrand is modest next to alpha the three
// terms are much larger than their sum, and log u carries a rounding error of
// up to half a unit in its last place, which alpha multiplies: at alpha = 140,
// z = -42.7 the sum is 16.7, the terms 148, -8 and -123, and that alone can
// cost 1.6e-14.  So log u is taken as k log 2 + log1p(m - 1), from u = m 2^k
// with m in [1/sqrt(2), sqrt(2)), where m - 1 is exact and |log1p(m - 1)| is at
// most 0.35; log 2 is split into a leading part whose low 21 bits are zero,
// so that k times it is exact, and the rest.  Each product and sum is then
// formed together with its rounding error (the products by fma, the sums by
// two_sum), and the errors are added back at the end.  What remains is the
// rounding of log1p, at most about 6e-17 alpha, and a unit or two in the last
// place of the sum.
double peak_exponent(double alpha, double u, double z) {
  const double ln2_hi = 6.93147180369123816490e-01;
  const double ln2_lo = 1.90821492927058770002e-10;   // log 2 - ln2_hi
  int k;
  double m = std::frexp(u, &k);   // u = m 2^k, m in [1/2, 1)
  if (m < 0.70710678118654752440) {
    m *= 2.0;
    --k;
  }
  const double log_hi = k * ln2_hi;   // exact
  const double log_lo = k * ln2_lo + std::log1p(m - 1.0);
  double err = 0.0;
  const double t0 = alpha * log_hi;
  err += mhn::fma(alpha, log_hi, -t0);
  const double t1 = alpha * log_lo;
  err += mhn::fma(alpha, log_lo, -t1);
  const double t2 = -(u * u);
  err -= mhn::fma(u, u, t2);
  const double t3 = z * u;
  err += mhn::fma(z, u, -t3);
  double e;
  double s = mhn::two_sum(t0, t1, e);
  err += e;
  s = mhn::two_sum(s, t2, e);
  err += e;
  s = mhn::two_sum(s, t3, e);
  err += e;
  // An infinite sum (the value overflows, for alpha beyond about 5e305) has no
  // rounding error to add back, and err would be an infinity of either sign.
  if (!std::isfinite(s)) return s;
  return s + err;
}

// log Psi[alpha/2, z] by quadrature in v = log u.
//
// In v the integrand is u^alpha exp(-u^2 + z u), u = e^v, whose logarithm has
// derivative alpha + u (z - 2u): for either sign of z it has a single peak, at
// the positive root u* of 2u^2 - z u - alpha = 0.  Writing u = u* e^w and
// using that equation to eliminate z,
//
//   Psi[alpha/2, z] = 2 u*^alpha exp(-u*^2 + z u*)
//                     * int exp(-alpha (e^w - 1 - w) - u*^2 (e^w - 1)^2) dw,
//
// over the whole line.  Both terms of the exponent are <= 0 and vanish at
// w = 0, so the integrand is 1 at the peak and is formed without subtracting
// anything large, however strong the tilt: z enters only through u*, and
// u* (e^w - 1) = u - u* is the offset on the u axis.  Under a strong negative
// tilt u* is about alpha / |z|, and the integrand tends to
// exp(-alpha (e^w - 1 - w)), up to a constant the density of
// w = log(G / alpha) for G = |z| u ~ Gamma(alpha): the Gamma limit of the
// law, on the logarithmic scale.  The integrand is also analytic, where in u
// the power u^(alpha-1) has unbounded derivatives at the origin for
// alpha < 2 and, for alpha just above 1, bends within about (alpha - 1) / |z|
// of it.
//
// The curvature of the exponent at the peak is c = alpha + 2 u*^2, so the
// integrand is a bell of width 1 / sqrt(c) there.  To the right of the peak
// the exponent falls by at least c w^2 / 2 (e^w - 1 - w >= w^2 / 2 and
// (e^w - 1)^2 >= w^2 for w >= 0) and by at least alpha (e^w - 1 - w).  The
// first reaches L = -log(tol) + 40 at sqrt(2 L / c); the second where
// e^w = A + w, A = 1 + L / alpha, at a root r = log(A + r) that is at most
// 1 + log A (e^(1 + log A) = e A >= A + 1 + log A), and so also at most
// log(A + 1 + log A).  d_right, the smaller of sqrt(2 L / c) and 1 + log A,
// ends the panel on the right, and the exponent falls faster beyond it.  To
// the left, at w = -d, it has fallen by
//
//   F(d) = alpha (d - 1 + e^-d) + u*^2 (1 - e^-d)^2,
//
// which grows only like alpha d once e^-d is small: that is the mass
// u^(alpha-1) keeps near the origin, all of it as alpha -> 0.  Its rate of
// growth, F'(d) = (1 - e^-d) (alpha + 2 u*^2 e^-d), is at least
// alpha (1 - e^-d), so everything to the left of -d is at most
// exp(-F(d)) / (alpha (1 - e^-d)).  The left side is therefore taken on
// panels [-2d, -d], outwards from d = sqrt(2 L / c), until that bound falls
// below exp(-L) of the sum, or until the head below u = x0 is reached, which
// is taken in closed form (below).
//
// The two panels next to the peak are then at most sqrt(2 L) = 10.8 bell
// widths long, about as many as those of the integration in u under a
// positive tilt (psi_integrate_core; at most 15.2 there): on a range that
// holds a narrow peak at one end, an adaptive rule can accept an estimate
// whose two embedded rules agree by chance while both are off.
//
// The same can happen on the right, where d_right, set by lower bounds on
// the fall, can lie well past the point where the integrand has died.  When
// alpha (e^w - 1 - w) is the larger term of the fall at d_right, as for a
// small shape or under a strong negative tilt, where u* is about
// alpha / |z|, the integrand stays within a factor e of 1 until e^w is
// about 1 + 1 / alpha, and falls to e^-L within a few units of w, a cliff
// that a single rule over the whole panel resolves poorly.  Taken whole at
// the tolerance of the other panels, that panel is off by 4e-13 of itself
// at (alpha, z) = (0.16, -1993), about 400 units in the last place of
// log Psi (8.8e-14), and log Psi by up to 4000 eps (1 + |log Psi|) on a
// scan of 4.8e4 points with alpha from 1e-4 to 1.  There the panel ends
// instead at log(A + 1 + log A), where that term alone has passed L, and
// is taken in two halves.  When u*^2 (e^w - 1)^2 is the larger term, it
// alone reaches L at w = log(1 + sqrt(L) / u*), which can lie well short of
// d_right: at (15.9, -1.69) d_right is 2.05, and the integrand is below
// e^-L from w of about 1.3 on.  Taken whole, the panel leaves log Psi off
// by 170 eps (1 + |log Psi|) there.  It ends at that point instead, where
// the fall is at least L since alpha (e^w - 1 - w) >= 0, and is taken in
// two halves; ended there but taken whole, it misses by up to 34 eps
// (1 + |log Psi|), at (0.058, -1.65).  Under a negative tilt the integral
// then takes fewer evaluations, 163 per call on average against 228 where
// u*^2 (e^w - 1)^2 is the larger term; above z = 1e3 it takes 31 more, 248
// against 217.  Against the power series of Psi and its asymptotic
// expansion in 1 / z^2, summed in arbitrary precision, log Psi is then
// within 4.4 eps (1 + |log Psi|) at 1.27e5 points with alpha from 1e-4 to
// 2000 and z from -1e5 to -1e-4, and within 3.2 at all but one of them, at
// (73.7, -23.4).
//
// With `relative` set, the logarithm of the integral over w alone is
// returned, log Psi - log 2 - (alpha log u* - u*^2 + z u*), a quantity of
// order one or log(1 / sqrt(c)) that is formed without the peak exponent
// (see psi_log_axis_relative).
double psi_log_axis(double alpha, double z, double tol,
                    bool relative = false) {
  const double inf = std::numeric_limits<double>::infinity();
  if (std::isinf(z)) {
    if (!relative) return (z > 0.0) ? inf : -inf;
    // Under an infinite negative tilt u* = 0 and the integrand is
    // exp(-alpha (e^w - 1 - w)), whose integral is Gamma(alpha) e^alpha /
    // alpha^alpha (substitute G = alpha e^w).  From alpha = 100 on its
    // logarithm is written with Stirling's series, whose leading terms cancel
    // against alpha log alpha - alpha, truncated after the alpha^-5 term,
    // which leaves an error below 1/(1680 alpha^7) < 1e-17.
    if (z > 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (alpha < 100.0) {
      return std::lgamma(alpha) + alpha - alpha * std::log(alpha);
    }
    const double r = 1.0 / alpha;
    const double r2 = r * r;
    const double stirling =
        r * (1.0 / 12.0 - r2 * (1.0 / 360.0 - r2 * (1.0 / 1260.0)));
    return 0.5 * std::log(2.0 * M_PI) - 0.5 * std::log(alpha) + stirling;
  }
  const double L = -std::log(tol) + 40.0;

  // u* = (z + sqrt(z^2 + 8 alpha)) / 4, for z <= 0 written as
  // alpha / ((|z| + sqrt(z^2 + 8 alpha)) / 2), which only adds; hypot keeps
  // z^2 from overflowing.
  const double root = std::hypot(z, 2.0 * M_SQRT2 * std::sqrt(alpha));
  const double den = 0.5 * std::fabs(z) + 0.5 * root;
  const double u_star = (z > 0.0) ? 0.25 * z + 0.25 * root : alpha / den;
  const double u2 = u_star * u_star;
  const double c = alpha + 2.0 * u2;

  // At the peak 2 u*^2 - z u* = alpha, so the exponent there is
  // alpha (log u* - 1) + u*^2, and c is the curvature (above).  Once c
  // overflows, that exponent -- and log Psi, smaller by only about
  // log(c) / 2 -- is within a factor of four of the largest double or beyond
  // it, and log Psi is returned as Inf.  Under a positive tilt this happens
  // from z of about 1.9e154; from 2.7e154 u*^2 overflows as well, and the
  // three terms below would form Inf - Inf.  The integral over w alone is then
  // the Laplace integral sqrt(2 pi / c): its skewness corrections are of
  // relative order 1 / c, below 1e-308 (the third and fourth derivatives of
  // the exponent in w at the peak, alpha + 6 u*^2 and alpha + 14 u*^2, are
  // at most 7 c).
  if (mhn::opaque_value(c) > std::numeric_limits<double>::max()) {
    if (!relative) return inf;
    return 0.5 * std::log(M_PI) - std::log(u_star)
           - 0.5 * std::log1p(0.5 * (alpha / u_star) / u_star);
  }

  // The logarithm of the integrand at its peak, alpha log u* - u*^2 + z u*
  // (peak_exponent).  Its three terms are taken at the same rounded u*: at the
  // peak the sum is stationary in u*, so the rounding of u* itself drops out,
  // whereas forming each term from its own expression for u* would leave
  // about eps alpha of it.  Only when u* is subnormal or has underflowed, for a
  // small shape under a strong negative tilt, are log u* and z u* formed
  // without it.  The exponent is then about alpha (log u* - 1) with |log u*|
  // above 708, and the rounding of log u* costs it about eps alpha |log u*|,
  // a unit or so in its last place, as the compensated sum leaves elsewhere.
  double log_u_star, z_u_star, log_peak;
  if (u_star >= std::numeric_limits<double>::min()) {
    log_u_star = std::log(u_star);
    z_u_star = z * u_star;
    log_peak = peak_exponent(alpha, u_star, z);
  } else {
    log_u_star = std::log(alpha) - std::log(den);
    z_u_star = -alpha * (std::fabs(z) / den);
    log_peak = alpha * log_u_star - u2 + z_u_star;
  }
  // Beyond alpha of about 5e305 the peak exponent, and log Psi with it, can
  // overflow; the integral over w does not need it.
  if (!relative && !std::isfinite(log_peak)) return log_peak;

  auto scaled = [alpha, u_star](double w) -> double {
    const double s = u_star * std::expm1(w);   // u - u*
    const double fall = alpha * mhn::expm1_minus_id(w) + s * s;
    // Also 0 where u* has underflowed and e^w overflowed (s = 0 * Inf), far
    // to the right, where alpha (e^w - 1 - w) is already past L.
    return (fall < 745.0) ? std::exp(-fall) : 0.0;
  };
  const double panel_tol = tol / 10.0;

  // Right of the peak: the panel ends where the larger term of the fall at
  // d_right alone has passed L, and is taken in two halves (see above).
  // log A = log(1 + L / alpha) is written so that L / alpha cannot
  // overflow, and log(A + 1 + log A) as log A + log1p((1 + log A) / A), so
  // that A need not be formed.
  const double log_A = std::log(L + alpha) - std::log(alpha);
  const double d_right = std::min(std::sqrt(2.0 * L / c), 1.0 + log_A);
  const double em1_right = std::expm1(d_right);
  double val;
  if (u2 * em1_right * em1_right < alpha * mhn::expm1_minus_id(d_right)) {
    const double d_cliff = std::min(
        d_right, log_A + std::log1p((1.0 + log_A) * std::exp(-log_A)));
    val = mhn::gk_unit<31>(scaled, 0.0, 0.5 * d_cliff, panel_tol) +
          mhn::gk_unit<31>(scaled, 0.5 * d_cliff, d_cliff, panel_tol);
  } else {
    const double d_bell =
        std::min(d_right, std::log1p(std::sqrt(L) / u_star));
    val = mhn::gk_unit<31>(scaled, 0.0, 0.5 * d_bell, panel_tol) +
          mhn::gk_unit<31>(scaled, 0.5 * d_bell, d_bell, panel_tol);
  }

  // The head.  Below x0, the positive root of x^2 + |z| x = 4.5e-9
  // (mhn::kHeadBudget), exp(-u^2 + z u) agrees with 1 + z u - u^2 to within
  // 1e-17 (the first omitted term is (u^2 + |z| u)^2 / 2 at most), and term
  // by term
  //
  //   int_0^x0 u^(alpha-1) (1 + z u - u^2) du
  //     = (x0^alpha / alpha) (1 + z x0 alpha / (alpha+1)
  //                             - x0^2 alpha / (alpha+2)),
  //
  // as in mhn::kernel_head.  With x0 = u* e^w0 its logarithm less that of the
  // integrand's peak is alpha w0 - log(alpha) + u*^2 - z u* + log(1 + ...).
  // When x0 >= u* the head is taken up to the peak itself (w0 = 0).
  const double budget = mhn::kHeadBudget;
  const double x0_den =
      0.5 * std::fabs(z) + 0.5 * std::hypot(z, 2.0 * std::sqrt(budget));
  const double log_x0 = std::log(budget) - std::log(x0_den);
  const double w0 = std::min(log_x0 - log_u_star, 0.0);

  // Left of the peak: [-d, 0], then [-2d, -d], [-4d, -2d], ... down to w0.
  double d = std::min(std::sqrt(2.0 * L / c), -w0);
  val += mhn::gk_unit<31>(scaled, -d, 0.0, panel_tol);
  while (d > 0.0 && d < -w0) {
    const double one_minus_e = -std::expm1(-d);   // 1 - e^-d
    const double s = u_star * one_minus_e;
    const double fall = alpha * mhn::expm1_minus_id(-d) + s * s;
    if (-fall - std::log(alpha * one_minus_e) < std::log(val) - L) break;
    const double d_next = std::min(2.0 * d, -w0);
    val += mhn::gk_unit<31>(scaled, -d_next, -d, panel_tol);
    d = d_next;
  }
  if (!(val > 0.0) || !std::isfinite(val)) {
    return relative ? std::numeric_limits<double>::quiet_NaN() : -inf;
  }

  double log_val = std::log(val);
  if (d >= -w0) {
    const double t = std::exp(w0);   // x0 / u*
    const double shape = 1.0 + z_u_star * t * alpha / (alpha + 1.0)
                             - u2 * t * t * alpha / (alpha + 2.0);
    const double log_head = alpha * w0 - std::log(alpha) + u2 - z_u_star
                            + std::log(shape);
    // log(val + exp(log_head)); for alpha -> 0 the head, of order 1 / alpha,
    // is most of the integral.
    log_val = std::max(log_val, log_head)
              + std::log1p(std::exp(-std::fabs(log_val - log_head)));
  }
  if (relative) return log_val;
  return std::log(2.0) + log_peak + log_val;
}

// log Psi[alpha/2, z], less z^2/4 when `shifted`, where z = gamma / sqrt(beta)
// itself overflows: |gamma| exceeds the largest double times sqrt(beta),
// which needs beta < 1.  log |z| = log |gamma| - log(beta) / 2 is still in
// range, and the law is then at its limit to double precision.
//
// Under a negative tilt, with w = |z| u,
//     Psi = 2 |z|^(-alpha) int_0^inf w^(alpha-1) e^(-w) e^(-w^2 / z^2) dw
//         = 2 Gamma(alpha) |z|^(-alpha) E[exp(-W^2 / z^2)],
// W ~ Gamma(alpha, 1).  The cumulants of W^2 are E W^2 = alpha (alpha + 1),
// Var W^2 = alpha (alpha + 1) (4 alpha + 6) and so on, each higher one larger
// by a factor of order alpha, so log E[exp(-W^2 / z^2)] is
// -alpha (alpha + 1) / z^2 to within a relative (2 alpha + 3) / z^2, below
// 1.2e-308 for every alpha.  The correction is formed as the product of
// alpha sqrt(beta) / |gamma| and (alpha + 1) sqrt(beta) / |gamma|, each at
// most about 1, so that nothing is squared.
//
// Under a positive tilt the kernel is normal about its peak m, with
// curvature 2 + (alpha - 1) / m^2 = 2 to within 2e-308 and a cubic term
// (alpha - 1) s^3 / (3 m^3) that is as small, so
//     log Psi - z^2/4 = log 2 + log(pi)/2 + (alpha - 1) log m - (m - z/2)^2
// to double precision.  With t = (alpha - 1) / z = (alpha - 1) sqrt(beta) /
// gamma, m - z/2 = 2 (alpha - 1) / (sqrt(z^2 + 8 (alpha - 1)) + z) is t to
// within a relative 2 |t| / z, and log m = log(z/2) + 2 t / z to within
// (2 t / z)^2, so the sum is
//     log 2 + log(pi)/2 + (alpha - 1) log(z/2) + t^2.
// Unshifted, log Psi exceeds z^2/4 > 8e615 and is +Inf.
double psi_infinite_tilt(double alpha, double beta, double gamma,
                         bool shifted) {
  const double log_abs_z = std::log(std::fabs(gamma)) - 0.5 * std::log(beta);
  const double sb = std::sqrt(beta);
  if (gamma < 0.0) {
    const double r0 = (alpha * sb) / -gamma;
    const double r1 = ((alpha + 1.0) * sb) / -gamma;
    // log Gamma(alpha) and alpha log|z| each overflow past alpha of about
    // 1e305 while their difference, about -alpha (1 + log(|z| / alpha)),
    // need not; Stirling's series, exact to double precision there, is then
    // written with the two logarithms together.
    double lead = std::lgamma(alpha) - alpha * log_abs_z;
    if (!std::isfinite(lead)) {
      lead = alpha * (std::log(alpha) - 1.0 - log_abs_z)
             - 0.5 * std::log(alpha) + 0.5 * std::log(2.0 * M_PI);
    }
    return M_LN2 + lead - r0 * r1;
  }
  if (!shifted) return std::numeric_limits<double>::infinity();
  const double t = ((alpha - 1.0) * sb) / gamma;
  return M_LN2 + 0.5 * std::log(M_PI) + (alpha - 1.0) * (log_abs_z - M_LN2)
         + t * t;
}

// Shared core of the integration in u under a positive tilt.  Returns
// log Psi minus `shift`, where `shift` is either 0 or z^2/4.  The shifted
// form removes z^2/4 from the peak exponent analytically and integrates in
// the offset from the peak.  The unshifted form integrates in u itself, and
// hands over to psi_log_axis for z <= 0, for z beyond 1e3, for a shape
// beyond 1e6, for alpha < 6 when the peak lies within sqrt(L0) of the
// origin, and at a shape so large that its panels overflow.
double psi_integrate_core(double alpha, double beta, double gamma, double tol,
                          bool shifted) {
  // Psi[alpha/2, z] = 2 int_0^inf u^(alpha-1) exp(-u^2 + z u) du with
  // z = gamma / sqrt(beta), so the quadrature depends on the tilt only through
  // z.
  const double z = gamma / std::sqrt(beta);
  if (std::isinf(z)) return psi_infinite_tilt(alpha, beta, gamma, shifted);

  // For z <= 0, and under a strong positive tilt or at a large shape, the
  // integral is taken in log u (psi_log_axis).  In u the log integrand is
  // formed from alpha log u, u^2 and z u, about z^2/4 and z^2/2 near the
  // peak, and their rounding reaches the integrand through exp.  It exceeds
  // the tolerance asked of a panel from z of a few thousand and from alpha
  // near 1e7, and the adaptive rule would then subdivide to its depth limit,
  // costing 15 ms at (2.5, 1, 1e6) and 0.8 ms at (1e7, 1, 5) against a few
  // microseconds in log u.  Past a nat,
  // at z of about 1.3e8, the panels cannot resolve the bell at all, and
  // further out the exponent passes the range of exp, where tanh_sinh (which
  // takes a panel for alpha < 1) stops at the first infinite value.  In
  // log u the exponent is formed relative to the peak, and no such terms
  // arise.
  if (!shifted && (!(z > 0.0) || z > 1.0e3 || alpha > 1.0e6)) {
    return psi_log_axis(alpha, z, tol);
  }

  // Past the peak m_u the log integrand falls, over a distance s, by
  //     s^2 + (alpha-1) [ s/m_u - log(1 + s/m_u) ],
  // z having cancelled through the stationarity condition
  // 2 m_u - z = (alpha-1) / m_u.  For alpha >= 1 the bracketed term is
  // non-negative, so the fall is at least s^2 whatever the shape, and the
  // upper limit M_u = m_u + sqrt(L0), L0 = -log(tol) + 40, leaves out a tail
  // below exp(-L0) of the peak.  For alpha < 1 the bracket gives some of that
  // back, less than a third of a nat once the peak lies beyond sqrt(L0) (see
  // below).
  const double m_u = integrand_mode(alpha, z);
  const double M_u = m_u + std::sqrt(-std::log(tol) + 40.0);

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
    // D = sqrt(z^2 + 8 (alpha-1)), without squaring z where that overflows;
    // the sum D + z is halved first for the same reason.
    const double D = mhn::root_discriminant(1.0, z, alpha - 1.0);
    double offset;
    if (z > 0.0 && D > 0.0) {
      offset = (alpha - 1.0) / (0.5 * D + 0.5 * z);
      // Take the offset of the double m_u itself where the two differ by
      // more than rounding at the size of the offset: when integrand_mode did
      // not take the root branch, and when the spacing of the doubles at m_u
      // is larger than that, as it is for a large z.
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
  // the range reaches min(sqrt(L) + 2, m) from it: a few units for a shape of
  // order one, so the origin -- and for alpha < 1 the u^(alpha-1) singularity
  // there -- is far outside it and needs no separate treatment, and up to m
  // itself for a large shape, where the integrand vanishes at the origin.
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
    // A tenth of the caller's tolerance per panel, as for the unshifted path
    // below, on panels mapped onto [-1, 1] (see mhn_quadrature.h).
    const double panel_tol = tol / 10.0;
    const double val_c = mhn::gk_unit<31>(centred, -reach, 0.0, panel_tol)
                       + mhn::gk_unit<31>(centred, 0.0, reach, panel_tol);
    if (!(val_c > 0.0) || !std::isfinite(val_c)) {
      return -std::numeric_limits<double>::infinity();
    }
    return std::log(2.0) + log_peak_shifted + std::log(val_c);
  }

  // Below x0 the integral is taken in closed form -- see mhn::kernel_head;
  // quadrature alone cannot reach the head of the range as alpha approaches
  // zero.
  const mhn::KernelHead head = mhn::kernel_head(alpha, 1.0, z, m_u, log_peak);

  // A positive tilt carries the peak out to about z/2, and once it lies
  // beyond sqrt(L0) the integrand is a bell about one unit wide there: for
  // alpha >= 1 the log integrand falls by at least s^2 over a distance s on
  // either side of the peak; for alpha < 1 the factor exp(-u^2 + z u) alone
  // does so on the left, since z - 2 m_u = (1 - alpha) / m_u > 0 (on the
  // right, which tanh_sinh takes for alpha < 1, the fall is short of s^2 by
  // (1 - alpha) [s/m_u - log(1 + s/m_u)], less than a third of a nat).  Its
  // curvature at the peak, 2 + (alpha-1) / m_u^2, lies between about 2 and 4.
  // So the range below m_u - sqrt(L0) holds at most exp(-L0) m_u log(m_u / x0)
  // of the peak value and is left out, and on each side of the peak a panel of
  // at most sqrt(L0) = 7.6 units remains.  The 31-point Kronrod rule
  // integrates such a panel to about 3e-14 without subdividing it.  That is
  // what makes the result safe: an adaptive rule accepts an estimate as soon
  // as its two embedded rules agree, and on a panel that holds a narrow bell
  // at one end they can agree by chance while both are off -- by 3.6e-4 in
  // log Psi for tanh_sinh over [x0, m_u] at alpha = 71.4, z = 189, and by
  // 1e-4 in log Psi (2e-4 in that panel) for the 15-point pair over
  // [m_u, m_u + 15] at alpha = 82.8, z = 965.  For a peak nearer the origin
  // than sqrt(L0) the power u^(alpha-1) at the origin and the bell share the
  // range [x0, m_u], and below alpha = 6 its unbounded derivative of order
  // ceil(alpha - 1) there costs the rule hundreds to thousands of units in
  // the last place (2679 at alpha = 3.77, z = 6.1e-5); the integral is then
  // taken in log u instead.  From alpha = 6 on the loss is within two units.
  const double left_cut = 2.0 * m_u - M_u;   // m_u - sqrt(L0)
  if (!(alpha >= 6.0 || left_cut > head.x0)) {
    return psi_log_axis(alpha, z, tol);
  }

  auto scaled_fn = [alpha, z, log_peak](double u) -> double {
    if (u <= 0.0) return 0.0;
    const double log_integrand = (alpha - 1.0) * std::log(u) - u * u + z * u;
    return std::exp(log_integrand - log_peak);
  };

  // Split at the peak.  The scaled integrand is 1 at m_u and falls away on
  // both sides, so each panel is monotone.  Each panel is asked for a tenth of
  // the caller's tolerance, so that their sum still meets it.
  const double panel_tol = tol / 10.0;
  double val = head.value;
  val += mhn::gk_unit<31>(scaled_fn, std::max(head.x0, left_cut), m_u,
                          panel_tol);

  // Right panel, [m_u, M_u], the right half of the bell, mapped onto [-1, 1]
  // before Boost's stopping rule sees it (see mhn_quadrature.h).  For
  // alpha < 1 tanh_sinh takes it, its nodes crowding towards the peak at the
  // left end.
  if (alpha < 1.0) {
    boost::math::quadrature::tanh_sinh<double> ts;
    val += ts.integrate(scaled_fn, m_u, M_u, panel_tol);
  } else {
    val += mhn::gk_unit<31>(scaled_fn, m_u, M_u, panel_tol);
  }

  // The log integrand is formed from terms of size alpha log u, whose
  // rounding, about eps alpha log u, passes the range of exp once alpha nears
  // 1e17: the panels then overflow.  The integrand in log u forms no such
  // terms, and takes over there.
  if (!(val > 0.0) || !std::isfinite(val)) {
    return psi_log_axis(alpha, z, tol);
  }
  return std::log(2.0) + log_peak_shifted + std::log(val);
}

}  // namespace

// log of the integral over w of exp(-alpha (e^w - 1 - w) - u*^2 (e^w - 1)^2),
// u* the peak of u^alpha exp(-u^2 + z u) (see psi_log_axis).
double psi_log_axis_relative(double alpha, double z, double tol) {
  return psi_log_axis(alpha, z, tol, /*relative=*/true);
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
