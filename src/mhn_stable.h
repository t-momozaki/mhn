// Numerically stable evaluations of closed forms that recur across the
// package.
//
// Each expression here has a textbook form that is written one way in the
// source papers and a mathematically identical form that behaves better in
// double precision.  Which one is safe depends on the sign of gamma, so the
// choice is made once, here, rather than at each call site.

#ifndef MHN_STABLE_H
#define MHN_STABLE_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

// Positive root of  2 beta x^2 - gamma x - c = 0,  i.e.
//
//     x = ( gamma + sqrt(gamma^2 + 8 beta c) ) / (4 beta),      beta > 0, c >= 0.
//
// This is the shape of the MHN mode (Sun et al. 2023, Lemma 3b, with
// c = alpha - 1), of the mode of the log-axis density (c = alpha), and of the
// integrand peak used by the quadrature paths.
//
// Written literally, the numerator subtracts two nearly equal quantities when
// gamma < 0: sqrt(gamma^2 + 8 beta c) -> |gamma| as |gamma| grows, so the
// leading digits cancel and the result loses all precision (it reaches exactly
// 0 near gamma = -1e8).  Multiplying through by the conjugate gives
//
//     x = 2 c / ( sqrt(gamma^2 + 8 beta c) - gamma ),
//
// which adds two positive quantities for gamma < 0 and is exact there.  The
// original form is the well-conditioned one for gamma >= 0, so each sign takes
// the branch that only ever adds.  Gao and Wang (2025, Eq. 9) give both forms
// for the same reason.
inline double positive_root(double beta, double gamma, double c) {
  const double disc = std::sqrt(gamma * gamma + 8.0 * beta * c);
  if (gamma >= 0.0) return (gamma + disc) / (4.0 * beta);
  return (c > 0.0) ? 2.0 * c / (disc - gamma) : 0.0;
}

// Optimal squared-scale of the sqrt-Gamma proposal in Algorithm 1 of
// Sun et al. (2023, Theorem 1b):
//
//     delta_opt = beta + ( gamma^2 - gamma sqrt(gamma^2 + 8 alpha beta) )
//                        / (4 alpha).
//
// Theorem 1a requires delta_opt to lie strictly inside (0, beta), and
// delta_opt -> 0+ as gamma grows.  Evaluated literally the bracket is a
// difference of two quantities that both grow like gamma^2, so for large
// gamma the result is dominated by rounding: it reaches 0, then turns
// negative, then snaps to beta.  A negative delta_opt is not merely
// inaccurate -- it makes the proposal's Gamma scale negative.
//
// Factoring gamma out and rationalising removes every subtraction:
//
//     delta_opt = beta (sqrt(D) - gamma) / (sqrt(D) + gamma),   D = gamma^2 + 8 alpha beta,
//
// and since (sqrt(D) - gamma)(sqrt(D) + gamma) = 8 alpha beta, this is
//
//     delta_opt = 8 alpha beta^2 / (sqrt(D) + gamma)^2          for gamma >= 0,
//     delta_opt = (sqrt(D) + |gamma|)^2 / (8 alpha)             for gamma < 0,
//
// each of which only ever adds.  Algorithm 1 is only used for gamma > 0; the
// second branch is kept so the helper is correct wherever it is called.
inline double sun_delta_opt(double alpha, double beta, double gamma) {
  const double disc = std::sqrt(gamma * gamma + 8.0 * alpha * beta);
  if (gamma >= 0.0) {
    const double denom = disc + gamma;
    return 8.0 * alpha * beta * beta / (denom * denom);
  }
  const double numer = disc - gamma;  // gamma < 0, so this adds
  return numer * numer / (8.0 * alpha);
}

// log(1 + exp(t)), evaluated without overflowing for large t.
inline double log1p_exp(double t) {
  return (t > 0.0) ? t + std::log1p(std::exp(-t)) : std::log1p(std::exp(t));
}

// log|exp(t) - 1|, evaluated without overflowing for large |t|.
//
// Forming expm1(t) first and taking its logarithm overflows to +Inf once
// t exceeds about 709, which is well inside the range reached by a steep
// envelope piece (slope times width).  Factoring the dominant exponential out
// of each branch keeps every intermediate in range:
//
//     t > 0:  exp(t) - 1 = exp(t) (1 - exp(-t))  ->  t + log1p(-exp(-t))
//     t < 0:  1 - exp(t)                         ->  log1p(-exp(t))
inline double log_abs_expm1(double t) {
  if (t > 0.0) return t + std::log1p(-std::exp(-t));
  if (t < 0.0) return std::log1p(-std::exp(t));
  return -std::numeric_limits<double>::infinity();
}

// log(1 + u (exp(t) - 1)) for u in [0, 1], the inverse-CDF core of an
// exponential envelope piece of log-slope t = slope * width.
//
// The direct form log1p(u * expm1(t)) overflows for t above about 709, and
// then returns +Inf for every u, so the piece proposes nothing inside its own
// support and contributes no accepted draws.  For t > 0 the dominant factor
// is u exp(t), and pulling it out leaves only bounded terms:
//
//     1 + u (e^t - 1) = u e^t (1 + ((1-u)/u) e^{-t}).
inline double log1p_u_expm1(double u, double t) {
  if (u <= 0.0) return 0.0;
  if (t <= 0.0) return std::log1p(u * std::expm1(t));
  const double log_u = std::log(u);
  return log_u + t + log1p_exp(-t + std::log1p(-u) - log_u);
}

// The Algorithm 1 tangency slope 2 beta mu_opt - gamma.
//
// With mu_opt = (gamma + sqrt(D)) / (4 beta) and D = gamma^2 + 8 (alpha-1) beta,
// this is (sqrt(D) - gamma) / 2, a difference of two quantities that both grow
// like gamma.  It reaches exactly 0 by gamma = 1e9, where the true value is
// 6e-9 -- and it is a divisor of the acceptance test and the argument of a
// logarithm, so zero there is fatal.  Since D - gamma^2 = 8 (alpha-1) beta,
//
//     (sqrt(D) - gamma) / 2 = 4 (alpha-1) beta / (sqrt(D) + gamma),
//
// which only ever adds for gamma > 0, the sign Algorithm 1 runs at.
inline double sun_tangency_slope(double alpha, double beta, double gamma) {
  const double disc = std::sqrt(gamma * gamma + 8.0 * (alpha - 1.0) * beta);
  if (gamma >= 0.0) return 4.0 * (alpha - 1.0) * beta / (disc + gamma);
  return 0.5 * (disc - gamma);
}

// How far past its peak the MHN kernel x^(alpha-1) exp(-beta x^2 + gamma x)
// must travel before its logarithm has fallen by L.
//
// Writing s for that distance and using the stationarity condition at the peak
// m, the drop is
//
//     (alpha-1) [ log(1 + s/m) - s/m ] - beta s^2,
//
// two non-positive terms when alpha > 1.  Either alone reaching L is enough,
// so the smaller of the two distances is a safe reach: the quadratic term
// needs sqrt(L / beta), and the bracket, which behaves like -(alpha-1) s/m
// once s exceeds m, needs L m / (alpha-1).  The second is much the shorter
// whenever a negative tilt squeezes the peak towards the origin -- at
// (alpha, beta, gamma) = (2.5, 1, -30) it is 2.0 against 7.8 -- and using only
// the Gaussian figure there leaves a range twenty times wider than the support,
// on which a fixed-order rule sees nothing.
//
// For alpha <= 1 the bracket has the opposite sign and cannot be relied on, so
// only the quadratic term counts.
inline double kernel_decay_reach(double alpha, double beta, double m, double L) {
  const double gaussian = std::sqrt(L / beta);
  if (alpha > 1.0 && m > 0.0) {
    return std::min(gaussian, L * m / (alpha - 1.0));
  }
  return gaussian;
}

// The same reach, measured from a point that is not the peak.
//
// At the peak the first-order terms cancel through the stationarity condition
// and the drop is governed by the curvature alone, which is what the reach
// above assumes.  Anchored anywhere else the linear term survives: the drop
// over a distance s is at least |L'(a)| s + beta s^2, so
//
//     s = 2 L / ( |L'(a)| + sqrt( L'(a)^2 + 4 beta L ) ),
//
// written so that it only ever adds.  At L'(a) = 0 this is exactly the
// Gaussian reach sqrt(L / beta); for a steep slope it is L / |L'(a)|, which is
// the point -- a strong tilt makes the kernel die within a fraction of the
// Gaussian width, and a range set by the width alone leaves the whole mass
// inside a sliver a fixed-order rule never samples.
inline double kernel_decay_reach_at(double beta, double slope, double L) {
  const double s = std::fabs(slope);
  return 2.0 * L / (s + std::sqrt(s * s + 4.0 * beta * L));
}

// The slope to hand kernel_decay_reach_at, given an anchor a.
//
// Over a distance s the log kernel changes by
//     (alpha-1) log(1 + s/a) + s (gamma - 2 beta a) - beta s^2,
// and the reach above needs a linear coefficient that under-states the decay,
// so that solving for a drop of L guarantees at least that much.
//
// For alpha > 1 the power term is bounded by (alpha-1) s / a, which folds into
// the exponential part to give exactly L'(a) -- zero at the peak, where the
// reach is then the Gaussian one.  For alpha <= 1 that bound points the other
// way: the term is negative and grows only logarithmically, so charging it
// linearly over-states the decay and truncates the range far too early.  At
// alpha = 0.05, gamma = -1000, an anchor of 1e-6 gave (alpha-1)/a = -9.5e5 and
// a reach of 8e-5, cutting off a tail that runs out to 1e-2.  Count only the
// exponential part there and let the power-law decay be a bonus.
inline double kernel_decay_slope(double alpha, double beta, double gamma,
                                 double a) {
  const double exponential = gamma - 2.0 * beta * a;
  if (alpha > 1.0 && a > 0.0) return (alpha - 1.0) / a + exponential;
  return exponential;
}

// Scaled complementary error function erfcx(x) = exp(x^2) erfc(x), for x >= 2,
// by the continued fraction
//
//     erfcx(x) = 1 / ( sqrt(pi) ( x + (1/2)/(x + 1/(x + (3/2)/(x + ...))) ) ).
//
// The point of the scaling is that erfcx decays only like 1/(x sqrt(pi)), so it
// stays in range where erfc itself underflows and where exp(x^2) overflows.
// Evaluated backwards from a fixed depth, which converges quickly for x >= 2.
inline double erfcx_cf(double x) {
  double f = 0.0;
  for (int k = 100; k >= 1; --k) f = (0.5 * k) / (x + f);
  return 1.0 / (std::sqrt(M_PI) * (x + f));
}

// The head of int_0^X x^(alpha-1) exp(-beta x^2 + gamma x) dx, in closed form,
// divided by exp(log_scale).
//
// The factor x^(alpha-1) carries mass all the way down to x = exp(-1/alpha),
// which for alpha below about 1e-3 is smaller than any representable double --
// at alpha = 1e-6 it is exp(-1e6).  No quadrature can place an abscissa there:
// tanh_sinh reaches about 1e-300 and stops, which accounts for only
// (1 - 1e-300^alpha)/alpha of the integral.  So the head has to be done in
// closed form rather than sampled.
//
// Take x0 short enough that exp(-beta x^2 + gamma x) agrees with its expansion
// 1 - beta x^2 + gamma x to within 1e-17 over [0, x0] -- the first discarded
// term is about (beta x^2 + |gamma| x)^2 / 2, which gives the bound below --
// and integrate term by term:
//
//   int_0^x0 x^(a-1) (1 - beta x^2 + gamma x) dx
//       = x0^a/a - beta x0^(a+2)/(a+2) + gamma x0^(a+1)/(a+1)
//       = (x0^a / a) [ 1 - a beta x0^2/(a+2) + gamma a x0/(a+1) ].
//
// Factoring the leading term out keeps the whole expression in log space, so it
// neither overflows for tiny alpha nor cancels.  For moderate alpha the head is
// a negligible sliver that the quadrature handled correctly anyway, and the two
// agree to rounding.
struct KernelHead {
  double value;  // the integral over [0, x0], divided by exp(log_scale)
  double x0;     // where the closed form stops and quadrature should start
};

inline KernelHead kernel_head(double alpha, double beta, double gamma,
                              double x_limit, double log_scale) {
  if (!(x_limit > 0.0) || !(alpha > 0.0)) return KernelHead{0.0, 0.0};
  const double budget = 4.5e-9;  // sqrt(2 * 1e-17)
  const double g = std::fabs(gamma);
  double x0 = 2.0 * budget / (g + std::sqrt(g * g + 4.0 * beta * budget));
  if (x0 > x_limit) x0 = x_limit;
  const double log_lead = alpha * std::log(x0) - std::log(alpha) - log_scale;
  const double shape = 1.0 - alpha * beta * x0 * x0 / (alpha + 2.0)
                           + gamma * alpha * x0 / (alpha + 1.0);
  const double v = (log_lead > -700.0 && shape > 0.0)
                       ? std::exp(log_lead) * shape
                       : 0.0;
  return KernelHead{v, x0};
}

}  // namespace mhn

#endif  // MHN_STABLE_H
