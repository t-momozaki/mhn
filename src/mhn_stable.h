// Numerically stable evaluations of closed forms that recur across the
// package.
//
// Each expression here has a textbook form that is written one way in the
// source papers and a mathematically identical form that behaves better in
// double precision.  Which one is safe depends on the sign of gamma, so the
// choice is made once, here, rather than at each call site.  The small
// building blocks several of them need -- log(1 - e^t), e^w - 1 - w, an
// exact sum of two doubles and an exact Newton step on the quadratic of the
// mode -- are kept here too, so that every caller forms them alike.

#ifndef MHN_STABLE_H
#define MHN_STABLE_H

#include <algorithm>
#include <cmath>
#include <limits>

#include "mhn_constants.h"
#include "mhn_fma.h"

namespace mhn {

// The largest finite double.
constexpr double kDoubleMax = std::numeric_limits<double>::max();

// The smallest sum of two products that the forms below take as it stands.
// A product under the smallest normal double, 2^-1022, is rounded to a
// multiple of 2^-1074, so a sum of two products such as gamma^2 + 8 beta c
// can be off by 2^-1074 whatever its size.  From 2^-1021 up that is at most
// half a unit in the last place of the sum; below, it takes digits off the
// sum, and all of them where both products underflow to 0.
constexpr double kProductSumMin = 0x1p-1021;

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
// printed form is the well-conditioned one for gamma >= 0, so each sign takes
// the branch that only ever adds.  Gao & Wang (2025) print both forms in
// their Eq. (9) (and Eq. (12) for the log axis), and the R code in their
// Appendix B picks between them on the sign of their beta, which is
// -gamma / sqrt(beta) here, the same way, except at gamma = 0, where
// Appendix B takes the conjugate form and this function the printed one
// (both exact there).
//
// Both forms square gamma, which overflows once |gamma| passes 1.34e154
// although the root itself is representable (about gamma / (2 beta) or
// c / |gamma|), and 8 beta c, 4 beta and 2 c can overflow on their own.  At
// the other end gamma^2 and 8 beta c lose digits below the smallest normal
// double, and their sum with them (kProductSumMin), although the root, about
// sqrt(c / (2 beta)) at a small beta, is far inside the range.  Where any of
// that would happen, the root is formed from ratios instead
// (positive_root_scaled below); everywhere else the forms above are used as
// they stand.
namespace detail {

// The root of positive_root without squaring anything.  With
// r = sqrt(8 beta c), the discriminant is sqrt(gamma^2 + r^2), and r itself
// is formed from sqrt(beta) sqrt(c), so that it overflows only where it
// exceeds the largest double.  Whichever of |gamma| and r is the larger is
// factored out:
//
//   * |gamma| >= r, with q = r / |gamma| <= 1 and S = 1 + sqrt(1 + q^2) in
//     [2, 1 + sqrt 2]:
//         gamma >= 0:  x = (gamma / (4 beta)) S,
//         gamma <  0:  x = 2 c / (|gamma| S) = (c / |gamma|) (2 / S);
//   * r > |gamma|, with p = |gamma| / r < 1 and T = p + sqrt(1 + p^2) in
//     [1, 1 + sqrt 2), since r / (4 beta) = 2 c / r = sqrt(c / (2 beta)):
//         gamma >= 0:  x = sqrt(c / (2 beta)) T,
//         gamma <  0:  x = sqrt(c / (2 beta)) / T.
//
// Each is a quantity within a factor 1 + sqrt 2 of the root times one of
// order one, so it overflows or underflows only where the root itself does.
// r can overflow on its own, where it exceeds the largest double, which
// needs beta and c both above 2.2e307.  It then exceeds |gamma| as well, and
// p, which can be as large as 1, is formed as
// (|gamma| / sqrt(beta)) / sqrt(8 c).
inline double positive_root_scaled(double beta, double gamma, double c) {
  const double g = std::fabs(gamma);
  const double r = 2.0 * kSqrt2 * (std::sqrt(beta) * std::sqrt(c));
  if (g >= r) {
    if (g == 0.0) return 0.0;
    const double q = r / g;
    const double S = 1.0 + std::sqrt(1.0 + q * q);
    if (gamma >= 0.0) return (0.25 * g / beta) * S;
    return (c / g) * (2.0 / S);
  }
  const double p = (opaque_value(r) <= kDoubleMax)
      ? g / r
      : (g / std::sqrt(beta)) / (2.0 * kSqrt2 * std::sqrt(c));
  const double T = p + std::sqrt(1.0 + p * p);
  const double base = std::sqrt(c) / (kSqrt2 * std::sqrt(beta));
  return (gamma >= 0.0) ? base * T : base / T;
}

}  // namespace detail

inline double positive_root(double beta, double gamma, double c) {
  const double d2 = gamma * gamma + 8.0 * beta * c;
  if (!(d2 >= kProductSumMin && opaque_value(d2) <= kDoubleMax &&
        beta <= 0.25 * kDoubleMax && c <= 0.5 * kDoubleMax)) {
    return detail::positive_root_scaled(beta, gamma, c);
  }
  const double disc = std::sqrt(d2);
  if (gamma >= 0.0) return (gamma + disc) / (4.0 * beta);
  return (c > 0.0) ? 2.0 * c / (disc - gamma) : 0.0;
}

// log of positive_root(beta, gamma, c), for beta > 0 and c > 0, finite also
// where the root itself lies beyond the largest double or below the
// smallest positive one: about gamma / (2 beta) for a positive tilt at a
// small beta, or c / |gamma| for a strong negative one.  Where the root is a
// normal double its logarithm is taken; otherwise (a root below the
// smallest normal double has lost digits, and its logarithm with them) the
// forms of detail::positive_root_scaled are taken on the log scale, as the
// log of the leading factor plus the log of the factor of order one.
inline double log_positive_root(double beta, double gamma, double c) {
  const double x = positive_root(beta, gamma, c);
  if (x >= std::numeric_limits<double>::min() &&
      opaque_value(x) <= kDoubleMax) {
    return std::log(x);
  }
  const double g = std::fabs(gamma);
  const double r = 2.0 * kSqrt2 * (std::sqrt(beta) * std::sqrt(c));
  if (g >= r) {
    const double q = r / g;
    const double S = 1.0 + std::sqrt(1.0 + q * q);
    if (gamma >= 0.0) return std::log(0.25 * g) - std::log(beta) + std::log(S);
    return std::log(c) - std::log(g) + std::log(2.0 / S);
  }
  const double p = g / r;
  const double log_T = std::log(p + std::sqrt(1.0 + p * p));
  const double log_base = 0.5 * (std::log(c) - kLn2 - std::log(beta));
  return (gamma >= 0.0) ? log_base + log_T : log_base - log_T;
}

// sqrt(gamma^2 + 8 beta c), the square root of the discriminant of
// 2 beta x^2 - gamma x - c = 0, for beta > 0 and c of either sign; NaN where
// it is negative.  As in positive_root, the literal form is used wherever
// gamma^2 + 8 beta c is finite and gamma^2 and 8 beta |c| are not both below
// kProductSumMin, and otherwise, with s = sqrt(8 beta |c|) formed from
// square roots, the same quantity as hypot(gamma, s) for c >= 0 and
// sqrt(|gamma| - s) sqrt(|gamma| + s) for c < 0, neither of which squares
// anything.
inline double root_discriminant(double beta, double gamma, double c) {
  const double d2 = gamma * gamma + 8.0 * beta * c;
  const bool tiny = gamma * gamma < kProductSumMin &&
                    8.0 * beta * std::fabs(c) < kProductSumMin;
  if (opaque_value(std::fabs(d2)) <= kDoubleMax && !tiny) return std::sqrt(d2);
  const double s = 2.0 * kSqrt2 * (std::sqrt(beta) * std::sqrt(std::fabs(c)));
  if (c >= 0.0) return std::hypot(gamma, s);
  const double g = std::fabs(gamma);
  return kSqrt2 * std::sqrt(g - s) * std::sqrt(0.5 * g + 0.5 * s);
}

// Interior maximum of the kernel t^(alpha-1) exp(-beta t^2 + gamma t) on
// t > 0, or 0 when there is none.  Its stationary points solve
// 2 beta t^2 - gamma t - (alpha-1) = 0, and which root is a maximum depends on
// the sign of alpha - 1 (Sun et al. 2023, Lemma 3):
//
//   * alpha > 1 (Lemma 3b): one positive root, the unique mode, for every
//     gamma -- positive_root with c = alpha - 1.
//   * alpha < 1 and gamma > 0 (Lemma 3c): the larger root
//     (gamma + sqrt(gamma^2 + 8 beta (alpha-1))) / (4 beta) is a local
//     maximum when the discriminant is positive, that is when
//     alpha > 1 - gamma^2 / (8 beta).  Lemma 3c prints >=, but at equality
//     the two stationary points merge into an inflection point and there is
//     no maximum.  positive_root assumes c >= 0, so the discriminant, which
//     can be negative in this case, is checked here.
//   * Otherwise the kernel decreases from its value at the origin, which is
//     unbounded for alpha < 1, and the result is 0: by Lemma 3d for
//     gamma < 0 with alpha <= 1 and for gamma > 0 with
//     alpha < 1 - gamma^2 / (8 beta); at gamma = 0 because the derivative of
//     the log kernel, -(2 beta t^2 + 1 - alpha) / t, is negative; and at
//     alpha = 1 - gamma^2 / (8 beta) because of the inflection point above.
//
// At alpha = 1 the same formula gives gamma / (2 beta) for gamma > 0, the mode
// of the truncated normal (Lemma 6b), and 0 otherwise.
//
// For alpha <= 1 the discriminant is formed literally wherever it is finite
// and 4 beta is, and gamma^2 and 8 beta (1 - alpha) are not both below
// kProductSumMin.  Otherwise (gamma^2 or 8 beta overflows, or both terms
// lose digits, as gamma^2 does at alpha = 1 for a gamma below 1.5e-154) the
// root is written with s = sqrt(8 beta (1 - alpha)) and q = s / gamma: a
// maximum exists when gamma > s, that is q < 1, and it lies at
// (gamma / (4 beta)) (1 + sqrt((1 - q)(1 + q))).
inline double kernel_interior_peak(double alpha, double beta, double gamma) {
  if (alpha > 1.0) return positive_root(beta, gamma, alpha - 1.0);
  const double disc = gamma * gamma + 8.0 * beta * (alpha - 1.0);
  const bool tiny = gamma * gamma < kProductSumMin &&
                    8.0 * beta * (1.0 - alpha) < kProductSumMin;
  if (opaque_value(std::fabs(disc)) <= kDoubleMax &&
      beta <= 0.25 * kDoubleMax && !tiny) {
    return (gamma > 0.0 && disc > 0.0)
               ? (gamma + std::sqrt(disc)) / (4.0 * beta)
               : 0.0;
  }
  if (!(gamma > 0.0)) return 0.0;
  const double s =
      2.0 * kSqrt2 * (std::sqrt(beta) * std::sqrt(1.0 - alpha));
  const double q = s / gamma;
  if (!(q < 1.0)) return 0.0;
  return (0.25 * gamma / beta) * (1.0 + std::sqrt((1.0 - q) * (1.0 + q)));
}

// Optimal rate delta of the sqrt-Gamma proposal (X = sqrt(T),
// T ~ Gamma(alpha/2, rate delta)) in Algorithm 1 of Sun et al. (2023,
// Theorem 1, optimal constants):
//
//     delta_opt = beta + ( gamma^2 - gamma sqrt(gamma^2 + 8 alpha beta) )
//                        / (4 alpha).
//
// The first part of Theorem 1 (printed "a.") requires delta to lie strictly
// inside (0, beta), and delta_opt -> 0+ as gamma grows.  Evaluated literally
// the bracket is a difference of two quantities that both grow like gamma^2,
// so for large gamma the result is dominated by rounding: it reaches 0, then
// turns negative, then snaps to beta.  A negative delta_opt is not merely
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
// second form is the one for gamma < 0.
//
// The first form is evaluated as 8 alpha r^2 with r = beta / (sqrt(D) + gamma).
// delta_opt is of order beta, but beta^2 is not: it overflows once beta
// exceeds about 1e154, which would make delta_opt +Inf, and below about
// 1e-154 it falls into the subnormal range, losing digits until it underflows
// to zero near 1e-162.  The ratio r is of order sqrt(beta), so r^2 stays in
// range wherever D does.
//
// D itself overflows once |gamma| passes 1.34e154 or 8 alpha beta passes the
// largest double, and the second form's square can overflow before D does.
// Written with p = |gamma| / sqrt(8 alpha beta) and
// sigma = p + sqrt(1 + p^2) >= 1, the two forms are
//
//     delta_opt = beta / sigma^2   for gamma >= 0,
//     delta_opt = beta sigma^2     for gamma < 0,
//
// and where the forms above would leave the double range these are used,
// with sqrt(8 alpha beta) formed from square roots and sqrt(1 + p^2) as a
// hypot, so that nothing squares a large number.
//
// This helper requires that sqrt(8 alpha beta) be finite, which it is not
// once alpha and beta both exceed 2.2e307, and that D be a normal double;
// outside that range its result is not accurate.
inline double sun_delta_opt(double alpha, double beta, double gamma) {
  const double d2 = gamma * gamma + 8.0 * alpha * beta;
  if (opaque_value(d2) <= kDoubleMax) {
    const double disc = std::sqrt(d2);
    if (gamma >= 0.0) {
      const double r = beta / (disc + gamma);
      return 8.0 * alpha * r * r;
    }
    const double numer = disc - gamma;  // gamma < 0, so this adds
    const double v = numer * numer / (8.0 * alpha);
    if (opaque_value(v) <= kDoubleMax) return v;
  }
  const double s = 2.0 * kSqrt2 * (std::sqrt(alpha) * std::sqrt(beta));
  const double p = std::fabs(gamma) / s;
  const double sigma = p + std::hypot(1.0, p);
  if (gamma >= 0.0) return (beta / sigma) / sigma;
  return (beta * sigma) * sigma;
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
//
// The t < 0 branch cancels as t -> 0: exp(t) carries a rounding of about
// eps, and 1 - exp(t) is about |t|, so the result keeps a relative accuracy
// of only about eps / |t|.  The envelope pieces that call this pass a
// log-slope times a width that is not small; log1m_exp below serves a
// t <= 0 that can come close to 0.
inline double log_abs_expm1(double t) {
  if (t > 0.0) return t + std::log1p(-std::exp(-t));
  if (t < 0.0) return std::log1p(-std::exp(t));
  return -std::numeric_limits<double>::infinity();
}

// log(1 - exp(t)) for t <= 0, taking whichever of two forms does not cancel
// (Maechler 2012, "Accurately computing log(1 - exp(-|a|))"; R's own
// R_Log1_Exp): log(-expm1(t)) for t > -log 2, where expm1 keeps the relative
// accuracy of 1 - exp(t) however close t is to 0, and log1p(-exp(t))
// otherwise, where exp(t) <= 1/2 and log1p keeps the accuracy of a small
// result.  log1p(-exp(t)) alone returns -Inf for t above about -5.6e-17,
// where exp(t) rounds to 1, and is visibly wrong somewhat below that
// (-36.737 against -36.841 at t = -1e-16), where the answer is an ordinary
// number.
inline double log1m_exp(double t) {
  return (t > -kLn2) ? std::log(-std::expm1(t)) : std::log1p(-std::exp(t));
}

// e^w - 1 - w.  Near w = 0 the difference expm1(w) - w cancels, losing about
// log2(2 / |w|) bits, so for |w| <= 1/2 the Taylor series sum_k w^k / k! is
// summed from k = 2 instead, by Horner's rule from k = 15 down; the terms
// past k = 15 add less than 1e-17 of the sum there.
inline double expm1_minus_id(double w) {
  if (std::fabs(w) > 0.5) return std::expm1(w) - w;
  double p = 1.0 / 1307674368000.0;   // 1/15!
  p = p * w + 1.0 / 87178291200.0;    // 1/14!
  p = p * w + 1.0 / 6227020800.0;
  p = p * w + 1.0 / 479001600.0;
  p = p * w + 1.0 / 39916800.0;
  p = p * w + 1.0 / 3628800.0;
  p = p * w + 1.0 / 362880.0;
  p = p * w + 1.0 / 40320.0;
  p = p * w + 1.0 / 5040.0;
  p = p * w + 1.0 / 720.0;
  p = p * w + 1.0 / 120.0;
  p = p * w + 1.0 / 24.0;
  p = p * w + 1.0 / 6.0;
  p = p * w + 0.5;                    // 1/2!
  return w * w * p;
}

// a + b = s + err exactly, s being the rounded sum a + b (Knuth's TwoSum,
// which needs no ordering of |a| and |b|).  The operands are read through
// opaque_value: a compiler that contracts across statements (clang with
// -ffp-contract=fast, and GCC by default where the target has a fused
// multiply-add) could otherwise fuse the product that formed a or b into
// the sum, so that s would round the exact product plus the other operand
// while err is formed as the error of a + b.  At extreme parameters that
// moved pmhn() by up to 0.999 and RTDR's placement of the mode by many
// spacings of the doubles.
inline double two_sum(double a, double b, double& err) {
  const double x = opaque_value(a);
  const double y = opaque_value(b);
  const double s = x + y;
  const double bb = s - x;
  err = (x - (s - bb)) + (y - bb);
  return s;
}

// One Newton step from x0 towards the positive root of
//
//     2 beta x^2 - gamma x - (c + c_lo) = 0,
//
// where c + c_lo is a constant held as a double c and the part c_lo that c
// leaves out.  x0 is taken to lie within a few units in its last place of
// the root, and x0 plus the step returned gives the root to about twice
// double precision.  The residual 2 beta x0^2 - gamma x0 - c - c_lo is
// formed without rounding error: each product is split exactly by fma, and
// the three leading terms 2 beta x0^2, -gamma x0 and -c, which cancel to a
// number of order ulp(x0) times the derivative, are summed exactly by
// two_sum; the rounding errors are added back at the end.  The derivative
// 4 beta x0 - gamma = sqrt(gamma^2 + 8 beta c) is a sum of two positive terms
// for gamma <= 0 and at least gamma for gamma > 0, so it is formed to full
// relative precision.  A product that leaves the double range gives a step
// that is NaN or infinite, which the callers handle.
inline double quadratic_root_step(double beta, double gamma, double c,
                                  double c_lo, double x0) {
  const double a = beta * x0,  a_lo = mhn::fma(beta, x0, -a);
  const double p = a * x0,     p_lo = mhn::fma(a, x0, -p);
  const double q = gamma * x0, q_lo = mhn::fma(gamma, x0, -q);
  double e1, e2;
  const double s1 = two_sum(2.0 * p, -q, e1);
  const double s2 = two_sum(s1, -c, e2);
  const double r =
      s2 + ((e1 + e2) + ((2.0 * (p_lo + a_lo * x0) - q_lo) - c_lo));
  return -r / (4.0 * a - gamma);
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

// log1p(u) - u, without the cancellation the literal form has for small u.
// std::log1p(u) is accurate to within about one rounding, an absolute error of
// about eps |u|, while the difference is about u^2/2.  The subtraction
// therefore keeps a relative accuracy of only about 2 eps / |u|: some eleven
// digits at u = 1e-5, seven at 1e-9, and six at the u = d/m = 2e-10 that the
// mode-centred region A ordinate of the RTDR sampler meets at
// gamma_norm = 1e10.  Below |u| = 1e-4 the Taylor expansion of log(1+u) with
// its linear term removed is used instead.  Truncated after u^4, its relative
// error is about 0.4 |u|^3, at most 4e-13.
inline double log1p_minus_u(double u) {
  if (std::abs(u) < 1e-4) {
    return -u * u * (0.5 - u / 3.0 + u * u / 4.0);
  }
  return std::log1p(u) - u;
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
//
// Where D would overflow, the same two forms are written with
// r = sqrt(8 (alpha-1) beta), formed from square roots, and
// p = |gamma| / r: (sqrt(D) - gamma) / 2 = (r/2) / (p + sqrt(1 + p^2)) for
// gamma >= 0, and (sqrt(D) + |gamma|) / 2 = (r/2) (p + sqrt(1 + p^2)) for
// gamma < 0.
//
// This helper requires that r and |gamma| / r be finite (r is not once
// alpha - 1 and beta both exceed 2.2e307) and that 8 (alpha-1) beta be a
// normal double; outside that range its result is not accurate.
inline double sun_tangency_slope(double alpha, double beta, double gamma) {
  const double d2 = gamma * gamma + 8.0 * (alpha - 1.0) * beta;
  if (opaque_value(d2) <= kDoubleMax) {
    const double disc = std::sqrt(d2);
    if (gamma >= 0.0) return 4.0 * (alpha - 1.0) * beta / (disc + gamma);
    return 0.5 * (disc - gamma);
  }
  const double half_r =
      kSqrt2 * std::sqrt(alpha - 1.0) * std::sqrt(beta);
  const double p = std::fabs(gamma) / (2.0 * half_r);
  const double sigma = p + std::hypot(1.0, p);
  return (gamma >= 0.0) ? half_r / sigma : half_r * sigma;
}

// How far past its peak the MHN kernel x^(alpha-1) exp(-beta x^2 + gamma x)
// must travel before its logarithm has fallen by L.
//
// Writing s for that distance, t = s/m, and using the stationarity condition
// at the peak m, the change in the logarithm is
//
//     (alpha-1) [ log(1 + t) - t ] - beta s^2,
//
// two non-positive terms when alpha > 1.  Either alone reaching -L is enough,
// so the smaller of the two distances is a safe reach.  The quadratic term
// needs |s| = sqrt(L / beta).  The bracket needs t - log(1 + t) = c, with
// c = L / (alpha-1), an equation with no closed-form root.  An upper bound for
// the root comes from
//
//     t - log(1 + t) >= t^2 / (2 (1 + |t|))      for every t > -1.
//
// For t >= 0 the difference of the two sides is 0 at t = 0 and has derivative
// t^2 / (2 (1 + t)^2) >= 0; for -1 < t < 0 the left side is the series
// t^2/2 + |t|^3/3 + |t|^4/4 + ..., which is at least t^2/2.  So the bracket
// term has fallen by at least L once |t| reaches c + sqrt(c^2 + 2c), the
// positive root of t^2 = 2c (1 + t).  For small c -- large alpha -- that is
// about sqrt(2c), because near the peak the bracket is quadratic, t^2/2;
// charging it linearly, as t, would stop at |t| = c, where it has fallen by
// far less than L (by about L c / 2).
//
// The bound is the starting point, not the answer: it overshoots the root by
// a factor that grows from 1 at small c to 2 at large c (at c = 5 it is 10.9
// against 7.1, where the bracket has fallen by 76 nats rather than 45), and
// the quadrature pays for every stretch of range on which the kernel is
// already negligible.  Newton's method on g(t) = t - log(1 + t) - c brings it
// down to the root.  For t > 0, g is increasing and convex, so from a point
// where g >= 0 each step moves down without passing the root, and every
// iterate is still a safe reach.  After four steps t agrees with the root to
// within about 1e-12 relative for every c from 1e-306 to 2e17.  The largest
// departures are just above t = 1e-4, where g is formed by the subtraction
// log1p(t) - t and carries a rounding of about eps / t relative; where t lies
// below the root it does so only by that error in g, a shortfall of at most
// about 2e-12 of L.
//
// The same t serves on the left of the peak, where the bracket falls faster:
// for 0 < t < 1, -t - log(1 - t) = t^2/2 + t^3/3 + ... is at least
// t - log(1 + t) = t^2/2 - t^3/3 + ..., and for t >= 1 the range is cut at
// the origin anyway.
//
// The bracket's figure is much the shorter whenever a negative tilt squeezes
// the peak towards the origin -- with L = 45, at (alpha, beta, gamma) =
// (2.5, 1, -30) it is 1.7 against 6.7 -- and using only the Gaussian figure
// there leaves a range far wider than the region that holds the mass, which a
// fixed-order rule then samples only sparsely.
//
// For alpha <= 1 the bracket has the opposite sign and cannot be relied on, so
// only the quadratic term counts.
//
// L / beta overflows for beta below L over the largest double, about 3e-307
// for the L = 58 of the distribution function, while its square root is
// still far inside the range; it is then taken as sqrt(L) / sqrt(beta).  An
// infinite Gaussian reach would leave m t as the reach, far wider than the
// law, and a quadrature over it could miss the law altogether.
inline double kernel_decay_reach(double alpha, double beta, double m, double L) {
  const double q = L / beta;
  const double gaussian =
      (opaque_value(q) <= kDoubleMax) ? std::sqrt(q)
                                       : std::sqrt(L) / std::sqrt(beta);
  if (alpha > 1.0 && m > 0.0) {
    const double c = L / (alpha - 1.0);
    double t = c + std::sqrt(c * c + 2.0 * c);
    // g(t) = t - log(1 + t) - c and g'(t) = t / (1 + t); log1p_minus_u keeps
    // g accurate for the small t that large alpha brings.
    for (int k = 0; k < 4; ++k) t -= (-log1p_minus_u(t) - c) * (1.0 + t) / t;
    return std::min(gaussian, m * t);
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
//
// L'(a)^2 overflows once the slope passes 1.34e154, and 4 beta L can too; the
// reach would then come out 0.  There, and where their sum has lost digits
// below kProductSumMin, the same quantity is written with everything halved,
// L / (|L'(a)|/2 + hypot(|L'(a)|/2, sqrt(beta L))), which squares nothing.
inline double kernel_decay_reach_at(double beta, double slope, double L) {
  const double s = std::fabs(slope);
  const double d2 = s * s + 4.0 * beta * L;
  if (d2 >= kProductSumMin && opaque_value(d2) <= kDoubleMax) {
    return 2.0 * L / (s + std::sqrt(d2));
  }
  const double hs = 0.5 * s;
  return L / (hs + std::hypot(hs, std::sqrt(beta) * std::sqrt(L)));
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
// alpha = 0.05, gamma = -1000, an anchor of 1e-6 would give
// (alpha-1)/a = -9.5e5 and, for L = 58, a reach of 6e-5, cutting off a tail
// that runs out to 1e-2.  Count only the exponential part there and let the power-law decay
// be a bonus.
inline double kernel_decay_slope(double alpha, double beta, double gamma,
                                 double a) {
  // 2 beta overflows for beta above about 9e307; beta a is then taken first.
  double exponential = gamma - 2.0 * beta * a;
  if (!std::isfinite(exponential)) exponential = gamma - 2.0 * (beta * a);
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
//
// The fraction is evaluated backwards from a depth of n levels, and the depth
// it needs falls quickly as x grows: n = 16 + 400 / x^2 levels, capped at
// 100, leave a truncation error of at most 6.3e-23 relative, at x = 2 where
// the cap applies (1e-25 at x = 4, 3e-30 at x = 10), against erfc evaluated
// at 300 bits at 401 points from 2 to 1e4 -- far below the rounding of the
// result.  The rounded results agree bit for bit with a fixed depth of 100
// levels at every one of 1.2e8 points tried between 2 and 1e150, while a
// large x takes 17 levels instead of 100.
inline double erfcx_cf(double x) {
  const double levels = 16.0 + 400.0 / (x * x);
  const int n = (levels < 100.0) ? static_cast<int>(std::ceil(levels)) : 100;
  double f = 0.0;
  for (int k = n; k >= 1; --k) f = (0.5 * k) / (x + f);
  return 1.0 / (std::sqrt(kPi) * (x + f));
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
// a negligible sliver that the quadrature would handle correctly anyway, and
// the two agree to rounding.
//
// x0 is the positive root of beta x^2 + |gamma| x = budget (kernel_head_x0),
// written so that it only adds.  Once |gamma| passes 1.34e154 its square
// overflows, and the root, about budget / |gamma|, is then formed with
// everything halved, as in kernel_decay_reach_at.  So it is where gamma^2
// and 4 beta budget both fall below the normal range, at a beta below about
// 1e-300 with a small gamma: once both underflow, the literal form gives
// 2 budget / |gamma|, far past the root sqrt(budget / beta), and the closed
// form would be used where its expansion does not hold.  The two
// corrections in the bracket are at most about budget in size, but their
// products alpha beta and gamma alpha, taken first, can overflow for a huge
// alpha; they are then regrouped as (beta x0) x0 and gamma x0, each at most
// budget.
struct KernelHead {
  double value;  // the integral over [0, x0], divided by exp(log_scale)
  double x0;     // where the closed form stops and quadrature should start
};

// The budget: the first discarded term, at most budget^2 / 2, is 1e-17.
constexpr double kHeadBudget = 4.5e-9;  // sqrt(2 * 1e-17)

inline double kernel_head_x0(double beta, double gamma) {
  const double budget = kHeadBudget;
  const double g = std::fabs(gamma);
  const double d2 = g * g + 4.0 * beta * budget;
  return (d2 >= kProductSumMin && opaque_value(d2) <= kDoubleMax)
      ? 2.0 * budget / (g + std::sqrt(d2))
      : budget / (0.5 * g + std::hypot(0.5 * g,
                                       std::sqrt(beta) * std::sqrt(budget)));
}

inline KernelHead kernel_head(double alpha, double beta, double gamma,
                              double x_limit, double log_scale) {
  if (!(x_limit > 0.0) || !(alpha > 0.0)) return KernelHead{0.0, 0.0};
  double x0 = kernel_head_x0(beta, gamma);
  if (x0 > x_limit) x0 = x_limit;
  const double log_lead = alpha * std::log(x0) - std::log(alpha) - log_scale;
  double shape = 1.0 - alpha * beta * x0 * x0 / (alpha + 2.0)
                     + gamma * alpha * x0 / (alpha + 1.0);
  if (!std::isfinite(shape)) {
    shape = 1.0 - (beta * x0) * x0 * (alpha / (alpha + 2.0))
                + (gamma * x0) * (alpha / (alpha + 1.0));
  }
  const double v = (log_lead > -700.0 && shape > 0.0)
                       ? std::exp(log_lead) * shape
                       : 0.0;
  return KernelHead{v, x0};
}

}  // namespace mhn

#endif  // MHN_STABLE_H
