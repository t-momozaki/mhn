// mhn_vanishing.h
//
// The law at a subnormal shape, alpha below the smallest normal double,
// 2^-1022 (about 2.2e-308).
//
// There half of alpha is not always a double -- half of 3 * 2^-1074 rounds to
// 2 * 2^-1074 -- and 1 / alpha passes the largest double once alpha is below
// about 5.6e-309.  The paths written for an ordinary shape go wrong with
// them: the series for the normalising constant and for the distribution
// function start from Gamma(alpha / 2), the square-root-Gamma case of
// gamma = 0 takes alpha / 2 as its shape, and the quadratures of the moments
// hold the mass by the origin, about 1 / alpha.
//
// At such an alpha, t^alpha = 1 to double precision at every positive double
// t (|alpha log t| < 1.7e-305), and the law takes a form in which alpha only
// multiplies.  With U = sqrt(beta) X, z = gamma / sqrt(beta) and the
// normalising integral of Sun et al. (2023), Lemma 1a,
//
//   int_0^Inf t^(alpha-1) exp(-t^2 + z t) dt = (1 + w) / alpha,
//   w = alpha K0(z),
//   K0(z) = int_0^1 (exp(-t^2 + z t) - 1) / t dt
//           + int_1^Inf exp(-t^2 + z t) / t dt,
//
// the density of X is alpha exp(-beta x^2 + gamma x) / (x (1 + w)), the upper
// tail P(X > x) is alpha J(sqrt(beta) x) / (1 + w) with
// J(u) = int_u^Inf exp(-t^2 + z t) / t dt, the lower tail P(X <= x) is
// (1 + alpha G(sqrt(beta) x)) / (1 + w) with G(u) = K0(z) - J(u)
// = log u + int_0^u (exp(-t^2 + z t) - 1) / t dt, and, by Lemma 2a,
// E(U^k) = alpha A_k / (1 + w) with A_k = int_0^Inf t^(k-1) exp(-t^2 + z t) dt.
// w is negative for z below about 0.3, where K0 is, and for a large z it is
// the mass at the peak of the kernel against the mass by the origin,
// 1 / alpha.
//
// These forms serve the whole law while w <= 2^64, up to a tilt z between
// about 55 and 56.3 as alpha falls.  Above it the mass by the origin is at
// most 2^-64 (5.4e-20) of the whole, and
// the density, the upper tail and the moments do not depend on alpha to
// double precision: the ordinary paths are right for them.  The lower tail
// far below the peak still does.  There it is (1 + alpha G(u)) / (1 + w), the
// mass by the origin plus the part of the peak below u, and where
// alpha G(u) is below 2^64 the first counts; the paths written for an
// ordinary shape lose it (at 2^-1074 they leave it out, and at
// 3 * 2^-1074 take it at 3/4 of its size).  That part of the lower tail is
// taken from the forms at every tilt (VanishingShape::lower_only).

#ifndef MHN_VANISHING_H
#define MHN_VANISHING_H

#include <limits>

namespace mhn {

// Whether alpha can be served here at all: only a subnormal shape can.  The
// full test, which also needs z and w, is VanishingShape::set.
inline bool vanishing_candidate(double alpha) {
  return alpha < std::numeric_limits<double>::min();
}

struct VanishingShape {
  double gamma = 0.0;
  double sqrt_beta = 1.0;
  double z = 0.0;
  double log_alpha = 0.0;
  double half_log_beta = 0.0;
  // log w where it is computed (z > 30), and -Inf below, where w is below
  // 1e-200 and log(1 + w) is 0 to double precision.
  double log_w = -std::numeric_limits<double>::infinity();
  double log1p_w = 0.0;
  // Past the tilt where the forms hand over to the ordinary paths (w > 2^64):
  // whether they still serve the lower tail far below the peak, and log u of
  // the point above which alpha G(u) exceeds 2^64 for certain, so that the
  // ordinary paths serve every u above it.
  bool lower_only = false;
  double log_u_far = 0.0;
  // log D(z) for z above 30 (see set_tilt), and NaN otherwise.
  double log_D_z = std::numeric_limits<double>::quiet_NaN();

  // Sets the triple, and returns whether the forms above serve the whole law
  // for it: a subnormal positive alpha, a finite z and w at most 2^64.  For a
  // subnormal alpha and a finite z with w above 2^64 it returns false and
  // sets lower_only.  It is set_tilt(beta, gamma) and then set_shape(alpha).
  bool set(double alpha, double beta, double gamma);
  // The part of set that depends on beta and gamma alone, log D among it.
  void set_tilt(double beta, double gamma);
  // The rest of set, for the beta and gamma of the last set_tilt: a caller
  // that holds several shapes at one beta and gamma evaluates log D once.
  bool set_shape(double alpha);

  // log f(x), for an x that is not NA or NaN.
  double log_density(double x) const;
  // P(X <= q) or P(X > q), on the log scale if log_p, for a q that is not NA
  // or NaN.
  double p_one(double q, bool lower_tail, bool log_p) const;
  // The quantile of p in the same forms; NA and NaN pass through, and a p
  // outside the range of a probability gives NaN.
  double q_one(double p, bool lower_tail, bool log_p) const;
  // log P(U > u) from log u, u = sqrt(beta) x, which need not be a double.
  double log_upper(double log_u) const;
  // log P(U <= u) from log u, for w > 1 and u below z / 2 + 5, which the
  // lower tail below one half lies within.
  double log_lower(double log_u) const;

  // Where lower_only is set, for the lower tail past the hand-over:
  // log P(X <= q) at a q > 0 for which alpha G(u) is at most 2^64 and
  // P(X <= q) at most one half, and false for any other q, which the
  // ordinary paths serve.
  bool log_lower_far(double q, double* log_F) const;
  // The quantile of p for a lower tail at most (1 + 2^64) / (1 + w) and one
  // half, and false for any other p, NA, NaN and a p outside the range of a
  // probability included, which the ordinary paths serve.
  bool q_lower_far(double p, bool lower_tail, bool log_p, double* x) const;

 private:
  // log u, u = sqrt(beta) x, from x > 0, and x from log u.
  double log_u_of(double x) const;
  double x_of(double log_u) const;
  // The x at which log P(X <= x) = log_F, solved on the lower tail for
  // log u up to log_u_top; 0 where P(X <= x) at the smallest positive double
  // already reaches it, and NaN where log_u_top does not bracket it.
  double lower_root(double log_F, double log_u_top) const;
};

struct VanishingMoments {
  double mean;
  double var;
  double skewness;
  double kurtosis;
};

// The four moments of the law that v holds (v.set must have returned true).
VanishingMoments vanishing_moments(const VanishingShape& v);

}  // namespace mhn

#endif  // MHN_VANISHING_H
