// Rejection samplers from Sun et al. (2023): Algorithm 1 (gamma > 0,
// alpha > 1; Normal or sqrt-Gamma proposal) and Algorithm 3 (AM-GM
// proposal; gamma < 0 in its heading, gamma <= 0 in Section 4.2 and
// Theorem 4).  Each setup is built once per parameter triple via
// build_sun_algoX(), then reused for repeated sample_sun_algoX() draws.

#ifndef MHN_SUN_H
#define MHN_SUN_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

// Gamma draws at small shapes, for the sqrt-Gamma special case of rmhn and
// for Algorithm 3.
//
// R::rgamma(shape, scale) draws a Gamma(shape, 1) variate and multiplies it by
// the scale.  For a small shape that variate falls below the smallest normal
// double, DBL_MIN = 2.2e-308, with probability DBL_MIN^shape / Gamma(shape + 1)
// to first order, and it then comes back subnormal, with only some of its
// significant digits, or as exactly 0: at shape 0.005, 2.9% of the draws fall
// below DBL_MIN and 2.4% of all draws are exactly 0.  A product below DBL_MIN
// is equally imprecise when the scale is below 1.  Both are the draws below
//
//     t_min = DBL_MIN max(1, scale).
//
// Given T < t_min, T has density proportional to t^(shape-1) exp(-rate t) on
// (0, t_min), rate = 1 / scale.  It is drawn as t_min S, with S = U^(1/shape)
// and U uniform on (0, 1), which has density proportional to t^(shape-1)
// there, and is accepted with probability exp(-c S), c = rate t_min; in
// logarithms it keeps every digit however far below the double range it lies.
// c is DBL_MIN max(1, rate), at most DBL_MIN DBL_MAX, about 4, so the
// acceptance is at least exp(-4).  For a rate below 5e291, c is below the unit
// round-off and exp(-c S) is 1 to double precision, so the test is skipped
// there.  Draws at or above t_min are used as they come; each replacement
// consumes one extra uniform (two per trial above that rate), so the random
// stream moves only after a draw that would otherwise have lost its digits.
// Sun et al. (2023, Section 4.2) point out the same limitation of rgamma at
// small shapes, citing Liu, Martin and Syring (2017).
inline double gamma_precision_floor(double scale) {
  return std::numeric_limits<double>::min() * std::max(1.0, scale);
}

// log T for T drawn from its conditional law below t_min (see above), for a
// Gamma law of the given shape and rate.  Defined out of line, in
// mhn_sun.cpp, so that this rarely taken branch does not weigh on the
// inlining of the samplers that call it.
double gamma_tail_log_draw(double t_min, double shape, double rate);

// Gamma draws at large shapes: the largest shape at which a sampler of this
// package takes its Gamma variate from R::rgamma, 2^59 (about 5.8e17).
//
// For a shape a >= 1, R::rgamma (Ahrens and Dieter's algorithm GD) forms its
// unit-scale variate as (s + t/2)^2, t a standard normal deviate, from
// s = sqrt(a - 1/2) rounded to double precision; the rest of the algorithm
// only decides whether that value is accepted.  s + t/2 has standard
// deviation 1/2, and the rounding of s, up to half a unit in its last place,
// or 2^-53 s, moves every draw alike: by up to 2^-52 sqrt(a) standard
// deviations.  (a - 1/2 itself is exact below a = 2^52 and off by at most
// 1/2 above, which moves the draws by at most 1/(2 sqrt(a)) standard
// deviations, below 1e-8.)  R::rgamma then multiplies that variate by the
// scale 1 / rate, which is rounded once, by a relative 2^-53 at most, and so
// scales every draw alike: that moves them by up to 2^-53 sqrt(a) standard
// deviations.  In Algorithm 3 the same rounding makes the rate of the
// proposal differ from the rate in the acceptance test, which tilts the
// accepted law by as much.  Together these shifts are below
// 3 * 2^-53 sqrt(a) standard deviations, plus at most 1e-8 from a - 1/2:
// about 2.5e-7 at a = 2^59, the whole of the common shift of the
// square-root-Gamma closed form.  Algorithm 3's own setup adds a common tilt
// of the same order: it rounds the proposal's shape alpha r, the power r and
// the rate, and forms the scale m as u / sqrt(beta), so that its proposal and
// its acceptance test no longer match exactly.  The accepted law is then the
// target times exp(D (x / m - 1)), where D, 0 in exact arithmetic, is of
// order 2^-53 alpha, and x / m - 1 has a standard deviation of about
// 1 / sqrt(alpha): again a shift of every draw alike, of order
// 2^-53 sqrt(a) standard deviations.  Over shapes from 2^58 to 2^59 the
// total common shift of Algorithm 3 measured up to about
// 4.5 * 2^-53 sqrt(a) standard deviations (3.6e-7, at a = 5.2e17), a
// coefficient that gives 3.8e-7 at a = 2^59.  Each total is a few times the
// 1.2e-7 that Algorithm 1 allows for the rounding of its Normal proposal's
// mean (see build_sun_algo1), and grows with sqrt(a), to a whole standard
// deviation near a = 2^103, close to the shape 2^104 at which the Gamma law
// becomes as narrow as the spacing of the doubles at its mean.
// The remaining roundings, those of each draw itself, each of a unit in the
// last place or less, do not act alike on every draw.  Beyond this shape the
// square-root-Gamma closed form of rmhn and Algorithm 3 are not used (see
// rebuild_cache in mhn_rmhn.cpp and build_sun_algo3); the sqrt-Gamma
// proposal of Algorithm 1 never reaches it, as it is used only up to
// alpha = 1e18, a shape of 5e17.
constexpr double kRgammaShapeMax = 576460752303423488.0;  // 2^59

// Setup for Sun et al. (2023) Algorithm 1 (gamma > 0, alpha > 1).
// Theorem 1 of that paper gives the envelope constants of its two proposals,
// and Theorem 2e bounds the acceptance below by 0.8 for alpha >= 4 when the
// better of the two is taken.  Above alpha = 1e18 the Normal proposal is
// taken whatever the constants (see build_sun_algo1), and under a small tilt
// it accepts about 1/sqrt(2) of its proposals (0.706 at alpha = 2e18,
// gamma / sqrt(beta) up to 1e6), below that bound; the law drawn is exact.
//
// The fields from mu_opt to sg_centre describe the sampler for Y = y_per_x X.
// y_per_x is 1 where the setup is formed with beta and gamma directly.
// Outside the bounds build_sun_algo1 checks for that, it is the power of two
// 2^k with beta / 4^k in [1, 4), and the setup is built for Y, which is
// MHN(alpha, beta / 4^k, gamma / 2^k).  Where even that setup would leave
// the double range, y_per_x is 1 again and the fields describe the Normal
// proposal centred at the mode, in X (see build_sun_algo1).
struct SunAlgo1Setup {
  enum Proposal { NORMAL, SQRT_GAMMA } chosen;
  double alpha, beta, gamma;
  double y_per_x;
  double mu_opt;          // optimal Normal-proposal mean
  // The optimal Normal-proposal mean is mu_opt + mu_lo, mu_lo being the part
  // that the double mu_opt leaves out.  mu_lo is used only where
  // offset_normal is set, and there the proposal's offset from that mean
  // enters the acceptance test directly (see sample_sun_algo1).
  double mu_lo;
  bool offset_normal;
  double delta_opt;       // optimal sqrt-Gamma-proposal rate
  double sigma;           // = 1 / sqrt(2*beta), Normal-proposal SD
  double slope;           // = 2 beta mu_opt - gamma, the tangency slope
  // Whether the Normal acceptance test is evaluated in its cancellation-free
  // form (large alpha); see sample_sun_algo1.
  bool stable_normal;
  // The sqrt-Gamma acceptance exponent -(beta - delta) X^2 + gamma X
  // - gamma^2 / (4 (beta - delta)) is the perfect square
  // -(beta - delta) (X - c)^2 with c = gamma / (2 (beta - delta)).
  double sg_curv;         // = beta - delta_opt
  double sg_centre;       // = c
  // <0 -> NORMAL, >=0 -> SQRT_GAMMA, except that the Normal proposal is
  // taken at any sign above alpha = 1e18.  NaN, as is delta_opt, where
  // build_sun_algo1 takes the Normal proposal without the comparison.
  double log_K1_minus_K2;
};

// Setup for Sun et al. (2023) Algorithm 3 (gamma <= 0, alpha > 0).
// Theorem 3 of that paper gives the envelope for any matching point m > 0
// and Theorem 4 (alpha > 1) its acceptance probability; the 1/sqrt(2)
// bound of Theorem 4c is proved at m_opt and m_star, not at the one-step
// m used here, and its proof in the Supplementary Material treats
// alpha >= 3 only.
struct SunAlgo3Setup {
  double alpha, beta, gamma_abs;
  double m;                  // matching point (after Newton refinement)
  double r;                  // = (beta*m + |gamma|) / (2*beta*m + |gamma|)
  double shape;              // = alpha * r
  // = m * (beta*m + |gamma|): both the rate of the Gamma proposal and the
  // coefficient of T in the acceptance test.
  double rate;
  // Below this value a Gamma(shape, rate) draw from R::rgamma has lost its
  // precision or underflowed to zero; see gamma_precision_floor.
  double t_min;
  // beta m / (beta m + |gamma|) and |gamma| / (beta m + |gamma|), and whether
  // the acceptance test is evaluated in its cancellation-free form (large
  // rate); see a3_log_acc_stable in mhn_sun.cpp.
  double acc_p, acc_q;
  bool stable_test;
  double m_init;             // diagnostic: m before Newton refinement
  bool used_inflex_heuristic;// true when alpha > 1.1 elaborate path is taken
};

// Build the Algorithm 1 setup.  Throws Rcpp::stop on alpha <= 1 or gamma <= 0,
// and where the mode of the law lies beyond the largest double (see
// build_sun_algo1).
//
// With normal_only set, the setup is that of the Normal proposal alone,
// N(m, 1/(2 beta)) centred at the mode m, without the comparison with the
// sqrt-Gamma proposal, and gamma may have either sign.  The Normal envelope
// needs only alpha > 1 and m, the positive root of
// 2 beta x^2 - gamma x - (alpha - 1) = 0: the log ratio of the target to the
// proposal is then (alpha - 1) [log(x/m) - (x - m)/m] plus a constant, at
// most 0 and equal to it at x = m, whatever the sign of gamma.  rmhn uses it
// for a negligible tilt at a shape too large for R's Gamma generator (see
// kRgammaShapeMax), where its acceptance rate is about 1/sqrt(2).
SunAlgo1Setup build_sun_algo1(double alpha, double beta, double gamma,
                              bool normal_only = false);

// Draw one sample using the Algorithm 1 setup.  retries_out, if non-null,
// is incremented by the number of accept/reject retries.  After max_retries
// rejected proposals the draw is returned as NaN.
double sample_sun_algo1(const SunAlgo1Setup& s, int* retries_out = nullptr,
                        int max_retries = 1000);

// Build the Algorithm 3 setup.  Throws Rcpp::stop on gamma > 0, and for a
// shape above kRgammaShapeMax, beyond which its Gamma proposal cannot be
// drawn precisely enough (see build_sun_algo3).
SunAlgo3Setup build_sun_algo3(double alpha, double beta, double gamma);

// Draw one sample using the Algorithm 3 setup, as sample_sun_algo1 does.
double sample_sun_algo3(const SunAlgo3Setup& s, int* retries_out = nullptr,
                        int max_retries = 1000);

// Stops the call with an error where the law of MHN(alpha, beta, gamma) lies
// at or beyond the largest double.  Its standard deviation, at most
// 1 / sqrt(2 beta) <= 3.2e161 wherever alpha >= 1, is then below 1e-140 times
// its mode, so that no double precision draw could represent it.  Used by
// the Algorithm 1 setup and by the closed forms of rmhn.
[[noreturn]] void stop_beyond_range(double alpha, double beta, double gamma);

}  // namespace mhn

#endif  // MHN_SUN_H
