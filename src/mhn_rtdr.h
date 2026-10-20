// Relaxed Transformed Density Rejection (RTDR) sampler for MHN, from
// Gao & Wang (2025).  Envelope construction is region-dependent; the
// region enum below documents the (alpha, gamma) split inherited from
// that paper.  Drawing uses a piecewise envelope of exponential and
// inverse-square (T_{-1/2}) pieces, sampled by sample_rtdr().

#ifndef MHN_RTDR_H
#define MHN_RTDR_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace mhn {

// Region classification for the Gao & Wang (2025) RTDR algorithm.
// Region A follows Section 3.1 (Theorem 3.1) of that paper, region BC
// Section 3.2 (Theorem 3.2) and region D Section 4 (Theorem 4.4).
enum RtdrRegion {
  REGION_A   = 0,   // alpha > 1, log-concave on f(x); T_0 tangent hat
  REGION_BC  = 1,   // alpha <= 1 apart from the region D corner, on g(y):
                    //   T_0 tangent hat if gamma <= 0 (log-concave), else
                    //   T_{-1/2} tangent hat (Theorem 3.2)
  REGION_D   = 2    // alpha < 1/2 and gamma > 2(1 - sqrt(1 - 2*alpha)),
                    //   inflection-point envelope
};

// Piece type for the piecewise envelope.
//
// Log-linear (T_0 / exponential) pieces, h(x) = exp(base_log_dens + slope*(x - a)):
//   PIECE_EXP_RIGHT (slope < 0):    semi-infinite tail [a, +inf)
//   PIECE_EXP_LEFT  (slope > 0):    semi-infinite tail (-inf, a]
//   PIECE_EXP_BOUNDED (any slope):  bounded interval [a, b], h(a) = base_log_dens
//   PIECE_PLATEAU   (slope unused): bounded interval [a, b], h(x) = base_log_dens
//   PIECE_SECANT (any slope):       bounded interval [a, b], same form as EXP_BOUNDED
//                                   (kept distinct for clarity in region D code)
//   PIECE_EXP_BOUNDED_RIGHT (slope > 0): bounded interval [a, b] anchored at
//                                   its right end,
//                                   h(x) = exp(base_log_dens + slope*(x - b)),
//                                   with `aux` = expm1(-slope*(b - a))
//
// A left tangent whose left end lies far from the mode (x = 0 in region A, the
// last secant breakpoint in region D) carries its mass within a few 1/slope of
// b.  Anchored at a, both base_log_dens and slope*(x - a) are of size
// slope*(b - a) and cancel to the O(1) value of log h, which then carries an
// error of order ulp(slope*(b - a)); anchored at b, the same line has no
// cancellation where its mass is.
//
// Inverse-square (T_{-1/2}) tangent pieces used by the region BC envelope
// when gamma > 0 (Gao & Wang 2025, Section 3.2); region D uses none.  Here
// `a` is the contact point t, `slope` is the half log-derivative L'(t)/2, and
// the hat is
//   h(x) = exp(base_log_dens) / (slope*(x - a) - 1)^2,
// so log h(x) = base_log_dens - 2*log|slope*(x - a) - 1|.  `b` holds the
// sampling anchor (ppl / ppr) and `aux` the sampling scale (om1 / om3); the
// left piece covers (-inf, pl] and the right piece [pr, +inf).
//   PIECE_TNEGHALF_LEFT  (slope > 0): (-inf, pl],  x = b - aux/u
//   PIECE_TNEGHALF_RIGHT (slope < 0): [pr, +inf),  x = b + aux/u
enum PieceType {
  PIECE_PLATEAU,
  PIECE_EXP_RIGHT,
  PIECE_EXP_LEFT,
  PIECE_EXP_BOUNDED,
  PIECE_SECANT,
  PIECE_TNEGHALF_LEFT,
  PIECE_TNEGHALF_RIGHT,
  PIECE_EXP_BOUNDED_RIGHT
};

// One envelope piece.  See PieceType for the geometry encoded by `type`.
struct EnvelopePiece {
  PieceType type;
  double a, b;             // bounds/anchors; b unused for exp tails, = ppl/ppr for T_{-1/2}
  double slope;            // exp rate (log-linear) or half log-derivative (T_{-1/2})
  double base_log_dens;    // log h at the anchor: b for EXP_BOUNDED_RIGHT, else a
  double log_area;         // log of the piece's contribution to the total area
  double aux = 0.0;        // T_{-1/2} sampling scale (om1/om3), or
                           // expm1(-slope*(b - a)) for EXP_BOUNDED_RIGHT
};

// Piecewise envelope for RTDR sampling.  Built once per (alpha, beta, gamma)
// triple by build_rtdr_envelope() and consumed repeatedly by sample_rtdr().
// All envelope quantities are stored in the beta=1 normalized coordinate
// x_norm = sqrt(beta) x.  Off the centred path (below) a draw is returned as
// x_norm / sqrt_beta.
struct RtdrEnvelope {
  RtdrRegion region;
  double alpha;           // beta=1 normalized
  double gamma_norm;      // = gamma / sqrt(beta)
  double sqrt_beta;       // scale-restoration factor
  double beta = 1.0;      // the caller's beta and gamma, as passed
  double gamma = 0.0;

  // Common envelope quantities, filled by the region-specific setups.  Every
  // position -- mode, t_l, t_r, p_l, p_r, y_star, y_break and the piece
  // bounds -- is in the sampling coordinate s = (x or y) - origin.
  double mode = 0.0;
  double log_dens_mode = 0.0;
  // When the raw log-density at the mode is so large that an O(1) difference
  // of two of its values has no significant digits left, the setup measures
  // the ordinate from the mode instead: log_dens_mode is then 0 and every
  // base_log_dens / log_area is relative to the peak.  sample_rtdr must
  // measure the target the same way.
  bool centred = false;
  // On the centred path the abscissa is measured from the mode as well:
  // origin is the double mode on the region's own axis, m in region A and
  // m_g = log(u_mode) in regions BC and D, so the sampling coordinate s is
  // resolved to full relative precision near the mode.  A relative error of
  // order 1e-16 in the double mode changes the law of s only by as much; its
  // location is supplied separately, as the mode in the caller's coordinate
  // carried to about twice double precision in x_mode_hi + x_mode_lo: the
  // mode of X in region A, and e^{m_g}/sqrt(beta) in regions BC and D, where
  // x = e^{m_g + s}/sqrt(beta).  A draw is returned as
  // x_mode_hi + (x_mode_lo + s/sqrt_beta) in region A and as
  // x_mode_hi + (x_mode_lo + x_mode_hi expm1(s)) in regions BC and D, the
  // exact draw rounded once.  Elsewhere origin and x_mode_hi/lo are 0.
  double origin = 0.0;
  double u_mode = 0.0;      // e^y at the mode (regions BC and D), beta = 1
  double x_mode_hi = 0.0, x_mode_lo = 0.0;
  double t_l = 0.0, t_r = 0.0;
  double p_l = 0.0, p_r = 0.0;
  double slope_l = 0.0, slope_r = 0.0;
  bool simplified = false;       // region a only
  // Region A only: the contact points were searched from next to the mode
  // (see rtdr_region_a_centred below) rather than walked in.
  bool quick_start = false;

  // Region D specific.
  // Number of secant breakpoints built (= y_break.size()): ceil(rho) for
  // envelope (14) of Gao & Wang (2025); for envelope (15) the package, like
  // the R code of their Appendix B, stops the chain at the first breakpoint
  // on or below the left tangent (Eq. (15) as printed keeps all ceil(rho),
  // ending at the dual point t_dot_l).  At most ceil(rho).
  int K_eff = 0;
  double y_star = 0.0;
  double rho = 0.0;
  std::vector<double> y_break;
  std::vector<double> log_dens_break;
  std::vector<double> alpha_k;
  bool has_left_tangent_d = false;
  // Piece table built by the region-specific setups and consumed by
  // sample_rtdr.  Pieces appear left-to-right by support range.
  // piece_log_area mirrors pieces[i].log_area and exists so that the
  // diagnostic dump, and the tests that read it, can inspect the areas
  // without unpacking the piece table; the sampling path does not use it,
  // selecting instead through piece_cum_area / piece_area_total below.
  std::vector<EnvelopePiece> pieces;
  std::vector<double> piece_log_area;

  // Running sum of the piece areas, normalised by the largest log-area so the
  // exponentials stay in range, built once per envelope so that selecting a
  // piece costs one uniform and a binary search rather than a pass over every
  // piece.
  std::vector<double> piece_cum_area;
  double piece_area_total = 0.0;
};

// Build the envelope for the given (alpha, beta, gamma) triple.
RtdrEnvelope build_rtdr_envelope(double alpha, double beta, double gamma);

// Where region A's setup (alpha > 1) starts its contact-point search, which
// sets most of its cost.  It starts next to the mode, at m +- inc (the start
// Gao & Wang 2025, Section 2, recommend), in two places: on the centred path,
// where the magnitudes of the three terms of the raw log-density at the mode,
// (alpha-1) log m, m^2 and gamma_norm m, sum past 0.01 / DBL_EPSILON, and at
// a strong negative tilt, gamma_norm < -1e10 sqrt(alpha).
// Elsewhere it walks the two points in from m/2 and m + max(1, m/2), at a
// cost that grows with log2 of the mode.  setup_region_a makes both tests
// with the same code, and the default method of rmhn (rebuild_cache in
// mhn_rmhn.cpp) calls these predicates to know the cost of RTDR's setup, so
// the two cannot disagree.  gamma_norm = gamma / sqrt(beta), as in
// RtdrEnvelope.  For -1e10 sqrt(alpha) <= gamma_norm <= 2 sqrt(alpha) the
// centred path starts at a shape of about 2.35e12 at the earliest (at
// gamma_norm = 2 sqrt(alpha)), and below 1e12 that sum stays under 0.42 of
// the switch.  It is at least alpha - 2, so from a shape of about 4.5e13 the
// path is centred at every tilt.
bool rtdr_region_a_centred(double alpha, double gamma_norm);
inline bool rtdr_region_a_strong_negative_tilt(double alpha, double gamma_norm) {
  return gamma_norm < -1e10 * std::sqrt(alpha);
}

// Draw one sample using the supplied envelope.  retries_out, if non-null,
// is incremented by the number of accept/reject retries.  After max_retries
// rejected proposals the draw is returned as NaN.  Implementation is in
// mhn_rtdr.cpp.
double sample_rtdr(const RtdrEnvelope& env, int* retries_out = nullptr,
                   int max_retries = 1000);

}  // namespace mhn

#endif  // MHN_RTDR_H
