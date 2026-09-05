// Relaxed Transformed Density Rejection (RTDR) sampler for MHN, from
// Gao & Wang (2025).  Envelope construction is region-dependent; the
// region enum below documents the (alpha, gamma) split inherited from
// that paper.  Drawing uses a piecewise log-linear envelope evaluated
// at runtime by sample_rtdr().

#ifndef MHN_RTDR_H
#define MHN_RTDR_H

#include <Rcpp.h>
#include <vector>

namespace mhn {

// Region classification for the Gao & Wang (2025) RTDR algorithm.
// The (alpha, gamma) cases follow Section 4 of that paper.
enum RtdrRegion {
  REGION_A   = 0,   // alpha > 1, log-concave on f(x); T_0 tangent hat
  REGION_BC  = 1,   // alpha <= 1 apart from the region D corner, on g(y):
                    //   T_0 tangent hat if gamma <= 0 (log-concave), else
                    //   T_{-1/2} tangent hat (Thm 3.2)
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
//
// Inverse-square (T_{-1/2}) tangent pieces used by the alpha < 1, gamma > 0
// envelope (Gao & Wang 2025, Section 3.2).  Here `a` is the contact point t,
// `slope` is the half log-derivative L'(t)/2, and the hat is
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
  PIECE_TNEGHALF_RIGHT
};

// One envelope piece.  See PieceType for the geometry encoded by `type`.
struct EnvelopePiece {
  PieceType type;
  double a, b;             // bounds/anchors; b unused for exp tails, = ppl/ppr for T_{-1/2}
  double slope;            // exp rate (log-linear) or half log-derivative (T_{-1/2})
  double base_log_dens;    // log h at the anchor point a
  double log_area;         // log of the piece's contribution to the total area
  double aux = 0.0;        // T_{-1/2} sampling scale (om1/om3); unused otherwise
};

// Piecewise envelope for RTDR sampling.  Built once per (alpha, beta, gamma)
// triple by build_rtdr_envelope() and consumed repeatedly by sample_rtdr().
// All envelope quantities are stored in the beta=1 normalized coordinate;
// sqrt_beta is the scale factor restored on output (X_real = X_norm / sqrt_beta).
struct RtdrEnvelope {
  RtdrRegion region;
  double alpha;           // beta=1 normalized
  double gamma_norm;      // = gamma / sqrt(beta)
  double sqrt_beta;       // scale-restoration factor

  // Common envelope quantities, filled by the region-specific setups.
  double mode = 0.0;
  double log_dens_mode = 0.0;
  double t_l = 0.0, t_r = 0.0;
  double p_l = 0.0, p_r = 0.0;
  double slope_l = 0.0, slope_r = 0.0;
  bool simplified = false;       // region a only

  // Region D specific.
  int K_eff = 0;
  double y_star = 0.0;
  double rho = 0.0;
  std::vector<double> y_break;
  std::vector<double> log_dens_break;
  std::vector<double> alpha_k;
  bool has_left_tangent_d = false;
  // Set when a region-D setup could not place a single secant breakpoint and
  // handed the triple to the region-BC construction instead.  Both envelopes
  // are valid for the same density, so this is recorded rather than warned
  // about; it is exposed for inspection because the theory says it should
  // never happen.
  bool fell_back_to_bc = false;

  // Piece table built by the region-specific setups and consumed by
  // sample_rtdr.  Pieces appear left-to-right by support range.
  // piece_log_area mirrors pieces[i].log_area and exists so that the
  // diagnostic dump, and the tests that read it, can inspect the areas
  // without unpacking the piece table; the sampling path does not use it,
  // selecting instead through piece_cum_area / piece_area_total below.
  std::vector<EnvelopePiece> pieces;
  std::vector<double> piece_log_area;

  // Running sum of the piece areas, normalised by the largest log-area so the
  // exponentials stay in range.  Selecting a piece needs only a uniform draw
  // scaled by piece_area_total and a search in this table, so the per-draw
  // cost is a search rather than a fresh pass over every piece.  Region D can
  // hold thousands of pieces, which is where that matters.
  std::vector<double> piece_cum_area;
  double piece_area_total = 0.0;
};

// Build the envelope for the given (alpha, beta, gamma) triple.
RtdrEnvelope build_rtdr_envelope(double alpha, double beta, double gamma);

// Draw one sample using the supplied envelope.  retries_out, if non-null,
// is incremented by the number of accept/reject retries.  Implementation
// is in mhn_rtdr.cpp.
double sample_rtdr(const RtdrEnvelope& env, int* retries_out = nullptr);

}  // namespace mhn

#endif  // MHN_RTDR_H
