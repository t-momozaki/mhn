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
  REGION_A   = 0,   // alpha >= 1, log-concave on f(x)
  REGION_BC  = 1,   // alpha < 1, T_{-1/2}-concave on g(y)
  REGION_D   = 2    // alpha < 1/2 and gamma > 2(1 - sqrt(1 - 2*alpha)),
                    //   inflection-point envelope
};

// Piece type for the piecewise envelope.  Sign convention on `slope`:
//   PIECE_EXP_RIGHT (slope < 0):    semi-infinite tail [a, +inf)
//   PIECE_EXP_LEFT  (slope > 0):    semi-infinite tail (-inf, a]
//   PIECE_EXP_BOUNDED (any slope):  bounded interval [a, b], h(a) = base_log_dens
//   PIECE_PLATEAU   (slope unused): bounded interval [a, b], h(x) = base_log_dens
//   PIECE_SECANT (any slope):       bounded interval [a, b], same form as EXP_BOUNDED
//                                   (kept distinct for clarity in region D code)
enum PieceType {
  PIECE_PLATEAU,
  PIECE_EXP_RIGHT,
  PIECE_EXP_LEFT,
  PIECE_EXP_BOUNDED,
  PIECE_SECANT
};

// One envelope piece.  See PieceType for the geometry encoded by `type`.
struct EnvelopePiece {
  PieceType type;
  double a, b;             // bounds; b is unused for semi-infinite tails
  double slope;            // exponential rate, signed per PieceType
  double base_log_dens;    // log h at the anchor point a
  double log_area;         // log of the piece's contribution to the total area
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

  // Common envelope quantities (filled by region-specific setups in Step 3.2).
  double mode = 0.0;
  double log_dens_mode = 0.0;
  double t_l = 0.0, t_r = 0.0;
  double p_l = 0.0, p_r = 0.0;
  double slope_l = 0.0, slope_r = 0.0;
  bool simplified = false;       // region a only

  // Region D specific (filled in Step 3.2).
  int K_eff = 0;
  double y_star = 0.0;
  double rho = 0.0;
  std::vector<double> y_break;
  std::vector<double> log_dens_break;
  std::vector<double> alpha_k;
  bool has_left_tangent_d = false;

  // Piece table built by the region-specific setups and consumed by
  // sample_rtdr.  Pieces appear left-to-right by support range; piece_log_area
  // duplicates the log-areas for convenient access (and is what dispatchers
  // typically need).  The two stay in lockstep (same indexing).
  std::vector<EnvelopePiece> pieces;
  std::vector<double> piece_log_area;
};

// Build the envelope for the given (alpha, beta, gamma) triple.
// Step 3.1 ships a stub that throws; region setups are filled in Step 3.2.
RtdrEnvelope build_rtdr_envelope(double alpha, double beta, double gamma);

// Draw one sample using the supplied envelope.  retries_out, if non-null,
// is incremented by the number of accept/reject retries.  Implementation
// in Step 3.2.
double sample_rtdr(const RtdrEnvelope& env, int* retries_out = nullptr);

}  // namespace mhn

#endif  // MHN_RTDR_H
