// Relaxed Transformed Density Rejection (RTDR) sampler for the
// Modified Half-Normal distribution, following Gao & Wang (2025).
// This file holds the common infrastructure (Newton contact-point
// search, piecewise envelope sampler) together with the
// region-specific setups (a, bc, d) and the dispatcher.

#include "mhn_rtdr.h"
#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_log_arith.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace {

// EnvelopePiece and PieceType live in namespace mhn (see mhn_rtdr.h) so
// that RtdrEnvelope can hold a std::vector<EnvelopePiece> directly.
using mhn::EnvelopePiece;
using mhn::PieceType;
using mhn::PIECE_PLATEAU;
using mhn::PIECE_EXP_RIGHT;
using mhn::PIECE_EXP_LEFT;
using mhn::PIECE_EXP_BOUNDED;
using mhn::PIECE_SECANT;

// ====================================================================
// Section 2: contact_point_newton — Eq. (8) iteration.
// ====================================================================

// Iteratively solve for a contact point t such that
//     log f(m) - log f(t) = delta
// where delta = 1 (log-concave) or delta = log 4 (T_{-1/2}-concave).
//
// Acceptable-range termination: returns early if the increment
// log f(m) - log f(t) lies in [accept_low, accept_high] (caller-supplied,
// e.g. [0.46, 2.49] for log-concave or [0.93, 1.99] for T_{-1/2}).
//
// Optional bounds (t_min, t_max): each Newton step is clamped to
// halfway-bisected position if it would exit the [t_min, t_max] interval.
// Defaults are ±inf.
//
// On non-convergence within max_iter, returns the current t (Region setups
// are responsible for downstream sanity checks like log_area > 0).
double contact_point_newton(double t_init,
                            double log_dens_mode_val,
                            double delta,
                            const std::function<double(double)>& log_dens_at,
                            const std::function<double(double)>& log_dens_deriv_at,
                            double accept_low,
                            double accept_high,
                            int max_iter = 30,
                            double tol = 1e-10,
                            double t_min = -std::numeric_limits<double>::infinity(),
                            double t_max =  std::numeric_limits<double>::infinity()) {
  double t = t_init;
  for (int iter = 0; iter < max_iter; ++iter) {
    const double ld = log_dens_at(t);
    const double increment = log_dens_mode_val - ld;
    if (increment >= accept_low && increment <= accept_high) return t;
    const double deriv = log_dens_deriv_at(t);
    if (!std::isfinite(deriv) || std::abs(deriv) < tol) break;
    double t_new = t + (increment - delta) / deriv;
    if (!std::isfinite(t_new)) break;
    // Clamp via bisection toward the violated bound to keep t in (t_min, t_max).
    if (t_new <= t_min) t_new = 0.5 * (t + t_min);
    if (t_new >= t_max) t_new = 0.5 * (t + t_max);
    if (std::abs(t_new - t) < tol) return t_new;
    t = t_new;
  }
  return t;
}

// ====================================================================
// Section 3: piecewise_envelope_sample — piece selection + per-piece
//   inverse-CDF sampling.
// ====================================================================

// Pick a piece index proportional to exp(piece.log_area).
// Uses log-space stable normalization (subtract max before exp).
int select_piece(const std::vector<EnvelopePiece>& pieces) {
  const std::size_t K = pieces.size();
  double m = -std::numeric_limits<double>::infinity();
  for (const auto& p : pieces) if (p.log_area > m) m = p.log_area;
  double total = 0.0;
  std::vector<double> cum(K);
  for (std::size_t i = 0; i < K; ++i) {
    total += std::exp(pieces[i].log_area - m);
    cum[i] = total;
  }
  const double u = R::runif(0.0, 1.0) * total;
  for (std::size_t i = 0; i < K; ++i) {
    if (u <= cum[i]) return static_cast<int>(i);
  }
  return static_cast<int>(K - 1);
}

// Inverse-CDF sampling within one piece.  Returns x (or y, depending on
// the region's coordinate system — caller handles the mapping).
double sample_within_piece(const EnvelopePiece& piece) {
  const double u = R::runif(0.0, 1.0);
  switch (piece.type) {
    case PIECE_PLATEAU:
      return piece.a + u * (piece.b - piece.a);
    case PIECE_EXP_RIGHT: {
      // h(x) = h(a) * exp(slope * (x - a)),  slope < 0,  x >= a
      // F(x) ∝ 1 - exp(slope * (x - a)),  x = a + log(1 - u) / slope
      return piece.a + std::log1p(-u) / piece.slope;
    }
    case PIECE_EXP_LEFT: {
      // h(x) = h(a) * exp(slope * (x - a)),  slope > 0,  x <= a
      // F(x) ∝ exp(slope * (x - a)),  x = a + log(u) / slope
      return piece.a + std::log(u) / piece.slope;
    }
    case PIECE_EXP_BOUNDED:
    case PIECE_SECANT: {
      // h(x) = h(a) * exp(slope * (x - a)) on [a, b].  Same form for both
      // EXP_BOUNDED and SECANT; the type tag is informational (region D
      // distinguishes secants from tangents in the piece table).
      // F^-1: x = a + log(1 + u * (exp(slope*(b-a)) - 1)) / slope.
      // Stable form via expm1 / log1p.
      const double width = piece.b - piece.a;
      const double z = std::expm1(piece.slope * width);
      return piece.a + std::log1p(u * z) / piece.slope;
    }
    default:
      Rcpp::stop("piecewise_envelope_sample: unknown PieceType.");
  }
}

// One draw from the piecewise envelope.  Caller owns the accept/reject
// decision (caller knows the target density on x).
double piecewise_envelope_sample(const std::vector<EnvelopePiece>& pieces) {
  const int idx = select_piece(pieces);
  return sample_within_piece(pieces[idx]);
}

// log h(x) on a single piece, evaluated at point x.  Used inside the
// accept/reject loop.
double log_piece_at(const EnvelopePiece& piece, double x) {
  switch (piece.type) {
    case PIECE_PLATEAU:
      return piece.base_log_dens;
    case PIECE_EXP_RIGHT:
    case PIECE_EXP_LEFT:
    case PIECE_EXP_BOUNDED:
    case PIECE_SECANT:
      // h(x) = exp(base_log_dens) * exp(slope * (x - a))
      return piece.base_log_dens + piece.slope * (x - piece.a);
  }
  return -std::numeric_limits<double>::infinity();
}

// Classify the (alpha, gamma_norm) point into one of the three RTDR regions.
mhn::RtdrRegion classify_region(double alpha, double gamma_norm) {
  if (alpha >= 1.0) return mhn::REGION_A;
  if (alpha >= 0.5) return mhn::REGION_BC;
  // alpha in (0, 1/2).  Threshold gamma_d = 2(1 - sqrt(1 - 2*alpha)).
  const double gamma_d = 2.0 * (1.0 - std::sqrt(1.0 - 2.0 * alpha));
  if (gamma_norm <= gamma_d) return mhn::REGION_BC;
  return mhn::REGION_D;
}

// Forward declarations of region setups (definitions below; stubbed in
// Step 3.2.1, filled in Steps 3.2.3 / 3.2.4 / 3.2.5).
void setup_region_a(mhn::RtdrEnvelope& env);
void setup_region_bc(mhn::RtdrEnvelope& env);
void setup_region_d(mhn::RtdrEnvelope& env);

// ====================================================================
// Region setup helpers — area/anchor formulas
// ====================================================================
// log h on [a, b] of an exponential piece h(x) = exp(base + slope * (x - a)):
//   integral over [a, b] = exp(base) * expm1(slope * (b - a)) / slope.
// In log-space, with width = b - a:
//   log_area = base + log(|expm1(slope * width)|) - log(|slope|).
// Both factors share the same sign so the absolute-value cancels into a
// real positive area.  We assume slope != 0 (caller guarantees).
double log_area_exp_bounded(double base, double slope, double width) {
  return base + std::log(std::abs(std::expm1(slope * width)))
              - std::log(std::abs(slope));
}

// Set base_log_dens and log_area on a piece, given anchor convention:
//   PIECE_PLATEAU:    base = mode_log_dens, log_area = base + log(b - a)
//   PIECE_EXP_BOUNDED [a, b]: anchor at piece.a; log_area via expm1 form
//   PIECE_EXP_LEFT  (-inf, a]: anchor at piece.a; log_area = base - log(slope)
//   PIECE_EXP_RIGHT [a, +inf): anchor at piece.a; log_area = base - log(-slope)
EnvelopePiece make_plateau(double a, double b, double base_log_dens) {
  EnvelopePiece p;
  p.type = PIECE_PLATEAU; p.a = a; p.b = b; p.slope = 0.0;
  p.base_log_dens = base_log_dens;
  p.log_area = base_log_dens + std::log(b - a);
  return p;
}
EnvelopePiece make_exp_bounded(double a, double b, double slope, double base_log_dens) {
  EnvelopePiece p;
  p.type = PIECE_EXP_BOUNDED; p.a = a; p.b = b; p.slope = slope;
  p.base_log_dens = base_log_dens;
  p.log_area = log_area_exp_bounded(base_log_dens, slope, b - a);
  return p;
}
EnvelopePiece make_exp_left(double a, double slope, double base_log_dens) {
  // (-inf, a], slope > 0
  EnvelopePiece p;
  p.type = PIECE_EXP_LEFT; p.a = a; p.b = 0.0; p.slope = slope;
  p.base_log_dens = base_log_dens;
  p.log_area = base_log_dens - std::log(slope);
  return p;
}
EnvelopePiece make_exp_right(double a, double slope, double base_log_dens) {
  // [a, +inf), slope < 0
  EnvelopePiece p;
  p.type = PIECE_EXP_RIGHT; p.a = a; p.b = 0.0; p.slope = slope;
  p.base_log_dens = base_log_dens;
  p.log_area = base_log_dens - std::log(-slope);
  return p;
}
EnvelopePiece make_secant(double a, double b, double slope, double base_log_dens) {
  // [a, b], slope free; same area formula as EXP_BOUNDED.
  EnvelopePiece p;
  p.type = PIECE_SECANT; p.a = a; p.b = b; p.slope = slope;
  p.base_log_dens = base_log_dens;
  p.log_area = log_area_exp_bounded(base_log_dens, slope, b - a);
  return p;
}

void finalize_pieces(mhn::RtdrEnvelope& env, std::vector<EnvelopePiece>& pieces) {
  env.pieces = pieces;
  env.piece_log_area.clear();
  env.piece_log_area.reserve(pieces.size());
  for (const auto& p : pieces) env.piece_log_area.push_back(p.log_area);
}

// ====================================================================
// Region A setup [alpha >= 1, log-concave on f(x)]
// ====================================================================
void setup_region_a(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // mode of log f(x) = (alpha-1) log x - x^2 + gamma_norm * x
  const double m = (gn + std::sqrt(gn * gn + 8.0 * (alpha - 1.0))) / 4.0;
  env.mode = m;

  auto log_f = [alpha, gn](double x) -> double {
    if (x <= 0.0) return -std::numeric_limits<double>::infinity();
    return (alpha - 1.0) * std::log(x) - x * x + gn * x;
  };
  auto dlog_f = [alpha, gn](double x) -> double {
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    return (alpha - 1.0) / x - 2.0 * x + gn;
  };
  const double log_f_mode = log_f(m);
  env.log_dens_mode = log_f_mode;

  // Simplified-envelope condition from Gao & Wang (2025) Theorem 3.1.
  const double e_const = std::exp(1.0);
  bool simplified = false;
  if (alpha >= 1.0 && alpha <= e_const) {
    const double inner = std::sqrt((e_const - 1.0) * (alpha - 1.0)) - (alpha - 1.0);
    if (inner >= 0.0) {
      const double gamma_threshold = std::sqrt(8.0 * inner);
      if (gn <= gamma_threshold) simplified = true;
    }
  }
  env.simplified = simplified;

  // Right contact point: t_r > m, slope_r = L'(t_r) < 0.
  // Initial heuristic: m + max(1, m/2).
  const double t_r_init = m + std::max(1.0, 0.5 * m);
  env.t_r = ::contact_point_newton(t_r_init, log_f_mode, 1.0, log_f, dlog_f,
                                   0.46, 2.49, /*max_iter=*/30, /*tol=*/1e-10,
                                   /*t_min=*/m, /*t_max=*/std::numeric_limits<double>::infinity());
  env.slope_r = dlog_f(env.t_r);
  if (env.slope_r >= 0.0) {
    Rcpp::stop("setup_region_a: right contact has non-negative slope (Newton failure).");
  }
  env.p_r = env.t_r + (log_f_mode - log_f(env.t_r)) / env.slope_r;

  std::vector<EnvelopePiece> pieces;

  if (!simplified) {
    // Left contact point: 0 < t_l < m, slope_l = L'(t_l) > 0.
    // Initial: m / 2 (kept strictly positive).
    const double t_l_init = std::max(0.5 * m, 1e-6);
    env.t_l = ::contact_point_newton(t_l_init, log_f_mode, 1.0, log_f, dlog_f,
                                     0.46, 2.49, /*max_iter=*/30, /*tol=*/1e-10,
                                     /*t_min=*/0.0, /*t_max=*/m);
    if (env.t_l <= 0.0) env.t_l = 0.5 * m;  // defensive
    env.slope_l = dlog_f(env.t_l);
    if (env.slope_l <= 0.0) {
      Rcpp::stop("setup_region_a: left contact has non-positive slope (Newton failure).");
    }
    env.p_l = env.t_l + (log_f_mode - log_f(env.t_l)) / env.slope_l;
    if (env.p_l <= 0.0 || env.p_l >= env.p_r) {
      Rcpp::stop("setup_region_a: invalid intersection points (p_l out of range).");
    }
    // Left tangent: [0, p_l], h(p_l) = f(m), so base at piece.a = 0
    // is log_f_mode - slope_l * p_l.
    pieces.push_back(make_exp_bounded(
        0.0, env.p_l, env.slope_l, log_f_mode - env.slope_l * env.p_l));
  }

  // Plateau: [plateau_a, p_r], h = f(m).
  const double plateau_a = simplified ? 0.0 : env.p_l;
  pieces.push_back(make_plateau(plateau_a, env.p_r, log_f_mode));

  // Right tangent: [p_r, +inf), slope_r < 0, h(p_r) = f(m).
  pieces.push_back(make_exp_right(env.p_r, env.slope_r, log_f_mode));

  finalize_pieces(env, pieces);
}

// ====================================================================
// Region BC setup [alpha < 1, T_{-1/2}-concave on g(y) = exp(alpha*y
//   - exp(2y) + gamma_norm * exp(y)), y in R]
// ====================================================================
void setup_region_bc(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // Mode in u = exp(y): 2 u^2 - gamma_norm u - alpha = 0,
  //   u_+ = (gamma_norm + sqrt(gamma_norm^2 + 8 alpha)) / 4 > 0.
  const double u_mode = (gn + std::sqrt(gn * gn + 8.0 * alpha)) / 4.0;
  const double m_g = std::log(u_mode);
  env.mode = m_g;

  auto log_g = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha * y - ey * ey + gn * ey;
  };
  auto dlog_g = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha - 2.0 * ey * ey + gn * ey;
  };
  const double log_g_mode = log_g(m_g);
  env.log_dens_mode = log_g_mode;
  env.simplified = false;  // BC envelope is always 3-piece

  // Left contact point: y_l < m_g, slope_l > 0.
  const double t_l_init = m_g - 1.0;
  env.t_l = ::contact_point_newton(t_l_init, log_g_mode, std::log(4.0),
                                   log_g, dlog_g, 0.93, 1.99,
                                   /*max_iter=*/30, /*tol=*/1e-10,
                                   /*t_min=*/-std::numeric_limits<double>::infinity(),
                                   /*t_max=*/m_g);
  env.slope_l = dlog_g(env.t_l);
  if (env.slope_l <= 0.0) {
    Rcpp::stop("setup_region_bc: left contact has non-positive slope (Newton failure).");
  }
  env.p_l = env.t_l + (log_g_mode - log_g(env.t_l)) / env.slope_l;

  // Right contact point: y_r > m_g, slope_r < 0.
  const double t_r_init = m_g + 1.0;
  env.t_r = ::contact_point_newton(t_r_init, log_g_mode, std::log(4.0),
                                   log_g, dlog_g, 0.93, 1.99,
                                   /*max_iter=*/30, /*tol=*/1e-10,
                                   /*t_min=*/m_g,
                                   /*t_max=*/std::numeric_limits<double>::infinity());
  env.slope_r = dlog_g(env.t_r);
  if (env.slope_r >= 0.0) {
    Rcpp::stop("setup_region_bc: right contact has non-negative slope (Newton failure).");
  }
  env.p_r = env.t_r + (log_g_mode - log_g(env.t_r)) / env.slope_r;

  if (env.p_l >= env.p_r) {
    Rcpp::stop("setup_region_bc: invalid intersection points (p_l >= p_r).");
  }

  std::vector<EnvelopePiece> pieces;
  pieces.push_back(make_exp_left(env.p_l, env.slope_l, log_g_mode));
  pieces.push_back(make_plateau(env.p_l, env.p_r, log_g_mode));
  pieces.push_back(make_exp_right(env.p_r, env.slope_r, log_g_mode));
  finalize_pieces(env, pieces);
}

// ====================================================================
// Region D setup [alpha < 1/2, gamma > 2(1 - sqrt(1 - 2 alpha)),
//   inflection-point envelope on T_{-1/2}-transformed g(y)].
//
// Implements both envelope 14 and envelope 15 of Gao & Wang (2025), with
// automatic switching on the criterion log(g(m_g)/g(y*)) > 2.49.
//   - envelope 14 (default): y_hat_star = y*, plateau spans (y*, p_r].
//     Used when the log-density drop from mode to inflection is moderate.
//   - envelope 15: y_hat_star = t_hat_l, an extra left-tangent piece is
//     inserted between the secant chain and the plateau.  Used when the
//     mode-to-inflection drop exceeds 2.49 (Gao & Wang 2025,
//     Theorem 4.4 -- required for the uniform 1/e acceptance bound).
// ====================================================================
void setup_region_d(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  if (gn <= 0.0) Rcpp::stop("setup_region_d: gamma_norm must be positive");

  // Mode (in y) and inflection point (in y).
  const double u_mode = (gn + std::sqrt(gn * gn + 8.0 * alpha)) / 4.0;
  const double m_g = std::log(u_mode);
  const double y_star = std::log(gn / 4.0);
  env.mode = m_g;
  env.y_star = y_star;

  auto log_g = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha * y - ey * ey + gn * ey;
  };
  auto dlog_g = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha - 2.0 * ey * ey + gn * ey;
  };
  const double log_g_mode = log_g(m_g);
  const double log_g_y_star = log_g(y_star);
  env.log_dens_mode = log_g_mode;
  env.simplified = false;

  // Envelope selection.  log(g(m_g)/g(y_star)) > 2.49 triggers env 15.
  const bool use_env_15 = (log_g_mode - log_g_y_star) > 2.49;

  // Compute rho according to the selected envelope.  env 15 also derives
  // the contact point t_l (in the log-concave half [y_star, m_g]) and the
  // tangent intersection p_l for the extra left-tangent piece.
  double rho;
  double t_l_d = 0.0, t_hat_l = 0.0, slope_l_d = 0.0;
  double log_g_t_l = 0.0, p_l_d = 0.0;
  if (!use_env_15) {
    // Envelope 14: y_hat_star = y_star.  Closed form rho = 3 gn^2 / 16.
    rho = log_g_y_star - alpha * y_star;
  } else {
    // Envelope 15: contact-point search inside [y_star, m_g] (log-concave
    // sub-region of g).  Initial value at the mid-point — the m_g - 1
    // BC heuristic can drop below y_star for extreme parameters.
    const double y_0 = std::log(gn / 2.0);   // dual point of -inf
    const double t_l_init = 0.5 * (y_star + m_g);
    t_l_d = ::contact_point_newton(t_l_init, log_g_mode, std::log(4.0),
                                   log_g, dlog_g, 0.93, 1.99,
                                   30, 1e-10, y_star, m_g);
    if (!(t_l_d < y_0)) {
      Rcpp::stop("setup_region_d (env 15): t_l >= y_0 (cannot compute t_hat_l).");
    }
    t_hat_l = std::log(gn / 2.0 - std::exp(t_l_d));
    slope_l_d = dlog_g(t_l_d);
    log_g_t_l = log_g(t_l_d);
    if (slope_l_d <= alpha) {
      // Gao & Wang (2025) Lemma 4.5 requires L'(t_l) > alpha for the
      // dual-point construction.
      Rcpp::stop("setup_region_d (env 15): L'(t_l) <= alpha "
                 "(Gao & Wang 2025 Lemma 4.5 violated).");
    }
    p_l_d = t_l_d + (log_g_mode - log_g_t_l) / slope_l_d;
    // rho = L_g(t_hat_l) - alpha * t_hat_l
    //     = closed form gn^2/4 - exp(2 t_l_d), but we use the direct form.
    rho = log_g(t_hat_l) - alpha * t_hat_l;
  }
  env.rho = rho;
  if (rho <= 0.0) Rcpp::stop("setup_region_d: rho <= 0 (unexpected)");

  const int K = static_cast<int>(std::ceil(rho));

  // K_eff guard: largest k <= K with gn^2 > 4 k (rho/K).  Theoretically
  // both envelopes always satisfy the guard at every k (rho < gn^2 / 4 in
  // both forms), but we keep the loop for defensive programming.
  int K_eff = K;
  while (K_eff >= 1 && gn * gn <= 4.0 * static_cast<double>(K_eff) * (rho / K)) {
    --K_eff;
  }

  if (K_eff < 1) {
    Rcpp::warning("setup_region_d: K_eff < 1, falling back to region BC.");
    env.region = mhn::REGION_BC;
    setup_region_bc(env);
    return;
  }

  // Breakpoints y_k built incrementally with early termination.
  //
  // Gao & Wang (2025) Theorem 4.4 sets K = ceil(rho) as the upper
  // bound on the number of secant pieces.  For extreme (alpha << 1,
  // gamma >> 1) corners K can be tens of thousands -- e.g.
  // (alpha=0.3, gamma=10000) gives K ~ 12000 -- and constructing
  // every piece dominates the setup cost.
  //
  // Gao & Wang (2025) Lemma 4.1 (proof in Appendix A.6) justifies
  // terminating the secant sequence as soon as the next y_k drops
  // below the left tangent line at t_l: any further secants would be
  // dominated by the tangent envelope and contribute nothing.  We
  // mirror that condition here for envelope (15) of Gao & Wang (2025)
  // -- the dual-point + left-tangent construction selected by the
  // `use_env_15` branch above.
  // Envelope 14 has no left tangent (the leftmost piece is the
  // alpha-slope exponential tail) and K stays small in practice, so it
  // falls through with no early termination.
  //
  // The break decision is made AFTER the current k is appended, so K_eff
  // is at least 1 and at most the original ceil(rho).
  const double drho = rho / static_cast<double>(K);
  std::vector<double> y_b;
  std::vector<double> log_g_b;
  y_b.reserve(64);
  log_g_b.reserve(64);
  y_b.push_back(0.0);          // 1-indexed: y_b[0] unused
  log_g_b.push_back(0.0);
  int K_actual = 0;
  for (int k = 1; k <= K_eff; ++k) {
    const double rk = 4.0 * static_cast<double>(k) * drho;
    const double inside = gn * gn - rk;
    if (inside <= 0.0) {
      Rcpp::stop("setup_region_d: breakpoint inside non-positive (K_eff guard).");
    }
    const double yk = std::log(2.0 * static_cast<double>(k) * drho
                               / (std::sqrt(inside) + gn));
    const double log_g_yk = log_g(yk);
    y_b.push_back(yk);
    log_g_b.push_back(log_g_yk);
    K_actual = k;
    if (use_env_15) {
      // Gao & Wang (2025) Lemma 4.1 (proof in Appendix A.6): break
      // once log g(y_k) drops below the tangent line at t_l_d.  Beyond
      // this point further secants are redundant.
      const double tangent = log_g_t_l + slope_l_d * (yk - t_l_d);
      if (log_g_yk <= tangent) break;
    }
  }
  K_eff = K_actual;
  env.K_eff = K_eff;
  env.y_break.assign(y_b.begin() + 1, y_b.end());
  env.log_dens_break.assign(log_g_b.begin() + 1, log_g_b.end());

  // Secant slopes alpha_k for k=2..K_eff.
  std::vector<double> a_k(K_eff + 1);
  for (int k = 2; k <= K_eff; ++k) {
    a_k[k] = (log_g_b[k] - log_g_b[k-1]) / (y_b[k] - y_b[k-1]);
  }
  if (K_eff >= 2) {
    env.alpha_k.assign(a_k.begin() + 2, a_k.end());
  } else {
    env.alpha_k.clear();
  }

  // Right contact point in the log-concave region (same as BC right tangent).
  const double t_r_init = m_g + 1.0;
  env.t_r = ::contact_point_newton(t_r_init, log_g_mode, std::log(4.0),
                                   log_g, dlog_g, 0.93, 1.99,
                                   30, 1e-10, m_g,
                                   std::numeric_limits<double>::infinity());
  env.slope_r = dlog_g(env.t_r);
  if (env.slope_r >= 0.0) {
    Rcpp::stop("setup_region_d: right contact slope >= 0 (Newton failure).");
  }
  env.p_r = env.t_r + (log_g_mode - log_g(env.t_r)) / env.slope_r;

  // Build piece table.
  std::vector<EnvelopePiece> pieces;

  // 1. Leftmost exponential tail (-inf, y_1] with slope alpha (the asymptotic
  //    slope of log g(y) as y -> -inf).
  pieces.push_back(make_exp_left(y_b[1], alpha, log_g_b[1]));

  // 2. Secants on (y_{k-1}, y_k], k=2..K_eff.
  for (int k = 2; k <= K_eff; ++k) {
    pieces.push_back(make_secant(y_b[k-1], y_b[k], a_k[k], log_g_b[k-1]));
  }

  if (use_env_15) {
    // 3a. Left-tangent piece (y_K = t_hat_l, p_l] for envelope 15.
    //     h(y) = log_g(t_l) + L'(t_l) * (y - t_l) on this interval.
    //     Anchor at piece.a = t_hat_l: base = log_g(t_l) + L'(t_l)*(t_hat_l - t_l).
    const double tangent_base = log_g_t_l + slope_l_d * (y_b[K_eff] - t_l_d);
    pieces.push_back(make_exp_bounded(y_b[K_eff], p_l_d, slope_l_d, tangent_base));
    // 4a. Plateau (p_l, p_r] at height log g(m_g).
    if (env.p_r > p_l_d) {
      pieces.push_back(make_plateau(p_l_d, env.p_r, log_g_mode));
    }
    env.t_l = t_l_d;
    env.p_l = p_l_d;
    env.slope_l = slope_l_d;
    env.has_left_tangent_d = true;
  } else {
    // 3b. Plateau (y_K = y_star, p_r] for envelope 14.
    if (env.p_r > y_b[K_eff]) {
      pieces.push_back(make_plateau(y_b[K_eff], env.p_r, log_g_mode));
    }
    env.has_left_tangent_d = false;
  }

  // 5. Right tangent (p_r, inf).
  pieces.push_back(make_exp_right(env.p_r, env.slope_r, log_g_mode));

  finalize_pieces(env, pieces);
}

// Target log-density for the accept/reject test, in the coordinate system
// matching the region's envelope.  Both forms drop the (constant) log of
// the normalizing constant; the constant cancels in the ratio.
double log_target_x(double x, double alpha, double gamma_norm) {
  // Region A target on x > 0: log f(x) = (alpha-1) log x - x^2 + gamma_norm * x
  if (x <= 0.0) return -std::numeric_limits<double>::infinity();
  return (alpha - 1.0) * std::log(x) - x * x + gamma_norm * x;
}

double log_target_y(double y, double alpha, double gamma_norm) {
  // Region BC/D target on y in R: log g(y) = alpha*y - exp(2y) + gamma_norm*exp(y)
  return alpha * y - std::exp(2.0 * y) + gamma_norm * std::exp(y);
}

// ====================================================================
// Section 4: Test-only fixtures for the Newton iterator.
// ====================================================================
// Three small log-concave functions used by .rtdr_contact_point_newton_test_cpp.
// Kept tiny and pure so unit tests can verify convergence in isolation.

// 1) Standard normal log-density (up to constant): log f(t) = -t^2 / 2
inline double tf_normal_log(double t)        { return -0.5 * t * t; }
inline double tf_normal_dlog(double t)       { return -t; }

// 2) Gamma-like (alpha=3, gamma=0) on x>0:
//    log f(t) = 2 log t - t^2,  L'(t) = 2/t - 2t
inline double tf_gamma_like_log(double t)    { return 2.0 * std::log(t) - t * t; }
inline double tf_gamma_like_dlog(double t)   { return 2.0 / t - 2.0 * t; }

// 3) T_{-1/2}-concave synthetic on y in R:
//    g(y) = exp(0.7 y - exp(2y)),  log g(y) = 0.7 y - exp(2y)
//    L'(y) = 0.7 - 2 exp(2y)
inline double tf_tneghalf_log(double y)      { return 0.7 * y - std::exp(2.0 * y); }
inline double tf_tneghalf_dlog(double y)     { return 0.7 - 2.0 * std::exp(2.0 * y); }

}  // namespace

// ====================================================================
// Section 5: Public API stubs for Step 3.2 to fill.
// ====================================================================

namespace mhn {

RtdrEnvelope build_rtdr_envelope(double alpha, double beta, double gamma) {
  RtdrEnvelope env;
  env.sqrt_beta = std::sqrt(beta);
  env.alpha = alpha;
  env.gamma_norm = gamma / env.sqrt_beta;
  env.region = ::classify_region(alpha, env.gamma_norm);
  switch (env.region) {
    case REGION_A:  ::setup_region_a(env);  break;
    case REGION_BC: ::setup_region_bc(env); break;
    case REGION_D:  ::setup_region_d(env);  break;
  }
  return env;
}

double sample_rtdr(const RtdrEnvelope& env, int* retries_out) {
  // Region A: sampling and accept/reject in x-space (x > 0).
  // Region BC/D: sampling and accept/reject in y-space (y in R), then
  //   x = exp(y) on output.
  // Final scale restoration: real_X = X' / sqrt(beta).
  const bool y_space = (env.region != REGION_A);
  const int max_retries = 1000;  // safety; theory guarantees acceptance prob >= 1/e
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    const int idx = ::select_piece(env.pieces);
    const double s = ::sample_within_piece(env.pieces[idx]);
    const double log_h = ::log_piece_at(env.pieces[idx], s);
    const double log_target = y_space
        ? ::log_target_y(s, env.alpha, env.gamma_norm)
        : ::log_target_x(s, env.alpha, env.gamma_norm);
    const double log_u = std::log(R::runif(0.0, 1.0));
    if (log_u + log_h <= log_target) {
      const double x_norm = y_space ? std::exp(s) : s;
      if (retries_out != nullptr) *retries_out += retries;
      return x_norm / env.sqrt_beta;
    }
    ++retries;
  }
  Rcpp::warning("sample_rtdr: max retries exceeded; returning last proposal.");
  if (retries_out != nullptr) *retries_out += retries;
  // Fall through: return last proposal scaled.  This should be statistically
  // rare under the theoretical 1/e acceptance lower bound.
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace mhn

// ====================================================================
// Section 6: Test-only Rcpp export for isolated Newton testing.
// ====================================================================
// concavity tag selects one of three fixed test fixtures defined above.

// [[Rcpp::export(.rmhn_rtdr_cpp)]]
Rcpp::NumericVector rmhn_rtdr_cpp(int n, double alpha, double beta, double gamma) {
  if (n < 0) Rcpp::stop("n must be non-negative");
  if (n == 0) return Rcpp::NumericVector(0);
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::RtdrEnvelope env = mhn::build_rtdr_envelope(alpha, beta, gamma);
  Rcpp::NumericVector out(n);
  int retries = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = mhn::sample_rtdr(env, &retries);
  }
  out.attr("rtdr_retries") = retries;
  out.attr("rtdr_region") = static_cast<int>(env.region);
  return out;
}

// [[Rcpp::export(.dump_rtdr_envelope_cpp)]]
Rcpp::List dump_rtdr_envelope_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::RtdrEnvelope env = mhn::build_rtdr_envelope(alpha, beta, gamma);
  // Convert pieces to a list of lists for R inspection.
  const std::size_t K = env.pieces.size();
  Rcpp::List pieces(K);
  for (std::size_t i = 0; i < K; ++i) {
    const auto& p = env.pieces[i];
    pieces[i] = Rcpp::List::create(
      Rcpp::Named("type")          = static_cast<int>(p.type),
      Rcpp::Named("a")             = p.a,
      Rcpp::Named("b")             = p.b,
      Rcpp::Named("slope")         = p.slope,
      Rcpp::Named("base_log_dens") = p.base_log_dens,
      Rcpp::Named("log_area")      = p.log_area
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("region")          = static_cast<int>(env.region),
    Rcpp::Named("alpha")           = env.alpha,
    Rcpp::Named("gamma_norm")      = env.gamma_norm,
    Rcpp::Named("sqrt_beta")       = env.sqrt_beta,
    Rcpp::Named("mode")            = env.mode,
    Rcpp::Named("log_dens_mode")   = env.log_dens_mode,
    Rcpp::Named("t_l")             = env.t_l,
    Rcpp::Named("t_r")             = env.t_r,
    Rcpp::Named("p_l")             = env.p_l,
    Rcpp::Named("p_r")             = env.p_r,
    Rcpp::Named("slope_l")         = env.slope_l,
    Rcpp::Named("slope_r")         = env.slope_r,
    Rcpp::Named("simplified")          = env.simplified,
    Rcpp::Named("has_left_tangent_d")  = env.has_left_tangent_d,
    Rcpp::Named("K_eff")               = env.K_eff,
    Rcpp::Named("y_star")          = env.y_star,
    Rcpp::Named("rho")             = env.rho,
    Rcpp::Named("y_break")         = env.y_break,
    Rcpp::Named("log_dens_break")  = env.log_dens_break,
    Rcpp::Named("alpha_k")         = env.alpha_k,
    Rcpp::Named("piece_log_area")  = env.piece_log_area,
    Rcpp::Named("pieces")          = pieces
  );
}

// [[Rcpp::export(.rtdr_contact_point_newton_test_cpp)]]
Rcpp::List rtdr_contact_point_newton_test_cpp(double t_init,
                                              double log_dens_mode_val,
                                              double delta,
                                              std::string concavity,
                                              int max_iter) {
  std::function<double(double)> ld, dld;
  double low, high;
  if (concavity == "log_concave_normal") {
    ld = tf_normal_log;     dld = tf_normal_dlog;
    low = 0.46; high = 2.49;
  } else if (concavity == "log_concave_gamma_like") {
    ld = tf_gamma_like_log; dld = tf_gamma_like_dlog;
    low = 0.46; high = 2.49;
  } else if (concavity == "tneghalf") {
    ld = tf_tneghalf_log;   dld = tf_tneghalf_dlog;
    low = 0.93; high = 1.99;
  } else {
    Rcpp::stop("Unknown concavity tag: %s", concavity.c_str());
  }
  const double t = contact_point_newton(t_init, log_dens_mode_val, delta,
                                        ld, dld, low, high, max_iter);
  return Rcpp::List::create(
    Rcpp::Named("t")          = t,
    Rcpp::Named("log_f_t")    = ld(t),
    Rcpp::Named("dlog_f_t")   = dld(t),
    Rcpp::Named("increment")  = log_dens_mode_val - ld(t)
  );
}
