// Relaxed Transformed Density Rejection (RTDR) sampler for the
// Modified Half-Normal distribution, following Gao & Wang (2025).
// This file holds the common infrastructure (Newton contact-point
// search, piecewise envelope sampler) together with the
// region-specific setups (a, bc, d) and the dispatcher.

#include "mhn_rtdr.h"
#include "mhn_check.h"
#include "mhn_stable.h"

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
using mhn::PIECE_TNEGHALF_LEFT;
using mhn::PIECE_TNEGHALF_RIGHT;

// Contact-point search: the iteration of Gao & Wang (2025), Eq. (8).
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
  // Bring the start inside the bounds before the first evaluation.  The
  // recommended start of Gao & Wang (2025), Eq. (8) is m +- sqrt(-2 delta /
  // L''(m)), which is far outside for a flat mode: at alpha = 1e-5 the offset
  // is 447, and the log-axis density carries exp(2y), so log g and its
  // derivative both came back non-finite on the very first step.  The loop then
  // broke out and returned the start unchanged.
  double t = std::min(std::max(t_init, t_min), t_max);
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

// Upper bound for a right-contact search on the log axis.
//
// The log-axis density carries exp(2y), which overflows past y = 354, and the
// derivative then comes back as NaN or -Inf.  Leaving the search unbounded let
// a single Newton step land there: for alpha = 1e-4 the derivative at the start
// is 1e-4, so the step is ten thousand.  The search returned a point whose
// slope was not a number, and because NaN fails every ordinary comparison the
// sign guard let it through -- the plateau piece was then dropped without any
// error, truncating the sampler's support so that it excluded the mode.
//
// contact_point_newton bisects *toward* this bound when a step would exceed it,
// so the bound must stay above the mode even when the mode itself is large.
double right_contact_bound(double m_g) {
  return m_g + std::min(40.0, std::max(1.0, 354.0 - m_g));
}

// log1p(u) - u, without the cancellation the literal form has for small u.
// The two agree to about ten digits at u = 1e-5 and to none at 1e-9, and the
// mode-centred region A ordinate below evaluates it at u = d/m, which is 2e-10
// at gamma_norm = 1e10.  The series is the Taylor expansion of log(1+u) with
// its linear term removed, which is what makes it exact where the subtraction
// is not.
double log1p_minus_u(double u) {
  if (std::abs(u) < 1e-4) {
    return -u * u * (0.5 - u / 3.0 + u * u / 4.0);
  }
  return std::log1p(u) - u;
}

// Piece selection and per-piece inverse-CDF sampling.
// Pick a piece index proportional to exp(piece.log_area), reading the
// cumulative table that finalize_pieces built once for this envelope.
// lower_bound returns the first entry with cum >= u, which is the index the
// equivalent linear scan would have chosen.
int select_piece(const mhn::RtdrEnvelope& env) {
  const std::vector<double>& cum = env.piece_cum_area;
  const double u = R::runif(0.0, 1.0) * env.piece_area_total;
  const auto it = std::lower_bound(cum.begin(), cum.end(), u);
  if (it == cum.end()) return static_cast<int>(cum.size()) - 1;
  return static_cast<int>(it - cum.begin());
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
      const double width = piece.b - piece.a;
      return piece.a
             + mhn::log1p_u_expm1(u, piece.slope * width) / piece.slope;
    }
    case PIECE_TNEGHALF_LEFT:
      // Inverse-square hat on (-inf, pl].  With b = ppl, aux = om1,
      // F^-1: x = ppl - om1 / u  (Gao & Wang 2025, Appendix B).
      return piece.b - piece.aux / u;
    case PIECE_TNEGHALF_RIGHT:
      // Inverse-square hat on [pr, +inf).  With b = ppr, aux = om3,
      // F^-1: x = ppr + om3 / u.
      return piece.b + piece.aux / u;
    default:
      Rcpp::stop("rmhn: internal invariant violated (unknown envelope piece). Please report this to the package maintainer.");
  }
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
    case PIECE_TNEGHALF_LEFT:
    case PIECE_TNEGHALF_RIGHT:
      // h(x) = exp(base_log_dens) / (slope*(x - a) - 1)^2,  slope = L'(t)/2
      return piece.base_log_dens
             - 2.0 * std::log(std::abs(piece.slope * (x - piece.a) - 1.0));
  }
  return -std::numeric_limits<double>::infinity();
}

// Report a failure to build the envelope in the caller's own terms.
//
// These conditions are all "the construction of Gao & Wang (2025) did not
// produce a usable envelope at this parameter triple".  Naming the internal
// routine and the sign of an intermediate slope tells a user nothing they can
// act on, so say what failed, at which parameters, and what to do about it.
// The region and stage are kept as a short tag for a bug report.
[[noreturn]] void envelope_failure(const mhn::RtdrEnvelope& env,
                                   const char* stage) {
  const double beta = env.sqrt_beta * env.sqrt_beta;
  const double gamma = env.gamma_norm * env.sqrt_beta;
  Rcpp::stop("rmhn: could not construct a sampling envelope for alpha = %g, "
             "beta = %g, gamma = %g. This is a defect; please report these "
             "parameter values to the package maintainer. [%s]",
             env.alpha, beta, gamma, stage);
}

// Report the one failure that is not a defect: a tilt (or a shape) so extreme
// that the whole contact structure -- the contact points, the plateau, the
// tangency -- falls inside one ulp of the mode, so the density's own shape has
// no representable neighbourhood to sample on.  This happens at gamma_norm
// ~ 1e16 in region A, where the sampling coordinate is x and the mode is
// gamma_norm/2 against a width of 1/sqrt 2; at alpha ~ 2e31 for the same
// reason with the mode driven by the shape instead; and at gamma_norm ~ 1e15
// in region D, where the coordinate is log x and the mode's width there is
// about 1.4/gamma_norm.  Before this check the caller saw the generic
// "this is a defect, please report it" message, or -- worse, on the older code
// -- a vector that was almost all NaN.  Neither is true or useful: nothing is
// wrong with the package, the numbers simply do not exist.
[[noreturn]] void envelope_unrepresentable(const mhn::RtdrEnvelope& env) {
  const double beta = env.sqrt_beta * env.sqrt_beta;
  const double gamma = env.gamma_norm * env.sqrt_beta;
  Rcpp::stop("rmhn: cannot sample at alpha = %g, beta = %g, gamma = %g: the "
             "distribution is narrower than the spacing of the double "
             "precision numbers around its own mode, so its shape has no "
             "representable neighbourhood and no sampler can reproduce it. "
             "At this tilt the law is a normal of mean mhn_mode(alpha, beta, "
             "gamma) and variance 1/(2*beta) to far below double precision; "
             "draw from that instead.",
             env.alpha, beta, gamma);
}

// Classify the (alpha, gamma_norm) point into one of the three RTDR regions.
mhn::RtdrRegion classify_region(double alpha, double gamma_norm) {
  // Region A needs an interior mode in x, and at alpha = 1 exactly there is
  // none: the mode formula returns 0 for a non-positive tilt, log f(0) is
  // -Inf, and the envelope degenerated to a plateau of height -Inf spanning
  // [0, Inf], so every proposal was rejected and the draws came back NaN.
  // The log-axis mode used by region BC is strictly positive for every
  // alpha > 0, and the concavity test there is on the tilt alone, so alpha = 1
  // belongs on that side of the boundary.
  if (alpha > 1.0) return mhn::REGION_A;
  if (alpha >= 0.5) return mhn::REGION_BC;
  // alpha in (0, 1/2).  Threshold gamma_d = 2(1 - sqrt(1 - 2*alpha)).
  const double gamma_d = 2.0 * (1.0 - std::sqrt(1.0 - 2.0 * alpha));
  if (gamma_norm <= gamma_d) return mhn::REGION_BC;
  return mhn::REGION_D;
}

// Forward declarations of the region setups; definitions follow below.
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
  return base + mhn::log_abs_expm1(slope * width)
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
// T_{-1/2} tangent pieces (Gao & Wang 2025, Section 3.2 / Appendix B).  The hat
// is h(y) = exp(logf_t) / (half_slope*(y - t) - 1)^2 with half_slope = L'(t)/2.
// gap = L(m) - L(t) is the mode-to-contact log-drop; the piece area om equals
// exp(-gap/2) / |half_slope|.  The (unnormalised) area shared with the plateau
// carries the common factor exp(logf_m), giving log_area = (logf_m + logf_t)/2
// - log|half_slope|.  b stores the sampling anchor (ppl/ppr), aux the scale om.
EnvelopePiece make_tneghalf_left(double t, double half_slope,
                                 double logf_t, double logf_m) {
  const double om1 = std::exp(-(logf_m - logf_t) / 2.0) / half_slope;  // >0
  const double pl = t + 1.0 / half_slope - om1;
  EnvelopePiece p;
  p.type = PIECE_TNEGHALF_LEFT; p.a = t; p.b = pl + om1;  // b = ppl
  p.slope = half_slope; p.base_log_dens = logf_t; p.aux = om1;
  p.log_area = 0.5 * (logf_m + logf_t) - std::log(half_slope);
  return p;
}
EnvelopePiece make_tneghalf_right(double t, double half_slope,
                                  double logf_t, double logf_m) {
  const double om3 = -std::exp(-(logf_m - logf_t) / 2.0) / half_slope;  // half_slope<0 -> >0
  const double pr = t + 1.0 / half_slope + om3;
  EnvelopePiece p;
  p.type = PIECE_TNEGHALF_RIGHT; p.a = t; p.b = pr - om3;  // b = ppr
  p.slope = half_slope; p.base_log_dens = logf_t; p.aux = om3;
  p.log_area = 0.5 * (logf_m + logf_t) - std::log(-half_slope);
  return p;
}

void finalize_pieces(mhn::RtdrEnvelope& env, std::vector<EnvelopePiece>& pieces) {
  const std::size_t K = pieces.size();
  env.piece_log_area.clear();
  env.piece_log_area.reserve(K);
  env.piece_cum_area.clear();
  env.piece_cum_area.reserve(K);

  // Normalise by the largest log-area before exponentiating, so that a very
  // peaked envelope cannot overflow the running sum.
  double max_log_area = -std::numeric_limits<double>::infinity();
  for (const auto& p : pieces) if (p.log_area > max_log_area) max_log_area = p.log_area;

  double total = 0.0;
  for (const auto& p : pieces) {
    env.piece_log_area.push_back(p.log_area);
    total += std::exp(p.log_area - max_log_area);
    env.piece_cum_area.push_back(total);
  }
  env.piece_area_total = total;
  env.pieces = std::move(pieces);
}

// ====================================================================
// Region A setup [alpha >= 1, log-concave on f(x)]
// ====================================================================
void setup_region_a(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // mode of log f(x) = (alpha-1) log x - x^2 + gamma_norm * x, i.e. the
  // positive root of 2x^2 - gamma_norm x - (alpha-1) = 0 at beta = 1.
  const double m = mhn::positive_root(1.0, gn, alpha - 1.0);
  env.mode = m;

  auto log_f_raw = [alpha, gn](double x) -> double {
    if (x <= 0.0) return -std::numeric_limits<double>::infinity();
    return (alpha - 1.0) * std::log(x) - x * x + gn * x;
  };
  auto dlog_f_raw = [alpha, gn](double x) -> double {
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    return (alpha - 1.0) / x - 2.0 * x + gn;
  };
  const double log_f_mode_raw = log_f_raw(m);

  // The same conditioning failure region BC had, in x rather than on the log
  // axis.  At the mode -x^2 + gamma_norm x is gamma_norm^2/4 - (x -
  // gamma_norm/2)^2, so every value of log f carries gamma_norm^2/4 -- 2.5e19
  // at gamma_norm = 1e10 -- while the construction reads O(1) differences out
  // of it: the contact drop (delta = 1), the tangent intersections p_l and
  // p_r, the piece areas, and the accept test.  Measured against 300-bit
  // arithmetic at alpha = 2, gamma_norm = 1e10, a drop whose true value is -25
  // comes back as -4096 from the plain difference.
  //
  // Past the switch, measure from the mode.  With d = x - m and the mode
  // equation 2m - gamma_norm = (alpha-1)/m,
  //     log f(x) - log f(m) = (alpha-1) [ log1p(d/m) - d/m ] - d^2
  //     (log f)'(x)         = -d [ (alpha-1)/(x m) + 2 ]
  // log1p(u) - u <= 0 for every u > -1 and alpha > 1 throughout region A, so
  // both terms carry the same sign and nothing cancels.  Below the switch the
  // arithmetic is untouched and the draw sequence is bit-identical.
  const bool centred =
      std::abs(log_f_mode_raw) * std::numeric_limits<double>::epsilon() > 0.01;
  auto log_f = [&](double x) -> double {
    if (!centred) return log_f_raw(x);
    if (x <= 0.0) return -std::numeric_limits<double>::infinity();
    const double d = x - m;
    return (alpha - 1.0) * ::log1p_minus_u(d / m) - d * d;
  };
  auto dlog_f = [&](double x) -> double {
    if (!centred) return dlog_f_raw(x);
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    return -(x - m) * ((alpha - 1.0) / (x * m) + 2.0);
  };
  const double log_f_mode = centred ? 0.0 : log_f_mode_raw;
  env.centred = centred;
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

  // Contact-point starts.  The heuristic m + max(1, m/2) is an offset of order
  // m, and past the switch the ordinate behaves like -d^2 near the mode, on
  // which one Newton step only halves d.  Reaching the acceptance band
  // therefore takes log2(m/2) steps -- 26 at gamma_norm = 1e8, 33 at 1e10, 40
  // at 1e12 -- against a budget of 30, so past gamma_norm ~ 2e9 the search
  // returns a point whose drop is hundreds rather than 1.  The hat still
  // dominates (log f is concave, so every tangent does) and the draws stay
  // correctly distributed, but the plateau widens to the contact points
  // instead of the mode: measured acceptance falls from 0.86 at gamma_norm
  // = 1e9 to 0.61 at 1e10, 0.076 at 1e11 and 0.0076 at 1e12, and the
  // 1000-retry budget then starts running out -- 8 NaN draws in 2e4 at
  // gamma_norm = 1e12, and 18495 of 2e4 at 1e14.
  //
  // Past the switch use the offset recommended by Gao & Wang (2025, Eq. 8),
  // m +- sqrt(-2 delta / L''(m)) with -L''(m) = (alpha-1)/m^2 + 2, which is
  // what region BC already uses and which lands inside the band on the first
  // evaluation at every tilt measured.  Below the switch the old start always
  // converges -- !centred bounds |log f(m)| by 0.01/epsilon, and log f(m)
  // = m^2 + (alpha-1)(log m - 1) >= m^2 for m >= e, so m <= 6.7e6 and the walk
  // in costs at most 23 steps -- and it is kept there so that the draw
  // sequence stays bit-identical.
  //
  // The offset is only used when it is representable: once ulp(m) exceeds it,
  // m + inc rounds back to m, the search starts where the derivative is zero
  // and returns the mode itself, and the slope guard below would stop with a
  // message telling the user they have found a defect.  They have not -- that
  // is exactly the point at which the density is narrower than one double
  // (gamma_norm ~ 6e15, or alpha ~ 2e31 on the shape axis) and its mode's
  // neighbourhood is no longer representable in x at all.  Fall back to the
  // old start there, so that corner behaves no worse than it does today.
  const double neg_ddlog_f_mode = ((alpha - 1.0) / m) / m + 2.0;
  const double inc = std::sqrt(2.0 / neg_ddlog_f_mode);
  // Half the offset is where the plateau's edges land, so once that is below
  // one ulp of the mode the plateau collapses to a point and p_l == p_r: the
  // density is narrower than the grid it would have to be drawn on.  Say so
  // rather than letting the generic intersection guard call it a defect.
  if (centred && !(m + 0.5 * inc > m)) {
    ::envelope_unrepresentable(env);
  }
  const bool inc_usable = centred && (m + inc > m) && (m - inc < m);

  // Right contact point: t_r > m, slope_r = L'(t_r) < 0.
  const double t_r_init = inc_usable ? m + inc : m + std::max(1.0, 0.5 * m);
  env.t_r = ::contact_point_newton(t_r_init, log_f_mode, 1.0, log_f, dlog_f,
                                   0.46, 2.49, /*max_iter=*/30, /*tol=*/1e-10,
                                   /*t_min=*/m, /*t_max=*/std::numeric_limits<double>::infinity());
  env.slope_r = dlog_f(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region A right contact");
  }
  env.p_r = env.t_r + (log_f_mode - log_f(env.t_r)) / env.slope_r;
  // p_r is the plateau's right edge and the right tail's anchor.  A non-finite
  // one reaches make_plateau's log(b - a) and poisons every piece area, and so
  // every draw, with no error raised; the simplified branch below skips the
  // p_l guard entirely and would otherwise check no endpoint at all.
  if (!std::isfinite(env.p_r)) {
    ::envelope_failure(env, "region A right intersection");
  }

  std::vector<EnvelopePiece> pieces;

  if (!simplified) {
    // Left contact point: 0 < t_l < m, slope_l = L'(t_l) > 0.
    // Start at half the mode.  An absolute floor of 1e-6 was applied here as
    // well, which put the start to the *right* of the mode whenever the mode
    // was smaller than 2e-6 -- and the search tests the acceptance band before
    // it clamps into [t_min, t_max], so it could return a point on the wrong
    // side with a non-positive slope.  rmhn(200, 4, 1, -1e7), where the mode is
    // 3e-7, failed outright on the default path.  m > 0 here because this
    // branch is only reached for alpha > 1.
    // m - min(inc, m/2) is the Eq. (8) start kept strictly inside (0, m): a
    // bare m - inc lands at or below zero whenever the mode is smaller than the
    // mode's own width, and the search tests the acceptance band before it
    // clamps into [t_min, t_max].
    const double t_l_init =
        inc_usable ? m - std::min(inc, 0.5 * m) : 0.5 * m;
    env.t_l = ::contact_point_newton(t_l_init, log_f_mode, 1.0, log_f, dlog_f,
                                     0.46, 2.49, /*max_iter=*/30, /*tol=*/1e-10,
                                     /*t_min=*/0.0, /*t_max=*/m);
    if (env.t_l <= 0.0) env.t_l = 0.5 * m;  // defensive
    env.slope_l = dlog_f(env.t_l);
    if (!(env.slope_l > 0.0) || !std::isfinite(env.slope_l)) {
      ::envelope_failure(env, "region A left contact");
    }
    env.p_l = env.t_l + (log_f_mode - log_f(env.t_l)) / env.slope_l;
    // Negated rather than written p_l <= 0 || p_l >= p_r so a non-finite
    // endpoint is caught: every ordinary comparison is false for NaN, and a
    // NaN p_l would otherwise walk past this guard and poison every piece area.
    if (!(env.p_l > 0.0) || !(env.p_l < env.p_r)) {
      ::envelope_failure(env, "region A intersection");
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
// Region BC setup [alpha <= 1, T_{-1/2}-concave on g(y) = exp(alpha*y
//   - exp(2y) + gamma_norm * exp(y)), y in R]
// ====================================================================
void setup_region_bc(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // Mode in u = exp(y): 2 u^2 - gamma_norm u - alpha = 0,
  //   u_+ = (gamma_norm + sqrt(gamma_norm^2 + 8 alpha)) / 4 > 0.
  const double u_mode = mhn::positive_root(1.0, gn, alpha);
  const double m_g = std::log(u_mode);
  env.mode = m_g;

  auto log_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha * y - ey * ey + gn * ey;
  };
  auto dlog_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha - 2.0 * ey * ey + gn * ey;
  };
  const double log_g_mode_raw = log_g_raw(m_g);

  // Near the mode, -e^{2y} + gamma_norm e^y = gamma_norm^2/4 - (e^y -
  // gamma_norm/2)^2, so every evaluation of log g carries the constant
  // gamma_norm^2/4 -- 2.5e23 at gamma_norm = 1e12.  Everything the envelope
  // reads out of log g is an O(1) *difference* of two such values: the contact
  // drop (delta = log 4), the tangent gaps in om1/om3, the between-piece area
  // ratios, and the accept/reject test.  Once ulp(log g(m_g)) approaches delta
  // those differences carry no digits -- at gamma_norm = 1e12 they are
  // quantised to multiples of 3.4e7 against a target of log 4 -- so t_l, om1
  // and p_l are decided by the last bit of exp().  That is the whole of the
  // platform split: the same triple built a NaN-bounded envelope on one libm
  // and tripped the p_l >= p_r guard on another, non-monotonically in gamma.
  //
  // Past that point, measure the ordinate from the mode.  With d = y - m_g,
  // E = expm1(d), and the mode equation 2 u_m^2 - gamma_norm u_m - alpha = 0,
  //     log g(y) - log g(m_g) = -( u_m^2 E^2 + alpha (E - d) )
  //     (log g)'(y)           = -E ( alpha + 2 u_m^2 e^d )
  // Both only ever add same-signed terms, so the drop stays exact at any tilt.
  // The switch is placed at ulp(log g(m_g)) = delta/100, i.e. gamma_norm
  // ~ 1.6e7, which measurement puts just below the tilt at which the plain
  // difference starts to bend the sampled distribution (a KS test of 2e5 draws
  // against the exact Gaussian limit is clean at gamma_norm = 2.5e7,
  // D = 0.0025, and fails at 4.0e7, D = 0.0054) and far above any ordinary
  // tilt -- the package's own test grid reaches gamma_norm = 1e6.  Below the
  // switch the arithmetic is untouched and the draw sequence is bit-identical.
  const bool centred =
      std::abs(log_g_mode_raw) * std::numeric_limits<double>::epsilon() >
      0.01 * ((gn > 0.0) ? std::log(4.0) : 1.0);
  const double um2 = u_mode * u_mode;
  auto log_g = [&](double y) -> double {
    if (!centred) return log_g_raw(y);
    const double d = y - m_g;
    const double E = std::expm1(d);
    return -(um2 * E * E + alpha * (E - d));
  };
  auto dlog_g = [&](double y) -> double {
    if (!centred) return dlog_g_raw(y);
    const double d = y - m_g;
    return -std::expm1(d) * (alpha + 2.0 * um2 * std::exp(d));
  };
  const double log_g_mode = centred ? 0.0 : log_g_mode_raw;
  env.centred = centred;
  env.u_mode = u_mode;
  env.log_dens_mode = log_g_mode;
  env.simplified = false;  // BC envelope is always 3-piece

  // g(y) is log-concave iff gamma_norm <= 0 (Gao & Wang 2025, Theorem 3.2).
  // For gamma_norm <= 0 use the T_0 (log-tangent) hat with acceptance band
  // [0.46, 2.49] and delta = 1.  For gamma_norm > 0 the density is only
  // T_{-1/2}-concave -- log g is convex on (-inf, log(gamma_norm/4)) -- so a
  // log-tangent hat would fail to dominate there; use the T_{-1/2}
  // (inverse-square) hat with band [0.93, 1.99] and delta = log 4.
  const bool tneghalf = (gn > 0.0);
  const double delta = tneghalf ? std::log(4.0) : 1.0;
  const double band_lo = tneghalf ? 0.93 : 0.46;
  const double band_hi = tneghalf ? 1.99 : 2.49;

  // Recommended contact-point start (Gao & Wang 2025, Eq. 8):
  // inc = sqrt(-2*delta / L''(m_g)), with L''(m_g) = -u_mode*(4 u_mode - gn) < 0.
  const double ddlog_g_mode = -u_mode * (4.0 * u_mode - gn);
  const double inc = std::sqrt(-2.0 * delta / ddlog_g_mode);

  // Left contact point: t_l < m_g, slope_l > 0.
  env.t_l = ::contact_point_newton(m_g - inc, log_g_mode, delta,
                                   log_g, dlog_g, band_lo, band_hi,
                                   /*max_iter=*/30, /*tol=*/1e-10,
                                   /*t_min=*/-std::numeric_limits<double>::infinity(),
                                   /*t_max=*/m_g);
  env.slope_l = dlog_g(env.t_l);
  if (!(env.slope_l > 0.0) || !std::isfinite(env.slope_l)) {
    ::envelope_failure(env, "region BC left contact");
  }
  const double log_g_tl = log_g(env.t_l);

  // Right contact point: t_r > m_g, slope_r < 0.
  env.t_r = ::contact_point_newton(m_g + inc, log_g_mode, delta,
                                   log_g, dlog_g, band_lo, band_hi,
                                   /*max_iter=*/200, /*tol=*/1e-10,
                                   /*t_min=*/m_g,
                                   /*t_max=*/::right_contact_bound(m_g));
  env.slope_r = dlog_g(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region BC right contact");
  }
  const double log_g_tr = log_g(env.t_r);

  std::vector<EnvelopePiece> pieces;
  if (tneghalf) {
    // T_{-1/2} inverse-square tangents; the tangent/plateau intersection
    // p_l = ppl - om1 (and p_r = ppr + om3) come from the built pieces.
    const EnvelopePiece left =
        make_tneghalf_left(env.t_l, env.slope_l / 2.0, log_g_tl, log_g_mode);
    const EnvelopePiece right =
        make_tneghalf_right(env.t_r, env.slope_r / 2.0, log_g_tr, log_g_mode);
    env.p_l = left.b - left.aux;    // ppl - om1
    env.p_r = right.b + right.aux;  // ppr + om3
    // Written !(p_l < p_r) rather than p_l >= p_r so a non-finite endpoint is
    // caught: `>=` is false for NaN, and a NaN p_l used to walk past this
    // guard into the plateau's width and poison every piece area.
    if (!(env.p_l < env.p_r)) {
      ::envelope_failure(env, "region BC intersection");
    }
    pieces.push_back(left);
    pieces.push_back(make_plateau(env.p_l, env.p_r, log_g_mode));
    pieces.push_back(right);
  } else {
    // T_0 log-tangents (g log-concave for gamma_norm <= 0).
    env.p_l = env.t_l + (log_g_mode - log_g_tl) / env.slope_l;
    env.p_r = env.t_r + (log_g_mode - log_g_tr) / env.slope_r;
    if (!(env.p_l < env.p_r)) {
      ::envelope_failure(env, "region BC intersection");
    }
    pieces.push_back(make_exp_left(env.p_l, env.slope_l, log_g_mode));
    pieces.push_back(make_plateau(env.p_l, env.p_r, log_g_mode));
    pieces.push_back(make_exp_right(env.p_r, env.slope_r, log_g_mode));
  }
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

  if (gn <= 0.0) ::envelope_failure(env, "region D tilt");

  // Mode (in y) and inflection point (in y).  y_star is where L'' changes
  // sign: L''(y) = e^y (gamma_norm - 4 e^y), so log g is convex to the left of
  // log(gamma_norm/4) and concave to the right -- the secant chain covers the
  // convex half and the tangent hat the concave half.
  const double disc = std::sqrt(gn * gn + 8.0 * alpha);   // S
  const double u_mode = mhn::positive_root(1.0, gn, alpha);
  const double m_g = std::log(u_mode);
  const double y_star = std::log(gn / 4.0);
  env.mode = m_g;
  env.y_star = y_star;

  auto log_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha * y - ey * ey + gn * ey;
  };
  auto dlog_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha - 2.0 * ey * ey + gn * ey;
  };
  const double log_g_mode_raw = log_g_raw(m_g);

  // Region D works on the same g(y) as region BC and carries the same defect
  // verbatim.  Near the mode -e^{2y} + gamma_norm e^y = gamma_norm^2/4 -
  // (e^y - gamma_norm/2)^2, so every value of log g carries gamma_norm^2/4 --
  // 2.5e19 at gamma_norm = 1e10, one ulp of which is 4096 -- while everything
  // the envelope reads out of it is an O(1) *difference* of two such values:
  // the contact drop (delta = log 4), the tangent/plateau intersections p_l
  // and p_r, the secant slopes, the between-piece area ratios exponentiated in
  // finalize_pieces, and the accept/reject test.  Measured: at gamma_norm
  // = 1e9 and above the numerator of p_l rounds to exactly zero, so p_l == t_l
  // bit for bit; at 1e10 p_r came out to the *left* of p_l, the plateau was
  // silently dropped, and the sampler's support no longer contained the mode
  // (sd/true = 0.0146, KS D = 1.000 against the exact Gaussian limit).
  //
  // The remedy is region BC's: with d = y - m_g, E = expm1(d), and the mode
  // equation 2 u_m^2 - gamma_norm u_m - alpha = 0,
  //     log g(y) - log g(m_g) = -( u_m^2 E^2 + alpha (E - d) )
  //     (log g)'(y)           = -E ( alpha + 2 u_m^2 e^d )
  // Both only ever add same-signed terms (E - d = expm1(d) - d >= 0), so the
  // drop stays exact at any tilt.  Region D always has gamma_norm > 0 and so
  // always uses the T_{-1/2} hat, delta = log 4.
  //
  // The switch is region BC's rule -- ulp(log g(m_g)) against a fraction of
  // delta -- but with a tenth of BC's fraction, which is a measured decision
  // rather than a stylistic one.  Six seeds of 2e5 draws each, KS against the
  // exact Gaussian limit on the raw path, give a combined p of 0.73 at
  // gamma_norm = 1e6, 0.56 at 5e6, 0.27 at 1.1e7, 0.41 at 1.6e7, then 0.0048
  // at 2e7 and 7e-21 at 3e7.  BC's own constant would put the switch at
  // gamma_norm 1.6e7, a factor of 1.25 below that onset; a tenth of it puts
  // the switch at 5.0e6, a factor of four in tilt and sixteen in ulp below it,
  // and still five times above the largest tilt the package's own test grid
  // reaches.  Region D earns the wider margin because it has strictly more
  // difference sites than BC's three-piece table -- the dual point, the secant
  // chain and the left-tangent/plateau join on top of the same contact search
  // and the same area ratios -- and because the cost of switching early is
  // nil: at gamma_norm = 5e6 the two constructions already agree to a part in
  // 1e4.  Below the switch every expression reverts to the exact pre-existing
  // one and the draw sequence is bit-identical.
  const double delta = std::log(4.0);
  const bool centred =
      std::abs(log_g_mode_raw) * std::numeric_limits<double>::epsilon() >
      0.001 * delta;
  const double um2 = u_mode * u_mode;
  auto log_g = [&](double y) -> double {
    if (!centred) return log_g_raw(y);
    const double d = y - m_g;
    const double E = std::expm1(d);
    return -(um2 * E * E + alpha * (E - d));
  };
  auto dlog_g = [&](double y) -> double {
    if (!centred) return dlog_g_raw(y);
    const double d = y - m_g;
    return -std::expm1(d) * (alpha + 2.0 * um2 * std::exp(d));
  };
  const double log_g_mode = centred ? 0.0 : log_g_mode_raw;
  const double log_g_y_star = log_g(y_star);
  env.centred = centred;
  env.u_mode = u_mode;
  env.log_dens_mode = log_g_mode;
  env.simplified = false;

  // Contact-point starts on the curvature scale (Gao & Wang 2025, Eq. 8), as
  // region BC uses, with L''(m_g) = -u_mode (4 u_mode - gamma_norm) < 0.  The
  // old starts -- the midpoint of [y_star, m_g] on the left and m_g + 1 on the
  // right -- sit a fixed O(1) distance from the mode on the log axis while the
  // contact points sit 2 sqrt(log 4) / gamma_norm away, 2.4e-10 at gamma_norm
  // = 1e10.  Newton halves that distance per step, so it needs log2(0.15
  // gamma_norm) steps on the left: 30.5 at 1e10, 37 at 1e12, against a budget
  // of 30.  Even in exact arithmetic the left search then returned a point
  // whose drop was 2.28 at 1e10 and 12731 at 1e12, outside the [0.93, 1.99]
  // band it is supposed to land in.  The right search has a budget of 200 but
  // is stopped instead by the *absolute* step tolerance, which exceeds the
  // whole contact structure once gamma_norm passes ~1.4e10; scale the
  // tolerance to the mode's own width so it means the same thing at every
  // tilt.  From the curvature start the first evaluation is already inside the
  // band, so neither limit is reached: measured acceptance is 0.87 flat from
  // gamma_norm = 1e6 to 1e14, against 0.24 at 1e10 and 0.27 at 1e12 before.
  const double ddlog_g_mode = -u_mode * (4.0 * u_mode - gn);
  const double inc = std::sqrt(-2.0 * delta / ddlog_g_mode);
  const double newton_tol = centred ? 1e-6 * inc : 1e-10;
  // The plateau's edges sit half an offset from the mode on the log axis.
  // Once that is inside one ulp of m_g -- gamma_norm ~ 1e15, where the mode is
  // near 34 and its width is 2e-15 against an ulp of 7e-15 -- the whole
  // contact structure collapses onto the mode and the search returns the mode
  // itself, with a zero slope and a non-positive dual point.  That is the
  // representational floor of the log axis, not a defect in the construction.
  if (centred && !(m_g - 0.5 * inc < m_g)) {
    ::envelope_unrepresentable(env);
  }

  // Envelope selection.  log(g(m_g)/g(y_star)) > 2.49 triggers env 15.
  // Centred, log_g_y_star is already the negated drop and log_g_mode is zero,
  // so the test reads the same in both coordinate systems.  It is well
  // conditioned either way: the drop is gamma_norm^2/16, not O(1).
  const bool use_env_15 = (log_g_mode - log_g_y_star) > 2.49;

  // Compute rho according to the selected envelope.  env 15 also derives
  // the contact point t_l (in the log-concave half [y_star, m_g]) and the
  // tangent intersection p_l for the extra left-tangent piece.
  double rho;
  double t_l_d = 0.0, t_hat_l = 0.0, slope_l_d = 0.0;
  double log_g_t_l = 0.0, p_l_d = 0.0;
  // rho is L_g(y_hat) - alpha y_hat.  Since L_g(y) - alpha y is exactly
  // -e^{2y} + gamma_norm e^y, that is
  //     rho = v (gamma_norm - v),      v = e^{y_hat} <= gamma_norm / 2,
  // a product of two positive factors.  Read instead as a difference of two
  // log-density values it subtracts two numbers of size gamma_norm^2/4 to
  // leave one of size gamma_norm: measured against the theoretical value it
  // came out 41x too large at gamma_norm = 1e10 and 6e5 times too large at
  // 1e14.  It also cannot be read off a *centred* log g at all -- that would
  // give rho ~ -gamma_norm^2/4, tripping the rho <= 0 guard below -- so the
  // closed form is mandatory once the switch fires rather than merely better.
  // For envelope 14, v = e^{y_star} = gamma_norm/4 and it reduces to the
  // documented 3 gamma_norm^2 / 16.
  if (!use_env_15) {
    // Envelope 14: y_hat_star = y_star.  Closed form rho = 3 gn^2 / 16.
    rho = centred ? (gn / 4.0) * (gn - gn / 4.0)
                  : (log_g_y_star - alpha * y_star);
  } else {
    // Envelope 15: contact-point search inside [y_star, m_g] (log-concave
    // sub-region of g).  contact_point_newton clamps the start into
    // [y_star, m_g] itself, so the flat-mode case the mid-point start was
    // guarding against is still covered by the curvature start.
    const double y_0 = std::log(gn / 2.0);   // dual point of -inf
    const double t_l_init = centred ? m_g - inc : 0.5 * (y_star + m_g);
    t_l_d = ::contact_point_newton(t_l_init, log_g_mode, delta,
                                   log_g, dlog_g, 0.93, 1.99,
                                   30, newton_tol, y_star, m_g);
    // The dual point.  v = e^{t_hat_l} = gamma_norm/2 - e^{t_l} equals
    // sqrt(log 4) = 1.1774 for *every* gamma_norm, so written literally it
    // subtracts two numbers of size gamma_norm/2 to leave an O(1) answer and
    // loses a digit per decade of tilt.  With d = t_l - m_g,
    // e^{t_l} = u_mode (1 + expm1(d)), and the conjugate form of the mode
    // equation gives gamma_norm/2 - u_mode = (gamma_norm - S)/4 =
    // -2 alpha / (gamma_norm + S) exactly, so
    //     v = -u_mode expm1(d) - 2 alpha / (gamma_norm + S)
    // adds two terms of known sign (expm1(d) < 0 here, and the subtrahend is
    // the tiny one).  The old guard t_l < log(gamma_norm/2) is exactly v > 0,
    // and written this way it also traps a non-finite t_l.
    const double v = centred
        ? (-u_mode * std::expm1(t_l_d - m_g) - 2.0 * alpha / (gn + disc))
        : (gn / 2.0 - std::exp(t_l_d));
    if (centred ? !(v > 0.0) : !(t_l_d < y_0)) {
      ::envelope_failure(env, "region D envelope 15 breakpoint");
    }
    t_hat_l = std::log(v);
    slope_l_d = dlog_g(t_l_d);
    log_g_t_l = log_g(t_l_d);
    if (slope_l_d <= alpha) {
      // Gao & Wang (2025) Lemma 4.5 requires L'(t_l) > alpha for the
      // dual-point construction.
      ::envelope_failure(env, "region D envelope 15 slope");
    }
    p_l_d = t_l_d + (log_g_mode - log_g_t_l) / slope_l_d;
    rho = centred ? v * (gn - v)
                  : (log_g_raw(t_hat_l) - alpha * t_hat_l);
  }
  env.rho = rho;
  if (rho <= 0.0) ::envelope_failure(env, "region D secant count");

  // K = ceil(rho) is the theoretical secant count.  For envelope 15 rho is
  // about sqrt(log 4) * gamma_norm, so ceil(rho) passes INT_MAX at gamma_norm
  // ~ 1.8e9 -- well inside the default dispatch, since rmhn's auto path sends
  // every alpha < 1 with gamma > 0 to RTDR.  Converting a double above INT_MAX
  // to int is undefined behaviour and the architectures disagree about it:
  // arm64 saturates to INT_MAX and stays in region D, x86-64 yields INT_MIN,
  // which fails the K_eff >= 1 test below and hands the triple silently to the
  // region BC envelope -- a T_{-1/2} tangent hat for a density that, by the
  // definition of this region, is not T_{-1/2}-concave anywhere (the minimum
  // of (L')^2/2 - L'' over the axis is exactly (2 alpha - 1)/2 < 0 for every
  // alpha < 1/2, at every tilt).  So the two platforms ran different
  // algorithms, one of them without a domination guarantee.  Nothing needs K
  // as an int; ceil is exact in double, so keeping the count there is
  // value-identical wherever the old conversion was defined.
  const double K = std::ceil(rho);

  // K_eff guard: largest k <= K with gn^2 > 4 k (rho/K).  Theoretically
  // both envelopes always satisfy the guard at every k -- gn^2 - 4 rho is
  // (gn - 2v)^2 for envelope 15 and gn^2/4 for envelope 14, both positive, so
  // the condition is already false at k = K and the loop exits at once -- but
  // we keep the loop for defensive programming.
  double K_eff = K;
  while (K_eff >= 1.0 && gn * gn <= 4.0 * K_eff * (rho / K)) {
    K_eff -= 1.0;
  }

  if (K_eff < 1.0) {
    // No secant breakpoint fits; the region-BC envelope is equally valid for
    // this density, so use it and record the fallback on the envelope.
    env.fell_back_to_bc = true;
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
  const double drho = rho / K;
  std::vector<double> y_b;
  std::vector<double> log_g_b;
  y_b.reserve(64);
  log_g_b.reserve(64);
  y_b.push_back(0.0);          // 1-indexed: y_b[0] unused
  log_g_b.push_back(0.0);
  // K_eff is a double because K can exceed INT_MAX, so the counter is a long
  // long and is capped.  Envelope 14 is only selected when the mode-to-
  // inflection drop is at most 2.49, which forces gamma_norm <= 6.31 and hence
  // ceil(rho) = ceil(3 gamma_norm^2/16) <= 8; envelope 15 breaks out by Lemma
  // 4.1 within a handful of steps at every parameter measured.  If neither
  // bound fires, truncating the chain would leave the left tangent spanning an
  // interval Lemma 4.1 does not cover, i.e. a hat that is not one, so fail
  // loudly rather than build it.
  const long long k_cap = 1000000LL;
  long long K_actual = 0;
  for (long long k = 1; static_cast<double>(k) <= K_eff && k <= k_cap; ++k) {
    const double rk = 4.0 * static_cast<double>(k) * drho;
    const double inside = gn * gn - rk;
    if (inside <= 0.0) {
      ::envelope_failure(env, "region D breakpoint");
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
  if (K_actual == k_cap && static_cast<double>(k_cap) < K_eff) {
    ::envelope_failure(env, "region D secant chain");
  }
  const int K_used = static_cast<int>(K_actual);
  env.K_eff = K_used;
  env.y_break.assign(y_b.begin() + 1, y_b.end());
  env.log_dens_break.assign(log_g_b.begin() + 1, log_g_b.end());

  // Secant slopes alpha_k for k=2..K_used.
  //
  // y_k above is built as the *small* root of u^2 - gamma_norm u + k drho = 0,
  // so at every breakpoint -u_k^2 + gamma_norm u_k = k drho exactly and
  // therefore
  //     L(y_k) - L(y_{k-1}) = alpha (y_k - y_{k-1}) + drho
  //     alpha_k             = alpha + drho / (y_k - y_{k-1}).
  // Written as a difference quotient of two ordinates this is the
  // worst-conditioned expression in the region: the true numerator is O(1)
  // while each ordinate is O(rho) raw and O(gamma_norm^2/4) once the ordinates
  // are measured from the mode.  Centred, the difference rounds to zero, and a
  // zero slope makes log_area_exp_bounded return base + log|expm1(0)| -
  // log(0) = NaN, which poisons piece_cum_area and hands std::lower_bound a
  // NaN-ordered array.  The identity never forms the difference at all, so it
  // is immune to the coordinate choice; it is gated only so that the draw
  // sequence below the switch stays bit-identical.
  std::vector<double> a_k(K_used + 1);
  for (int k = 2; k <= K_used; ++k) {
    a_k[k] = centred
        ? (alpha + drho / (y_b[k] - y_b[k-1]))
        : ((log_g_b[k] - log_g_b[k-1]) / (y_b[k] - y_b[k-1]));
  }
  if (K_used >= 2) {
    env.alpha_k.assign(a_k.begin() + 2, a_k.end());
  } else {
    env.alpha_k.clear();
  }

  // Right contact point in the log-concave region (same as BC right tangent).
  // The curvature start and the scaled tolerance matter here for the same
  // reason they do on the left: from m_g + 1 the search needs 35 halvings at
  // gamma_norm = 1e10 and the absolute 1e-10 step tolerance ends it earlier
  // still, returning a contact whose drop was 16.2 at gamma_norm = 1e11 and
  // 1526 at 1e12 against the target log 4.  That widens the plateau by a
  // factor of 33 at 1e12; the draws stay correctly distributed, because a
  // wider hat over a correct accept test is still a correct sampler, so the
  // KS invariant does not see it -- only the acceptance rate does, and it fell
  // to 0.086 against the 1/e floor Gao & Wang (2025) prove.
  const double t_r_init = centred ? m_g + inc : m_g + 1.0;
  env.t_r = ::contact_point_newton(t_r_init, log_g_mode, delta,
                                   log_g, dlog_g, 0.93, 1.99,
                                   200, newton_tol, m_g,
                                   ::right_contact_bound(m_g));
  env.slope_r = dlog_g(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region D right contact");
  }
  env.p_r = env.t_r + (log_g_mode - log_g(env.t_r)) / env.slope_r;

  // Build piece table.
  std::vector<EnvelopePiece> pieces;

  // 1. Leftmost exponential tail (-inf, y_1] with slope alpha (the asymptotic
  //    slope of log g(y) as y -> -inf).
  pieces.push_back(make_exp_left(y_b[1], alpha, log_g_b[1]));

  // 2. Secants on (y_{k-1}, y_k], k=2..K_used.
  for (int k = 2; k <= K_used; ++k) {
    pieces.push_back(make_secant(y_b[k-1], y_b[k], a_k[k], log_g_b[k-1]));
  }

  if (use_env_15) {
    // 3a. Left-tangent piece (y_K = t_hat_l, p_l] for envelope 15.
    //     h(y) = log_g(t_l) + L'(t_l) * (y - t_l) on this interval.
    //     Anchor at piece.a = t_hat_l: base = log_g(t_l) + L'(t_l)*(t_hat_l - t_l).
    // The plateau's left edge must lie to the right of the last breakpoint, or
    // the piece is empty and its area formula below is a log of a
    // non-positive number.
    if (!(y_b[K_used] < p_l_d)) {
      ::envelope_failure(env, "region D tangent interval");
    }
    const double tangent_base = log_g_t_l + slope_l_d * (y_b[K_used] - t_l_d);
    EnvelopePiece tangent =
        make_exp_bounded(y_b[K_used], p_l_d, slope_l_d, tangent_base);
    if (centred) {
      // By construction of p_l this tangent reaches the mode height there, so
      // its area is exp(log g(m_g)) (1 - exp(-slope * width)) / slope.  The
      // generic helper instead adds log|expm1(slope*width)| to a base of the
      // same magnitude and opposite sign -- both about 5.3e11 at gamma_norm
      // = 1e10 -- and loses ulp(slope*width) out of a piece carrying a fifth
      // of the envelope area.  This form has no cancellation.
      tangent.log_area =
          log_g_mode - std::log(slope_l_d)
          + std::log1p(-std::exp(-slope_l_d * (p_l_d - y_b[K_used])));
    }
    pieces.push_back(tangent);
    // 4a. Plateau (p_l, p_r] at height log g(m_g).  p_l < m_g < p_r holds
    // mathematically, so an empty plateau means the construction failed;
    // dropping it silently leaves a support that does not contain the mode,
    // which is exactly what produced sd/true = 0.0146 at gamma_norm = 1e10.
    // Written !(p_l < p_r) so a non-finite endpoint is caught too, as region
    // BC's guard already is.
    if (!(p_l_d < env.p_r)) {
      ::envelope_failure(env, "region D intersection");
    }
    pieces.push_back(make_plateau(p_l_d, env.p_r, log_g_mode));
    env.t_l = t_l_d;
    env.p_l = p_l_d;
    env.slope_l = slope_l_d;
    env.has_left_tangent_d = true;
  } else {
    // 3b. Plateau (y_K = y_star, p_r] for envelope 14.
    if (!(y_b[K_used] < env.p_r)) {
      ::envelope_failure(env, "region D intersection");
    }
    pieces.push_back(make_plateau(y_b[K_used], env.p_r, log_g_mode));
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

// The same target measured from the mode, for an envelope whose ordinate is
// centred (see setup_region_bc).  Written this way the accept/reject test
// compares two O(1) numbers; in the raw form both sides are gamma_norm^2/4, so
// at a large tilt log_u + log_h == log_h exactly and the uniform draw stops
// influencing acceptance at all.
// Region A's target measured from the mode, for the same reason.
double log_target_x_centred(double x, double alpha, double m) {
  if (x <= 0.0) return -std::numeric_limits<double>::infinity();
  const double d = x - m;
  return (alpha - 1.0) * ::log1p_minus_u(d / m) - d * d;
}

double log_target_y_centred(double y, double alpha, double m_g, double u_mode) {
  const double d = y - m_g;
  const double E = std::expm1(d);
  return -(u_mode * u_mode * E * E + alpha * (E - d));
}

// ====================================================================
// Fixtures used only by the contact-point unit test below.
// ====================================================================
// Three small targets whose contact points can be written down, so that
// tests/testthat/test-rtdr-newton.R can check the search in isolation from any
// envelope: a normal log-density, a gamma-like kernel on t > 0, and a
// T_{-1/2}-concave function on the log axis.

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
// Public API.
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
    const int idx = ::select_piece(env);
    const double s = ::sample_within_piece(env.pieces[idx]);
    const double log_h = ::log_piece_at(env.pieces[idx], s);
    // The target must be measured in whatever coordinate the envelope was
    // built in.  If they disagree the two sides of the test are again of size
    // gamma_norm^2/4, log_u + log_h == log_h exactly, and the uniform draw
    // stops influencing acceptance at all -- the sampler then draws from its
    // own hat.
    const double log_target =
        env.centred
          ? (y_space ? ::log_target_y_centred(s, env.alpha, env.mode, env.u_mode)
                     : ::log_target_x_centred(s, env.alpha, env.mode))
          : (y_space ? ::log_target_y(s, env.alpha, env.gamma_norm)
                     : ::log_target_x(s, env.alpha, env.gamma_norm));
    const double log_u = std::log(R::runif(0.0, 1.0));
    if (log_u + log_h <= log_target) {
      const double x_norm = y_space ? std::exp(s) : s;
      if (retries_out != nullptr) *retries_out += retries;
      return x_norm / env.sqrt_beta;
    }
    ++retries;
  }
  // Retry budget exhausted, which the 1/e acceptance lower bound of
  // Gao & Wang (2025) makes vanishingly unlikely.  Report it as NaN and leave
  // the reporting to the caller: raising an R warning here would leave this
  // frame by a long jump under options(warn = 2), skipping the destructors of
  // the envelope's vectors.
  if (retries_out != nullptr) *retries_out += retries;
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace mhn

// ====================================================================
// Entry points used by the test suite; not part of the public API.
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

// [[Rcpp::export(.dump_rtdr_envelope_cpp, rng = false)]]
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
      Rcpp::Named("log_area")      = p.log_area,
      Rcpp::Named("aux")           = p.aux
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("region")          = static_cast<int>(env.region),
    Rcpp::Named("alpha")           = env.alpha,
    Rcpp::Named("gamma_norm")      = env.gamma_norm,
    Rcpp::Named("sqrt_beta")       = env.sqrt_beta,
    Rcpp::Named("mode")            = env.mode,
    Rcpp::Named("log_dens_mode")   = env.log_dens_mode,
    // Exposed so that a test rebuilding log f or log g for a domination check
    // knows which ordinate the piece table is measured in: past the switch
    // every base_log_dens and log_area is relative to the peak, and comparing
    // them against a raw log-density silently compares the wrong two things.
    Rcpp::Named("centred")         = env.centred,
    Rcpp::Named("u_mode")          = env.u_mode,
    Rcpp::Named("t_l")             = env.t_l,
    Rcpp::Named("t_r")             = env.t_r,
    Rcpp::Named("p_l")             = env.p_l,
    Rcpp::Named("p_r")             = env.p_r,
    Rcpp::Named("slope_l")         = env.slope_l,
    Rcpp::Named("slope_r")         = env.slope_r,
    Rcpp::Named("simplified")          = env.simplified,
    Rcpp::Named("has_left_tangent_d")  = env.has_left_tangent_d,
    Rcpp::Named("fell_back_to_bc")     = env.fell_back_to_bc,
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

// [[Rcpp::export(.rtdr_contact_point_newton_test_cpp, rng = false)]]
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
