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
using mhn::PIECE_EXP_BOUNDED_RIGHT;

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
// Defaults are +/-inf.
//
// The band is not guaranteed.  If it is not reached within max_iter, or the
// derivative or the step degenerates, the last iterate is returned, and no
// caller tests its drop against the band; the callers check only the sign and
// finiteness of the slope there and the order of the tangent intersections.
// A tangent to a concave transformed density lies above it whatever the
// contact point, so the hat stays valid and only the bound on the rejection
// constant that the band buys (Gao & Wang 2025, Eq. (7)) is lost.  The case
// met in practice is region A's left contact for alpha just above 1 and a
// small tilt: log f falls towards x = 0 only like (alpha - 1) log(m / x),
// each step is cut back to a halving of t, and the band may lie below the
// smallest positive double as well.  The rejection constant measured there is
// about 1.1.
//
// The callables are template parameters so that the region setups' lambdas are
// called directly; a std::function would allocate for each of them.
template <class LogDens, class DLogDens>
double contact_point_newton(double t_init,
                            double log_dens_mode_val,
                            double delta,
                            const LogDens& log_dens_at,
                            const DLogDens& log_dens_deriv_at,
                            double accept_low,
                            double accept_high,
                            int max_iter = 30,
                            double tol = 1e-10,
                            double t_min = -std::numeric_limits<double>::infinity(),
                            double t_max =  std::numeric_limits<double>::infinity()) {
  // Bring the start inside the bounds before the first evaluation.  The start
  // Gao & Wang (2025, Section 2) recommend for their Eq. (8) iteration is
  // m +- sqrt(-2 delta / L''(m)), which is far outside for a flat mode: at
  // alpha = 1e-5 the offset is 447 under a negative tilt (316 at gamma = 0,
  // about 372 under a small positive one), and the log-axis density carries
  // exp(2y), so without the clamp log g and its derivative would both be
  // non-finite on the very first step, and the loop would break out and
  // return the start unchanged.
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
// derivative then comes back as NaN or -Inf.  Unbounded, a single Newton step
// can land there: for alpha = 1e-4 the derivative at the start is 1e-4, so the
// step is ten thousand.  The search would then return a point whose slope is
// not a number, which fails every ordinary comparison, so a sign test written
// as slope >= 0 would let it through, the plateau piece would be dropped
// without any error, and the sampler's support would exclude the mode.
//
// contact_point_newton bisects *toward* this bound when a step would exceed it,
// so the bound must stay above the mode even when the mode itself is large.
// The bound is placed on the log axis and returned in the sampling coordinate
// y - origin; measured from the mode, the density's u_mode^2 (e^{y - m_g} - 1)^2
// term overflows at the same y.
double right_contact_bound(double m_g, double origin) {
  return (m_g - origin) + std::min(40.0, std::max(1.0, 354.0 - m_g));
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
// the region's coordinate system -- caller handles the mapping).
double sample_within_piece(const EnvelopePiece& piece) {
  const double u = R::runif(0.0, 1.0);
  switch (piece.type) {
    case PIECE_PLATEAU:
      return piece.a + u * (piece.b - piece.a);
    case PIECE_EXP_RIGHT: {
      // h(x) = h(a) * exp(slope * (x - a)),  slope < 0,  x >= a
      // F(x) is proportional to 1 - exp(slope * (x - a)), so
      // x = a + log(1 - u) / slope
      return piece.a + std::log1p(-u) / piece.slope;
    }
    case PIECE_EXP_LEFT: {
      // h(x) = h(a) * exp(slope * (x - a)),  slope > 0,  x <= a
      // F(x) is proportional to exp(slope * (x - a)), so
      // x = a + log(u) / slope
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
    case PIECE_EXP_BOUNDED_RIGHT:
      // h(x) = h(b) * exp(slope * (x - b)) on [a, b], slope > 0, with the mass
      // measured from b: 1 - F(x) = (1 - exp(slope (x - b))) / (1 - exp(-t)),
      // t = slope (b - a), so x = b + log1p(u * expm1(-t)) / slope, with
      // aux = expm1(-t) stored at setup.
      return piece.b + std::log1p(u * piece.aux) / piece.slope;
    case PIECE_TNEGHALF_LEFT:
      // Inverse-square hat on (-inf, pl].  With b = ppl, aux = om1,
      // F^-1: x = ppl - om1 / u; Gao & Wang (2025, Appendix B) print
      // ppl - om1^2 / U with U uniform on (0, om1), and u = U / om1.
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
    case PIECE_EXP_BOUNDED_RIGHT:
      return piece.base_log_dens + piece.slope * (x - piece.b);
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
  Rcpp::stop("rmhn: could not construct a sampling envelope for alpha = %g, "
             "beta = %g, gamma = %g. This is a defect; please report these "
             "parameter values to the package maintainer. [%s]",
             env.alpha, env.beta, env.gamma, stage);
}

// The failure that is not a defect: the law is too narrow for the grid of
// doubles its draws must lie on.  On the centred path every region returns the
// exact draw rounded once to the nearest double (see RtdrEnvelope::origin), so
// the only limit is that grid itself, and the sampler stops once the law's
// standard deviation at its mode, in the caller's coordinate, is no larger
// than the spacing of the doubles there: even the law rounded to double then
// has only a few atoms and says nothing about its shape.  With beta = 1 this
// happens from a tilt of 2^53 = 9.0e15, where the mode passes 2^52 and the
// doubles are a unit apart against a standard deviation of 1/sqrt(2).  On the
// shape axis the doubles at the mode are half a unit apart for alpha from
// 2^103 = 1.0e31, against a standard deviation that differs from 1/2 by a
// relative 1/(8 m), so there the outcome turns on the last bit of the
// computed standard deviation; from 2^105 = 4.1e31 they are a unit apart.
// That is at a small tilt.  A negative tilt widens the law relative to its
// mode, from 1/sqrt(2 alpha) towards the 1/sqrt(alpha - 1) of its Gamma
// limit, and under a strong one the stop comes between alpha - 1 = 2^104 and
// 2^106 (8.1e31), depending on where the mode falls between two powers of
// two; from 2^106 on it comes at any tilt (see setup_region_a).
//
// The message offers the law's normal limit.  With
//     log f(m + d) - log f(m) = (alpha - 1) [log1p(d/m) - d/m] - beta d^2,
// the normal of variance 1 / (-(log f)''(m)) = 1 / (2 beta + (alpha - 1)/m^2)
// -- 1/(2 beta) on the tilt axis and about 1/(4 beta) on the shape axis --
// leaves out a term (alpha - 1) d^3 / (3 m^3) in the log-density.  Wherever
// the stop fires the relative width sd/m is below 2.2e-16, and that term is
// at most (k^3/3) sd/m at k standard deviations: negligible on the tilt axis,
// and of order 1e-15 over +-3 standard deviations on the shape axis.
[[noreturn]] void envelope_unrepresentable(const mhn::RtdrEnvelope& env) {
  Rcpp::stop("rmhn: cannot sample at alpha = %g, beta = %g, gamma = %g: the "
             "distribution is narrower than the spacing of the double "
             "precision numbers around its mode (its standard deviation is no "
             "larger than one spacing), so draws could take only a few "
             "distinct values and could not represent its shape. Near its "
             "mode the law there is, to within about 1e-14 in log-density, a "
             "normal with mean m = mhn_mode(alpha, beta, gamma) and variance "
             "1/(2*beta + (alpha - 1)/m^2); use that instead.",
             env.alpha, env.beta, env.gamma);
}

// The other failure that is not a defect.  Every region builds its envelope
// for sqrt(beta) X, and under a strong negative tilt the law of sqrt(beta) X
// sits near alpha / |gamma_norm|, gamma_norm = gamma / sqrt(beta), with its
// mode near (alpha - 1) / |gamma_norm| in region A.  The envelope cannot be
// built once gamma_norm itself is -Inf (|gamma| / sqrt(beta) past the largest
// double), once the mode underflows (alpha - 1 in region A, or alpha in
// region BC, too small for the tilt), or, in region A, once the left
// tangent's slope, a few times |gamma_norm|, passes the largest double --
// although X itself may be well inside the double range.  In each case the
// factor exp(-beta x^2) of the density is 1 wherever the law has its mass:
// beta x^2 is about alpha^2 / gamma_norm^2, below 1e-500 for every alpha that
// reaches this point (from alpha - 1 = 2^106 on, region A refuses first,
// through envelope_unrepresentable, and the slope passes the largest double
// only for an alpha of a few).  So the law of X is, to double precision,
// Gamma(alpha, rate -gamma), which the message offers instead.
[[noreturn]] void envelope_out_of_range(const mhn::RtdrEnvelope& env) {
  Rcpp::stop("rmhn: cannot sample at alpha = %g, beta = %g, gamma = %g: "
             "gamma/sqrt(beta) is so far below zero, for this alpha, that "
             "the sampler, which works with sqrt(beta)*x, cannot place the "
             "distribution within the range of double precision numbers. To "
             "double precision the distribution is then a gamma distribution "
             "with shape alpha and rate -gamma; use rgamma(n, shape = alpha, "
             "rate = -gamma) instead.",
             env.alpha, env.beta, env.gamma);
}

// Spacing of the doubles just above m.
double grid_spacing(double m) {
  return std::nextafter(m, std::numeric_limits<double>::infinity()) - m;
}

// The mode in the caller's coordinate, as an unevaluated sum hi + lo, for the
// centred path.  It is the positive root of
//     2 beta x^2 - gamma x - c = 0,
// with c = alpha - 1 for the mode of X (region A) and c = alpha for
// e^{m_g}/sqrt(beta) (regions BC and D); c_lo is the part of c that the double
// c leaves out.  x0 = (mode in the beta = 1 coordinate)/sqrt(beta) is within a
// few units in the last place of the root, the rounding of gamma_norm and of
// sqrt(beta) included, and one Newton step from it on the residual, formed
// without rounding error from the caller's own beta, gamma and alpha
// (mhn::quadratic_root_step), gives the root to about twice double
// precision.  The products that step forms are the beta = 1 quantities m^2
// and gamma_norm m, reached through beta x0 = sqrt(beta) m, so none of them
// overflows at any parameters the representability test lets through;
// beyond it an overflow or a non-finite x0 yields a NaN or infinite hi,
// which that test refuses.
void caller_axis_mode(double beta, double gamma, double c, double c_lo,
                      double x0, double& hi, double& lo) {
  const double step = mhn::quadratic_root_step(beta, gamma, c, c_lo, x0);
  hi = x0 + step;
  lo = step - (hi - x0);
}

// The representability test behind envelope_unrepresentable, on the centred
// path, once env.x_mode_hi is set.  m is the mode in the beta = 1 coordinate
// (u_mode in regions BC and D, which differs from the mode of x by a relative
// 2/gamma_norm^2 there), and m / sqrt(2 m^2 + alpha - 1) is the standard
// deviation 1/sqrt(2 + (alpha - 1)/m^2) of the normal limit written so that
// it neither overflows nor divides by an underflowed m^2.  Written as a
// negated comparison, so that a NaN mode or spacing is refused too.
void check_representable(const mhn::RtdrEnvelope& env, double m) {
  const double sd = m / std::sqrt(2.0 * m * m + (env.alpha - 1.0))
                    / env.sqrt_beta;
  if (!(::grid_spacing(env.x_mode_hi) < sd)) ::envelope_unrepresentable(env);
}

// Classify the (alpha, gamma_norm) point into one of the three RTDR regions.
mhn::RtdrRegion classify_region(double alpha, double gamma_norm) {
  // Region A builds its hat around an interior mode in x.  At alpha = 1 with
  // gamma <= 0 the density's maximum is at the boundary x = 0 (f(0) = 1
  // there; the mode formula returns 0), so region A's construction does not
  // apply.  With gamma > 0, alpha = 1 has the interior mode gamma / (2 beta)
  // but goes to region BC as well, which Theorem 3.2 of Gao & Wang (2025)
  // (T_{-1/2}-concavity for lambda >= 1/2) covers.  The mode of g on the log
  // axis lies at a finite y for every alpha > 0 (its x = e^y is positive),
  // and the concavity test there is on the tilt alone, so alpha = 1 belongs
  // on that side of the boundary.  Gao & Wang themselves put lambda = 1 with
  // lambda > 1 (their Section 3.1, and Appendix B, `if(lambda >= 1)`), that
  // is in region A; the package departs from them here.
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
// Region setup helpers -- area/anchor formulas
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
// [a, b], slope > 0, anchored at b: h(x) = exp(base + slope (x - b)), whose
// integral is exp(base) (1 - exp(-t)) / slope with t = slope (b - a).  No
// term is larger than log h itself, however far a lies from b.
EnvelopePiece make_exp_bounded_right(double a, double b, double slope,
                                     double base_log_dens) {
  EnvelopePiece p;
  p.type = PIECE_EXP_BOUNDED_RIGHT; p.a = a; p.b = b; p.slope = slope;
  p.base_log_dens = base_log_dens;
  const double t = slope * (b - a);
  p.aux = std::expm1(-t);
  p.log_area = base_log_dens + mhn::log_abs_expm1(-t) - std::log(slope);
  return p;
}
// T_{-1/2} tangent pieces: the tangent hat of Gao & Wang (2025, Section 2,
// Eq. (4)) under the T_{-1/2} of their Eq. (6), valid by the
// T_{-1/2}-concavity of Section 3.2 (Theorem 3.2); coded as in their
// Appendix B.  The hat is h(y) = exp(logf_t) / (half_slope*(y - t) - 1)^2
// with half_slope = L'(t)/2.
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
// Region A setup [alpha > 1, log-concave on f(x)]
// ====================================================================
// The switch to the centred path of region A (see setup_region_a), at the
// mode m of the raw log-density log f(x) = (alpha-1) log x - x^2 +
// gamma_norm x at beta = 1.  The three terms at m each carry a rounding error
// of order DBL_EPSILON times their own size, so the error of the raw log f(m)
// is of order DBL_EPSILON times the sum of their magnitudes, and the centred
// path is taken once that sum passes 0.01 / DBL_EPSILON.  Shared by
// setup_region_a and mhn::rtdr_region_a_centred, so that the default method
// of rmhn (rebuild_cache in mhn_rmhn.cpp) makes the same test, on the same
// doubles, as the setup.
bool region_a_centred_at(double alpha, double gn, double m) {
  const double size =
      std::abs((alpha - 1.0) * std::log(m)) + m * m + std::abs(gn * m);
  return !(size * std::numeric_limits<double>::epsilon() <= 0.01);
}

void setup_region_a(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // mode of log f(x) = (alpha-1) log x - x^2 + gamma_norm * x, i.e. the
  // positive root of 2x^2 - gamma_norm x - (alpha-1) = 0 at beta = 1.
  const double c = alpha - 1.0;
  const double m = mhn::positive_root(1.0, gn, c);
  if (!(m > 0.0)) {
    // positive_root forms the root without squaring gamma_norm, so for
    // alpha > 1 it is 0 only where gamma_norm is -Inf (the caller's
    // gamma / sqrt(beta) is past the largest double) or where the root,
    // about (alpha - 1) / |gamma_norm|, lies below the smallest positive
    // double.  The law of X is then Gamma(alpha, rate |gamma|) to double
    // precision.  Its relative width -- the standard deviation
    // 1/sqrt(2 + (alpha - 1)/m^2) of its normal limit, over the mode m -- is
    // below 1/sqrt(alpha - 1), and the doubles at any positive x are more
    // than x 2^-53 apart, so from alpha - 1 = 2^106 on the law is narrower
    // than the spacing wherever its mode lies.  Below that the envelope
    // cannot be placed (see envelope_out_of_range).
    if (!(c < std::ldexp(1.0, 106))) ::envelope_unrepresentable(env);
    ::envelope_out_of_range(env);
  }

  auto log_f_raw = [alpha, gn](double x) -> double {
    if (x <= 0.0) return -std::numeric_limits<double>::infinity();
    return (alpha - 1.0) * std::log(x) - x * x + gn * x;
  };
  auto dlog_f_raw = [alpha, gn](double x) -> double {
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    return (alpha - 1.0) / x - 2.0 * x + gn;
  };
  const double log_f_mode_raw = log_f_raw(m);

  // Region BC's conditioning problem (see setup_region_bc) arises here too, in
  // x rather than on the log axis.  At the mode -x^2 + gamma_norm x is
  // gamma_norm^2/4 - (x - gamma_norm/2)^2, so every value of log f carries
  // gamma_norm^2/4 -- 2.5e19 at gamma_norm = 1e10 -- while the construction
  // reads O(1) differences out of it: the contact drop (delta = 1), the
  // tangent intersections p_l and p_r, the piece areas, and the accept test.
  // Measured against 300-bit arithmetic at alpha = 2, gamma_norm = 1e10, a
  // drop whose true value is -25 comes back as -4096 from the plain
  // difference.
  //
  // Past the switch, where the rounding error of log f(m) can exceed a
  // hundredth of delta, measure from the mode, in the ordinate and in the
  // abscissa alike: the envelope is built, stored and sampled in s = x - m.
  // That error is of order the machine epsilon times the sum of the
  // magnitudes of the three terms of log f(m), (alpha-1) log m, m^2 and
  // gamma_norm m, and the switch is placed where that product passes 0.01
  // (see region_a_centred_at).  |log f(m)| itself does not measure the
  // error: at a tilt near -alpha/e the mode is near e, where terms of order
  // alpha cancel, and a small log f(m) still carries an error of order
  // alpha DBL_EPSILON.
  // Tested on |log f(m)| instead, that pocket would keep the raw
  // expressions, whose draws drift visibly from a shape of about 3e15 (the
  // mean is 0.75 standard deviations off at 1e17).  The sum of the
  // magnitudes is at least alpha - 2 at every tilt, so the path is centred
  // whatever the tilt from a shape of 0.01 / DBL_EPSILON + 2, about 4.5e13,
  // and in the pocket from about 2.3e13.
  // With the mode equation 2m - gamma_norm = (alpha-1)/m,
  //     log f(m + s) - log f(m) = (alpha-1) [ log1p(s/m) - s/m ] - s^2
  //     (log f)'(m + s)         = -s [ (alpha-1)/((m + s) m) + 2 ]
  // log1p(u) - u <= 0 for every u > -1 and alpha > 1 throughout region A, so
  // both terms carry the same sign and nothing cancels.  In x itself the
  // doubles near the mode are ulp(m) apart, a tenth of a standard deviation
  // from a tilt of 1.1e15 or a shape of 1.6e29, and a sampler that proposes x
  // and accepts or rejects at the double a proposal rounds to gives each double
  // only approximately the law's mass on the cell that rounds to it; in s the
  // accept test sees the proposal at full precision.  The pair of formulas is
  // exact for the law whose mode is the double m -- the tilt 2m - (alpha-1)/m
  // -- and a relative error of order 1e-16 in m changes its log-density by a
  // relative 1e-16 over the range where the mass is, so s is a draw of X - mode
  // in the beta = 1 coordinate, and it is placed at the mode computed in the
  // caller's coordinate (see RtdrEnvelope::origin).  Below the switch the raw
  // expressions are used as they stand.
  //
  // The switch is written as a negated comparison so that a term that has
  // overflowed to Inf, as gamma_norm m does once a positive gamma_norm passes
  // about 1.9e154, or a NaN takes the centred path and meets the
  // representability test below rather than a construction built on NaN.
  const bool centred = ::region_a_centred_at(alpha, gn, m);
  const double origin = centred ? m : 0.0;
  const double mode_s = m - origin;          // the mode in s: 0 when centred
  const double zero_s = centred ? -m : 0.0;  // x = 0 in s
  // The curvature term (alpha-1)/m^2 of -(log f)''(m).  It passes the largest
  // double once the mode is below about 7e-155 sqrt(alpha - 1), that is at
  // gamma_norm below about -1.3e154 sqrt(alpha - 1): m is then about
  // (alpha-1)/|gamma_norm|, and (alpha-1)/m^2 about gamma_norm^2/(alpha-1).
  // There the terms of the derivative are regrouped.  Centred it is formed as
  // -2 (s/(m + s)) ((alpha-1)/2 / m) - 2 s: near the mode s/(m + s) is of
  // order 1/sqrt(alpha - 1), and (alpha-1)/2 / m is about |gamma_norm|/2,
  // at most half the largest double, so neither factor overflows, nor does
  // their product while the true slope is finite.  (Grouped as
  // (s/m) ((alpha-1)/(m + s)), the second factor exceeds |gamma_norm| by a
  // relative s/m left of the mode and overflows once |gamma_norm| is within
  // that factor of the largest double.)  Raw it is formed as
  // (alpha - 1 + gamma_norm x)/x - 2x, since (alpha-1)/x alone overflows for
  // x below (alpha-1)/DBL_MAX, which half the mode reaches once |gamma_norm|
  // passes DBL_MAX/2; the numerator cancels near the mode as
  // (alpha-1)/x + gamma_norm does, to the same absolute error.
  const double curv = (c / m) / m;
  const bool curv_overflows =
      !(mhn::opaque_value(curv) <= std::numeric_limits<double>::max());
  auto log_f = [&](double s) -> double {
    if (!centred) return log_f_raw(s);
    if (s <= -m) return -std::numeric_limits<double>::infinity();
    return (alpha - 1.0) * mhn::log1p_minus_u(s / m) - s * s;
  };
  auto dlog_f = [&](double s) -> double {
    if (!centred) {
      if (curv_overflows && s > 0.0) return (c + gn * s) / s - 2.0 * s;
      return dlog_f_raw(s);
    }
    if (s <= -m) return std::numeric_limits<double>::infinity();
    if (curv_overflows) {
      return -2.0 * ((s / (m + s)) * ((0.5 * c) / m)) - 2.0 * s;
    }
    return -s * ((alpha - 1.0) / ((m + s) * m) + 2.0);
  };
  const double log_f_mode = centred ? 0.0 : log_f_mode_raw;
  env.mode = mode_s;
  env.origin = origin;
  env.centred = centred;
  env.log_dens_mode = log_f_mode;

  // Simplified-envelope condition, stated with envelope (10) in Gao & Wang
  // (2025), Section 3.1; their Theorem 3.1 validates (10) and (11).
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

  // Contact-point starts.  Gao & Wang (2025, Section 2) recommend starting
  // their Eq. (8) iteration at m +- inc, inc = sqrt(-2 delta / L''(m)) with
  // -L''(m) = (alpha-1)/m^2 + 2, as region BC does, and that start is paired
  // here with a step tolerance on the same scale, 1e-6 inc.  The fixed starts
  // m + max(1, m/2) on the right and m/2 on the left, with an absolute step
  // tolerance of 1e-10, are used except in two places where they fail:
  //
  //  * Past the switch.  The right start's offset is of order m, and there
  //    the ordinate behaves like -d^2 near the mode, on which one Newton step
  //    only halves d.  Reaching the band takes about log2(m/2) steps -- 31 at
  //    gamma_norm = 1e10, 38 at 1e12 -- against a budget of 30, and the point
  //    returned has a drop of hundreds rather than 1.  The hat still
  //    dominates (log f is concave, so every tangent does), but the plateau
  //    widens to the contact points, and the acceptance falls by an order of
  //    magnitude for each factor of ten in gamma_norm past about 1e10, until
  //    the 1000-retry budget runs out.  Below the switch m^2, one of the three terms of
  //    log f(m), is at most 0.01/epsilon, so m <= 6.7e6 and the walk in
  //    costs at most 23 steps.
  //  * At a strong negative tilt, gamma_norm < -1e10 sqrt(alpha).  The law is
  //    then close to a gamma law of rate -gamma_norm, whose standard
  //    deviation sqrt(alpha)/(-gamma_norm) is below the absolute tolerance
  //    1e-10, so a Newton step can be shorter than the tolerance while the
  //    iterate is still many standard deviations from the band, and the
  //    search stops there: measured acceptance is 0.086 at alpha = 1000,
  //    gamma_norm = -1e15, and 0.003, with NaN draws, at alpha = 1e6,
  //    gamma_norm = -1e18.  Further out -- from gamma_norm of about -1e25 at
  //    alpha = 1.01 and -1e27 at alpha = 2 -- the Newton step from the right
  //    start, a unit away, is lost to rounding and cut back to a halving, and
  //    thirty of them leave the contact point 2^-30 from the mode.  The
  //    tangent intersection p_r then cancels to 0, left of the mode, and the
  //    hat fails to cover the law.  Measured over alpha from 1 + 1e-8 to 1e12
  //    and gamma_norm from -1e2 to -1e30, the contact points from the fixed
  //    starts first leave the band where that standard deviation is about
  //    1.05e-10, just above the threshold, near alpha = 36 and
  //    gamma_norm = -5.7e10: the left contact's drop there is 2.53 to 2.68,
  //    against the band's 2.49.  The hat, a tangent to a concave log f, still
  //    dominates there, and the acceptance measured is 0.82, against 0.88
  //    just past the switch.
  //
  // The default method of rmhn reads which start is taken from
  // mhn::rtdr_region_a_centred and mhn::rtdr_region_a_strong_negative_tilt
  // (mhn_rtdr.h), which make the two tests above: the walk in from the fixed
  // starts, about log2 m Newton steps, makes this setup dearer than that of
  // Sun et al.'s (2023) Algorithm 3 (gamma < 0) from a shape of about 1e3.
  // env.quick_start records the outcome.
  const double neg_ddlog_f_mode = curv + 2.0;
  // Where curv overflows the same offset is m sqrt(2 / (2 m^2 + alpha - 1)).
  const double inc = curv_overflows
      ? m * std::sqrt(2.0 / (2.0 * m * m + c))
      : std::sqrt(2.0 / neg_ddlog_f_mode);
  const bool inc_usable =
      centred || mhn::rtdr_region_a_strong_negative_tilt(alpha, gn);
  env.quick_start = inc_usable;
  const double newton_tol = inc_usable ? 1e-6 * inc : 1e-10;
  if (centred) {
    ::caller_axis_mode(env.beta, env.gamma, c, (alpha - c) - 1.0,
                       m / env.sqrt_beta, env.x_mode_hi, env.x_mode_lo);
    ::check_representable(env, m);
  }

  // Right contact point: t_r > m, slope_r = L'(t_r) < 0.
  const double t_r_init = inc_usable ? mode_s + inc : m + std::max(1.0, 0.5 * m);
  env.t_r = ::contact_point_newton(t_r_init, log_f_mode, 1.0, log_f, dlog_f,
                                   0.46, 2.49, /*max_iter=*/30, newton_tol,
                                   /*t_min=*/mode_s, /*t_max=*/std::numeric_limits<double>::infinity());
  env.slope_r = dlog_f(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region A right contact");
  }
  env.p_r = env.t_r + (log_f_mode - log_f(env.t_r)) / env.slope_r;
  // p_r is the plateau's right edge and the right tail's anchor.  For an
  // exact tangent it lies right of the mode; one at or left of it means the
  // intersection has been lost to rounding, and the tail, which starts at the
  // mode's height, would then pass below the target at the mode.  A
  // non-finite p_r reaches make_plateau's log(b - a) and poisons every piece
  // area, and so every draw.  Neither would raise an error anywhere else: the
  // simplified branch below checks no other endpoint.
  if (!(env.p_r > mode_s) || !std::isfinite(env.p_r)) {
    ::envelope_failure(env, "region A right intersection");
  }

  std::vector<EnvelopePiece> pieces;
  pieces.reserve(3);

  if (!simplified) {
    // Left contact point: 0 < t_l < m, slope_l = L'(t_l) > 0.
    // The fixed start is half the mode, inside (0, m) at any scale.  A start
    // with an absolute floor, max(m/2, 1e-6), would lie right of the mode
    // whenever the mode is below 2e-6 (it is 3e-7 at alpha = 4,
    // gamma_norm = -1e7), and since the search tests the acceptance band
    // before it clamps into [t_min, t_max], it could return a point on the
    // wrong side with a non-positive slope.  The Eq. (8) start is
    // m - min(inc, m/2), kept strictly inside (0, m) for the same reason: a
    // bare m - inc lands at or below zero whenever the mode is smaller than
    // the mode's own width.
    const double t_l_init =
        inc_usable ? mode_s - std::min(inc, 0.5 * m) : 0.5 * m;
    env.t_l = ::contact_point_newton(t_l_init, log_f_mode, 1.0, log_f, dlog_f,
                                     0.46, 2.49, /*max_iter=*/30, newton_tol,
                                     /*t_min=*/zero_s, /*t_max=*/mode_s);
    if (env.t_l <= zero_s) env.t_l = zero_s + 0.5 * m;  // defensive
    env.slope_l = dlog_f(env.t_l);
    // Left of the mode the slope is about (alpha-1)/x - |gamma_norm|, and at
    // the contact point of a small alpha - 1 that is a few times |gamma_norm|
    // -- from 1.2 to 7 times over the acceptance band at alpha = 3 -- so it
    // can pass DBL_MAX once |gamma_norm| is within that factor of it: the
    // edge of the range of doubles, not a defect.
    if (curv_overflows && std::isinf(env.slope_l) && env.slope_l > 0.0) {
      ::envelope_out_of_range(env);
    }
    if (!(env.slope_l > 0.0) || !std::isfinite(env.slope_l)) {
      ::envelope_failure(env, "region A left contact");
    }
    env.p_l = env.t_l + (log_f_mode - log_f(env.t_l)) / env.slope_l;
    // An exact tangent meets the mode's height inside (0, m), and with
    // p_r > m that keeps the plateau non-empty.  Negated rather than written
    // p_l <= 0 || p_l >= m so a non-finite endpoint is caught: every ordinary
    // comparison is false for NaN, and a NaN p_l would otherwise walk past
    // this guard and poison every piece area.
    if (!(env.p_l > zero_s) || !(env.p_l < mode_s)) {
      ::envelope_failure(env, "region A intersection");
    }
    // Left tangent: [0, p_l] in x, h(p_l) = f(m).  Its mass lies within a few
    // 1/slope_l of p_l, while its left end x = 0 lies m away, which past the
    // switch is millions of standard deviations (m/sd = sqrt(2 m^2 + alpha - 1)
    // in the beta = 1 coordinate), so there it is anchored at p_l (see
    // PIECE_EXP_BOUNDED_RIGHT).  Anchored at 0 instead, the base is
    // log_f_mode - slope_l * p_l.
    if (centred) {
      pieces.push_back(make_exp_bounded_right(zero_s, env.p_l, env.slope_l,
                                              log_f_mode));
    } else {
      pieces.push_back(make_exp_bounded(
          0.0, env.p_l, env.slope_l, log_f_mode - env.slope_l * env.p_l));
    }
  }

  // Plateau: [plateau_a, p_r], h = f(m).
  const double plateau_a = simplified ? zero_s : env.p_l;
  pieces.push_back(make_plateau(plateau_a, env.p_r, log_f_mode));

  // Right tangent: [p_r, +inf), slope_r < 0, h(p_r) = f(m).
  pieces.push_back(make_exp_right(env.p_r, env.slope_r, log_f_mode));

  finalize_pieces(env, pieces);
}

// ====================================================================
// Region BC setup [alpha <= 1, T_{-1/2}-concave on g(y) = exp(alpha*y
//   - exp(2y) + gamma_norm * exp(y)), y in R; for gamma_norm <= 0, where g is
//   log-concave, the hat is built from T_0 tangents instead]
// ====================================================================
void setup_region_bc(mhn::RtdrEnvelope& env) {
  const double alpha = env.alpha;
  const double gn = env.gamma_norm;

  // Mode in u = exp(y): 2 u^2 - gamma_norm u - alpha = 0,
  //   u_+ = (gamma_norm + sqrt(gamma_norm^2 + 8 alpha)) / 4 > 0.
  // positive_root forms it without squaring gamma_norm, and it is 0 only
  // where gamma_norm is -Inf or the root, about alpha / |gamma_norm|, lies
  // below the smallest positive double (see envelope_out_of_range).
  const double u_mode = mhn::positive_root(1.0, gn, alpha);
  if (!(u_mode > 0.0)) ::envelope_out_of_range(env);
  const double m_g = std::log(u_mode);

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
  // and p_l are decided by the last bit of exp(), and the outcome depends on
  // the platform's libm: the same triple can build a NaN-bounded envelope on
  // one and trip the p_l >= p_r guard on another, non-monotonically in gamma.
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
  // switch the raw expressions are used as they stand.
  //
  // Past the switch the abscissa is measured from the mode as well: the
  // envelope is built, stored and sampled in s = y - m_g.  In y itself the
  // doubles near m_g ~ log(gamma_norm / 2) are ulp(m_g) apart against a law
  // about 1.4 / gamma_norm wide, so a draw exp(y) could take only values
  // x ulp(log x) apart -- 0.018 at gamma_norm = 1e13, eighteen spacings of the
  // doubles at the mode and 2.5% of the standard deviation -- and past about
  // 3e14 the contact structure falls inside one ulp of m_g.  In s everything
  // near the mode keeps full relative precision and the formulas above hold
  // with d = s.  They are exact for the law whose log-axis mode is the double
  // u_mode, and a relative error of order 1e-16 in u_mode only rescales s by
  // as much, so s is a draw of Y - m_g; it is placed at the mode computed in
  // the caller's coordinate (see RtdrEnvelope::origin), which makes the draw
  // the exact draw rounded once.
  //
  // The switch is a negated comparison so that a raw peak that has overflowed
  // to NaN (a positive gamma_norm past about 1.3e154) takes the centred path
  // and meets the representability test below.
  const bool centred =
      !(std::abs(log_g_mode_raw) * std::numeric_limits<double>::epsilon() <=
        0.01 * ((gn > 0.0) ? std::log(4.0) : 1.0));
  const double origin = centred ? m_g : 0.0;
  const double mode_s = m_g - origin;       // the mode in s: 0 when centred
  const double um2 = u_mode * u_mode;
  auto log_g = [&](double s) -> double {
    if (!centred) return log_g_raw(s);
    const double E = std::expm1(s);
    return -(um2 * E * E + alpha * (E - s));
  };
  auto dlog_g = [&](double s) -> double {
    if (!centred) return dlog_g_raw(s);
    return -std::expm1(s) * (alpha + 2.0 * um2 * std::exp(s));
  };
  const double log_g_mode = centred ? 0.0 : log_g_mode_raw;
  env.mode = mode_s;
  env.origin = origin;
  env.centred = centred;
  env.u_mode = u_mode;
  env.log_dens_mode = log_g_mode;
  env.simplified = false;  // BC envelope is always 3-piece

  // g(y) is log-concave iff gamma_norm <= 0 (Gao & Wang 2025, Lemma 3.3,
  // whose beta is -gamma_norm).  For gamma_norm <= 0 use the T_0
  // (log-tangent) hat with acceptance band [0.46, 2.49] and delta = 1.  For
  // gamma_norm > 0 the density is only T_{-1/2}-concave (their Theorem 3.2)
  // -- log g is convex on (-inf, log(gamma_norm/4)) -- so a log-tangent hat
  // would fail to dominate there; use the T_{-1/2} (inverse-square) hat with
  // band [0.93, 1.99] and delta = log 4.
  const bool tneghalf = (gn > 0.0);
  const double delta = tneghalf ? std::log(4.0) : 1.0;
  const double band_lo = tneghalf ? 0.93 : 0.46;
  const double band_hi = tneghalf ? 1.99 : 2.49;

  // Recommended contact-point start (Gao & Wang 2025, Section 2, stated after
  // their Eq. (8)):
  // inc = sqrt(-2*delta / L''(m_g)), with L''(m_g) = -u_mode*(4 u_mode - gn) < 0.
  const double ddlog_g_mode = -u_mode * (4.0 * u_mode - gn);
  const double inc = std::sqrt(-2.0 * delta / ddlog_g_mode);
  // On the centred path a draw is the exact draw rounded once, so the only
  // limit left is the grid of doubles at the mode, in the caller's coordinate
  // (see check_representable).  With gamma_norm <= 0 the centred path is not
  // reached: for any positive double u_mode, |log g(m_g)| is below 750 there
  // (alpha |m_g| <= 745, and u_mode^2 + |gamma_norm| u_mode <= alpha by the
  // mode equation), and a u_mode of 0 has been refused above.
  if (centred) {
    if (!(gn > 0.0)) ::envelope_failure(env, "region BC mode");
    ::caller_axis_mode(env.beta, env.gamma, alpha, 0.0,
                       u_mode / env.sqrt_beta, env.x_mode_hi, env.x_mode_lo);
    ::check_representable(env, u_mode);
  }
  // A step tolerance on the scale of the contact offset, as in region D, so
  // that it means the same at every tilt.
  const double newton_tol = centred ? 1e-6 * inc : 1e-10;

  // Left contact point: t_l < m_g, slope_l > 0.
  env.t_l = ::contact_point_newton(mode_s - inc, log_g_mode, delta,
                                   log_g, dlog_g, band_lo, band_hi,
                                   /*max_iter=*/30, newton_tol,
                                   /*t_min=*/-std::numeric_limits<double>::infinity(),
                                   /*t_max=*/mode_s);
  env.slope_l = dlog_g(env.t_l);
  if (!(env.slope_l > 0.0) || !std::isfinite(env.slope_l)) {
    ::envelope_failure(env, "region BC left contact");
  }
  const double log_g_tl = log_g(env.t_l);

  // Right contact point: t_r > m_g, slope_r < 0.
  env.t_r = ::contact_point_newton(mode_s + inc, log_g_mode, delta,
                                   log_g, dlog_g, band_lo, band_hi,
                                   /*max_iter=*/200, newton_tol,
                                   /*t_min=*/mode_s,
                                   /*t_max=*/::right_contact_bound(m_g, origin));
  env.slope_r = dlog_g(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region BC right contact");
  }
  const double log_g_tr = log_g(env.t_r);

  std::vector<EnvelopePiece> pieces;
  pieces.reserve(3);
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
    // caught: `>=` is false for NaN, and a NaN p_l would otherwise walk past
    // this guard into the plateau's width and poison every piece area.
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
    // For a tiny alpha the left search starts between 1 / sqrt(alpha) (at
    // gamma_norm = 0) and sqrt(2 / alpha) (for gamma_norm well below
    // -sqrt(alpha)) left of the mode -- 1e20 to 1.4e20 at alpha = 1e-40 --
    // and stops there at once, the
    // slope alpha being below its step tolerance.  The tangent there is exact,
    // but the plain p_l adds back to t_l a term of the same size, so it
    // carries an error of about one spacing of the doubles at t_l, 16384 at
    // alpha = 1e-40, against a true p_l between half a unit and a unit left of
    // the mode and a p_r within 40 right of it.  An error e in p_l moves the
    // hat by alpha e in log-density, about 2e-16 sqrt(alpha) here, which no
    // draw can detect, but once it carries p_l past p_r the plateau is lost.
    // There p_l is formed instead with the terms alpha (m_g - t_l) cancelled:
    // from L(y) = alpha y - e^{2y} + gamma_norm e^y,
    //   p_l - m_g = [ (m_g - t_l)(2 e^{2 t_l} - gamma_norm e^{t_l})
    //                 + (gamma_norm u_m - u_m^2)
    //                 - (gamma_norm e^{t_l} - e^{2 t_l}) ] / L'(t_l),
    // in which, with e^{t_l} = 0, only the middle term, between -alpha and
    // -alpha/2 by the mode equation, is left.  (m_g - t_l) e^{t_l} tends to 0
    // as t_l falls to -Inf, where the search starts once inc overflows, for
    // alpha below about 1e-308.  Only the uncentred path reaches this branch
    // (see above), so m_g and t_l are both on the log axis.
    if (!(env.p_l < env.p_r) && !centred) {
      const double et = std::exp(env.t_l);
      const double lever = (et > 0.0) ? (m_g - env.t_l) * et : 0.0;
      env.p_l = m_g + (lever * (2.0 * et - gn) + (gn * u_mode - um2)
                       - (gn * et - et * et)) / env.slope_l;
    }
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
//   inflection-point envelope on log g(y); every piece is log-linear].
//
// Implements both envelope (14) and envelope (15) of Gao & Wang (2025), with
// automatic switching on the criterion log(g(m_g)/g(y*)) > 2.49.  In their
// notation the secant chain ends at \hat y^*, which is y* for envelope (14)
// and the dual point \dot t_l of the left contact point t_l (their Lemma 4.5,
// Theorem 4.3) for envelope (15); the code writes the latter t_dot_l.
//   - envelope (14) (default): \hat y^* = y*, plateau spans (y*, p_r] (the
//     paper's [y*, p_r], half open here so that the pieces do not overlap).
//     Used when the log-density drop from mode to inflection is moderate.
//   - envelope (15): \hat y^* = t_dot_l, an extra left-tangent piece is
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
  const double disc = mhn::root_discriminant(1.0, gn, alpha);   // S
  const double u_mode = mhn::positive_root(1.0, gn, alpha);
  const double m_g = std::log(u_mode);
  const double y_star = std::log(gn / 4.0);

  auto log_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha * y - ey * ey + gn * ey;
  };
  auto dlog_g_raw = [alpha, gn](double y) -> double {
    const double ey = std::exp(y);
    return alpha - 2.0 * ey * ey + gn * ey;
  };
  const double log_g_mode_raw = log_g_raw(m_g);

  // Region D works on the same g(y) as region BC and has the same
  // conditioning problem.  Near the mode -e^{2y} + gamma_norm e^y =
  // gamma_norm^2/4 - (e^y - gamma_norm/2)^2, so every value of log g carries
  // gamma_norm^2/4 -- 2.5e19 at gamma_norm = 1e10, one ulp of which is 4096 --
  // while everything the envelope reads out of it is an O(1) *difference* of
  // two such values: the contact drop (delta = log 4), the tangent/plateau
  // intersections p_l and p_r, the secant slopes, the between-piece area
  // ratios exponentiated in finalize_pieces, and the accept/reject test.
  // Measured on the raw expressions: at gamma_norm = 1e9 and above the
  // numerator of p_l rounds to exactly zero, so p_l == t_l bit for bit, and at
  // 1e10 p_r lands to the *left* of p_l; without the plateau the sampler's
  // support would not contain the mode (sd/true = 0.0146, KS D = 1.000 against
  // the exact Gaussian limit).
  //
  // The remedy is region BC's: with d = y - m_g, E = expm1(d), and the mode
  // equation 2 u_m^2 - gamma_norm u_m - alpha = 0,
  //     log g(y) - log g(m_g) = -( u_m^2 E^2 + alpha (E - d) )
  //     (log g)'(y)           = -E ( alpha + 2 u_m^2 e^d )
  // Both only ever add same-signed terms (E - d = expm1(d) - d >= 0), so the
  // drop stays exact at any tilt.  Region D's hat is log-linear (T_0)
  // throughout, as in Gao & Wang (2025), Eqs. (14)-(15).  Its contact points
  // are searched with delta = log 4 and accepted in [0.93, 1.99], which lies
  // inside the T_0 band [0.46, 2.49] of their Eq. (7) that Eqs. (14)-(15)
  // call for.
  //
  // The switch is region BC's rule -- ulp(log g(m_g)) against a fraction of
  // delta -- but with a tenth of BC's fraction.  On the raw path the draws
  // drift from the exact Gaussian limit from gamma_norm of about 2e7, and
  // BC's own constant would put the switch only a factor of 1.25 below that;
  // a tenth of it puts the switch at 5.0e6, a factor of four in tilt and
  // sixteen in ulp below it, and still five times above the largest tilt the
  // package's own test grid reaches.  Region D earns the wider margin because it has strictly more
  // difference sites than BC's three-piece table -- the dual point, the secant
  // chain and the left-tangent/plateau join on top of the same contact search
  // and the same area ratios -- and because the cost of switching early is
  // nil: at gamma_norm = 5e6 the two constructions already agree to a part in
  // 1e4.  Below the switch the raw expressions are used as they stand.
  //
  // Past the switch the abscissa is measured from the mode too, s = y - m_g,
  // for region BC's reason (see setup_region_bc): every position below --
  // contact points, dual point, breakpoints, piece bounds -- is in s.  The
  // switch is negated for region BC's reason as well.
  const double delta = std::log(4.0);
  const bool centred =
      !(std::abs(log_g_mode_raw) * std::numeric_limits<double>::epsilon() <=
        0.001 * delta);
  const double origin = centred ? m_g : 0.0;
  const double mode_s = m_g - origin;       // the mode in s: 0 when centred
  // y* - m_g = log((gamma_norm/4) / u_mode), formed as one quotient.
  const double y_star_s = centred ? std::log((gn / 4.0) / u_mode) : y_star;
  const double um2 = u_mode * u_mode;
  auto log_g = [&](double s) -> double {
    if (!centred) return log_g_raw(s);
    const double E = std::expm1(s);
    return -(um2 * E * E + alpha * (E - s));
  };
  auto dlog_g = [&](double s) -> double {
    if (!centred) return dlog_g_raw(s);
    return -std::expm1(s) * (alpha + 2.0 * um2 * std::exp(s));
  };
  const double log_g_mode = centred ? 0.0 : log_g_mode_raw;
  const double log_g_y_star = log_g(y_star_s);
  env.mode = mode_s;
  env.origin = origin;
  env.y_star = y_star_s;
  env.centred = centred;
  env.u_mode = u_mode;
  env.log_dens_mode = log_g_mode;
  env.simplified = false;

  // Contact-point starts on the curvature scale recommended by Gao & Wang
  // (2025, Section 2) for their Eq. (8) iteration, as region BC uses, with
  // L''(m_g) = -u_mode (4 u_mode - gamma_norm) < 0, past the switch.  The
  // fixed starts used below it -- the midpoint of [y_star, m_g] on the left
  // and m_g + 1 on the right -- sit a fixed O(1) distance from the mode on the
  // log axis while the contact points sit 2 sqrt(log 4) / gamma_norm away,
  // 2.4e-10 at gamma_norm = 1e10.  Newton halves that distance per step, so
  // from them it would need log2(0.15 gamma_norm) steps on the left: 30.5 at
  // 1e10, 37 at 1e12, against a budget of 30.  Even in exact arithmetic the
  // left search would then return a point whose drop is 2.28 at 1e10 and
  // 12731 at 1e12, outside the [0.93, 1.99] band it is supposed to land in.
  // The right search has a budget of 200 but would be stopped instead by the
  // *absolute* step tolerance, which exceeds the whole contact structure once
  // gamma_norm passes ~1.4e10; past the switch the tolerance is scaled to the
  // mode's own width so that it means the same thing at every tilt.  From the
  // curvature start the first evaluation is already inside the band, so
  // neither limit is reached: measured acceptance is 0.87 flat from
  // gamma_norm = 1e6 to 1e14, against 0.24 at 1e10 and 0.27 at 1e12 from the
  // fixed starts.
  const double ddlog_g_mode = -u_mode * (4.0 * u_mode - gn);
  const double inc = std::sqrt(-2.0 * delta / ddlog_g_mode);
  const double newton_tol = centred ? 1e-6 * inc : 1e-10;
  // In s the construction has no floor of its own; the grid of the draws, the
  // doubles at the mode in the caller's coordinate, is the limit, as in region
  // BC.  Region D has gamma_norm > 0.
  if (centred) {
    ::caller_axis_mode(env.beta, env.gamma, alpha, 0.0,
                       u_mode / env.sqrt_beta, env.x_mode_hi, env.x_mode_lo);
    ::check_representable(env, u_mode);
  }

  // Envelope selection.  log(g(m_g)/g(y_star)) > 2.49 selects envelope (15).
  // Centred, log_g_y_star is already the negated drop and log_g_mode is zero,
  // so the test reads the same in both coordinate systems.  It is well
  // conditioned either way: the drop is gamma_norm^2/16, not O(1).
  const bool use_env_15 = (log_g_mode - log_g_y_star) > 2.49;

  // Compute rho according to the selected envelope.  Envelope (15) also
  // derives the contact point t_l (in the log-concave half [y_star, m_g]) and
  // the tangent intersection p_l for the extra left-tangent piece.
  double rho;
  double t_l_d = 0.0, slope_l_d = 0.0;
  double log_g_t_l = 0.0, p_l_d = 0.0;
  // rho is L_g(\hat y^*) - alpha \hat y^*.  Since L_g(y) - alpha y is exactly
  // -e^{2y} + gamma_norm e^y, that is
  //     rho = v (gamma_norm - v),      v = e^{\hat y^*} <= gamma_norm / 2,
  // a product of two positive factors.  Read instead as a difference of two
  // log-density values it subtracts two numbers of size gamma_norm^2/4 to
  // leave one of size gamma_norm: measured against the theoretical value it
  // comes out 41x too large at gamma_norm = 1e10 and 6e5 times too large at
  // 1e14.  It also cannot be read off a *centred* log g at all -- that would
  // give rho ~ -gamma_norm^2/4, tripping the rho <= 0 guard below -- so the
  // closed form is mandatory once the switch fires rather than merely better.
  // For envelope (14), v = e^{y_star} = gamma_norm/4 and it reduces to
  // 3 gamma_norm^2 / 16.
  if (!use_env_15) {
    // Envelope (14): \hat y^* = y_star.  Closed form rho = 3 gn^2 / 16.
    rho = centred ? (gn / 4.0) * (gn - gn / 4.0)
                  : (log_g_y_star - alpha * y_star);
  } else {
    // Envelope (15): contact-point search inside [y_star, m_g] (log-concave
    // sub-region of g).  contact_point_newton clamps the start into
    // [y_star, m_g] itself, so the flat-mode case the mid-point start guards
    // against is covered by the curvature start too.
    const double y_0 = std::log(gn / 2.0);   // dual point of -inf
    const double t_l_init = centred ? mode_s - inc : 0.5 * (y_star + m_g);
    t_l_d = ::contact_point_newton(t_l_init, log_g_mode, delta,
                                   log_g, dlog_g, 0.93, 1.99,
                                   30, newton_tol, y_star_s, mode_s);
    // The dual point \dot t_l of t_l (Gao & Wang 2025, Lemma 4.5 and
    // Section 4, where it is log(-beta/2 - e^{t_l}) with their beta equal to
    // -gamma_norm).  v = e^{t_dot_l} = gamma_norm/2 - e^{t_l} is
    // sqrt(drop - alpha (E - d)) - 2 alpha/(gamma_norm + S), with drop the
    // contact drop the search accepted and E = expm1(d) below: about the
    // square root of that drop, 1.1774 once the curvature start lands in the
    // band, and between about 0.91 and 1.39 over the region.  So it is O(1) at
    // every tilt, and written literally it subtracts two numbers of size
    // gamma_norm/2 and loses a digit per decade of tilt.  With d = t_l - m_g,
    // e^{t_l} = u_mode (1 + expm1(d)), and the conjugate form of the mode
    // equation gives gamma_norm/2 - u_mode = (gamma_norm - S)/4 =
    // -2 alpha / (gamma_norm + S) exactly, so
    //     v = -u_mode expm1(d) - 2 alpha / (gamma_norm + S)
    // adds two terms of known sign (expm1(d) < 0 here, and the subtrahend is
    // the tiny one).  The uncentred guard t_l < log(gamma_norm/2) is exactly
    // v > 0, and written this way it also traps a non-finite t_l.
    const double v = centred
        ? (-u_mode * std::expm1(t_l_d - mode_s) - 2.0 * alpha / (gn + disc))
        : (gn / 2.0 - std::exp(t_l_d));
    if (centred ? !(v > 0.0) : !(t_l_d < y_0)) {
      ::envelope_failure(env, "region D envelope 15 breakpoint");
    }
    slope_l_d = dlog_g(t_l_d);
    log_g_t_l = log_g(t_l_d);
    if (slope_l_d <= alpha) {
      // Gao & Wang (2025) Lemma 4.5 requires L'(t_l) > alpha for the
      // dual-point construction.
      ::envelope_failure(env, "region D envelope 15 slope");
    }
    p_l_d = t_l_d + (log_g_mode - log_g_t_l) / slope_l_d;
    if (centred) {
      rho = v * (gn - v);
    } else {
      const double t_dot_l = std::log(v);
      rho = log_g_raw(t_dot_l) - alpha * t_dot_l;
    }
  }
  env.rho = rho;
  if (rho <= 0.0) ::envelope_failure(env, "region D secant count");

  // K = ceil(rho) is Gao & Wang's number of breakpoints y_1 < ... < y_K,
  // between which lie K - 1 secants.  For envelope (15) rho is about
  // sqrt(log 4) * gamma_norm, so ceil(rho) passes INT_MAX at gamma_norm
  // ~ 1.8e9 -- well inside the default dispatch, since rmhn's auto path sends
  // every alpha < 1 with gamma > 0 to RTDR.  Converting a double above INT_MAX
  // to int is undefined behaviour, which arm64 and x86-64 resolve differently
  // (INT_MAX and INT_MIN), the second failing the K_max >= 1 test below.
  // Nothing needs K as an int, so the count is kept in double, where ceil is
  // exact.
  const double K = std::ceil(rho);

  // K_max: the largest k <= K with gn^2 > 4 k (rho/K), i.e. the last index at
  // which the breakpoint formula below has a real root.  Theoretically both
  // envelopes always satisfy the guard at every k -- gn^2 - 4 rho is
  // (gn - 2v)^2 for envelope (15) and gn^2/4 for envelope (14), both positive,
  // so the condition is already false at k = K and the loop exits at once --
  // but we keep the loop for defensive programming.  The number of
  // breakpoints actually built, env.K_eff, is at most K_max.
  double K_max = K;
  while (K_max >= 1.0 && gn * gn <= 4.0 * K_max * (rho / K)) {
    K_max -= 1.0;
  }

  if (K_max < 1.0) {
    // Unreachable in exact arithmetic: gn^2 - 4 rho is gn^2/4 for envelope (14)
    // and (gn - 2v)^2 >= gn^2/4 for envelope (15) (see the guard above).  No
    // envelope of this region can be built without a breakpoint, and the hat
    // of region BC need not dominate a region-D density (the minimum over y
    // of (L')^2/2 - L'' is negative here and tends to (2 alpha - 1)/2 as the
    // tilt grows), so the call stops rather than draw from it.
    ::envelope_failure(env, "region D breakpoint count");
  }

  // Breakpoints y_k built incrementally with early termination.
  //
  // Gao & Wang (2025, Section 4, the definition of K following Eq. (15))
  // take K = ceil(rho) breakpoints, so K - 1 secants and the exponential
  // piece left of y_1; their Theorem 4.4 validates that envelope.  For extreme (alpha << 1, gamma >> 1) corners K can be tens
  // of thousands -- e.g. (alpha=0.3, gamma=10000) gives K ~ 12000 -- and
  // constructing every piece dominates the setup cost.
  //
  // When the left tangent at t_l is steeper than alpha it also envelops
  // log g from its left intersection with log g up to y* (Gao & Wang
  // 2025, Section 4, paragraph before Lemma 4.5); their Appendix B code
  // therefore stops the secant sequence at the first y_k at which log g
  // lies on or below that tangent: the tangent then dominates log g from
  // y_k to y*, so the hat stays valid.  The secants dropped lie below the
  // tangent, so stopping early can only enlarge the hat; Theorem 4.4 is
  // proved for the full chain of ceil(rho) breakpoints, and for this one
  // the bound is checked empirically.  We mirror that condition here for
  // envelope (15) of Gao & Wang (2025) -- the dual-point + left-tangent
  // construction selected by the `use_env_15` branch above.
  // Envelope (14) has no left tangent (the leftmost piece is the
  // alpha-slope exponential tail) and K stays small in practice, so it
  // falls through with no early termination.
  //
  // The break decision is made AFTER the current k is appended, so K_used,
  // the number of breakpoints built, is at least 1 and at most K_max <=
  // ceil(rho).
  const double drho = rho / K;
  std::vector<double> y_b;
  std::vector<double> log_g_b;
  y_b.reserve(64);
  log_g_b.reserve(64);
  y_b.push_back(0.0);          // 1-indexed: y_b[0] unused
  log_g_b.push_back(0.0);
  // K_max is a double because K can exceed INT_MAX, so the counter is a long
  // long and is capped.  Envelope (14) is only selected when the mode-to-
  // inflection drop is at most 2.49, which forces gamma_norm <= 6.31 and hence
  // ceil(rho) = ceil(3 gamma_norm^2/16) <= 8; envelope (15) stops at the
  // tangent test below within a handful of steps at every parameter measured.
  // If neither bound fires, truncating the chain would leave the left tangent
  // spanning an interval on which the tangent has not been shown to dominate
  // log g, i.e. a hat that is not one, so fail loudly rather than build it.
  const long long k_cap = 1000000LL;
  long long K_actual = 0;
  for (long long k = 1; static_cast<double>(k) <= K_max && k <= k_cap; ++k) {
    const double rk = 4.0 * static_cast<double>(k) * drho;
    const double inside = gn * gn - rk;
    if (inside <= 0.0) {
      ::envelope_failure(env, "region D breakpoint");
    }
    // Placed on the log axis, then moved to s.  Every breakpoint lies at or
    // left of the dual point, where e^y <= v < 1.5, while e^{m_g} is about
    // gamma_norm/2 > 1e6 on the centred path, so log g there is about
    // -gamma_norm^2/4 below the mode: those pieces carry no mass in double,
    // and the log-axis resolution of their position is ample.
    const double yk = std::log(2.0 * static_cast<double>(k) * drho
                               / (std::sqrt(inside) + gn)) - origin;
    const double log_g_yk = log_g(yk);
    y_b.push_back(yk);
    log_g_b.push_back(log_g_yk);
    K_actual = k;
    if (use_env_15) {
      // Gao & Wang (2025) Appendix B: break once log g(y_k) is on or
      // below the tangent line at t_l_d.  Beyond this point further
      // secants are not needed for a valid hat (they would only lower it).
      const double tangent = log_g_t_l + slope_l_d * (yk - t_l_d);
      if (log_g_yk <= tangent) break;
    }
  }
  if (K_actual == k_cap && static_cast<double>(k_cap) < K_max) {
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
  // is immune to the coordinate choice; it is gated only so that below the
  // switch every expression is the uncentred one.
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
  // reason they do on the left: from m_g + 1 the search would need 35
  // halvings at gamma_norm = 1e10, and the absolute 1e-10 step tolerance ends
  // it earlier still, returning a contact whose drop is 16.2 at
  // gamma_norm = 1e11 and 1526 at 1e12 against the target log 4.  That widens
  // the plateau by a factor of 33 at 1e12; the draws stay correctly
  // distributed, because a wider hat over a correct accept test is still a
  // correct sampler, so the KS invariant does not see it -- only the
  // acceptance rate does, and it falls to 0.086 against the 1/e that
  // Gao & Wang (2025, Theorem 4.4) prove for the full secant chain.
  const double t_r_init = centred ? mode_s + inc : m_g + 1.0;
  env.t_r = ::contact_point_newton(t_r_init, log_g_mode, delta,
                                   log_g, dlog_g, 0.93, 1.99,
                                   200, newton_tol, mode_s,
                                   ::right_contact_bound(m_g, origin));
  env.slope_r = dlog_g(env.t_r);
  if (!(env.slope_r < 0.0) || !std::isfinite(env.slope_r)) {
    ::envelope_failure(env, "region D right contact");
  }
  env.p_r = env.t_r + (log_g_mode - log_g(env.t_r)) / env.slope_r;

  // Build piece table: the exponential tail, K_used - 1 secants, at most one
  // left tangent, the plateau and the right tail.
  std::vector<EnvelopePiece> pieces;
  pieces.reserve(static_cast<std::size_t>(K_used) + 3);

  // 1. Leftmost exponential tail (-inf, y_1] with slope alpha (the asymptotic
  //    slope of log g(y) as y -> -inf).
  pieces.push_back(make_exp_left(y_b[1], alpha, log_g_b[1]));

  // 2. Secants on (y_{k-1}, y_k], k=2..K_used.
  for (int k = 2; k <= K_used; ++k) {
    pieces.push_back(make_secant(y_b[k-1], y_b[k], a_k[k], log_g_b[k-1]));
  }

  if (use_env_15) {
    // 3a. Left-tangent piece (y_K, p_l] for envelope (15), where
    //     y_K = y_b[K_used] is the last breakpoint kept: y_K <= t_dot_l, with
    //     equality only when the chain runs to its full length ceil(rho),
    //     and strictly below whenever the early-termination test above
    //     stopped it sooner -- at k = 1 for every large tilt, where y_K is
    //     about -log(gamma_norm), some 2 log(gamma_norm) left of the mode.
    //     h(y) = log_g(t_l) + L'(t_l) * (y - t_l) on this interval.
    //     Gao & Wang (2025) print this interval as (y*, p_l) in Eq. (15),
    //     which leaves (t_dot_l, y*] uncovered because their secants end at
    //     t_dot_l < y*.  Their Theorem 4.3 has (\dot t_l, p_l], and their
    //     Appendix B code, like this one, starts the piece at the last
    //     breakpoint kept.
    // The plateau's left edge must lie to the right of the last breakpoint, or
    // the piece is empty and its area formula below is a log of a
    // non-positive number.
    if (!(y_b[K_used] < p_l_d)) {
      ::envelope_failure(env, "region D tangent interval");
    }
    if (centred) {
      // Past the switch the piece is written as Gao & Wang (2025, Eq. (15))
      // do, g(m) exp{L'(t_l) (y - p_l)}: anchored at p_l, where by construction
      // the tangent reaches the mode's height log_g_mode.  Its mass lies within
      // a few 1/L'(t_l) of p_l, and anchored at y_K instead both the base and
      // L'(t_l) (y - y_K) would be of size L'(t_l) (p_l - y_K) -- 6e13 at
      // gamma_norm = 1e12 -- cancelling in the hat's value and in its inverse
      // CDF to leave errors of order 1e-2 in log h, enough for the hat to dip
      // below the target (see PIECE_EXP_BOUNDED_RIGHT).
      pieces.push_back(make_exp_bounded_right(y_b[K_used], p_l_d, slope_l_d,
                                              log_g_mode));
    } else {
      // Anchored at piece.a = y_K: base = log_g(t_l) + L'(t_l)*(y_K - t_l).
      const double tangent_base =
          log_g_t_l + slope_l_d * (y_b[K_used] - t_l_d);
      pieces.push_back(
          make_exp_bounded(y_b[K_used], p_l_d, slope_l_d, tangent_base));
    }
    // 4a. Plateau (p_l, p_r] at height log g(m_g).  p_l < m_g < p_r holds
    // mathematically, so an empty plateau means the construction failed;
    // dropping it silently would leave a support that does not contain the
    // mode (sd/true = 0.0146 at gamma_norm = 1e10 on the raw expressions, see
    // above).  Written !(p_l < p_r) so a non-finite endpoint is caught too, as
    // in region BC's guard.
    if (!(p_l_d < env.p_r)) {
      ::envelope_failure(env, "region D intersection");
    }
    pieces.push_back(make_plateau(p_l_d, env.p_r, log_g_mode));
    env.t_l = t_l_d;
    env.p_l = p_l_d;
    env.slope_l = slope_l_d;
    env.has_left_tangent_d = true;
  } else {
    // 3b. Plateau (y_K = y_star, p_r] for envelope (14) (the paper's
    //     [y*, p_r], half open so that the last secant keeps y_K).
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
// influencing acceptance at all.  On the log axis the abscissa is centred as
// well, so the argument is s = y - m_g itself.
double log_target_y_centred(double s, double alpha, double u_mode) {
  const double E = std::expm1(s);
  return -(u_mode * u_mode * E * E + alpha * (E - s));
}

// Region A's target measured from the mode, for the same reason, at the
// offset s = x - m from its mode m (see setup_region_a).
double log_target_x_centred(double s, double alpha, double m) {
  if (s <= -m) return -std::numeric_limits<double>::infinity();
  return (alpha - 1.0) * mhn::log1p_minus_u(s / m) - s * s;
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

bool rtdr_region_a_centred(double alpha, double gamma_norm) {
  const double m = mhn::positive_root(1.0, gamma_norm, alpha - 1.0);
  return ::region_a_centred_at(alpha, gamma_norm, m);
}

RtdrEnvelope build_rtdr_envelope(double alpha, double beta, double gamma) {
  RtdrEnvelope env;
  env.sqrt_beta = std::sqrt(beta);
  env.alpha = alpha;
  env.beta = beta;
  env.gamma = gamma;
  env.gamma_norm = gamma / env.sqrt_beta;
  env.region = ::classify_region(alpha, env.gamma_norm);
  switch (env.region) {
    case REGION_A:  ::setup_region_a(env);  break;
    case REGION_BC: ::setup_region_bc(env); break;
    case REGION_D:  ::setup_region_d(env);  break;
  }
  return env;
}

double sample_rtdr(const RtdrEnvelope& env, int* retries_out, int max_retries) {
  // Region A: sampling and accept/reject in x-space (x > 0).
  // Region BC/D: sampling and accept/reject in y-space (y in R), then
  //   x = exp(y) on output.
  // Scale restoration: X = X' / sqrt(beta).  On the centred path the offset s
  // from the mode is placed at the mode in the caller's coordinate instead
  // (see RtdrEnvelope::origin), so that the draw rounds once, to the double
  // nearest the exact draw.
  // The retry budget, 1000 proposals unless the caller sets another, is a
  // safety net: the acceptance probability is at least 1/e (Gao & Wang 2025)
  // when the contact points lie in the band of their Eq. (7), which the
  // search does not guarantee (see contact_point_newton), and, in region D,
  // when the secant chain is the full one, which the package truncates (see
  // setup_region_d).
  const bool y_space = (env.region != REGION_A);
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    const int idx = ::select_piece(env);
    const EnvelopePiece& piece = env.pieces[idx];
    double s = ::sample_within_piece(piece);
    // A draw from the leftmost piece of regions BC and D,
    // h(y) = h(a) exp(slope (y - a)) on the log axis (region D's exponential
    // tail in Gao & Wang 2025, Eqs. (14)-(15)), is a + log(u) / slope, which
    // passes the largest double, and is -Inf, where |log u| exceeds
    // slope DBL_MAX: for a shape near 2^-1022 (slope alpha) at u below e^-4.
    // There log g(y) - log h(y) = (alpha - slope) y + phi(e^y) + const, with
    // phi(t) = gamma_norm t - t^2, and phi(e^y) is 0 to double precision
    // below y = -745.  So the draws past y0 = -DBL_MAX are accepted with
    // probability (slope / alpha) g(y0) / h(y0) on average: their density is
    // proportional to exp(slope y), and that of the ratio to
    // exp((alpha - slope) y).  The draw is placed at y0 for the test, with
    // that factor, and is 0 at every beta: e^y / sqrt(beta) is below the
    // smallest positive double for y <= y0, sqrt(beta) being at least
    // 2.2e-162.  The factor is 1 wherever this can happen: region D's tail
    // has slope alpha by construction, and region BC's left contact lies so
    // far left at such a shape (between -sqrt(2 / alpha) and
    // -1 / sqrt(alpha)) that its slope is alpha to the bit; it is kept so that
    // the test stays exact for any tail of a smaller slope.  Taken as -Inf,
    // both sides of the test would be -Inf and every such draw accepted, so
    // that 0 would come out too often.
    double log_factor = 0.0;
    const bool past_doubles =
        (piece.type == PIECE_EXP_LEFT &&
         s == -std::numeric_limits<double>::infinity());
    if (past_doubles) {
      s = -std::numeric_limits<double>::max();
      log_factor = std::min(0.0, std::log(piece.slope / env.alpha));
    }
    const double log_h = ::log_piece_at(piece, s);
    // The target must be measured in whatever coordinate the envelope was
    // built in.  If they disagree the two sides of the test are again of size
    // gamma_norm^2/4, log_u + log_h == log_h exactly, and the uniform draw
    // stops influencing acceptance at all -- the sampler then draws from its
    // own hat.
    const double log_target =
        env.centred
          ? (y_space ? ::log_target_y_centred(s, env.alpha, env.u_mode)
                     : ::log_target_x_centred(s, env.alpha, env.origin))
          : (y_space ? ::log_target_y(s, env.alpha, env.gamma_norm)
                     : ::log_target_x(s, env.alpha, env.gamma_norm));
    const double log_u = std::log(R::runif(0.0, 1.0));
    if (log_u + log_h <= log_target + log_factor) {
      if (retries_out != nullptr) *retries_out += retries;
      if (past_doubles) return 0.0;
      if (!env.centred) {
        if (!y_space) return s / env.sqrt_beta;
        // e^s is the draw of sqrt(beta) X.  Under a strong negative tilt at a
        // small shape it can lie below the smallest normal double, where it
        // keeps only some of its digits or none, while X itself, larger by
        // 1 / sqrt(beta) for beta < 1, need not: at (1e-3, 1e-300, -100)
        // seven in ten draws of sqrt(beta) X fall below it.  There X is
        // formed as exp(s - log sqrt(beta)).
        const double u = std::exp(s);
        return (u >= std::numeric_limits<double>::min() || env.sqrt_beta >= 1.0)
            ? u / env.sqrt_beta
            : std::exp(s - std::log(env.sqrt_beta));
      }
      return y_space
          ? env.x_mode_hi + (env.x_mode_lo + env.x_mode_hi * std::expm1(s))
          : env.x_mode_hi + (env.x_mode_lo + s / env.sqrt_beta);
    }
    ++retries;
  }
  // Retry budget exhausted.  At the acceptance probability of 1/e that
  // Gao & Wang (2025) guarantee for contact points in the band of their
  // Eq. (7) (in region D, for the full secant chain, which the package
  // truncates; see setup_region_d), 1000 proposals all fail with probability
  // below 1e-199; outside the band the guarantee is lost (see
  // contact_point_newton).  Report it as NaN and leave the reporting to the
  // caller: raising an R warning here would leave this frame by a long jump
  // under options(warn = 2), skipping the destructors of the envelope's
  // vectors.
  if (retries_out != nullptr) *retries_out += retries;
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace mhn

// ====================================================================
// Entry points used by the test suite; not part of the public API.
// ====================================================================

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

// The envelope for inspection.  Positions (mode, t_l, t_r, p_l, p_r, y_star,
// y_break) are reported on the region's own axis, sqrt(beta) x in region A and
// y = log(sqrt(beta) x) in regions BC and D, by adding `origin` back; the
// piece table is reported as stored, in the sampling coordinate
// s = (x or y) - origin, which keeps the resolution those numbers lose once
// origin is added.
// [[Rcpp::export(.dump_rtdr_envelope_cpp, rng = false)]]
Rcpp::List dump_rtdr_envelope_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::RtdrEnvelope env = mhn::build_rtdr_envelope(alpha, beta, gamma);
  const double o = env.origin;
  Rcpp::NumericVector y_break(env.y_break.begin(), env.y_break.end());
  y_break = y_break + o;
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
    Rcpp::Named("mode")            = env.mode + o,
    Rcpp::Named("log_dens_mode")   = env.log_dens_mode,
    // Exposed so that a test rebuilding log f or log g for a domination check
    // knows which ordinate the piece table is measured in: past the switch
    // every base_log_dens and log_area is relative to the peak, and comparing
    // them against a raw log-density silently compares the wrong two things.
    Rcpp::Named("centred")         = env.centred,
    Rcpp::Named("origin")          = o,
    Rcpp::Named("u_mode")          = env.u_mode,
    Rcpp::Named("x_mode_hi")       = env.x_mode_hi,
    Rcpp::Named("x_mode_lo")       = env.x_mode_lo,
    Rcpp::Named("t_l")             = env.t_l + o,
    Rcpp::Named("t_r")             = env.t_r + o,
    Rcpp::Named("p_l")             = env.p_l + o,
    Rcpp::Named("p_r")             = env.p_r + o,
    Rcpp::Named("slope_l")         = env.slope_l,
    Rcpp::Named("slope_r")         = env.slope_r,
    Rcpp::Named("simplified")          = env.simplified,
    // Exposed so that a test can hold the default method of rmhn to the
    // setup's own choice of start (see rtdr_region_a_centred in mhn_rtdr.h).
    Rcpp::Named("quick_start")         = env.quick_start,
    Rcpp::Named("has_left_tangent_d")  = env.has_left_tangent_d,
    Rcpp::Named("K_eff")               = env.K_eff,
    Rcpp::Named("y_star")          = env.y_star + o,
    Rcpp::Named("rho")             = env.rho,
    Rcpp::Named("y_break")         = y_break,
    Rcpp::Named("log_dens_break")  = env.log_dens_break,
    Rcpp::Named("alpha_k")         = env.alpha_k,
    Rcpp::Named("piece_log_area")  = env.piece_log_area,
    Rcpp::Named("pieces")          = pieces
  );
}

// The concavity tag selects one of the three fixed test fixtures defined above.
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
