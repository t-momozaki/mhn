// Quantile function (qmhn) -- vectorized over p and (alpha, beta, gamma).
//
// Dispatch:
//   * gamma  = 0   -> sqrt-Gamma special case (Sun et al. 2023, Lemma 6a):
//                     x = sqrt(y / beta), where P(alpha/2, y) = p and P is
//                     the regularised lower incomplete gamma function.
//   * alpha  = 1   -> truncated normal inverse (Sun et al. 2023, Lemma 6b).
//                     Both are taken at those values exactly (see
//                     mhn_special_cases.h, is_sqrt_gamma_exact).
//   * otherwise    -> Boost TOMS 748 root-finder (bisection where it cannot
//                     be trusted or does not converge) on F(x) - p_target = 0,
//                     using CdfState from mhn_pmhn_state.h to evaluate F, or
//                     for lower.tail = FALSE on
//                     log p_target - log P(X > x) = 0 (unless the survival
//                     probability is so close to 1 that the logarithm of
//                     the upper tail no longer resolves 1 - p, see
//                     upper_route_limit, where F(x) = 1 - p is solved, or
//                     the survival root polished on log F, see
//                     polishes_on_lower_tail); on the upper tail for a
//                     lower-tail p within lower_route_limit of 1 (at
//                     least 1e-7, set from the error bounds of the two
//                     tails); and for a lower-tail
//                     log p below log(DBL_MIN) on log F(x) - log p = 0.
//
// The search bracket starts at [s sqrt(eps), s] with s = E(X) + 8 sqrt(Var(X)).
// Its lower end is halved (down to the denormal floor) while F there is still
// above the target, and its upper end doubled (up to the largest double)
// while F there is below it.  Where TOMS 748 cannot be trusted with the
// bracket, or does not converge on it, the bracket is bisected instead.

#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_pmhn_state.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_param_slots.h"
#include "mhn_warn.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <algorithm>
#include <cmath>
#include <exception>
#include <cstdint>
#include <limits>

// Defined in mhn_moments.cpp, at file scope beside the Rcpp exports of the
// mean and the variance, whose doubles it returns; the quantile bracket wants
// those guarded forms, not a second copy of the algebra.
void mhn_mean_var_cpp(double alpha, double beta, double gamma,
                      double log_psi_alpha, double* mean, double* var);

namespace {

// What the root find needs of a parameter triple besides its distribution
// function, cached so that all elements sharing the triple reuse it.
struct QuantileMetadata {
  bool valid = false;
  // The law's standard deviation.  It is kept as it is and not as a variance:
  // for beta below about 2.8e-309 its square passes the largest double.
  double sd = 0.0;
  double initial_x_hi = 1.0;
  // The root is wanted to a relative 2^-46 of x, and in the peak frame
  // (mhn_cdf_peak.h), where the law can be narrower than that, also to within
  // this absolute width, 2^-46 of its standard deviation; +Inf elsewhere.
  double abs_tol = std::numeric_limits<double>::infinity();
};

// The metadata is derived from a CdfState and is only valid for the triple
// that state holds, so the two are cached together.
struct QuantileState {
  mhn::CdfState state;
  QuantileMetadata meta;
  void recompute(double a, double b, double g) {
    state.recompute(a, b, g);
    meta.valid = false;  // recomputed on demand, general case only
  }
};

// The mean and variance of the law for the triple the state holds, and the
// scale s = E(X) + 8 sd(X) from which the search bracket is set (see
// solve_quantile_general).
void recompute_metadata(const mhn::CdfState& state, QuantileMetadata& meta) {
  const double alpha = state.alpha;
  const double beta = state.beta;
  const double gamma = state.gamma;
  meta.abs_tol = std::numeric_limits<double>::infinity();
  if (state.peak_frame) {
    // The law is concentrated about the peak v* of the density of log X,
    // with standard deviation 1 / sqrt(2 beta + (alpha - 1) / v*^2) to a
    // relative 1 / sqrt(alpha) or 1 / z; that is scale enough for the
    // bracket, and it is also the width the root must be found within.  For
    // alpha >= 1 it is formed as v* / hypot(sqrt(alpha - 1), sqrt(2 beta) v*),
    // which neither overflows nor underflows where v* is tiny or huge; where
    // v* lies beyond the largest double it is 1 / sqrt(2 beta).
    const double v = state.peak_ok ? state.peak.v : R_NaN;
    double sd;
    if (!(mhn::opaque_value(v) <= std::numeric_limits<double>::max())) {
      sd = M_SQRT1_2 / std::sqrt(beta);
    } else if (alpha >= 1.0) {
      sd = v / std::hypot(std::sqrt(alpha - 1.0),
                          M_SQRT2 * std::sqrt(beta) * v);
    } else {
      sd = 1.0 / std::sqrt(2.0 * beta + ((alpha - 1.0) / v) / v);
    }
    meta.sd = sd;
    double x_hi = v + 8.0 * sd;
    if (!(mhn::opaque_value(x_hi) <= std::numeric_limits<double>::max())) {
      x_hi = std::numeric_limits<double>::max();
    }
    if (!(x_hi > 0.0)) x_hi = 1.0;
    meta.initial_x_hi = x_hi;
    meta.abs_tol = std::ldexp(sd, -46);
    if (!(meta.abs_tol > 0.0)) {
      meta.abs_tol = std::numeric_limits<double>::infinity();
    }
    meta.valid = true;
    return;
  }
  // Both the mean and the variance are taken from the guarded helpers of
  // mhn_moments.cpp rather than rebuilt here.  Forming
  // exp(log Psi(alpha+1) - log Psi(alpha)) directly subtracts two quantities
  // that are each about gamma^2/(4 beta) for a large positive tilt -- 2.5e9
  // at gamma = 1e5 -- and the mean switches to a quadrature past the point
  // where that costs more than the working tolerance, while the variance
  // switches on the measured cancellation in the Lemma 2c form.  This is only
  // the bracket for the root find, but a bracket built from a mean that is 31
  // percent wrong is a slow start at best.  The helper forms the Lemma 2a
  // ratio once for both, and takes its denominator log Psi[alpha/2, z] from
  // the state, which holds it for the general case below z = 1e3 (NaN
  // otherwise, when it is evaluated there if needed).
  double mu, v;
  mhn_mean_var_cpp(alpha, beta, gamma, state.log_psi, &mu, &v);
  meta.sd = std::sqrt(v);
  double x_hi = mu + 8.0 * std::sqrt(v);
  if (!(x_hi > 0.0) || !std::isfinite(x_hi)) {
    x_hi = std::max(mu, 1.0);
    if (!(x_hi > 0.0)) x_hi = 1.0;
  }
  meta.initial_x_hi = x_hi;
  meta.valid = true;
}

// Within this distance of 1 an upper-tail probability p is always solved for
// on the lower tail, F(x) = 1 - p (directly, or by polishing a root of the
// survival equation on it, see polish_on_lower_tail), and a lower-tail one
// on the upper tail (see upper_log_prob).  Near the lower end of the law
// log P(X > x) is about -F(x), and it carries an absolute error of at least
// a few units of eps (see kUpperLogErrorScale in mhn_pmhn_state.h), which
// is all of F once F is near 1e-16; near the upper end F is resolved only
// to the rounding of 1 in the same way.
const double kSurvivalSolveLimit = 1e-7;

// The survival root that starts the polish is found only to a relative
// 2^(1 - kPolishStartDigits) = 1.2e-7 of x.  Where the polish is used the
// survival equation resolves d = 1 - p only to upper_log_error / d of
// itself, a bound of 1e-9 or more there, and x to that divided by
// x f(x) / F(x); near its root it is flat to within that noise, and a root
// finder asked for 2^-46 of x spends most of its evaluations of the upper
// tail there.  The Newton steps on log F take the root on to full accuracy
// either way: from 1.2e-7 one step leaves a relative error of C times its
// square (see polish_on_lower_tail), and the next confirms it.  Measured
// per value, this halves the cost of the polish at small tilts.
const int kPolishStartDigits = 24;

// Within kSurvivalSolveLimit of 1, the survival equation whose root starts
// the polish is solved for 1 - p no smaller than this many times
// upper_log_error.  The log upper tail resolves such a 1 - p to within 1/16
// of itself, so the root lies where F is within that of 1 - p; for a
// smaller 1 - p it could lie anywhere F is of the order of that error, and
// the search for it could run down to the smallest double.
const double kSurvivalStartErrors = 16.0;

// An upper-tail probability p = 1 - d is solved for on log P(X > x) = log p
// only while d is at least this large.  The logarithm of the upper tail
// carries an absolute error of up to state.upper_log_error, a relative error
// of that divided by d in d, and the quantile moves by that relative error
// times F / (x f(x)) of itself and F / (f(x) sd) of the law's standard
// deviation.  Solved on it, the quantile of d = 1e-6 would be 9e-5 of a
// standard deviation off at (9e5, 1, -500), where those terms reach 7e6,
// and 3e-6 off at (0.5, 1, 900), where log Psi is 2e5.  Below the limit that
// relative error could pass kUpperLogRelTol, and the lower-tail equation
// F(x) = d is used instead, solved directly or by polishing the survival
// root on it (see polishes_on_lower_tail), with d formed exactly: as
// -expm1(log p) on the log scale, and on the linear scale by a subtraction
// that Sterbenz's lemma makes exact for p >= 1/2.  The limit is at most 1/2
// (d below 1/2 keeps the subtraction exact), and at least
// kSurvivalSolveLimit.  For the special cases and the peak frame, whose log
// upper tail keeps its relative accuracy near 0, it is kSurvivalSolveLimit.
double upper_route_limit(const mhn::CdfState& state) {
  const double limit = state.upper_log_error / mhn::kUpperLogRelTol;
  if (!(limit <= 0.5)) return 0.5;
  return (limit > kSurvivalSolveLimit) ? limit : kSurvivalSolveLimit;
}

// A lower-tail p = 1 - d is solved for on the upper tail,
// log P(X > x) = log d, while d is below this limit, and on F(x) = p above
// it: the counterpart of upper_route_limit, asking of the equation what
// log_cdf_general asks of log F (see there).  Below the peak of the density
// of log X, F carries an absolute error of up to about lower_log_error, a
// relative error of that divided by d in d, which for a small shape, where
// F is close to 1 over most of that range, can pass kUpperLogRelTol long
// before d reaches kSurvivalSolveLimit: the quantile of a lower-tail
// log p = log1p(-1e-6) at (1e-8, 1, 0.5) was 1.3e-7 of itself off.  The
// upper tail resolves d to a relative error of up to upper_log_error.  So d
// is taken on the upper tail where lower_log_error / d could pass
// kUpperLogRelTol and 2 upper_log_error d / (1 - d) < lower_log_error, as
// in log_cdf_general, the factor 2 leaving the lower tail where the two
// bounds, each about twice the errors measured, do not tell which is the
// better: d below lower_log_error / kUpperLogRelTol and below r / (1 + r),
// r = lower_log_error / (2 upper_log_error).  The limit is at most 1/2,
// where 1 - p is exact for a p on the linear scale (Sterbenz's lemma), and
// at least kSurvivalSolveLimit.  The special cases and the peak frame,
// whose lower_log_error is 0, keep kSurvivalSolveLimit.
double lower_route_limit(const mhn::CdfState& state) {
  const double lower = state.lower_log_error;
  const double upper = state.upper_log_error;
  if (!(lower > 0.0) || !(upper > 0.0)) return kSurvivalSolveLimit;
  const double r = lower / (2.0 * upper);
  double limit = std::min(lower / mhn::kUpperLogRelTol, r / (1.0 + r));
  if (!(limit <= 0.5)) limit = 0.5;
  return (limit > kSurvivalSolveLimit) ? limit : kSurvivalSolveLimit;
}

// Whether an upper-tail p within upper_route_limit of 1 is solved on the
// upper tail and then polished on the lower one (see polish_on_lower_tail)
// rather than solved on the lower tail.  Where the
// lower tail comes from the Lemma 1b series it costs more than the
// upper-tail quadrature: up to about 100 microseconds for |z| <= 20
// (see kCdfSeriesMaxAbsZ in mhn_cdf_series.h) and, for alpha < 1 past it,
// milliseconds (25 ms at (0.01, 1, 500), 9 ms at (0.3, 1, 300)), against
// about 1 to 5 microseconds.  A lower-tail solve spends some 30 evaluations
// of it where the polish spends a few, and past |z| = 20 the polish
// takes the lower-tail quadrature instead wherever it agrees with the upper
// tail (CdfState::log_lower_resolving).  Where the lower tail comes from a
// quadrature it costs about as much as the upper tail, and the lower-tail
// solve took from about 0.2 to 1 times as long as the survival solve at the
// shapes and tilts measured: 0.2 to 0.4 at (2, 1, -30), 0.4 to 0.8 at
// alpha = 412 and about 1 at alpha = 5e4.
bool polishes_on_lower_tail(const mhn::CdfState& state) {
  return state.kind == mhn::CDF_GENERAL && !state.peak_frame &&
         !state.shifted && state.series_first;
}

// True when p, of either tail, is within limit of 1: on the log scale when
// log p > -limit, which differs from 1 - p < limit only by a relative
// limit / 2 in the limit.  p is not NA or NaN.
bool near_one(double p, bool log_p, double limit) {
  return log_p ? (p > -limit) : (p > 1.0 - limit);
}

// Convert a user-provided (p, lower_tail, log_p) to a lower-tail probability
// p_low in [0, 1].  Returns NA_REAL for NA inputs and NaN for out-of-range.
double to_lower_tail_prob(double p, bool lower_tail, bool log_p) {
  if (R_IsNA(p)) return NA_REAL;
  if (R_IsNaN(p)) return R_NaN;
  double p_low;
  if (log_p) {
    if (p > 0.0) return R_NaN;
    if (p == R_NegInf) p_low = 0.0;
    else               p_low = std::exp(p);
  } else {
    if (p < 0.0 || p > 1.0) return R_NaN;
    p_low = p;
  }
  // For an upper-tail log p near 0, 1 - exp(p) cancels and -expm1(p) does
  // not; this value is used for such a p (see upper_route_limit).
  if (!lower_tail) p_low = log_p ? -std::expm1(p) : 1.0 - p_low;
  if (p_low < 0.0) p_low = 0.0;
  if (p_low > 1.0) p_low = 1.0;
  return p_low;
}

// The upper-tail probability on the log scale, or NaN when the lower-tail
// equation is to be solved instead: for a lower-tail p not within
// lower_limit, the state's lower_route_limit, of 1, and for an upper-tail
// one within limit, the state's upper_route_limit, of 1.  Keeping it in
// logs is the whole point: an upper-tail p of 1e-30 is perfectly
// representable, and only the conversion to 1 - p destroys it.
//
// A lower-tail p within lower_limit of 1 is solved for on the upper tail as
// well, at least within kSurvivalSolveLimit.  Near the upper end of the
// law F is resolved only to the
// rounding of 1, eps, and the root of F(x) = p moves by eps / f(x), all of
// the distance to the true quantile once 1 - p is of order eps; the upper
// tail is formed there directly, to its own relative accuracy, and 1 - p is
// exact: by Sterbenz's lemma on the linear scale, where p >= 1/2, and as
// -expm1(p) to a unit in its last place on the log scale.
double upper_log_prob(double p, bool lower_tail, bool log_p, double limit,
                      double lower_limit) {
  if (R_IsNA(p) || R_IsNaN(p)) return R_NaN;
  if (lower_tail) {
    if (log_p) {
      if (p > 0.0 || !near_one(p, true, lower_limit)) return R_NaN;
      return std::log(-std::expm1(p));
    }
    if (p < 0.0 || p > 1.0 || !near_one(p, false, lower_limit)) {
      return R_NaN;
    }
    return std::log(1.0 - p);
  }
  if (log_p) {
    if (p > 0.0 || near_one(p, true, limit)) return R_NaN;
    return p;
  }
  if (p < 0.0 || p > 1.0 || near_one(p, false, limit)) return R_NaN;
  return std::log(p);
}

// A lower-tail log probability too small for exp() of it to hold as a normal
// double, or NaN otherwise.  Below log(DBL_MIN), about -708.4, exp(p) keeps
// fewer significant bits the smaller it is, and below about -745.1 it is 0,
// while the quantile itself can still be a normal double: at (2, 1, 3),
// F(x) is about 0.0198 x^2 near the origin, so log p = -1400 has its quantile
// near 7e-304.  Such a p is solved for on log F instead.  For an upper-tail
// log p within limit of 0 the lower-tail log probability is
// log(-expm1(p)), which is below log(DBL_MIN) for |p| below DBL_MIN.
double lower_log_prob(double p, bool lower_tail, bool log_p, double limit) {
  static const double log_dbl_min =
      std::log(std::numeric_limits<double>::min());
  if (!log_p || R_IsNA(p) || R_IsNaN(p)) return R_NaN;
  double lp = p;
  if (!lower_tail) {
    if (p > 0.0 || !near_one(p, true, limit)) return R_NaN;
    lp = std::log(-std::expm1(p));
  }
  if (!(lp < log_dbl_min) || lp == R_NegInf) return R_NaN;
  return lp;
}

// True when no double lies strictly between a < b, so that a bracket [a, b]
// cannot be narrowed any further.
bool adjacent_doubles(double a, double b) {
  return !(a < b) || std::nextafter(a, b) == b;
}

// The bracket [x_lo, x_hi] of the root of an increasing f, with
// f(x_lo) <= 0 <= f(x_hi), found from a starting pair by halving the lower
// end and doubling the upper one; or, where there is no such bracket to hand
// a root finder, the value the quantile takes instead (done).
struct QuantileBracket {
  double x_lo, x_hi, f_lo, f_hi;
  bool done;
  double value;
};

template <class F>
QuantileBracket bracket_quantile(const F& f_at, double x_lo, double x_hi,
                                 const mhn::CdfState& state) {
  double f_lo = f_at(x_lo);
  double f_hi = f_at(x_hi);

  // If the lower end is not below the target, halve it until it is.  For
  // alpha < 1 the distribution function behaves like x^alpha near the origin,
  // so the quantile of a small p sits at about p^(1/alpha) in units of the
  // scale: of order 1e-40 at alpha = 0.1 and p = 1e-4.  Each abandoned end
  // still has f > 0, so it becomes the upper end, and the root finder starts
  // from a bracket a factor of two wide rather than one spanning every decade
  // the search went down.  The loop may run from any double down to the
  // denormal floor, about 2100 halvings; each is one CDF evaluation, and it
  // stops as soon as the bracket straddles the root.
  int shrink = 0;
  while (f_lo > 0.0 && x_lo > std::numeric_limits<double>::denorm_min() &&
         shrink < 2100) {
    const double x_try = x_lo * 0.5;
    double f_try;
    try {
      f_try = f_at(x_try);
    } catch (const std::exception&) {
      // For alpha < 1 the kernel x^(alpha-1) grows without bound towards the
      // origin, and far enough down the quadrature can no longer represent it.
      // That point is below the smallest quantile this arithmetic can express,
      // so stop here and report the endpoint reached.
      break;
    }
    if (!std::isfinite(f_try)) break;
    x_hi = x_lo;
    f_hi = f_lo;
    x_lo = x_try;
    f_lo = f_try;
    ++shrink;
  }

  // Likewise double the upper end until it is above the target, moving the
  // lower end up behind it, up to the largest double: at most about 2100
  // doublings from the smallest positive double, one CDF evaluation each,
  // and a handful unless the target lies far out in the upper tail.  A root
  // beyond the largest double is a quantile beyond the range of doubles,
  // returned as +Inf, as qgamma returns a quantile that overflows.
  const double x_max = std::numeric_limits<double>::max();
  int expand = 0;
  while (f_hi < 0.0 && expand < 2200) {
    if (x_hi >= x_max) {
      return QuantileBracket{x_lo, x_hi, f_lo, f_hi, true, R_PosInf};
    }
    x_lo = x_hi;
    f_lo = f_hi;
    x_hi = (x_hi <= 0.5 * x_max) ? 2.0 * x_hi : x_max;
    f_hi = f_at(x_hi);
    ++expand;
  }
  // f_hi may be +Inf, where the survival probability is 0 in double
  // precision or its logarithm below the double range: a legitimate upper end.
  if (f_hi < 0.0 || ISNAN(f_hi) || ISNAN(f_lo)) {
    // There is no root to report.  NA, as for any value that cannot be
    // determined, rather than the NaN reserved for an invalid probability.
    mhn::warn("qmhn: could not bracket the quantile for alpha=%g, "
              "beta=%g, gamma=%g; returning NA",
              state.alpha, state.beta, state.gamma);
    return QuantileBracket{x_lo, x_hi, f_lo, f_hi, true, NA_REAL};
  }
  if (f_lo > 0.0) {
    // The root lies below the lowest point the search reached.  Where that is
    // the smallest positive double the quantile is 0 in double precision, and
    // is reported so, as qgamma and the special cases here report a quantile
    // that underflows.  Otherwise the search stopped because the distribution
    // function could no longer be evaluated further down, and the point
    // reached is the best bound available.
    return QuantileBracket{
        x_lo, x_hi, f_lo, f_hi, true,
        (x_lo <= std::numeric_limits<double>::denorm_min()) ? 0.0 : x_lo};
  }
  return QuantileBracket{x_lo, x_hi, f_lo, f_hi, false,
                         std::numeric_limits<double>::quiet_NaN()};
}

// The root of the equation of solve_quantile_general taken in t = log x,
// log F(e^t) = log_target, or the survival form when upper is set, by
// TOMS 748 from the bracket [b_lo, b_hi] of x, to a relative 2^(1 - digits)
// in x (an absolute one in t) or to adjacent doubles in t; NaN where the
// bracket does not hold in t or the root finder does not converge.
double solve_in_log_x(const mhn::CdfState& state, double b_lo, double b_hi,
                      bool upper, double log_p_up, double log_target,
                      int digits, std::uintmax_t iter_limit) {
  using boost::math::tools::toms748_solve;
  auto g_at = [&state, log_target, upper, log_p_up](double t) -> double {
    const double x = std::exp(t);
    double g;
    if (upper) g = log_p_up - state.log_upper_tail(x);
    else g = state.log_cdf_general(x) - log_target;
    if (g > 1e10) return 1e10;
    if (g < -1e10) return -1e10;
    return g;
  };
  const double t_tol = std::ldexp(1.0, 1 - digits);
  auto t_close = [t_tol](double a, double b) -> bool {
    return std::fabs(b - a) <= t_tol || adjacent_doubles(a, b);
  };
  const double t_lo = std::log(b_lo);
  const double t_hi = std::log(b_hi);
  const double g_lo = g_at(t_lo);
  const double g_hi = g_at(t_hi);
  if (std::isfinite(g_lo) && std::isfinite(g_hi) && g_lo <= 0.0 &&
      g_hi >= 0.0 && std::isfinite(t_lo) && std::isfinite(t_hi)) {
    std::uintmax_t max_iter = iter_limit;
    const std::pair<double, double> root_t =
        toms748_solve(g_at, t_lo, t_hi, g_lo, g_hi, t_close, max_iter);
    if (std::isfinite(root_t.first) && std::isfinite(root_t.second) &&
        (max_iter < iter_limit || t_close(root_t.first, root_t.second))) {
      return std::exp(0.5 * (root_t.first + root_t.second));
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// The root is wanted to a relative 2^(1 - kRootDigits) = 2^-46, about 14
// digits, unless the caller asks for less.
const int kRootDigits = std::numeric_limits<double>::digits - 6;

// Solve F(x) = p_low; or, when log_p_up is not NaN, the survival equation
// P(X > x) = exp(log_p_up); or, when log_p_low is not NaN, log F(x) =
// log_p_low.  The survival form exists because converting an upper-tail
// probability by p_low = 1 - p destroys it for small p: below about 1e-16
// that subtraction returns exactly 1, where base R's qgamma still resolves the
// answer at 1e-30.  The log form serves the lower-tail log probabilities that
// lower_log_prob picks out.  The root is found to a relative 2^(1 - digits).
double solve_quantile_general(const mhn::CdfState& state, double p_low,
                              const QuantileMetadata& meta,
                              double log_p_up = R_NaN,
                              double log_p_low = R_NaN,
                              int digits = kRootDigits) {
  const bool upper = !ISNAN(log_p_up);
  const bool log_lower = !upper && !ISNAN(log_p_low);
  if (upper) {
    if (log_p_up >= 0.0) return 0.0;            // survival 1
    if (log_p_up == R_NegInf) return R_PosInf;
  } else if (!log_lower) {
    if (p_low <= 0.0) return 0.0;
    if (p_low >= 1.0) return R_PosInf;
  }

  // The bracket is set on the scale of the distribution, s = E(X) + 8 sd(X),
  // as [s sqrt(eps), s].  X / s has a law that depends on alpha and
  // z = gamma / sqrt(beta) alone (Sun et al. 2023, Lemmas 1a and 2a), so the
  // search then takes the same steps at every beta; an absolute left end of
  // sqrt(eps) lies above every quantile once the scale falls below about
  // 1e-16.
  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double x_hi = meta.initial_x_hi;
  double x_lo = x_hi * sqrt_eps;

  // Increasing in x in every form: for the upper tail the survival falls, so
  // its negated log rises.
  auto f_at = [&state, p_low, upper, log_p_up, log_lower,
               log_p_low](double x) -> double {
    if (upper) return log_p_up - state.log_upper_tail(x);
    if (log_lower) return state.log_cdf_general(x) - log_p_low;
    return state.cdf_linear(x) - p_low;
  };
  // The root finders see it kept within +-1e10, since only its sign locates
  // the root: a log tail far from the target can be of the order of the
  // largest double (-2.5e307 in the upper tail of a law narrower than the
  // spacing of the doubles), and the interpolations of TOMS 748 multiply such
  // values by differences of abscissae, which would overflow to a NaN.  The
  // values the three forms take on a bracket set as below lie inside those
  // bounds at every shape up to 1e6.
  auto clamp = [](double f) -> double {
    return (f > 1e10) ? 1e10 : ((f < -1e10) ? -1e10 : f);
  };
  auto f_root = [&f_at, &clamp](double x) -> double { return clamp(f_at(x)); };

  const QuantileBracket br = bracket_quantile(f_at, x_lo, x_hi, state);
  if (br.done) return br.value;
  x_lo = br.x_lo;
  x_hi = br.x_hi;
  const double f_lo = br.f_lo;
  const double f_hi = br.f_hi;

  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;
  const std::uintmax_t iter_limit = 100;
  std::uintmax_t max_iter = iter_limit;
  // The root is wanted to a relative 2^(1 - digits), 2^-46 unless the caller
  // asks for less.  Below the smallest normal double,
  // DBL_MIN, the doubles are evenly spaced denorm_min apart, so only an
  // absolute resolution is attainable there, and a quantile below about
  // 3.5e-310 cannot meet the relative tolerance even with its bracket one
  // spacing wide.  The bracket is therefore also accepted once no double lies
  // strictly inside it.  Above DBL_MIN the relative test is always met first,
  // except in the peak frame, where the law can be narrower than 2^-46 of x
  // and the bracket must also be within 2^-46 of its standard deviation
  // (meta.abs_tol); a law narrower than the spacing of the doubles is then
  // resolved to adjacent doubles.
  eps_tolerance<double> rel_tol(digits);
  const double abs_tol = meta.abs_tol;
  auto tol = [&rel_tol, abs_tol](double a, double b) -> bool {
    return (rel_tol(a, b) && b - a <= abs_tol) || adjacent_doubles(a, b);
  };
  // The midpoint, formed so that it cannot overflow near the largest double.
  auto middle = [](double a, double b) -> double {
    const double sum = a + b;
    return std::isfinite(sum) ? 0.5 * sum : a + 0.5 * (b - a);
  };
  // The value returned for a final bracket [a, b].  Outside the peak frame,
  // and inside it for a bracket within abs_tol, the midpoint.  A bracket of
  // adjacent doubles wider than that holds the quantile q* somewhere
  // inside, and the nearer of a and b is returned.  Where the law is wider
  // than the spacing of the doubles, f is smooth across [a, b] and q* is
  // placed by linear interpolation of f.  Where it is narrower, f rises as a
  // step and says nothing about where q* lies; but the law there is normal
  // about v* with standard deviation sd to within relative corrections below
  // 1e-3 of the offset, so q* = v* + sd t with t the normal quantile of the
  // lower-tail probability, and q* - a = (v* - a) + sd t, with v* - a from
  // the peak frame's exact offset.
  const double sd_frame = meta.sd;
  auto finish = [&](double a, double b) -> double {
    if (!std::isfinite(abs_tol) || !(b - a > abs_tol)) return middle(a, b);
    double from_a;   // (q* - a) / (b - a)
    if (b - a > sd_frame && state.peak_ok) {
      double t;
      if (upper) t = R::qnorm(log_p_up, 0.0, 1.0, 0, 1);
      else if (log_lower) t = R::qnorm(log_p_low, 0.0, 1.0, 1, 1);
      else t = R::qnorm(p_low, 0.0, 1.0, 1, 0);
      from_a = (sd_frame * t - state.peak.offset(a)) / (b - a);
    } else {
      const double fa = f_root(a);
      const double fb = f_root(b);
      if (!(fb > fa)) return middle(a, b);
      from_a = -fa / (fb - fa);
    }
    if (ISNAN(from_a)) return middle(a, b);
    return (from_a <= 0.5) ? a : b;
  };

  // Plain bisection of [a, b], f(a) <= 0 <= f(b), to the tolerance: at most
  // about 60 halvings from a bracket a factor of two wide, one CDF
  // evaluation each.  It serves where TOMS 748 does not converge -- on a
  // distribution function that rises as a step within one spacing of the
  // doubles, which interpolation cannot follow -- and where the bracket lies
  // near the largest double, where the products formed by its
  // interpolations would overflow.
  auto bisect = [&f_root, &tol, &finish, &state](double a,
                                                 double b) -> double {
    for (int it = 0; it < 2200 && !tol(a, b); ++it) {
      const double c =
          (mhn::opaque_value(a + b) <= std::numeric_limits<double>::max())
              ? 0.5 * (a + b)
              : a + 0.5 * (b - a);
      if (!(c > a && c < b)) break;
      const double fc = f_root(c);
      if (ISNAN(fc)) {
        mhn::warn("qmhn: the distribution function could not be "
                  "evaluated for alpha=%g, beta=%g, gamma=%g; "
                  "returning NA", state.alpha, state.beta, state.gamma);
        return NA_REAL;
      }
      if (fc < 0.0) a = c; else b = c;
    }
    return finish(a, b);
  };

  if (x_hi > 1e298) return bisect(x_lo, x_hi);

  std::pair<double, double> root = toms748_solve(
      f_root, x_lo, x_hi, clamp(f_lo), clamp(f_hi), tol, max_iter);
  const bool root_ok = std::isfinite(root.first) &&
                       std::isfinite(root.second) &&
                       root.first >= x_lo && root.second <= x_hi &&
                       root.first <= root.second;
  if (root_ok && (max_iter < iter_limit || tol(root.first, root.second))) {
    return finish(root.first, root.second);
  }

  // TOMS 748 used up its iterations without meeting the tolerance.  This
  // happens for large alpha and a tiny lower-tail p: F(x) behaves like
  // x^alpha near the origin, so f is flat at -p across most of the bracket,
  // and interpolation on it narrows the bracket only slowly.  In t = log x
  // the same equation, log F(e^t) = log p, is close to linear with slope
  // alpha, and the root finder closes on it in a few steps.  It is solved
  // there from the bracket reached so far, to the same relative accuracy in
  // x, which is an absolute accuracy in t -- or, from |t| = 128 on, where the
  // doubles near t are spaced more widely than that, to adjacent doubles.
  // Failing that, and in the peak frame, where adjacent doubles in t can lie
  // many of the law's standard deviations apart in x, the bracket reached so
  // far is bisected.
  const double b_lo = root_ok ? root.first : x_lo;
  const double b_hi = root_ok ? root.second : x_hi;
  if (std::isfinite(abs_tol)) return bisect(b_lo, b_hi);
  const double log_target = upper ? 0.0
                            : (log_lower ? log_p_low : std::log(p_low));
  const double x = solve_in_log_x(state, b_lo, b_hi, upper, log_p_up,
                                  log_target, digits, iter_limit);
  return ISNAN(x) ? bisect(b_lo, b_hi) : x;
}

// log F(x) as the polish takes it: log_lower_resolving below x_split where
// the series would cost milliseconds, and log_cdf_general otherwise.
double polish_log_cdf(const mhn::CdfState& state, double x) {
  return (state.series_costly && x <= state.x_split)
      ? state.log_lower_resolving(x, state.log_upper_tail(x))
      : state.log_cdf_general(x);
}

// Newton's method on log F(x) = log_d in t = log x, from x, the root of the
// survival equation for the upper-tail probability 1 - d.  The survival
// root carries a relative error of up to upper_log_error / d in d (see
// upper_route_limit): at most about 7e-3 for d >= kSurvivalSolveLimit,
// since upper_log_error is at most about 7.4e-10 for the states polished
// (at |z| = 832, where the series stops serving, with beta at either end of
// the double range), and d is at most 1/2.  In t, g(t) = log F(e^t) - log_d
// has slope x f(x) / F(x), and a step lands within C times the square of
// the error it starts from, C = |g''| / (2 g').  Next to the origin, where
// F is about a power of x, g is close to linear and C small; where the law
// is close to normal, C = (x / sd) |w + phi(w) / Phi(w)| / 2 at the
// standardised quantile w, which increases with w to 0.4 x / sd at
// d = 1/2, and x / sd is below about 710 for the states polished, whose
// peak lies below 500 scale units.  Steps are taken until one is below
// 5e-9, which leaves a relative error below about 7e-15 (C at most 284),
// within the tolerance of the root finder, or no longer moves x.
//
// With from_afar set, x can lie far from the root: for d below
// kSurvivalSolveLimit x is the root of the survival equation for 1 - d or
// for a larger d' (kSurvivalStartErrors), so log F(x) can differ from log_d
// by up to 1/16 or by log(d' / d), and the caller also sets it for a second
// try where the bounded steps fail.  The steps are then not bounded, and up
// to eight are taken.  g is increasing, and for gamma <= 0 also concave:
// the density of log X, proportional to exp(alpha t - beta e^(2t) +
// gamma e^t), is then log-concave, and so is its distribution function.
// Newton's method on it lands at or below the root after the first step and
// then rises to it without passing it.  For gamma > 0 the density of log X
// is log-concave only above x = gamma / (4 beta).  Next to the origin F is
// about c x^alpha, and g close to the line alpha t + log c - log_d.  A step
// to the left that would leave the positive doubles stops at the smallest
// of them, denorm_min: if F there is still at least d, the root lies below
// it and 0 is returned, as bracket_quantile reports such a quantile, and
// otherwise the steps go on from there.  Returns NaN where a step cannot be
// formed, leaves the positive doubles otherwise or does not settle, and the
// caller then solves the lower-tail equation instead.
double polish_on_lower_tail(const mhn::CdfState& state, double x,
                            double log_d, bool from_afar) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    return from_afar ? std::numeric_limits<double>::quiet_NaN() : x;
  }
  const int max_steps = from_afar ? 8 : 4;
  for (int k = 0; k < max_steps; ++k) {
    const double log_x = std::log(x);
    const double log_F = polish_log_cdf(state, x);
    // The ordinary log density, as dmhn forms it for these states (no shift
    // and no peak frame), to the rounding of its terms: a relative error
    // far below what a Newton step needs.
    const double log_f = state.log_prefactor + (state.alpha - 1.0) * log_x -
                         state.beta * x * x + state.gamma * x;
    const double step = (log_F - log_d) * std::exp(log_F - log_f - log_x);
    if (!std::isfinite(step) || (!from_afar && std::fabs(step) > 0.1)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    // A step that no longer moves x also ends the search: below
    // .Machine$double.xmin the doubles are 4.9e-324 apart, and the root is
    // resolved only to that spacing (a relative 2e-4 at x = 2.5e-320).
    double x_new = x * std::exp(-step);
    const double x_min = std::numeric_limits<double>::denorm_min();
    if (from_afar && step > 0.0 && !(x_new >= x_min)) {
      if (x == x_min || polish_log_cdf(state, x_min) >= log_d) return 0.0;
      x_new = x_min;
    }
    if (!(x_new > 0.0) || !std::isfinite(x_new)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::fabs(step) <= 5e-9 || x_new == x) return x_new;
    x = x_new;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// The quantile of one probability p of the general case, for the state and
// metadata of its triple: the root of the equation the routing above
// chooses, polished on the lower tail where polishes_on_lower_tail asks.
double quantile_general_one(const mhn::CdfState& state,
                            const QuantileMetadata& meta, double p,
                            bool lower_tail, bool log_p) {
  const double p_low = to_lower_tail_prob(p, lower_tail, log_p);
  // R_IsNA / R_IsNaN must be tested in this order: NA_REAL also passes
  // R_IsNaN, but only R_IsNA distinguishes it from a generic NaN.
  if (R_IsNA(p_low)) return NA_REAL;
  if (R_IsNaN(p_low)) return R_NaN;
  // A subnormal alpha past the hand-over of mhn_vanishing.h: a lower tail
  // small enough that its quantile lies where F depends on alpha is solved
  // in the forms there.
  double x_far;
  if (state.vanishing_lower &&
      state.vanishing.q_lower_far(p, lower_tail, log_p, &x_far)) {
    return x_far;
  }
  const double limit = upper_route_limit(state);
  if (!lower_tail && p_low > 0.0 && near_one(p, log_p, limit) &&
      polishes_on_lower_tail(state)) {
    // Within kSurvivalSolveLimit of 1 the survival root only starts the
    // polish, and a failed search or polish falls back on the lower-tail
    // solve; the survival equation is then solved for a 1 - p it resolves.
    const bool from_afar = near_one(p, log_p, kSurvivalSolveLimit);
    const double d_start = kSurvivalStartErrors * state.upper_log_error;
    const double log_p_up = (from_afar && p_low < d_start)
        ? std::log1p(-d_start) : (log_p ? p : std::log(p));
    double x = solve_quantile_general(state, p_low, meta, log_p_up, R_NaN,
                                      kPolishStartDigits);
    if (ISNAN(x) && !from_afar) return x;
    if (from_afar && x == 0.0) {
      // The survival root lies below the smallest positive double, which
      // happens for a small shape, where F is about c x^alpha.  If F there
      // is still at least 1 - p, the quantile lies below it too and is 0 in
      // double precision, as bracket_quantile reports it; otherwise the
      // polish starts from that double.
      const double x_min = std::numeric_limits<double>::denorm_min();
      if (polish_log_cdf(state, x_min) >= std::log(p_low)) return 0.0;
      x = x_min;
    }
    // Where the bounded steps fail, as for a subnormal root, which the
    // survival root resolves only to the spacing of the doubles there, so
    // that the first step exceeds the bound, the unbounded ones are tried
    // from the same x before the lower-tail equation is solved from
    // scratch: for a series state that solve halves its bracket hundreds of
    // times, one sum of the series each, which can take seconds.
    double polished =
        polish_on_lower_tail(state, x, std::log(p_low), from_afar);
    if (ISNAN(polished) && !from_afar) {
      polished = polish_on_lower_tail(state, x, std::log(p_low), true);
    }
    if (!ISNAN(polished)) return polished;
  }
  return solve_quantile_general(
      state, p_low, meta,
      upper_log_prob(p, lower_tail, log_p, limit, lower_route_limit(state)),
      lower_log_prob(p, lower_tail, log_p, limit));
}

}  // namespace

// [[Rcpp::export(.qmhn_cpp, rng = false)]]
Rcpp::NumericVector qmhn_cpp(Rcpp::NumericVector p,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool lower_tail,
                             bool log_p) {
  const R_xlen_t np = p.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // A zero-length argument gives numeric(0), as in qnorm.
  if (np == 0 || na == 0 || nb == 0 || ng == 0) return Rcpp::NumericVector(0);

  // Invalid parameters are handled element by element, as in dmhn_cpp.
  const bool all_valid = mhn::all_params_valid(alpha, beta, gamma);

  // Fast path: all parameters scalar.
  if (all_valid && na == 1 && nb == 1 && ng == 1) {
    const double a = alpha[0], b = beta[0], g = gamma[0];
    // A subnormal alpha goes through the state, which takes the forms of
    // mhn_vanishing.h where they apply, before the special cases.
    if (!mhn::vanishing_candidate(a)) {
      if (mhn::is_sqrt_gamma_exact(g) && !mhn::sqrt_gamma_in_peak_frame(a)) {
        return mhn::qmhn_sqrt_gamma(p, a, b, lower_tail, log_p);
      }
      if (mhn::is_truncated_normal_exact(a)) {
        return mhn::qmhn_truncated_normal(p, b, g, lower_tail, log_p);
      }
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    if (state.kind == mhn::CDF_VANISHING) {
      Rcpp::NumericVector out(np);
      for (R_xlen_t i = 0; i < np; ++i) {
        Rcpp::checkUserInterrupt();
        out[i] = state.vanishing.q_one(p[i], lower_tail, log_p);
      }
      return out;
    }
    QuantileMetadata meta;
    recompute_metadata(state, meta);
    Rcpp::NumericVector out(np);
    for (R_xlen_t i = 0; i < np; ++i) {
      // Each element is a bracketed root find over the CDF, and where that is
      // evaluated by quadrature an element can take milliseconds.
      Rcpp::checkUserInterrupt();
      out[i] = quantile_general_one(state, meta, p[i], lower_tail, log_p);
    }
    return out;
  }

  const R_xlen_t n = std::max({np, na, nb, ng});
  Rcpp::NumericVector out(n);

  // One cache slot per distinct triple in the recycling cycle.  The metadata
  // travels with the state it was derived from, so a slot carries both.
  mhn::ParamSlots<QuantileState> slots(na, nb, ng, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double pi = p[i % np];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!all_valid && !mhn::params_valid(a, b, g)) {
      out[i] = mhn::invalid_result(pi, a, b, g);
      continue;
    }
    QuantileState& qs = slots.at(i, a, b, g);
    mhn::CdfState& state = qs.state;
    QuantileMetadata& meta = qs.meta;

    // The special cases go to the per-element quantiles that the all-scalar
    // path loops over, with the constants the slot already holds for this
    // triple, and honour lower.tail / log.p themselves.
    if (state.kind == mhn::CDF_SQRT_GAMMA) {
      out[i] = mhn::qmhn_sqrt_gamma_one(pi, state.sg_shape, state.sg_sqrt_beta,
                                        lower_tail, log_p);
      continue;
    }
    if (state.kind == mhn::CDF_TRUNCATED_NORMAL) {
      out[i] = mhn::qmhn_truncated_normal_one(pi, state.tn, lower_tail, log_p);
      continue;
    }
    if (state.kind == mhn::CDF_VANISHING) {
      Rcpp::checkUserInterrupt();
      out[i] = state.vanishing.q_one(pi, lower_tail, log_p);
      continue;
    }

    const double p_low = to_lower_tail_prob(pi, lower_tail, log_p);
    if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }
    if (!meta.valid) recompute_metadata(state, meta);
    // Each element runs a bracketed root find over the CDF, so this loop is
    // the slowest in the package; check for an interrupt every element.
    Rcpp::checkUserInterrupt();
    out[i] = quantile_general_one(state, meta, pi, lower_tail, log_p);
  }
  return out;
}
