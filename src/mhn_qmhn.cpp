// Quantile function (qmhn) -- vectorized over p and (alpha, beta, gamma).
//
// Dispatch:
//   * gamma  ~ 0   -> sqrt-Gamma special case (Sun et al. 2023, Lemma 6a):
//                     sqrt(qgamma(p, alpha/2, scale = 1/beta)).
//   * alpha  ~ 1   -> truncated normal inverse (Sun et al. 2023, Lemma 6b).
//   * otherwise    -> Boost TOMS 748 root-finder on F(x) - p_target = 0,
//                     using CdfState from mhn_pmhn_state.h to evaluate F.
//
// Search bracket starts at [sqrt(eps), E(X) + 8 sqrt(Var(X))] (per spec
// §3.3); the upper end is doubled up to 30 times if F(x_hi) < p_target.

#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_pmhn_state.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_param_slots.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <algorithm>
#include <cmath>
#include <exception>
#include <cstdint>
#include <limits>

// Defined in mhn_moments.cpp.  Declared rather than included because both are
// Rcpp exports at file scope; the quantile bracket wants the guarded forms, not
// a second copy of the algebra.
double mhn_mean_cpp(double alpha, double beta, double gamma);
double mhn_var_cpp(double alpha, double beta, double gamma);

namespace {

// Cached mean / variance for the general path so the bracket can be reused
// across all elements that share a parameter triple.
struct QuantileMetadata {
  bool valid = false;
  double mean = 0.0;
  double var = 0.0;
  double initial_x_hi = 1.0;
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

// The caller has just built a CdfState for the same triple, and for the
// general case that state already holds log Psi(alpha).  Evaluating it a second
// time here doubled the setup cost of every distinct triple, and log Psi is the
// expensive part -- a series or a quadrature, not a closed form.
void recompute_metadata(const mhn::CdfState& state, QuantileMetadata& meta) {
  const double alpha = state.alpha;
  const double beta = state.beta;
  const double gamma = state.gamma;
  // Both the mean and the variance are taken from the exported helpers rather
  // than rebuilt here.  Forming exp(log Psi(alpha+1) - log Psi(alpha)) directly
  // subtracts two quantities that are each about gamma^2/(4 beta) for a large
  // positive tilt -- 2.5e9 at gamma = 1e5 -- and mhn_mean_cpp already switches
  // to a quadrature mean past the point where that costs more than the working
  // tolerance, while mhn_var_cpp switches on the measured cancellation in the
  // Lemma 2c form.  This is only the bracket for the root find, but a bracket
  // built from a mean that is 31 percent wrong is a slow start at best.
  const double mu = mhn_mean_cpp(alpha, beta, gamma);
  const double v = mhn_var_cpp(alpha, beta, gamma);
  meta.mean = mu;
  meta.var = v;
  double x_hi = mu + 8.0 * std::sqrt(v);
  if (!(x_hi > 0.0) || !std::isfinite(x_hi)) {
    x_hi = std::max(mu, 1.0);
    if (!(x_hi > 0.0)) x_hi = 1.0;
  }
  meta.initial_x_hi = x_hi;
  meta.valid = true;
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
  if (!lower_tail) p_low = 1.0 - p_low;
  if (p_low < 0.0) p_low = 0.0;
  if (p_low > 1.0) p_low = 1.0;
  return p_low;
}

// The upper-tail probability on the log scale, or NaN when the caller asked
// for the lower tail.  Keeping it in logs is the whole point: an upper-tail p
// of 1e-30 is perfectly representable, and only the conversion to 1 - p
// destroys it.
double upper_log_prob(double p, bool lower_tail, bool log_p) {
  if (lower_tail) return R_NaN;
  if (R_IsNA(p) || R_IsNaN(p)) return R_NaN;
  if (log_p) return (p > 0.0) ? R_NaN : p;
  if (p < 0.0 || p > 1.0) return R_NaN;
  return std::log(p);
}

// Solve F(x) = p_low, or -- when log_p_up is finite -- the survival equation
// P(X > x) = exp(log_p_up).  The second form exists because converting an
// upper-tail probability by p_low = 1 - p destroys it for small p: below about
// 1e-16 that subtraction returns exactly 1 and the quantile came back as Inf,
// where base R's qgamma still resolves the answer at 1e-30.
double solve_quantile_general(const mhn::CdfState& state, double p_low,
                              const QuantileMetadata& meta,
                              double log_p_up = R_NaN) {
  const bool upper = !ISNAN(log_p_up);
  if (!upper) {
    if (p_low <= 0.0) return 0.0;
    if (p_low >= 1.0) return R_PosInf;
  } else {
    if (log_p_up >= 0.0) return 0.0;            // survival 1
    if (log_p_up == R_NegInf) return R_PosInf;
  }

  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double x_lo = sqrt_eps;
  double x_hi = std::max(meta.initial_x_hi, x_lo * 2.0);

  // Increasing in x either way: for the upper tail the survival falls, so its
  // negated log rises.
  auto f_at = [&state, p_low, upper, log_p_up](double x) -> double {
    if (upper) return log_p_up - state.log_upper_tail(x);
    return state.cdf_linear(x) - p_low;
  };

  double f_lo = f_at(x_lo);
  double f_hi = f_at(x_hi);

  // If the lower endpoint isn't below the target, shrink it geometrically.
  // For alpha < 1 the distribution function behaves like x^alpha near the
  // origin, so the quantile of a small p sits at about p^(1/alpha): at
  // alpha = 0.1 and p = 1e-4 that is of order 1e-40, far below anything
  // 30 halvings of sqrt(eps) can reach.  Stopping there and returning the
  // bracket end gave answers wrong by tens of orders of magnitude -- and
  // silently, since the caller could not tell a converged root from a
  // surrendered bracket.  Allow the search to run down to the denormal floor
  // instead; each step is one CDF evaluation and the loop stops as soon as it
  // straddles the root.
  int shrink = 0;
  while (f_lo > 0.0 && x_lo > std::numeric_limits<double>::denorm_min() &&
         shrink < 1100) {
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
    x_lo = x_try;
    f_lo = f_try;
    ++shrink;
  }

  // Expand x_hi until we bracket the root.
  int expand = 0;
  while (f_hi < 0.0 && expand < 30 && std::isfinite(x_hi)) {
    x_hi *= 2.0;
    f_hi = f_at(x_hi);
    ++expand;
  }
  if (f_hi < 0.0 || !std::isfinite(f_hi)) {
    // Returning x_hi here passed off the untouched initial bracket as the
    // answer: the value came back identical for every p, with no warning.
    // There is no root to report, so say so.
    Rcpp::warning("qmhn: could not bracket the quantile for alpha=%g, "
                  "beta=%g, gamma=%g; returning NA",
                  state.alpha, state.beta, state.gamma);
    return R_NaN;
  }
  if (f_lo > 0.0) {
    return x_lo;  // root is below the smallest bracket the shrink loop reached
  }

  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;
  const int digits = std::numeric_limits<double>::digits - 6;  // ~14 digits
  std::uintmax_t max_iter = 100;
  eps_tolerance<double> tol(digits);
  std::pair<double, double> root = toms748_solve(f_at, x_lo, x_hi,
                                                 f_lo, f_hi, tol, max_iter);
  return 0.5 * (root.first + root.second);
}

}  // namespace

// [[Rcpp::export(.qmhn_cpp, rng = false)]]
Rcpp::NumericVector qmhn_cpp(Rcpp::NumericVector p,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool lower_tail,
                             bool log_p) {
  mhn::check_params_vector(alpha, beta, gamma);

  const R_xlen_t np = p.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  if (np == 0) return Rcpp::NumericVector(0);

  // Fast path: all parameters scalar.
  if (na == 1 && nb == 1 && ng == 1) {
    const double a = alpha[0], b = beta[0], g = gamma[0];
    if (mhn::is_sqrt_gamma(g, b)) {
      return mhn::qmhn_sqrt_gamma(p, a, b, lower_tail, log_p);
    }
    if (mhn::is_truncated_normal(a)) {
      return mhn::qmhn_truncated_normal(p, b, g, lower_tail, log_p);
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    QuantileMetadata meta;
    recompute_metadata(state, meta);
    Rcpp::NumericVector out(np);
    for (R_xlen_t i = 0; i < np; ++i) {
      const double p_low = to_lower_tail_prob(p[i], lower_tail, log_p);
      // R_IsNA / R_IsNaN must be tested in this order: NA_REAL also passes
      // R_IsNaN, but only R_IsNA distinguishes it from a generic NaN.
      if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
      if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }
      out[i] = solve_quantile_general(state, p_low, meta,
                                      upper_log_prob(p[i], lower_tail, log_p));
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

    QuantileState& qs = slots.at(i, a, b, g);
    mhn::CdfState& state = qs.state;
    QuantileMetadata& meta = qs.meta;

    if (state.kind == mhn::CDF_SQRT_GAMMA) {
      // Inline element-wise sqrt-gamma quantile via the helper that already
      // honours lower.tail / log.p.
      Rcpp::NumericVector p_one(1);
      p_one[0] = pi;
      Rcpp::NumericVector q_one = mhn::qmhn_sqrt_gamma(p_one, a, b,
                                                       lower_tail, log_p);
      out[i] = q_one[0];
      continue;
    }
    if (state.kind == mhn::CDF_TRUNCATED_NORMAL) {
      Rcpp::NumericVector p_one(1);
      p_one[0] = pi;
      Rcpp::NumericVector q_one = mhn::qmhn_truncated_normal(p_one, b, g,
                                                             lower_tail, log_p);
      out[i] = q_one[0];
      continue;
    }

    const double p_low = to_lower_tail_prob(pi, lower_tail, log_p);
    if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }

    if (!meta.valid) recompute_metadata(state, meta);
    // Each element runs a bracketed root find over the CDF, so this loop is
    // the slowest in the package; check for an interrupt every element.
    Rcpp::checkUserInterrupt();
    out[i] = solve_quantile_general(state, p_low, meta,
                                    upper_log_prob(pi, lower_tail, log_p));
  }
  return out;
}
