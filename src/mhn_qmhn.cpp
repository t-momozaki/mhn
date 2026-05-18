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

#include <Rcpp.h>
#include <Rmath.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {

// Cached mean / variance for the general path so the bracket can be reused
// across all elements that share a parameter triple.
struct QuantileMetadata {
  bool valid = false;
  double mean = 0.0;
  double var = 0.0;
  double initial_x_hi = 1.0;
};

void recompute_metadata(double alpha, double beta, double gamma,
                        QuantileMetadata& meta) {
  const double log_psi_a = mhn::mhn_log_normalizing_const(alpha, beta, gamma,
                                                          -1.0);
  const double log_psi_a1 = mhn::mhn_log_normalizing_const(alpha + 1.0, beta,
                                                           gamma, -1.0);
  const double mu = std::exp(log_psi_a1 - 0.5 * std::log(beta) - log_psi_a);
  const double v_raw = alpha / (2.0 * beta) + mu * (gamma / (2.0 * beta) - mu);
  const double v = std::max(0.0, v_raw);
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

double solve_quantile_general(const mhn::CdfState& state, double p_low,
                              const QuantileMetadata& meta) {
  if (p_low <= 0.0) return 0.0;
  if (p_low >= 1.0) return R_PosInf;

  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double x_lo = sqrt_eps;
  double x_hi = std::max(meta.initial_x_hi, x_lo * 2.0);

  auto f_at = [&state, p_low](double x) -> double {
    return state.cdf_linear(x) - p_low;
  };

  double f_lo = f_at(x_lo);
  double f_hi = f_at(x_hi);

  // If the lower endpoint isn't below the target (can happen for alpha < 1
  // where F rises rapidly near 0), shrink it geometrically a few times.
  int shrink = 0;
  while (f_lo > 0.0 && x_lo > std::numeric_limits<double>::min() && shrink < 30) {
    x_lo *= 0.5;
    f_lo = f_at(x_lo);
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
    return x_hi;  // best-effort fallback
  }
  if (f_lo > 0.0) {
    return x_lo;  // root is below x_lo; sqrt(eps) is the best we can offer
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

// [[Rcpp::export(.qmhn_cpp)]]
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
    if (mhn::is_sqrt_gamma(g)) {
      return mhn::qmhn_sqrt_gamma(p, a, b, lower_tail, log_p);
    }
    if (mhn::is_truncated_normal(a)) {
      return mhn::qmhn_truncated_normal(p, b, g, lower_tail, log_p);
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    QuantileMetadata meta;
    recompute_metadata(a, b, g, meta);
    Rcpp::NumericVector out(np);
    for (R_xlen_t i = 0; i < np; ++i) {
      const double p_low = to_lower_tail_prob(p[i], lower_tail, log_p);
      // R_IsNA / R_IsNaN must be tested in this order: NA_REAL also passes
      // R_IsNaN, but only R_IsNA distinguishes it from a generic NaN.
      if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
      if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }
      out[i] = solve_quantile_general(state, p_low, meta);
    }
    return out;
  }

  const R_xlen_t n = std::max({np, na, nb, ng});
  Rcpp::NumericVector out(n);

  mhn::CdfState state;
  QuantileMetadata meta;
  double prev_a = std::numeric_limits<double>::quiet_NaN();
  double prev_b = std::numeric_limits<double>::quiet_NaN();
  double prev_g = std::numeric_limits<double>::quiet_NaN();
  bool primed = false;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double pi = p[i % np];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!primed || a != prev_a || b != prev_b || g != prev_g) {
      state.recompute(a, b, g);
      meta.valid = false;  // metadata recomputed only when needed (general case)
      prev_a = a; prev_b = b; prev_g = g;
      primed = true;
    }

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

    if (!meta.valid) recompute_metadata(a, b, g, meta);
    out[i] = solve_quantile_general(state, p_low, meta);
  }
  return out;
}
