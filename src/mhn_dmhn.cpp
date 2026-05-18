// MHN density f(x | alpha, beta, gamma), computed in log space via the
// Fox-Wright Psi normalizing constant.  Special cases (gamma = 0 ->
// sqrt-Gamma; alpha = 1 -> truncated normal; both -> half-normal) are
// detected and dispatched to mhn_special_cases.cpp for closed-form
// evaluation.  Vectorized over both x and (alpha, beta, gamma) with the
// ParamCache pattern -- mirrored by mhn_pmhn.cpp.
//
// R-side wrapper: dmhn() in mhn/R/dmhn.R.

#include "mhn_check.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>

namespace {

// Special-case kind detected for a given (alpha, gamma) pair.
enum SpecialKind {
  KIND_GENERAL = 0,
  KIND_SQRT_GAMMA = 1,      // gamma == 0
  KIND_TRUNCATED_NORMAL = 2 // alpha == 1
};

// Per-parameter cached state used by the vectorized loop.
struct ParamCache {
  // General-case constant: log(2) + (alpha/2) log(beta) - log Psi
  double log_const = 0.0;
  // Truncated-normal pre-computed quantities
  double tn_mu = 0.0;
  double tn_sigma = 0.0;
  double tn_log_norm = 0.0;
  // Boundary value at x=0 for the sqrt-Gamma branch
  double sqg_log_zero = R_NegInf;
  SpecialKind kind = KIND_GENERAL;

  void recompute(double alpha, double beta, double gamma) {
    if (mhn::is_sqrt_gamma(gamma)) {
      kind = KIND_SQRT_GAMMA;
      if (alpha > 1.0) {
        sqg_log_zero = R_NegInf;
      } else if (alpha < 1.0) {
        sqg_log_zero = R_PosInf;
      } else {
        sqg_log_zero = std::log(2.0) + 0.5 * std::log(beta) - 0.5 * std::log(M_PI);
      }
    } else if (mhn::is_truncated_normal(alpha)) {
      kind = KIND_TRUNCATED_NORMAL;
      tn_mu = gamma / (2.0 * beta);
      tn_sigma = 1.0 / std::sqrt(2.0 * beta);
      tn_log_norm = R::pnorm(tn_mu / tn_sigma, 0.0, 1.0,
                             /*lower_tail=*/1, /*log_p=*/1);
    } else {
      kind = KIND_GENERAL;
      const double log_nc = mhn::mhn_log_normalizing_const(alpha, beta, gamma, -1.0);
      log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;
    }
  }

  double log_density_at(double xi, double alpha, double beta, double gamma) const {
    if (Rcpp::NumericVector::is_na(xi)) return NA_REAL;
    switch (kind) {
      case KIND_SQRT_GAMMA:
        if (xi > 0.0) {
          return std::log(2.0) + std::log(xi)
                 + R::dgamma(xi * xi, alpha / 2.0,
                             /*scale=*/1.0 / beta, /*log=*/1);
        }
        if (xi == 0.0) return sqg_log_zero;
        return R_NegInf;
      case KIND_TRUNCATED_NORMAL:
        if (xi >= 0.0) {
          return R::dnorm(xi, tn_mu, tn_sigma, /*log=*/1) - tn_log_norm;
        }
        return R_NegInf;
      case KIND_GENERAL:
      default:
        if (xi > 0.0) {
          return log_const + (alpha - 1.0) * std::log(xi)
                 - beta * xi * xi + gamma * xi;
        }
        if (xi == 0.0) {
          if (alpha > 1.0) return R_NegInf;
          if (alpha < 1.0) return R_PosInf;
          return log_const;  // alpha == 1 is normally caught by KIND_TRUNCATED_NORMAL
        }
        return R_NegInf;
    }
  }
};

// Fast path when all of alpha, beta, gamma are scalars.
// Identical layout to the pre-vectorization implementation, kept hot
// because this is the dominant call shape in practice (MCMC density
// evaluations with fixed parameters).
Rcpp::NumericVector dmhn_scalar_path(const Rcpp::NumericVector& x,
                                     double alpha, double beta, double gamma,
                                     bool log_p) {
  if (mhn::is_sqrt_gamma(gamma)) {
    return mhn::dmhn_sqrt_gamma(x, alpha, beta, log_p);
  }
  if (mhn::is_truncated_normal(alpha)) {
    return mhn::dmhn_truncated_normal(x, beta, gamma, log_p);
  }

  const double log_nc = mhn::mhn_log_normalizing_const(alpha, beta, gamma, -1.0);
  const double log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      log_f[i] = NA_REAL;
    } else if (xi > 0.0) {
      log_f[i] = log_const + (alpha - 1.0) * std::log(xi)
                 - beta * xi * xi + gamma * xi;
    } else if (xi == 0.0) {
      if (alpha > 1.0) log_f[i] = R_NegInf;
      else if (alpha < 1.0) log_f[i] = R_PosInf;
      // alpha == 1 is handled by the truncated-normal dispatch above.
    }
  }

  if (!log_p) {
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}

}  // namespace

// [[Rcpp::export(.dmhn_cpp)]]
Rcpp::NumericVector dmhn_cpp(Rcpp::NumericVector x,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool log_p) {
  mhn::check_params_vector(alpha, beta, gamma);

  const R_xlen_t nx = x.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // Empty x with non-empty params returns numeric(0), matching R conventions.
  if (nx == 0) return Rcpp::NumericVector(0);

  // Fast path: all parameters scalar.
  if (na == 1 && nb == 1 && ng == 1) {
    return dmhn_scalar_path(x, alpha[0], beta[0], gamma[0], log_p);
  }

  // Vectorized path. Recycle each input to length n = max(nx, na, nb, ng).
  const R_xlen_t n = std::max({nx, na, nb, ng});
  Rcpp::NumericVector log_f(n);

  // Cache the most recently used parameter triple to avoid recomputing
  // log Psi when consecutive elements share parameters (common in
  // grouped-parameter sweeps).
  ParamCache cache;
  double prev_a = std::numeric_limits<double>::quiet_NaN();
  double prev_b = std::numeric_limits<double>::quiet_NaN();
  double prev_g = std::numeric_limits<double>::quiet_NaN();
  bool primed = false;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i % nx];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!primed || a != prev_a || b != prev_b || g != prev_g) {
      cache.recompute(a, b, g);
      prev_a = a; prev_b = b; prev_g = g;
      primed = true;
    }

    log_f[i] = cache.log_density_at(xi, a, b, g);
  }

  if (!log_p) {
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}
