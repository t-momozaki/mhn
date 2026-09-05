// MHN density f(x | alpha, beta, gamma), computed in log space via the
// Fox-Wright Psi normalizing constant.  Special cases (gamma = 0 ->
// sqrt-Gamma; alpha = 1 -> truncated normal; both -> half-normal) are
// detected and dispatched to mhn_special_cases.cpp for closed-form
// evaluation.  Vectorized over both x and (alpha, beta, gamma) with the
// ParamCache pattern -- mirrored by mhn_pmhn.cpp.
//
// R-side wrapper: dmhn() in mhn/R/dmhn.R.

#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_param_slots.h"

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

// The general-case log density.
//
// Both the scalar fast path and the recycling loop need this, and when each
// carried its own copy they were free to drift apart -- which is how pmhn came
// to return two different answers for the same truncated-normal tail.  One
// definition, used by both.
struct GeneralDensity {
  double log_const = 0.0;
  double kernel_centre = 0.0;
  bool completed_square = false;

  void set(double alpha, double beta, double gamma) {
    // For a large positive tilt, log Psi and the kernel exponent are each about
    // gamma^2/(4 beta) while the log density is of order one, so assembling the
    // two directly leaves an absolute error of about eps gamma^2/(4 beta).
    // That passes the working tolerance at |z| = gamma/sqrt(beta) around 2e4
    // and reaches 5 percent by 3e7.  Past a threshold well inside the safe
    // range, take log Psi with that term removed and complete the square in the
    // kernel, so it cancels analytically instead of numerically.  Below the
    // threshold the ordinary assembly is used unchanged, and is the more
    // accurate of the two there.
    const double z = gamma / std::sqrt(beta);
    completed_square = (z > 1.0e3);
    if (completed_square) {
      const double log_nc =
          mhn::psi_integrate_shifted(alpha, beta, gamma, mhn::mhn_eps());
      log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;
      kernel_centre = gamma / (2.0 * beta);
    } else {
      const double log_nc =
          mhn::mhn_log_normalizing_const(alpha, beta, gamma, -1.0);
      log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;
      kernel_centre = 0.0;
    }
  }

  double at(double xi, double alpha, double beta, double gamma) const {
    // The density vanishes at infinity, as every base R density does.  Left to
    // the general expression it came back as NaN: (alpha-1) log(x) and gamma x
    // are +Inf while -beta x^2 is -Inf.
    if (xi == R_PosInf) return R_NegInf;
    if (xi > 0.0) {
      if (completed_square) {
        const double d = xi - kernel_centre;
        return log_const + (alpha - 1.0) * std::log(xi) - beta * d * d;
      }
      return log_const + (alpha - 1.0) * std::log(xi)
             - beta * xi * xi + gamma * xi;
    }
    if (xi == 0.0) {
      if (alpha > 1.0) return R_NegInf;
      if (alpha < 1.0) return R_PosInf;
      return log_const;  // alpha == 1 is normally intercepted as a special case
    }
    return R_NegInf;
  }
};

// Per-parameter cached state used by the vectorized loop.
struct ParamCache {
  // Truncated-normal pre-computed quantities
  double tn_mu = 0.0;
  double tn_sigma = 0.0;
  double tn_log_norm = 0.0;
  // Boundary value at x=0 for the sqrt-Gamma branch
  double sqg_log_zero = R_NegInf;
  SpecialKind kind = KIND_GENERAL;
  GeneralDensity general;

  void recompute(double alpha, double beta, double gamma) {
    if (mhn::is_sqrt_gamma(gamma, beta)) {
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
      general.set(alpha, beta, gamma);
    }
  }

  double log_density_at(double xi, double alpha, double beta, double gamma) const {
    // NA and NaN are distinct in the base R d/p/q contract, and Rcpp's
    // is_na does not separate them -- it is ISNAN -- so each is tested
    // with the R_Is* predicate that isolates it.
    if (R_IsNA(xi)) return NA_REAL;
    if (R_IsNaN(xi)) return R_NaN;
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
        return general.at(xi, alpha, beta, gamma);
    }
  }
};

// Fast path when all of alpha, beta, gamma are scalars.
// Kept separate from the recycling loop because this is the dominant call
// shape in practice: MCMC density evaluations at fixed parameters, where
// the modular indexing and the cache check would be pure overhead.
Rcpp::NumericVector dmhn_scalar_path(const Rcpp::NumericVector& x,
                                     double alpha, double beta, double gamma,
                                     bool log_p) {
  if (mhn::is_sqrt_gamma(gamma, beta)) {
    return mhn::dmhn_sqrt_gamma(x, alpha, beta, log_p);
  }
  if (mhn::is_truncated_normal(alpha)) {
    return mhn::dmhn_truncated_normal(x, beta, gamma, log_p);
  }

  GeneralDensity general;
  general.set(alpha, beta, gamma);

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (ISNAN(xi)) {
      log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
    } else {
      log_f[i] = general.at(xi, alpha, beta, gamma);
    }
  }

  if (!log_p) {
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}

}  // namespace

// [[Rcpp::export(.dmhn_cpp, rng = false)]]
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

  // One cache slot per distinct triple in the recycling cycle, so that
  // log Psi is computed once per triple rather than once per element.
  mhn::ParamSlots<ParamCache> slots(na, nb, ng, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i % nx];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    log_f[i] = slots.at(i, a, b, g).log_density_at(xi, a, b, g);
    if ((i & 255) == 0) Rcpp::checkUserInterrupt();
  }

  if (!log_p) {
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}
