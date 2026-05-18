// CDF (pmhn) -- vectorized over q and (alpha, beta, gamma).
//
// Dispatch (shared CdfState in mhn_pmhn_state.h):
//   * gamma  ~ 0   -> sqrt-Gamma special case (Sun et al. 2023, Lemma 6a):
//                     F(q) = pgamma(q^2; alpha/2, scale = 1/beta).
//   * alpha  ~ 1   -> truncated normal special case (Sun et al. 2023,
//                     Lemma 6b).
//   * default      -> Sun et al. (2023) Lemma 1(b) series, truncated at
//                     the Lemma 10(d) K = max(K1, K2).  The series
//                     returns NaN when its double-precision
//                     alternating-sign cancellation guard fires (see
//                     mhn_cdf_series.cpp), in which case the Boost.Math
//                     integration of the unnormalised density (Sun
//                     et al. 2023 Lemma 11) is used as the runtime
//                     fallback.
//
// Recycling, ParamCache, and the {q, alpha, beta, gamma} -> length max(...)
// rule mirror mhn_dmhn.cpp.

#include "mhn_cdf_integrate.h"
#include "mhn_cdf_series.h"
#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_pmhn_state.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

// Translate F in [0, 1] to the requested (lower.tail, log.p) form.
double tail_log_transform(double F, bool lower_tail, bool log_p) {
  if (Rcpp::NumericVector::is_na(F)) return NA_REAL;
  if (log_p) {
    if (lower_tail) {
      if (F <= 0.0) return R_NegInf;
      if (F >= 1.0) return 0.0;
      return std::log(F);
    }
    if (F <= 0.0) return 0.0;
    if (F >= 1.0) return R_NegInf;
    return std::log1p(-F);
  }
  return lower_tail ? F : (1.0 - F);
}

}  // namespace

// [[Rcpp::export(.pmhn_cpp)]]
Rcpp::NumericVector pmhn_cpp(Rcpp::NumericVector q,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool lower_tail,
                             bool log_p) {
  mhn::check_params_vector(alpha, beta, gamma);

  const R_xlen_t nq = q.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  if (nq == 0) return Rcpp::NumericVector(0);

  // Fast path: all parameters scalar.  Dispatches special cases directly
  // to R::pgamma / R::pnorm with the user's flags intact, preserving
  // pnorm/pgamma tail precision.
  if (na == 1 && nb == 1 && ng == 1) {
    const double a = alpha[0], b = beta[0], g = gamma[0];
    if (mhn::is_sqrt_gamma(g)) {
      return mhn::pmhn_sqrt_gamma(q, a, b, lower_tail, log_p);
    }
    if (mhn::is_truncated_normal(a)) {
      return mhn::pmhn_truncated_normal(q, b, g, lower_tail, log_p);
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    Rcpp::NumericVector out(nq);
    for (R_xlen_t i = 0; i < nq; ++i) {
      const double F = state.cdf_linear(q[i]);
      out[i] = tail_log_transform(F, lower_tail, log_p);
    }
    return out;
  }

  const R_xlen_t n = std::max({nq, na, nb, ng});
  Rcpp::NumericVector out(n);

  mhn::CdfState state;
  double prev_a = std::numeric_limits<double>::quiet_NaN();
  double prev_b = std::numeric_limits<double>::quiet_NaN();
  double prev_g = std::numeric_limits<double>::quiet_NaN();
  bool primed = false;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double qi = q[i % nq];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!primed || a != prev_a || b != prev_b || g != prev_g) {
      state.recompute(a, b, g);
      prev_a = a; prev_b = b; prev_g = g;
      primed = true;
    }

    const double F = state.cdf_linear(qi);
    out[i] = tail_log_transform(F, lower_tail, log_p);
  }
  return out;
}

// Diagnostic hook for inst/audits/cdf_series_accuracy.R: evaluates
// the general-case CDF via a single, caller-specified code path,
// bypassing both the special-case shortcuts and the automatic
// series-to-integration fallback used by `pmhn()`.  Not exposed to
// package users; intended for benchmarks that compare the Sun et al.
// (2023) Lemma 1b series against the Boost.Math integration fallback
// at controlled (alpha, beta, gamma, q) points.
//
// method = "series"    -> mhn::log_cdf_series; returns NaN when its
//                         catastrophic-cancellation guard fires.
// method = "integrate" -> mhn::log_cdf_integrate (Boost.Math), with
//                         the log-prefactor added by the caller.
//
// Boundary handling matches the dispatcher: NA -> NA;
// q <= 0 -> 0; q == +Inf -> 1.
//
// [[Rcpp::export(.pmhn_force_cpp)]]
double pmhn_force_cpp(double q, double alpha, double beta, double gamma,
                      std::string method) {
  if (!(alpha > 0.0)) Rcpp::stop("alpha must be positive");
  if (!(beta > 0.0)) Rcpp::stop("beta must be positive");
  if (method != "series" && method != "integrate") {
    Rcpp::stop("method must be \"series\" or \"integrate\"");
  }
  if (Rcpp::NumericVector::is_na(q)) return NA_REAL;
  if (q <= 0.0) return 0.0;
  if (q == R_PosInf) return 1.0;

  const double tol = mhn::mhn_eps();
  const double log_psi =
      mhn::mhn_log_normalizing_const(alpha, beta, gamma, -1.0);

  double log_F;
  if (method == "series") {
    log_F = mhn::log_cdf_series(alpha, beta, gamma, q, log_psi, tol);
    if (!R_finite(log_F)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  } else {  // "integrate"
    const double log_prefactor =
        std::log(2.0) + 0.5 * alpha * std::log(beta) - log_psi;
    const double log_I = mhn::log_cdf_integrate(alpha, beta, gamma, q, tol);
    log_F = log_prefactor + log_I;
  }

  double F = std::exp(log_F);
  if (F < 0.0) F = 0.0;
  if (F > 1.0) F = 1.0;
  return F;
}
