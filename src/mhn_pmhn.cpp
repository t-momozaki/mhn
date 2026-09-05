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
#include "mhn_param_slots.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

// Evaluate one element in the requested (lower.tail, log.p) form.
//
// The upper tail is taken from the state's own survival evaluation rather than
// derived as 1 - F.  Deriving it loses the answer once F rounds to 1, which
// happens as soon as the survival probability drops below about 1e-16: at
// alpha = 2.5, beta = 1, gamma = 1 that is q = 8, where the true value is
// 5e-18 and 1 - F is exactly 0.
double eval_tail(const mhn::CdfState& state, double qi,
                 bool lower_tail, bool log_p) {
  // NA and NaN are distinct in the base R d/p/q contract, and
  // Rcpp::NumericVector::is_na does not separate them -- it is ISNAN.
  // qmhn already tested them apart; these two collapsed NaN to NA.
  if (R_IsNA(qi)) return NA_REAL;
  if (R_IsNaN(qi)) return R_NaN;
  if (!lower_tail) {
    const double log_Q = state.log_upper_tail(qi);
    if (Rcpp::NumericVector::is_na(log_Q)) return NA_REAL;
    return log_p ? log_Q : std::exp(log_Q);
  }
  const double F = state.cdf_linear(qi);
  if (Rcpp::NumericVector::is_na(F)) return NA_REAL;
  if (!log_p) return F;
  if (F <= 0.0) return R_NegInf;
  if (F >= 1.0) return 0.0;
  return std::log(F);
}

}  // namespace

// [[Rcpp::export(.pmhn_cpp, rng = false)]]
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
    if (mhn::is_sqrt_gamma(g, b)) {
      return mhn::pmhn_sqrt_gamma(q, a, b, lower_tail, log_p);
    }
    if (mhn::is_truncated_normal(a)) {
      return mhn::pmhn_truncated_normal(q, b, g, lower_tail, log_p);
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    Rcpp::NumericVector out(nq);
    for (R_xlen_t i = 0; i < nq; ++i) {
      out[i] = eval_tail(state, q[i], lower_tail, log_p);
      if ((i & 255) == 0) Rcpp::checkUserInterrupt();
    }
    return out;
  }

  const R_xlen_t n = std::max({nq, na, nb, ng});
  Rcpp::NumericVector out(n);

  // One cache slot per distinct triple in the recycling cycle: under
  // recycling consecutive elements rarely share a triple, so a single
  // most-recently-used slot missed on every element.
  mhn::ParamSlots<mhn::CdfState> slots(na, nb, ng, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double qi = q[i % nq];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    out[i] = eval_tail(slots.at(i, a, b, g), qi, lower_tail, log_p);
    if ((i & 255) == 0) Rcpp::checkUserInterrupt();
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
// [[Rcpp::export(.pmhn_force_cpp, rng = false)]]
double pmhn_force_cpp(double q, double alpha, double beta, double gamma,
                      std::string method) {
  if (!(alpha > 0.0)) Rcpp::stop("alpha must be positive");
  if (!(beta > 0.0)) Rcpp::stop("beta must be positive");
  if (method != "series" && method != "integrate") {
    Rcpp::stop("method must be \"series\" or \"integrate\"");
  }
  // NA and NaN are distinct in the base R d/p/q contract, and
  // Rcpp::NumericVector::is_na does not separate them -- it is ISNAN.
  if (R_IsNA(q)) return NA_REAL;
  if (R_IsNaN(q)) return R_NaN;
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
