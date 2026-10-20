// CDF (pmhn) -- vectorized over q and (alpha, beta, gamma).
//
// Dispatch (shared CdfState in mhn_pmhn_state.h):
//   * gamma  = 0   -> sqrt-Gamma special case (Sun et al. 2023, Lemma 6a):
//                     F(q) = P(alpha/2, beta q^2), P the regularised lower
//                     incomplete gamma function.
//   * alpha  = 1   -> truncated normal special case (Sun et al. 2023,
//                     Lemma 6b).  Both are taken at those values exactly
//                     (see mhn_special_cases.h, is_sqrt_gamma_exact).
//   * default      -> each tail is computed directly on its own side of the
//                     peak of the density of log X, the positive root of
//                     2 beta x^2 - gamma x - alpha = 0, which lies at or
//                     above the median.  The upper tail P(X > q) is the
//                     Boost.Math quadrature of the density over [q, Inf),
//                     except on the log scale below the median (see
//                     below).  Above the peak the lower tail is one minus
//                     that upper tail, which is then at most one half, so
//                     nothing cancels.  Up to the peak, for alpha < 1 with
//                     gamma/sqrt(beta) <= 1e3 and for |gamma|/sqrt(beta) <= 20
//                     with alpha <= 100 (mhn_cdf_series.h), the lower tail is
//                     the Sun et al. (2023) Lemma 1b series, truncated at
//                     the Supplementary Lemma 10(d) K = max(K1, K2).  The
//                     series returns NaN when its alternating-sign
//                     cancellation guard fires (gamma < 0), when its last
//                     index 2K + 1 exceeds the ceiling of 1e6, and when
//                     beta q^2 is
//                     below the smallest normal double (see
//                     mhn_cdf_series.cpp); F is then, and outside those
//                     ranges always, taken by Boost.Math quadrature of the
//                     unnormalised density x^(alpha-1) exp(-beta x^2 +
//                     gamma x) over [0, q].  Below the median, where the
//                     upper-tail quadrature exceeds one half and F is at
//                     most one half, log P(X > q) for log.p = TRUE is about
//                     -F, and the log of that quadrature carries an absolute
//                     error set by the size of the terms it is assembled
//                     from; where that error would pass 1e-9 of -F it is
//                     taken as log1p(-F) instead
//                     (CdfState::log_upper_general).
//                     For gamma/sqrt(beta) > 1e3 both quadratures are used
//                     with gamma^2/(4 beta) removed from them and from Psi.
//   * alpha > 1e6 or gamma/sqrt(beta) > 1e6, or a shape or tilt m above
//     1e3 with m (1 + |log v*|) > 1e7 (v* the peak below; use_peak_frame)
//                  -> both tails are measured from
//                     the peak of the density of log X in exact offsets
//                     (mhn_cdf_peak.h), where the terms the forms above
//                     assemble the probability from would carry more than
//                     about 2e-9 of rounding, and where the law can be
//                     narrower than the spacing of the doubles.
//
// Recycling (one cached CdfState per distinct parameter triple, see
// mhn_param_slots.h) and the {q, alpha, beta, gamma} -> length max(...) rule
// mirror mhn_dmhn.cpp.

#include "mhn_cdf_integrate.h"
#include "mhn_cdf_series.h"
#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_pmhn_state.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_param_slots.h"
#include "mhn_warn.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

// Evaluate one element in the requested (lower.tail, log.p) form.
//
// The upper tail is taken from the state's own survival evaluation rather than
// derived here as 1 - F.  Deriving it loses the answer once F rounds to 1,
// which happens as soon as the survival probability drops below about 1e-16:
// at alpha = 2.5, beta = 1, gamma = 1 that is q = 7, where the true value is
// 5.8e-19 and 1 - F is exactly 0.  For the same reason log F is taken on the
// log scale throughout rather than as the log of F, which would be 0 once F
// rounds to 1 (log F at alpha = 2.5, gamma = 1, q = 10 is -9.7e-40) and -Inf
// once F underflows (at alpha = 2, gamma = 3, q = 1e-200 it is -924.96), and
// log P(X > q) likewise, which is about -F where F is small and is taken
// from log F below the median wherever the upper tail would not resolve it
// (CdfState::log_upper_general).
double eval_tail(const mhn::CdfState& state, double qi,
                 bool lower_tail, bool log_p) {
  // NA and NaN are distinct in the base R d/p/q contract, and
  // Rcpp::NumericVector::is_na does not separate them -- it is ISNAN -- so
  // each is tested with the R_Is* predicate that isolates it.
  if (R_IsNA(qi)) return NA_REAL;
  if (R_IsNaN(qi)) return R_NaN;
  if (state.kind == mhn::CDF_VANISHING) {
    return state.vanishing.p_one(qi, lower_tail, log_p);
  }
  // The special cases take the caller's flags themselves, through the same
  // per-element functions the all-scalar path loops over, so a recycled
  // parameter gives the same double as a scalar one.
  if (state.kind == mhn::CDF_SQRT_GAMMA) {
    return mhn::pmhn_sqrt_gamma_one(qi, state.sg_shape, state.sg_sqrt_beta,
                                    lower_tail, log_p);
  }
  if (state.kind == mhn::CDF_TRUNCATED_NORMAL) {
    return mhn::pmhn_truncated_normal_one(qi, state.tn, lower_tail, log_p);
  }
  if (!lower_tail) {
    if (log_p) {
      if (qi <= 0.0) return 0.0;
      if (qi == R_PosInf) return R_NegInf;
      const double log_Q = state.log_upper_general(qi);
      if (Rcpp::NumericVector::is_na(log_Q)) return NA_REAL;
      return log_Q;
    }
    const double log_Q = state.log_upper_tail(qi);
    if (Rcpp::NumericVector::is_na(log_Q)) return NA_REAL;
    return std::exp(log_Q);
  }
  if (log_p) {
    if (qi <= 0.0) return R_NegInf;
    if (qi == R_PosInf) return 0.0;
    const double log_F = state.log_cdf_general(qi);
    if (Rcpp::NumericVector::is_na(log_F)) return NA_REAL;
    return log_F;
  }
  const double F = state.cdf_linear(qi);
  if (Rcpp::NumericVector::is_na(F)) return NA_REAL;
  return F;
}

// Warns, once per call, of a probability that could not be evaluated: an NA
// for a q that is not NA or NaN, at valid parameters, means that the
// distribution function could not be formed for that triple (see CdfState),
// which the warning names.
// A pass of its own over the result, so that the loops that compute it carry
// no test.
void report_failure(const Rcpp::NumericVector& out,
                    const Rcpp::NumericVector& q,
                    const Rcpp::NumericVector& alpha,
                    const Rcpp::NumericVector& beta,
                    const Rcpp::NumericVector& gamma) {
  const R_xlen_t n = out.size();
  const double* value = out.begin();
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!(std::isnan(value[i]) && R_IsNA(value[i]))) continue;
    const double a = alpha[i % alpha.size()];
    const double b = beta[i % beta.size()];
    const double g = gamma[i % gamma.size()];
    if (!ISNAN(q[i % q.size()]) && mhn::params_valid(a, b, g)) {
      mhn::warn("pmhn: the distribution function could not be evaluated "
                "for alpha=%g, beta=%g, gamma=%g; returning NA", a, b, g);
      return;
    }
  }
}

}  // namespace

// [[Rcpp::export(.pmhn_cpp, rng = false)]]
Rcpp::NumericVector pmhn_cpp(Rcpp::NumericVector q,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool lower_tail,
                             bool log_p) {
  const R_xlen_t nq = q.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // A zero-length argument gives numeric(0), as in pnorm.
  if (nq == 0 || na == 0 || nb == 0 || ng == 0) return Rcpp::NumericVector(0);

  // Invalid parameters are handled element by element, as in dmhn_cpp.
  const bool all_valid = mhn::all_params_valid(alpha, beta, gamma);

  // Fast path: all parameters scalar.  The special cases go to the
  // per-element functions of mhn_special_cases.cpp, which take the user's
  // lower.tail / log.p themselves; the general case loops over eval_tail.
  if (all_valid && na == 1 && nb == 1 && ng == 1) {
    const double a = alpha[0], b = beta[0], g = gamma[0];
    // A subnormal alpha goes through the state, which takes the forms of
    // mhn_vanishing.h where they apply, before the special cases.
    if (!mhn::vanishing_candidate(a)) {
      if (mhn::is_sqrt_gamma_exact(g) && !mhn::sqrt_gamma_in_peak_frame(a)) {
        return mhn::pmhn_sqrt_gamma(q, a, b, lower_tail, log_p);
      }
      if (mhn::is_truncated_normal_exact(a)) {
        return mhn::pmhn_truncated_normal(q, b, g, lower_tail, log_p);
      }
    }
    mhn::CdfState state;
    state.recompute(a, b, g);
    Rcpp::NumericVector out(nq);
    for (R_xlen_t i = 0; i < nq; ++i) {
      out[i] = eval_tail(state, q[i], lower_tail, log_p);
      if ((i & 255) == 0) Rcpp::checkUserInterrupt();
    }
    report_failure(out, q, alpha, beta, gamma);
    return out;
  }

  const R_xlen_t n = std::max({nq, na, nb, ng});
  Rcpp::NumericVector out(n);

  // One cache slot per distinct triple in the recycling cycle: under
  // recycling consecutive elements rarely share a triple, so a single
  // most-recently-used slot would miss on every element.
  mhn::ParamSlots<mhn::CdfState> slots(na, nb, ng, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double qi = q[i % nq];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!all_valid && !mhn::params_valid(a, b, g)) {
      out[i] = mhn::invalid_result(qi, a, b, g);
      continue;
    }
    out[i] = eval_tail(slots.at(i, a, b, g), qi, lower_tail, log_p);
    if ((i & 255) == 0) Rcpp::checkUserInterrupt();
  }
  report_failure(out, q, alpha, beta, gamma);
  return out;
}

// Diagnostic hook for the tests and for accuracy audits: evaluates
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
