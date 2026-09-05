// CDF series for the MHN distribution (Sun et al. 2023, Lemma 1b).
//
// After simplifying beta^(alpha/2) * beta^(-(alpha+i)/2) = beta^(-i/2), the
// Lemma 1b expression reduces to
//
//   F(x) = (1 / Psi[alpha/2, z]) * sum_{i=0}^{infty} z^i / i!
//                                  * Gamma(s_i) * P(s_i, y)
//
// where z = gamma / sqrt(beta), s_i = (alpha + i) / 2, y = beta x^2, and
// P(s, y) is the regularized lower incomplete gamma function
// (R::pgamma(y, s, scale=1, lower=TRUE)).
// (Equivalent to the paper's non-regularized form gamma(s_i, y) via the
//  identity gamma(s, y) = Gamma(s) * P(s, y).)
//
// Numerical strategy:
//   * Work in log space.  log T_i = i log|z| + lgamma(s_i)
//                                   + log P(s_i, y) - lgamma(i + 1).
//   * gamma >= 0 (z >= 0): every T_i is non-negative => single log_sum_exp.
//   * gamma <  0 (z <  0): T_i alternates sign; collect log|T_i| separately
//                          for even / odd i, then subtract in log space as
//                          log_S_pos + log1p(-exp(log_S_neg - log_S_pos)),
//                          under the cancellation guard described below.
//   * Truncation: the upper bound is taken from Sun et al. (2023)
//     Supplementary, Lemma 10(d).  Because the i-th CDF term satisfies
//     |T_{2k}| <= A(k) and |T_{2k+1}| <= |B(k)|, the same K = max(K1, K2)
//     that bounds Psi's truncation residual at epsilon also bounds the
//     Lemma 1b CDF residual at epsilon / Psi.  (Note: that per-term
//     bound is a package-internal derivation; it is not stated in Sun
//     et al. (2023) but follows immediately from comparing the
//     Lemma 1(b) CDF term log T_i = i log|z| + lgamma(s_i)
//     + log P(s_i, y) - lgamma(i+1) to the Psi series terms A(k) and
//     B(k) defined in Lemma 10.)  No patience heuristic, no fixed-
//     magnitude threshold -- the K is constructive and depends only on
//     (alpha, |z|, tol).
//
// Cancellation guard (gamma < 0):
//   The log_sum_exp accumulators for log_S_pos and log_S_neg each
//   carry relative precision ~ 2^-52 (= eps_d) per term, and over the
//   i_max terms of the sum that roundoff accumulates as a random walk,
//   ~ eps_d * sqrt(i_max).  Their difference S_pos - S_neg therefore
//   has relative error roughly
//
//      rel_err = 2 * eps_d * sqrt(i_max)
//                / (1 - exp(log_S_neg - log_S_pos)),
//
//   so to keep the final F within the user's tolerance tol_eff we
//   require 1 - exp(log_S_neg - log_S_pos)
//     >= 2 * eps_d * sqrt(i_max) / tol_eff.
//   When the inequality is violated, the routine returns NaN so the
//   dispatcher can fall back to the Boost.Math integration of the
//   unnormalised density.  This formalises Sun et al.'s qualitative
//   observation (Supplementary Section 1) that the Lemma 10 series
//   accumulates lgamma rounding errors that can exceed Psi for
//   gamma < 0, by quantifying "too many digits lost" in terms of
//   double precision's known floor.
//   inst/audits/series_breakdown_audit.R reproduces the measurement that
//   fixes the threshold, against an arbitrary-precision reference.

#include "mhn_cdf_series.h"
#include "mhn_constants.h"
#include "mhn_log_arith.h"
#include "mhn_psi.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mhn {

namespace {

// Lemma 10's q parameter chosen to match the existing psi_series
// implementation (mhn_psi_series.cpp).  Same q keeps the truncation
// rule consistent across the two routines that share K1/K2.
constexpr double Q_LEMMA10 = 0.5;

// 2 * 2^-52 -- twice the double-precision unit roundoff.  This is only
// the constant numerator of the cancellation-loss formula derived in
// the file header comment; the guard in log_cdf_series multiplies it
// by sqrt(i_max) / tol_eff to get the minimum allowable
// 1 - exp(log_S_neg - log_S_pos).
constexpr double DOUBLE_EPS_TIMES_TWO = 2.0 * 2.220446049250313e-16;

// Hard ceiling on the number of series terms.  The Lemma 10(d) truncation
// length grows like z^2, so a large |gamma| asks for a summation whose cost
// and memory are quadratic in the tilt: |z| = 1e4 alone would request ~5e7
// terms, several hundred megabytes of accumulator and tens of seconds.
// Quadrature evaluates the same probability in a bounded number of function
// calls, so past this ceiling it is the better method on every axis.
// Exceeding it returns the NaN sentinel, which routes the caller to
// quadrature exactly as a cancellation failure would.
constexpr long I_MAX_GUARD = 1000000L;

}  // namespace

double log_cdf_series(double alpha, double beta, double gamma,
                      double x, double log_psi, double tol) {
  if (!(x > 0.0)) return -std::numeric_limits<double>::infinity();
  if (gamma == 0.0) {
    // Caller should have dispatched to sqrt-Gamma special case.  Guard
    // against accidental entry: log T_0 = lgamma(alpha/2) + log P,
    // higher terms vanish (z = 0).
    const double y = beta * x * x;
    const double log_p0 = R::pgamma(y, alpha / 2.0, 1.0, 1, 1);
    return R::lgammafn(alpha / 2.0) + log_p0 - log_psi;
  }

  const double y = beta * x * x;
  const double abs_z = std::fabs(gamma) / std::sqrt(beta);
  const double z2 = abs_z * abs_z;
  const double log_abs_z = std::log(abs_z);
  const bool gamma_negative = (gamma < 0.0);
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();

  // Sun et al. (2023) Supplementary Lemma 10(d) truncation bound,
  // reused for the Lemma 1b CDF series (|T_{2k}| <= A(k),
  // |T_{2k+1}| <= |B(k)|).
  const double log_eps_quarter = std::log(tol_eff / 4.0);
  const long C1 = lemma10_C(alpha,        z2, Q_LEMMA10,
                            6.0 * Q_LEMMA10,  4.0 * Q_LEMMA10);
  const long C2 = lemma10_C(alpha + 1.0,  z2, Q_LEMMA10,
                            10.0 * Q_LEMMA10, 12.0 * Q_LEMMA10);
  const long K1 = lemma10_K(alpha / 2.0,        abs_z, z2, C1,
                            Q_LEMMA10, log_eps_quarter, /*is_A=*/true);
  // Both branches take a = alpha/2: the odd CDF term i = 2k+1 has
  // s_i = alpha/2 + k + 1/2, and lemma10_K's B branch adds that half itself.
  // Passing (alpha+1)/2 counted it twice and over-stated K2.
  const long K2 = lemma10_K(alpha / 2.0, abs_z, z2, C2,
                            Q_LEMMA10, log_eps_quarter, /*is_A=*/false);
  const long i_max = 2L * std::max(K1, K2) + 1L;
  if (i_max > I_MAX_GUARD) {
    return std::numeric_limits<double>::quiet_NaN();  // -> quadrature
  }

  std::vector<double> log_pos;
  std::vector<double> log_neg;
  log_pos.reserve(static_cast<std::size_t>(i_max + 1L));
  if (gamma_negative) log_neg.reserve(static_cast<std::size_t>(i_max / 2L + 1L));

  bool used_any = false;
  for (long i = 0; i <= i_max; ++i) {
    const double s_i = (alpha + static_cast<double>(i)) / 2.0;
    const double log_P = R::pgamma(y, s_i, /*scale=*/1.0,
                                   /*lower_tail=*/1, /*log_p=*/1);
    if (!R_finite(log_P) || log_P == R_NegInf) {
      // P(s, y) underflowed to zero in double precision; later terms
      // (with even larger s) will too, so subsequent contributions to
      // the sum are bounded by the Lemma 10 residual budget and may
      // be skipped without affecting the truncation error.
      continue;
    }
    const double log_T_i = static_cast<double>(i) * log_abs_z
                           + R::lgammafn(s_i)
                           + log_P
                           - R::lgammafn(static_cast<double>(i) + 1.0);
    if (!R_finite(log_T_i)) continue;
    used_any = true;
    if (gamma_negative && ((i & 1L) == 1L)) log_neg.push_back(log_T_i);
    else                                     log_pos.push_back(log_T_i);

    if ((i & 1023L) == 0L) Rcpp::checkUserInterrupt();
  }

  if (!used_any || log_pos.empty()) {
    return -std::numeric_limits<double>::infinity();
  }

  const double log_S_pos = log_sum_exp(log_pos);
  double log_sum;
  if (log_neg.empty()) {
    log_sum = log_S_pos;
  } else {
    const double log_S_neg = log_sum_exp(log_neg);
    // Double-precision cancellation guard.  Require
    //   1 - exp(log_S_neg - log_S_pos) >= 2 * 2^-52 * sqrt(i_max) / tol_eff
    // so the relative error in S_pos - S_neg stays within the user's
    // tolerance.  Testing only for total wipeout (log_S_neg >= log_S_pos)
    // is not enough: the series returns silently wrong values well before
    // that, as inst/audits/series_breakdown_audit.R demonstrates.
    //
    // The sqrt(i_max) is what each accumulator actually carries.  Charging one
    // unit of roundoff assumes a single rounding, but the two sums run over
    // i_max terms and their errors accumulate as a random walk.  Without it the
    // guard admitted results some twenty times worse than the tolerance asked
    // for: at (alpha, gamma) = (10, -100) the dispatcher returned 0.542194283
    // where an independent quadrature gives 0.542194428, a relative error of
    // 2.7e-7 against a requested 1.5e-8.  The package's own integration path
    // matches the reference there, so the cost of firing earlier is only that
    // more calls take it.
    const double min_one_minus_ratio =
        DOUBLE_EPS_TIMES_TWO * std::sqrt(static_cast<double>(i_max)) / tol_eff;
    const double one_minus_ratio = -std::expm1(log_S_neg - log_S_pos);
    if (!(one_minus_ratio > min_one_minus_ratio)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    log_sum = log_S_pos + std::log1p(-std::exp(log_S_neg - log_S_pos));
  }
  return log_sum - log_psi;
}

}  // namespace mhn
