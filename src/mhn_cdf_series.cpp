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
//     Supplementary Lemma 10(d).  Because the i-th CDF term satisfies
//     |T_{2k}| <= A(k) and |T_{2k+1}| <= |B(k)|, the same K = max(K1, K2)
//     that bounds Psi's truncation residual at epsilon also bounds the
//     Lemma 1b CDF residual at epsilon / Psi.  (Note: that per-term
//     bound is a package-internal derivation; it is not stated in Sun
//     et al. (2023) but follows immediately from comparing the
//     Lemma 1b CDF term log T_i = i log|z| + lgamma(s_i)
//     + log P(s_i, y) - lgamma(i+1) to the Psi series terms A(k) and
//     B(k) defined in Lemma 10.)  No patience heuristic, no fixed-
//     magnitude threshold -- the K is constructive and depends only on
//     (alpha, |z|, tol).  For gamma > 0, whose terms are all positive, the
//     sum also stops where the terms fall below machine precision relative
//     to Psi, if that comes first; the truncation residual in F is then below
//     the larger of epsilon / Psi and 2^-53 F (see log_cdf_series).
//
// Cancellation guard (gamma < 0):
//   The i-th term is |T_i| = exp(l_i), where l_i is the sum of the four
//   logarithms i log|z|, lgamma(s_i), log P(s_i, y) and -lgamma(i + 1).  Each
//   is computed to within a few units of rounding of its own magnitude, so
//   l_i carries an absolute error of about eps_d L_i, eps_d = 2^-52 and L_i
//   the sum of the four magnitudes, and |T_i| that relative error.  Summed,
//   the error of S_pos - S_neg is about eps_d L_bar (S_pos + S_neg), where
//   L_bar = sum |T_i| L_i / sum |T_i| over both signs, and relative to the
//   difference it is at most
//
//      rel_err = 2 * eps_d * L_bar / (1 - exp(log_S_neg - log_S_pos)).
//
//   The log-sum-exp accumulators add a rounding of about eps_d |log S_pos|
//   relative, and |log S_pos| is at most about L_bar, so they change this by
//   no more than a small factor.  To keep F within the working tolerance
//   tol_eff the series is therefore used only when
//     1 - exp(log_S_neg - log_S_pos) >= 2 * eps_d * L_bar / tol_eff;
//   otherwise the routine returns NaN and the dispatcher falls back to the
//   Boost.Math integration of the unnormalised density.
//
//   L_bar grows with alpha through lgamma(s_i), and with |z| through
//   i log|z| and lgamma(i + 1), each taken at the largest terms.  These lie
//   near the largest term of the series of Psi at |z|, whose term ratio is
//   about |z| sqrt(s_i) / (i + 1): at i = (z^2 + sqrt(z^4 + 8 alpha z^2)) / 4,
//   about z^2 / 2 when z^2 >> alpha and about |z| sqrt(alpha / 2) when
//   alpha >> z^2 -- or below it when beta x^2 is small, since P(s_i, beta x^2)
//   then cuts the terms off earlier.  A bound that charges only the
//   number of terms misses the first: at alpha = 2000, gamma = -0.2 it lets
//   through series values in error by 8 times the tolerance.  This quantifies
//   the qualitative observation of Sun et al. (2023, Supplementary
//   Section 1) that for gamma < 0 the lgamma rounding in the Lemma 10 series
//   becomes significant against Psi.

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

// 2 * 2^-52, the constant factor of the cancellation-loss bound derived in
// the file header comment; the guard in log_cdf_series multiplies it by
// L_bar / tol_eff to get the minimum allowable
// 1 - exp(log_S_neg - log_S_pos).
constexpr double DOUBLE_EPS_TIMES_TWO = 2.0 * 2.220446049250313e-16;

// Hard ceiling on the last index of the series.  The series runs to index
// 2K + 1, which grows like 1.44 z^2 (K = max(K1, K2), the Lemma 10(d)
// length), so a large |gamma| asks for a summation whose cost and memory are
// quadratic in the tilt: |z| = 1e4 alone would request about 1.4e8 terms,
// over a gigabyte of accumulator and tens of seconds.
// Quadrature evaluates the same probability in a bounded number of function
// calls, so past this ceiling it is the better method on every axis.
// Exceeding it returns the NaN sentinel, which routes the caller to
// quadrature exactly as a cancellation failure would.
constexpr long I_MAX_GUARD = 1000000L;

// The length max(K1, K2) of Lemma 10(d), in terms of which the series is
// summed up to index 2 max(K1, K2) + 1, for z = |gamma| / sqrt(beta) and the
// logarithm of a quarter of the tolerance; NaN where that index would pass
// the ceiling I_MAX_GUARD, and the series is refused whatever x.  C1 and C2
// receive the critical indices of Lemma 10(a, b).
//
// Two lower bounds on the last index 2K + 1 decide, before K is formed, that
// the ceiling would refuse the series.
//
//  * z^2: the Lemma 10(a) index C_1 is at least z^2/2 - 2.  It is the
//    positive root (z^2 - 3 + sqrt((z^2 - 3)^2 + 4 alpha z^2 - 8)) / 4 of
//    its quadratic (the form of the proof, Supplementary Section 2.2.1,
//    with q = 1/2; the statement of Lemma 10(a) prints alpha x for
//    alpha x^2 -- see ERRATA in mhn_psi_series.cpp), rounded up, and so
//    at least
//    (z^2 - 3 + sqrt((z^2 - 3)^2 - 8)) / 4, which is z^2/2 - 2 or more once
//    z^2 >= 6; below that C_1 >= 1 exceeds z^2/2 - 2.  Hence
//    i_max = 2 max(K1, K2) + 1 is at least z^2 - 3.
//
//  * alpha: for 1 <= k < C_1 the ratio A(k+1)/A(k) =
//    (alpha/2 + k) z^2 / ((2k + 1)(2k + 2)) exceeds q = 1/2, since C_1 is
//    the larger root of the quadratic on which it equals q, and when
//    C_1 > 1 the smaller root lies below 1/2.  Hence A(C_1) >=
//    A(1) 2^-(C_1 - 1), and the length K_1 = C_1 + ceil(log2(4 A(C_1) /
//    tol)) that lemma10_K takes (an upper bound on the K_1 of
//    Lemma 10(d)) is at least 1 + log2(4 A(1) / tol), with A(1) =
//    Gamma(alpha/2 + 1) z^2 / 2.  This grows like (alpha/2) log2(alpha/2)
//    and exceeds the ceiling at every z once alpha is above about 7.3e4.
//
// The lengths are formed as a long, which is 32 bits on some platforms, and
// left to grow they overflow it: the second term of K_1 reaches 2^31 at
// alpha of about 1.7e8.  With both bounds met, alpha is below about 7.5e4
// (unless z is 0, when every length is 1) and z^2 below 2e6, and the
// lengths stay below about 2.2e6.
double series_length(double alpha, double abs_z, double log_eps_quarter,
                     double* C1_out, double* C2_out) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double z2 = abs_z * abs_z;
  if (!(z2 <= 2.0 * static_cast<double>(I_MAX_GUARD))) return nan;
  const double log_A1 =
      R::lgammafn(alpha / 2.0 + 1.0) + 2.0 * std::log(abs_z) - M_LN2;
  if ((log_A1 - log_eps_quarter) / M_LN2 >
      0.5 * static_cast<double>(I_MAX_GUARD)) {
    return nan;
  }

  // Sun et al. (2023) Supplementary Lemma 10(d) truncation bound,
  // reused for the Lemma 1b CDF series (|T_{2k}| <= A(k),
  // |T_{2k+1}| <= |B(k)|).
  const double C1 = lemma10_C(alpha,        z2, Q_LEMMA10,
                              6.0 * Q_LEMMA10,  4.0 * Q_LEMMA10);
  const double C2 = lemma10_C(alpha + 1.0,  z2, Q_LEMMA10,
                              10.0 * Q_LEMMA10, 12.0 * Q_LEMMA10);
  const double K1 = lemma10_K(alpha / 2.0,        abs_z, z2, C1,
                              Q_LEMMA10, log_eps_quarter, /*is_A=*/true);
  // Both branches take a = alpha/2: the odd CDF term i = 2k+1 has
  // s_i = alpha/2 + k + 1/2, and lemma10_K's B branch adds that half itself,
  // so passing (alpha+1)/2 would count it twice and over-state K2.
  const double K2 = lemma10_K(alpha / 2.0, abs_z, z2, C2,
                              Q_LEMMA10, log_eps_quarter, /*is_A=*/false);
  // i_max = 2 max(K1, K2) + 1 exceeds the ceiling exactly when max(K1, K2)
  // exceeds (I_MAX_GUARD - 1) / 2.  The lengths are integer-valued doubles
  // (see lemma10_C), tested before any is converted, so that none can
  // overflow an integer type on a platform where long has 32 bits.
  const double K_max = std::max(K1, K2);
  if (!(K_max <= static_cast<double>((I_MAX_GUARD - 1L) / 2L))) return nan;
  if (C1_out != nullptr) *C1_out = C1;
  if (C2_out != nullptr) *C2_out = C2;
  return K_max;
}

}  // namespace

bool cdf_series_fits(double alpha, double beta, double gamma, double tol) {
  if (gamma == 0.0) return true;
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();
  const double abs_z = std::fabs(gamma) / std::sqrt(beta);
  return !std::isnan(series_length(alpha, abs_z, std::log(tol_eff / 4.0),
                                   nullptr, nullptr));
}

double log_cdf_series(double alpha, double beta, double gamma,
                      double x, double log_psi, double tol) {
  if (!(x > 0.0)) return -std::numeric_limits<double>::infinity();
  const double y = beta * x * x;
  // Below the smallest normal double y is subnormal and carries fewer
  // significant bits the smaller it is (a relative 5e-4 at y = 1e-320), and
  // P(s_i, y), which behaves like y^s_i, inherits that.  Such an x lies far
  // inside the stretch the quadrature takes in closed form.
  if (!(y >= std::numeric_limits<double>::min())) {
    return std::numeric_limits<double>::quiet_NaN();  // -> quadrature
  }
  if (gamma == 0.0) {
    // Caller should have dispatched to sqrt-Gamma special case.  Guard
    // against accidental entry: log T_0 = lgamma(alpha/2) + log P,
    // higher terms vanish (z = 0).
    const double log_p0 = R::pgamma(y, alpha / 2.0, 1.0, 1, 1);
    return R::lgammafn(alpha / 2.0) + log_p0 - log_psi;
  }
  const double abs_z = std::fabs(gamma) / std::sqrt(beta);
  const double z2 = abs_z * abs_z;
  const double log_abs_z = std::log(abs_z);
  const bool gamma_negative = (gamma < 0.0);
  const double tol_eff = (tol > 0.0) ? tol : mhn_eps();
  const double log_eps_quarter = std::log(tol_eff / 4.0);
  double C1, C2;
  const double K_max = series_length(alpha, abs_z, log_eps_quarter, &C1, &C2);
  if (std::isnan(K_max)) {
    return std::numeric_limits<double>::quiet_NaN();  // -> quadrature
  }
  const long i_max = 2L * static_cast<long>(K_max) + 1L;

  // For gamma > 0 the sum can stop earlier without changing its value beyond
  // rounding.  Every term is positive, and P(s_i, y) falls as i grows, so with
  // I the last index kept, the terms dropped are at most P(s_I, y) times the
  // corresponding terms A(k), B(k) of Psi, while those kept are at least
  // P(s_I, y) times theirs.  Once A(k) and B(k) fall below 2^-52 Psi / 4 the
  // dropped part is therefore below about 2^-53 of the part kept, however
  // small F is.  Truncating there as well as at tol / 4 matters once Psi is
  // large -- at alpha = 1000, gamma = 1 it cuts the sum from about 7700 terms
  // to about 130 -- and never lengthens it.  The length ceiling above is still
  // applied to the Lemma 10(d) length, so the hand-off to quadrature does not
  // move.
  long i_last = i_max;
  if (!gamma_negative) {
    const double log_rel_quarter =
        std::log(std::numeric_limits<double>::epsilon() / 4.0) + log_psi;
    if (log_rel_quarter > log_eps_quarter) {
      // With the same C and a larger threshold these lengths are no larger
      // than K1 and K2, which the ceiling has already bounded.
      const double K1_rel = lemma10_K(alpha / 2.0, abs_z, z2, C1, Q_LEMMA10,
                                      log_rel_quarter, /*is_A=*/true);
      const double K2_rel = lemma10_K(alpha / 2.0, abs_z, z2, C2, Q_LEMMA10,
                                      log_rel_quarter, /*is_A=*/false);
      i_last = std::min(
          i_max, 2L * static_cast<long>(std::max(K1_rel, K2_rel)) + 1L);
    }
  }

  std::vector<double> log_pos;
  std::vector<double> log_neg;
  log_pos.reserve(static_cast<std::size_t>(i_last + 1L));
  if (gamma_negative) log_neg.reserve(static_cast<std::size_t>(i_max / 2L + 1L));

  // For gamma < 0, the |T_i|-weighted mean L_bar of the rounding magnitudes
  // L_i that the cancellation guard charges (see the header).  The weights
  // exp(log T_i - w_max) are taken relative to the largest log T_i seen so
  // far, and the two sums are rescaled whenever that maximum moves, so no
  // weight overflows.
  double w_max = -std::numeric_limits<double>::infinity();
  double w_sum = 0.0;
  double wl_sum = 0.0;

  bool used_any = false;
  for (long i = 0; i <= i_last; ++i) {
    const double s_i = (alpha + static_cast<double>(i)) / 2.0;
    const double log_P = R::pgamma(y, s_i, /*scale=*/1.0,
                                   /*lower_tail=*/1, /*log_p=*/1);
    if (!std::isfinite(log_P)) {
      // P(s, y) underflowed to zero in double precision; later terms
      // (with even larger s) will too, so subsequent contributions to
      // the sum are bounded by the Lemma 10 residual budget and may
      // be skipped without affecting the truncation error.
      continue;
    }
    const double log_z_part = static_cast<double>(i) * log_abs_z;
    const double lgamma_s = R::lgammafn(s_i);
    const double lgamma_i = R::lgammafn(static_cast<double>(i) + 1.0);
    const double log_T_i = log_z_part + lgamma_s + log_P - lgamma_i;
    if (!std::isfinite(log_T_i)) continue;
    used_any = true;
    if (gamma_negative && ((i & 1L) == 1L)) log_neg.push_back(log_T_i);
    else                                     log_pos.push_back(log_T_i);

    if (gamma_negative) {
      const double L_i = std::fabs(log_z_part) + std::fabs(lgamma_s)
                         + std::fabs(log_P) + std::fabs(lgamma_i);
      if (log_T_i > w_max) {
        const double shrink = std::exp(w_max - log_T_i);  // 0 on the first term
        w_sum = w_sum * shrink + 1.0;
        wl_sum = wl_sum * shrink + L_i;
        w_max = log_T_i;
      } else {
        const double w = std::exp(log_T_i - w_max);
        w_sum += w;
        wl_sum += w * L_i;
      }
    }

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
    // Double-precision cancellation guard (see the header).  Require
    //   1 - exp(log_S_neg - log_S_pos) >= 2 * 2^-52 * L_bar / tol_eff
    // so that the relative error of S_pos - S_neg stays within the working
    // tolerance.  Testing only for total wipeout (log_S_neg >= log_S_pos)
    // is not enough: the series returns silently wrong values well before
    // that.  Where the guard fires, the quadrature fallback answers instead.
    const double L_bar = wl_sum / w_sum;
    const double min_one_minus_ratio =
        DOUBLE_EPS_TIMES_TWO * L_bar / tol_eff;
    const double one_minus_ratio = -std::expm1(log_S_neg - log_S_pos);
    if (!(one_minus_ratio > min_one_minus_ratio)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    log_sum = log_S_pos + std::log1p(-std::exp(log_S_neg - log_S_pos));
  }
  return log_sum - log_psi;
}

}  // namespace mhn
