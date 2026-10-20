// Series expansion for the Fox-Wright Psi function with gamma > 0
// (Sun et al. 2023, Supplementary Lemma 10).
//
//   Psi[alpha/2, z] = sum_k A(k) + sum_k B(k),  z = gamma / sqrt(beta) > 0
//   A(k) = Gamma(alpha/2 + k) z^(2k) / (2k)!
//   B(k) = Gamma((alpha+1)/2 + k) z^(2k+1) / (2k+1)!
//
// Lemma 10 splits each sum at a critical index (parts a, b) past which
// term ratios drop below a chosen q in (0, 1), then bounds the
// truncation length K (part d) so the remainder is below the requested
// tolerance.

#include "mhn_psi.h"
#include "mhn_log_arith.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace mhn {

namespace {

// Hard ceiling on the truncation length K returned by lemma10_K.  The
// dispatcher sends only 0 < z <= 20 with alpha <= 100 here, where K is at most
// about 700.  K grows like z^2, and by about log Gamma(alpha/2) / log 2 for a
// large shape (see psi_series), so a direct call passes the ceiling at z = 20
// near alpha = 1.3e5.  Past it the series declines, and the caller turns to
// the numerical-integration path instead.
constexpr long K_MAX_GUARD = 1000000L;

}  // namespace

// Critical index C from Sun et al. (2023) Supplementary Lemma 10(a, b).
// Declared in mhn_psi.h so the CDF series (mhn_cdf_series.cpp) can
// reuse the same Lemma 10 truncation analysis -- since |T_{2k}| <= A(k)
// and |T_{2k+1}| <= |B(k)| for the Lemma 1b CDF series, the K
// returned here also bounds the CDF series remainder.
//
// ERRATA: the discriminant uses alpha_adj * z^2, NOT alpha_adj * z.
// The statement of Supplementary Lemma 10 of Sun et al. (2023) contains a
// typesetting error here; its proof (Supplementary Section 2.2) uses z^2,
// which is the form implemented below.
double lemma10_C(double alpha_adj, double z2, double q,
                 double c1, double c2) {
  const double denom = 8.0 * q;
  const double t = c1 - z2;
  const double disc = t * t - denom * (c2 - alpha_adj * z2);
  if (disc > 0.0) {
    const double cval = (-t + std::sqrt(disc)) / denom;
    return std::max(std::ceil(cval), 1.0);
  }
  return 1.0;
}

namespace {

// log A(C) (is_A) or log B(C), the term at the index C of Lemma 10(a, b).
double lemma10_log_term(double a, double z, double z2, double C, bool is_A) {
  if (is_A) {
    // log A(C) = lgamma(a + C) + C * log(z^2) - lgamma(2C + 1)
    return R::lgammafn(a + C)
           + C * std::log(z2)
           - R::lgammafn(2.0 * C + 1.0);
  }
  // log B(C) = lgamma(a + 1/2 + C) + (2C + 1) * log(z) - lgamma(2C + 2)
  return R::lgammafn(a + 0.5 + C)
         + (2.0 * C + 1.0) * std::log(z)
         - R::lgammafn(2.0 * C + 2.0);
}

// K = C + ceil(max(0, (log T(C) - log_threshold) / -log(q))), from the term
// log T(C) at C: past C the terms fall by at least the factor q, so from K on
// each is below exp(log_threshold).
double lemma10_K_from_term(double C, double log_T_C, double q,
                           double log_threshold) {
  double extra = 0.0;
  if (log_T_C > log_threshold) {
    extra = std::ceil((log_T_C - log_threshold) / (-std::log(q)));
    if (!(extra > 0.0)) extra = 0.0;
  }
  return C + extra;
}

}  // namespace

// Upper bound K for Lemma 10(d) truncation, derived from the geometric decay
// of Lemma 10(a, b), the term ratio below q past C:
//   K = C + ceil(max(0, (log T(C) - log eps_quarter) / -log(q))).
double lemma10_K(double a, double z, double z2, double C, double q,
                 double log_eps_quarter, bool is_A) {
  return lemma10_K_from_term(C, lemma10_log_term(a, z, z2, C, is_A), q,
                             log_eps_quarter);
}

namespace {

// Append to `log_terms` log A(k) or log B(k), k = 0..K: start_log followed
// by its running sum with the log term ratios, and then further terms for as
// long as the last one appended exceeds exp(log_negligible).
//   start_log    : log of the k = 0 term (lgamma(a) for A, lgamma(a+1/2) + log(z) for B)
//   ratio_offset : 0 for A-series (log(a + k - 1) ...), 0.5 for B-series
//   denom_shift  : 0 for A-series ((2k-1)(2k)), 1 for B-series ((2k)(2k+1))
// For k >= 1 the log ratio is:
//   log_ratios_A(k) = log(a + k - 1) + log(z^2) - log(2k - 1) - log(2k)        (denom_shift=0)
//   log_ratios_B(k) = log(a + 0.5 + k - 1) + log(z^2) - log(2k) - log(2k+1)    (denom_shift=1)
void append_log_terms(std::vector<double>& log_terms, long K, double a,
                      double log_z2, double start_log, double ratio_offset,
                      int denom_shift, double log_negligible) {
  log_terms.push_back(start_log);
  double acc = start_log;
  for (long k = 1; k <= K || (acc > log_negligible && k <= K_MAX_GUARD);
       ++k) {
    const double dk = static_cast<double>(k);
    // dk - 1 is exact, so the factor a + offset + (k - 1) never passes through
    // a + k.  Grouped as (a + k) - 1 it would lose the digits of a below the
    // rounding of k: at k = 1 of the A series (offset 0) the factor is a
    // itself, and (a + 1) - 1 keeps only the leading bits of a tiny shape --
    // none at all once a < 1.1e-16 -- and every later A term inherits the
    // error.
    const double ratio = std::log(a + ratio_offset + (dk - 1.0))
                         + log_z2
                         - std::log(2.0 * dk - 1.0 + denom_shift)
                         - std::log(2.0 * dk + denom_shift);
    acc += ratio;
    log_terms.push_back(acc);
    if (k % 1024 == 0) Rcpp::checkUserInterrupt();
  }
}

}  // namespace

double psi_series(double alpha, double beta, double gamma, double tol) {
  const double a = alpha / 2.0;
  const double z = gamma / std::sqrt(beta);
  const double z2 = z * z;
  const double log_z2 = std::log(z2);
  const double q = 0.5;

  // Critical indices from Sun et al. (2023) Supplementary Lemma 10(a, b).
  // ERRATA:
  // discriminant uses alpha_adj * z^2 (see lemma10_C comment above).
  const double C1 = lemma10_C(alpha,       z2, q, /*c1=*/6.0 * q,  /*c2=*/4.0 * q);
  const double C2 = lemma10_C(alpha + 1.0, z2, q, /*c1=*/10.0 * q, /*c2=*/12.0 * q);

  // Truncation lengths (Lemma 10(d)).
  const double log_eps_quarter = std::log(tol / 4.0);
  const double log_A_C1 = lemma10_log_term(a, z, z2, C1, /*is_A=*/true);
  const double log_B_C2 = lemma10_log_term(a, z, z2, C2, /*is_A=*/false);
  double K1 = lemma10_K_from_term(C1, log_A_C1, q, log_eps_quarter);
  double K2 = lemma10_K_from_term(C2, log_B_C2, q, log_eps_quarter);

  if (!(K1 <= K_MAX_GUARD && K2 <= K_MAX_GUARD)) {
    // The Lemma 10(d) truncation length grows like z^2, and since it bounds
    // the absolute truncation error it also grows with the size of the terms,
    // by about log Gamma(alpha/2) / log 2 for a large shape.  Past the ceiling
    // the summation is neither affordable nor necessary: quadrature evaluates
    // the same quantity in a bounded number of function calls.  Report the
    // refusal so the caller can hand over, rather than raising.
    return std::numeric_limits<double>::quiet_NaN();
  }

  // The same bound taken relative to the terms kept.  All terms are positive,
  // and past C the terms of each series fall at least by the factor q = 1/2
  // (Lemma 10(a, b)).  So once A(k) is below 2^-56 of the larger of A(0) and
  // A(C1), every later A term is too, and each is below a quarter of the
  // spacing of the doubles at the running sum it is added to (the sum below
  // takes the A terms first, then the B terms): adding it leaves that sum
  // unchanged to the bit, and so does leaving it out.  The same holds for the
  // B terms against the larger of A(0), A(C1), B(0) and B(C2).  The absolute
  // bound tol / 4 of Lemma 10(d) grows with the size of the terms, by about
  // log Gamma(alpha/2) / log 2 terms, and this one does not: at alpha = 1000
  // it keeps a few dozen terms past C where the absolute bound keeps
  // thousands.  The value returned is the same double either way.
  // The terms at 0 and at C are formed here once, and serve both this bound
  // and the sums below.
  const double margin = 56.0 * M_LN2;
  const double log_A0 = R::lgammafn(a);
  const double log_B0 = R::lgammafn(a + 0.5) + std::log(z);
  const double log_A_lb = std::max(log_A0, log_A_C1);
  const double log_B_lb = std::max({log_A_lb, log_B0, log_B_C2});
  if (log_A_lb - margin > log_eps_quarter) {
    K1 = std::min(K1, lemma10_K_from_term(C1, log_A_C1, q, log_A_lb - margin));
  }
  if (log_B_lb - margin > log_eps_quarter) {
    K2 = std::min(K2, lemma10_K_from_term(C2, log_B_C2, q, log_B_lb - margin));
  }

  // The absolute bound can also stop a series too early for the double it
  // returns.  When the terms are of order one and fall well below the
  // factor q per step, K can come out barely past C, and the first term left
  // out, below tol / 4 but not below the spacing of the doubles at the sum,
  // is lost: at alpha = 2.3, z = 0.0024 the Lemma 10(d) length keeps B(0)
  // and B(1), just below tol / 4, and leaves out B(2), 3e-15 of Psi, 14
  // units in the last place of log Psi.  So each series is continued past K
  // while its last term is above 2^-56 of the lower bound of its sum used
  // above.  K is at least C, so every later term is smaller than that one by
  // at least a factor q per step, and all of them together are no larger than
  // it: the part left out is below 2^-56 of the sum.  Where K already reaches
  // that point, as it does unless K lies close to C, no term is added.
  //
  // One buffer holds both series, the A terms and then the B terms,
  // log Psi = log( sum_k A(k) + sum_k B(k) ) = log_sum_exp(combined pool).
  // Pooling both series into one buffer lets a single max-and-sum scan do
  // the work of two separate log-sum-exps and a combining step.
  // Both lengths are at most the ceiling, 1e6, and fit any integer type.
  const long K1l = static_cast<long>(K1);
  const long K2l = static_cast<long>(K2);
  std::vector<double> log_terms;
  log_terms.reserve(static_cast<std::size_t>(K1l + K2l) + 8);

  // A-series: log A(0) = lgamma(a); ratio_offset = 0, denom_shift = 0.
  append_log_terms(log_terms, K1l, a, log_z2, log_A0, 0.0, 0,
                   log_A_lb - margin);
  // B-series: log B(0) = lgamma(a + 1/2) + log(z); ratio_offset = 0.5, denom_shift = 1.
  append_log_terms(log_terms, K2l, a, log_z2, log_B0, 0.5, 1,
                   log_B_lb - margin);

  return log_sum_exp(log_terms);
}

}  // namespace mhn
