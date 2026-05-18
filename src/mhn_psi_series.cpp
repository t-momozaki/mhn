// Series expansion for the Fox-Wright Psi function with gamma > 0
// (Sun et al. 2023 Supplementary, Lemma 10).
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

// Hard ceiling on the truncation length K returned by lemma10_K.  For
// realistic |z| the analytical bound is in the low thousands; values
// approaching this guard imply pathological parameters that should be
// routed to the numerical-integration path instead.
constexpr long K_MAX_GUARD = 1000000L;

}  // namespace

// Critical index C from Sun et al. (2023) Supplementary, Lemma 10(a, b).
// Declared in mhn_psi.h so the CDF series (mhn_cdf_series.cpp) can
// reuse the same Lemma 10 truncation analysis -- since |T_{2k}| <= A(k)
// and |T_{2k+1}| <= |B(k)| for the Lemma 1(b) CDF series, the K
// returned here also bounds the CDF series remainder.
//
// ERRATA: the discriminant uses alpha_adj * z^2, NOT alpha_adj * z.
// The Sun et al. (2023) main-text statement of Lemma 10 contains a
// typesetting error here; the supplementary proof uses z^2, which is
// the form implemented below.
long lemma10_C(double alpha_adj, double z2, double q,
               double c1, double c2) {
  const double denom = 8.0 * q;
  const double t = c1 - z2;
  const double disc = t * t - denom * (c2 - alpha_adj * z2);
  if (disc > 0.0) {
    const double cval = (-t + std::sqrt(disc)) / denom;
    return std::max<long>(static_cast<long>(std::ceil(cval)), 1L);
  }
  return 1L;
}

// Upper bound K for Lemma 10(d) truncation, derived from the geometric decay
// from Lemma 10(a, c).  K = C + ceil(max(0, (log T(C) - log eps_quarter) / -log(q))).
long lemma10_K(double a, double z, double z2, long C, double q,
               double log_eps_quarter, bool is_A) {
  double log_T_C;
  if (is_A) {
    // log A(C) = lgamma(a + C) + C * log(z^2) - lgamma(2C + 1)
    log_T_C = R::lgammafn(a + C)
              + static_cast<double>(C) * std::log(z2)
              - R::lgammafn(2.0 * C + 1.0);
  } else {
    // log B(C) = lgamma(a + 1/2 + C) + (2C + 1) * log(z) - lgamma(2C + 2)
    log_T_C = R::lgammafn(a + 0.5 + C)
              + (2.0 * C + 1.0) * std::log(z)
              - R::lgammafn(2.0 * C + 2.0);
  }
  long extra = 0L;
  if (log_T_C > log_eps_quarter) {
    extra = static_cast<long>(
        std::ceil((log_T_C - log_eps_quarter) / (-std::log(q))));
    if (extra < 0) extra = 0;
  }
  return C + extra;
}

namespace {

// Fill a contiguous slice of `log_terms` (length K + 1, starting at `start_idx`)
// with the cumulative log of the series term ratios, mirroring the R vectorized
// form via cumsum.
//   start_log    : log of the k = 0 term (lgamma(a) for A, lgamma(a+1/2) + log(z) for B)
//   ratio_offset : 0 for A-series (log(a + k - 1) ...), 0.5 for B-series
//   denom_shift  : 0 for A-series ((2k-1)(2k)), 1 for B-series ((2k)(2k+1))
// In vectorized form k_A = 1..K1, where the term is:
//   log_ratios_A(k) = log(a + k - 1) + log(z^2) - log(2k - 1) - log(2k)        (denom_shift=0)
//   log_ratios_B(k) = log(a + 0.5 + k - 1) + log(z^2) - log(2k) - log(2k+1)    (denom_shift=1)
void fill_log_terms(std::vector<double>& log_terms, std::size_t start_idx,
                    long K, double a, double log_z2,
                    double start_log, double ratio_offset, int denom_shift) {
  log_terms[start_idx] = start_log;
  double acc = start_log;
  for (long k = 1; k <= K; ++k) {
    const double dk = static_cast<double>(k);
    const double ratio = std::log(a + ratio_offset + dk - 1.0)
                         + log_z2
                         - std::log(2.0 * dk - 1.0 + denom_shift)
                         - std::log(2.0 * dk + denom_shift);
    acc += ratio;
    log_terms[start_idx + static_cast<std::size_t>(k)] = acc;
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

  // Critical indices from Sun et al. (2023) Lemma 10(a, b).  ERRATA:
  // discriminant uses alpha_adj * z^2 (see lemma10_C comment above).
  const long C1 = lemma10_C(alpha,        z2, q, /*c1=*/6.0 * q,  /*c2=*/4.0 * q);
  const long C2 = lemma10_C(alpha + 1.0,  z2, q, /*c1=*/10.0 * q, /*c2=*/12.0 * q);

  // Truncation lengths (Lemma 10 d).
  const double log_eps_quarter = std::log(tol / 4.0);
  const long K1 = lemma10_K(a, z, z2, C1, q, log_eps_quarter, /*is_A=*/true);
  const long K2 = lemma10_K(a, z, z2, C2, q, log_eps_quarter, /*is_A=*/false);

  if (K1 > K_MAX_GUARD || K2 > K_MAX_GUARD) {
    Rcpp::stop("psi_series: truncation length exceeds K_MAX_GUARD; "
               "consider tightening tol or using closed-form/integration paths");
  }

  // Single allocation holding both series:
  //   [0 .. K1]                  -> A-series (K1 + 1 entries)
  //   [K1 + 1 .. K1 + K2 + 1]    -> B-series (K2 + 1 entries)
  // log Psi = log( sum_k A(k) + sum_k B(k) ) = log_sum_exp(combined pool).
  // This fuses what was two log_sum_exp + one log_add_exp into a single
  // max + sum scan.
  const std::size_t n_A = static_cast<std::size_t>(K1) + 1;
  const std::size_t n_B = static_cast<std::size_t>(K2) + 1;
  std::vector<double> log_terms(n_A + n_B);

  // A-series: log A(0) = lgamma(a); ratio_offset = 0, denom_shift = 0.
  fill_log_terms(log_terms, /*start_idx=*/0, K1, a, log_z2,
                 R::lgammafn(a), 0.0, 0);
  // B-series: log B(0) = lgamma(a + 1/2) + log(z); ratio_offset = 0.5, denom_shift = 1.
  fill_log_terms(log_terms, /*start_idx=*/n_A, K2, a, log_z2,
                 R::lgammafn(a + 0.5) + std::log(z), 0.5, 1);

  return log_sum_exp(log_terms);
}

}  // namespace mhn
