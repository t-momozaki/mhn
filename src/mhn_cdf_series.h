// CDF of MHN(alpha, beta, gamma) evaluated by the Sun et al. (2023)
// Lemma 1b series.  See log_cdf_series() below for the truncation /
// cancellation strategy.

#ifndef MHN_CDF_SERIES_H
#define MHN_CDF_SERIES_H

namespace mhn {

// Sun et al. (2023) Lemma 1b series-based CDF for MHN(alpha, beta, gamma).
//
// Returns log F(x | alpha, beta, gamma), given a precomputed
//   log_psi = log Psi[alpha/2, gamma/sqrt(beta)].
// Caller must guarantee x > 0 (q <= 0 and q = +Inf are handled at the
// dispatcher).
//
// `tol` is the relative truncation tolerance for the series; pass a
// non-positive value to use sqrt(.Machine$double.eps) (mhn_eps()).
//
// Returns NaN to signal that the series could not produce a usable
// result (e.g. catastrophic cancellation when gamma < 0); the caller
// should then fall back to log_cdf_integrate().
double log_cdf_series(double alpha, double beta, double gamma,
                      double x, double log_psi, double tol);

}  // namespace mhn

#endif  // MHN_CDF_SERIES_H
