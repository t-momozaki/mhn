// Numerical-integration fallback for the MHN CDF, used by the pmhn
// dispatcher when the Sun et al. (2023) Lemma 1b series is unsuitable.
// Built on Boost.Math quadrature with peak-normalized integrands.

#ifndef MHN_CDF_INTEGRATE_H
#define MHN_CDF_INTEGRATE_H

namespace mhn {

// Numerical-integration fallback for the MHN CDF, used when the Sun
// et al. (2023) Lemma 1b series is unsuitable (log_cdf_series returns
// NaN when the alternating-sign cancellation guard for gamma < 0
// fires).
//
// Computes log integral_0^x g(t) dt where
//   g(t) = t^(alpha-1) exp(-beta t^2 + gamma t)
// is the unnormalized MHN density kernel.  The caller adds
// log(2) + (alpha/2) log(beta) - log Psi to recover log F(x).
//
// Implementation uses peak normalization to keep the magnitudes near unity
// and switches between Gauss-Kronrod (alpha >= 1, smooth integrand) and
// tanh-sinh (alpha < 1, x^(alpha-1) endpoint singularity at t = 0).
//
// `tol` is the relative tolerance; pass a non-positive value to use
// sqrt(.Machine$double.eps).
double log_cdf_integrate(double alpha, double beta, double gamma,
                         double x, double tol);

// The same integral with gamma^2/(4 beta) removed, for the caller that pairs it
// with psi_integrate_shifted.  For a large positive tilt log Psi and this
// integral are each about gamma^2/(4 beta) -- 2.5e13 at gamma = 1e7 -- while
// their difference, log F, is of order one, so assembling F from the unshifted
// pair loses it.  Removing the term from both analytically leaves nothing large.
double log_cdf_integrate_shifted(double alpha, double beta, double gamma,
                                 double x, double tol);

// The complement: log integral_x^inf g(t) dt, over the same kernel and with
// the same caller-supplied prefactor, giving log(1 - F(x)).
//
// Deriving the upper tail as 1 - F loses it entirely once F rounds to 1, which
// happens as soon as the survival probability falls below about 1e-16 -- at
// alpha = 2.5, beta = 1, gamma = 1 that is already q = 8.  Integrating the
// other way keeps the full range, exactly as base R's distribution functions
// do with their lower.tail argument.
double log_ccdf_integrate(double alpha, double beta, double gamma,
                          double x, double tol);

// The complement with gamma^2/(4 beta) removed, to pair with the shifted
// normalising constant exactly as log_cdf_integrate_shifted does.
double log_ccdf_integrate_shifted(double alpha, double beta, double gamma,
                                  double x, double tol);

}  // namespace mhn

#endif  // MHN_CDF_INTEGRATE_H
