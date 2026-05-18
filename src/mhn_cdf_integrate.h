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

}  // namespace mhn

#endif  // MHN_CDF_INTEGRATE_H
