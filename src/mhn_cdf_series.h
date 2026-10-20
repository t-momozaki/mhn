// CDF of MHN(alpha, beta, gamma) evaluated by the Sun et al. (2023)
// Lemma 1b series.  See log_cdf_series() below for the truncation /
// cancellation strategy.

#ifndef MHN_CDF_SERIES_H
#define MHN_CDF_SERIES_H

namespace mhn {

// Where pmhn and qmhn try the series before the quadrature: for a shape
// below 1 at every tilt, and otherwise for z = gamma / sqrt(beta) with
// |z| <= kCdfSeriesMaxAbsZ and alpha <= kCdfSeriesMaxAlpha, the bounds within
// which the normalising constant takes the Lemma 10 series for a positive
// tilt (mhn_psi_dispatch.cpp).  For z above 1e3 the shifted quadrature
// serves instead of either (see log_lower_general in mhn_pmhn_state.h).
//
// The choice is measured, below the peak of the density of log X, where the
// series is tried.  Over shapes from 0.05 to 1e4 and tilts up to 1e3 in size
// the quadrature (mhn_cdf_integrate.cpp) costs 0.02 to 13 us per
// evaluation.  The series needs about z^2 terms, and more for a large shape,
// each with an incomplete gamma function: 26 to 98 us at |z| = 20 with
// alpha <= 100, 0.8 to 10 ms at |z| = 100 and 7 to 18 ms at 300.  From
// alpha = 100 on it costs 8 to 250 us at a positive tilt up to 20, and 0.1 to
// 9 ms at a negative one, where its cancellation guard then refuses it
// almost always.  Against an independent reference (the integral of the law
// of log X measured from its peak, at 4e5 points with shapes from 0.01 to
// 1e4), outside |z| <= 20 with alpha <= 100 the series was off by up to
// 2.5e-9.  For a shape of 1 or more the quadrature was off by up to 1.1e-9
// for alpha from 1.01 to 1.3, under a tilt of either sign, up to 6e-9 at
// scattered points for alpha within 0.01 of 1, and by 3e-10 or less
// elsewhere, all within the working tolerance, and there it replaces the
// series.  For a shape below 1 the quadrature, tanh-sinh against the
// power t^(alpha - 1) near the origin, missed by more than the working
// tolerance at scattered points, about one in 4000 under a negative tilt (by
// up to 3.5e-7) and one in 35000 under a positive one, where the series was
// within 3.1e-10; so below 1 the series is tried at every tilt, at its cost.
// Inside the region the series is kept at every shape: it costs 1.2 to 500
// times as much as the quadrature, but for a small shape it is the more
// accurate, to about 1e-13 against up to 1.4e-9 and, at scattered points,
// more than the working tolerance, and it gives the values of the paper's
// formula.
constexpr double kCdfSeriesMaxAbsZ = 20.0;
constexpr double kCdfSeriesMaxAlpha = 100.0;

// Whether the truncation length of the series for this triple, which does
// not depend on x, stays within the ceiling of log_cdf_series.  Where it
// does not, from |z| of about 830 for a shape below 1, log_cdf_series
// returns NaN at every x, and callers need not try it.
bool cdf_series_fits(double alpha, double beta, double gamma, double tol);

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
