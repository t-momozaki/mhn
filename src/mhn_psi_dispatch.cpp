// Dispatcher for log Psi[alpha/2, gamma/sqrt(beta)].  Routes each call
// to the cheapest accurate evaluation method for the parameter region:
// the closed forms of Sun et al. (2023) Supplementary Lemma 9(c), the
// series of Supplementary Lemma 10, or numerical integration.

#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>

namespace mhn {

double resolve_psi_tol(double tol) {
  return (tol > 0.0) ? tol : mhn_eps();
}

double mhn_log_normalizing_const(double alpha, double beta, double gamma,
                                 double tol) {
  if (is_sqrt_gamma_exact(gamma)) {
    // Psi[alpha/2, 0] = Gamma(alpha/2).
    return R::lgammafn(alpha / 2.0);
  }
  // Each closed form holds at its own value exactly, and is taken only there.
  // Psi is a smooth function of both arguments: d log Psi[alpha/2, z] /
  // d alpha = E[log U] and d log Psi / dz = E[U], with U = sqrt(beta) X, so
  // returning the alpha = 1 value for a shape delta away would be off by
  // about delta E[log U], -delta (log|z| + 0.577) under a strong negative
  // tilt, 2e-7 at z = -1e6 for delta = 1.4e-8; and returning Gamma(alpha/2)
  // for a tilt z would be off by about z E[U].  The moment formulas divide
  // Psi at alpha + 1 by Psi at alpha, so such an offset in one of the two
  // would pass straight into the mean, and dmhn would pair it with the
  // kernel at the true parameters.  The general paths below hold their
  // accuracy right up to these points: within a few units in the last place
  // of log Psi (or of 1, where |log Psi| < 1) for shapes within sqrt(eps) of
  // 1 or 2 and for tilts within sqrt(eps) of 0.  (rmhn does treat a shape
  // within sqrt(eps) of 1 as 1, and such a tilt as 0, but it replaces the
  // whole law by the special case when it does, and forms no Psi.)
  if (alpha == 1.0) {
    // alpha = 1 (a = 1/2): closed form for all gamma
    // (Sun et al. 2023, Supplementary Lemma 9(c)).
    return psi_alpha1(gamma, beta);
  }
  if (alpha == 2.0 && gamma >= 0.0) {
    // alpha = 2 (a = 1), gamma >= 0: closed form via the ratio identity of
    // Sun et al. (2023) Supplementary Lemma 9(c).
    return psi_alpha2(gamma, beta);
  }
  // Psi[alpha/2, z] = 2 int_0^Inf u^(alpha-1) exp(-u^2 + z u) du depends on the
  // tilt only through z = gamma / sqrt(beta), and the quadrature serves either
  // sign of z.  The series of Sun et al. (2023) Supplementary Lemma 10 is used
  // for 0 < z <= 20 with alpha <= 100, and the quadrature everywhere else.
  //
  // The series has only positive terms, and below z = 20 it is as accurate as
  // the quadrature: both are within a few units in the last place of log Psi
  // on the shapes and tilts checked (the series within 4 eps (1 + |log Psi|)
  // at 23000 points against its sum in arbitrary precision).  But the length
  // of the series grows in two ways.  Its terms peak near index z^2/4, so
  // their number grows like z^2 and the rounding of their running sum with
  // it: past z of about 15 the quadrature is faster whatever the shape, and
  // above 20 it is also the more accurate, within about one unit in the last
  // place of log Psi on scans of several million points up to z = 1e4, where
  // the series is off by tens to hundreds of them near z = 1e3.  And the
  // Lemma 10(d) truncation length bounds an absolute error, so it also grows
  // with the size of the terms, by about log Gamma(alpha/2) / log 2: from
  // alpha = 100 on the series is the slower of the two at every z, and it is
  // no more accurate (16 units in the last place at alpha = 1e5, z = 10,
  // against one).  Under a negative tilt, where only the quadrature serves,
  // it is within 4.4 eps (1 + |log Psi|) at 1.27e5 points, and within 3.2 at
  // all but one of them (see psi_log_axis in mhn_psi_integrate.cpp).  The
  // series also starts from Gamma(alpha / 2), which is infinite where
  // alpha / 2 underflows to 0, at the smallest positive double; that shape
  // goes to the quadrature as well.
  if (gamma < 0.0 || alpha > 100.0 || gamma / std::sqrt(beta) > 20.0 ||
      !(0.5 * alpha > 0.0)) {
    return psi_integrate(alpha, beta, gamma, resolve_psi_tol(tol));
  }
  // 0 < z <= 20 and alpha <= 100: the Lemma 10 series, whose truncation length
  // there is at most about 700.
  const double series = psi_series(alpha, beta, gamma, resolve_psi_tol(tol));
  if (std::isnan(series)) {
    return psi_integrate(alpha, beta, gamma, resolve_psi_tol(tol));
  }
  return series;
}

}  // namespace mhn

// ---------------------------------------------------------------------------
// R-visible exports.  Each accepts a sentinel tol = -1.0 meaning "use the
// default, sqrt(.Machine$double.eps)".
// ---------------------------------------------------------------------------

// [[Rcpp::export(.mhn_log_normalizing_const, rng = false)]]
double mhn_log_normalizing_const_R(double alpha, double beta, double gamma,
                                   double tol = -1.0) {
  return mhn::mhn_log_normalizing_const(alpha, beta, gamma, tol);
}

// [[Rcpp::export(.psi_series, rng = false)]]
double psi_series_R(double alpha, double beta, double gamma,
                    double tol = -1.0) {
  return mhn::psi_series(alpha, beta, gamma, mhn::resolve_psi_tol(tol));
}

// [[Rcpp::export(.psi_integrate, rng = false)]]
double psi_integrate_R(double alpha, double beta, double gamma,
                       double tol = -1.0) {
  return mhn::psi_integrate(alpha, beta, gamma, mhn::resolve_psi_tol(tol));
}
