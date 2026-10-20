// Fox-Wright Psi[alpha/2, gamma/sqrt(beta)] evaluation interfaces.
// Four numerical paths are exposed -- closed form for alpha = 1 and
// alpha = 2 (Sun et al. 2023, Supplementary Lemma 9(c)), series
// expansion for gamma > 0 (Supplementary Lemma 10), and numerical
// integration of the integral representation for either sign of gamma --
// plus a dispatcher that picks the right one: the series for
// 0 < gamma/sqrt(beta) <= 20 with alpha <= 100, the integration elsewhere.

#ifndef MHN_PSI_H
#define MHN_PSI_H

namespace mhn {

// Closed-form Psi for alpha = 1 (Sun et al. 2023, Supplementary Lemma 9(c)).
// Returns log Psi[1/2, gamma/sqrt(beta)], stable for all real gamma.
double psi_alpha1(double gamma, double beta);

// Closed-form Psi for alpha = 2 with gamma >= 0, derived from the ratio
// identity in Sun et al. (2023) Supplementary Lemma 9(c).
double psi_alpha2(double gamma, double beta);

// Series expansion for gamma > 0
// (Sun et al. 2023, Supplementary Lemma 10).
// Returns log Psi[alpha/2, gamma/sqrt(beta)].
double psi_series(double alpha, double beta, double gamma, double tol);

// Sun et al. (2023) Supplementary Lemma 10(a, b) -- the index C (C_1 for
// A(k), when c1 = 6q, c2 = 4q; C_2 for B(k), when c1 = 10q, c2 = 12q) from
// which on the term ratio A(k+1)/A(k), resp. B(k+1)/B(k), is at most q, so
// the terms decay at least geometrically from there (Lemma 10(c) states that
// they are strictly decreasing).  Used both by psi_series and by the CDF
// series in mhn_cdf_series.cpp.
//
// ERRATA: the discriminant uses alpha_adj * z^2 (the form in the proof,
// Supplementary Section 2.2), not alpha_adj * z as printed in the statement
// of the lemma (Supplementary Section 1).
//
// The index, like the length below, is an integer held in a double: it is
// about z^2 / 2, which no integer type need hold (a 32-bit long overflows at
// |z| of about 6.6e4), and the callers test it against their length ceilings
// before converting it.
double lemma10_C(double alpha_adj, double z2, double q,
                 double c1, double c2);

// Sun et al. (2023) Supplementary Lemma 10(d) -- smallest integer
// k >= C for which A(k) (is_A = true) or |B(k)| (is_A = false) is
// guaranteed to be below the prescribed log_eps_quarter, derived from
// the geometric decay rate q established in Lemma 10(a, b).  The
// pair (K_1, K_2) returned for the A and B series bounds the
// truncation error of Psi[alpha/2, z] -- and, via the relation
// |T_{2k}| <= A(k), |T_{2k+1}| <= |B(k)|, also bounds the truncation
// error of the Sun et al. (2023) Lemma 1b CDF series.
// (Note: that per-term bound is a package-internal derivation; it is
//  not stated in Sun et al. (2023), but follows immediately from
//  comparing the Lemma 1b CDF term to A(k) and B(k) defined in
//  Lemma 10.)
// Truncation length for the Lemma 10(d) bound.  `a` is the Fox-Wright first
// index alpha/2 for both branches; the B branch adds its own +1/2.  Returned
// as an integer-valued double, for the reason given at lemma10_C.
double lemma10_K(double a, double z, double z2, double C, double q,
                 double log_eps_quarter, bool is_A);

// Numerical integration of
//   Psi[alpha/2, z] = 2 int_0^Inf u^(alpha-1) exp(-u^2 + z u) du,
// z = gamma / sqrt(beta), for either sign of gamma; the dispatcher uses it
// for gamma < 0, for z > 20 and for alpha > 100.  See mhn_psi_integrate.cpp.
// Returns log Psi[alpha/2, gamma/sqrt(beta)].
double psi_integrate(double alpha, double beta, double gamma, double tol);

// log Psi[alpha/2, z] - z^2/4 with z = gamma / sqrt(beta), for z > 0.  See the
// definition in mhn_psi_integrate.cpp: the subtracted term is the one that
// cancels against the kernel exponent, and removing it from both analytically
// is what keeps the log density accurate at a large tilt.
double psi_integrate_shifted(double alpha, double beta, double gamma,
                             double tol);

// The integral over the whole line of
//   exp(-alpha (e^w - 1 - w) - u*^2 (e^w - 1)^2),
// u* the positive root of 2u^2 - z u - alpha = 0, the peak of the integrand
// u^alpha exp(-u^2 + z u) of Psi in log u, with u = u* e^w; so that
// Psi[alpha/2, z] = 2 u*^alpha exp(-u*^2 + z u*) times this integral.  Returns
// its logarithm, formed without the peak exponent, which for a large shape or
// tilt is far larger than the result.  For z = -Inf (u* = 0) it is
// log(Gamma(alpha) e^alpha / alpha^alpha); z = +Inf returns NaN.
double psi_log_axis_relative(double alpha, double z, double tol);

// Dispatcher: log Psi[alpha/2, gamma/sqrt(beta)] for all (alpha, beta, gamma).
double mhn_log_normalizing_const(double alpha, double beta, double gamma,
                                 double tol);

// Resolves a sentinel tol (any non-positive value) to sqrt(eps).
double resolve_psi_tol(double tol);

}  // namespace mhn

#endif  // MHN_PSI_H
