// Fox-Wright Psi[alpha/2, gamma/sqrt(beta)] evaluation interfaces.
// Four numerical paths are exposed -- closed form for alpha = 1 and
// alpha = 2 (Sun et al. 2023 Supplementary, Lemma 9c), series
// expansion for gamma > 0 (Lemma 10), and numerical integration for
// gamma < 0 (Lemma 11) -- plus a dispatcher that picks the right one.

#ifndef MHN_PSI_H
#define MHN_PSI_H

namespace mhn {

// Closed-form Psi for alpha = 1 (Sun et al. 2023 Supplementary, Lemma 9(c)).
// Returns log Psi[1/2, gamma/sqrt(beta)], stable for all real gamma.
double psi_alpha1(double gamma, double beta);

// Closed-form Psi for alpha = 2 with gamma >= 0, derived from the ratio
// identity in Sun et al. (2023) Lemma 9(c).
double psi_alpha2(double gamma, double beta);

// Series expansion for gamma > 0
// (Sun et al. 2023 Supplementary, Lemma 10).
// Returns log Psi[alpha/2, gamma/sqrt(beta)].
double psi_series(double alpha, double beta, double gamma, double tol);

// Sun et al. (2023) Supplementary, Lemma 10(a, b) -- smallest integer
// k beyond which A(k) (when c1 = 6q, c2 = 4q) or B(k) (when c1 = 10q,
// c2 = 12q) is strictly decreasing in k.  Used both by psi_series and
// by the CDF series in mhn_cdf_series.cpp.
//
// ERRATA: the discriminant uses alpha_adj * z^2 (the form proved in
// the supplementary), not alpha_adj * z as in the main-text statement.
long lemma10_C(double alpha_adj, double z2, double q,
               double c1, double c2);

// Sun et al. (2023) Supplementary, Lemma 10(d) -- smallest integer
// k >= C for which A(k) (is_A = true) or |B(k)| (is_A = false) is
// guaranteed to be below the prescribed log_eps_quarter, derived from
// the geometric decay rate q established in Lemma 10(a, b, c).  The
// pair (K_1, K_2) returned for the A and B series bounds the
// truncation error of Psi[alpha/2, z] -- and, via the relation
// |T_{2k}| <= A(k), |T_{2k+1}| <= |B(k)|, also bounds the truncation
// error of the Sun et al. (2023) Lemma 1(b) CDF series.
// (Note: that per-term bound is a package-internal derivation; it is
//  not stated in Sun et al. (2023), but follows immediately from
//  comparing the Lemma 1(b) CDF term to A(k) and B(k) defined in
//  Lemma 10.)
long lemma10_K(double a, double z, double z2, long C, double q,
               double log_eps_quarter, bool is_A);

// Numerical integration for gamma < 0
// (Sun et al. 2023 Supplementary, Lemma 11).
// Returns log Psi[alpha/2, gamma/sqrt(beta)].
double psi_integrate(double alpha, double beta, double gamma, double tol);

// Dispatcher: log Psi[alpha/2, gamma/sqrt(beta)] for all (alpha, beta, gamma).
// Mirrors R-side `.mhn_log_normalizing_const`.
double mhn_log_normalizing_const(double alpha, double beta, double gamma,
                                 double tol);

// Resolves a sentinel tol (any non-positive value) to sqrt(eps).
double resolve_psi_tol(double tol);

}  // namespace mhn

#endif  // MHN_PSI_H
