// Special-case detection predicates and dispatch helpers for the three
// reductions of MHN(alpha, beta, gamma) noted in Sun et al. (2023)
// Lemma 6: gamma = 0 (sqrt-Gamma), alpha = 1 (truncated normal), and
// alpha = 1 with gamma = 0 (half-normal).  Density, CDF, and quantile
// helpers all live here so the dispatcher in dmhn / pmhn / qmhn only
// needs is_sqrt_gamma() / is_truncated_normal() to route the call.

#ifndef MHN_SPECIAL_CASES_H
#define MHN_SPECIAL_CASES_H

#include <Rcpp.h>

namespace mhn {

bool is_sqrt_gamma(double gamma);
bool is_truncated_normal(double alpha);

// log = false: density; log = true: log-density.
Rcpp::NumericVector dmhn_sqrt_gamma(const Rcpp::NumericVector& x,
                                    double alpha, double beta, bool log_p);

Rcpp::NumericVector dmhn_truncated_normal(const Rcpp::NumericVector& x,
                                          double beta, double gamma,
                                          bool log_p);

// ---------------------------------------------------------------------------
// CDF and quantile helpers for the special-case dispatch.  Each returns the
// value(s) honouring R's lower.tail / log.p conventions, so the dispatcher
// can pass the user's flags straight through to R::pgamma / R::pnorm.
// ---------------------------------------------------------------------------

// CDF for gamma == 0 (sqrt-Gamma): F(q) = pgamma(q^2, alpha/2, scale = 1/beta).
//   q <= 0  -> 0  (or -Inf if log.p)
//   q = Inf -> 1  (or 0 if log.p)
Rcpp::NumericVector pmhn_sqrt_gamma(const Rcpp::NumericVector& q,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p);

// CDF for alpha == 1 (truncated normal, Lemma 6b).
//   mu = gamma/(2 beta), sigma = 1/sqrt(2 beta);
//   F(q) = (Phi((q - mu)/sigma) - Phi(-mu/sigma)) / Phi(mu/sigma).
Rcpp::NumericVector pmhn_truncated_normal(const Rcpp::NumericVector& q,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p);

// Quantile for gamma == 0: qmhn(p) = sqrt(qgamma(p, alpha/2, scale=1/beta)).
//   p = 0 -> 0, p = 1 -> +Inf.
Rcpp::NumericVector qmhn_sqrt_gamma(const Rcpp::NumericVector& p,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p);

// Quantile for alpha == 1 (truncated normal inverse).
Rcpp::NumericVector qmhn_truncated_normal(const Rcpp::NumericVector& p,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p);

}  // namespace mhn

#endif  // MHN_SPECIAL_CASES_H
