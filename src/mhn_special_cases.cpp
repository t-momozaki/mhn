// Special-case detection and density dispatch.
// Mirrors R-side mhn/R/mhn_special_cases.R.

#include "mhn_special_cases.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

bool is_sqrt_gamma(double gamma) {
  return std::fabs(gamma) < mhn_eps();
}

bool is_truncated_normal(double alpha) {
  return std::fabs(alpha - 1.0) < mhn_eps();
}

Rcpp::NumericVector dmhn_sqrt_gamma(const Rcpp::NumericVector& x,
                                    double alpha, double beta, bool log_p) {
  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  // x = 0 boundary value (depends on alpha) -- precomputed.
  double log_zero;
  if (alpha > 1.0) {
    log_zero = R_NegInf;            // density 0
  } else if (alpha < 1.0) {
    log_zero = R_PosInf;            // density +Inf
  } else {
    // alpha == 1 (half-normal): density at 0 is 2 * sqrt(beta / pi)
    log_zero = std::log(2.0) + 0.5 * std::log(beta) - 0.5 * std::log(M_PI);
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      log_f[i] = NA_REAL;
    } else if (xi > 0.0) {
      // log_f = log(2) + log(x) + dgamma(x^2, alpha/2, rate=beta, log=TRUE)
      log_f[i] = std::log(2.0) + std::log(xi)
                 + R::dgamma(xi * xi, alpha / 2.0,
                             /*scale=*/1.0 / beta, /*log=*/1);
    } else if (xi == 0.0) {
      log_f[i] = log_zero;
    }
    // xi < 0 falls through to R_NegInf
  }

  if (!log_p) {
    // exp(NA) = NA, exp(-Inf) = 0
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}

Rcpp::NumericVector dmhn_truncated_normal(const Rcpp::NumericVector& x,
                                          double beta, double gamma,
                                          bool log_p) {
  const double mu = gamma / (2.0 * beta);
  const double sigma = 1.0 / std::sqrt(2.0 * beta);
  const double log_norm_const = R::pnorm(mu / sigma, 0.0, 1.0,
                                         /*lower_tail=*/1, /*log_p=*/1);

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      log_f[i] = NA_REAL;
    } else if (xi >= 0.0) {
      // Evaluate at xi = 0 as well; the truncated-normal density is
      // finite and positive at the lower boundary.
      log_f[i] = R::dnorm(xi, mu, sigma, /*log=*/1) - log_norm_const;
    }
    // xi < 0 falls through to R_NegInf
  }

  if (!log_p) {
    for (R_xlen_t i = 0; i < n; ++i) log_f[i] = std::exp(log_f[i]);
  }
  return log_f;
}

// ---------------------------------------------------------------------------
// CDF and quantile helpers for the special-case dispatch.
//
// Each accepts R's lower.tail / log.p flags directly and forwards them to
// R::pgamma / R::pnorm / R::qgamma / R::qnorm where possible.  The dispatcher
// in mhn_pmhn.cpp / mhn_qmhn.cpp is therefore free to pass user flags
// through without intermediate clipping, preserving precision in the tails.
// ---------------------------------------------------------------------------

Rcpp::NumericVector pmhn_sqrt_gamma(const Rcpp::NumericVector& q,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p) {
  const double shape = alpha / 2.0;
  const double scale = 1.0 / beta;
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double qi = q[i];
    if (Rcpp::NumericVector::is_na(qi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (qi <= 0.0) {
      // F = 0
      if (log_p) out[i] = lower_tail ? R_NegInf : 0.0;
      else       out[i] = lower_tail ? 0.0 : 1.0;
      continue;
    }
    if (qi == R_PosInf) {
      // F = 1
      if (log_p) out[i] = lower_tail ? 0.0 : R_NegInf;
      else       out[i] = lower_tail ? 1.0 : 0.0;
      continue;
    }
    out[i] = R::pgamma(qi * qi, shape, scale,
                       lower_tail ? 1 : 0, log_p ? 1 : 0);
  }
  return out;
}

Rcpp::NumericVector pmhn_truncated_normal(const Rcpp::NumericVector& q,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p) {
  const double mu = gamma / (2.0 * beta);
  const double sigma = 1.0 / std::sqrt(2.0 * beta);
  // P(Z > -mu/sigma) on log scale: denominator of the truncated-normal CDF.
  const double log_denom = R::pnorm(mu / sigma, 0.0, 1.0,
                                    /*lower_tail=*/1, /*log_p=*/1);

  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double qi = q[i];
    if (Rcpp::NumericVector::is_na(qi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (qi <= 0.0) {
      if (log_p) out[i] = lower_tail ? R_NegInf : 0.0;
      else       out[i] = lower_tail ? 0.0 : 1.0;
      continue;
    }
    if (qi == R_PosInf) {
      if (log_p) out[i] = lower_tail ? 0.0 : R_NegInf;
      else       out[i] = lower_tail ? 1.0 : 0.0;
      continue;
    }
    const double z_q = (qi - mu) / sigma;
    // Upper-tail of the underlying normal, evaluated in log space:
    //   log P(Z > z_q) - log P(Z > -mu/sigma)
    // = log Q where Q = P(X > qi).
    const double log_upper_num = R::pnorm(z_q, 0.0, 1.0,
                                          /*lower_tail=*/0, /*log_p=*/1);
    const double log_upper = log_upper_num - log_denom;  // log P(X > qi)

    // Translate to the requested (lower.tail, log.p) combination.
    if (log_p) {
      if (lower_tail) {
        // log(1 - P(X > qi)) = log1p(-exp(log_upper))
        if (log_upper >= 0.0) out[i] = R_NegInf;
        else                  out[i] = std::log1p(-std::exp(log_upper));
      } else {
        out[i] = log_upper;
      }
    } else {
      const double upper = std::exp(log_upper);
      if (lower_tail) out[i] = 1.0 - upper;
      else            out[i] = upper;
    }
  }
  return out;
}

namespace {

// Convert a user-provided (p, lower_tail, log_p) to a lower-tail probability
// p_low in [0, 1].  Returns NA_REAL for true NA, R_NaN for NaN or
// out-of-range inputs (matches R's qnorm / qgamma conventions).
double to_lower_tail_prob(double p, bool lower_tail, bool log_p) {
  if (R_IsNA(p)) return NA_REAL;
  if (R_IsNaN(p)) return R_NaN;
  double p_low;
  if (log_p) {
    if (p > 0.0) return R_NaN;            // log p > 0 -> p > 1
    if (p == R_NegInf) p_low = 0.0;
    else               p_low = std::exp(p);
  } else {
    if (p < 0.0 || p > 1.0) return R_NaN;
    p_low = p;
  }
  if (!lower_tail) p_low = 1.0 - p_low;
  // Numerical clamp.
  if (p_low < 0.0) p_low = 0.0;
  if (p_low > 1.0) p_low = 1.0;
  return p_low;
}

}  // namespace

Rcpp::NumericVector qmhn_sqrt_gamma(const Rcpp::NumericVector& p,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p) {
  const double shape = alpha / 2.0;
  const double scale = 1.0 / beta;
  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double p_low = to_lower_tail_prob(p[i], lower_tail, log_p);
    if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }
    if (p_low <= 0.0) { out[i] = 0.0; continue; }
    if (p_low >= 1.0) { out[i] = R_PosInf; continue; }
    const double q_g = R::qgamma(p_low, shape, scale,
                                 /*lower_tail=*/1, /*log_p=*/0);
    out[i] = std::sqrt(std::max(q_g, 0.0));
  }
  return out;
}

Rcpp::NumericVector qmhn_truncated_normal(const Rcpp::NumericVector& p,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p) {
  const double mu = gamma / (2.0 * beta);
  const double sigma = 1.0 / std::sqrt(2.0 * beta);
  const double Phi_low = R::pnorm(-mu / sigma, 0.0, 1.0, 1, 0);
  const double Phi_up  = 1.0 - Phi_low;                  // = Phi(mu/sigma)

  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double p_low = to_lower_tail_prob(p[i], lower_tail, log_p);
    if (R_IsNA(p_low))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(p_low)) { out[i] = R_NaN;   continue; }
    if (p_low <= 0.0) { out[i] = 0.0; continue; }
    if (p_low >= 1.0) { out[i] = R_PosInf; continue; }
    // F(x) = (Phi((x - mu)/sigma) - Phi_low) / Phi_up
    // Solve for x:
    //   (x - mu)/sigma = qnorm( Phi_low + p_low * Phi_up )
    // For p_low close to 1, evaluate upper tail to preserve precision:
    //   1 - (Phi_low + p_low * Phi_up) = Phi_up * (1 - p_low)
    double x;
    if (p_low < 0.5) {
      const double q_arg = Phi_low + p_low * Phi_up;
      x = mu + sigma * R::qnorm(q_arg, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
    } else {
      const double q_arg_up = Phi_up * (1.0 - p_low);
      x = mu + sigma * R::qnorm(q_arg_up, 0.0, 1.0,
                                /*lower_tail=*/0, /*log_p=*/0);
    }
    out[i] = std::max(x, 0.0);
  }
  return out;
}

}  // namespace mhn

// ---------------------------------------------------------------------------
// R-visible exports
// ---------------------------------------------------------------------------

// [[Rcpp::export(.is_sqrt_gamma)]]
bool is_sqrt_gamma_R(double gamma) {
  return mhn::is_sqrt_gamma(gamma);
}

// [[Rcpp::export(.is_truncated_normal)]]
bool is_truncated_normal_R(double alpha) {
  return mhn::is_truncated_normal(alpha);
}
