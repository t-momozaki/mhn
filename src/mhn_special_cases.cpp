// Special-case detection and density dispatch.

#include "mhn_special_cases.h"
#include "mhn_stable.h"
#include "mhn_constants.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mhn {

bool is_sqrt_gamma(double gamma, double beta) {
  // Test the scale-free tilt Delta = gamma / sqrt(beta).  Testing |gamma|
  // alone silently discarded a fully significant tilt whenever beta was
  // small: at beta = 4e-18 a gamma of 1e-8 is Delta = 5, an ordinary value,
  // yet it was routed to the sqrt-Gamma closed form.
  return std::fabs(gamma) < mhn_eps() * std::sqrt(beta);
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
    // NA and NaN are distinct in the base R d/p/q contract, and Rcpp's
    // is_na does not separate them -- it is ISNAN -- so NaN was
    // silently promoted to NA on these paths.
    if (ISNAN(xi)) {
      log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
    } else if (xi == R_PosInf) {
      // Zero at infinity, as for every base R density; the expression below
      // gives log(Inf) + (-Inf) = NaN.
      log_f[i] = R_NegInf;
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
    if (ISNAN(xi)) {
      log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
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
    if (ISNAN(qi)) {
      out[i] = R_IsNA(qi) ? NA_REAL : R_NaN;
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
    if (ISNAN(qi)) {
      out[i] = R_IsNA(qi) ? NA_REAL : R_NaN;
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
    const double z_0 = -mu / sigma;

    // Compute whichever tail is the small one directly, and derive the other
    // from it.  Always working from the upper tail loses the lower tail
    // entirely once it falls below the rounding of 1: P(X > q) is then 1 to
    // machine precision, and 1 - that is exactly 0 rather than the true value
    // (1.3e-17 at alpha = 1, gamma = 18, q = 3, for instance).
    double log_lower, log_upper;
    if (z_q <= 0.0) {
      // qi below the untruncated mean: the lower tail is the small one, and it
      // is a difference of two lower-tail normals, both small and ordered.
      const double lo_q = R::pnorm(z_q, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
      const double lo_0 = R::pnorm(z_0, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
      const double d = lo_0 - lo_q;                       // <= 0
      log_lower = (d < 0.0) ? lo_q + std::log1p(-std::exp(d)) : R_NegInf;
      log_lower -= log_denom;
      log_upper = (log_lower < 0.0) ? std::log1p(-std::exp(log_lower)) : R_NegInf;
    } else if (z_0 >= 3.0) {
      // Both ordinates are far into the right tail, where log Phi_upper(z) is
      // about -z^2/2.  Subtracting the two logs then cancels: at
      // (1, 1e-4, -1000) each is -2.5e9 while the answer is -0.01, so the
      // result carried an absolute error of 2.5e9 * eps = 5.5e-7 -- thirty
      // times the working tolerance, and growing as z^2.
      //
      // Writing Phi_upper(z) = exp(-z^2/2) erfcx(z/sqrt2) / 2 removes it:
      //
      //   log Phi_up(z_q) - log Phi_up(z_0)
      //       = -(z_q^2 - z_0^2)/2 + log erfcx(z_q/sqrt2) - log erfcx(z_0/sqrt2)
      //
      // and z_q^2 - z_0^2 factors as d (2 z_0 + d) with d = z_q - z_0 = q/sigma,
      // which only ever adds.  erfcx decays like 1/(z sqrt(pi)), so the two
      // remaining logarithms are of order one and their difference is exact.
      const double d = qi / sigma;                        // = z_q - z_0
      const double quad = 0.5 * d * (2.0 * z_0 + d);      // (z_q^2 - z_0^2)/2
      log_upper = -quad
                  + std::log(mhn::erfcx_cf(z_q * M_SQRT1_2))
                  - std::log(mhn::erfcx_cf(z_0 * M_SQRT1_2));
      log_lower = (log_upper < 0.0) ? std::log1p(-std::exp(log_upper)) : R_NegInf;
    } else {
      const double log_upper_num = R::pnorm(z_q, 0.0, 1.0,
                                            /*lower_tail=*/0, /*log_p=*/1);
      log_upper = log_upper_num - log_denom;              // log P(X > qi)
      log_lower = (log_upper < 0.0) ? std::log1p(-std::exp(log_upper)) : R_NegInf;
    }

    // Translate to the requested (lower.tail, log.p) combination.
    if (log_p) {
      out[i] = lower_tail ? log_lower : log_upper;
    } else {
      out[i] = std::exp(lower_tail ? log_lower : log_upper);
    }
  }
  return out;
}

Rcpp::NumericVector qmhn_sqrt_gamma(const Rcpp::NumericVector& p,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p) {
  const double shape = alpha / 2.0;
  const double scale = 1.0 / beta;
  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double pi = p[i];
    if (R_IsNA(pi))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(pi)) { out[i] = R_NaN;   continue; }
    if (!log_p && (pi < 0.0 || pi > 1.0)) { out[i] = R_NaN; continue; }
    if (log_p && pi > 0.0)                { out[i] = R_NaN; continue; }
    // Hand the tail flags to qgamma rather than converting first.  Turning an
    // upper-tail p into 1 - p loses it below about 1e-16 -- the quantile came
    // back as Inf at 1e-20, where qgamma itself still resolves 1e-30 -- and
    // exponentiating a log-scale p throws away its whole point.
    const double q_g = R::qgamma(pi, shape, scale,
                                 lower_tail ? 1 : 0, log_p ? 1 : 0);
    out[i] = std::sqrt(std::max(q_g, 0.0));
  }
  return out;
}

Rcpp::NumericVector qmhn_truncated_normal(const Rcpp::NumericVector& p,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p) {
  const double mu = gamma / (2.0 * beta);
  const double sigma = 1.0 / std::sqrt(2.0 * beta);

  // Work with the standardised lower bound z0 = -mu/sigma and with
  // t = x/sigma >= 0, so that the truncation identity reads
  //
  //     Phibar(t + z0) = Phibar(z0) * (1 - p),
  //
  // and hence  t + z0 = qnorm(log Phibar(z0) + log1p(-p), upper, log).
  //
  // Two things are avoided by writing it this way.  The truncation mass was
  // previously formed as 1 - Phi(z0), which underflows to exactly 0 once z0
  // passes about 8.3 -- reached at gamma = -12 with beta = 1 -- after which
  // every quantile came back as Inf.  Taking it from pnorm's own upper tail on
  // the log scale keeps it representable however far out z0 goes.  And the
  // answer is assembled as sigma * t rather than mu + sigma * (t + z0): those
  // are equal, since sigma * z0 = -mu, but the second subtracts two nearly
  // equal quantities once mu is far below zero.
  const double z0 = -mu / sigma;
  const double log_mass = R::pnorm(z0, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/1);

  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double pi = p[i];
    if (R_IsNA(pi))  { out[i] = NA_REAL; continue; }
    if (R_IsNaN(pi)) { out[i] = R_NaN;   continue; }
    if (!log_p && (pi < 0.0 || pi > 1.0)) { out[i] = R_NaN; continue; }
    if (log_p && pi > 0.0)                { out[i] = R_NaN; continue; }

    // What the identity needs is log(1 - p_lower), which is the log of the
    // upper-tail probability.  When the caller supplied that directly -- as an
    // upper tail, on either scale -- use it as given rather than round-tripping
    // through 1 - p, which would destroy it below about 1e-16.
    double log_upper;
    if (!lower_tail) {
      log_upper = log_p ? pi : std::log(pi);
    } else {
      const double p_low = log_p ? std::exp(pi) : pi;
      if (p_low <= 0.0) { out[i] = 0.0; continue; }
      if (p_low >= 1.0) { out[i] = R_PosInf; continue; }
      log_upper = std::log1p(-p_low);
    }
    if (log_upper >= 0.0)      { out[i] = 0.0; continue; }
    if (log_upper == R_NegInf) { out[i] = R_PosInf; continue; }

    // Solve for the offset t = z - z0 directly.  Assembling it as
    // sigma * (z - z0) avoids the mu + sigma z cancellation, but only moves it:
    // at gamma = -1000, beta = 1e-4 both z and z0 are about 70710 while their
    // difference is 1.4e-7, so the subtraction keeps four digits of the answer.
    //
    // The defining equation, in the erfcx form that pmhn_truncated_normal uses
    // for the same reason,
    //
    //     -t (2 z0 + t)/2 + log erfcx((z0+t)/sqrt2) - log erfcx(z0/sqrt2)
    //         = log(1 - p),
    //
    // contains no such subtraction: to leading order t = -log(1-p)/z0, and
    // Newton on it converges in a few steps.  d/dt of the left side is
    // -2(z0+t) + sqrt(2/pi)/erfcx((z0+t)/sqrt2).
    double t;
    if (z0 >= 3.0) {
      const double e0 = mhn::erfcx_cf(z0 * M_SQRT1_2);
      const double log_e0 = std::log(e0);
      double d = -log_upper / z0;                 // leading order, > 0
      if (!(d > 0.0)) d = 0.0;
      for (int it = 0; it < 60; ++it) {
        const double ea = mhn::erfcx_cf((z0 + d) * M_SQRT1_2);
        const double G = -d * (2.0 * z0 + d) * 0.5 + std::log(ea) - log_e0
                         - log_upper;
        const double dG = -2.0 * (z0 + d)
                          + M_SQRT2 / (std::sqrt(M_PI) * ea);
        if (!std::isfinite(G) || !std::isfinite(dG) || dG == 0.0) break;
        const double step = G / dG;
        if (!std::isfinite(step)) break;
        d -= step;
        if (!(d > 0.0)) d = 0.0;
        if (std::fabs(step) <= 1e-15 * (1.0 + d)) break;
      }
      t = d;
    } else {
      const double z = R::qnorm(log_mass + log_upper, 0.0, 1.0,
                                /*lower_tail=*/0, /*log_p=*/1);
      t = z - z0;
    }
    out[i] = std::max(sigma * t, 0.0);
  }
  return out;
}

}  // namespace mhn

// ---------------------------------------------------------------------------
// R-visible exports
// ---------------------------------------------------------------------------

// [[Rcpp::export(.is_sqrt_gamma, rng = false)]]
bool is_sqrt_gamma_R(double gamma, double beta = 1.0) {
  return mhn::is_sqrt_gamma(gamma, beta);
}

// [[Rcpp::export(.is_truncated_normal, rng = false)]]
bool is_truncated_normal_R(double alpha) {
  return mhn::is_truncated_normal(alpha);
}
