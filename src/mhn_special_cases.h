// Special-case detection predicates and dispatch helpers for the three
// reductions of MHN(alpha, beta, gamma) noted in Sun et al. (2023)
// Lemma 6: gamma = 0 (sqrt-Gamma), alpha = 1 (truncated normal), and
// alpha = 1 with gamma = 0 (half-normal).  Density, CDF, and quantile
// helpers all live here, and the dispatcher in dmhn / pmhn / qmhn routes a
// call to them with is_sqrt_gamma_exact() / is_truncated_normal_exact(),
// after giving a subnormal shape to mhn_vanishing.h (vanishing_candidate)
// and leaving gamma = 0 above alpha = 1e6 to the peak frame
// (sqrt_gamma_in_peak_frame).
//
// Each CDF and quantile is written once, as a function of one element and of
// constants that depend only on the parameter triple.  The all-scalar fast
// paths loop over those functions, and the recycling loops call them with the
// constants cached per triple, so the two call shapes return the same double.

#ifndef MHN_SPECIAL_CASES_H
#define MHN_SPECIAL_CASES_H

#include <Rcpp.h>
#include <Rmath.h>
#include <cfloat>
#include <cmath>
#include <limits>

#include "mhn_fma.h"

namespace mhn {

// dmhn, pmhn, qmhn and the normalising constant take the special cases at
// gamma = 0 and alpha = 1 exactly, and nowhere else.  With U = sqrt(beta) X,
// d log Psi[alpha/2, z] / d alpha = E log U and d log Psi / dz = E U
// (differentiate Psi = 2 int u^(alpha-1) exp(-u^2 + z u) du under the
// integral), so at alpha = 1 + delta the log density differs from the
// truncated normal's by about delta (log(sqrt(beta) x) - E log U), and at a
// tilt z = gamma / sqrt(beta) from the square root of a Gamma law's by about
// z (sqrt(beta) x - E U).  Taken as the special case, a delta of 1.5e-8
// would put the density 1e-5 off next to x = 1e-300, where
// log(sqrt(beta) x) = -691, and a z of 1.5e-8 the upper tail 8e-8 off at
// sqrt(beta) x = 6.  The general forms keep their accuracy up to both
// values.
inline bool is_sqrt_gamma_exact(double gamma) { return gamma == 0.0; }
inline bool is_truncated_normal_exact(double alpha) { return alpha == 1.0; }

// rmhn draws from the special-case law also for a shape within sqrt(eps) of
// 1 and a tilt within sqrt(eps) of 0 (see closed_form_triple in
// mhn_rmhn.cpp).  If X ~ MHN(alpha, beta, gamma) then sqrt(beta) X ~
// MHN(alpha, 1, gamma / sqrt(beta)) (change variables in the density of Sun
// et al. 2023, Lemma 1a), so the shape depends on gamma only through the
// scale-free tilt Delta = gamma / sqrt(beta), and the "is the tilt
// negligible?" test has to be made on Delta, not on gamma.
// beta defaults to 1, for which the two coincide.
bool is_sqrt_gamma(double gamma, double beta = 1.0);
bool is_truncated_normal(double alpha);

// ---------------------------------------------------------------------------
// gamma = 0: the square root of a gamma variable.
//
// By Sun et al. (2023), Lemma 6a, X^2 ~ Gamma(alpha/2, rate beta), so
// Y = beta X^2 ~ Gamma(alpha/2, 1).  Everything below is written in the
// scale-free variable t = sqrt(beta) x, with Y = t^2:
//
//   log f(x) = log 2 + (1/2) log beta - lgamma(alpha/2)
//              + (alpha - 1) log t - t^2,
//   F(x)     = P(alpha/2, t^2),
//
// with P the regularised lower incomplete gamma function.  The square is never
// formed as x * x: that underflows once x falls below 1.5e-154 and is exactly 0
// below 1e-162, where f and F are still ordinary numbers (f(1e-170) is
// 2/sqrt(pi) at the defaults).  t^2 leaves the normal range only where
// beta x^2 itself is below DBL_MIN, and there exp(-t^2) = 1 and
// P(a, t^2) = t^(2a) / Gamma(a + 1) hold exactly in double precision, so those
// elements are taken from log t directly.
//
// These forms serve up to alpha = 1e6, at every beta, since nothing in them
// depends on the size of beta.  Above it d, p and q take gamma = 0 in the
// peak frame of the general case (mhn_cdf_peak.h, sqrt_gamma_in_peak_frame),
// as for any other tilt: the Gamma law of t^2 has a relative width of
// sqrt(2 / alpha), and t^2, rounded to within about 3 eps relative, is off by
// up to 3 eps sqrt(alpha / 2) of its standard deviations, 1.5e-8 at
// alpha = 1e15 and 4.7e-6 at alpha = 1e20.  The law is then close to normal,
// so within three standard deviations of the peak that moves F by at most
// phi(0), 0.4, times the shift and log f by at most three times it: F by up
// to 6e-9 and log f by up to 4.5e-8 at alpha = 1e15, and by up to 1.9e-6 and
// 1.4e-5 at alpha = 1e20 (the errors measured reach about half of these at a
// beta that is not a power of 2, where beta x^2 rounds twice, and less at
// one that is).  It puts them anywhere at all once the law is narrower than
// the spacing of the doubles (alpha of about 1e31).
// ---------------------------------------------------------------------------

// Per-triple constants of the gamma = 0 log density.
struct SqrtGammaDensity {
  double alpha = 1.0;
  double shape = 0.5;           // alpha / 2
  double sqrt_beta = 1.0;
  double log_sqrt_beta = 0.0;
  double log_c_dgamma = 0.0;    // log 2 + (1/2) log beta
  double log_c = 0.0;           // log_c_dgamma - lgamma(alpha/2), alpha <= 20
  double log_zero = 0.0;        // log f(0)
  bool use_dgamma = false;

  void set(double a, double beta) {
    alpha = a;
    shape = 0.5 * a;
    // The constants of beta are kept while beta stays the same, as it does
    // when alpha alone is recycled.
    if (!(beta == beta_)) {
      beta_ = beta;
      sqrt_beta = std::sqrt(beta);
      log_sqrt_beta = 0.5 * std::log(beta);
      log_c_dgamma = M_LN2 + log_sqrt_beta;
    }
    // The closed form adds terms of size about (alpha/2) log(alpha/2) that
    // cancel to a log density of order one, so it loses about eps times that
    // much in log f: 4e-15 at alpha = 20, 4e-14 at 100, 6e-13 at 1e3.  R's
    // dgamma evaluates the same density through a saddle-point expansion
    // (its dpois_raw), which has no such cancellation, and is kept for the
    // larger shapes; up to alpha = 20 the two are equally accurate.  The
    // larger shapes need lgamma(alpha/2) only below beta x^2 = DBL_MIN, so it
    // is formed there, and a triple met at a single element does not pay for
    // it otherwise.
    use_dgamma = (a > 20.0);
    log_c = use_dgamma ? R_NaN : log_c_dgamma - R::lgammafn(shape);
    // f(0) is 0 for alpha > 1, +Inf for alpha < 1, and 2 sqrt(beta / pi) for
    // the half-normal.
    if (a > 1.0)      log_zero = R_NegInf;
    else if (a < 1.0) log_zero = R_PosInf;
    else              log_zero = log_c_dgamma - 0.5 * std::log(M_PI);
  }

  // log f(x) for x not NaN.
  double log_at(double x) const {
    if (x > 0.0) {
      const double t = sqrt_beta * x;
      if (t == R_PosInf) return R_NegInf;  // f vanishes at infinity
      const double y = t * t;
      if (use_dgamma && y >= DBL_MIN) {
        // f(x) = 2 beta x g(beta x^2) with g the Gamma(alpha/2, 1) density.
        return log_c_dgamma + std::log(t) + R::dgamma(y, shape, 1.0, 1);
      }
      const double log_t = (t >= DBL_MIN) ? std::log(t)
                                          : log_sqrt_beta + std::log(x);
      const double c = use_dgamma ? log_c_dgamma - R::lgammafn(shape) : log_c;
      return c + (alpha - 1.0) * log_t - y;
    }
    if (x == 0.0) return log_zero;
    return R_NegInf;
  }

 private:
  double beta_ = std::numeric_limits<double>::quiet_NaN();
};

// ---------------------------------------------------------------------------
// alpha = 1: the truncated normal.
//
// By Sun et al. (2023), Lemma 6b, X is N(mu, sigma^2) restricted to
// (0, Inf), with mu = gamma/(2 beta) and sigma = 1/sqrt(2 beta).  In the
// standardised variable, with z0 = -mu/sigma the truncation point and
// d = x/sigma the distance above it,
//
//   F(x) = (Phi(z0 + d) - Phi(z0)) / Phibar(z0) = lambda(z0) I(d),
//   I(d) = int_0^d exp(-z0 s - s^2/2) ds,
//
// where lambda(z) = phi(z)/Phibar(z) is the normal hazard, since
// phi(z0 + s) = phi(z0) exp(-z0 s - s^2/2).  The second form carries no
// subtraction, and is what makes a small lower tail accurate.
// ---------------------------------------------------------------------------

// Beyond this |z0| a truncated normal is measured from gamma and beta rather
// than from mu = gamma / (2 beta) (see TruncNormal::zq and the density).
const double kTnLargeZ0 = 1e3;

// Per-triple constants of the truncated normal.
struct TruncNormal {
  double mu = 0.0;
  double sigma = 1.0;
  double z0 = 0.0;
  // Below z0 = 3 the CDF works with log_mass = log Phibar(z0) =
  // log Phi(mu/sigma), the log mass kept by the truncation; from z0 = 3 on
  // it works with erfcx0 = erfcx(z0/sqrt2) instead.  log_mass is formed up to
  // z0 = kTnQnormMaxZ0, where the quantile still uses it, and erfcx0 from
  // z0 = 3 on; each is NaN elsewhere.
  double log_mass = 0.0;
  double erfcx0 = 0.0;
  // Phibar(z0) itself, beside log_mass, for the quantile's linear form.
  double mass = 0.0;
  // The parameters themselves, for the offsets below.
  double beta = 1.0;
  double gamma = 0.0;
  // Where z0 = -gamma / sqrt(2 beta) is +Inf the law is Exponential(rate
  // -gamma) to double precision (see set); the tails and quantiles are then
  // taken from it.
  bool exp_limit = false;
  // Where (q - mu) / sigma is formed from 2 beta q - gamma by fma (see zq),
  // and where it is formed as written.
  bool fma_offset = false;
  bool plain_offset = true;

  void set(double beta, double gamma);

  // The standardised position (q - mu) / sigma.  Under a strong positive tilt,
  // z0 < -kTnLargeZ0, q and mu are close and mu carries the rounding of
  // gamma / (2 beta), half a spacing of the doubles at mu, which is a
  // fraction eps |z0| / 2 of sigma: the difference is then taken as
  // (beta q - gamma / 2) (2 sigma), since 1 / (2 beta sigma) = sigma, with
  // beta q - gamma / 2 rounded once (fma).  That is half of 2 beta q - gamma
  // rounded once, so the two forms agree to the bit wherever 2 beta is a
  // double, and this one needs no 2 beta, which overflows for beta above
  // about 9e307.  The same form serves where mu is beyond the largest double
  // under a positive tilt, and z0 + q / sigma where it is under a negative
  // one.
  double zq(double q) const {
    if (plain_offset) return (q - mu) / sigma;
    if (fma_offset) return mhn::fma(beta, q, -0.5 * gamma) * (2.0 * sigma);
    return z0 + q / sigma;
  }

  // log erfcx(z0/sqrt2), for z0 >= 3; log lambda(z0), lambda(z0) itself (0
  // where it underflows) and log Phi(z0).  Each is formed the first time an
  // element needs it and then kept, so a triple shared by many elements pays
  // for it once and a triple met at a single element pays only for what that
  // element uses.
  double log_erfcx0() const {
    if (!have_log_erfcx0_) {
      log_erfcx0_ = std::log(erfcx0);
      have_log_erfcx0_ = true;
    }
    return log_erfcx0_;
  }
  double log_hazard() const {
    if (!have_hazard_) make_hazard();
    return log_hazard_;
  }
  double hazard() const {
    if (!have_hazard_) make_hazard();
    return hazard_;
  }
  double log_lower0() const {
    if (!have_lower0_) {
      log_lower0_ = R::pnorm(z0, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
      have_lower0_ = true;
    }
    return log_lower0_;
  }

 private:
  void make_hazard() const;
  mutable bool have_log_erfcx0_ = false;
  mutable double log_erfcx0_ = 0.0;
  mutable bool have_hazard_ = false;
  mutable double log_hazard_ = 0.0;
  mutable double hazard_ = 0.0;
  mutable bool have_lower0_ = false;
  mutable double log_lower0_ = 0.0;
};

// ---------------------------------------------------------------------------
// Up to this z0 the quantile of the truncated normal is taken through qnorm
// (see qmhn_truncated_normal_one).
const double kTnQnormMaxZ0 = 30.0;

// The log density of the truncated normal, shared by the all-scalar and the
// recycled paths of dmhn.  Up to |z0| = kTnLargeZ0 it is
// dnorm(x, mu, sigma, log) - log Phi(mu / sigma).  Beyond, each of those two
// logarithms is about -z0^2 / 2 under a negative tilt, and their difference
// carries eps z0^2 / 2 of rounding (5e-7 at z0 = 7e4); it is then taken as
// log lambda(z0) - log sigma - d (z0 + d / 2), d = x / sigma, the ratio
// phi(z0 + d) / (sigma Phibar(z0)) written with the normal hazard
// lambda(z0) = phi(z0) / Phibar(z0) (see TruncNormal).  Under a positive tilt
// the offset from mu is taken from zq instead, and where z0 = +Inf the density
// is the Exponential one, log(-gamma) + gamma x.
struct TruncNormalDensity {
  TruncNormal tn;
  // 0: dnorm - log Phi(mu / sigma); 1: through the hazard (z0 > kTnLargeZ0);
  // 2: through zq (z0 < -kTnLargeZ0, or mu beyond the doubles under a
  // positive tilt); 3: the Exponential limit (z0 = +Inf).
  int form = 0;
  double log_norm = 0.0;   // log Phi(mu / sigma), for form 0
  double log_sigma = 0.0;
  double log_hazard = 0.0;
  void set(double beta, double gamma);
  // log f(x) for x not NaN.
  double log_at(double x) const;
};

// Density, CDF and quantile helpers.  Each honours R's lower.tail / log.p
// conventions.  The vector forms serve the all-scalar-parameter fast paths;
// the *_one forms are the per-element kernels both paths share.
// ---------------------------------------------------------------------------

// Takes a log density to the density in place.  An NA or NaN is left as it
// is: the exponential of NA is NA or NaN as the platform has it (R's ?NaN),
// and NA and NaN are distinct in the d/p/q contract.
inline void exp_in_place(Rcpp::NumericVector& log_f) {
  for (double& v : log_f) {
    if (!ISNAN(v)) v = std::exp(v);
  }
}

// log = false: density; log = true: log-density.
Rcpp::NumericVector dmhn_sqrt_gamma(const Rcpp::NumericVector& x,
                                    double alpha, double beta, bool log_p);

Rcpp::NumericVector dmhn_truncated_normal(const Rcpp::NumericVector& x,
                                          double beta, double gamma,
                                          bool log_p);

// CDF for gamma == 0 (sqrt-Gamma): F(q) = P(alpha/2, beta q^2).
//   q <= 0  -> 0  (or -Inf if log.p)
//   q = Inf -> 1  (or 0 if log.p)
// shape = alpha/2 and sqrt_beta = sqrt(beta) are the per-triple constants.
double pmhn_sqrt_gamma_one(double q, double shape, double sqrt_beta,
                           bool lower_tail, bool log_p);
Rcpp::NumericVector pmhn_sqrt_gamma(const Rcpp::NumericVector& q,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p);

// CDF for alpha == 1 (truncated normal, Lemma 6b).
//   F(q) = (Phi((q - mu)/sigma) - Phi(-mu/sigma)) / Phi(mu/sigma).
double pmhn_truncated_normal_one(double q, const TruncNormal& tn,
                                 bool lower_tail, bool log_p);
Rcpp::NumericVector pmhn_truncated_normal(const Rcpp::NumericVector& q,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p);

// Quantile for gamma == 0: qmhn(p) = sqrt(qgamma(p, alpha/2, 1)) / sqrt(beta).
//   p = 0 -> 0, p = 1 -> +Inf.
double qmhn_sqrt_gamma_one(double p, double shape, double sqrt_beta,
                           bool lower_tail, bool log_p);
Rcpp::NumericVector qmhn_sqrt_gamma(const Rcpp::NumericVector& p,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p);

// Quantile for alpha == 1 (truncated normal inverse).
double qmhn_truncated_normal_one(double p, const TruncNormal& tn,
                                 bool lower_tail, bool log_p);
Rcpp::NumericVector qmhn_truncated_normal(const Rcpp::NumericVector& p,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p);

}  // namespace mhn

#endif  // MHN_SPECIAL_CASES_H
