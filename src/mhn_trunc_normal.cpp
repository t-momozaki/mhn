// Implementation of the one-sided truncated-normal sampler declared in
// mhn_trunc_normal.h, following Robert (1995).

#include "mhn_trunc_normal.h"
#include "mhn_fma.h"

#include <Rcpp.h>
#include <cmath>

namespace mhn {

double rtnorm_robert_zero(double mu, double sigma, double mu_lo) {
  // Standardized lower bound: a = (0 - mu) / sigma.
  const double a = -mu / sigma;

  // Both loops below accept with probability bounded well away from zero
  // (>= 0.5 in the first, >= ~0.76 in the second), so they are unbounded only
  // in principle.  The interrupt check costs nothing at this stride and keeps
  // the call escapable if an RNG ever misbehaves.
  long iter = 0;

  if (a <= 0.0) {
    // mu >= 0: rejection from the untruncated N(mu, sigma).
    // Acceptance probability is 1 - Phi(a) = Phi(mu/sigma) >= 0.5.
    //
    // R::rnorm(mu, sigma) returns mu + sigma Z rounded, Z = norm_rand().
    // With mu_lo the draw is mu + (mu_lo + sigma Z), the inner sum formed in
    // a single rounding by fma: the exact draw from the mean mu + mu_lo,
    // rounded once more to the grid of doubles at mu, whose spacing is far
    // above any error left in mu_lo + sigma Z.  Z is drawn the same way in
    // both branches.
    while (true) {
      const double x = (mu_lo == 0.0)
                         ? R::rnorm(mu, sigma)
                         : mu + mhn::fma(sigma, R::norm_rand(), mu_lo);
      if (x >= 0.0) return x;
      if ((++iter & 1023L) == 0L) Rcpp::checkUserInterrupt();
    }
  }

  // mu < 0: translated-exponential proposal of Robert (1995), Section 2.1
  // (Proposition 2.3).  Its rate alpha_star = (a + sqrt(a^2 + 4)) / 2
  // maximises the acceptance probability.  From a of about 2e8, a^2 + 4
  // rounds to a^2, whose square root rounds to a, so that the expression
  // gives a itself; a^2 overflows from 1.3e154, and from a = 1e150 on
  // alpha_star is taken as a directly.  The proposal is then accepted at
  // once, as y - alpha_star = e is about 1 / a.
  const double alpha_star =
    (a < 1e150) ? 0.5 * (a + std::sqrt(a * a + 4.0)) : a;
  // R::rexp takes the scale (= 1/rate), so pass 1/alpha_star.
  const double inv_rate = 1.0 / alpha_star;
  while (true) {
    const double excess = R::rexp(inv_rate);   // e >= 0, the amount past a
    const double y = a + excess;
    const double diff = y - alpha_star;
    const double log_u = std::log(R::runif(0.0, 1.0));
    if (log_u < -0.5 * diff * diff) {
      // The draw is mu + sigma * y with y = a + e and a = -mu / sigma, so the
      // mu cancels exactly: mu + sigma (a + e) = sigma e.  The literal
      // mu + sigma * y would subtract two nearly equal quantities once mu is
      // far below zero; for gamma near -1e8 that gives negative values --
      // outside the support -- and puts 2e5 draws on 13 distinct doubles.
      return sigma * excess;
    }
    if ((++iter & 1023L) == 0L) Rcpp::checkUserInterrupt();
  }
}

}  // namespace mhn

// [[Rcpp::export(.rtnorm_robert_cpp)]]
Rcpp::NumericVector rtnorm_robert_cpp(int n, double mu, double sigma) {
  if (n < 0) Rcpp::stop("n must be non-negative");
  if (n == 0) return Rcpp::NumericVector(0);
  if (!(sigma > 0.0)) Rcpp::stop("sigma must be positive");
  if (!std::isfinite(mu)) Rcpp::stop("mu must be finite");
  if (!std::isfinite(mu / sigma)) Rcpp::stop("mu / sigma must be finite");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = mhn::rtnorm_robert_zero(mu, sigma);
  }
  return out;
}
