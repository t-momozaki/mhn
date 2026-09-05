// Implementation of the one-sided truncated-normal sampler declared in
// mhn_trunc_normal.h, following Robert (1995).

#include "mhn_trunc_normal.h"

#include <Rcpp.h>
#include <cmath>

namespace mhn {

double rtnorm_robert_zero(double mu, double sigma) {
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
    while (true) {
      const double x = R::rnorm(mu, sigma);
      if (x >= 0.0) return x;
      if ((++iter & 1023L) == 0L) Rcpp::checkUserInterrupt();
    }
  }

  // mu < 0: exponential proposal from Robert (1995), Section 3.  The
  // optimal rate alpha_star = (a + sqrt(a^2 + 4)) / 2 maximises the
  // acceptance probability for the displaced-exponential proposal.
  const double alpha_star = 0.5 * (a + std::sqrt(a * a + 4.0));
  // R::rexp takes the scale (= 1/rate), so pass 1/alpha_star.
  const double inv_rate = 1.0 / alpha_star;
  while (true) {
    const double excess = R::rexp(inv_rate);   // e >= 0, the amount past a
    const double y = a + excess;
    const double diff = y - alpha_star;
    const double log_u = std::log(R::runif(0.0, 1.0));
    if (log_u < -0.5 * diff * diff) {
      // The draw is mu + sigma * y with y = a + e and a = -mu / sigma, so the
      // mu cancels exactly: mu + sigma (a + e) = sigma e.  Returning the
      // literal mu + sigma * y instead subtracts two nearly equal quantities
      // once mu is far below zero, which for gamma near -1e8 produced negative
      // values -- outside the support -- and collapsed 2e5 draws onto 13
      // distinct doubles.
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
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = mhn::rtnorm_robert_zero(mu, sigma);
  }
  return out;
}
