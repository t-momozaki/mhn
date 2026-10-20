// One-sided truncated-normal sampling on [0, infinity).  Used by the
// alpha = 1 special-case branch of rmhn (Sun et al. 2023, Lemma 6b).

#ifndef MHN_TRUNC_NORMAL_H
#define MHN_TRUNC_NORMAL_H

#include <Rcpp.h>

namespace mhn {

// Sample one variate from N(mu + mu_lo, sigma^2) truncated to [0, infinity).
// Implements Robert (1995): rejection from the untruncated normal when
// mu >= 0, exponential proposal otherwise.  mu_lo, if not 0, is the part of
// the mean that the double mu leaves out, at most half a unit in its last
// place; a draw from the untruncated normal is then placed from mu + mu_lo
// (see mhn_trunc_normal.cpp).  Caller must guarantee sigma > 0, mu and
// mu / sigma finite, and mu_lo = 0 whenever mu < 0.  Consumes R RNG state.
double rtnorm_robert_zero(double mu, double sigma, double mu_lo = 0.0);

}  // namespace mhn

#endif  // MHN_TRUNC_NORMAL_H
