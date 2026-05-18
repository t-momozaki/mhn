// One-sided truncated-normal sampling on [0, infinity).  Used by the
// alpha = 1 special-case branch of rmhn (Sun et al. 2023, Lemma 6b)
// and by the rare alpha = 1 fallback inside the RTDR setup.

#ifndef MHN_TRUNC_NORMAL_H
#define MHN_TRUNC_NORMAL_H

#include <Rcpp.h>

namespace mhn {

// Sample one variate from N(mu, sigma^2) truncated to [0, infinity).
// Implements Robert (1995): rejection from the untruncated normal when
// mu >= 0, exponential proposal otherwise.  Caller must guarantee
// sigma > 0 and mu finite.  Consumes R RNG state.
double rtnorm_robert_zero(double mu, double sigma);

}  // namespace mhn

#endif  // MHN_TRUNC_NORMAL_H
