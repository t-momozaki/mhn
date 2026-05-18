// rmhn() dispatcher.  Selects between the RTDR sampler (Gao & Wang
// 2025, mhn_rtdr.cpp) and the Sun et al. (2023) Algorithm 1 / 3 paths
// (mhn_sun.cpp), with closed-form shortcuts for the three special
// cases identified in Sun Lemma 6.  Vectorized over parameter inputs
// with the ParamCacheRmhn re-use pattern; auto-method decision rules
// follow spec §8.7.3 and are spelled out below.
//
// R-side wrapper: rmhn() in mhn/R/rmhn.R.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_special_cases.h"
#include "mhn_rtdr.h"
#include "mhn_sun.h"
#include "mhn_trunc_normal.h"

namespace {

enum class Kind {
  NONE,
  HALF_NORMAL,
  SQRT_GAMMA,
  TRUNC_NORMAL,
  GENERAL_RTDR,
  GENERAL_SUN_A1,
  GENERAL_SUN_A3
};

struct ParamCacheRmhn {
  Kind kind = Kind::NONE;
  double prev_a = std::numeric_limits<double>::quiet_NaN();
  double prev_b = std::numeric_limits<double>::quiet_NaN();
  double prev_g = std::numeric_limits<double>::quiet_NaN();

  // Special-case scratchpad: only the field matching `kind` is meaningful.
  double sigma = 0.0;     // HALF_NORMAL, TRUNC_NORMAL: 1/sqrt(2 beta)
  double tn_mu = 0.0;     // TRUNC_NORMAL: gamma / (2 beta)
  double sg_shape = 0.0;  // SQRT_GAMMA: alpha / 2
  double sg_scale = 0.0;  // SQRT_GAMMA: 1 / beta  (R::rgamma takes scale)

  // General samplers (one populated per kind)
  mhn::RtdrEnvelope    rtdr;
  mhn::SunAlgo1Setup   sun_a1;
  mhn::SunAlgo3Setup   sun_a3;

  // Reserved for Step 3.7 diagnostic instrumentation.
  int last_retries = 0;

  bool needs_rebuild(double a, double b, double g) const {
    return kind == Kind::NONE || a != prev_a || b != prev_b || g != prev_g;
  }
};

// The setup (rebuild_cache) does NOT consume R RNG state. This invariant
// is critical for the vectorization-reproducibility test (spec Sec.12.7),
// where rmhn(n, alpha=v) must equal vapply(v, function(a) rmhn(1, a)).
// build_rtdr_envelope / build_sun_algo1 / build_sun_algo3 are all
// deterministic Newton/closed-form computations (Step 3.1-3.3).
//
// samples_per_setup: estimated number of samples drawn per setup, used by
// the auto path's gamma<0 dispatch (Step 3.7.3).  For scalar params this
// equals n; for fully vectorized params it can be 1.  Other paths ignore
// this parameter.
void rebuild_cache(ParamCacheRmhn& c, double a, double b, double g,
                   const std::string& method,
                   R_xlen_t samples_per_setup) {
  const bool is_sg = mhn::is_sqrt_gamma(g);
  const bool is_tn = mhn::is_truncated_normal(a);

  if (is_sg && is_tn) {
    c.kind = Kind::HALF_NORMAL;
    c.sigma = 1.0 / std::sqrt(2.0 * b);
  } else if (is_sg) {
    c.kind = Kind::SQRT_GAMMA;
    c.sg_shape = a / 2.0;
    c.sg_scale = 1.0 / b;
  } else if (is_tn) {
    c.kind = Kind::TRUNC_NORMAL;
    c.tn_mu = g / (2.0 * b);
    c.sigma = 1.0 / std::sqrt(2.0 * b);
  } else if (method == "rtdr") {
    c.kind = Kind::GENERAL_RTDR;
    c.rtdr = mhn::build_rtdr_envelope(a, b, g);
  } else if (method == "sun") {
    if (g > 0.0 && a > 1.0) {
      c.kind = Kind::GENERAL_SUN_A1;
      c.sun_a1 = mhn::build_sun_algo1(a, b, g);
    } else if (g <= 0.0) {
      c.kind = Kind::GENERAL_SUN_A3;
      c.sun_a3 = mhn::build_sun_algo3(a, b, g);
    } else {
      // alpha < 1 && gamma > 0 with method="sun": should have been
      // caught by prescan_sun_compat. Defensive stop.
      Rcpp::stop("internal: method=\"sun\" with alpha<1 and gamma>0 "
                 "should have been caught by pre-scan");
    }
  } else {
    // method == "auto".  Decision rules from spec Sec 8.7.3, benchmarked
    // in mhn/inst/benchmarks/auto_dispatch.R (iter=50, seven n_per_call
    // patterns 1/5/10/25/50/100/10000).
    //
    // Reaching this branch implies neither special case applies, so
    // |gamma| >= MHN_EPS (intercepted as SQRT_GAMMA above) and
    // |alpha - 1| >= MHN_EPS (intercepted as TRUNC_NORMAL/HALF_NORMAL).
    // Hence the simple sign tests below need no epsilon guards.
    //
    //   gamma > 0 && alpha > 1 -> Sun A1.  Closed-form setup and
    //                              Sun et al. (2023) Theorem 2e gives
    //                              acceptance >= 0.8 for alpha >= 4;
    //                              wins uniformly across all n.
    //   gamma > 0 && alpha < 1 -> RTDR.  Sun A2 is not implemented in
    //                              this package; spec Sec 3.4 [2]
    //                              confirms this region for RTDR.
    //   gamma < 0              -> n-dependent.  The per-cell winner
    //                              flips from Sun A3 (small n) to RTDR
    //                              (large n) at n_per_call ~ 10-25.
    //                              Sun A3 has lighter setup but a
    //                              heavier per-sample cost (rgamma +
    //                              (X/m)^(1/r)) than RTDR's piecewise
    //                              inverse CDF.  We use 25 as the
    //                              cutoff.
    if (g > 0.0) {
      if (a > 1.0) {
        c.kind = Kind::GENERAL_SUN_A1;
        c.sun_a1 = mhn::build_sun_algo1(a, b, g);
      } else {
        c.kind = Kind::GENERAL_RTDR;
        c.rtdr = mhn::build_rtdr_envelope(a, b, g);
      }
    } else {  // g < 0.0
      if (samples_per_setup >= 25) {
        c.kind = Kind::GENERAL_RTDR;
        c.rtdr = mhn::build_rtdr_envelope(a, b, g);
      } else {
        c.kind = Kind::GENERAL_SUN_A3;
        c.sun_a3 = mhn::build_sun_algo3(a, b, g);
      }
    }
  }

  c.prev_a = a;
  c.prev_b = b;
  c.prev_g = g;
}

double sample_one(ParamCacheRmhn& c) {
  switch (c.kind) {
    case Kind::HALF_NORMAL:
      return std::abs(R::rnorm(0.0, c.sigma));
    case Kind::SQRT_GAMMA:
      return std::sqrt(R::rgamma(c.sg_shape, c.sg_scale));
    case Kind::TRUNC_NORMAL:
      return mhn::rtnorm_robert_zero(c.tn_mu, c.sigma);
    case Kind::GENERAL_RTDR:
      return mhn::sample_rtdr(c.rtdr, nullptr);
    case Kind::GENERAL_SUN_A1:
      return mhn::sample_sun_algo1(c.sun_a1, nullptr);
    case Kind::GENERAL_SUN_A3:
      return mhn::sample_sun_algo3(c.sun_a3, nullptr);
    default:
      Rcpp::stop("internal: ParamCacheRmhn not initialized");
  }
}

// Pre-scan in method="sun" mode. Reject (alpha<1 && gamma>0) before
// the main loop consumes any RNG state. NA / non-finite gamma elements
// are skipped; they emit NA in the main loop.
void prescan_sun_compat(const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& gamma,
                        R_xlen_t n) {
  const R_xlen_t na = alpha.size();
  const R_xlen_t ng = gamma.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double a = alpha[i % na];
    const double g = gamma[i % ng];
    if (Rcpp::NumericVector::is_na(a) || Rcpp::NumericVector::is_na(g)) continue;
    if (!std::isfinite(g)) continue;
    if (a < 1.0 && g > 0.0) {
      Rcpp::stop("Sun method not available for alpha<1 and gamma>0; "
                 "use method=\"rtdr\"");
    }
  }
}

}  // anonymous namespace

// [[Rcpp::export(.rmhn_cpp)]]
Rcpp::NumericVector rmhn_cpp(int n,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             std::string method) {
  if (n < 0) Rcpp::stop("n must be non-negative");
  if (n == 0) return Rcpp::NumericVector(0);

  if (method != "auto" && method != "rtdr" && method != "sun") {
    Rcpp::stop("'method' must be one of \"auto\", \"rtdr\", \"sun\"");
  }

  mhn::check_params_vector_allow_na(alpha, beta, gamma);

  const R_xlen_t nn = static_cast<R_xlen_t>(n);

  if (method == "sun") {
    prescan_sun_compat(alpha, gamma, nn);
  }

  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // Estimate of samples drawn per cache rebuild, used by the auto path's
  // gamma<0 dispatch (Step 3.7.3).  For scalar params (na=nb=ng=1) this
  // equals n; for fully vectorized params (max(na,nb,ng) >= n) it is 1.
  // Partial vectorisation lands somewhere in between.
  const R_xlen_t L_param = std::max({na, nb, ng});
  const R_xlen_t samples_per_setup = (L_param >= nn) ? 1 : (nn / L_param);

  Rcpp::NumericVector out(nn);
  ParamCacheRmhn cache;

  for (R_xlen_t i = 0; i < nn; ++i) {
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (Rcpp::NumericVector::is_na(a) ||
        Rcpp::NumericVector::is_na(b) ||
        Rcpp::NumericVector::is_na(g) ||
        !std::isfinite(g)) {
      out[i] = NA_REAL;
      continue;
    }

    if (cache.needs_rebuild(a, b, g)) {
      rebuild_cache(cache, a, b, g, method, samples_per_setup);
    }
    out[i] = sample_one(cache);
  }

  return out;
}
