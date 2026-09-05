// rmhn() dispatcher.  Selects between the RTDR sampler (Gao & Wang
// 2025, mhn_rtdr.cpp) and the Sun et al. (2023) Algorithm 1 / 3 paths
// (mhn_sun.cpp), with closed-form shortcuts for the three special
// cases identified in Sun Lemma 6.  Vectorized over parameter inputs
// with the ParamCacheRmhn re-use pattern; the auto-method decision rules
// are spelled out below.
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

  bool needs_rebuild(double a, double b, double g) const {
    return kind == Kind::NONE || a != prev_a || b != prev_b || g != prev_g;
  }
};

// The setup (rebuild_cache) does NOT consume R RNG state. This invariant
// is what makes a vectorised call reproduce the element-wise loop:
// rmhn(n, alpha = v) equals vapply(v, function(a) rmhn(1, a), 0) under one
// seed.  build_rtdr_envelope / build_sun_algo1 / build_sun_algo3 are all
// deterministic Newton or closed-form computations.
//
// samples_per_setup: estimated number of samples drawn per setup, used by
// the auto path's gamma<0 dispatch.  For scalar parameters this equals n;
// for fully vectorised parameters it can be 1.  Other paths ignore it.
void rebuild_cache(ParamCacheRmhn& c, double a, double b, double g,
                   const std::string& method,
                   R_xlen_t samples_per_setup) {
  const bool is_sg = mhn::is_sqrt_gamma(g, b);
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
      Rcpp::stop("rmhn: internal invariant violated (an unsupported alpha < 1 "
                 "with gamma > 0 reached the Sun sampler). Please report this "
                 "to the package maintainer.");
    }
  } else {
    // method == "auto".  The thresholds below are benchmarked rather than
    // theoretical; inst/benchmarks/auto_dispatch.R reproduces the
    // measurement (50 iterations over seven per-call sample counts,
    // 1/5/10/25/50/100/10000).
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
    //                              this package, and RTDR's acceptance
    //                              bound holds uniformly in this region.
    //   gamma < 0              -> n-dependent.  For small n (Gibbs) Sun
    //                              A3's lighter setup wins; for large n
    //                              its heavier per-sample cost (rgamma +
    //                              (X/m)^(1/r)) loses to RTDR's piecewise
    //                              inverse CDF, so the winner flips at
    //                              n ~ 10-25 and we use 25 as the cutoff.
    //                              Exception: for alpha >= 10 Sun A3's
    //                              per-proposal cost falls below RTDR's
    //                              even for large n (RTDR builds a fuller
    //                              envelope as the mode sharpens), so Sun
    //                              A3 wins in both regimes there; the
    //                              alpha >= 10 carve-out keeps large-n
    //                              dispatch on the measured optimum
    //                              (benchmarked in auto_dispatch.R across
    //                              three independent runs: it cuts the
    //                              worst-case regret from ~11% to ~4%).
    if (g > 0.0) {
      if (a > 1.0) {
        c.kind = Kind::GENERAL_SUN_A1;
        c.sun_a1 = mhn::build_sun_algo1(a, b, g);
      } else {
        c.kind = Kind::GENERAL_RTDR;
        c.rtdr = mhn::build_rtdr_envelope(a, b, g);
      }
    } else {  // g < 0.0
      // The Algorithm 3 -> RTDR crossover moves right as the shape shrinks:
      // it sits near n = 25 for alpha around 0.8 and above, but by
      // alpha = 0.01 it has moved out to n ~ 100, because Algorithm 3's
      // setup gets relatively cheaper as the density spikes at the origin.
      // Switching at 25 there would cost up to 10% for 25 <= n < 100, so
      // very small shapes wait longer.
      const R_xlen_t n_switch = (a < 0.1) ? 100 : 25;
      if (samples_per_setup >= n_switch && a < 10.0) {
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
      Rcpp::stop("rmhn: internal invariant violated (sampler state not "
                 "initialised). Please report this to the package maintainer.");
  }
}

// Pre-scan in method="sun" mode. Reject (alpha<1 && gamma>0) before
// the main loop consumes any RNG state. NA / non-finite gamma elements
// are skipped; they emit NA in the main loop.
//
// The rows the dispatcher answers with a closed form are skipped too.  Testing
// the bare floats rejected parameters the sun path handles perfectly well:
// rebuild_cache intercepts |alpha - 1| < sqrt(eps) as the truncated normal and
// a negligible scale-free tilt as sqrt-Gamma before it ever reaches a Sun
// algorithm, so rmhn(3, 1 - 1e-9, 1, 2, method = "sun") raised "not available
// for alpha<1 and gamma>0" for a call that method = "auto" answers in closed
// form.  Mirroring the dispatcher's own predicates keeps the two in step.
void prescan_sun_compat(const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& beta,
                        const Rcpp::NumericVector& gamma,
                        R_xlen_t n) {
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];
    if (Rcpp::NumericVector::is_na(a) || Rcpp::NumericVector::is_na(g)) continue;
    if (!std::isfinite(g)) continue;
    if (Rcpp::NumericVector::is_na(b) || !std::isfinite(b)) continue;
    if (mhn::is_sqrt_gamma(g, b) || mhn::is_truncated_normal(a)) continue;
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
    prescan_sun_compat(alpha, beta, gamma, nn);
  }

  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // Estimate of samples drawn per cache rebuild, used by the auto path's
  // gamma<0 dispatch.  For scalar params (na=nb=ng=1) this
  // equals n; for fully vectorized params (max(na,nb,ng) >= n) it is 1.
  // Partial vectorisation lands somewhere in between.
  const R_xlen_t L_param = std::max({na, nb, ng});
  const R_xlen_t samples_per_setup = (L_param >= nn) ? 1 : (nn / L_param);

  Rcpp::NumericVector out(nn);
  R_xlen_t n_exhausted = 0;

  // The cache is scoped so that it is destroyed before the warning below is
  // raised.  Under options(warn = 2) an R warning becomes an error and leaves
  // the C++ frame by a long jump, which would not run the destructors of the
  // vectors held inside the cached RTDR envelope.
  {
    ParamCacheRmhn cache;

    for (R_xlen_t i = 0; i < nn; ++i) {
      const double a = alpha[i % na];
      const double b = beta[i % nb];
      const double g = gamma[i % ng];

      // The documented contract is an NA draw for any non-finite parameter,
      // not only for gamma: an infinite alpha or beta used to reach the
      // sampler and come back as something the caller could not distinguish
      // from a draw.  R_IsNA covers NA specifically; !isfinite covers NaN and
      // the infinities.
      if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(g)) {
        out[i] = NA_REAL;
        continue;
      }

      if (cache.needs_rebuild(a, b, g)) {
        rebuild_cache(cache, a, b, g, method, samples_per_setup);
      }
      out[i] = sample_one(cache);
      if ((i & 1023) == 0) Rcpp::checkUserInterrupt();
      // A rejection sampler that exhausts its retry budget returns NaN; no
      // successful draw ever does.  Count them and report once, rather than
      // emitting one warning per failed draw.
      // R_IsNA, not Rcpp::NumericVector::is_na: for doubles the latter is
      // ISNAN, so the conjunction was identically false and the count never
      // rose -- rmhn returned all-NaN vectors in complete silence.
      if (ISNAN(out[i]) && !R_IsNA(out[i])) ++n_exhausted;
    }
  }

  if (n_exhausted > 0) {
    Rcpp::warning("rmhn: the rejection sampler exhausted its retry budget for "
                  "%lld of %lld draws, which were returned as NaN.",
                  static_cast<long long>(n_exhausted),
                  static_cast<long long>(nn));
  }

  return out;
}
