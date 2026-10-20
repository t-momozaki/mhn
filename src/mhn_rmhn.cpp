// rmhn() dispatcher.  Selects between the RTDR sampler (Gao & Wang
// 2025, mhn_rtdr.cpp) and the Sun et al. (2023) Algorithm 1 / 3 paths
// (mhn_sun.cpp), with closed-form shortcuts for the three special
// cases identified in Sun et al. (2023) Lemma 6.  Vectorized over parameter
// inputs with the ParamCacheRmhn re-use pattern; the auto-method decision
// rules are spelled out below.
//
// R-side wrapper: rmhn() in mhn/R/rmhn.R.

#include <Rcpp.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_fma.h"
#include "mhn_param_slots.h"
#include "mhn_special_cases.h"
#include "mhn_rtdr.h"
#include "mhn_stable.h"
#include "mhn_sun.h"
#include "mhn_trunc_normal.h"
#include "mhn_vanishing.h"
#include "mhn_warn.h"

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

  // Proposals a rejection sampler may reject for one draw before it returns
  // NaN (see rmhn_cpp).
  int max_retries = 1000;

  // Special-case scratchpad: only the field matching `kind` is meaningful.
  double sigma = 0.0;     // HALF_NORMAL, TRUNC_NORMAL: 1/sqrt(2 beta)
  double tn_mu = 0.0;     // TRUNC_NORMAL: gamma / (2 beta)
  double tn_mu_lo = 0.0;  // TRUNC_NORMAL: the part of the mean tn_mu leaves out
  double tn_exp_rate = 0.0;  // TRUNC_NORMAL: -gamma where gamma/(2 beta) = -Inf
  // SQRT_GAMMA draws X = sqrt(T) / sg_root with T ~ Gamma(sg_shape, sg_rate):
  // sg_rate = beta and sg_root = 1, or, where T of rate beta would leave the
  // double range, sg_rate = beta / 4^k in [1, 4) and sg_root = 2^k (see
  // build_sampler).
  double sg_shape = 0.0;     // alpha / 2
  double sg_rate = 0.0;
  double sg_scale = 0.0;     // 1 / sg_rate  (R::rgamma takes the scale)
  double sg_root = 0.0;
  double sg_log_root = 0.0;  // log(sg_root)
  double sg_t_min = 0.0;     // see mhn::gamma_precision_floor

  // General samplers (one populated per kind)
  mhn::RtdrEnvelope    rtdr;
  mhn::SunAlgo1Setup   sun_a1;
  mhn::SunAlgo3Setup   sun_a3;

  // A subnormal shape (see rebuild_cache): the probability of a draw below
  // the smallest positive double, and sqrt(beta); -1 for any other shape.
  double vanishing_q = -1.0;
  double vanishing_sqrt_beta = 1.0;
  // The part of the forms that depends on beta and gamma alone, log D among
  // it, for the beta and gamma of the last subnormal shape (if vanishing_set),
  // and whether the sampler above was last built for the proxy shape at that
  // beta and gamma, so that another subnormal shape there reuses both.
  mhn::VanishingShape vanishing;
  bool vanishing_set = false;
  double vanishing_b = 0.0;
  double vanishing_g = 0.0;
  bool proxy_built = false;

  bool needs_rebuild(double a, double b, double g) const {
    // A subnormal shape whose draws are all 0 builds no sampler and leaves
    // kind NONE, but its vanishing_q of 1 marks the setup as made.
    return (kind == Kind::NONE && vanishing_q < 0.0) || a != prev_a ||
           b != prev_b || g != prev_g;
  }
};

// 1 / sqrt(2 beta), the standard deviation of the kernel exp(-beta x^2).
// 2 beta passes the largest double once beta exceeds DBL_MAX / 2, and there
// sqrt(1/2) is divided by sqrt(beta) instead.
double normal_sd(double b) {
  return (b <= 0.5 * std::numeric_limits<double>::max())
           ? 1.0 / std::sqrt(2.0 * b)
           : 0.70710678118654752 / std::sqrt(b);
}

// Whether the square-root-Gamma closed form serves a negligible tilt: it
// draws a Gamma variate of shape alpha / 2 from R::rgamma, which is used only
// up to the shape mhn::kRgammaShapeMax (see mhn_sun.h).
bool sqrt_gamma_closed_form(double a) {
  return !(0.5 * a > mhn::kRgammaShapeMax);
}

// Whether the log-axis mode of a negative tilt is representable.  For
// alpha <= 1 RTDR works on the log axis, y = log(sqrt(beta) X), and the mode
// of that law lies at y = log u, u the positive root of
// 2 u^2 - gamma_norm u - alpha = 0, gamma_norm = gamma / sqrt(beta), about
// alpha / |gamma_norm| under a strong negative tilt.  For a shape below about
// 2.5e-324 |gamma_norm| (half the smallest positive double times
// |gamma_norm|) u rounds to 0, and RTDR cannot place the law (see
// envelope_out_of_range in mhn_rtdr.cpp).  All but a fraction of about
// 1e3 alpha of the mass of X then rounds to 0, and Algorithm 3 samples it.
// A finite |gamma_norm| is below 1.8e308, so u can round to 0 only for a
// shape below 4.5e-16, and only below 1e-15 is it formed.  From that shape
// on an infinite gamma_norm is left to RTDR, which refuses it and names the
// Gamma law (envelope_out_of_range).
bool rtdr_places_mode(double a, double b, double g) {
  return !(a < 1e-15) || mhn::positive_root(1.0, g / std::sqrt(b), a) > 0.0;
}

// Whether the triple is answered by a closed form rather than a rejection
// sampler, with the predicates build_sampler uses.  A subnormal shape, which
// rebuild_cache hands to the mixture or answers with 0, is counted here as
// for its own shape; that only sets how many setups are kept.  They take a shape within
// sqrt(eps) of 1 as 1 and a tilt z = gamma / sqrt(beta) within sqrt(eps) of
// 0 as 0, as d, p and q do not (see is_sqrt_gamma_exact in
// mhn_special_cases.h).  The law drawn there is the special case, and the
// one requested has density ratio 1 + delta (log U - E log U) or
// 1 + z (U - E U) to it at first order, U = sqrt(beta) X, so the two are a
// total variation distance of about |delta| E|log U - E log U| / 2 or
// |z| E|U - E U| / 2 apart.  E|log U - E log U| under the truncated normal
// is at most 0.983, its limit under a strong negative tilt, and E|U - E U|
// under the square root of a Gamma law at most the standard deviation of U,
// below 1/sqrt(2).  Where both are within sqrt(eps) the half-normal is drawn,
// and the ratio is 1 + delta L + z V with L = log U - E log U and
// V = U - E U, a distance of E|delta L + z V| / 2, largest at
// |delta| = |z| = sqrt(eps) of one sign, where E|L + V| = 1.150 under the
// half-normal.  So the distance is below 8.6e-9, and below 7.4e-9 where only
// one of the two is within sqrt(eps): no sample of feasible size can tell the
// laws apart.  The general samplers draw the exact law there too, but of Sun
// et al.'s samplers only Algorithm 2, which is not
// implemented, covers alpha < 1 with gamma > 0, a case method = "sun" would
// then refuse next to alpha = 1 from below under a positive tilt, and next
// to gamma = 0 from above for any alpha < 1.
bool closed_form_triple(double a, double b, double g) {
  return mhn::is_truncated_normal(a) ||
         (mhn::is_sqrt_gamma(g, b) && sqrt_gamma_closed_form(a));
}

// The setup (rebuild_cache) does NOT consume R RNG state. This invariant
// is what makes a vectorised call reproduce the element-wise loop:
// rmhn(n, alpha = v) equals vapply(v, function(a) rmhn(1, a), 0) under one
// seed.  build_rtdr_envelope / build_sun_algo1 / build_sun_algo3 are all
// deterministic closed-form, bisection or Newton computations.
//
// samples_per_setup: estimated number of samples drawn per setup, which
// the auto path's dispatch reads (see rmhn_cpp).  For scalar parameters
// it equals n; for fully vectorised parameters it is 1.  Other paths ignore it.
void build_sampler(ParamCacheRmhn& c, double a, double b, double g,
                   const std::string& method,
                   R_xlen_t samples_per_setup) {
  const bool is_sg = mhn::is_sqrt_gamma(g, b);
  const bool is_tn = mhn::is_truncated_normal(a);

  if (is_sg && is_tn) {
    c.kind = Kind::HALF_NORMAL;
    c.sigma = normal_sd(b);
  } else if (is_sg && !sqrt_gamma_closed_form(a)) {
    // A negligible tilt at a shape above 2^60, where R::rgamma could shift
    // every Gamma variate of shape alpha / 2 by more than 2.5e-7 of its
    // standard deviation (see mhn::kRgammaShapeMax).  The law is then close
    // to normal about its mode, and it is drawn by a rejection sampler that
    // places each draw from the mode computed to about twice double
    // precision: RTDR under method = "rtdr", and otherwise the Normal
    // proposal of Sun et al.'s Algorithm 1 centred at the mode, which is
    // valid for a tilt of either sign (see build_sun_algo1) and accepts about
    // 1/sqrt(2) of its proposals here, with a setup far cheaper than RTDR's.
    if (method == "rtdr") {
      c.kind = Kind::GENERAL_RTDR;
      c.rtdr = mhn::build_rtdr_envelope(a, b, g);
    } else {
      c.kind = Kind::GENERAL_SUN_A1;
      c.sun_a1 = mhn::build_sun_algo1(a, b, g, /*normal_only=*/true);
    }
  } else if (is_sg) {
    c.kind = Kind::SQRT_GAMMA;
    c.sg_shape = a / 2.0;
    // T of rate beta has mean alpha / (2 beta), and R::rgamma takes the scale
    // 1 / beta.  Once either passes 1e300, T can pass the largest double while
    // X = sqrt(T) is far inside the range.  There T is drawn at the rate
    // beta / 4^k, which lies in [1, 4), and its square root divided by 2^k.
    // Scaling by a power of two is exact and commutes with the roundings in
    // R::rgamma and in the square root, so that, but for draws replaced below
    // t_min, the draws are those at rate beta wherever those stay inside the
    // double range.
    if (b >= 1e-300 && c.sg_shape <= 1e300 * b) {
      c.sg_rate = b;
      c.sg_root = 1.0;
    } else {
      const int e = std::ilogb(b);
      const int k = (e - (e < 0)) / 2;  // floor(e / 2)
      c.sg_rate = std::ldexp(b, -2 * k);
      c.sg_root = std::ldexp(1.0, k);
    }
    c.sg_scale = 1.0 / c.sg_rate;
    c.sg_log_root = std::log(c.sg_root);
    c.sg_t_min = mhn::gamma_precision_floor(c.sg_scale);
    // sqrt(alpha / (2 beta)), the root mean square of X, beyond the largest
    // double needs alpha / (2 beta) above 3.2e616 and so alpha above 3e293:
    // the law then lies beyond the largest double, with a standard deviation
    // of about 1 / sqrt(2 alpha) times its mode.
    if (!std::isfinite(std::sqrt(c.sg_shape * c.sg_scale) / c.sg_root)) {
      mhn::stop_beyond_range(a, b, g);
    }
  } else if (is_tn) {
    c.kind = Kind::TRUNC_NORMAL;
    // (gamma / 2) / beta is gamma / (2 beta) rounded once, without forming
    // 2 beta.  Where it overflows, the standard deviation 1 / sqrt(2 beta),
    // at most 3.2e161, is below 1e-146 times |gamma| / (2 beta).  For
    // gamma < 0 the law, N(gamma / (2 beta), 1 / (2 beta)) truncated to
    // x > 0, then lies that many standard deviations to the right of its
    // mean, and is exponential with rate -gamma to double precision; for
    // gamma > 0 it lies beyond the largest double, and the call stops.
    c.tn_mu = (0.5 * g) / b;
    c.sigma = normal_sd(b);
    c.tn_exp_rate = (g < 0.0 && !std::isfinite(c.tn_mu)) ? -g : 0.0;
    if (g > 0.0 && !std::isfinite(c.tn_mu)) mhn::stop_beyond_range(a, b, g);
    // The rounding of the mean, up to half a unit in the last place of tn_mu,
    // moves every draw alike, by up to 2^-53 tn_mu / sigma standard
    // deviations: below 1.2e-7 while tn_mu <= 2^30 sigma, but 0.08 at a tilt
    // gamma / sqrt(beta) of 1e15.  Beyond 2^30 sigma the part of the mean
    // that tn_mu leaves out is kept too, and the draws are placed from the sum
    // (see mhn::rtnorm_robert_zero), as Algorithm 1 places its Normal
    // proposals (mhn_sun.cpp).  The residual 0.5 gamma - beta tn_mu of the
    // rounded quotient is a double, which fma forms exactly; 0.5 gamma is a
    // normal double here, as tn_mu > 2^30 sigma needs
    // gamma > 2^30 sqrt(2 beta), above 3e-153 for every beta.
    c.tn_mu_lo = (c.tn_mu > 1073741824.0 * c.sigma)
                   ? mhn::fma(-c.tn_mu, b, 0.5 * g) / b : 0.0;
  } else if (method == "rtdr") {
    // Where the log-axis mode cannot be placed (see rtdr_places_mode), RTDR
    // cannot be built, and Algorithm 3, which samples the law exactly there,
    // serves instead, as on the default path.
    if (g < 0.0 && !rtdr_places_mode(a, b, g)) {
      c.kind = Kind::GENERAL_SUN_A3;
      c.sun_a3 = mhn::build_sun_algo3(a, b, g);
    } else {
      c.kind = Kind::GENERAL_RTDR;
      c.rtdr = mhn::build_rtdr_envelope(a, b, g);
    }
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
    // theoretical: they are crossovers of the costs of the two samplers,
    // measured on an Apple M1 Ultra by two paired timing races of the two
    // samplers and by three independent runs of
    // inst/benchmarks/auto_dispatch.R.  The choice changes only the speed:
    // every sampler draws the exact law.  The second race times shapes up to 1e14, tilts up to 3 sqrt(alpha) and
    // vectorised calls with one setup per draw; the thresholds from
    // alpha = 1e3 on, and those for a positive tilt above sqrt(alpha), rest
    // on it.  Each sampler costs one setup per parameter triple plus a cost
    // per draw, so which is faster depends on S = samples_per_setup: the
    // lighter setup wins at small S, the cheaper draw at large S.
    //
    // Reaching this branch implies neither special case applies, so
    // |gamma| >= mhn_eps() * sqrt(beta) (intercepted as SQRT_GAMMA above) and
    // |alpha - 1| >= mhn_eps() (intercepted as TRUNC_NORMAL/HALF_NORMAL).
    // Hence the simple sign tests below need no epsilon guards.
    //
    // For alpha > 1 the cost of RTDR's setup depends on where its region-A
    // construction starts the search for its two contact points (see
    // mhn::rtdr_region_a_centred).  On its centred path, from a shape of
    // 2.4e12 to 3.1e12 at moderate tilts, of about 2.3e13 near
    // z = gamma / sqrt(beta) = -alpha/e and of about 4.5e13 at any tilt, and
    // at a tilt z below -1e10 sqrt(alpha) it starts next to the mode.
    // Elsewhere it walks them in from m/2 and m + max(1, m/2), at a cost
    // growing with log2 of the mode m, and from a shape of about 1e3 that walk
    // makes RTDR's setup dearer than that of Sun et al.'s Algorithm 3.
    //
    //   gamma > 0, alpha < 1 -> RTDR.  Sun et al.'s Algorithm 2 is not
    //       implemented in this package, and RTDR's acceptance is at least
    //       1/e (Gao & Wang 2025) where its contact points lie in the band
    //       of their Eq. (7), which the search does not always reach (see
    //       contact_point_newton in mhn_rtdr.cpp; in region D the bound is
    //       also proved only for the full secant chain, see
    //       setup_region_d).
    //   gamma > 0, alpha > 1 -> RTDR if z <= 2 sqrt(alpha) and S reaches
    //       the threshold below; Algorithm 1 otherwise.  Algorithm 1's
    //       closed-form setup is the lighter, so Algorithm 1 wins or ties
    //       S = 1 at every shape and tilt (RTDR takes 1.6 to 13 times as
    //       long in a vectorised call).  Algorithm 1's cost per draw depends
    //       on z / sqrt(alpha): up to 1 (to within 2% from alpha = 3 on) it
    //       proposes from the square root of a Gamma law, an rgamma draw per
    //       proposal, and RTDR's cheaper draw wins from a lower S; above it
    //       Algorithm 1 proposes from a Normal law that accepts 0.81 to 0.92
    //       of its proposals up to z / sqrt(alpha) = 2, RTDR's draw is at
    //       most 14% cheaper (34% at alpha = 1e7), and the crossing comes at
    //       a larger S.  Beyond 2 sqrt(alpha) RTDR gains at most 3.4% at
    //       S = 1000 (at 2.5 sqrt(alpha)), or 13% at alpha = 1e7, where
    //       Algorithm 1's acceptance test is dear.  The thresholds, for
    //       z <= sqrt(alpha) / sqrt(alpha) < z <= 2 sqrt(alpha):
    //         alpha < 10           5 / 25
    //         10 <= alpha < 1e3   25 / 50
    //         1e3 <= alpha < 1e8  50 / 100
    //         alpha >= 1e8       100 / 1000
    //       and 5 / 25 wherever RTDR's setup starts next to the mode.  Each
    //       is a round S near the crossings timed in its band; above it RTDR
    //       is up to 72% faster.
    //   gamma < 0, alpha < 1 -> RTDR if S reaches 25 (alpha < 0.1) or 5;
    //       Algorithm 3 otherwise.  Algorithm 3 draws a Gamma variate and a
    //       fractional power per proposal, so RTDR's draw is the cheaper at
    //       every shape (RTDR is 42-61% faster at S = 1e4), but here
    //       Algorithm 3's setup is the lighter: Algorithm 3 wins S = 1, by
    //       11-23% at alpha = 0.01 and by at most 2% for alpha from 0.3 to
    //       0.8, and at alpha = 0.01 it still wins S = 10, by up to 11%.
    //   gamma < 0, alpha > 1 -> RTDR if its setup starts next to the mode,
    //       or if alpha < 1e3, or if S reaches 10 (alpha < 1e5), 25
    //       (alpha < 1e8) or 50; Algorithm 3 otherwise.  Algorithm 3 locates
    //       its matching point by bisection and Newton steps, and RTDR wins
    //       at every S for alpha from 1.5 to 100 (by up to about 6% at
    //       S = 1; at alpha = 1.1 the two are within 2%), and wherever its
    //       setup starts next to the mode (by 3.5-5% at S = 1, and 16-33% in
    //       a vectorised call).  Where it walks in, from alpha = 1e3,
    //       Algorithm 3 wins S = 1: RTDR takes up to 1.4 times as long as
    //       Algorithm 3 there, and up to 5.2 times in a vectorised call.
    //       RTDR is faster from S = 5 or 10 at alpha = 1e3, 10 or 25 at
    //       1e4, 25 from 1e5 to 1e7 and 50 from 1e8 on, give or take a call
    //       shape in a few cells.
    const double z = g / std::sqrt(b);
    const double sa = std::sqrt(a);
    if (g > 0.0) {
      // No timing covers shapes above 1e14, and RTDR cannot represent the
      // law from a shape of about 1e32 on (see envelope_unrepresentable in
      // mhn_rtdr.cpp), while Algorithm 1 draws it, so from
      // mhn::kRgammaShapeMax on Algorithm 1 serves at any S.  The tilt test is
      // formed as gamma / sqrt(beta) <= 2 sqrt(alpha), not as
      // gamma^2 <= 4 alpha beta, whose sides can overflow; where
      // gamma / sqrt(beta) itself overflows, the test fails and Algorithm 1
      // serves.
      bool use_rtdr = (a < 1.0);
      if (!use_rtdr && a <= mhn::kRgammaShapeMax && z <= 2.0 * sa) {
        const bool mild = (z <= sa);
        R_xlen_t n_switch = (a < 10.0)  ? (mild ? 5 : 25)
                          : (a < 1e3)   ? (mild ? 25 : 50)
                          : (a < 1e8)   ? (mild ? 50 : 100)
                          :               (mild ? 100 : 1000);
        // The centred path needs a shape above 1e12 at every tilt this
        // branch reaches (see mhn::rtdr_region_a_centred), so the test is
        // made only there, and only where it can change the choice.
        if (samples_per_setup < n_switch && a > 1e12 &&
            samples_per_setup >= (mild ? 5 : 25) &&
            mhn::rtdr_region_a_centred(a, z)) {
          n_switch = mild ? 5 : 25;
        }
        use_rtdr = (samples_per_setup >= n_switch);
      }
      if (use_rtdr) {
        c.kind = Kind::GENERAL_RTDR;
        c.rtdr = mhn::build_rtdr_envelope(a, b, g);
      } else {
        c.kind = Kind::GENERAL_SUN_A1;
        c.sun_a1 = mhn::build_sun_algo1(a, b, g);
      }
    } else {  // g < 0.0
      // Three guards take precedence over S.
      //
      // Past a tilt |gamma| / sqrt(beta) of 1e150 Algorithm 3 proposes the
      // Gamma(alpha, rate |gamma|) limit of the law itself and accepts
      // practically every proposal (see build_sun_algo3), at a setup cost no
      // envelope construction can undercut, so it is used there at any S.
      //
      // Above the shape mhn::kRgammaShapeMax Algorithm 3 is not available:
      // R::rgamma could shift its Gamma proposals by more than 2.5e-7 of a
      // standard deviation (see mhn_sun.h).  RTDR, which places its draws from
      // the mode computed to about twice double precision, serves there at
      // any S and tilt.
      //
      // Where the mode of log(sqrt(beta) X), on whose axis RTDR works for
      // alpha <= 1, lies at the log of a number that rounds to 0 (see
      // rtdr_places_mode), RTDR cannot place the law, and Algorithm 3 is
      // used instead.
      const R_xlen_t n_switch = (a < 0.1) ? 25 : (a < 1.0) ? 5
                              : (a < 1e3) ? 1 : (a < 1e5) ? 10
                              : (a < 1e8) ? 25 : 50;
      // Whether RTDR's region-A setup starts next to the mode, where RTDR
      // wins at any S (see mhn::rtdr_region_a_centred).  It is tested only
      // where it can change the choice, and the centred path, which needs a
      // shape above 1e12 at every tilt not already past -1e10 sqrt(alpha),
      // only there.
      const bool quick_start =
          samples_per_setup < n_switch && a > 1.0 &&
          (mhn::rtdr_region_a_strong_negative_tilt(a, z) ||
           (a > 1e12 && mhn::rtdr_region_a_centred(a, z)));
      if (a > mhn::kRgammaShapeMax ||
          ((samples_per_setup >= n_switch || quick_start) &&
           -z <= 1e150 && rtdr_places_mode(a, b, g))) {
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

// The shape from which a draw for a subnormal shape is taken where it is not
// 0 (see rebuild_cache).
constexpr double kVanishingProxyShape = 0x1p-60;

// A subnormal shape alpha is drawn from its law as mhn_vanishing.h gives it,
// a mixture of the origin, with weight q = 1 / (1 + w), and the peak.  A draw
// from the origin's part lies below the smallest positive double (it exceeds
// it with probability below alpha 1117 < 2.5e-305) and is 0.  The peak's
// part has density proportional to (exp(-t^2 + z t) - [t < 1]) / t,
// t = sqrt(beta) x.  A draw comes from it only where q falls below 1 in double
// precision, which needs w above about 5.6e-17 (2^-54) and so z above about
// 52: there its mass below t = 1 is below exp(-600) of it.  On t >= 1 the law
// of the shape 2^-60 has the peak's density times t^(2^-60), which varies by
// under 4e-18 over the peak, z / 2 +- 40.  So a draw is taken from that law,
// by the sampler the dispatch chooses for it, until it reaches t = 1; its own
// mass below t = 1, (2^60 + I(1)) / (2^60 + K0(z)), is below exp(-600)
// there, and the first draw almost always does.  Where q is 1 in double
// precision every draw is 0, and none consumes a random number (q is then 1,
// which sample_one reads as that).  log D and the sampler of the proxy depend
// on beta and gamma alone, and a cache that last served a subnormal shape at
// the same beta and gamma keeps them: vectorised shapes at one tilt cost a
// setup each otherwise, seven times the draws.  The choice
// between 0 and the peak compares one uniform variate with q, and carries
// the resolution of the generator, 2^-32 for R's default: within 2^-32 of
// either end the rarer outcome is never drawn, an error in probability below
// 2.3e-10.  The paths for an ordinary shape do not serve here: they draw 0
// far more often than q.
void rebuild_cache(ParamCacheRmhn& c, double a, double b, double g,
                   const std::string& method,
                   R_xlen_t samples_per_setup) {
  c.vanishing_q = -1.0;
  if (mhn::vanishing_candidate(a)) {
    if (!(c.vanishing_set && b == c.vanishing_b && g == c.vanishing_g)) {
      c.vanishing_set = false;  // until set_tilt has returned
      c.vanishing.set_tilt(b, g);
      c.vanishing_set = true;
      c.vanishing_b = b;
      c.vanishing_g = g;
      c.proxy_built = false;
    }
    mhn::VanishingShape& v = c.vanishing;
    if (v.set_shape(a) || v.lower_only) {
      const double q = std::exp(-v.log1p_w);
      if (q < 1.0 && !c.proxy_built) {
        // A refusal of the law (one narrower than the spacing of the doubles
        // at its mode, say) names the shape asked for, not the one drawn.
        try {
          build_sampler(c, kVanishingProxyShape, b, g, method,
                        samples_per_setup);
        } catch (Rcpp::exception& e) {
          std::string msg = e.what();
          const std::string drawn = tfm::format("%g", kVanishingProxyShape);
          const std::size_t at = msg.find(drawn);
          if (at != std::string::npos) {
            msg.replace(at, drawn.size(), tfm::format("%g", a));
          }
          Rcpp::stop(msg);
        }
        c.proxy_built = true;
      }
      c.vanishing_q = (q < 1.0) ? q : 1.0;
      c.vanishing_sqrt_beta = v.sqrt_beta;
      c.prev_a = a;
      c.prev_b = b;
      c.prev_g = g;
      return;
    }
  }
  c.proxy_built = false;
  build_sampler(c, a, b, g, method, samples_per_setup);
}

double sample_kind(ParamCacheRmhn& c);

double sample_one(ParamCacheRmhn& c) {
  if (c.vanishing_q >= 0.0) {
    if (c.vanishing_q >= 1.0) return 0.0;
    if (unif_rand() < c.vanishing_q) return 0.0;
    // A draw below t = 1 is rejected like a proposal, within the same number
    // of tries, max_retries, that the sampler of the proxy has for each.
    for (int i = 0; i < c.max_retries; ++i) {
      const double x = sample_kind(c);
      if (ISNAN(x) || c.vanishing_sqrt_beta * x >= 1.0) return x;
    }
    return R_NaN;
  }
  return sample_kind(c);
}

double sample_kind(ParamCacheRmhn& c) {
  switch (c.kind) {
    case Kind::HALF_NORMAL:
      return std::abs(R::rnorm(0.0, c.sigma));
    case Kind::SQRT_GAMMA: {
      // X = sqrt(T) with T ~ Gamma(alpha/2, rate beta) (Sun et al. 2023,
      // Lemma 6a), or sqrt(T) / 2^k with T of rate beta / 4^k (see
      // build_sampler).  A T that has lost its precision or underflowed to
      // zero is drawn again from its conditional law below t_min, in
      // logarithms, so that X keeps its digits (see
      // mhn::gamma_precision_floor).
      const double T = R::rgamma(c.sg_shape, c.sg_scale);
      if (T >= c.sg_t_min) return std::sqrt(T) / c.sg_root;
      return std::exp(0.5 * mhn::gamma_tail_log_draw(c.sg_t_min, c.sg_shape,
                                                     c.sg_rate)
                      - c.sg_log_root);
    }
    case Kind::TRUNC_NORMAL:
      if (c.tn_exp_rate > 0.0) return R::exp_rand() / c.tn_exp_rate;
      return mhn::rtnorm_robert_zero(c.tn_mu, c.sigma, c.tn_mu_lo);
    case Kind::GENERAL_RTDR:
      return mhn::sample_rtdr(c.rtdr, nullptr, c.max_retries);
    case Kind::GENERAL_SUN_A1:
      return mhn::sample_sun_algo1(c.sun_a1, nullptr, c.max_retries);
    case Kind::GENERAL_SUN_A3:
      return mhn::sample_sun_algo3(c.sun_a3, nullptr, c.max_retries);
    default:
      Rcpp::stop("rmhn: internal invariant violated (sampler state not "
                 "initialised). Please report this to the package maintainer.");
  }
}

// Pre-scan in method="sun" mode. Reject (alpha<1 && gamma>0) before
// the main loop consumes any RNG state. NA / non-finite gamma elements
// are skipped; they emit NA in the main loop.
//
// The rows the dispatcher answers with a closed form are skipped too, with the
// dispatcher's own predicates.  build_sampler intercepts
// |alpha - 1| < sqrt(eps) as the truncated normal and a negligible
// scale-free tilt as sqrt-Gamma before any Sun algorithm is reached, so
// rmhn(3, 1 - 1e-9, 1, 2, method = "sun") is answered in closed form, as with
// method = "auto"; a test of the bare alpha < 1 && gamma > 0 would refuse it.
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
    if (!mhn::params_valid(a, b, g)) continue;
    if (mhn::is_sqrt_gamma(g, b) || mhn::is_truncated_normal(a)) continue;
    if (a < 1.0 && g > 0.0) {
      Rcpp::stop("Sun method not available for alpha<1 and gamma>0; "
                 "use method=\"rtdr\"");
    }
  }
}

// Whether any of the first `period` recycled triples is answered by a
// rejection sampler rather than a closed form, with the dispatcher's own
// predicates (see build_sampler).  Invalid triples draw NaN and count as
// neither.
bool any_rejection_triple(const Rcpp::NumericVector& alpha,
                          const Rcpp::NumericVector& beta,
                          const Rcpp::NumericVector& gamma,
                          R_xlen_t period) {
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();
  for (R_xlen_t k = 0; k < period; ++k) {
    const double a = alpha[k % na];
    const double b = beta[k % nb];
    const double g = gamma[k % ng];
    if (!mhn::params_valid(a, b, g)) continue;
    if (!closed_form_triple(a, b, g)) return true;
  }
  return false;
}

}  // anonymous namespace

// max_retries is the number of proposals a rejection sampler may reject for
// one draw before it returns that draw as NaN.  rmhn() leaves it at 1000, far
// more than any supported law needs; the argument exists so that the handling
// of an exhausted budget can be exercised.
// [[Rcpp::export(.rmhn_cpp)]]
Rcpp::NumericVector rmhn_cpp(double n,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             std::string method,
                             int max_retries = 1000) {
  // rmhn() passes n as a double: NA for an n that is NA, NULL or a string
  // that is not a number.  A fractional n is truncated, as rnorm() does.
  if (ISNAN(n) || n > static_cast<double>(INT_MAX)) {
    Rcpp::stop("n must be a finite, non-negative number no larger than "
               ".Machine$integer.max");
  }
  if (n < 0.0) Rcpp::stop("n must be non-negative");
  const R_xlen_t nn = static_cast<R_xlen_t>(n);
  if (nn == 0) return Rcpp::NumericVector(0);

  if (method != "auto" && method != "rtdr" && method != "sun") {
    Rcpp::stop("'method' must be one of \"auto\", \"rtdr\", \"sun\"");
  }
  if (max_retries < 0) Rcpp::stop("'max_retries' must be non-negative");

  // A zero-length parameter draws NA, with a warning, as in rnorm().
  if (alpha.size() == 0 || beta.size() == 0 || gamma.size() == 0) {
    mhn::warn("NAs produced");
    return Rcpp::NumericVector(nn, NA_REAL);
  }

  if (method == "sun") {
    prescan_sun_compat(alpha, beta, gamma, nn);
  }

  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  const R_xlen_t L_param = std::max({na, nb, ng});

  // Under recycling the triple at element i depends only on i mod P, with
  // P = lcm(na, nb, ng), so the same P setups serve the whole call.  A single
  // most-recently-used setup would be rebuilt at every element whenever the
  // recycled triples alternate, as in rmhn(n, alpha = c(a1, a2)), and the
  // setup of a rejection sampler can cost as much as many draws.  So when every
  // residue recurs several times (4 P <= n) and some residue needs a
  // rejection sampler, the setups are kept one per residue, as in the density
  // evaluators (mhn_param_slots.h), each created on its first use; otherwise
  // a single setup is kept, as a closed form rebuilds in less time than a
  // kept slot takes to create and reach.  A residue whose triple equals the
  // one just drawn from uses that setup rather than building its own, so
  // parameters grouped into runs of equal triples still cost one setup per
  // run.  Which setup serves a draw cannot change the draw: the setups
  // consume no random numbers (see rebuild_cache).  P is capped at 256, past
  // which recycling_period returns 257 and a single setup is kept.
  const R_xlen_t max_slots = 256;
  const R_xlen_t period_found = mhn::recycling_period(na, nb, ng, max_slots);
  const R_xlen_t period =
    (period_found >= 2 && period_found <= max_slots && 4 * period_found <= nn &&
     any_rejection_triple(alpha, beta, gamma, period_found))
      ? period_found : 1;

  // Samples drawn per setup, which the auto path's dispatch reads.  For
  // scalar parameters it is n.  Where the setups are kept one per residue of
  // the recycling period, each serves n / period draws.  Otherwise, with
  // vectorised parameters, consecutive draws may need different setups, and
  // it is 1, as it is for fully vectorised parameters (max(na, nb, ng) >= n).
  // Runs of equal consecutive triples share a setup but are not counted.
  const R_xlen_t samples_per_setup =
    (L_param == 1) ? nn : (period > 1 ? nn / period : 1);

  Rcpp::NumericVector out(nn);
  R_xlen_t n_exhausted = 0;
  R_xlen_t n_invalid = 0;

  // The caches are scoped so that they are destroyed before the warning below
  // is raised.  Under options(warn = 2) an R warning becomes an error, which
  // with an Rcpp older than 1.0.10 leaves the C++ frame by a long jump (see
  // mhn_warn.h) and would not run the destructors of the vectors held inside
  // the cached RTDR envelopes.
  {
    ParamCacheRmhn single;
    single.max_retries = max_retries;
    // slot_of[k] is the setup of residue k once it exists.  Each residue
    // creates at most one setup, so reserving `period` places keeps `slots`
    // from reallocating and the pointers into it valid; with period 1 nothing
    // is allocated.
    std::vector<ParamCacheRmhn*> slot_of(
      static_cast<std::size_t>(period > 1 ? period : 0), nullptr);
    std::vector<ParamCacheRmhn> slots;
    if (period > 1) slots.reserve(static_cast<std::size_t>(period));
    ParamCacheRmhn* last = nullptr;

    for (R_xlen_t i = 0; i < nn; ++i) {
      const double a = alpha[i % na];
      const double b = beta[i % nb];
      const double g = gamma[i % ng];

      // An element whose parameters are NA, NaN or outside the parameter
      // space is NaN, as in rnorm(), and draws no random numbers.
      if (!mhn::params_valid(a, b, g)) {
        out[i] = R_NaN;
        ++n_invalid;
        continue;
      }

      ParamCacheRmhn* cache = &single;
      if (period == 1) {
        if (single.needs_rebuild(a, b, g)) {
          rebuild_cache(single, a, b, g, method, samples_per_setup);
        }
      } else {
        ParamCacheRmhn*& slot = slot_of[static_cast<std::size_t>(i % period)];
        if (slot != nullptr && !slot->needs_rebuild(a, b, g)) {
          cache = slot;
        } else if (last != nullptr && !last->needs_rebuild(a, b, g)) {
          cache = last;
        } else {
          if (slot == nullptr) {
            slots.emplace_back();
            slot = &slots.back();
            slot->max_retries = max_retries;
          }
          rebuild_cache(*slot, a, b, g, method, samples_per_setup);
          cache = slot;
        }
        last = cache;
      }
      out[i] = sample_one(*cache);
      if ((i & 1023) == 0) Rcpp::checkUserInterrupt();
      // A rejection sampler that exhausts its retry budget returns NaN; no
      // successful draw ever does.  Count them and report once, rather than
      // emitting one warning per failed draw.
      if (ISNAN(out[i])) ++n_exhausted;
    }
  }

  if (n_invalid > 0) mhn::warn("NAs produced");

  if (n_exhausted > 0) {
    mhn::warn("rmhn: the rejection sampler exhausted its retry budget for "
              "%lld of %lld draws, which were returned as NaN.",
              static_cast<long long>(n_exhausted),
              static_cast<long long>(nn));
  }

  return out;
}
