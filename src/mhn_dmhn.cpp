// MHN density f(x | alpha, beta, gamma), computed in log space via the
// Fox-Wright Psi normalizing constant, or, for a shape or a tilt
// gamma / sqrt(beta) above 1e6 or wherever the terms of that assembly pass
// about 1e7 in size, from the peak of the density of log X
// (mhn_cdf_peak.h, use_peak_frame).  Special cases (gamma = 0 ->
// sqrt-Gamma; alpha = 1 -> truncated normal; both -> half-normal) are
// detected and dispatched to mhn_special_cases.cpp for closed-form
// evaluation.  Vectorized over both x and (alpha, beta, gamma) with the
// ParamCache pattern -- mirrored by mhn_pmhn.cpp.
//
// R-side wrapper: dmhn() in mhn/R/dmhn.R.

#include "mhn_cdf_peak.h"
#include "mhn_check.h"
#include "mhn_constants.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_param_slots.h"
#include "mhn_vanishing.h"
#include "mhn_warn.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Special-case kind detected for a given (alpha, gamma) pair.
enum SpecialKind {
  KIND_GENERAL = 0,
  KIND_SQRT_GAMMA = 1,      // gamma == 0
  KIND_TRUNCATED_NORMAL = 2, // alpha == 1
  KIND_VANISHING = 3         // a subnormal alpha (mhn_vanishing.h)
};

// The general-case log density.
//
// Both the scalar fast path and the recycling loop need this.  It is defined
// once, so that the two call shapes cannot drift apart and return different
// values for the same parameters.
struct GeneralDensity {
  double log_const = 0.0;
  double kernel_centre = 0.0;
  // False for the completed square and for the peak frame below, so that the
  // ordinary assembly is reached through a single test.
  bool plain = true;
  bool completed_square = false;
  // For a shape or a tilt above 1e6, and wherever the terms the two forms
  // below add pass about 1e7 in size (mhn_cdf_peak.h, use_peak_frame), the
  // log density is measured from the peak of the density of log X instead:
  // the rounding of those terms would otherwise pass about 2e-9, and the
  // working tolerance further out, 1e-7 at alpha = 1e8.  peak_ok is false
  // where that frame could not be set up, and the density is then NA, with a
  // warning.
  bool peak_frame = false;
  bool peak_ok = false;
  mhn::PeakFrame peak;

  void set(double alpha, double beta, double gamma) {
    peak_frame = mhn::use_peak_frame(alpha, beta, gamma);
    if (peak_frame) {
      plain = false;
      completed_square = false;
      peak_ok = peak.set(alpha, beta, gamma, mhn::mhn_eps());
      return;
    }
    // For a large positive tilt, log Psi and the kernel exponent are each about
    // gamma^2/(4 beta) while the log density is of order one, so assembling the
    // two directly leaves an absolute error of about eps gamma^2/(4 beta).
    // That passes the working tolerance at |z| = gamma/sqrt(beta) around 2e4
    // and reaches 5 percent by 3e7.  Past a threshold well inside the safe
    // range, take log Psi with that term removed and complete the square in the
    // kernel, so it cancels analytically instead of numerically.  Below the
    // threshold the ordinary assembly is used, with log Psi from
    // mhn_log_normalizing_const: the closed form of Sun et al. (2023,
    // Supplementary Lemma 9(c)) at alpha = 2 with gamma >= 0, the series of
    // Supplementary Lemma 10 for 0 < z <= 20 with alpha <= 100, and
    // quadrature otherwise.  The error of the assembly there is at most about
    // eps z^2/4, some 6e-11 at z = 1e3.
    const double z = gamma / std::sqrt(beta);
    completed_square = (z > 1.0e3);
    plain = !completed_square;
    if (completed_square) {
      const double log_nc =
          mhn::psi_integrate_shifted(alpha, beta, gamma, mhn::mhn_eps());
      log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;
      // 2 beta overflows for beta above half the largest double; half of
      // gamma is then divided by beta instead.  The test is made on beta (see
      // opaque_value in mhn_constants.h).
      kernel_centre = (beta <= 0.5 * std::numeric_limits<double>::max())
                          ? gamma / (2.0 * beta)
                          : (0.5 * gamma) / beta;
    } else {
      const double log_nc =
          mhn::mhn_log_normalizing_const(alpha, beta, gamma, -1.0);
      log_const = std::log(2.0) + (alpha / 2.0) * std::log(beta) - log_nc;
      kernel_centre = 0.0;
    }
  }

  double at(double xi, double alpha, double beta, double gamma) const {
    // The density vanishes at infinity, as every base R density does.  Left to
    // the general expression it would be NaN: (alpha-1) log(x) and gamma x
    // are +Inf while -beta x^2 is -Inf.
    if (xi == R_PosInf) return R_NegInf;
    if (xi > 0.0) {
      if (!plain) {
        if (peak_frame) return peak_ok ? peak.log_density(xi) : NA_REAL;
        const double d = xi - kernel_centre;
        return log_const + (alpha - 1.0) * std::log(xi) - beta * d * d;
      }
      return log_const + (alpha - 1.0) * std::log(xi)
             - beta * xi * xi + gamma * xi;
    }
    if (xi == 0.0) {
      if (alpha > 1.0) return R_NegInf;
      if (alpha < 1.0) return R_PosInf;
      return log_const;  // alpha = 1 is taken as the truncated normal
    }
    return R_NegInf;
  }
};

// Per-parameter cached state used by the vectorized loop.
struct ParamCache {
  // Truncated-normal constants, shared with the all-scalar path
  mhn::TruncNormalDensity truncated_normal;
  // sqrt-Gamma constants, shared with the all-scalar path
  mhn::SqrtGammaDensity sqrt_gamma;
  SpecialKind kind = KIND_GENERAL;
  GeneralDensity general;
  mhn::VanishingShape vanishing;

  void recompute(double alpha, double beta, double gamma) {
    // A subnormal alpha is served by the forms of mhn_vanishing.h wherever
    // they apply, at gamma = 0 too.  gamma = 0 above alpha = 1e6 is taken in
    // the peak frame (see mhn::SqrtGammaDensity).  The special cases serve
    // gamma = 0 and alpha = 1 exactly (see mhn::is_sqrt_gamma_exact).
    if (mhn::vanishing_candidate(alpha) &&
        vanishing.set(alpha, beta, gamma)) {
      kind = KIND_VANISHING;
    } else if (mhn::is_sqrt_gamma_exact(gamma) &&
               !mhn::sqrt_gamma_in_peak_frame(alpha)) {
      kind = KIND_SQRT_GAMMA;
      sqrt_gamma.set(alpha, beta);
    } else if (mhn::is_truncated_normal_exact(alpha)) {
      kind = KIND_TRUNCATED_NORMAL;
      truncated_normal.set(beta, gamma);
    } else {
      kind = KIND_GENERAL;
      general.set(alpha, beta, gamma);
    }
  }

  double log_density_at(double xi, double alpha, double beta, double gamma) const {
    // NA and NaN are distinct in the base R d/p/q contract, and Rcpp's
    // is_na does not separate them -- it is ISNAN -- so each is tested
    // with the R_Is* predicate that isolates it.
    if (R_IsNA(xi)) return NA_REAL;
    if (R_IsNaN(xi)) return R_NaN;
    switch (kind) {
      case KIND_SQRT_GAMMA:
        return sqrt_gamma.log_at(xi);
      case KIND_TRUNCATED_NORMAL:
        return truncated_normal.log_at(xi);
      case KIND_VANISHING:
        return vanishing.log_density(xi);
      case KIND_GENERAL:
      default:
        return general.at(xi, alpha, beta, gamma);
    }
  }
};

// Warns of a density that could not be evaluated (an NA for an x that is not
// NA or NaN), once per call, naming the first such triple.
void report_failure(double alpha, double beta, double gamma) {
  mhn::warn("dmhn: the density could not be evaluated for alpha=%g, "
            "beta=%g, gamma=%g; returning NA", alpha, beta, gamma);
}

// Fast path when all of alpha, beta, gamma are scalars.
// Kept separate from the recycling loop because this is the dominant call
// shape in practice: MCMC density evaluations at fixed parameters, where
// the modular indexing and the cache check would be pure overhead.
Rcpp::NumericVector dmhn_scalar_path(const Rcpp::NumericVector& x,
                                     double alpha, double beta, double gamma,
                                     bool log_p) {
  if (mhn::vanishing_candidate(alpha)) {
    mhn::VanishingShape vanishing;
    if (vanishing.set(alpha, beta, gamma)) {
      const R_xlen_t n = x.size();
      Rcpp::NumericVector out(n);
      for (R_xlen_t i = 0; i < n; ++i) {
        const double xi = x[i];
        if (ISNAN(xi)) {
          out[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
        } else {
          const double l = vanishing.log_density(xi);
          out[i] = log_p ? l : std::exp(l);
        }
      }
      return out;
    }
  }
  if (mhn::is_sqrt_gamma_exact(gamma) &&
      !mhn::sqrt_gamma_in_peak_frame(alpha)) {
    return mhn::dmhn_sqrt_gamma(x, alpha, beta, log_p);
  }
  if (mhn::is_truncated_normal_exact(alpha)) {
    return mhn::dmhn_truncated_normal(x, beta, gamma, log_p);
  }

  GeneralDensity general;
  general.set(alpha, beta, gamma);

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  if (general.completed_square) {
    // The completed square, with its constants in locals; the same
    // arithmetic as GeneralDensity::at, which the recycled path shares.
    const double log_const = general.log_const;
    const double centre = general.kernel_centre;
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (ISNAN(xi)) {
        log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
      } else if (xi > 0.0 && xi < R_PosInf) {
        const double d = xi - centre;
        log_f[i] = log_const + (alpha - 1.0) * std::log(xi) - beta * d * d;
      } else {
        log_f[i] = general.at(xi, alpha, beta, gamma);
      }
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (ISNAN(xi)) {
        log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
      } else {
        log_f[i] = general.at(xi, alpha, beta, gamma);
      }
    }
  }
  // Only the peak frame gives NA for an x that is not NA or NaN: at every x
  // where it could not be set up, and next to a peak it cannot resolve
  // (PeakFrame::unresolved).  Each such NA is looked for, as in the recycled
  // path below.
  if (general.peak_frame) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (R_IsNA(log_f[i]) && !ISNAN(x[i])) {
        report_failure(alpha, beta, gamma);
        break;
      }
    }
  }

  if (!log_p) mhn::exp_in_place(log_f);
  return log_f;
}

}  // namespace

// [[Rcpp::export(.dmhn_cpp, rng = false)]]
Rcpp::NumericVector dmhn_cpp(Rcpp::NumericVector x,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector beta,
                             Rcpp::NumericVector gamma,
                             bool log_p) {
  const R_xlen_t nx = x.size();
  const R_xlen_t na = alpha.size();
  const R_xlen_t nb = beta.size();
  const R_xlen_t ng = gamma.size();

  // A zero-length argument gives numeric(0), as in dnorm.
  if (nx == 0 || na == 0 || nb == 0 || ng == 0) return Rcpp::NumericVector(0);

  // An element whose parameters lie outside the parameter space is NaN, and
  // one with an NA or NaN parameter NA or NaN (mhn::invalid_result); the R
  // wrapper warns of the NaN.  When every element is valid, the loops below
  // carry no test for it.
  const bool all_valid = mhn::all_params_valid(alpha, beta, gamma);

  // Fast path: all parameters scalar.
  if (all_valid && na == 1 && nb == 1 && ng == 1) {
    return dmhn_scalar_path(x, alpha[0], beta[0], gamma[0], log_p);
  }

  // Vectorized path. Recycle each input to length n = max(nx, na, nb, ng).
  const R_xlen_t n = std::max({nx, na, nb, ng});
  Rcpp::NumericVector log_f(n);

  // One cache slot per distinct triple in the recycling cycle, so that
  // log Psi is computed once per triple rather than once per element.
  mhn::ParamSlots<ParamCache> slots(na, nb, ng, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i % nx];
    const double a = alpha[i % na];
    const double b = beta[i % nb];
    const double g = gamma[i % ng];

    if (!all_valid && !mhn::params_valid(a, b, g)) {
      log_f[i] = mhn::invalid_result(xi, a, b, g);
      continue;
    }
    log_f[i] = slots.at(i, a, b, g).log_density_at(xi, a, b, g);
    if ((i & 255) == 0) Rcpp::checkUserInterrupt();
  }
  // An NA for an x that is not NA or NaN, at valid parameters, is a density
  // that could not be formed.  Looked for in a pass of its own, so that the
  // loop above carries no test, and first as any NaN at all, a reduction
  // without branches.
  const double* value = log_f.begin();
  bool any_nan = false;
  for (R_xlen_t i = 0; i < n; ++i) any_nan |= std::isnan(value[i]);
  if (any_nan) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double a = alpha[i % na];
      const double b = beta[i % nb];
      const double g = gamma[i % ng];
      if (std::isnan(value[i]) && R_IsNA(value[i]) && !ISNAN(x[i % nx]) &&
          mhn::params_valid(a, b, g)) {
        report_failure(a, b, g);
        break;
      }
    }
  }

  if (!log_p) mhn::exp_in_place(log_f);
  return log_f;
}
