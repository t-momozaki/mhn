#ifndef MHN_PMHN_STATE_H
#define MHN_PMHN_STATE_H

// Shared evaluation state for the MHN CDF, used by both pmhn (linear CDF
// evaluation per element) and qmhn (the same routine called repeatedly
// inside the root finder).
//
// Defined inline in a header so that both translation units can link
// without pulling in a dedicated .cpp file.

#include "mhn_cdf_integrate.h"
#include "mhn_cdf_peak.h"
#include "mhn_cdf_series.h"
#include "mhn_constants.h"
#include "mhn_psi.h"
#include "mhn_special_cases.h"
#include "mhn_vanishing.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace mhn {

// The logarithm of the upper-tail quadrature, log_upper_tail in the ordinary
// general case, is a sum: log 2 + (alpha/2) log beta, -log Psi (with
// gamma^2/(4 beta) removed where the pair is shifted), and the logarithm of
// an integral normalised at an anchor, whose own terms, alpha log x at the
// abscissae placed in log x and the tilt terms, reach about
// m (1 + |log v*|), m the larger of alpha and z = gamma / sqrt(beta) and v*
// the peak of the density of log X.  It carries an absolute error of eps
// times the size of those terms or so.  Measured against 160-bit references
// below the median, at 1088 points on 139 triples covering each term up to
// its size at the edge of the ordinary range (about 7e6), the error stays
// below 1.8 eps (1 + |log Psi| + |(alpha/2) log beta| + m (1 + |log v*|));
// kUpperLogErrorScale times that sum is the bound used, together with the
// error of log Psi itself (kLogPsiErrorScale).
const double kUpperLogErrorScale = 4.0;

// log Psi comes from a quadrature, a series or a closed form (see
// mhn_log_normalizing_const), each accurate to a few units in its last
// place, and its error passes whole into log_upper_tail.  Against
// arbitrary-precision references (the power series of Psi, or its
// asymptotic expansion in 1 / z) it is within 4.0 eps (1 + |log Psi|) at
// 23000 points with 0 < z <= 20 and alpha <= 100, where the Lemma 10 series
// serves, and within 1.8 eps at 5000 points with z from 20 to 1e6, or alpha
// from 100 to 1e4 and z up to 20.  Under a negative tilt, where the
// quadrature in log u serves, it is within 4.4 eps (1 + |log Psi|) at
// 1.27e5 points with alpha from 1e-4 to 2000 and z from -1e5 to -1e-4 (see
// psi_log_axis in mhn_psi_integrate.cpp), and at each of them within 0.18
// of upper_log_error, the whole bound.  The bound charges this many times
// eps (1 + |log Psi|), about twice the largest of these, log Psi taken as
// the one used.  At the points of the calibration above the error of
// log Psi is part of the 1.8 eps measured, so there it is charged twice.
const double kLogPsiErrorScale = 8.0;

// Near the lower end of the law log P(X > q) is about -F(q), and the bound
// above is a relative error of that bound divided by F in it.  pmhn and qmhn
// use log_upper_tail there only while that relative error is at most this;
// elsewhere they work from the lower tail.  It is of the order of the
// accuracy of the ordinary forms, about 2e-9 where their terms are largest.
const double kUpperLogRelTol = 1e-9;

enum CdfKind {
  CDF_GENERAL = 0,
  CDF_SQRT_GAMMA = 1,
  CDF_TRUNCATED_NORMAL = 2,
  CDF_VANISHING = 3  // a subnormal alpha (mhn_vanishing.h)
};

struct CdfState {
  CdfKind kind = CDF_GENERAL;
  double alpha = 1.0;
  double beta = 1.0;
  double gamma = 0.0;
  // General-case state:
  double log_psi = 0.0;
  double log_prefactor = 0.0;       // log(2) + (alpha/2) log(beta) - log Psi
  // The same with gamma^2/(4 beta) removed, paired with the shifted integral.
  bool shifted = false;
  double log_prefactor_shifted = 0.0;
  // Whether the lower tail below x_split is tried with the Lemma 1b series
  // first (see log_lower_general).
  bool series_first = false;
  // Whether that series is tried outside |z| <= kCdfSeriesMaxAbsZ, as it is
  // for alpha < 1 at every tilt short of the shift: it needs about z^2 terms
  // there and costs from about 0.1 ms at |z| = 20 to milliseconds from
  // |z| = 100 on (see kCdfSeriesMaxAbsZ and log_lower_resolving).  Where its
  // length passes the ceiling of log_cdf_series (cdf_series_fits), which it
  // does at every x once |z| is about 830 or more, it would be refused
  // whatever x, and neither flag is set.
  bool series_costly = false;
  // An absolute bound on the error of log_upper_tail (see
  // kUpperLogErrorScale); 0 for the special cases and the peak frame, whose
  // log upper tail keeps its relative accuracy near 0, and NaN where the
  // state could not be formed.
  double upper_log_error = 0.0;
  // An estimate of the absolute error of log F from the lower tail, for
  // log_cdf_general: the terms that cancel in it, log Psi,
  // (alpha / 2) log beta and (alpha - 1) log x near x_split, charged as in
  // upper_log_error, whose last term is m (1 + |log v*|) with m the larger
  // of alpha and z, for the abscissae the upper-tail quadrature places in
  // log x; here it is |alpha - 1| (1 + |log v*|).  So this estimate is not
  // above upper_log_error for alpha >= 1/2 or z >= 1 - alpha, and above it
  // only for a smaller shape at a small or negative tilt (8.5 times at
  // (3e-6, 1e-300, -1e-146)), where the lower tail is the one to doubt.
  double lower_log_error = 0.0;
  // Above this point F is taken as 1 - P(X > q) (see cdf_linear).
  double x_split = 0.0;
  // For a shape or a tilt z = gamma / sqrt(beta) above 1e6, and wherever the
  // terms of the forms below pass about 1e7 in size (use_peak_frame), both
  // tails are measured from the peak of the density of log X
  // (mhn_cdf_peak.h), and none of the general-case state above is formed.
  // peak_ok is false where that frame could not be set up, and every
  // evaluation is then NA.
  bool peak_frame = false;
  bool peak_ok = false;
  PeakFrame peak;
  // Truncated-normal state (alpha == 1):
  TruncNormal tn;
  // sqrt-Gamma state (gamma == 0):
  double sg_shape = 0.5;            // alpha / 2
  double sg_sqrt_beta = 1.0;
  // Subnormal-shape state:
  VanishingShape vanishing;
  // A subnormal alpha past the tilt where the forms of mhn_vanishing.h hand
  // over: the general case serves the triple, except for the lower tail far
  // below the peak, where it depends on alpha and the forms still serve it
  // (see log_lower_vanishing).
  bool vanishing_lower = false;

  void recompute(double a, double b, double g) {
    alpha = a; beta = b; gamma = g;
    peak_frame = false;
    upper_log_error = 0.0;
    lower_log_error = 0.0;
    vanishing_lower = false;
    // A subnormal alpha is served by the forms of mhn_vanishing.h wherever
    // they apply, at gamma = 0 too.
    if (vanishing_candidate(a)) {
      if (vanishing.set(a, b, g)) {
        kind = CDF_VANISHING;
        return;
      }
      vanishing_lower = vanishing.lower_only;
    }
    // gamma = 0 above alpha = 1e6 is taken in the peak frame (see
    // SqrtGammaDensity).  The special cases serve gamma = 0 and alpha = 1
    // exactly (see is_sqrt_gamma_exact).
    if (is_sqrt_gamma_exact(g) && !sqrt_gamma_in_peak_frame(a)) {
      kind = CDF_SQRT_GAMMA;
      sg_shape = 0.5 * a;
      sg_sqrt_beta = std::sqrt(b);
      return;
    }
    if (is_truncated_normal_exact(a)) {
      kind = CDF_TRUNCATED_NORMAL;
      tn.set(b, g);
      return;
    }
    kind = CDF_GENERAL;
    peak_frame = use_peak_frame(a, b, g);
    if (peak_frame) {
      shifted = false;
      log_psi = R_NaN;
      log_prefactor = R_NaN;
      log_prefactor_shifted = R_NaN;
      peak_ok = peak.set(a, b, g, mhn_eps());
      x_split = peak_ok ? peak.v : R_NaN;
      return;
    }

    // For a large positive tilt, log Psi and the log of the integral over
    // [0, x] are each about gamma^2/(4 beta) while the log probability they
    // combine into is of order one.  Assembled directly the result carries
    // their rounding: at gamma = 1e7 both are 2.5e13, and pmhn at the mode
    // would come out 1.1e-3 away from 0.5, and, where the quadrature rounds
    // differently, far enough away to be clamped to 1.  Past a threshold
    // well inside the safe range, take both with that term removed.  Only the
    // pair that is used is formed, as in dmhn: the unshifted log Psi is never
    // read when the shifted pair is, and at such a tilt it is the most
    // expensive part of the setup.
    const double z = g / std::sqrt(b);
    shifted = (z > 1.0e3);
    series_first = (a < 1.0 || (std::fabs(z) <= kCdfSeriesMaxAbsZ &&
                                a <= kCdfSeriesMaxAlpha));
    series_costly = !shifted && a < 1.0 && std::fabs(z) > kCdfSeriesMaxAbsZ;
    if (series_costly && !cdf_series_fits(a, b, g, mhn_eps())) {
      series_first = false;
      series_costly = false;
    }
    if (shifted) {
      log_psi = R_NaN;
      log_prefactor = R_NaN;
      log_prefactor_shifted = std::log(2.0) + 0.5 * a * std::log(b)
                              - psi_integrate_shifted(a, b, g, mhn_eps());
    } else {
      log_psi = mhn_log_normalizing_const(a, b, g, -1.0);
      log_prefactor = std::log(2.0) + 0.5 * a * std::log(b) - log_psi;
    }
    // The quadrature for log Psi returns -Inf when it fails.  The prefactor
    // is then +Inf, and set against a finite integral it would read as F = 1
    // (log F = +Inf).  The state is marked NaN instead, so that every
    // evaluation for this triple returns NA.
    if (!std::isfinite(shifted ? log_prefactor_shifted : log_prefactor)) {
      log_psi = R_NaN;
      log_prefactor = R_NaN;
      log_prefactor_shifted = R_NaN;
    }
    // The peak of the density of log X, the positive root of
    // 2 beta t^2 - gamma t - alpha = 0 (see cdf_linear).
    x_split = positive_root(b, g, a);
    // The terms of log_upper_tail (see kUpperLogErrorScale), log Psi taken
    // as the one paired with the integral, and the error of log Psi itself
    // (kLogPsiErrorScale).  NaN where log Psi failed.
    const double half_a_log_b = 0.5 * a * std::log(b);
    const double log_psi_used =
        shifted ? M_LN2 + half_a_log_b - log_prefactor_shifted : log_psi;
    upper_log_error = DBL_EPSILON *
        (kUpperLogErrorScale *
             (1.0 + std::fabs(log_psi_used) + std::fabs(half_a_log_b) +
              std::max(a, z) * (1.0 + std::fabs(std::log(x_split)))) +
         kLogPsiErrorScale * (1.0 + std::fabs(log_psi_used)));
    lower_log_error = DBL_EPSILON *
        (kUpperLogErrorScale *
             (1.0 + std::fabs(log_psi_used) + std::fabs(half_a_log_b) +
              std::fabs(a - 1.0) * (1.0 + std::fabs(std::log(x_split)))) +
         kLogPsiErrorScale * (1.0 + std::fabs(log_psi_used)));
  }

  // Returns F(qi) in [0, 1].  Boundary handling: qi <= 0 -> 0,
  // qi = +Inf -> 1, NA -> NA_REAL.
  double cdf_linear(double qi) const {
    // NA and NaN are distinct in the base R d/p/q contract, and
    // Rcpp::NumericVector::is_na does not separate them -- it is ISNAN.
    if (R_IsNA(qi)) return NA_REAL;
    if (R_IsNaN(qi)) return R_NaN;
    if (qi <= 0.0) return 0.0;
    if (qi == R_PosInf) return 1.0;
    switch (kind) {
      // The special cases have one implementation each, in
      // mhn_special_cases.cpp, shared with the all-scalar paths of pmhn and
      // qmhn.
      case CDF_SQRT_GAMMA:
        return pmhn_sqrt_gamma_one(qi, sg_shape, sg_sqrt_beta, true, false);
      case CDF_TRUNCATED_NORMAL:
        return pmhn_truncated_normal_one(qi, tn, true, false);
      case CDF_GENERAL:
      default: {
        double log_F_far;
        if (log_lower_vanishing(qi, &log_F_far)) return std::exp(log_F_far);
        if (peak_frame) return peak_ok ? peak.cdf(qi) : R_NaN;
        // Each tail is computed directly on its own side of x_split, the peak
        // of the density of log X.  That peak lies at or above the median
        // (the law of log X is skewed to the left), so above it F is
        // 1 - P(X > qi) with an upper tail of at most one half from its own
        // quadrature: nothing cancels, F + P(X > qi) is 1 to rounding, and F
        // rises wherever the upper tail falls.  Taken from the lower tail
        // there, F would approach 1 through the cancelling sum of the
        // gamma < 0 series, whose rounding near F = 1 is of the order of the
        // tolerance.  The upper-tail quadrature is also the cheaper of the two
        // there.
        if (qi > x_split) return -std::expm1(log_upper_tail(qi));
        double F = std::exp(log_lower_general(qi));
        if (F < 0.0) F = 0.0;
        if (F > 1.0) F = 1.0;
        return F;
      }
    }
  }

  // log P(X > qi), computed from the upper tail rather than derived as
  // log(1 - F).  The derived form loses the survival probability entirely once
  // F rounds to 1, which happens as soon as it falls below about 1e-16; base R
  // exposes lower.tail precisely so that its distribution functions do not have
  // that limitation, and this is what lets pmhn and qmhn match.
  double log_upper_tail(double qi) const {
    if (R_IsNA(qi)) return NA_REAL;
    if (R_IsNaN(qi)) return R_NaN;
    if (qi <= 0.0) return 0.0;                 // log 1
    if (qi == R_PosInf) return R_NegInf;
    switch (kind) {
      case CDF_SQRT_GAMMA:
        return pmhn_sqrt_gamma_one(qi, sg_shape, sg_sqrt_beta, false, true);
      case CDF_TRUNCATED_NORMAL:
        return pmhn_truncated_normal_one(qi, tn, false, true);
      case CDF_GENERAL:
      default: {
        if (peak_frame) return peak_ok ? peak.log_upper(qi) : R_NaN;
        const double tol = mhn_eps();
        double log_Q = shifted
            ? log_prefactor_shifted
                  + log_ccdf_integrate_shifted(alpha, beta, gamma, qi, tol)
            : log_prefactor + log_ccdf_integrate(alpha, beta, gamma, qi, tol);
        if (log_Q > 0.0) log_Q = 0.0;
        return log_Q;
      }
    }
  }

  // General case only: log P(X > qi) for finite qi > 0, as pmhn returns it
  // for lower.tail = FALSE, log.p = TRUE, the counterpart of log_cdf_general.
  // It is log_upper_tail, except up to x_split where that upper tail exceeds
  // one half, below the median, and is too close to 1 for its logarithm to
  // keep its relative accuracy.  There log P(X > qi) = log(1 - F) is about
  // -F, while log_upper_tail carries an absolute error of up to
  // upper_log_error, a relative error of that divided by F: at F = 1e-10 it
  // is off by 2e-4 at (3, 1, 10) and by 9 percent at (7730, 1.8e-4,
  // 1.7576), and at F = 0.05 by 2.6e-8 at (2e4, 1e300, 5e151), whose terms
  // reach 7e6.  Wherever that relative error could pass kUpperLogRelTol and
  // F is at most one half, it is log1p(-F), formed from log F without
  // cancellation, which keeps the relative accuracy of F.  Between the
  // median and x_split, where F exceeds one half, 1 - F is the smaller
  // tail, and log1p(-F) would carry the relative error of F multiplied by
  // F / (1 - F), which for a small shape can pass 1e3 just below x_split
  // (F = 0.99914 at (1e-4, 1, -5)).
  //
  // The upper tail is evaluated first, and the lower tail only where the
  // test asks for it (see log_lower_resolving).  Where the lower tail comes
  // from the series it costs more than the upper-tail quadrature, about 4
  // times as much at (2.5, 1, 1.5), and at a triple whose terms are of order
  // one the test asks for it only below F of about 1e-5.  A NaN from the
  // upper tail is passed on, as pmhn passes it on for log.p = FALSE.
  // P(X > qi) itself, on the linear scale, keeps log_upper_tail at every
  // qi: there its absolute error is about that of 1 - F.
  //
  // For a subnormal alpha past the hand-over of mhn_vanishing.h, where F is
  // below one half and depends on alpha, log(1 - F) is formed from the F of
  // its forms, in the peak frame too: log_Q from the ordinary paths is off
  // by up to the mass by the origin, which may be all of F there.
  double log_upper_general(double qi) const {
    const double log_Q = log_upper_tail(qi);
    if (qi > x_split || !(log_Q > -M_LN2)) return log_Q;
    double log_F_far;
    if (log_lower_vanishing(qi, &log_F_far)) {
      return (log_F_far <= -M_LN2) ? mhn::log1m_exp(log_F_far) : log_Q;
    }
    if (peak_frame) return log_Q;
    if (-log_Q * kUpperLogRelTol >= upper_log_error) return log_Q;
    const double log_F = log_lower_resolving_ordinary(qi, log_Q);
    return (log_F <= -M_LN2) ? mhn::log1m_exp(log_F) : log_Q;
  }

  // General case, not the peak frame, qi <= x_split: log F(qi), given
  // log_Q = log_upper_tail(qi), for a caller that needs F where log_Q is too
  // close to 0 to resolve it to kUpperLogRelTol (log_upper_general, and
  // qmhn's polish of an upper-tail root).  It is log_lower_general, except
  // where that would sum the series at a cost of milliseconds
  // (series_costly: 9 ms at (0.3, 1, 300), 25 ms at (0.01, 1, 500)).  There
  // the lower-tail quadrature, which costs microseconds, is taken instead
  // wherever log(1 - F) from it agrees with log_Q to within twice
  // upper_log_error.  Below alpha = 1 that quadrature can miss by more than
  // the working tolerance at scattered points (see kCdfSeriesMaxAbsZ in
  // mhn_cdf_series.h), which is why log_lower_general does not take it
  // first.  A miss that exceeds the error bound of log_Q fails the
  // comparison, and the series is summed after all; one that passes leaves
  // log(1 - F) within three times that bound.  On random triples with
  // alpha < 1 and |z| from 20 to 1000 (to 3e4 under a negative tilt) the
  // quadrature was within 8e-11 of F at the 572 points checked where it is
  // used, against an arbitrary-precision reference; the misses found, up to
  // 1.1e-7, lay at larger F, where log_Q itself is used.
  double log_lower_resolving(double qi, double log_Q) const {
    double log_F_far;
    if (log_lower_vanishing(qi, &log_F_far)) return log_F_far;
    return log_lower_resolving_ordinary(qi, log_Q);
  }

  // A subnormal alpha past the hand-over of mhn_vanishing.h
  // (vanishing_lower): log F at qi from its forms where alpha G(u) is at
  // most 2^64 and F depends on alpha, and false where the general case
  // serves qi.  Every general-case entry to the lower tail asks it first.
  bool log_lower_vanishing(double qi, double* log_F) const {
    return vanishing_lower && vanishing.log_lower_far(qi, log_F);
  }

  // log_lower_resolving for a qi the forms of mhn_vanishing.h do not serve.
  double log_lower_resolving_ordinary(double qi, double log_Q) const {
    if (series_costly) {
      const double log_F =
          log_prefactor + log_cdf_integrate(alpha, beta, gamma, qi, mhn_eps());
      if (std::isfinite(log_F) && log_F <= -M_LN2 &&
          std::fabs(mhn::log1m_exp(log_F) - log_Q) <=
              2.0 * upper_log_error) {
        return log_F;
      }
    }
    return log_lower_general(qi);
  }

  // General case only: log F(qi) for finite qi > 0, as cdf_linear forms it
  // but without leaving the log scale, so that it stays accurate where F
  // itself underflows or rounds to 1.  pmhn returns it for lower.tail = TRUE,
  // log.p = TRUE, and qmhn solves on it when F spans too many decades for a
  // root finder working on F directly, and when the target log probability is
  // below log(DBL_MIN).
  double log_cdf_general(double qi) const {
    double log_F_far;
    if (log_lower_vanishing(qi, &log_F_far)) return log_F_far;
    if (peak_frame) return peak_ok ? peak.log_cdf(qi) : R_NaN;
    if (qi > x_split) {
      return mhn::log1m_exp(log_upper_tail(qi));
    }
    const double log_F = log_lower_general(qi);
    // The counterpart of log_upper_general.  Below x_split F can still
    // exceed one half: between the median and x_split, and for a small shape
    // over most of the range below the peak, where F is 1 less a tail of
    // alpha times a factor that is moderate under a small tilt (it grows
    // like exp(z^2 / 4) under a large one).  log F from the lower tail is the
    // sum of terms of the size of log Psi, (alpha / 2) log beta and
    // (alpha - 1) log x that cancel, and carries an absolute error of up to
    // about lower_log_error, a relative error of
    // that over |log F|, about the upper tail Q: it was 0 for -1.5e-98 at
    // (1e-100, 1, 5), q = 2.5, and 2.2e-6 off at (1e-10, 1, 0.5), q = 0.25.
    // Wherever that could pass kUpperLogRelTol and Q is below one half,
    // log F is log(1 - Q) from the upper tail, the value pmhn gives for
    // lower.tail = FALSE there.  An absolute error d in log Q is then an
    // absolute error of at most Q d / (1 - Q) in log F, and a relative error
    // below 2 d, whatever the size of Q.  The upper tail is taken only where
    // that, with d = upper_log_error and Q / (1 - Q) = exp(-log F) - 1 from
    // the lower tail, is below half lower_log_error, the lower tail's.  Both
    // are estimates of the size of rounding, not measured bounds, and where
    // they are about equal neither tail is clearly the more accurate, and the
    // lower tail stays: at a large shape the term |alpha - 1| log x is the
    // larger part of both, and between the median and x_split, where Q is
    // close to one half, the two tails were each the nearer to a reference
    // at about half the points.  The factor 2 keeps those points; any factor
    // from 1.5 to 4 decides alike on the points scanned, and 1 takes the
    // upper tail at 82 more of them.  The upper tail is evaluated only where
    // it can be taken.  Gated on upper_log_error alone, the switch took the
    // upper tail at a large tilt and an extreme beta, moving log F by
    // 6.7e-10 at (0.3, 1e300, 1e154), q at its peak, where
    // the lower tail agrees with its value at beta = 1 to 6e-14; gated on
    // the first test alone, it took it between the median and x_split of a
    // large shape, where lower_log_error passes about 7e-10 (|log Psi| about
    // 2.6e5 at beta = 1; alpha from 1e3 to 3e3 at beta = 1e300 or 1e-300)
    // and Q is close to one half, without gaining accuracy.  A log Q of -Inf
    // is a failed quadrature (log_ccdf_integrate), not an empty upper tail,
    // whose logarithm is formed finite, and leaves the lower tail in place.
    if (log_F > -M_LN2 && -log_F * kUpperLogRelTol < lower_log_error &&
        2.0 * upper_log_error * std::expm1(-log_F) < lower_log_error) {
      const double log_Q = log_upper_tail(qi);
      if (log_Q > R_NegInf && log_Q < -M_LN2) return mhn::log1m_exp(log_Q);
    }
    // Rounding can put log F slightly above 0; a NaN, from a failed
    // evaluation, is passed on rather than read as F = 1.
    return (log_F > 0.0) ? 0.0 : log_F;
  }

  // General case only: log F(qi) from the lower tail, for finite qi > 0.
  //
  // For alpha < 1, and for |z| <= kCdfSeriesMaxAbsZ with
  // alpha <= kCdfSeriesMaxAlpha (see mhn_cdf_series.h), the Sun et al. (2023)
  // Lemma 1b series is tried first (the paper's CDF formula, which the
  // package truncates at the length Lemma 10(d) gives for the Psi series).
  // The series returns NaN when its alternating-sign cancellation guard
  // fires (gamma < 0), when its last index 2K + 1 exceeds the ceiling of
  // 1e6, or when beta qi^2 is
  // below the smallest normal double; in each case, and outside those
  // ranges, the result is taken from the Boost.Math integration of the
  // unnormalised density.  See src/mhn_cdf_series.cpp for the cancellation
  // analysis.
  double log_lower_general(double qi) const {
    const double tol = mhn_eps();
    if (shifted) {
      // Where the shift is needed the series is not tried.  Its index
      // ceiling would refuse it at every such tilt -- the series runs to
      // index 2K + 1, which grows like 1.44 z^2 and passes the ceiling below
      // z = 900 -- and it returns log_sum - log Psi, which for a large
      // positive tilt are each about gamma^2/(4 beta), so the difference
      // would carry their rounding, as an unshifted quadrature does.
      return log_prefactor_shifted
             + log_cdf_integrate_shifted(alpha, beta, gamma, qi, tol);
    }
    if (series_first) {
      const double log_F =
          log_cdf_series(alpha, beta, gamma, qi, log_psi, tol);
      if (std::isfinite(log_F)) return log_F;
    }
    return log_prefactor + log_cdf_integrate(alpha, beta, gamma, qi, tol);
  }
};

}  // namespace mhn

#endif  // MHN_PMHN_STATE_H
