// The density and both tails of the MHN law measured from the peak of the
// density of log X, for a large shape or a strong positive tilt.
//
// The ordinary assembly of log f(x) and log F(x) adds log 2 + (alpha/2)
// log(beta) - log Psi to the log kernel (alpha - 1) log x - beta x^2 +
// gamma x, and log F adds the log of a quadrature of the kernel normalised at
// an anchor.  Each of these is a sum of terms of size up to about alpha
// log(x) or gamma^2 / (4 beta), and their rounding, eps times that, stays in
// the result, which is of order one: dmhn and pmhn are off by 1e-7 at
// alpha = 1e8, and by more as gamma / sqrt(beta) grows.  (Removing
// gamma^2 / (4 beta) analytically, as the completed square does for
// z = gamma / sqrt(beta) > 1e3, removes the largest of them but leaves the
// others.)  The quadratures, besides, place their abscissae at x itself or at
// log x, which round to the doubles next to the peak, and those are a fixed
// fraction of the law's width apart: about z 2^-53 of a standard deviation
// on the x axis, and z log(z) 2^-53 on the log axis.  Where the law is
// narrower than the spacing of the doubles no abscissa lies inside it.
//
// Here everything is measured from the peak v* of x g(x), g(x) =
// x^(alpha-1) exp(-beta x^2 + gamma x), the positive root of
// 2 beta x^2 - gamma x - alpha = 0, and in the variable w = log(x / v*).
// With u* = sqrt(beta) v*, the stationarity condition removes gamma and
//
//   log(x g(x)) - log(v* g(v*)) = log phi(w)
//     = -alpha (e^w - 1 - w) - u*^2 (e^w - 1)^2
//     = alpha (log(1 + rho) - rho) - beta d^2,   d = x - v*, rho = d / v*,
//
// both terms at most 0 and of the size of the result, so that
//
//   F(x) = int_{-inf}^{w_x} phi / V,   P(X > x) = int_{w_x}^{inf} phi / V,
//   f(x) = exp(log phi(w_x) - log(x / v*)) / (v* V),   V = int phi,
//
// with w_x = log(x / v*).  V is the integral over w of psi_log_axis_relative
// (mhn_psi_integrate.cpp).  The tail of phi beyond x is integrated in the
// offset sigma = w - w_x from x itself, where t = x e^sigma, and with
// Q(t) = 2 beta t^2 - gamma t - alpha, whose negative is the slope of
// log phi,
//
//   log phi(w_x + sigma) - log phi(w_x)
//     = -alpha (e^sigma - 1 - sigma) - Q(x) (e^sigma - 1) - beta (t - x)^2,
//
// each term again at most 0.  Q(x) = 2 d (beta x + alpha / (2 v*)) is formed
// from d, which only adds.  Nothing is rounded onto the doubles near the
// peak: d is (x - v) - s, with v a double next to v* and s = v* - v found
// from the residual of the quadratic at v, formed exactly as a sum of
// products (one Newton step, with its second-order term), and the offsets
// sigma are exact.
//
// The rounding the ordinary assembly leaves grows with the size of its terms,
// about m (1 + |log v*|) with m the larger of alpha and z, and this frame is
// used wherever that passes 1e7, where the ordinary assembly would lose more
// than about 2e-9, and for any alpha > 1e6 or z > 1e6.  At beta = 1 that is
// alpha > 1e6 or z above about 7.3e5; at beta = 1e-300 or 1e300 it is from a
// shape or a tilt of about 3e4.  Below the switch the two agree to about
// 2e-9.

#ifndef MHN_CDF_PEAK_H
#define MHN_CDF_PEAK_H

namespace mhn {

// Above either the density and both tails are taken from PeakFrame.
constexpr double kPeakFrameMinZ = 1.0e6;
constexpr double kPeakFrameMinAlpha = 1.0e6;
// Below these the ordinary assembly is used: the larger m of alpha and z at
// or below kPeakFrameScaleMin, or m (1 + |log v*|) at or below
// kPeakFrameMaxTerm (see use_peak_frame).
constexpr double kPeakFrameScaleMin = 1.0e3;
constexpr double kPeakFrameMaxTerm = 1.0e7;

// Whether (alpha, beta, gamma), in the general case, is taken from PeakFrame.
bool use_peak_frame(double alpha, double beta, double gamma);

// Whether gamma = 0 (exactly, is_sqrt_gamma_exact) is taken from
// PeakFrame rather than from the sqrt-Gamma closed forms.  Those are written
// in the scale-free t = sqrt(beta) x and lose nothing to the size of beta, so
// only the shape decides (see SqrtGammaDensity in mhn_special_cases.h).
inline bool sqrt_gamma_in_peak_frame(double alpha) {
  return alpha > kPeakFrameMinAlpha;
}

struct PeakFrame {
  double alpha = 1.0;
  double beta = 1.0;
  double gamma = 0.0;
  double v = 1.0;       // v* to a few spacings of the doubles, or +Inf
                        //   beyond the range
  double s = 0.0;       // v* - v, below a few spacings of the doubles at v,
  double s_lo = 0.0;    //   as the unrounded pair s + s_lo
  double s_err = 0.0;   // a bound on the error of s + s_lo
  double sd = 1.0;      // v* / sqrt(alpha + 2 beta v*^2), the law's width
  double log_v = 0.0;   // log v*
  double log_V = 0.0;   // log of the integral of phi over the line

  // Returns false where the frame cannot be set up (then nothing else may be
  // called).
  bool set(double a, double b, double g, double tol);

  // x - v*, for finite x > 0.
  double offset(double x) const;

  // log f(x) for finite x > 0.
  double log_density(double x) const;

  // log F(x) (lower = true, for x at or below v*) or log P(X > x)
  // (lower = false, for x at or above v*), for finite x > 0; NaN where the
  // integral could not be formed.
  double log_tail(double x, bool lower) const;

  // Whether an x at offset d from v* lies so close to it that the error of
  // s leaves the density and the tails there unresolved (both are then NA).
  bool unresolved(double d) const;

  // F(x), log P(X > x) and log F(x) for finite x > 0.  Each tail is taken
  // directly on its own side of v*, and the other as its complement: below
  // v* F is at most about one half, and above it P(X > x) is, so 1 minus
  // either loses nothing, and F + P(X > x) = 1 to rounding.
  double cdf(double x) const;
  double log_upper(double x) const;
  double log_cdf(double x) const;

 private:
  // rho = d / v* and log(1 + rho) - rho, for d = offset(x).
  void terms(double x, double d, double* rho, double* lmu) const;
};

}  // namespace mhn

#endif  // MHN_CDF_PEAK_H
