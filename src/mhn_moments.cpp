// Moment functions for the MHN distribution.  All formulas come from
// Sun et al. (2023) Lemma 2 (mean, variance) and the recurrence in
// Lemma 2b (higher raw moments).

#include "mhn_check.h"
#include "mhn_psi.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <vector>

namespace {

// Raw moments E(X^0..X^k_max) via the Sun et al. (2023) Lemma 2b
// recurrence.
std::vector<double> mhn_raw_moments(double alpha, double beta, double gamma,
                                    int k_max) {
  const double log_psi_num = mhn::mhn_log_normalizing_const(alpha + 1.0,
                                                            beta, gamma, -1.0);
  const double log_psi_den = mhn::mhn_log_normalizing_const(alpha, beta,
                                                            gamma, -1.0);
  const double mu = std::exp(log_psi_num - 0.5 * std::log(beta) - log_psi_den);

  std::vector<double> moments(static_cast<size_t>(k_max + 1), 0.0);
  moments[0] = 1.0;
  if (k_max >= 1) moments[1] = mu;

  for (int k = 0; k <= k_max - 2; ++k) {
    moments[static_cast<size_t>(k + 2)] =
        (alpha + k) / (2.0 * beta) * moments[static_cast<size_t>(k)]
        + gamma / (2.0 * beta) * moments[static_cast<size_t>(k + 1)];
  }
  return moments;
}

// Central moments E[(X - mu)^k] for k = 2, 3, 4, by quadrature of the
// unnormalised kernel.
//
// Forming these from raw moments -- m3 - 3 mu m2 + 2 mu^3, and the fourth-order
// analogue -- cancels catastrophically as soon as sigma is small next to mu,
// which is exactly what a large tilt produces.  At (alpha, beta, gamma) =
// (3, 1, 1000) the mean is about 500 and the standard deviation about 0.7, so
// mu^4 is 6e10 while the fourth central moment is under 1: every significant
// digit is lost, and the reported excess kurtosis came out at -701 for a
// quantity that cannot fall below -2.
//
// Quadrature of (x - mu)^k against the kernel has no such subtraction.  The
// kernel is peak-normalised so it never overflows, and the normalising
// constant cancels between numerator and denominator, so Psi is not needed.
struct CentralMoments {
  double mu, c2, c3, c4;
  bool ok;
};

CentralMoments central_moments(double alpha, double beta, double gamma) {
  CentralMoments out{0.0, 0.0, 0.0, 0.0, false};

  // Interior peak of the kernel, where one exists.
  double m = 0.0;
  if (alpha > 1.0) {
    m = mhn::positive_root(beta, gamma, alpha - 1.0);
  } else if (gamma > 0.0) {
    const double disc = gamma * gamma + 8.0 * beta * (alpha - 1.0);
    if (disc > 0.0) m = (gamma + std::sqrt(disc)) / (4.0 * beta);
  }
  if (!(m > 0.0) || !std::isfinite(m)) {
    // No interior peak: alpha <= 1 with a non-positive tilt, where the kernel
    // decreases from an unbounded value at the origin.  The law tends to
    // Gamma(alpha, |gamma|) as the tilt grows, so its scale is alpha/|gamma|
    // rather than anything Gaussian, and the raw-moment expansion this used to
    // fall back on is exactly where it fails worst: excess kurtosis at
    // (1, 1, -3000) came out as -58485 for a true value of 6.
    const double g_abs = std::fabs(gamma);
    const double ref = alpha / (g_abs + std::sqrt(2.0 * beta * alpha));
    if (!(ref > 0.0) || !std::isfinite(ref)) return out;
    const double Ldrop = 45.0;
    // Distance from the origin at which beta R^2 + |gamma| R reaches Ldrop.
    const double R = 2.0 * Ldrop /
        (g_abs + std::sqrt(g_abs * g_abs + 4.0 * beta * Ldrop));
    const double hi0 = std::max(R, 40.0 * ref);
    const double peak0 = (alpha - 1.0) * std::log(ref) - beta * ref * ref + gamma * ref;
    auto w0 = [alpha, beta, gamma, peak0](double x) -> double {
      if (x <= 0.0) return 0.0;
      const double e = (alpha - 1.0) * std::log(x) - beta * x * x + gamma * x - peak0;
      return (e < -700.0) ? 0.0 : std::exp(e);
    };
    const double tol0 = std::sqrt(std::numeric_limits<double>::epsilon());
    boost::math::quadrature::tanh_sinh<double> ts0;
    auto I0 = [&](const std::function<double(double)>& f) {
      return ts0.integrate(f, 0.0, hi0, tol0);
    };
    // The normalising integral is the one that reaches below the smallest
    // representable double for small alpha; the moment integrands carry an
    // extra factor of x and do not.  Under-computing it inflated every central
    // moment by the shortfall -- the variance at (1e-6, 1, -1) came out 1582
    // times too large.
    const mhn::KernelHead head0 =
        mhn::kernel_head(alpha, beta, gamma, hi0, peak0);
    double n0 = head0.value;
    if (head0.x0 < hi0) {
      n0 += ts0.integrate([&](double x) { return w0(x); }, head0.x0, hi0, tol0);
    }
    if (!(n0 > 0.0) || !std::isfinite(n0)) return out;
    const double mu0 = I0([&](double x) { return x * w0(x); }) / n0;
    out.mu = mu0;
    out.c2 = I0([&](double x) { const double d = x - mu0; return d * d * w0(x); }) / n0;
    out.c3 = I0([&](double x) { const double d = x - mu0; return d * d * d * w0(x); }) / n0;
    out.c4 = I0([&](double x) { const double d = x - mu0; const double d2 = d * d;
                                return d2 * d2 * w0(x); }) / n0;
    out.ok = std::isfinite(out.mu) && out.mu > 0.0 &&
             std::isfinite(out.c2) && out.c2 > 0.0 &&
             std::isfinite(out.c3) && std::isfinite(out.c4);
    return out;
  }

  // Work in s = x - m.  Using the stationarity condition at the peak,
  // gamma - 2 beta m = -(alpha-1)/m, the log kernel relative to its peak is
  //
  //     L(m + s) - L(m) = (alpha-1) [ log(1 + s/m) - s/m ] - beta s^2,
  //
  // in which nothing large appears.  Evaluated as written, L(x) - L(m)
  // subtracts beta x^2 and gamma x, both of order 1e13 at gamma = 1e7, to
  // leave a quantity of order 1 -- so the exponent carried an absolute error
  // of 5e-3 and the mean came out 556 adrift.
  //
  // The reach must follow the decay the kernel actually has.  Taking it from
  // the curvature alone assumes a Gaussian fall-off, but a negative tilt
  // squeezes the peak towards the origin and the linear term takes over: at
  // (1.01, 1, -100) the curvature figure is some four hundred times too wide,
  // and integrating over that left excess kurtosis at -0.56 where the true
  // value is +5.93.
  const double Ldrop = 45.0;
  const double reach = mhn::kernel_decay_reach(alpha, beta, m, Ldrop);
  const double s_lo = std::max(-m, -reach);
  const double s_hi = reach;
  if (!(s_hi > s_lo)) return out;

  auto w = [alpha, beta, m](double s) -> double {
    const double t = s / m;
    if (t <= -1.0) return 0.0;
    const double e = (alpha - 1.0) * (std::log1p(t) - t) - beta * s * s;
    return (e < -700.0) ? 0.0 : std::exp(e);
  };

  const double tol = std::sqrt(std::numeric_limits<double>::epsilon());
  auto integrate = [&](const std::function<double(double)>& f) -> double {
    using boost::math::quadrature::gauss_kronrod;
    // Split at the peak, s = 0, so each panel is monotone.
    return gauss_kronrod<double, 31>::integrate(f, s_lo, 0.0, 15, tol)
         + gauss_kronrod<double, 31>::integrate(f, 0.0, s_hi, 15, tol);
  };

  const double s0 = integrate([&](double s) { return w(s); });
  if (!(s0 > 0.0) || !std::isfinite(s0)) return out;

  // Mean as an offset from the peak, which is the accurate way to hold it.
  const double d_mu = integrate([&](double s) { return s * w(s); }) / s0;
  out.mu = m + d_mu;
  out.c2 = integrate([&](double s) { const double d = s - d_mu; return d * d * w(s); }) / s0;
  out.c3 = integrate([&](double s) { const double d = s - d_mu; return d * d * d * w(s); }) / s0;
  out.c4 = integrate([&](double s) { const double d = s - d_mu; const double d2 = d * d;
                                     return d2 * d2 * w(s); }) / s0;
  out.ok = std::isfinite(out.mu) && out.mu > 0.0 &&
           std::isfinite(out.c2) && out.c2 > 0.0 &&
           std::isfinite(out.c3) && std::isfinite(out.c4);
  return out;
}

}  // namespace

// [[Rcpp::export(.mhn_mean_cpp, rng = false)]]
double mhn_mean_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);

  // E[X] is a ratio of two Fox-Wright values, and for a strong tilt both of
  // their logarithms are about z^2/4 while their difference is only about
  // log z.  Subtracting them therefore costs an absolute error of roughly
  // (z^2/4) times the unit roundoff in the exponent: harmless at moderate
  // tilt, but 31% in the answer by gamma = 1e8.  Past the point where that
  // error would exceed the working tolerance, take the mean from the same
  // quadrature the central moments use, which never forms the difference.
  // The cancellation is one-sided.  For gamma > 0 both logarithms grow like
  // z^2/4 and the difference stays about log z, so subtracting them costs
  // (z^2/4) times the unit roundoff.  For gamma < 0 there is no such growth --
  // log Psi is only O(log|z|) there -- and the closed form is exact, so
  // switching to quadrature on that side traded a correct answer for a
  // truncation error of 1.4e-3.
  const double z = gamma / std::sqrt(beta);
  if (z > 0.0) {
    const double subtraction_error =
        0.5 * z * z * std::numeric_limits<double>::epsilon();
    if (subtraction_error > 1e-10) {
      const CentralMoments m = central_moments(alpha, beta, gamma);
      if (m.ok && m.mu > 0.0) return m.mu;
    }
  }

  const double log_psi_num = mhn::mhn_log_normalizing_const(alpha + 1.0,
                                                            beta, gamma, -1.0);
  const double log_psi_den = mhn::mhn_log_normalizing_const(alpha, beta,
                                                            gamma, -1.0);
  return std::exp(log_psi_num - 0.5 * std::log(beta) - log_psi_den);
}

// [[Rcpp::export(.mhn_var_cpp, rng = false)]]
double mhn_var_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  // The Lemma 2c form alpha/(2 beta) + mu (gamma/(2 beta) - mu) is exact
  // algebra and, where its terms do not cancel, more accurate than any
  // quadrature.  Under a strong negative tilt they do cancel: the first two are
  // each about alpha/2 while the answer is about alpha/gamma^2, so the form
  // amplifies the rounding of its own terms by gamma^2/(2 beta alpha).  The
  // mean is accurate to 1e-15 relative throughout, and yet the variance came
  // out 17 percent high at gamma = -1e7 and clamped to exactly 0 by -1e8.
  //
  // So measure the cancellation and switch only when it bites.  With the
  // largest term T, the closed form carries an absolute error of about eps T,
  // hence a relative error of eps T / v; asking that to stay under the working
  // tolerance is asking for v / T above 1e-8.  Below that, take the second
  // central moment from the same quadrature the third and fourth already use,
  // which forms no such difference.
  const double mu = mhn_mean_cpp(alpha, beta, gamma);
  const double t1 = alpha / (2.0 * beta);
  const double t2 = mu * gamma / (2.0 * beta);
  const double t3 = mu * mu;
  const double v = t1 + t2 - t3;
  const double largest = std::max(t1, std::max(std::fabs(t2), t3));
  if (v > 0.0 && largest > 0.0 && v / largest > 1e-8) return v;

  const CentralMoments m = central_moments(alpha, beta, gamma);
  if (m.ok) return m.c2;
  return std::max(0.0, v);
}

// [[Rcpp::export(.mhn_skewness_cpp, rng = false)]]
double mhn_skewness_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const CentralMoments m = central_moments(alpha, beta, gamma);
  if (m.ok) return m.c3 / std::pow(m.c2, 1.5);

  // Fall back on the raw-moment expansion where the quadrature could not be
  // set up.  It is the exact algebra, and accurate whenever sigma is not small
  // against the mean.
  const auto r = mhn_raw_moments(alpha, beta, gamma, 3);
  const double mu = r[1];
  const double sigma2 = std::max(0.0, r[2] - mu * mu);
  const double central3 = r[3] - 3.0 * mu * r[2] + 2.0 * mu * mu * mu;
  return central3 / std::pow(sigma2, 1.5);
}

// [[Rcpp::export(.mhn_kurtosis_cpp, rng = false)]]
double mhn_kurtosis_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  const CentralMoments m = central_moments(alpha, beta, gamma);
  if (m.ok) return m.c4 / (m.c2 * m.c2) - 3.0;

  const auto r = mhn_raw_moments(alpha, beta, gamma, 4);
  const double mu = r[1];
  const double sigma2 = std::max(0.0, r[2] - mu * mu);
  const double central4 = r[4] - 4.0 * mu * r[3]
                          + 6.0 * mu * mu * r[2] - 3.0 * mu * mu * mu * mu;
  return central4 / (sigma2 * sigma2) - 3.0;
}
