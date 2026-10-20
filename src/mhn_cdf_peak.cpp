// The density and the tails of MHN(alpha, beta, gamma) measured from the peak
// of the density of log X; see mhn_cdf_peak.h for the formulas and for why
// they are needed at a large shape or a strong positive tilt.

#include "mhn_cdf_peak.h"
#include "mhn_constants.h"
#include "mhn_fma.h"
#include "mhn_psi.h"
#include "mhn_quadrature.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// beta d^2, taken as (beta d) d where d^2 alone would overflow.
inline double beta_d2(double beta, double d) {
  const double d2 = d * d;
  return (mhn::opaque_value(d2) <= kDoubleMax) ? beta * d2 : (beta * d) * d;
}

// log of int_0^inf exp(-a t - c t^2 / 2) dt from log a and log c (c > 0),
// for a slope and a curvature too large to be formed as doubles: with
// y = a / sqrt(2 c) it is log(sqrt(pi / (2 c)) erfcx(y)), and 1 / a to within
// 1 / (2 y^2) relative once y passes 1e8.
double log_quadratic_tail_from_logs(double log_a, double log_c) {
  const double log_r = 0.5 * (log_c - M_LN2);        // r = sqrt(c / 2)
  const double log_y = log_a - M_LN2 - log_r;         // y = (a / 2) / r
  if (log_y > std::log(1e8)) return -log_a;
  const double y = std::exp(log_y);
  const double log_erfcx = (y >= 2.0)
      ? std::log(mhn::erfcx_cf(y))
      : y * y + M_LN2 + R::pnorm(-M_SQRT2 * y, 0.0, 1.0, 1, 1);
  return 0.5 * std::log(M_PI) - M_LN2 - log_r + log_erfcx;
}

}  // namespace

bool use_peak_frame(double alpha, double beta, double gamma) {
  const double z = gamma / std::sqrt(beta);
  if (alpha > kPeakFrameMinAlpha || z > kPeakFrameMinZ) return true;
  // The terms of the ordinary assembly are of size up to about
  // m (1 + |log v*|), m the larger of alpha and z: alpha log x,
  // (alpha / 2) log beta and log Psi for the shape, and for the tilt the
  // quadrature abscissae in log x near the peak, whose rounding is
  // eps |log v*| of v* or z |log v*| 2^-53 of a standard deviation.  At
  // beta = 1 that is below 1e7 up to m = 1e6 for a shape and 7.3e5 for a
  // tilt, but at beta = 1e-300 or 1e300, where |log v*| is about 350, it
  // passes 1e7 from m of about 3e4.  Past 1e7 the frame is taken, which
  // keeps the rounding the ordinary assembly leaves in the result to about
  // eps 1e7 times a small factor, 2e-9.  The test is made only for m above
  // 1e3, since |log v*| is below 1e4 for every v* and it cannot hold below,
  // so an ordinary triple does not pay for log v*.
  const double m = std::max(alpha, z);
  if (!(m > kPeakFrameScaleMin)) return false;
  const double log_v = mhn::log_positive_root(beta, gamma, alpha);
  return m * (1.0 + std::fabs(log_v)) > kPeakFrameMaxTerm;
}

bool PeakFrame::set(double a, double b, double g, double tol) {
  alpha = a;
  beta = b;
  gamma = g;
  v = mhn::positive_root(b, g, a);
  if (!(v > 0.0)) return false;
  if (mhn::opaque_value(v) <= kDoubleMax) {
    // One Newton step on Q(x) = 2 beta x^2 - gamma x - alpha from v,
    // s = -Q(v) / Q'(v), and its second-order term.  Q(v) is a difference
    // of terms of size alpha or
    // gamma v that cancel to a few spacings of the doubles at v times Q'(v),
    // and it is formed divided by 2 v, as
    //   Q(v) / (2 v) = (beta v - gamma / 2) - alpha / (2 v),
    // without rounding its terms: beta v = p + e1 exactly (fma gives the
    // error of the product), p - gamma / 2 = s1 + t1 exactly (two_sum), and
    // alpha / v = q1 + q2 + r2 / v with r1 = alpha - q1 v and
    // r2 = r1 - q2 v exact (fma, q2 = r1 / v rounded).  The six terms are
    // summed with the rounding error of each addition carried along
    // (Neumaier's summation) and kept as an unrounded pair.  s1 and q1 / 2
    // agree to within a few eps of the largest of the terms p, gamma / 2 and
    // alpha / v (their difference is Q(v) / (2 v) less the small terms), and
    // are added first: that difference is exact where they lie within a
    // factor of two of each other (Sterbenz's lemma), and otherwise both are
    // of that small size already.  Every later partial sum is then of that
    // size too, and each addition error of order eps^2 times the largest
    // term, where adding t1 or e1 to s1 first would leave one of order eps
    // times it.  Dividing by 2 v rather than multiplying through keeps every
    // term within the range of the doubles at every v.
    // Q'(v) / 2 = 2 beta v - gamma / 2 is kept as a
    // pair as well, and the step is formed relative to v, as
    //   u = s / v = -(Q(v) / (2 v)) / (Q'(v) / 2),
    // of the order of eps: Q'(v) / (2 v) itself would overflow under a strong
    // negative tilt, where v is about alpha / |gamma| and gamma / (2 v) can
    // pass the largest double.  u is corrected by the residual of the
    // division (fma), and s + s_lo = u v as an unrounded pair.  Q being
    // quadratic, Q(v + s) = Q(v) + Q'(v) s + 2 beta s^2 exactly, so the root
    // lies a further -2 beta s^2 / Q'(v) from v + s, a relative
    // 2 beta s / Q'(v) < |s| / v of s (Q'(v) > 2 beta v at the root), that
    // is up to a few eps of it; that term is added to s_lo.  What is left is
    // the rounding of the sum of the addition errors, at most 6 eps times the
    // sum of their sizes (0 where every term is exact, as at
    // v = gamma / (2 beta) exactly), and of r2 / v, and of order eps^2 of
    // the step (s_err); each is of order eps^3 v.
    const double eps = std::numeric_limits<double>::epsilon();
    const double p = b * v;
    const double e1 = mhn::fma(b, v, -p);
    double t1;
    const double s1 = mhn::two_sum(p, -0.5 * g, t1);
    const double q1 = a / v;
    const double r1 = mhn::fma(-q1, v, a);
    const double q2 = r1 / v;
    const double r2 = mhn::fma(-q2, v, r1);
    const double r2_term = 0.5 * (r2 / v);
    const double terms[6] = {s1, -0.5 * q1, t1, e1, -0.5 * q2, -r2_term};
    double sum = 0.0;
    double comp = 0.0;
    double comp_size = 0.0;
    for (double term : terms) {
      double er;
      sum = mhn::two_sum(sum, term, er);
      comp += er;
      comp_size += std::fabs(er);
    }
    double qv_lo;                                  // Q(v) / (2 v)
    const double qv_hi = mhn::two_sum(sum, comp, qv_lo);
    // Q'(v) / 2 = 2 beta v - gamma / 2 = (p + e1) + (s1 + t1 + e1), the
    // pair dq_hi + dq_lo, positive at the positive root.
    double dq_lo;
    const double dq_hi = mhn::two_sum(p, s1, dq_lo);
    dq_lo += t1 + 2.0 * e1;
    const double u = -qv_hi / dq_hi;
    const double u_lo =
        (mhn::fma(-u, dq_hi, -qv_hi) - qv_lo - u * dq_lo) / dq_hi;
    s = u * v;
    s_lo = mhn::fma(u, v, -s) + u_lo * v;
    // -2 beta s^2 / Q'(v) = -(beta s / (Q'(v) / 2)) s, the first factor
    // below |s| / v.
    s_lo -= ((b * s) / dq_hi) * s;
    // Not a step of a few spacings of the doubles at v: the frame cannot be
    // placed.
    if (!std::isfinite(s) || !std::isfinite(s_lo) || !(dq_hi > 0.0)) {
      return false;
    }
    s_err = ((6.0 * eps * comp_size + eps * std::fabs(r2_term)) / dq_hi) * v
            + 8.0 * eps * eps * std::fabs(s);
    if (!(s_err < std::numeric_limits<double>::infinity())) {
      s_err = std::numeric_limits<double>::infinity();
    }
    sd = v / std::hypot(std::sqrt(a), M_SQRT2 * std::sqrt(b) * v);
    log_v = std::log(v) + std::log1p(s / v);
  } else {
    // v* lies beyond the largest double, which under a positive tilt with
    // gamma^2 far above 8 alpha beta puts it at gamma / (2 beta) + alpha /
    // gamma to double precision; offset() then forms x - v* from those two.
    // Otherwise (a shape so large against a beta so small that
    // sqrt(alpha / (2 beta)) overflows) the frame is not set up.
    if (!(g > 0.0) ||
        !(M_SQRT2 * 2.0 * std::sqrt(b) * std::sqrt(a) <= 1e-8 * g)) {
      return false;
    }
    s = 0.0;
    s_lo = 0.0;
    s_err = 0.0;
    sd = M_SQRT1_2 / std::sqrt(b);
    log_v = mhn::log_positive_root(b, g, a);
    if (!std::isfinite(log_v)) return false;
  }
  const double z = g / std::sqrt(b);
  const double log_u = 0.5 * std::log(b) + log_v;
  if (z == kInf || log_u > 300.0) {
    // The integral of phi is then sqrt(2 pi / c) with c = alpha + 2 u*^2, to
    // within relative corrections of order 1 / u*^2 < 1e-260; u* itself may
    // lie beyond the largest double, so it is taken from its logarithm.
    const double r = std::exp(-2.0 * log_u) * (0.5 * a);   // alpha / (2 u*^2)
    log_V = 0.5 * std::log(M_PI) - log_u - 0.5 * std::log1p(r);
  } else {
    log_V = mhn::psi_log_axis_relative(a, z, tol);
  }
  return std::isfinite(log_V);
}

double PeakFrame::offset(double x) const {
  if (mhn::opaque_value(v) <= kDoubleMax) return ((x - v) - s) - s_lo;
  // x - gamma / (2 beta) = (beta x - gamma / 2) / beta, whose numerator is
  // one rounding by fma, less v* - gamma / (2 beta) = alpha / (2 beta v*),
  // which is alpha / gamma here.
  return mhn::fma(beta, x, -0.5 * gamma) / beta - alpha / gamma;
}

void PeakFrame::terms(double x, double d, double* rho, double* lmu) const {
  // rho = d / v*, and log(1 + rho) - rho: from log1p_minus_u, which keeps its
  // relative accuracy for a small rho, near the peak, and from log(x / v*)
  // far from it, where 1 + rho itself would carry the rounding of rho.
  if (mhn::opaque_value(v) <= kDoubleMax) {
    *rho = d / v;
    if (!std::isfinite(*rho)) {
      // x / v* beyond the largest double: log(1 + rho) - rho is -Inf.
      *lmu = -kInf;
      return;
    }
    if (std::fabs(*rho) < 0.5) {
      *lmu = mhn::log1p_minus_u(*rho);
    } else {
      // log(x / v*) = log(x / v) - log1p(s / v), s / v of order eps.  Where
      // x / v leaves the normal range (a subnormal x, or v* above 1e24 at
      // x = 1e-300) its logarithm is taken as log x - log v*, which then
      // loses nothing: the difference is at least 708 in size, and neither
      // term is above 745.
      const double ratio = x / v;
      *lmu = ((ratio >= DBL_MIN) ? (std::log(ratio) - s / v)
                                 : (std::log(x) - log_v)) - *rho;
    }
    return;
  }
  // Beyond the largest double every finite x lies below v*, and d / v* is
  // formed on the log scale, or as x / v* - 1 where d itself overflows.
  const double log_ratio = std::log(x) - log_v;
  *rho = std::isfinite(d) ? -std::exp(std::log(-d) - log_v)
                          : std::expm1(log_ratio);
  *lmu = (std::fabs(*rho) < 0.5) ? mhn::log1p_minus_u(*rho)
                                 : log_ratio - *rho;
}

bool PeakFrame::unresolved(double d) const {
  // The offset d carries the error of s, at most s_err.  Within 45 standard
  // deviations of the peak that moves the tails by about s_err / sd and the
  // log density by |d| s_err / sd^2, and the value is resolved where s_err
  // is below 1e-11 sd.  Beyond, the tails are below exp(-1000) and their
  // logarithms and the log density are about -d^2 / (2 sd^2), which it moves
  // by a relative 2 s_err / |d|; the value is resolved where s_err is below
  // 1e-11 |d|, which also keeps the true offset, within s_err of d, beyond
  // 45 standard deviations.  (A d that lies beyond 45 sd but within s_err of
  // the peak would otherwise be read as far out in a tail.)
  //
  // The addition errors of set() are of order eps^2 times the largest term
  // of Q(v) / (2 v), so s_err is of order eps^3 v, of order eps^2 times the
  // spacing of the doubles at v*.  The first test fails only for a law
  // narrower than that spacing by a factor of order 1e19: at a tilt beyond
  // about 1e35 (the spacing is about eps z / sqrt2 standard deviations), or
  // a shape beyond about 1e70 (eps sqrt(2 alpha)).  The second can fail only
  // there too, since it needs s_err above 45e-11 sd, and then only for a
  // double within about 1e11 s_err, some 1e-20 of a spacing, of the peak,
  // as where the peak falls on a double exactly.
  const double ad = std::fabs(d);
  return s_err > 1e-11 * ((ad < 45.0 * sd) ? sd : ad);
}

double PeakFrame::log_density(double x) const {
  const double d = offset(x);
  if (unresolved(d)) return NA_REAL;
  double rho, lmu;
  terms(x, d, &rho, &lmu);
  // x / v* beyond the largest double, where beta d^2, which dominates the
  // other terms, is too.
  if (!std::isfinite(rho)) return -kInf;
  // log g(x) - log g(v*) = (alpha - 1) (log(1 + rho) - rho) - rho - beta d^2,
  // and f(x) = g(x) / (v* g(v*) V).
  return (alpha - 1.0) * lmu - rho - beta_d2(beta, d) - log_v - log_V;
}

double PeakFrame::log_tail(double x, bool lower) const {
  const double d = offset(x);
  if (unresolved(d)) return NA_REAL;
  double rho, lmu;
  terms(x, d, &rho, &lmu);
  const double log_phi_x = alpha * lmu - beta_d2(beta, d);
  if (std::isnan(log_phi_x)) return NA_REAL;
  if (log_phi_x == -kInf) return -kInf;

  // 1 / (2 v*), 0 when v* lies beyond the largest double.
  const double half_inv_v =
      (mhn::opaque_value(v) <= kDoubleMax) ? 0.5 / v : 0.0;
  // Q(t) = 2 (t - v*) (beta t + alpha / (2 v*)), the negated slope of
  // log phi at t (by the stationarity condition beta v* - gamma / 2 =
  // alpha / (2 v*)); it is negative below v* and positive above.
  auto slope_at = [this, half_inv_v](double t_minus_v, double t) -> double {
    return 2.0 * t_minus_v * (beta * t + alpha * half_inv_v);
  };
  const double qx = slope_at(d, x);
  const double dir = lower ? -1.0 : 1.0;
  // The curvature of -log phi in sigma at x, x (4 beta x - gamma); positive
  // at and above v*, and possibly negative far below a peak set by a
  // positive tilt, where the slope then carries the fall.
  const double half_curv = mhn::fma(beta, x, mhn::fma(beta, x, -0.5 * gamma));
  const double curv = 2.0 * x * half_curv;
  const double a = std::fabs(qx);

  const double L = -std::log(mhn_eps()) + 40.0;
  double h0;
  if (curv > 0.0 && std::isfinite(curv)) {
    const double ha = 0.5 * a;
    h0 = L / (ha + std::hypot(ha, std::sqrt(0.5 * L * curv)));
  } else {
    h0 = L / a;
  }
  if (!(h0 > 0.0) || !std::isfinite(a) || !std::isfinite(curv) ||
      h0 < 1e-280) {
    // A slope or curvature beyond the double range, where the tail falls by
    // L within a stretch so short that its cubic term is nothing: the tail
    // is that of exp(-a sigma - curv sigma^2 / 2), from their logarithms.
    // (At the peak itself the slope is 0, log a = -Inf, and the tail is the
    // half Gaussian sqrt(pi / (2 curv)).)
    if (!(a > 0.0) && !(curv > 0.0)) return NA_REAL;
    const double log_a = std::isfinite(a)
        ? std::log(a)
        : M_LN2 + std::log(std::fabs(d))
              + std::log(beta * x + alpha * half_inv_v);
    double log_I;
    if (curv > 0.0) {
      const double log_c = std::isfinite(curv)
          ? std::log(curv)
          : M_LN2 + std::log(x) + std::log(half_curv);
      log_I = log_quadratic_tail_from_logs(log_a, log_c);
    } else {
      log_I = -log_a;
    }
    return log_phi_x + log_I - log_V;
  }

  // The kernel exp(log phi(w_x + dir tau) - log phi(w_x)) in tau >= 0.
  auto kernel = [this, x, qx, dir](double tau) -> double {
    const double sigma = dir * tau;
    const double e = std::expm1(sigma);
    const double delta = x * e;               // t - x
    const double k = -alpha * mhn::expm1_minus_id(sigma) - qx * e
                     - beta_d2(beta, delta);
    return (k > -745.0) ? std::exp(k) : 0.0;
  };

  // Below t = x0, the positive root of beta t^2 + |gamma| t = 4.5e-9
  // (mhn::kernel_head_x0), exp(-beta t^2 + gamma t) agrees with
  // 1 - beta t^2 + gamma t to within 1e-17, and the lower tail from t_e = x
  // e^sigma_e down to the origin is, with u = t / t_e,
  //
  //   exp(K(sigma_e)) int_0^1 u^(alpha-1) (1 + gamma t_e (u - 1)
  //                                        - beta t_e^2 (u^2 - 1)) du
  //     = exp(K(sigma_e)) / alpha (1 - gamma t_e / (alpha + 1)
  //                                + 2 beta t_e^2 / (alpha + 2)),
  //
  // which serves where the slope of log phi has fallen to about alpha and a
  // small shape would otherwise need panels out to sigma of order L / alpha.
  const double x0 = lower ? mhn::kernel_head_x0(beta, gamma) : 0.0;
  auto log_head = [this](double k_e, double t_e) -> double {
    const double shape = 1.0 - gamma * t_e / (alpha + 1.0)
                         + 2.0 * (beta * t_e) * t_e / (alpha + 2.0);
    return k_e - std::log(alpha) + std::log(shape);
  };

  // Panels [0, h0], [h0, 3 h0], [3 h0, 7 h0], ... outwards from x, each
  // integrated by the 31-point Gauss-Kronrod rule (mhn_quadrature.h), until
  // what lies beyond is below exp(-L) of the sum.  Past the end tau_e of a
  // panel the slope of log phi keeps at least the size it has there: above
  // v* |Q(t)| grows with t, and below it |Q(t)| is concave in t and equal to
  // alpha at t = 0, so it is at least the smaller of alpha and its value at
  // the panel's end.  The rest is then at most exp(K(tau_e)) / m_e.
  const double panel_tol = mhn_eps() / 10.0;
  double val = 0.0;
  double lo = 0.0;
  double w = std::min(h0, 1.0);
  for (int k = 0; k < 2200; ++k) {
    if (lower) {
      const double e_lo = std::expm1(-lo);
      const double t_lo = x + x * e_lo;
      if (t_lo <= x0) {
        const double k_lo = -alpha * mhn::expm1_minus_id(-lo) - qx * e_lo
                            - beta_d2(beta, x * e_lo);
        const double lh = log_head(k_lo, t_lo);
        if (!(val > 0.0)) return log_phi_x + lh - log_V;
        const double lv = std::log(val);
        const double hi_part = std::max(lv, lh);
        return log_phi_x + hi_part
               + std::log1p(std::exp(-std::fabs(lv - lh))) - log_V;
      }
    }
    const double hi = lo + w;
    val += mhn::gk_unit<31>(kernel, lo, hi, panel_tol);
    const double sigma = dir * hi;
    const double e = std::expm1(sigma);
    const double delta = x * e;
    const double k_end = -alpha * mhn::expm1_minus_id(sigma) - qx * e
                         - beta_d2(beta, delta);
    double m_end = std::fabs(slope_at(d + delta, x + delta));
    if (lower) m_end = std::min(alpha, m_end);
    if (val > 0.0 && (k_end == -kInf ||
                      k_end - std::log(m_end) < std::log(val) - L)) {
      return log_phi_x + std::log(val) - log_V;
    }
    lo = hi;
    w *= 2.0;
    if (!std::isfinite(lo)) break;
  }
  return NA_REAL;
}

double PeakFrame::cdf(double x) const {
  if (offset(x) < 0.0) {
    const double F = std::exp(log_tail(x, true));
    return (F > 1.0) ? 1.0 : F;
  }
  double log_Q = log_tail(x, false);
  if (log_Q > 0.0) log_Q = 0.0;
  return -std::expm1(log_Q);
}

double PeakFrame::log_upper(double x) const {
  if (offset(x) < 0.0) {
    double log_F = log_tail(x, true);
    if (log_F > 0.0) log_F = 0.0;
    return mhn::log1m_exp(log_F);
  }
  const double log_Q = log_tail(x, false);
  return (log_Q > 0.0) ? 0.0 : log_Q;
}

double PeakFrame::log_cdf(double x) const {
  if (offset(x) < 0.0) {
    const double log_F = log_tail(x, true);
    return (log_F > 0.0) ? 0.0 : log_F;
  }
  double log_Q = log_tail(x, false);
  if (log_Q > 0.0) log_Q = 0.0;
  return mhn::log1m_exp(log_Q);
}

}  // namespace mhn
