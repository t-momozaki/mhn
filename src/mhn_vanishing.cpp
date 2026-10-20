// The law at a subnormal shape (see mhn_vanishing.h).

#include "mhn_vanishing.h"
#include "mhn_constants.h"
#include "mhn_fma.h"
#include "mhn_quadrature.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <boost/math/tools/toms748_solve.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mhn {

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// The forms are used while w is at most 2^64 (see mhn_vanishing.h).
const double kMaxLogW = 64.0 * M_LN2;

// Below w = 2^-20 the mass by the origin outweighs the peak's a million
// times, and the central moments are built from the raw ones: a A_1^2 / A_2
// and the like are then about w, so nothing cancels (see vanishing_moments).
const double kRawMomentMaxLogW = -20.0 * M_LN2;

// Below z = 30 w is under 1e-200 (K0(30) is about 6e96 and alpha below
// 2.3e-308), and log(1 + w) is 0 to double precision.
const double kNegligibleWMaxZ = 30.0;

// The width of the first panel on one side of s0 (dir = +1 or -1): 1/8, or
// less where g changes by more than 4 across that, so that the panel resolves
// the integrand where it starts.  Next to u under a tilt such as
// z = -1e300 the integrand of J falls by e^-4 within 1e-200 of log u.
template <class G>
double first_width(const G& g, double s0, double dir) {
  const double g0 = g(s0);
  double width = 0.125;
  for (int i = 0; i < 400 && width > 1e-300; ++i) {
    const double d = g(s0 + dir * width) - g0;
    if (std::fabs(d) <= 4.0) break;
    width *= 0.125;
  }
  return width;
}

// log of int_a^b exp(g(s)) ds, for a smooth g on the log axis s = log t.  a
// may be -Inf and b +Inf.  g_max is the largest value of g on [a, b], by
// which the integrand is scaled, and s_start the point from which the panels
// are laid out: the first on either side is as wide as first_width allows,
// at most 1/8, and each is 1.5 times as wide as the last, so that the panels
// are fine where the integrand changes fastest and few elsewhere.  A side
// ends at its limit, or where g has fallen 60 below g_max.  On the right
// every integrand here falls faster than exponentially, and on the left
// those of A_k and I fall at least like e^s, so what is left out there is
// below e^-60 of the integral.  That of J is flat on the left, about 1 from
// 1 / (1 + |z|) down to u; it is cut only below a peak at least e^60 higher,
// and what is left out, at most |log u| e^-60 against a peak of width
// about 1 / z, is below 1e-22 of the integral.  Each panel is integrated by
// the adaptive 31-point Gauss-Kronrod rule to the package's working
// tolerance.
template <class G>
double log_integral(const G& g, double a, double b, double s_start,
                    double g_max) {
  if (!(g_max > -kInf)) return g_max;
  const double s0 = std::min(std::max(s_start, a), b);
  auto h = [&g, g_max](double s) -> double {
    return std::exp(g(s) - g_max);
  };
  const double tol = mhn_eps();
  double sum = 0.0;
  double lo = s0;
  double width = first_width(g, s0, 1.0);
  for (int i = 0; i < 400 && lo < b; ++i) {
    const double hi = std::min(lo + width, b);
    sum += gk_unit<31>(h, lo, hi, tol);
    if (hi >= b || g(hi) - g_max < -60.0) break;
    lo = hi;
    width *= 1.5;
  }
  double hi = s0;
  width = first_width(g, s0, -1.0);
  for (int i = 0; i < 400 && hi > a; ++i) {
    const double lo2 = std::max(hi - width, a);
    sum += gk_unit<31>(h, lo2, hi, tol);
    if (lo2 <= a || g(lo2) - g_max < -60.0) break;
    hi = lo2;
    width *= 1.5;
  }
  return g_max + std::log(sum);
}

// -t^2 + z t at t = exp(s), as t (z - t), which is 0 rather than NaN as t
// underflows and -Inf rather than NaN as it overflows.
inline double kernel_exponent(double s, double z) {
  const double t = std::exp(s);
  return t * (z - t);
}

// log J(u), J(u) = int_u^Inf exp(-t^2 + z t) / t dt, from log u.  On the log
// axis the integrand is exp(-t^2 + z t) itself.  It is largest at t = z / 2
// when that lies above u, and at u otherwise.  For z <= 0 it is flat, about
// 1, up to t of about 1 / (1 + |z|), and the panels start there, where it
// begins to fall, or at u if u is above that point.
//
// Where the integrand falls from u on, it falls at the rate
// lambda = -g'(log u) = u (2u - z) on the log axis.  Past lambda = 1e11 the
// distance over which it falls, 1 / lambda, is too short for the log axis to
// resolve near log u, whose spacing is about 1e-13, and the integral is taken
// from its expansion at the end point instead (Watson's lemma):
// int_0^Inf exp(g(log u + r) - g(log u)) dr = 1/lambda + g''/lambda^3 + ...,
// with g'' = u (z - 4u), the next term below 1e-22 of the sum there.
double log_J(double log_u, double z) {
  auto g = [z](double s) { return kernel_exponent(s, z); };
  if (z > 0.0) {
    const double s_peak = std::log(0.5 * z);
    if (log_u < s_peak) return log_integral(g, log_u, kInf, s_peak, g(s_peak));
  }
  const double u = std::exp(log_u);
  const double lambda = u * (2.0 * u - z);
  if (lambda > 1e11) {
    const double g2 = u * (z - 4.0 * u);
    return g(log_u) - std::log(lambda) + std::log1p(g2 / lambda / lambda);
  }
  if (z > 0.0) return log_integral(g, log_u, kInf, log_u, g(log_u));
  const double s_fall = -std::log1p(-z);
  return log_integral(g, log_u, kInf, std::max(log_u, s_fall), g(log_u));
}

// log I(u), I(u) = int_0^u (exp(-t^2 + z t) - 1) / t dt, from log u, for
// z > 0 and u below z, where the integrand is positive.  G(u) = log u + I(u)
// gives the lower tail (see VanishingShape::log_lower).  On the log axis the
// integrand is expm1(t (z - t)), formed as phi + log(1 - exp(-phi)) so that
// it neither overflows nor cancels.  It rises like z t from the origin to its
// largest point, at u or at z / 2 below u, where the panels start.
double log_I(double log_u, double z) {
  auto g = [z](double s) {
    const double phi = kernel_exponent(s, z);
    return phi + log1m_exp(-phi);
  };
  const double s_top = std::min(log_u, std::log(0.5 * z));
  return log_integral(g, -kInf, log_u, s_top, g(s_top));
}

// log A_k, A_k = int_0^Inf t^(k-1) exp(-t^2 + z t) dt.  On the log axis the
// integrand is exp(k s - t^2 + z t), largest where 2 t^2 - z t - k = 0, at
// t = (z + sqrt(z^2 + 8 k)) / 4, written for z <= 0 as
// 2 k / (sqrt(z^2 + 8 k) - z) so that it does not cancel, and with hypot so
// that z^2 cannot overflow.
double log_A(int k, double z) {
  const double kk = static_cast<double>(k);
  const double r = std::hypot(z, std::sqrt(8.0 * kk));
  const double t_peak = (z > 0.0) ? 0.25 * (z + r) : 2.0 * kk / (r - z);
  auto g = [z, kk](double s) { return kk * s + kernel_exponent(s, z); };
  const double s_peak = std::log(t_peak);
  return log_integral(g, -kInf, kInf, s_peak, g(s_peak));
}

// log D(z) for z > 30, the only tilts for which it is needed (see set_tilt),
// D(z) = int_0^Inf exp(-t^2 + z t) (1 - exp(-z t)) / t dt, so that
// K0(z) = D(z) - Euler's constant / 2: K0(0) = -Euler's constant / 2, from
// int_0^Inf t^(alpha-1) exp(-t^2) dt = Gamma(alpha/2) / 2
// = 1/alpha - 0.5772.../2 + O(alpha), and K0(z) - K0(0) is
// int_0^Inf exp(-t^2) (exp(z t) - 1) / t dt.
//
// With m = z / 2, D is taken from its expansion about the peak,
//   D = exp(m^2) (sqrt(pi) / m) (1 + sum_{k>=1} (2k-1)!! / (2 m^2)^k),
// the series of 1 / (m + d) over d ~ N(0, 1/2).  What it leaves out, the
// factor 1 - exp(-z t) and the parts below t = m / 2, is below exp(-m^2 / 4)
// of D, 3.7e-25 at z = 30, and its terms fall below 1e-17 of the sum within
// twelve terms there and seven at m = 50, far before they would grow again,
// near k = m^2.  m^2 is split exactly into hi + lo (fma), and the small
// terms are added to lo first, so that the result is rounded once at the
// size of m^2.  Its error is that rounding, half a unit in the last place,
// plus the rounding of the small terms, a few units in the last place of
// log m: below 1/64 of a unit of log D each for z > 30, and in all below
// 0.15 units, so that the whole is within 0.65 units.  Against D at 256
// bits (this expansion, which agrees with D's Taylor series
// sum_n z^n Gamma(n/2) / (2 n!) to within 2e-48 from z = 30 to 99) it is
// about half a unit off at most on scans of z from 30 to 1e150.  Added left
// to right, as m^2 + 0.5 log(pi) - log m + log1p(sum), it was up to 1.2
// units off at 308 points from z = 30 to 99 where the compiler fuses
// m^2 + 0.5 log(pi) into one rounding, and 1.4 where it does not.  A
// quadrature on the log axis, which this replaced, cost a setup of rmhn()
// seven times its own where the tilt changes from draw to draw, and from z
// of about 1e15 on, where the peak is finer than the spacing of the doubles
// near log m, it came out -Inf.
double log_D(double z) {
  const double m = 0.5 * z;
  if (!(m * m < kInf)) return kInf;
  const double r = 0.5 / (m * m);
  double term = 1.0;
  double sum = 0.0;
  for (int k = 1; k <= 60; ++k) {
    term *= (2.0 * k - 1.0) * r;
    sum += term;
    if (term <= 1e-17 * sum) break;
  }
  const double hi = m * m;
  const double lo = mhn::fma(m, m, -hi);
  return hi + (lo + (0.5 * std::log(M_PI) - std::log(m) + std::log1p(sum)));
}

// log P(X <= x) at the quantile x of p under the flags, and NaN where p is NA
// or NaN or outside the range of a probability.
double lower_log_target(double p, bool lower_tail, bool log_p) {
  if (ISNAN(p)) return R_NaN;
  if (log_p) {
    if (p > 0.0) return R_NaN;
    return lower_tail ? p : log1m_exp(p);
  }
  if (p < 0.0 || p > 1.0) return R_NaN;
  return lower_tail ? std::log(p) : std::log1p(-p);
}

}  // namespace

bool VanishingShape::set(double a, double b, double g) {
  lower_only = false;
  if (!(a > 0.0 && vanishing_candidate(a))) return false;
  set_tilt(b, g);
  return set_shape(a);
}

void VanishingShape::set_tilt(double b, double g) {
  gamma = g;
  sqrt_beta = std::sqrt(b);
  z = g / sqrt_beta;
  half_log_beta = 0.5 * std::log(b);
  // K0 = D - 0.2886..., and D is above 6e96 for z above 30, so
  // log K0 = log D to double precision.  Past z of about 2.7e154 log D passes
  // the largest double, and log w and log(1 + w) are Inf.
  log_D_z = (std::isfinite(z) && z > kNegligibleWMaxZ)
      ? log_D(z) : std::numeric_limits<double>::quiet_NaN();
  // Past the hand-over, alpha G(u) is above 2^64 from the point where
  // t (z - t) = 900 on: from z = 60 on that is the smaller root
  // 1800 / (z (1 + sqrt(1 - 3600 / z^2))), below z / 2, where I(u) is at
  // least e^899 / 1800 (t (z - 2t) is at most 1800 up to it, so the integrand
  // stays above e^899 over the last 1/1800 of log u), and alpha I(u) at
  // least e^147 for every subnormal alpha; below z = 60 it is z / 2 + 5,
  // where alpha G(u) is w times 1 - 7e-13.
  log_u_far = (z >= 60.0)
      ? std::log(1800.0 / z) -
            std::log1p(std::sqrt((1.0 - 60.0 / z) * (1.0 + 60.0 / z)))
      : std::log(0.5 * z + 5.0);
}

bool VanishingShape::set_shape(double a) {
  lower_only = false;
  if (!(a > 0.0 && vanishing_candidate(a)) || !std::isfinite(z)) return false;
  log_alpha = std::log(a);
  if (z <= kNegligibleWMaxZ) {
    log_w = -kInf;
    log1p_w = 0.0;
    return true;
  }
  // A log D of -Inf or NaN would be a failed evaluation, not a negligible w,
  // and is not taken as one.
  if (!(log_D_z > -kInf)) return false;
  log_w = log_alpha + log_D_z;
  log1p_w = log1p_exp(log_w);
  if (log_w <= kMaxLogW) return true;
  lower_only = true;  // past the hand-over
  return false;
}

// log u from x > 0: the log of the product where that is a normal double,
// which keeps log u to the rounding of u, and otherwise the sum of the logs,
// whose rounding grows with |log x| (5.7e-14 at x = 1e150).
double VanishingShape::log_u_of(double x) const {
  const double u = sqrt_beta * x;
  if (u >= std::numeric_limits<double>::min() && u < kInf) return std::log(u);
  return half_log_beta + std::log(x);
}

// x from s = log u, by the same rule.
double VanishingShape::x_of(double s) const {
  const double u = std::exp(s);
  if (u >= std::numeric_limits<double>::min() && u < kInf) {
    const double x = u / sqrt_beta;
    if (x >= std::numeric_limits<double>::min() && x < kInf) return x;
  }
  return std::exp(s - half_log_beta);
}

double VanishingShape::log_density(double x) const {
  if (x == kInf) return -kInf;
  if (x > 0.0) {
    // beta x^2 - gamma x = u (u - z) with u = sqrt(beta) x, which cannot
    // form Inf - Inf; where u underflows to 0 the term is -gamma x.
    const double u = sqrt_beta * x;
    const double q = (u > 0.0) ? u * (u - z) : -gamma * x;
    return log_alpha - std::log(x) - q - log1p_w;
  }
  if (x == 0.0) return kInf;  // alpha < 1
  return -kInf;
}

double VanishingShape::log_upper(double log_u) const {
  const double l = log_alpha + log_J(log_u, z) - log1p_w;
  return (l > 0.0) ? 0.0 : l;
}

// The lower tail is (1 + alpha G(u)) / (1 + w) with G(u) = log u + I(u).
// |log u| is below 1117 for any u = sqrt(beta) x with beta and x positive
// doubles, so alpha |log u| is below 2.5e-305, and 1 + alpha G(u) is
// 1 + alpha I(u) to double precision.
double VanishingShape::log_lower(double log_u) const {
  const double l = log1p_exp(log_alpha + log_I(log_u, z)) - log1p_w;
  return (l > 0.0) ? 0.0 : l;
}

double VanishingShape::p_one(double q, bool lower_tail, bool log_p) const {
  if (q <= 0.0) {
    return lower_tail ? (log_p ? -kInf : 0.0) : (log_p ? 0.0 : 1.0);
  }
  if (q == kInf) {
    return lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? -kInf : 0.0);
  }
  const double log_u = log_u_of(q);
  bool above_half = false;
  if (log_w > 0.0 && log_u < std::log(0.5 * z)) {
    // Below one half, where w > 1 lets the lower tail fall, it is taken
    // from its own form: log Q is the sum of three terms of about 750 and
    // carries their rounding, 1e-13, which 1 - Q would pass on to a small F
    // whole.  From z / 2 on F is above one half (the peak's part below its
    // centre is the larger, as 1 / t is), and below it the lower tail is
    // tried first, so that one tail is integrated at most points.
    const double log_F = log_lower(log_u);
    if (log_F <= -M_LN2) {
      if (lower_tail) return log_p ? log_F : std::exp(log_F);
      return log_p ? log1m_exp(log_F) : -std::expm1(log_F);
    }
    above_half = true;
  }
  double log_Q = log_upper(log_u);  // log P(X > q)
  // Where the lower tail found F above one half, Q is kept at most one half:
  // the two forms differ by their rounding there (2e-14 near the median), and
  // a Q above one half just past the switch would let F step back.
  if (above_half && log_Q > -M_LN2) log_Q = -M_LN2;
  if (!lower_tail) return log_p ? log_Q : std::exp(log_Q);
  if (log_p) return log1m_exp(log_Q);
  // +0 rather than -0 where Q rounds to 1.
  const double F = -std::expm1(log_Q);
  return (F > 0.0) ? F : 0.0;
}

double VanishingShape::q_one(double p, bool lower_tail, bool log_p) const {
  if (R_IsNA(p)) return NA_REAL;
  if (ISNAN(p)) return R_NaN;
  // log Q, Q = P(X > x) at the quantile x.
  double log_Q;
  if (log_p) {
    if (p > 0.0) return R_NaN;
    log_Q = lower_tail ? log1m_exp(p) : p;
  } else {
    if (p < 0.0 || p > 1.0) return R_NaN;
    log_Q = lower_tail ? std::log1p(-p) : std::log(p);
  }
  if (log_Q == -kInf) return kInf;
  if (log_Q >= 0.0) return 0.0;
  // A lower tail below one half, where w > 1 lets it fall that low, is
  // solved on its own form, for the reason p_one takes it from there.
  if (log_w > 0.0) {
    const double log_F = lower_log_target(p, lower_tail, log_p);
    if (log_F < -M_LN2) {
      const double x = lower_root(log_F, std::log(0.5 * z + 5.0));
      if (!ISNAN(x)) return x;
    }
  }

  // Solve log J(u) = target for s = log u; log J falls as s rises.
  const double target = log_Q + log1p_w - log_alpha;
  const double zz = z;
  auto f = [zz, target](double s) { return log_J(s, zz) - target; };
  const double s_min =
      std::log(std::numeric_limits<double>::denorm_min()) + half_log_beta;
  const double s_max =
      std::log(std::numeric_limits<double>::max()) + half_log_beta;
  const double f_min = f(s_min);
  // P(X > x) at the smallest positive double is already at most Q: the
  // quantile lies below that double and is 0 in double precision.
  if (!(f_min > 0.0)) return 0.0;
  // Step up from the point where the kernel starts to fall, in steps that
  // double, until log J is below the target.
  const double s_scale =
      (zz > 0.0) ? std::log(0.5 * zz + 1.0) : -std::log1p(-zz);
  double a = s_min, fa = f_min;
  double b = std::max(s_scale, s_min + 1.0);
  double fb = f(b);
  double step = 1.0;
  while (fb > 0.0) {
    a = b;
    fa = fb;
    if (b >= s_max) return kInf;
    b = std::min(b + step, s_max);
    step *= 2.0;
    fb = f(b);
  }
  if (fb == 0.0) return x_of(b);
  if (fb == -kInf) {
    // log J is -Inf where -t^2 + z t is; close in on its last finite value.
    double lo = a, hi = b;
    for (int i = 0; i < 200 && fb == -kInf; ++i) {
      const double mid = 0.5 * (lo + hi);
      const double fm = f(mid);
      if (fm > 0.0) {
        lo = mid;
        a = mid;
        fa = fm;
      } else if (fm == -kInf) {
        hi = mid;
      } else {
        b = mid;
        fb = fm;
      }
    }
  }
  if (fb == -kInf) return x_of(a);
  std::uintmax_t max_iter = 100;
  auto close = [](double l, double r) { return r - l <= 1e-14; };
  const std::pair<double, double> root =
      boost::math::tools::toms748_solve(f, a, b, fa, fb, close, max_iter);
  return x_of(0.5 * (root.first + root.second));
}

// log(1 + alpha I(u)) = log_F + log(1 + w), solved for s = log u as
// log I(u) = log(expm1(log_F + log(1 + w))) - log alpha; I rises with u.
double VanishingShape::lower_root(double log_F, double log_u_top) const {
  if (log_F == -kInf) return 0.0;
  const double T = log_F + log1p_w;
  // P(X <= x) is at least 1 / (1 + w) at every x > 0, to double precision
  // (less a relative alpha |log u|, below 2.5e-305).
  if (!(T > 0.0)) return 0.0;
  const double target = T + log1m_exp(-T) - log_alpha;
  const double zz = z;
  auto f = [zz, target](double s) { return log_I(s, zz) - target; };
  const double s_min =
      std::log(std::numeric_limits<double>::denorm_min()) + half_log_beta;
  if (!(s_min < log_u_top)) return R_NaN;
  double b = log_u_top;
  double fb = f(b);
  if (!(fb > 0.0)) return R_NaN;
  // Step down from the top, in steps that double, until log I is below the
  // target: the root lies within a few units of log u of the top unless the
  // target is small, and a bracket as wide as the doubles' range (1090 at
  // beta = 1e-300) can use up the iterations of the solver before it closes.
  double a = b, fa = fb;
  double step = 1.0;
  while (fa >= 0.0) {
    // P(X <= x) at the smallest positive double already reaches log_F: the
    // quantile lies below that double and is 0 in double precision.
    if (a <= s_min) return 0.0;
    b = a;
    fb = fa;
    a = std::max(a - step, s_min);
    step *= 2.0;
    fa = f(a);
  }
  if (!(fa < 0.0)) return R_NaN;
  std::uintmax_t max_iter = 100;
  auto close = [](double l, double r) { return r - l <= 1e-14; };
  const std::pair<double, double> root =
      boost::math::tools::toms748_solve(f, a, b, fa, fb, close, max_iter);
  return x_of(0.5 * (root.first + root.second));
}

bool VanishingShape::log_lower_far(double q, double* log_F) const {
  if (!lower_only || !(q > 0.0 && q < kInf)) return false;
  const double log_u = log_u_of(q);
  if (!(log_u < log_u_far)) return false;
  const double l = log_alpha + log_I(log_u, z);
  if (!(l <= kMaxLogW)) return false;
  const double lF = log1p_exp(l) - log1p_w;
  // Above one half the mass by the origin is below 2^-63 of F, and the
  // ordinary paths serve it; here log F would be the difference of two terms
  // of about 44 just past the hand-over, and lose its relative accuracy.
  if (!(lF <= -M_LN2)) return false;
  *log_F = lF;
  return true;
}

bool VanishingShape::q_lower_far(double p, bool lower_tail, bool log_p,
                                 double* x) const {
  if (!lower_only) return false;
  const double log_F = lower_log_target(p, lower_tail, log_p);
  // The largest lower tail log_lower_far serves: at alpha G(u) = 2^64, or
  // one half.
  if (!(log_F <= std::min(log1p_exp(kMaxLogW) - log1p_w, -M_LN2))) {
    return false;
  }
  const double r = lower_root(log_F, log_u_far);
  if (ISNAN(r)) return false;
  *x = r;
  return true;
}

VanishingMoments vanishing_moments(const VanishingShape& v) {
  VanishingMoments out;
  const double z = v.z;
  if (v.log_w < kRawMomentMaxLogW) {
    // E U^k = a A_k, a = alpha / (1 + w).  The central moments are a times
    //   nu_2 = A_2 - a A_1^2,
    //   nu_3 = A_3 - 3 a A_1 A_2 + 2 a^2 A_1^3,
    //   nu_4 = A_4 - 4 a A_1 A_3 + 6 a^2 A_1^2 A_2 - 3 a^3 A_1^4,
    // each written as A_k times one minus terms of the size of w, and every
    // product is formed on the log scale, where nothing over- or underflows.
    double lA[5];
    for (int k = 1; k <= 4; ++k) lA[k] = log_A(k, z);
    const double la = v.log_alpha - v.log1p_w;
    const double n2 = 1.0 - std::exp(la + 2.0 * lA[1] - lA[2]);
    const double n3 = 1.0 - 3.0 * std::exp(la + lA[1] + lA[2] - lA[3])
                      + 2.0 * std::exp(2.0 * la + 3.0 * lA[1] - lA[3]);
    const double n4 = 1.0 - 4.0 * std::exp(la + lA[1] + lA[3] - lA[4])
                      + 6.0 * std::exp(2.0 * la + 2.0 * lA[1] + lA[2] - lA[4])
                      - 3.0 * std::exp(3.0 * la + 4.0 * lA[1] - lA[4]);
    const double l2 = lA[2] + std::log(n2);  // log nu_2
    out.mean = std::exp(la + lA[1] - v.half_log_beta);
    out.var = std::exp(la + l2 - 2.0 * v.half_log_beta);
    out.skewness = std::exp(lA[3] + std::log(n3) - 1.5 * l2 - 0.5 * la);
    out.kurtosis = std::exp(lA[4] + std::log(n4) - 2.0 * l2 - la) - 3.0;
    return out;
  }
  // w >= 2^-20, which needs z above about 52: the mass at the peak counts.
  // The law is the mixture of the origin, with weight q = 1 / (1 + w), and
  // the peak, with weight 1 - q, whose mean m1 and central moments v_j are
  // integrated about t = z / 2, from t = 0.5, where the kernel is under
  // exp(-600) of its peak (z / 2 is above 26 here), to z / 2 + 35, where it
  // is under exp(-1200).  (The spread of the part by the origin, of the
  // order of alpha, adds nothing at double precision.)
  const double m = 0.5 * z;
  const double lo = 0.5;
  const double hi = m + 35.0;
  const double tol = mhn_eps();
  auto integrate = [lo, hi, tol](const auto& f) {
    double s = 0.0;
    for (double a = lo; a < hi; a += 5.0) {
      s += gk_unit<31>(f, a, std::min(a + 5.0, hi), tol);
    }
    return s;
  };
  const double K = integrate([m](double t) {
    const double d = t - m;
    return std::exp(-d * d) / t;
  });
  const double m1 = integrate([m](double t) {
    const double d = t - m;
    return std::exp(-d * d);
  }) / K;
  auto central = [&integrate, m, m1, K](int j) {
    return integrate([m, m1, j](double t) {
      const double d = t - m;
      return std::pow(t - m1, j) * std::exp(-d * d) / t;
    }) / K;
  };
  const double v2 = central(2), v3 = central(3), v4 = central(4);
  const double q = std::exp(-v.log1p_w);
  const double p = -std::expm1(-v.log1p_w);
  const double c = p * m1;
  const double d = q * m1;  // m1 - c
  const double mu2 = q * c * c + p * (v2 + d * d);
  const double mu3 = -q * c * c * c + p * (v3 + 3.0 * v2 * d + d * d * d);
  const double mu4 =
      q * c * c * c * c
      + p * (v4 + 4.0 * v3 * d + 6.0 * v2 * d * d + d * d * d * d);
  out.mean = c / v.sqrt_beta;
  out.var = mu2 / v.sqrt_beta / v.sqrt_beta;
  out.skewness = mu3 / std::pow(mu2, 1.5);
  out.kurtosis = mu4 / (mu2 * mu2) - 3.0;
  return out;
}

}  // namespace mhn

// log D(z) of the subnormal-shape forms, for the test of its accuracy.
// [[Rcpp::export(.vanishing_log_D, rng = false)]]
double vanishing_log_D(double z) { return mhn::log_D(z); }
