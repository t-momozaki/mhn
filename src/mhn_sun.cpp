// Sun et al. (2023) Algorithm 1 (gamma > 0, alpha > 1) and Algorithm 3
// (gamma <= 0, alpha > 0) for sampling from the Modified Half-Normal.
//
// Algorithm 2 (gamma > 0, any alpha > 0; the only one of Sun et al.'s
// samplers that covers alpha < 1 with gamma > 0) is intentionally NOT
// implemented.  It is inapplicable for gamma < 0 (negative mixing
// weights), and Sun et al. (2023, Section 4.1.2) note that its running
// time becomes significant when alpha or gamma^2 / beta is large: the
// truncation index M-dagger = max{[alpha], [gamma^2 / (eps1^2 beta)]} of its
// mixture weights (their Lemma 8) grows with alpha and with gamma^2 / beta,
// so its cost has no bound uniform in the parameters.
// For alpha < 1, gamma > 0 the dispatcher in rmhn_cpp routes to the
// Gao & Wang (2025) RTDR sampler in mhn_rtdr.cpp instead.

#include "mhn_fma.h"
#include "mhn_sun.h"
#include "mhn_check.h"
#include "mhn_stable.h"

#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// =====================================================================
// Algorithm 1 helpers
// =====================================================================

// Stirling remainder r(z) = lgamma(z) - [(z - 1/2) log z - z + log(2 pi) / 2],
// by its asymptotic series
//
//     r(z) = 1/(12 z) - 1/(360 z^3) + 1/(1260 z^5) - 1/(1680 z^7) + ...
//
// It is called only with z = alpha/2 > 5000 (see log_K1_minus_K2), where the
// first omitted term, 1/(1188 z^9), is below 1e-36; the difference above,
// formed directly, would lose about eps z log z.
double stirling_remainder(double z) {
  const double w = 1.0 / (z * z);
  return (1.0 / 12.0 - w * (1.0 / 360.0 - w * (1.0 / 1260.0 - w / 1680.0))) / z;
}

// log K_1 - log K_2 for Algorithm 1 at the optimal constants mu_opt and
// delta_opt, as a function of alpha and Delta = gamma / sqrt(beta).
//
// K_1 and K_2 are the envelope constants of Sun et al. (2023, Theorem 1).
// Psi is a common factor of both and cancels, and Theorem 1 also states that
// at the optimal constants both depend on (beta, gamma) only through Delta.
// Written with beta and gamma, the logarithm of their ratio is
//
//     log 2 + log(pi)/2 + (alpha-1) log( sqrt(beta) (alpha-1) / (2 beta mu - gamma) )
//       - (alpha-1) + beta mu^2 - (alpha/2) log beta - lgamma(alpha/2)
//       - gamma^2 / (4 (beta - delta)) + (alpha/2) log delta,
//
// a sum of terms of size alpha log alpha and gamma^2 / beta that cancel to a
// value of order one.  Each carries a rounding of about eps times its size,
// some 8 at alpha = 1e15 and 1e4 at alpha = 1e18, so there the sign that picks
// the proposal is noise.  The cancellation can be done in closed form.
// With s1 = sqrt(8 (alpha-1)), s0 = sqrt(8 alpha), t = Delta / s1,
// t' = Delta / s0 and h(t) = sqrt(1 + t^2), the optima give
// 2 beta mu - gamma = sqrt(beta) s1 (h(t) - t) / 2 and
// delta = beta (h(t') - t')^2, and with h(t) + t = exp(asinh t) every term
// of size alpha log alpha collects into
//
//     log K_1 - log K_2 = log(2)/2 + 1/2 + ((alpha-1)/2) log1p(-1/alpha)
//                         - asinh(t) + alpha asinh(q)
//                         - Delta / (s1 h(t) + s0 h(t')) - r(alpha/2),
//
// where r is the Stirling remainder above and asinh(q) = asinh(t) - asinh(t'),
// that is
//
//     q = t h(t') - t' h(t) = (t^2 - t'^2) / (t h(t') + t' h(t)),
//
// with t^2 - t'^2 = Delta^2 / (8 alpha (alpha-1)).  Every term is then of
// order one or computed without a subtraction of nearly equal quantities:
// t - t' = 8 Delta / (s0 s1 (s0 + s1)) since s0^2 - s1^2 = 8, and q is
// evaluated divided through by t, with rho = t'/t = s1 / s0, so that neither
// t^2 nor t h(t') overflows at a large tilt.
double log_K1_minus_K2_stable(double alpha, double tilt) {
  const double s1 = std::sqrt(8.0 * (alpha - 1.0));
  const double s0 = std::sqrt(8.0 * alpha);
  const double t   = tilt / s1;
  const double t_p = tilt / s0;
  const double h   = std::hypot(1.0, t);
  const double h_p = std::hypot(1.0, t_p);
  const double rho = s1 / s0;
  const double t_minus_t_p = 8.0 * tilt / (s0 * s1 * (s0 + s1));
  const double q = t_minus_t_p * (1.0 + rho) / (h_p + rho * h);
  return 0.5 * std::log(2.0) + 0.5
       + 0.5 * (alpha - 1.0) * std::log1p(-1.0 / alpha)
       - std::asinh(t) + alpha * std::asinh(q)
       - tilt / (s1 * h + s0 * h_p)
       - stirling_remainder(0.5 * alpha);
}

// The same quantity with only the terms of size gamma^2 / beta cancelled.
// With D = gamma^2 + 8 (alpha-1) beta and E = gamma^2 + 8 alpha beta, the
// optima give mu_opt = (gamma + sqrt(D)) / (4 beta) and beta - delta_opt
// = 2 beta gamma / (sqrt(E) + gamma), so
//
//     beta mu_opt^2 - gamma^2 / (4 (beta - delta_opt))
//       = -gamma / (sqrt(D) + sqrt(E)) + (alpha-1)/2,
//
// using sqrt(D) - sqrt(E) = -8 beta / (sqrt(D) + sqrt(E)).  The terms of size
// alpha log alpha and alpha log beta are left as they are.  Against a 700-bit
// evaluation their rounding is at most 3e-11 for alpha up to 1e4 and beta in
// [1e-3, 1e4], and 3e-10 at beta = 1e-170 or 1e160; the sign they decide
// only matters when K_1 and K_2 are that close, where either proposal
// serves.  This form costs about a third of the one above, and the setup is
// built once per distinct triple, so it is used up to alpha = 1e4.
double log_K1_minus_K2_direct(double alpha, double beta, double gamma) {
  const double D = std::sqrt(gamma * gamma + 8.0 * (alpha - 1.0) * beta);
  const double E = std::sqrt(gamma * gamma + 8.0 * alpha * beta);
  const double big = -gamma / (D + E) + 0.5 * (alpha - 1.0);
  const double slope = mhn::sun_tangency_slope(alpha, beta, gamma);
  const double delta = mhn::sun_delta_opt(alpha, beta, gamma);
  return std::log(2.0) + 0.5 * std::log(M_PI)
       + (alpha - 1.0) * (0.5 * std::log(beta) + std::log(alpha - 1.0)
                          - std::log(slope))
       - (alpha - 1.0)
       - 0.5 * alpha * std::log(beta) - std::lgamma(0.5 * alpha)
       + 0.5 * alpha * std::log(delta)
       + big;
}

// The proposal selector of Algorithm 1, log K_1 - log K_2.
// Sign convention: < 0 -> use NORMAL proposal, >= 0 -> use SQRT_GAMMA.
double log_K1_minus_K2(double alpha, double beta, double gamma) {
  return (alpha <= 1e4) ? log_K1_minus_K2_direct(alpha, beta, gamma)
                        : log_K1_minus_K2_stable(alpha, gamma / std::sqrt(beta));
}

// =====================================================================
// Algorithm 3 helpers
// =====================================================================

// Right inflection point of the unnormalized MHN density f(x) for alpha > 1.
// Inflection: f''(x) = 0 <=> L'(x)^2 + L''(x) = 0 (since f > 0), where
//   L(x) = (alpha-1) log x - beta x^2 + gamma x.
// At the mode L'(x_mode) = 0, so F(x_mode) = L''(x_mode) < 0.  As x -> inf,
// F(x) -> +inf, so a single root exists in (x_mode, inf).  Bisection.
//
// Every length here is a multiple of the width 1/sqrt(-L''(x_mode)) or of
// x_mode itself, so the root found at (alpha, beta, gamma) is 1/sqrt(beta)
// times the root at (alpha, 1, gamma/sqrt(beta)), as the scale identity
// requires.  The root lies one such width past the mode exactly in the
// Gamma limit beta -> 0 (where it is (alpha - 1 + sqrt(alpha - 1)) / |gamma|)
// and asymptotically for large alpha, where f is nearly normal; a scan of
// alpha > 1.1 and tilts from -1e12 to 1e4 puts it at most 1.33 widths past
// the mode, so two widths bracket it, and the doubling loop would cover any
// shortfall.  The bisection stops at a width relative to the root; an
// absolute width would stop at the scale of the unit and return a point that
// is not a root once the root itself is below that width.
double right_inflection_point(double alpha, double beta, double gamma_signed,
                              double x_mode) {
  if (!(alpha > 1.0)) {
    Rcpp::stop("rmhn: internal invariant violated (inflection point requested "
               "for alpha <= 1). Please report this to the package maintainer.");
  }
  auto F = [alpha, beta, gamma_signed](double x) -> double {
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    const double Lp  = (alpha - 1.0) / x - 2.0 * beta * x + gamma_signed;
    const double Lpp = -(alpha - 1.0) / (x * x) - 2.0 * beta;
    return Lp * Lp + Lpp;
  };
  double lo = x_mode * (1.0 + 1e-8);
  if (F(lo) >= 0.0) {
    // F(x_mode (1 + d)) < 0 needs d < 1 / sqrt(alpha - 1 + 2 beta x_mode^2),
    // which the offset above exceeds once alpha - 1 passes about 1e16: the
    // root is then less than 1e-8 x_mode past the mode, and the bracket
    // starts a few units in the last place above it instead.
    lo = x_mode * (1.0 + 4.0 * std::numeric_limits<double>::epsilon());
  }
  const double width =
    1.0 / std::sqrt((alpha - 1.0) / (x_mode * x_mode) + 2.0 * beta);
  double hi = x_mode + 2.0 * width;
  int doublings = 0;
  while (F(hi) < 0.0 && doublings < 50) {
    hi *= 2.0;
    ++doublings;
  }
  if (F(hi) < 0.0) {
    Rcpp::stop("rmhn: could not locate the inflection point needed by the "
               "Sun Algorithm 3 setup for alpha = %g, beta = %g, gamma = %g. "
               "This is a defect; please report these parameter values to the "
               "package maintainer.", alpha, beta, gamma_signed);
  }
  for (int it = 0; it < 100; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (hi - lo < 1e-12 * mid) return mid;
    if (F(mid) < 0.0) lo = mid;
    else              hi = mid;
  }
  return 0.5 * (lo + hi);
}

// Initial matching point for Algorithm 3.  Both branches and the
// alpha <= 1.1 / alpha > 1.1 cutoff are taken from the m_init definition
// in Sun et al. (2023) Section 4.2 (the equation for m_init following
// Theorem 4; Theorem 4 itself is only about the acceptance probability
// A_neg(m)).  For alpha > 1.1 we use the inflection-point heuristic (a
// lambda-weighted mix of the mode and the right-inflection point); for
// alpha <= 1.1 we fall back to the simple closed form
// m = alpha^2 / (1 + alpha), rescaled as below.
//
// Sun et al. print this branch without beta, and we read it as the value for
// beta = 1.  The matching point carries the units of x, and sqrt(beta) X is
// MHN(alpha, 1, gamma/sqrt(beta)), so it has to be divided by sqrt(beta) to be
// the same point.  Taken literally it would make Algorithm 3's acceptance
// depend on beta at a fixed tilt, which the scale identity rules out.  The
// alpha > 1.1 branch needs no rescaling: the mode and the inflection point are
// computed in x, and the inflection bisection works in multiples of the
// density's own width, so it is scale-equivariant as well.  (build_sun_algo3
// calls this in a unit of length in which the kernel's terms are of moderate
// size; see there.)
// `used_inflex_out` reports which path was taken (for diagnostics).
double m_init_algo3(double alpha, double beta, double gamma,
                    bool* used_inflex_out) {
  if (alpha > 1.1) {
    // X_mode for the kernel x^(alpha-1) exp(gamma x - beta x^2); the same
    // closed form as mhn_mode_cpp.  Valid for alpha > 1 (else the mode sits
    // at the boundary).  Algorithm 3 runs with gamma <= 0, which is exactly
    // where the literal form cancels, so use the stable one.
    const double x_mode = mhn::positive_root(beta, gamma, alpha - 1.0);
    const double x_inflex = ::right_inflection_point(alpha, beta, gamma, x_mode);
    const double x_left = 2.0 * x_mode - x_inflex;
    if (x_left > 0.0) {
      // lambda = f~(x_left) / (f~(x_left) + f~(x_inflex)) via log-space sigmoid.
      auto log_f_tilde = [alpha, beta, gamma](double x) -> double {
        return (alpha - 1.0) * std::log(x) + gamma * x - beta * x * x;
      };
      const double diff = log_f_tilde(x_left) - log_f_tilde(x_inflex);
      const double lambda = 1.0 / (1.0 + std::exp(-diff));
      if (used_inflex_out != nullptr) *used_inflex_out = true;
      return 1.5 * lambda * x_mode + (1.0 - 1.5 * lambda) * x_inflex;
    }
    // x_left <= 0 lies outside the support, where f_MHN is 0, so Sun et al.'s
    // lambda is 0 and m_init = x_inflex, which is used (log f~ without the
    // indicator is not finite there).
    if (used_inflex_out != nullptr) *used_inflex_out = true;
    return x_inflex;
  }
  if (used_inflex_out != nullptr) *used_inflex_out = false;
  return alpha * alpha / ((1.0 + alpha) * std::sqrt(beta));
}

// The digamma and trigamma functions enter l'(m) and l''(m) below at the shape
// x = alpha (beta m + |gamma|) / (2 beta m + |gamma|) of the Gamma proposal,
// and only multiplied by f = alpha beta |gamma| / (2 beta m + |gamma|)^2,
// which is x times beta |gamma| / ((beta m + |gamma|)(2 beta m + |gamma|)).
// Near x = 0,
//     psi(x)  = -1/x - euler + (pi^2/6) x + O(x^2),
//     psi1(x) = 1/x^2 + pi^2/6 - 2 zeta(3) x + O(x^2),
// so f psi(x) and f^2 psi1(x) stay of the size of f / x however small x is,
// while psi1(x) alone passes the largest double below x = 7.46e-155, the
// reciprocal of the square root of the largest double (and psi(x) below
// 5.6e-309), where Boost's functions raise an overflow error.
// Below kTinyShape, far above those points and where the omitted terms are
// below 1e-300 of the leading ones, the two products are formed from these
// expansions with the singular terms cancelled against f.  That happens only
// for a shape alpha below 2e-150; everywhere else Boost's values are used.
constexpr double kTinyShape = 1e-150;
constexpr double kEulerGamma = 0.57721566490153286061;
constexpr double kZeta2 = 1.64493406684822643647;   // pi^2 / 6
constexpr double kZeta3 = 1.20205690315959428540;

// First derivative l'(m) of the log-acceptance for Sun Algorithm 3,
// from Sun et al. (2023, Supplementary Section 2.17).  Maximising
// l(m) = log(A_neg(m, alpha, beta, gamma)) over m gives the matching
// point that maximises the rejection-sampler acceptance probability.
// Theorem 4 is stated for alpha > 1, but its part (a), a ratio of
// integrals, and l, l', l'' hold for every alpha > 0, where the package also
// uses them.
//
// gamma_abs = |gamma| (gamma <= 0 is required for Algorithm 3).
double sun_a3_l_prime(double m, double alpha, double beta, double gamma_abs) {
  const double A = beta * m + gamma_abs;        // beta m + |gamma|
  const double B = 2.0 * beta * m + gamma_abs;  // 2 beta m + |gamma|
  const double C = beta * m * m + m * gamma_abs; // beta m^2 + m |gamma|
  const double shape_term = alpha * A / B;
  const double factor = alpha * beta * gamma_abs / (B * B);
  if (shape_term < kTinyShape) {
    // factor * psi(shape_term), with factor / shape_term = beta |gamma| / (A B).
    const double f_psi = -(beta * gamma_abs / (A * B))
                         + factor * (kZeta2 * shape_term - kEulerGamma);
    return (f_psi - factor * std::log(C)) + 2.0 * beta / B - beta / A;
  }
  return factor * (boost::math::digamma(shape_term) - std::log(C))
         + 2.0 * beta / B - beta / A;
}

// Second derivative l''(m), the derivative of l'(m) above and, like it, valid
// for every alpha > 0 (see sun_a3_l_prime); used as the Hessian in
// Newton-Raphson for m_recommend.  Sun et al. (2023,
// Supplementary Section 2.17) print its first term as
//   +4 alpha beta^2 |gamma| / (2 beta m + |gamma|)^3 [psi(.) - log(beta m^2 + m |gamma|)],
// without a minus sign.  Differentiating the factor
// alpha beta |gamma| / (2 beta m + |gamma|)^2 of l'(m) gives
// -4 alpha beta^2 |gamma| / (2 beta m + |gamma|)^3, so the sign is missing in
// the print; T1 below carries it.  The remaining terms are as printed.
double sun_a3_l_double_prime(double m, double alpha, double beta,
                             double gamma_abs) {
  const double A = beta * m + gamma_abs;
  const double B = 2.0 * beta * m + gamma_abs;
  const double C = beta * m * m + m * gamma_abs;
  const double shape_term = alpha * A / B;
  // Trailing terms of l''(m): -1/(m + |gamma|/(2 beta))^2
  //                            + 1/(m + |gamma|/beta)^2,
  // equivalent to -(2 beta / B)^2 + (beta / A)^2.
  const double T3 = -(2.0 * beta / B) * (2.0 * beta / B)
                    + (beta / A) * (beta / A);
  if (shape_term < kTinyShape) {
    // The same terms with f psi(x) and f^2 psi1(x) from their expansions at
    // x = 0 (see kTinyShape); -4 alpha beta^2 |gamma| / B^3 is -4 beta / B
    // times f.
    const double f = alpha * beta * gamma_abs / (B * B);
    const double f_over_x = beta * gamma_abs / (A * B);
    const double f_psi = -f_over_x + f * (kZeta2 * shape_term - kEulerGamma);
    const double f2_psi1 = f_over_x * f_over_x
                           + f * f * (kZeta2 - 2.0 * kZeta3 * shape_term);
    const double T1 = (-4.0 * beta / B) * (f_psi - f * std::log(C));
    const double T2 = -f2_psi1 - f * B / C;
    return T1 + T2 + T3;
  }
  const double psi  = boost::math::digamma(shape_term);
  const double psi1 = boost::math::trigamma(shape_term);
  const double diff = psi - std::log(C);
  const double T1 = (-4.0 * alpha * beta * beta * gamma_abs / (B * B * B))
                    * diff;
  const double f2 = alpha * beta * gamma_abs / (B * B);
  const double T2 = f2 * (-f2 * psi1 - B / C);
  return T1 + T2 + T3;
}

// One Newton-Raphson step on l(m): the single step from m_init that Sun et
// al. (2023, Section 4.2) recommend (their m_recomend), with the step of the
// iteration of their Supplementary Section 2.17:
//   m_recommend = m_init - l'(m_init) / l''(m_init).
// Used at every alpha > 0, as l' and l'' are (see sun_a3_l_prime).
// Falls back to m_init when l''(m_init) is non-negative (the Newton
// step would move the wrong way) or when the step produces a
// non-finite or non-positive m.
double newton_refine_m(double m_init, double alpha, double beta,
                       double gamma_abs) {
  const double lp  = sun_a3_l_prime(m_init, alpha, beta, gamma_abs);
  const double lpp = sun_a3_l_double_prime(m_init, alpha, beta, gamma_abs);
  if (!std::isfinite(lp) || !std::isfinite(lpp) || !(lpp < 0.0)) {
    return m_init;
  }
  const double m_new = m_init - lp / lpp;
  return (std::isfinite(m_new) && m_new > 0.0) ? m_new : m_init;
}

// The Algorithm 3 acceptance exponent (Algorithm 3, step 5),
//
//     R (X/m)^(1/r) - beta X^2 - |gamma| X,      R = m (beta m + |gamma|),
//
// written without cancellation.  With w = X/m, p = beta m / (beta m + |gamma|)
// and q = 1 - p = |gamma| / (beta m + |gamma|), one has 1/r = 1 + p,
// beta m^2 = p R and |gamma| m = q R, so the exponent is
//
//     R (w^(1+p) - p w^2 - q w) = R w g(v),   g(v) = (1+v)^p - 1 - p v,
//
// with v = w - 1.  By Bernoulli's inequality g <= 0, with equality only at the
// matching point v = 0.  Its three terms are each of size R and cancel to
// order one, so as written they carry a rounding of about eps R: harmless
// for moderate R, but R grows like alpha, and near alpha = 1e16 that rounding
// reaches order one and distorts the accepted law.  Here g is formed from v,
// which X - m gives with a single rounding relative to itself:
//   * for |v| < 1e-3 by its binomial series,
//       g = -(p q / 2) v^2 [1 - (1+q) v/3 + (1+q)(2+q) v^2/12
//                           - (1+q)(2+q)(3+q) v^3/60 + ...],
//     whose first omitted term is below 4e-13 of the sum there;
//   * otherwise as expm1(p log1p(v)) - p v when p <= 1/2, and as
//     (1+v) expm1(-q log1p(v)) + q v when p > 1/2, in either case two terms
//     of size min(p, q) |v| whose difference is p q v^2 / 2 to leading order,
//     so the relative rounding is about eps / |v| <= 1e-12.  The first form
//     is also taken for v <= -1/2 when p > 1/2, where it stays finite above
//     X = 0 and the second need not: p > 1/2 at a rate above 1e6 needs
//     beta m > |gamma|, and so a shape above about 1e6, at which the
//     proposal puts no mass at v <= -1/2.  For p <= 1/2 the two terms of
//     the first form differ by at least a seventh of the larger at
//     v <= -1/2, so that it loses at most three bits there.
// At X = 0 (v = -1) every term of R (w^(1+p) - p w^2 - q w) vanishes, so the
// exponent is 0 whatever p is, including a p that has underflowed to 0.
// p and q are both passed so that neither is formed as 1 minus the other.
double a3_log_acc_stable(double v, double p, double q, double rate) {
  if (v <= -1.0) return 0.0;
  double g;
  if (std::abs(v) < 1e-3) {
    g = -0.5 * p * q * v * v
        * (1.0 - v * (1.0 + q) / 3.0
                 * (1.0 - v * (2.0 + q) / 4.0 * (1.0 - v * (3.0 + q) / 5.0)));
  } else if (p <= 0.5 || v <= -0.5) {
    g = std::expm1(p * std::log1p(v)) - p * v;
  } else {
    g = (1.0 + v) * std::expm1(-q * std::log1p(v)) + q * v;
  }
  return rate * (1.0 + v) * g;
}

}  // namespace

namespace mhn {

double gamma_tail_log_draw(double t_min, double shape, double rate) {
  const double c = rate * t_min;
  if (c <= 0.5 * std::numeric_limits<double>::epsilon()) {
    return std::log(t_min) + std::log(R::runif(0.0, 1.0)) / shape;
  }
  // T = t_min S is accepted with probability exp(-c S); see mhn_sun.h.
  while (true) {
    const double log_s = std::log(R::runif(0.0, 1.0)) / shape;
    if (std::log(R::runif(0.0, 1.0)) <= -c * std::exp(log_s)) {
      return std::log(t_min) + log_s;
    }
  }
}

void stop_beyond_range(double alpha, double beta, double gamma) {
  Rcpp::stop("rmhn: cannot sample at alpha = %g, beta = %g, gamma = %g: the "
             "distribution lies at or beyond the largest double precision "
             "number, and its standard deviation is below 1e-140 times its "
             "mode, so no double precision draw could represent it.",
             alpha, beta, gamma);
}

}  // namespace mhn

namespace {

// Raised for a shape above kRgammaShapeMax, where the Gamma variate that
// Algorithm 3 proposes would carry a shift from R::rgamma's rounding larger
// than 2.5e-7 of a standard deviation (see mhn_sun.h).  The RTDR sampler
// places its draws from the mode computed to about twice double precision,
// and rmhn's default method uses it there.
[[noreturn]] void stop_algo3_shape(double alpha, double beta, double gamma) {
  Rcpp::stop("rmhn: Sun Algorithm 3 is not used at alpha = %g, beta = %g, "
             "gamma = %g: above a shape of 2^59 (about 5.8e17) the Gamma "
             "variates it proposes, drawn by R's rgamma, can be shifted by "
             "more than 1e-7 of the distribution's standard deviation. The "
             "default method and method = \"rtdr\" sample this distribution "
             "without that shift.", alpha, beta, gamma);
}

// The part of the positive root of
//
//     2 beta x^2 - gamma x - (alpha - 1) = 0,
//
// the Algorithm 1 Normal-proposal mean, that the double x0 = mu_opt leaves
// out: one Newton step on the residual at x0, formed without rounding error
// from the caller's own alpha, beta and gamma (mhn::quadratic_root_step),
// with alpha - 1 split by Fast2Sum into c and c_lo (exact for alpha > 1).
// The RTDR sampler places its draws the same way (mhn_rtdr.cpp).  0 is
// returned where a product leaves the double range.  That needs beta x0^2,
// the square of the mode of sqrt(beta) X, above the largest double, and so
// a mode of sqrt(beta) X above 1.3e154 against a standard deviation below
// 1/sqrt(2): there the part left out cannot move a draw.
double algo1_mode_lo(double alpha, double beta, double gamma, double x0) {
  const double c = alpha - 1.0;
  const double c_lo = (alpha - c) - 1.0;
  const double lo = mhn::quadratic_root_step(beta, gamma, c, c_lo, x0);
  return std::isfinite(lo) ? lo : 0.0;
}

// The Algorithm 1 setup in the caller's coordinates, with y_per_x = 1.  Every
// quantity is formed from gamma^2, alpha beta, beta and gamma / beta, which
// build_sun_algo1 keeps well inside the double range before calling this.
// With normal_only set, the Normal proposal is taken without the comparison
// with the sqrt-Gamma one, and gamma may have either sign (see
// build_sun_algo1 in mhn_sun.h).
mhn::SunAlgo1Setup build_algo1_direct(double alpha, double beta, double gamma,
                                      bool normal_only) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  mhn::SunAlgo1Setup s;
  s.alpha = alpha;
  s.beta  = beta;
  s.gamma = gamma;
  s.y_per_x = 1.0;
  s.mu_lo = 0.0;
  s.offset_normal = false;
  // Optimal constants of Sun et al. (2023) Theorem 1 (the part printed as a
  // second "a."; their supplementary proof calls it part (b)).  Both optima
  // are evaluated through the stable forms in mhn_stable.h; written literally
  // they lose precision for large |gamma|, and delta_opt can leave the
  // (0, beta) interval that the first part of Theorem 1 (printed "a.")
  // requires.
  s.mu_opt    = mhn::positive_root(beta, gamma, alpha - 1.0);
  s.delta_opt = normal_only ? nan : mhn::sun_delta_opt(alpha, beta, gamma);
  s.sigma = 1.0 / std::sqrt(2.0 * beta);

  // Proposal selection: K_2 > K_1 (i.e. log K_1 - log K_2 < 0) -> Normal,
  // otherwise -> sqrt-Gamma.  K_1 / K_2 are the per-proposal rejection
  // constants for the Normal and sqrt-Gamma envelopes respectively
  // (Sun et al. 2023, Section 4.1.1 and Theorem 1).
  //
  // Above alpha = 1e18 the Normal proposal is taken whatever the sign.  Both
  // proposals are exact, so the choice only affects speed, but a sqrt-Gamma
  // draw sqrt(T) is placed by R::rgamma, which builds T from
  // sqrt(alpha/2 - 1/2) rounded to double precision.  That moves the
  // proposal by up to half a unit in the last place of the mode, and the
  // acceptance test cannot see the shift.  At a small tilt the mode is about
  // sqrt(2 alpha) standard deviations of the law, 2^30 at alpha = 1e18, where
  // the shift is below 1.6e-7 of a standard deviation, and below 2.4e-7 with
  // the one from the rounded scale 1 / delta_opt, which makes the proposal's
  // rate differ from the delta_opt of the acceptance test (see
  // kRgammaShapeMax in mhn_sun.h); at alpha = 3e28 it is 0.03, and at still
  // larger shapes the proposal takes a single value, which the test can
  // reject every time (at alpha = 1e100 and a unit tilt, for one).  The
  // Normal draw is placed from its mean to double precision there (see
  // below).
  s.log_K1_minus_K2 = normal_only ? nan : ::log_K1_minus_K2(alpha, beta, gamma);
  s.chosen = (normal_only || s.log_K1_minus_K2 < 0.0 || alpha > 1e18)
             ? mhn::SunAlgo1Setup::NORMAL
             : mhn::SunAlgo1Setup::SQRT_GAMMA;

  // Only the chosen proposal's acceptance constants are needed.
  s.slope = 0.0;
  s.stable_normal = false;
  s.sg_curv = 0.0;
  s.sg_centre = 0.0;
  if (s.chosen == mhn::SunAlgo1Setup::NORMAL) {
    // The Normal acceptance exponent as printed carries a rounding of about
    // eps (alpha - 1), below 1e-9 up to alpha - 1 = 1e6; there it is kept,
    // with the slope computed once here, since the log it needs is cheaper
    // than the log1p of the cancellation-free form (see sample_sun_algo1).
    s.slope = mhn::sun_tangency_slope(alpha, beta, gamma);
    s.stable_normal = (alpha - 1.0 > 1e6);
    // A draw mu_opt + sigma Z is off by the rounding of mu_opt, up to half a
    // unit in its last place, 2^-53 mu_opt / sigma standard deviations of the
    // proposal: below 1.2e-7 while mu_opt / sigma is below 2^30, but 0.050 at
    // (alpha, beta, gamma) = (2, 7, 1e15 sqrt(7)), a tilt gamma / sqrt(beta)
    // of 1e15.  Past 2^30 the draw is placed from mu_opt + mu_lo instead (see
    // sample_sun_algo1).
    if (s.mu_opt > 1073741824.0 * s.sigma) {
      s.offset_normal = true;
      s.mu_lo = ::algo1_mode_lo(alpha, beta, gamma, s.mu_opt);
    }
  } else {
    // The sqrt-Gamma proposal X is accepted when
    //   log U < -(beta - delta) X^2 + gamma X - gamma^2 / (4 (beta - delta))
    // (Algorithm 1, step 16).  The right side is -(beta - delta) (X - c)^2
    // with c = gamma / (2 (beta - delta)), and expanded, its three terms are
    // each of size gamma c and cancel to order one: at (alpha, beta, gamma) =
    // (1e16, 1, 1e10) their rounding is in the thousands.  With
    // E = sqrt(gamma^2 + 8 alpha beta), delta_opt = beta (E - gamma) / (E + gamma),
    // so beta - delta_opt = 2 beta gamma / (E + gamma), and
    // c = (E + gamma) / (4 beta).  Formed as a difference, beta - delta_opt
    // would carry the rounding of delta_opt magnified by
    // beta / (beta - delta_opt), about sqrt(2 alpha) / Delta; the quotient is
    // computed as 2 gamma r with r = beta / (E + gamma), as in sun_delta_opt,
    // so that beta gamma cannot overflow.
    const double E = std::sqrt(gamma * gamma + 8.0 * alpha * beta);
    s.sg_curv   = 2.0 * gamma * (beta / (E + gamma));
    s.sg_centre = (E + gamma) / (4.0 * beta);
  }
  return s;
}

// Whether a setup from build_algo1_direct can be used as it stands: the
// chosen proposal's constants are finite and positive, the sqrt-Gamma rate is
// a normal double, with all its digits, and the Gamma variate T = X^2 of that
// proposal, of mean alpha / (2 delta_opt), stays far inside the double range.
bool algo1_setup_usable(const mhn::SunAlgo1Setup& s) {
  if (s.chosen == mhn::SunAlgo1Setup::NORMAL) {
    return std::isfinite(s.mu_opt) && s.mu_opt > 0.0 &&
           std::isfinite(s.sigma) && s.sigma > 0.0 && std::isfinite(s.slope);
  }
  return std::isfinite(s.delta_opt) &&
         s.delta_opt >= std::numeric_limits<double>::min() &&
         0.5 * s.alpha / s.delta_opt <= 1e300 &&
         std::isfinite(s.sg_curv) && std::isfinite(s.sg_centre);
}

// The Algorithm 1 setup, in the caller's coordinates, where the law is normal
// to far below double precision (see build_sun_algo1): the Normal proposal
// N(m, 1/(2 beta)) centred at the mode m, the positive root of
// 2 beta x^2 - gamma x - (alpha - 1) = 0.  The root is written as
// q + sqrt(q^2 + (alpha - 1) / (2 beta)) with q = gamma / (4 beta), and
// formed as q + hypot(q, sqrt((alpha - 1) / 2) / sqrt(beta)); each of the two
// arguments of the hypot is at most the root, so that nothing overflows
// unless the root itself does.  sigma is sqrt(1/2) / sqrt(beta), finite for
// every beta.  The Normal proposal needs no comparison with the sqrt-Gamma
// one: with its acceptance test it samples the law exactly for any
// alpha > 1, and the proposal choice only affects speed.
mhn::SunAlgo1Setup build_algo1_normal_limit(double alpha, double beta,
                                            double gamma) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  mhn::SunAlgo1Setup s;
  s.alpha = alpha;
  s.beta  = beta;
  s.gamma = gamma;
  s.y_per_x = 1.0;
  s.chosen = mhn::SunAlgo1Setup::NORMAL;
  const double c = alpha - 1.0;
  const double q = (0.25 * gamma) / beta;
  s.mu_opt = q + std::hypot(q, std::sqrt(0.5 * c) / std::sqrt(beta));
  s.mu_lo = 0.0;
  s.offset_normal = true;
  s.sigma = 0.70710678118654752 / std::sqrt(beta);
  s.slope = c / s.mu_opt;
  s.stable_normal = true;
  s.delta_opt = nan;
  s.sg_curv = 0.0;
  s.sg_centre = 0.0;
  s.log_K1_minus_K2 = nan;
  return s;
}

}  // namespace

namespace mhn {

SunAlgo1Setup build_sun_algo1(double alpha, double beta, double gamma,
                              bool normal_only) {
  if (!(alpha > 1.0)) {
    Rcpp::stop("Sun Algorithm 1 requires alpha > 1 (got %g).", alpha);
  }
  if (!normal_only && !(gamma > 0.0)) {
    Rcpp::stop("Sun Algorithm 1 requires gamma > 0 (got %g).", gamma);
  }
  if (!(beta > 0.0)) {
    Rcpp::stop("Sun Algorithm 1 requires beta > 0 (got %g).", beta);
  }
  // The setup squares gamma, forms alpha beta, and divides beta, gamma and
  // their combinations by beta; within these bounds all of them stay inside
  // the double range with all their digits, and so does the mode, at most
  // |gamma| / (2 beta) + sqrt(alpha / (2 beta)).  The sqrt-Gamma rate
  // 8 alpha (beta / (sqrt(gamma^2 + 8 alpha beta) + gamma))^2 can still fall
  // below the normal range, or the square alpha / (2 beta) of the mode pass
  // the largest double, which the check catches.
  const double gamma_abs = std::fabs(gamma);
  if (gamma_abs <= 1e150 && beta >= 1e-300 && beta <= 1e300 &&
      alpha * beta <= 1e300 && gamma_abs <= 1e300 * beta) {
    const SunAlgo1Setup s = ::build_algo1_direct(alpha, beta, gamma,
                                                 normal_only);
    if (::algo1_setup_usable(s)) return s;
  }
  // Otherwise the setup is built for Y = 2^k X, with 2^k the power of two for
  // which beta_y = beta / 4^k lies in [1, 4).  Y is MHN(alpha, beta_y,
  // gamma_y) with gamma_y = gamma / 2^k, and each draw of Y is divided by 2^k.
  // Scaling by a power of two is exact, and it commutes with the rounding of
  // every product, quotient and square root in the setup and in the draw
  // (R::rnorm and R::rgamma included), so the draws are those the direct
  // setup gives wherever its arithmetic keeps to normal doubles, the
  // logarithms that choose the proposal aside; and a Normal draw placed from
  // mu_opt + mu_lo in Y is placed so in X.  That setup squares gamma_y and
  // forms 8 alpha beta_y; while their sum is finite, every constant it forms
  // is finite, and the same check applies to it.
  const int e = std::ilogb(beta);
  const int k = (e - (e < 0)) / 2;  // floor(e / 2)
  const double y_per_x = std::ldexp(1.0, k);
  const double beta_y  = std::ldexp(beta, -2 * k);
  const double gamma_y = std::ldexp(gamma, -k);
  if (std::isfinite(gamma_y * gamma_y + 8.0 * alpha * beta_y)) {
    SunAlgo1Setup s = ::build_algo1_direct(alpha, beta_y, gamma_y, normal_only);
    if (::algo1_setup_usable(s) && std::isfinite(s.mu_opt / y_per_x) &&
        std::isfinite(s.sg_centre / y_per_x)) {
      s.alpha = alpha;
      s.beta  = beta;
      s.gamma = gamma;
      s.y_per_x = y_per_x;
      return s;
    }
  }
  // What is left is a law whose mode in Y is above 5e152, where
  // gamma_y^2 + 8 alpha beta_y overflows, so that one of its terms passes half
  // the largest double: gamma_y above 9.4e153, and the mode of Y above
  // gamma_y / (2 beta_y) > 1.1e153, or alpha above 2.8e306, and the mode of Y
  // above sqrt((alpha - 1) / (2 beta_y)) > 5.9e152.  Or else a law whose
  // mode, in X, is at the top of the double range or beyond it.  The
  // log-density of Y has curvature at least 2 beta_y >= 2, so its standard
  // deviation is at most 1/sqrt(2), below 1.2e-153 times the mode in the
  // first case: far below one unit in its last place, so that every draw is
  // the mode as computed.  The Normal proposal centred there is used, in the
  // caller's coordinates.  In the second case the standard deviation of X,
  // at most 1/sqrt(2 beta), is below 1e-146 times the mode, as
  // beta >= 4.9e-324, and the call stops if the mode itself overflows.
  const SunAlgo1Setup s = ::build_algo1_normal_limit(alpha, beta, gamma);
  if (!std::isfinite(s.mu_opt)) stop_beyond_range(alpha, beta, gamma);
  return s;
}

double sample_sun_algo1(const SunAlgo1Setup& s, int* retries_out,
                        int max_retries) {
  // X below is the variable the setup describes: the draw itself, or the
  // draw times the power of two y_per_x at the extremes handled in
  // build_sun_algo1.  to_x returns the draw, the division being exact.
  auto to_x = [&s](double X) {
    return (s.y_per_x == 1.0) ? X : X / s.y_per_x;
  };
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    if (s.chosen == SunAlgo1Setup::NORMAL) {
      // X ~ N(mu, sigma^2) with mu = mu_opt + mu_lo; reject X <= 0.  Where
      // the draw is placed from mu_opt + mu_lo, z = sigma Z is the proposal's
      // offset from mu and X the proposal rounded to double precision, with
      // mu_lo + z formed in one rounding; otherwise X is R::rnorm's draw.
      double z = 0.0, X;
      if (s.offset_normal) {
        const double std_normal = R::norm_rand();
        z = s.sigma * std_normal;
        X = s.mu_opt + mhn::fma(s.sigma, std_normal, s.mu_lo);
      } else {
        X = R::rnorm(s.mu_opt, s.sigma);
      }
      if (X > 0.0) {
        // Step 6 of Algorithm 1: log U < (alpha-1) log(X/mu)
        // + (2 beta mu - gamma) (mu - X).  mu = mu_opt solves
        // 2 beta mu^2 - gamma mu - (alpha-1) = 0, so 2 beta mu - gamma is
        // (alpha-1) / mu and the right side is (alpha-1) [log(1+u) - u] with
        // u = (X - mu) / mu.  As two terms, each of size (alpha-1) |u|, it
        // cancels to order one; and log(X/mu) carries an absolute rounding of
        // about eps, which the factor alpha - 1 magnifies to order one from
        // alpha of about 1e15.  X - mu carries a single rounding relative to
        // itself, and log1p_minus_u keeps the rest.  That form is used once
        // alpha - 1 exceeds 1e6, and the form as printed at or below that
        // (see build_algo1_direct).  Where the draw is placed from
        // mu_opt + mu_lo, u = z / mu is the proposal's own offset, free of
        // the rounding of X, and the test is that of the proposal before its
        // rounding to double precision.
        const double log_u = std::log(R::runif(0.0, 1.0));
        const double log_acc =
          s.offset_normal
            ? (s.alpha - 1.0) * mhn::log1p_minus_u(z / s.mu_opt)
          : s.stable_normal
            ? (s.alpha - 1.0) * mhn::log1p_minus_u((X - s.mu_opt) / s.mu_opt)
            : (s.alpha - 1.0) * std::log(X / s.mu_opt)
                + s.slope * (s.mu_opt - X);
        if (log_u <= log_acc) {
          if (retries_out != nullptr) *retries_out += retries;
          return to_x(X);
        }
      }
    } else {
      // T ~ Gamma(alpha/2, rate = delta_opt); X = sqrt(T).
      const double T = R::rgamma(0.5 * s.alpha, 1.0 / s.delta_opt);
      const double X = std::sqrt(T);
      // Step 16 as the perfect square formed in build_sun_algo1.
      const double d = X - s.sg_centre;
      const double log_u = std::log(R::runif(0.0, 1.0));
      const double log_acc = -s.sg_curv * d * d;
      if (log_u <= log_acc) {
        if (retries_out != nullptr) *retries_out += retries;
        return to_x(X);
      }
    }
    ++retries;
  }
  // Retry budget exhausted.  Report failure as NaN and let the caller decide
  // how to surface it; raising an R warning from here would leave this frame
  // by a long jump under options(warn = 2).
  if (retries_out != nullptr) *retries_out += retries;
  return std::numeric_limits<double>::quiet_NaN();
}

SunAlgo3Setup build_sun_algo3(double alpha, double beta, double gamma) {
  if (gamma > 0.0) {
    Rcpp::stop("Sun Algorithm 3 requires gamma <= 0 (got %g).", gamma);
  }
  if (!(alpha > 0.0)) {
    Rcpp::stop("Sun Algorithm 3 requires alpha > 0 (got %g).", alpha);
  }
  if (!(beta > 0.0)) {
    Rcpp::stop("Sun Algorithm 3 requires beta > 0 (got %g).", beta);
  }
  // The Gamma proposal has shape alpha r, with r in [1/2, 1) (1/2 at
  // gamma = 0; it can round to 1 in double), and so at most alpha; beyond
  // kRgammaShapeMax R::rgamma cannot draw it closely enough.
  if (alpha > kRgammaShapeMax) ::stop_algo3_shape(alpha, beta, gamma);

  SunAlgo3Setup s;
  s.alpha     = alpha;
  s.beta      = beta;
  s.gamma_abs = std::abs(gamma);

  // m_init heuristic + a single Newton-Raphson refinement step on l(m)
  // (Sun et al. 2023, Section 4.2 and Supplementary Section 2.17), used at
  // every alpha > 0 although Theorem 4 is stated for alpha > 1 (see
  // sun_a3_l_prime).
  //
  // Both are computed for the kernel at beta = 1, that is for
  // MHN(alpha, 1, Delta) with Delta = gamma / sqrt(beta), and the matching
  // point u found there is then divided by sqrt(beta).  sqrt(beta) X is
  // MHN(alpha, 1, Delta), and the acceptance probability A_neg(m) of
  // Theorem 4a depends on (beta, gamma, m) only through Delta and
  // sqrt(beta) m, so this is the same point.  Working at beta = 1 also keeps
  // l'(m) and l''(m), whose terms carry beta^2 and (2 beta m + |gamma|)^3,
  // inside the double range: at beta beyond about 1e154 those overflow and the
  // Newton step would be skipped, and below 1e-154 beta^2 underflows and the
  // step would use a wrong l''.  r, the shape and the rate are functions of
  // (alpha, Delta, u) alone, and they are formed from u directly.
  //
  // At beta = 1 the setup still squares |Delta| (in the mode, and as
  // (2 u + |Delta|)^2 in l'(u) and l''(u)) and, for alpha > 1.1, divides by
  // u^2 with u = sqrt(beta) m of order alpha / |Delta|; both leave the double
  // range from |Delta| of about 1e154, and |Delta| itself overflows once
  // |gamma| / sqrt(beta) passes the largest double.  Past |Delta| = 1e150 the
  // unit of length is therefore the scale of the Gamma(alpha, rate |gamma|)
  // law that the MHN law approaches, alpha / |gamma| in x and
  // w = alpha / |Delta| in u: in v = u / w the kernel is
  // v^(alpha-1) exp(-w^2 v^2 - alpha v), whose terms are all of order alpha,
  // and m, r, the rate and p and q below are formed from v, never from
  // |Delta|.  There w v is far below |Delta|, so
  // r = (w^2 v + alpha) / (2 w^2 v + alpha) is 1 to double precision, and
  // Algorithm 3 proposes the Gamma limit itself and accepts it, whatever the
  // matching point.  If w^2 underflows to 0, l'' is 0 and the Newton step is
  // skipped.  For alpha <= 1.1 the matching point is the mean of the limit,
  // v = 1, which keeps the rate at about alpha: Sun et al.'s closed form,
  // alpha^2 / (1 + alpha), divided by sqrt(beta) as in m_init_algo3, would
  // put it alpha |Delta| / (1 + alpha) Gamma scales out, where the rate, of
  // order |Delta|, reaches the top of the double range, and it cannot be
  // formed at all once w underflows to 0.
  //
  // The rate of the Gamma proposal, u (u + |Delta|) or v (w^2 v + alpha),
  // has to be a normal double: the sampler divides by it (in t_min and in the
  // conditional draw below t_min), and a rate that has lost its digits or
  // underflowed no longer matches the matching point it was formed from.  At
  // a tiny shape the closed form alpha^2 / (1 + alpha) and the Newton step
  // can both put it far below the normal range (at a unit tilt the closed
  // form does so for alpha below about 1e-154), and a matching point that
  // underflows to 0 leaves no proposal at all.  Theorem 3 of Sun et al. gives
  // a valid envelope for any matching point m > 0, so the point is raised,
  // where necessary, to the one at which the rate is kRateFloor = 2^-1021,
  // twice the smallest normal double; wherever the rate is at least that
  // already, nothing changes.  That point is about kRateFloor / |Delta|.
  // Rounded to a subnormal double, as it is from |Delta| of about 2 on, it
  // still keeps the rate above 2/3 of kRateFloor, a normal double.  From
  // |Delta| of about 1.8e16 on it lies below half the smallest positive
  // double, 2^-1074, and would round to 0; that double is taken instead, and
  // it already puts the rate, 2^-1074 |Delta|, above kRateFloor.  The law of
  // X then puts all but a fraction of order 1e3 alpha of its mass on the
  // values that round to 0, and so do the proposals.
  constexpr double kRateFloor = 2.0 * std::numeric_limits<double>::min();
  const double sqrt_beta = std::sqrt(beta);
  const double tilt_abs = s.gamma_abs / sqrt_beta;
  if (tilt_abs <= 1e150) {
    // The positive root of u (u + |Delta|) = kRateFloor, in the form that
    // does not cancel, or the smallest positive double where that root is
    // smaller still.
    const double u_floor = std::max(
      2.0 * kRateFloor
      / (tilt_abs + std::sqrt(tilt_abs * tilt_abs + 4.0 * kRateFloor)),
      std::numeric_limits<double>::denorm_min());
    double u_init =
      ::m_init_algo3(alpha, 1.0, -tilt_abs, &s.used_inflex_heuristic);
    if (!(u_init >= u_floor)) u_init = u_floor;
    double u = ::newton_refine_m(u_init, alpha, 1.0, tilt_abs);
    if (!(u >= u_floor)) u = u_init;
    s.m_init = u_init / sqrt_beta;
    s.m      = u / sqrt_beta;
    const double u_g = u + tilt_abs;
    s.r      = u_g / (2.0 * u + tilt_abs);
    s.rate   = u * u_g;
    s.acc_p  = u / u_g;
    s.acc_q  = tilt_abs / u_g;
  } else {
    const double w = alpha / tilt_abs;  // 0 if tilt_abs has overflowed
    const double w2 = w * w;
    double v_init, v;
    if (alpha > 1.1) {
      v_init = ::m_init_algo3(alpha, w2, -alpha, &s.used_inflex_heuristic);
      v = ::newton_refine_m(v_init, alpha, w2, alpha);
    } else {
      // The rate v (w^2 v + alpha) is at least alpha v, so v = 1 keeps it at
      // or above kRateFloor except for a subnormal alpha.
      s.used_inflex_heuristic = false;
      v_init = v = (alpha >= kRateFloor) ? 1.0 : kRateFloor / alpha;
    }
    const double scale = alpha / s.gamma_abs;
    s.m_init = scale * v_init;
    s.m      = scale * v;
    // u + |Delta| = a_v / w, with a_v = w^2 v + alpha.
    const double a_v = w2 * v + alpha;
    s.r      = a_v / (w2 * v + a_v);
    s.rate   = v * a_v;
    s.acc_p  = w2 * v / a_v;
    s.acc_q  = alpha / a_v;
  }
  // m is finite for every shape up to kRgammaShapeMax: in the first branch it
  // is at most about sqrt(alpha / (2 beta)) + 1 / sqrt(beta), below 1e171 for
  // any beta >= 4.9e-324, and in the second it is of order
  // alpha / |gamma| < 1e30.  So is the whole law of X.
  if (!std::isfinite(s.m)) {
    Rcpp::stop("rmhn: internal invariant violated (the Sun Algorithm 3 "
               "matching point is not finite at alpha = %g, beta = %g, "
               "gamma = %g). Please report this to the package maintainer.",
               alpha, beta, gamma);
  }
  s.shape        = alpha * s.r;
  s.t_min        = mhn::gamma_precision_floor(1.0 / s.rate);
  // The acceptance exponent as printed carries a rounding of a few eps times
  // the rate, below 1e-9 up to a rate of 1e6 (a shape of order 1e6); there it
  // is kept, as it needs no transcendental function beyond the draw itself.
  s.stable_test  = (s.rate > 1e6);
  return s;
}

double sample_sun_algo3(const SunAlgo3Setup& s, int* retries_out,
                        int max_retries) {
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    // T ~ Gamma(shape, rate); R::rgamma takes (shape, scale=1/rate).
    double T = R::rgamma(s.shape, 1.0 / s.rate);
    double X;
    if (T >= s.t_min) {
      X = s.m * std::pow(T, s.r);
    } else {
      // T has lost its precision or underflowed to zero, and X = m T^r would
      // inherit that; draw it again from its conditional law below t_min, in
      // logarithms (see gamma_precision_floor in mhn_sun.h).  X is formed in
      // logarithms too: T^r can lie below the double range while m T^r, with
      // a large m, does not.  T itself enters only the printed acceptance
      // test, used up to a rate of 1e6, as rate T < DBL_MIN max(1, rate),
      // which is negligible there, so its underflow is harmless.
      const double log_t = mhn::gamma_tail_log_draw(s.t_min, s.shape, s.rate);
      T = std::exp(log_t);
      X = std::exp(std::log(s.m) + s.r * log_t);
    }
    // Acceptance (Algorithm 3, step 5): log U <= rate (X/m)^(1/r)
    // - beta X^2 - |gamma| X.  Since X = m T^r (step 3), (X/m)^(1/r) is T.
    const double log_u   = std::log(R::runif(0.0, 1.0));
    const double log_acc =
      s.stable_test
        ? ::a3_log_acc_stable((X - s.m) / s.m, s.acc_p, s.acc_q, s.rate)
        : s.rate * T - s.beta * X * X - s.gamma_abs * X;
    if (log_u <= log_acc) {
      if (retries_out != nullptr) *retries_out += retries;
      return X;
    }
    ++retries;
  }
  // Retry budget exhausted; see the note in sample_sun_algo1.
  if (retries_out != nullptr) *retries_out += retries;
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace mhn

// =====================================================================
// Entry points used by the test suite to exercise each sampler in
// isolation, bypassing the rmhn() dispatcher.  Not part of the public API.
// =====================================================================

// [[Rcpp::export(.rmhn_sun_algo1_cpp)]]
Rcpp::NumericVector rmhn_sun_algo1_cpp(int n, double alpha, double beta, double gamma) {
  if (n < 0) Rcpp::stop("n must be non-negative");
  if (n == 0) return Rcpp::NumericVector(0);
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo1Setup s = mhn::build_sun_algo1(alpha, beta, gamma);
  Rcpp::NumericVector out(n);
  int retries = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = mhn::sample_sun_algo1(s, &retries);
  }
  out.attr("sun_retries")  = retries;
  out.attr("sun_proposal") = (s.chosen == mhn::SunAlgo1Setup::NORMAL)
                             ? "normal" : "sqrt_gamma";
  return out;
}

// [[Rcpp::export(.dump_sun_algo1_cpp, rng = false)]]
Rcpp::List dump_sun_algo1_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo1Setup s = mhn::build_sun_algo1(alpha, beta, gamma);
  // The constants of the law of X, whichever variable the setup describes.
  // y_per_x is 2^k, and each conversion is exact unless its result is
  // subnormal.
  const int k = std::ilogb(s.y_per_x);
  return Rcpp::List::create(
    Rcpp::Named("chosen")           = (s.chosen == mhn::SunAlgo1Setup::NORMAL)
                                      ? "normal" : "sqrt_gamma",
    Rcpp::Named("mu_opt")           = std::ldexp(s.mu_opt, -k),
    Rcpp::Named("delta_opt")        = std::ldexp(s.delta_opt, 2 * k),
    Rcpp::Named("sigma")            = std::ldexp(s.sigma, -k),
    Rcpp::Named("log_K1_minus_K2")  = s.log_K1_minus_K2
  );
}

// [[Rcpp::export(.rmhn_sun_algo3_cpp)]]
Rcpp::NumericVector rmhn_sun_algo3_cpp(int n, double alpha, double beta, double gamma) {
  if (n < 0) Rcpp::stop("n must be non-negative");
  if (n == 0) return Rcpp::NumericVector(0);
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo3Setup s = mhn::build_sun_algo3(alpha, beta, gamma);
  Rcpp::NumericVector out(n);
  int retries = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = mhn::sample_sun_algo3(s, &retries);
  }
  out.attr("sun_retries")   = retries;
  out.attr("sun_m")         = s.m;
  out.attr("sun_r")         = s.r;
  out.attr("sun_used_inflex") = s.used_inflex_heuristic;
  return out;
}

// [[Rcpp::export(.dump_sun_algo3_cpp, rng = false)]]
Rcpp::List dump_sun_algo3_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo3Setup s = mhn::build_sun_algo3(alpha, beta, gamma);
  return Rcpp::List::create(
    Rcpp::Named("m")                     = s.m,
    Rcpp::Named("m_init")                = s.m_init,
    Rcpp::Named("r")                     = s.r,
    Rcpp::Named("shape")                 = s.shape,
    Rcpp::Named("rate")                  = s.rate,
    Rcpp::Named("used_inflex_heuristic") = s.used_inflex_heuristic
  );
}
