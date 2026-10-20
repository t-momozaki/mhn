// Adaptive Gauss-Kronrod quadrature whose accuracy does not depend on the
// width of the interval.
//
// Boost's gauss_kronrod<Real, N>::integrate(f, a, b, max_depth, tol) maps
// [a, b] onto [-1, 1] and bisects while the Kronrod-Gauss difference of the
// mapped integral K exceeds tol * (b - a)/2 * |K|.  The factor (b - a)/2 means
// the relative accuracy it demands is tol times the half-width of the panel,
// not tol.  A wide panel is accepted early, so the answer is less accurate
// than asked.  A narrow one is refined further than asked, and one whose
// half-width is below 2 eps / tol -- 3e-8 for tol = sqrt(eps) -- can never
// pass, because the difference is floored at 2 eps |K|.  It then bisects all
// the way to max_depth: 2^16 - 1 panels of N evaluations each, about a million
// for the 15-point rule and two million for the 31-point one.
//
// The panels in this package are set by the scale of the kernel: 1 / sqrt(beta)
// on the x axis, the distance over which a strongly tilted kernel dies, which
// shrinks like 1 / |z|, and that distance again on the log axis.  The
// quantities computed from them -- log Psi[alpha/2, z], F(q / sqrt(beta)),
// the skewness and kurtosis, and the mean and variance in units of
// 1 / sqrt(beta) and 1 / beta -- depend on beta and gamma only through
// z = gamma / sqrt(beta) (Sun et al. 2023, Lemmas 1a and 2a).  Integrated on
// the raw interval both their accuracy and their cost would move with beta,
// and with the tilt.  Handing Boost the interval [-1, 1] itself,
//
//     int_a^b f(x) dx = (h/2) int_{-1}^{1} f(c + (h/2) v) dv,
//     h = b - a,  c = (a + b) / 2,
//
// makes the factor (b - a)/2 in its test exactly 1, so every panel meets the
// relative standard tol that the caller asked for, whatever its width, and the
// result is the same at every scale.

#ifndef MHN_QUADRATURE_H
#define MHN_QUADRATURE_H

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace mhn {

// int_a^b f(x) dx by the N-point Gauss-Kronrod rule, bisected adaptively to a
// depth of at most 15 on the interval of unit half-width, [-1, 1].  An empty
// or reversed interval contributes 0.
template <unsigned N, class F>
double gk_unit(const F& f, double a, double b, double tol) {
  const double h = b - a;
  if (!(h > 0.0)) return 0.0;
  const double half = 0.5 * h;
  const double mid = a + half;
  auto on_unit = [&f, mid, half](double v) -> double {
    return f(mid + half * v);
  };
  return half * boost::math::quadrature::gauss_kronrod<double, N>::integrate(
                    on_unit, -1.0, 1.0, /*max_depth=*/15, tol);
}

}  // namespace mhn

#endif  // MHN_QUADRATURE_H
