// Sun et al. (2023) Algorithm 1 (gamma > 0, alpha > 1) and Algorithm 3
// (gamma <= 0, alpha > 0) for sampling from the Modified Half-Normal.
//
// Algorithm 2 (gamma > 0, alpha <= 1) is intentionally NOT implemented.
// It breaks for gamma < 0 (negative mixing weights) and is highly
// inefficient for alpha < 1; for that parameter regime the dispatcher
// in rmhn_cpp routes to the Gao & Wang (2025) RTDR sampler in
// mhn_rtdr.cpp instead.

#include "mhn_sun.h"
#include "mhn_check.h"

#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <cmath>
#include <limits>

namespace {

// =====================================================================
// Algorithm 1 helpers
// =====================================================================

// log K_1 - log K_2 with Psi cancelled (it appears as a common factor in
// both K_1(mu) and K_2(delta), per theory_sun.md sec 5.1 / Theorem 1a).
// Sign convention: < 0 -> use NORMAL proposal, >= 0 -> use SQRT_GAMMA.
double log_K1_minus_K2(double alpha, double beta, double gamma,
                       double mu, double delta) {
  // log K_1 = log(2 sqrt pi)
  //          + (alpha-1) [ 0.5 log(beta(alpha-1)) - log(2 beta mu - gamma) ]
  //          - (alpha-1) + beta mu^2
  //          - log Psi    (cancels)
  const double log_K1 = std::log(2.0) + 0.5 * std::log(M_PI)
    + (alpha - 1.0) *
      ( 0.5 * std::log(beta * (alpha - 1.0)) - std::log(2.0 * beta * mu - gamma) )
    - (alpha - 1.0) + beta * mu * mu;
  // log K_2 = (alpha/2) log beta + log Gamma(alpha/2)
  //          + gamma^2 / (4 (beta - delta))
  //          - (alpha/2) log delta
  //          - log Psi    (cancels)
  const double log_K2 = 0.5 * alpha * std::log(beta) + std::lgamma(0.5 * alpha)
    + gamma * gamma / (4.0 * (beta - delta))
    - 0.5 * alpha * std::log(delta);
  return log_K1 - log_K2;
}

// =====================================================================
// Algorithm 3 helpers
// =====================================================================

// Right inflection point of the unnormalized MHN density f(x) for alpha > 1.
// Inflection: f''(x) = 0 <=> L'(x)^2 + L''(x) = 0 (since f > 0), where
//   L(x) = (alpha-1) log x - beta x^2 + gamma x.
// At the mode L'(x_mode) = 0, so F(x_mode) = L''(x_mode) < 0.  As x -> inf,
// F(x) -> +inf, so a single root exists in (x_mode, inf).  Bisection.
double right_inflection_point(double alpha, double beta, double gamma_signed,
                              double x_mode) {
  if (!(alpha > 1.0)) {
    Rcpp::stop("right_inflection_point requires alpha > 1.");
  }
  auto F = [alpha, beta, gamma_signed](double x) -> double {
    if (x <= 0.0) return std::numeric_limits<double>::infinity();
    const double Lp  = (alpha - 1.0) / x - 2.0 * beta * x + gamma_signed;
    const double Lpp = -(alpha - 1.0) / (x * x) - 2.0 * beta;
    return Lp * Lp + Lpp;
  };
  double lo = x_mode * (1.0 + 1e-8);
  if (F(lo) >= 0.0) {
    // Defensive: should not happen analytically (F(x_mode) = L''(x_mode) < 0),
    // but very small x_mode + numerical noise can push it positive.  Walk
    // back toward x_mode.
    lo = x_mode + std::numeric_limits<double>::epsilon();
  }
  double hi = x_mode + 5.0 / std::sqrt(2.0 * beta) + 1.0;
  int doublings = 0;
  while (F(hi) < 0.0 && doublings < 50) {
    hi *= 2.0;
    ++doublings;
  }
  if (F(hi) < 0.0) {
    Rcpp::stop("right_inflection_point: failed to bracket.");
  }
  for (int it = 0; it < 100; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (hi - lo < 1e-12 * (1.0 + std::abs(mid))) return mid;
    if (F(mid) < 0.0) lo = mid;
    else              hi = mid;
  }
  return 0.5 * (lo + hi);
}

// Initial matching point for Algorithm 3.  Both branches and the
// alpha <= 1.1 / alpha > 1.1 cutoff are taken verbatim from the
// m_init definition in Sun et al. (2023) Section 4.2 (the equation
// for m_init following Theorem 4; Theorem 4 itself only states the
// acceptance-probability bounds).  For alpha > 1.1 we use the
// inflection-point heuristic (a lambda-weighted mix of the mode and
// the right-inflection point); for alpha <= 1.1 we fall back to the
// simple closed form m = alpha^2 / (1 + alpha).
// `used_inflex_out` reports which path was taken (for diagnostics).
double m_init_algo3(double alpha, double beta, double gamma,
                    bool* used_inflex_out) {
  if (alpha > 1.1) {
    // X_mode for f(x) ∝ x^(α-1) exp(γx - βx²); same closed form as
    // mhn_mode_cpp.  Valid for alpha > 1 (else mode at boundary).
    const double x_mode = (gamma + std::sqrt(gamma * gamma
                                             + 8.0 * beta * (alpha - 1.0)))
                          / (4.0 * beta);
    const double x_inflex = ::right_inflection_point(alpha, beta, gamma, x_mode);
    const double x_left = 2.0 * x_mode - x_inflex;
    if (x_left > 0.0) {
      // lambda = f̃(x_left) / (f̃(x_left) + f̃(x_inflex)) via log-space sigmoid.
      auto log_f_tilde = [alpha, beta, gamma](double x) -> double {
        return (alpha - 1.0) * std::log(x) + gamma * x - beta * x * x;
      };
      const double diff = log_f_tilde(x_left) - log_f_tilde(x_inflex);
      const double lambda = 1.0 / (1.0 + std::exp(-diff));
      if (used_inflex_out != nullptr) *used_inflex_out = true;
      return 1.5 * lambda * x_mode + (1.0 - 1.5 * lambda) * x_inflex;
    }
    // x_left <= 0: f̃(x_left) undefined; fall back to simple heuristic.
  }
  if (used_inflex_out != nullptr) *used_inflex_out = false;
  return alpha * alpha / (1.0 + alpha);
}

// First derivative l'(m) of the log-acceptance for Sun Algorithm 3,
// from Sun et al. (2023) Supplementary Section 2.17.  Maximising
// l(m) = log(A_neg(m, alpha, beta, gamma)) over m gives the matching
// point that maximises the rejection-sampler acceptance probability.
//
// gamma_abs = |gamma| (gamma <= 0 is required for Algorithm 3).
double sun_a3_l_prime(double m, double alpha, double beta, double gamma_abs) {
  const double A = beta * m + gamma_abs;        // beta m + |gamma|
  const double B = 2.0 * beta * m + gamma_abs;  // 2 beta m + |gamma|
  const double C = beta * m * m + m * gamma_abs; // beta m^2 + m |gamma|
  const double shape_term = alpha * A / B;
  const double factor = alpha * beta * gamma_abs / (B * B);
  return factor * (boost::math::digamma(shape_term) - std::log(C))
         + 2.0 * beta / B - beta / A;
}

// Second derivative l''(m), also from Sun et al. (2023) Supplementary
// Section 2.17.  Used as the Hessian in Newton-Raphson for m_recommend.
double sun_a3_l_double_prime(double m, double alpha, double beta,
                             double gamma_abs) {
  const double A = beta * m + gamma_abs;
  const double B = 2.0 * beta * m + gamma_abs;
  const double C = beta * m * m + m * gamma_abs;
  const double shape_term = alpha * A / B;
  const double psi  = boost::math::digamma(shape_term);
  const double psi1 = boost::math::trigamma(shape_term);
  const double diff = psi - std::log(C);
  const double T1 = (-4.0 * alpha * beta * beta * gamma_abs / (B * B * B))
                    * diff;
  const double f2 = alpha * beta * gamma_abs / (B * B);
  const double T2 = f2 * (-f2 * psi1 - B / C);
  // Trailing terms of l''(m): -1/(m + |gamma|/(2 beta))^2
  //                            + 1/(m + |gamma|/beta)^2,
  // equivalent to -(2 beta / B)^2 + (beta / A)^2.
  const double T3 = -(2.0 * beta / B) * (2.0 * beta / B)
                    + (beta / A) * (beta / A);
  return T1 + T2 + T3;
}

// One Newton-Raphson step on l(m), as in Sun et al. (2023)
// Supplementary Section 2.17:
//   m_recommend = m_init - l'(m_init) / l''(m_init).
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

}  // namespace

namespace mhn {

SunAlgo1Setup build_sun_algo1(double alpha, double beta, double gamma) {
  if (!(alpha > 1.0)) {
    Rcpp::stop("Sun Algorithm 1 requires alpha > 1 (got %g).", alpha);
  }
  if (!(gamma > 0.0)) {
    Rcpp::stop("Sun Algorithm 1 requires gamma > 0 (got %g).", gamma);
  }
  if (!(beta > 0.0)) {
    Rcpp::stop("Sun Algorithm 1 requires beta > 0 (got %g).", beta);
  }

  SunAlgo1Setup s;
  s.alpha = alpha;
  s.beta  = beta;
  s.gamma = gamma;
  // mu_opt = (gamma + sqrt(gamma^2 + 8 (alpha-1) beta)) / (4 beta)
  s.mu_opt = (gamma + std::sqrt(gamma * gamma + 8.0 * (alpha - 1.0) * beta))
             / (4.0 * beta);
  // delta_opt = beta + (gamma^2 - gamma sqrt(gamma^2 + 8 alpha beta)) / (4 alpha)
  s.delta_opt = beta + (gamma * gamma - gamma *
                        std::sqrt(gamma * gamma + 8.0 * alpha * beta))
                       / (4.0 * alpha);
  s.sigma = 1.0 / std::sqrt(2.0 * beta);

  // Proposal selection: K_2 > K_1 (i.e. log K_1 - log K_2 < 0) -> Normal,
  // otherwise -> sqrt-Gamma.  K_1 / K_2 are the per-proposal rejection
  // constants for the Normal and sqrt-Gamma envelopes respectively
  // (Sun et al. 2023 Section 5.1, Theorem 1).
  s.log_K1_minus_K2 = ::log_K1_minus_K2(alpha, beta, gamma,
                                        s.mu_opt, s.delta_opt);
  s.chosen = (s.log_K1_minus_K2 < 0.0)
             ? SunAlgo1Setup::NORMAL
             : SunAlgo1Setup::SQRT_GAMMA;
  return s;
}

double sample_sun_algo1(const SunAlgo1Setup& s, int* retries_out) {
  const int max_retries = 1000;
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    if (s.chosen == SunAlgo1Setup::NORMAL) {
      // X ~ N(mu_opt, sigma); reject X <= 0.
      const double X = R::rnorm(s.mu_opt, s.sigma);
      if (X > 0.0) {
        // log U < (alpha-1) log(X/mu) + (2 beta mu - gamma) (mu - X)
        const double log_u = std::log(R::runif(0.0, 1.0));
        const double log_acc =
          (s.alpha - 1.0) * std::log(X / s.mu_opt)
          + (2.0 * s.beta * s.mu_opt - s.gamma) * (s.mu_opt - X);
        if (log_u <= log_acc) {
          if (retries_out != nullptr) *retries_out += retries;
          return X;
        }
      }
    } else {
      // T ~ Gamma(alpha/2, rate = delta_opt); X = sqrt(T).
      const double T = R::rgamma(0.5 * s.alpha, 1.0 / s.delta_opt);
      const double X = std::sqrt(T);
      // log U < -(beta - delta) X^2 + gamma X - gamma^2 / (4 (beta - delta))
      const double bd = s.beta - s.delta_opt;
      const double log_u = std::log(R::runif(0.0, 1.0));
      const double log_acc =
        -bd * X * X + s.gamma * X - s.gamma * s.gamma / (4.0 * bd);
      if (log_u <= log_acc) {
        if (retries_out != nullptr) *retries_out += retries;
        return X;
      }
    }
    ++retries;
  }
  Rcpp::warning("sample_sun_algo1: max retries (%d) exceeded.", max_retries);
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

  SunAlgo3Setup s;
  s.alpha     = alpha;
  s.beta      = beta;
  s.gamma_abs = std::abs(gamma);

  // m_init heuristic + a single Newton-Raphson refinement step on l(m)
  // (Sun et al. 2023 Section 5.3 / Supplementary Section 2.17).
  s.m_init = ::m_init_algo3(alpha, beta, gamma, &s.used_inflex_heuristic);
  s.m      = ::newton_refine_m(s.m_init, alpha, beta, s.gamma_abs);

  const double bm_g = beta * s.m + s.gamma_abs;
  s.r            = bm_g / (2.0 * beta * s.m + s.gamma_abs);
  s.shape        = alpha * s.r;
  s.rate         = s.m * bm_g;
  s.m_betam_gam  = s.m * bm_g;
  return s;
}

double sample_sun_algo3(const SunAlgo3Setup& s, int* retries_out) {
  const int max_retries = 1000;
  int retries = 0;
  for (int iter = 0; iter < max_retries; ++iter) {
    // T ~ Gamma(shape, rate); R::rgamma takes (shape, scale=1/rate).
    const double T = R::rgamma(s.shape, 1.0 / s.rate);
    const double X = s.m * std::pow(T, s.r);
    // Acceptance: log U <= m_betam_gam (X/m)^(1/r) - beta X^2 - |gamma| X
    const double pow_term = std::pow(X / s.m, 1.0 / s.r);
    const double log_u   = std::log(R::runif(0.0, 1.0));
    const double log_acc = s.m_betam_gam * pow_term
                           - s.beta * X * X
                           - s.gamma_abs * X;
    if (log_u <= log_acc) {
      if (retries_out != nullptr) *retries_out += retries;
      return X;
    }
    ++retries;
  }
  Rcpp::warning("sample_sun_algo3: max retries (%d) exceeded.", max_retries);
  if (retries_out != nullptr) *retries_out += retries;
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace mhn

// =====================================================================
// Test-only Rcpp exports (for region-by-region validation ahead of the
// full rmhn() dispatcher in Step 3.4).
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

// [[Rcpp::export(.dump_sun_algo1_cpp)]]
Rcpp::List dump_sun_algo1_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo1Setup s = mhn::build_sun_algo1(alpha, beta, gamma);
  return Rcpp::List::create(
    Rcpp::Named("chosen")           = (s.chosen == mhn::SunAlgo1Setup::NORMAL)
                                      ? "normal" : "sqrt_gamma",
    Rcpp::Named("mu_opt")           = s.mu_opt,
    Rcpp::Named("delta_opt")        = s.delta_opt,
    Rcpp::Named("sigma")            = s.sigma,
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

// [[Rcpp::export(.dump_sun_algo3_cpp)]]
Rcpp::List dump_sun_algo3_cpp(double alpha, double beta, double gamma) {
  mhn::check_params_scalar(alpha, beta, gamma);
  mhn::SunAlgo3Setup s = mhn::build_sun_algo3(alpha, beta, gamma);
  return Rcpp::List::create(
    Rcpp::Named("m")                     = s.m,
    Rcpp::Named("m_init")                = s.m_init,
    Rcpp::Named("r")                     = s.r,
    Rcpp::Named("shape")                 = s.shape,
    Rcpp::Named("rate")                  = s.rate,
    Rcpp::Named("m_betam_gam")           = s.m_betam_gam,
    Rcpp::Named("used_inflex_heuristic") = s.used_inflex_heuristic
  );
}
