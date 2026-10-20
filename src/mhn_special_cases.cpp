// Special-case detection and density dispatch.

#include "mhn_special_cases.h"
#include "mhn_constants.h"
#include "mhn_fma.h"
#include "mhn_stable.h"
#include "mhn_warn.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace mhn {

namespace {

// sqrt(DBL_MIN): below this t, t * t is no longer a normal double.
const double kSqrtDblMin = std::sqrt(DBL_MIN);

const double kSqrtPi = 1.772453850905516027298;

// The scaled complementary error function erfcx(x) = exp(x^2) erfc(x) for
// x >= 3/sqrt2 = 2.12, which is where the truncated normal needs it (z0 >= 3,
// x = z/sqrt2).  It is formed as g(x) / (sqrt(pi) x), where
// g(x) = sqrt(pi) x erfcx(x) lies between 0.91 and 1 there and tends to 1 as
// x grows, and is returned as erfcx(x) / c, the division by c folded into the
// final one.
//
// g is a polynomial in each of two variables, both mapped onto t in [-1, 1]:
//
//   * x < 5:  y = 1/x in [1/5, sqrt2/3], degree 15;
//   * x >= 5: w = 1/x^2 in [0, 1/25], degree 13, which reaches w = 0, so that
//     it holds however large x is (at w = 0 it gives 1, the limit of g).
//
// Each is the leading part of the Chebyshev expansion of g in its variable,
// computed in 256-bit arithmetic from the values of g at 80 Chebyshev points
// and rewritten in powers of t; the Chebyshev terms left out sum to 2.0e-18
// and 1.1e-19 relative to g.  In powers of t the coefficients fall off
// quickly and their absolute values sum to less than 1.1 g, so the sum cannot
// cancel, and it is formed by Estrin's scheme, whose dependent chain is four
// multiply-adds long rather than one per degree.  Against
// exp(x^2) erfc(x) in 256-bit arithmetic the result is within 5e-16
// relative, and at most 3.2 units in the last place, over [2.12, 1e150].
const double kErfcxY[16] = {
    0.95121243654685694,    -0.034516122650552233,  -0.0037204829155527762,
    0.00086239768220414798, -6.2357051040285333e-05, -7.5667576509566674e-06,
    3.0027365808900129e-06, -4.2705819765177315e-07, 5.898926663092715e-09,
    1.2743451072891351e-08, -3.5052621279689111e-09, 4.7064643227808368e-10,
    6.1646862484239256e-12, -2.1646376530488549e-11, 6.1749482945017917e-12,
    -7.9345181860728888e-13};
const double kErfcxW[14] = {
    0.99028596471731922,    -0.0094412182246203045, 0.0002604872169531841,
    -1.1575329006190239e-05, 6.9701161691256636e-07, -5.230561185139081e-08,
    4.6562054844583167e-09, -4.7600452406842527e-10, 5.462341717282794e-11,
    -6.9177105127465571e-12, 9.5058682038065619e-13, -1.4112311954679222e-13,
    2.5157956328818864e-14, -4.3386588686686501e-15};

double erfcx_tail(double x, double c = 1.0) {
  const double y = 1.0 / x;
  double g;
  if (x < 5.0) {
    // t = (2 y - (1/5 + sqrt2/3)) / (sqrt2/3 - 1/5).
    const double t = 7.3690740086820554 * y - 2.473814801736411;
    const double t2 = t * t, t4 = t2 * t2, t8 = t4 * t4;
    const double* m = kErfcxY;
    const double q0 = (m[0] + m[1] * t) + (m[2] + m[3] * t) * t2;
    const double q1 = (m[4] + m[5] * t) + (m[6] + m[7] * t) * t2;
    const double q2 = (m[8] + m[9] * t) + (m[10] + m[11] * t) * t2;
    const double q3 = (m[12] + m[13] * t) + (m[14] + m[15] * t) * t2;
    g = (q0 + q1 * t4) + (q2 + q3 * t4) * t8;
  } else {
    // t = 50 w - 1.
    const double t = 50.0 * (y * y) - 1.0;
    const double t2 = t * t, t4 = t2 * t2, t8 = t4 * t4;
    const double* m = kErfcxW;
    const double q0 = (m[0] + m[1] * t) + (m[2] + m[3] * t) * t2;
    const double q1 = (m[4] + m[5] * t) + (m[6] + m[7] * t) * t2;
    const double q2 = (m[8] + m[9] * t) + (m[10] + m[11] * t) * t2;
    const double q3 = m[12] + m[13] * t;
    g = (q0 + q1 * t4) + (q2 + q3 * t4) * t8;
  }
  return g * y / (kSqrtPi * c);
}

// ---------------------------------------------------------------------------
// The integral I(d) = int_0^d exp(-z0 s - s^2/2) ds of the truncated normal
// (see mhn_special_cases.h), as J(d) = I(d)/d.
//
// The generating function of the Hermite polynomials,
// exp(z s - s^2/2) = sum_k He_k(z) s^k / k!, taken at s = -u gives
//
//   exp(-z0 u - u^2/2) = sum_k He_k(z0) (-u)^k / k!,
//
// and integrating over u in [0, d] and dividing by d,
//
//   J(d) = sum_k b_k / (k+1)!,     I'(d) = exp(-z0 d - d^2/2) = sum_k b_k / k!,
//
// with b_k = (-d)^k He_k(z0).  The recurrence He_{k+1}(z) = z He_k(z) -
// k He_{k-1}(z) becomes b_0 = 1, b_1 = -z0 d, b_{k+1} = -z0 d b_k - k d^2
// b_{k-1}, so no power or polynomial is formed explicitly.
//
// |He_k(z)| is bounded by the Hermite polynomial with all coefficients made
// positive, whose generating function is exp(|z| s + s^2/2), so the terms are
// bounded by those of exp(|z0| d + d^2/2) and the sum cancels by at most
// exp(2 h), h = d (|z0| + d).  The callers keep h near 1/4 or below, where
// that factor is 1.6 and twenty terms reach the last bit.  A single term can
// vanish at a zero of He_k (He_2(1) = 0), so the sum stops only after two
// consecutive negligible terms; the recurrence then bounds every later term by
// h / (k + 2) times the larger of the two before it.
// ---------------------------------------------------------------------------
const int kTnSeriesTerms = 40;

struct InverseFactorials {
  double v[kTnSeriesTerms + 2];
  InverseFactorials() {
    v[0] = 1.0;
    for (int k = 1; k < kTnSeriesTerms + 2; ++k) v[k] = v[k - 1] / k;
  }
};
const InverseFactorials kInvFact;

// Returns J(d); stores I'(d) in *slope when slope is not null.
double tn_mean_kernel(double d, double z0, double* slope) {
  const double a = -z0 * d;
  const double d2 = d * d;
  double b_prev = 1.0;
  double b = a;
  double J = 1.0 + 0.5 * a;
  double e = 1.0 + a;
  double last = std::fabs(0.5 * a);
  for (int k = 1; k < kTnSeriesTerms; ++k) {
    const double b_next = a * b - k * d2 * b_prev;
    b_prev = b;
    b = b_next;
    const double term = b * kInvFact.v[k + 2];
    J += term;
    e += b * kInvFact.v[k + 1];
    const double size = std::fabs(term);
    if (size <= 1e-17 * J && last <= 1e-17 * J) break;
    last = size;
  }
  if (slope != nullptr) *slope = e;
  return J;
}

// Where the lower tail of the truncated normal is taken from lambda(z0) I(d).
// The alternative forms get the log tail ratio, about -h, as the difference of
// two log tail probabilities of size L (of order one, or z0^2/2 for z0 far
// below 0), and so lose about eps L / h of it relative; at h = 1/4 that is a
// few eps L, and the series above needs about twenty terms.
const double kTnSmall = 0.25;

}  // namespace

bool is_sqrt_gamma(double gamma, double beta) {
  // Test the scale-free tilt Delta = gamma / sqrt(beta).  A test on |gamma|
  // alone would discard a fully significant tilt whenever beta is small: at
  // beta = 4e-18 a gamma of 1e-8 is Delta = 5, an ordinary value, which that
  // test would send to the sqrt-Gamma closed form.
  return std::fabs(gamma) < mhn_eps() * std::sqrt(beta);
}

bool is_truncated_normal(double alpha) {
  return std::fabs(alpha - 1.0) < mhn_eps();
}

void TruncNormal::set(double beta, double gamma) {
  this->beta = beta;
  this->gamma = gamma;
  // 2 beta overflows for beta above half the largest double, and
  // gamma / (2 beta) where beta is small against |gamma|, while sigma and
  // z0 = -gamma / sqrt(2 beta) are still doubles.  sigma and mu are then
  // formed from beta itself, and z0 as -gamma sigma, which is -mu / sigma
  // since 2 beta sigma^2 = 1.  2 beta is not formed where it would overflow:
  // the test is made on beta (see opaque_value in mhn_constants.h).
  if (beta <= 0.5 * DBL_MAX) {
    mu = gamma / (2.0 * beta);
    sigma = 1.0 / std::sqrt(2.0 * beta);
  } else {
    sigma = M_SQRT1_2 / std::sqrt(beta);
    mu = (0.5 * gamma) / beta;
  }
  z0 = -mu / sigma;
  if (!std::isfinite(mu)) z0 = -gamma * sigma;
  // z0 = +Inf only where gamma / sqrt(beta) itself passes the largest
  // double under a negative tilt.  The law is then Exponential(rate -gamma)
  // times exp(-beta x^2), and wherever the Exponential's log tail -rate x is
  // a double, beta x^2 is below a relative sqrt(beta) / rate < 1e-308 of it.
  exp_limit = (z0 == R_PosInf);
  fma_offset = (z0 < -kTnLargeZ0 || (!std::isfinite(mu) && gamma > 0.0));
  plain_offset = !fma_offset && std::isfinite(mu);
  if (z0 >= 3.0) {
    erfcx0 = exp_limit ? 0.0 : erfcx_tail(z0 * M_SQRT1_2);
  } else {
    erfcx0 = R_NaN;
  }
  if (z0 <= kTnQnormMaxZ0) {
    // P(Z > z0) on the log scale: the denominator of the truncated-normal CDF.
    log_mass = R::pnorm(mu / sigma, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
    mass = R::pnorm(mu / sigma, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
    if (!std::isfinite(mu)) {
      log_mass = R::pnorm(-z0, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
      mass = R::pnorm(-z0, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
    }
  } else {
    log_mass = R_NaN;
    mass = R_NaN;
  }
  have_log_erfcx0_ = false;
  have_hazard_ = false;
  have_lower0_ = false;
}

void TruncNormal::make_hazard() const {
  if (z0 >= 3.0) {
    // Phibar(z) = exp(-z^2/2) erfcx(z/sqrt2) / 2, so
    // lambda(z) = phi(z)/Phibar(z) = sqrt(2/pi) / erfcx(z/sqrt2), with no
    // z^2/2 to cancel however far out z0 lies.
    log_hazard_ = 0.5 * std::log(2.0 / M_PI) - log_erfcx0();
  } else {
    log_hazard_ = R::dnorm(z0, 0.0, 1.0, /*log=*/1) - log_mass;
  }
  // exp(log lambda) carries eps |log lambda| of rounding, 5e-14 at
  // z0 = 1e100; past z0 = kTnLargeZ0 lambda is formed as the quotient
  // itself, to a unit or two in its last place.
  hazard_ = (z0 > kTnLargeZ0) ? std::sqrt(2.0 / M_PI) / erfcx0
                              : std::exp(log_hazard_);
  have_hazard_ = true;
}

Rcpp::NumericVector dmhn_sqrt_gamma(const Rcpp::NumericVector& x,
                                    double alpha, double beta, bool log_p) {
  SqrtGammaDensity sg;
  sg.set(alpha, beta);

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    // NA and NaN are distinct in the base R d/p/q contract, and Rcpp's
    // is_na does not separate them -- it is ISNAN -- so each is tested
    // with the R_Is* predicate that isolates it.
    if (ISNAN(xi)) {
      log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
    } else {
      log_f[i] = sg.log_at(xi);
    }
  }

  if (!log_p) exp_in_place(log_f);
  return log_f;
}

void TruncNormalDensity::set(double beta, double gamma) {
  tn.set(beta, gamma);
  log_norm = R::pnorm(tn.mu / tn.sigma, 0.0, 1.0, /*lower_tail=*/1,
                      /*log_p=*/1);
  log_sigma = std::log(tn.sigma);
  log_hazard = (tn.z0 > kTnLargeZ0 && !tn.exp_limit) ? tn.log_hazard() : R_NaN;
  if (tn.exp_limit) {
    form = 3;
  } else if (tn.z0 > kTnLargeZ0) {
    form = 1;
  } else if (!tn.plain_offset) {
    form = 2;
  } else {
    form = 0;
  }
}

double TruncNormalDensity::log_at(double x) const {
  // Evaluated at x = 0 as well; the truncated-normal density is finite and
  // positive at the lower boundary.
  if (x >= 0.0) {
    if (form == 0) return R::dnorm(x, tn.mu, tn.sigma, /*log=*/1) - log_norm;
    if (form == 1) {
      const double d = x / tn.sigma;
      return log_hazard - log_sigma - d * (tn.z0 + 0.5 * d);
    }
    if (form == 2) {
      return R::dnorm(tn.zq(x), 0.0, 1.0, /*log=*/1) - log_sigma
             - tn.log_mass;
    }
    if (x == R_PosInf) return R_NegInf;
    return std::log(-tn.gamma) + tn.gamma * x;
  }
  return R_NegInf;
}

Rcpp::NumericVector dmhn_truncated_normal(const Rcpp::NumericVector& x,
                                          double beta, double gamma,
                                          bool log_p) {
  TruncNormalDensity density;
  density.set(beta, gamma);

  const R_xlen_t n = x.size();
  Rcpp::NumericVector log_f(n, R_NegInf);

  if (density.form == 0) {
    // The direct form, with its constants held in locals, as the loop the
    // recycled path shares element by element through log_at.
    const double mu = density.tn.mu;
    const double sigma = density.tn.sigma;
    const double log_norm_const = density.log_norm;
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (ISNAN(xi)) {
        log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
      } else if (xi >= 0.0) {
        log_f[i] = R::dnorm(xi, mu, sigma, /*log=*/1) - log_norm_const;
      }
      // xi < 0 falls through to R_NegInf
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (ISNAN(xi)) {
        log_f[i] = R_IsNA(xi) ? NA_REAL : R_NaN;
      } else {
        log_f[i] = density.log_at(xi);
      }
    }
  }

  if (!log_p) exp_in_place(log_f);
  return log_f;
}

// ---------------------------------------------------------------------------
// CDF and quantile helpers for the special-case dispatch.
//
// Each accepts R's lower.tail / log.p flags directly and forms the requested
// tail itself, handing the flags on to R::pgamma / R::pnorm / R::qgamma /
// R::qnorm where those give that tail accurately (R::qgamma is given the
// other tail of a log probability close to 0; see qmhn_sqrt_gamma_one).  The
// dispatcher in mhn_pmhn.cpp / mhn_qmhn.cpp is therefore free to pass user
// flags through without intermediate clipping, preserving precision in the
// tails.
// ---------------------------------------------------------------------------

double pmhn_sqrt_gamma_one(double q, double shape, double sqrt_beta,
                           bool lower_tail, bool log_p) {
  if (ISNAN(q)) return R_IsNA(q) ? NA_REAL : R_NaN;
  if (q <= 0.0) {
    // F = 0
    if (log_p) return lower_tail ? R_NegInf : 0.0;
    return lower_tail ? 0.0 : 1.0;
  }
  if (q == R_PosInf) {
    // F = 1
    if (log_p) return lower_tail ? 0.0 : R_NegInf;
    return lower_tail ? 1.0 : 0.0;
  }
  const double t = sqrt_beta * q;
  if (t >= kSqrtDblMin) {
    return R::pgamma(t * t, shape, /*scale=*/1.0,
                     lower_tail ? 1 : 0, log_p ? 1 : 0);
  }
  // beta q^2 < DBL_MIN: exp(-t^2) = 1 in double precision, and the series
  // P(a, y) = y^a e^(-y) sum_n y^n / Gamma(a + n + 1) reduces to its first
  // term, y^a / Gamma(a + 1).  Its logarithm is 2 a log t.
  const double log_t = (t >= DBL_MIN) ? std::log(t)
                                      : std::log(sqrt_beta) + std::log(q);
  const double log_P = 2.0 * shape * log_t - R::lgamma1p(shape);
  if (lower_tail) {
    if (log_p) return log_P;
    // exp(log P) would carry the rounding of log P, eps |log P| relative
    // (8e-14 at P = 1e-170); t^(2a) through pow carries about one unit.
    // With t below 1.5e-154, P is a normal double only for a < 1, where
    // 1/Gamma(a + 1) = exp(-lgamma1p(a)) lies in [1, 1.13] and is exact to a
    // unit or two.
    if (t >= DBL_MIN) {
      return std::pow(t, 2.0 * shape) * std::exp(-R::lgamma1p(shape));
    }
    return std::exp(log_P);
  }
  return log_p ? mhn::log1m_exp(log_P) : -std::expm1(log_P);
}

Rcpp::NumericVector pmhn_sqrt_gamma(const Rcpp::NumericVector& q,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p) {
  const double shape = alpha / 2.0;
  const double sqrt_beta = std::sqrt(beta);
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = pmhn_sqrt_gamma_one(q[i], shape, sqrt_beta, lower_tail, log_p);
  }
  return out;
}

namespace {

// The body of pmhn_truncated_normal_one, instantiated for an offset formed
// as written (kPlain: (q - mu) / sigma) and for the others (TruncNormal::zq,
// and the Exponential limit), so that the ordinary case carries no test for
// the others.
template <bool kPlain>
double pmhn_tn_impl(double q, const TruncNormal& tn, bool lower_tail,
                    bool log_p) {
  if (ISNAN(q)) return R_IsNA(q) ? NA_REAL : R_NaN;
  if (q <= 0.0) {
    if (log_p) return lower_tail ? R_NegInf : 0.0;
    return lower_tail ? 0.0 : 1.0;
  }
  if (q == R_PosInf) {
    if (log_p) return lower_tail ? 0.0 : R_NegInf;
    return lower_tail ? 1.0 : 0.0;
  }
  const double z_0 = tn.z0;
  const double d = q / tn.sigma;                          // = z_q - z_0

  // Near the truncation point, h = d (|z0| + d) <= 1/4, F = lambda(z0) d J(d)
  // is a product of positive factors.  The forms further down obtain the same
  // F as the difference of two tail probabilities of order one, and lose it
  // once it is small: they put pmhn(1e-13, 1, 1, -3) 2.7e-3 low and return 0
  // at q = 1e-16.  The upper tail on the log scale, log(1 - F), needs F just as
  // much.  The upper tail on the linear scale, 1 - F, needs no such care: in
  // this region it is at least 0.6, and the forms further down give it to
  // 2.5e-15 relative, at worst just below z0 = 3.
  if ((lower_tail || log_p) && d * (std::fabs(z_0) + d) <= kTnSmall) {
    // Here F is at most about 0.4, so the upper tail 1 - F loses nothing.
    const double J = tn_mean_kernel(d, z_0, nullptr);
    double F = R_NaN;
    if (tn.hazard() >= DBL_MIN && d >= DBL_MIN) F = tn.hazard() * d * J;
    if (F >= DBL_MIN) {
      if (lower_tail) return log_p ? std::log(F) : F;
      return log_p ? std::log1p(-F) : 1.0 - F;
    }
    // F, lambda(z0) or d below the normal range: assemble log F instead.
    const double log_d = (d >= DBL_MIN) ? std::log(d)
                                        : std::log(q) - std::log(tn.sigma);
    const double log_F = tn.log_hazard() + log_d + std::log(J);
    if (lower_tail) return log_p ? log_F : std::exp(log_F);
    return log_p ? mhn::log1m_exp(log_F) : 1.0 - std::exp(log_F);
  }

  const double z_q = kPlain ? (q - tn.mu) / tn.sigma : tn.zq(q);

  // Under a strong positive tilt, z0 < -kTnLargeZ0, z_q carries an absolute
  // rounding of about eps |z0|, and once d is no larger than that z_q rounds
  // onto z0 and the difference of the two lower-tail normals below is 0,
  // although F is not.  There the lower tail is taken from the ratio
  //
  //   Phi(z_q) / Phi(z0) = exp(r),
  //   r = d (|z0| - d / 2) + log(erfcx(-z_q/sqrt2) / erfcx(-z0/sqrt2)),
  //
  // since Phi(z) = exp(-z^2/2) erfcx(-z/sqrt2) / 2 for z < 0 and
  // z_q^2 - z0^2 = d (2 z0 + d).  d = q / sigma enters it exactly, and the
  // erfcx ratio, whose logarithm changes by only about 1/|z| per unit of z
  // out there, does not feel the rounding of z_q.  Then
  //
  //   log F = log Phi(z0) + r + log(1 - exp(-r)) - log Phibar(z0).
  //
  // For d <= |z0| / 2, z_q <= z0 / 2 and log Phi(z0) + r = log Phi(z_q) is at
  // least z0^2 / 8 in size, so the rounding of its terms, each at most about
  // z0^2 / 2, is a few eps of it.  Where log Phi(z0) is below the largest
  // negative double (|z0| beyond about 1.9e154) Phi(z0) is nothing next to
  // Phi(z_q), and the form below gives log Phi(z_q) - log Phibar(z0), with
  // the same few eps.
  if (!kPlain && lower_tail && z_0 < -kTnLargeZ0 && d <= -0.5 * z_0 &&
      tn.log_lower0() > R_NegInf) {
    const double r = d * (-z_0 - 0.5 * d) +
        std::log(erfcx_tail(-z_q * M_SQRT1_2,
                            erfcx_tail(-z_0 * M_SQRT1_2)));
    double log_lower = (r > 0.0)
        ? tn.log_lower0() + r + mhn::log1m_exp(-r) - tn.log_mass
        : R_NegInf;
    if (log_lower > 0.0) log_lower = 0.0;
    return log_p ? log_lower : std::exp(log_lower);
  }

  // Below the untruncated mean (z_q <= 0) the lower tail can be far below the
  // rounding of 1, and is computed directly: derived as 1 - P(X > q) it is
  // exactly 0 there (against a true 1.3e-17 at alpha = 1, gamma = 18, q = 3,
  // for instance).  Every other request is taken from the upper tail, which is
  // at least 1/2 when z_q <= 0 and so loses nothing there, and the lower tail
  // above the mean as its complement.
  if (z_q <= 0.0 && lower_tail) {
    // The lower tail as a difference of two lower-tail normals, both small
    // and ordered.
    const double lo_q = R::pnorm(z_q, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
    const double lo_0 = tn.log_lower0();
    const double gap = lo_0 - lo_q;                       // <= 0
    double log_lower = (gap < 0.0) ? lo_q + mhn::log1m_exp(gap) : R_NegInf;
    log_lower -= tn.log_mass;
    if (log_lower > 0.0) log_lower = 0.0;
    return log_p ? log_lower : std::exp(log_lower);
  }

  if (z_0 >= 3.0) {
    if (!kPlain && tn.exp_limit) {
      // z0 = +Inf: Exponential(rate -gamma), P(X > q) = exp(gamma q).
      const double log_upper = tn.gamma * q;
      if (!lower_tail) return log_p ? log_upper : std::exp(log_upper);
      return log_p ? mhn::log1m_exp(log_upper) : -std::expm1(log_upper);
    }
    // Both ordinates are far into the right tail, where log Phi_upper(z) is
    // about -z^2/2.  Subtracting the two logs then cancels: at
    // (1, 1e-4, -1000) each is -2.5e9 while the answer at q = 1e-3 is -1, so
    // the result would carry an absolute error of 2.5e9 * eps = 5.5e-7 --
    // thirty times the working tolerance, and growing as z^2.
    //
    // Writing Phi_upper(z) = exp(-z^2/2) erfcx(z/sqrt2) / 2 removes it:
    //
    //   P(X > q) = Phi_up(z_q) / Phi_up(z_0)
    //            = exp(-(z_q^2 - z_0^2)/2) erfcx(z_q/sqrt2) / erfcx(z_0/sqrt2)
    //
    // and z_q^2 - z_0^2 factors as d (2 z_0 + d) with d = z_q - z_0 = q/sigma,
    // which only ever adds.  The ratio of the two erfcx values lies in (0, 1],
    // since erfcx decreases, and is taken as one quotient, so nothing of the
    // size of z^2 is ever formed and subtracted.
    //
    // The lower tail reaches this point only outside the series region, where
    // h = d (z_0 + d) > 1/4, and so (z_q^2 - z_0^2)/2 >= h/2 and
    // P(X > q) <= exp(-1/8) = 0.88.  The lower tail is then at least 0.12, and
    // 1 minus the upper tail loses at most three bits of it.
    //
    // quad = (z_q^2 - z_0^2)/2 is formed as d (z_0 + d/2), which is
    // 0.5 d (2 z_0 + d) to the bit but does not overflow with 2 z_0 past
    // 9e307.  A NaN, which a failed ratio would be, is passed on rather than
    // read as one of the bounds.
    const double quad = d * (z_0 + 0.5 * d);
    const double ratio = erfcx_tail(z_q * M_SQRT1_2, tn.erfcx0);
    if (!lower_tail && log_p) {
      const double log_upper = std::log(ratio) - quad;
      return (log_upper > 0.0) ? 0.0 : log_upper;
    }
    double upper = std::exp(-quad) * ratio;
    if (upper > 1.0) upper = 1.0;
    if (!lower_tail) return upper;
    return log_p ? std::log1p(-upper) : 1.0 - upper;
  }

  // log P(X > q).
  double log_upper = R::pnorm(z_q, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/1)
                     - tn.log_mass;
  if (log_upper > 0.0) log_upper = 0.0;
  if (!lower_tail) return log_p ? log_upper : std::exp(log_upper);
  return log_p ? mhn::log1m_exp(log_upper) : -std::expm1(log_upper);
}

}  // namespace

double pmhn_truncated_normal_one(double q, const TruncNormal& tn,
                                 bool lower_tail, bool log_p) {
  return tn.plain_offset ? pmhn_tn_impl<true>(q, tn, lower_tail, log_p)
                         : pmhn_tn_impl<false>(q, tn, lower_tail, log_p);
}

Rcpp::NumericVector pmhn_truncated_normal(const Rcpp::NumericVector& q,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p) {
  TruncNormal tn;
  tn.set(beta, gamma);
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  if (tn.plain_offset) {
    for (R_xlen_t i = 0; i < n; ++i) {
      out[i] = pmhn_tn_impl<true>(q[i], tn, lower_tail, log_p);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      out[i] = pmhn_tn_impl<false>(q[i], tn, lower_tail, log_p);
    }
  }
  return out;
}

double qmhn_sqrt_gamma_one(double p, double shape, double sqrt_beta,
                           bool lower_tail, bool log_p) {
  if (R_IsNA(p))  return NA_REAL;
  if (R_IsNaN(p)) return R_NaN;
  if (!log_p && (p < 0.0 || p > 1.0)) return R_NaN;
  if (log_p && p > 0.0)                return R_NaN;
  // Hand the tail flags to qgamma rather than converting first.  Turning an
  // upper-tail p into 1 - p loses it below about 1e-16 -- the quantile would
  // be Inf at 1e-20, where qgamma itself still resolves 1e-30 -- and
  // exponentiating a log-scale p throws away its whole point.  One case is
  // converted the other way: a log probability above log(1/2), of either
  // tail, is passed as the log of the other tail, log(1 - exp(p)).  R's
  // qgamma loses the far tail when it is given as the log of the near one,
  // close to 0.  As a lower-tail log probability, at shape 5e-11 and
  // log p = log(1 - 1e-300) it returns 119.2, where log P(Y > 119.2) is
  // -147.8 rather than -690.8, and at shape 5 it is 0.4 percent off; as an
  // upper-tail one, at shape 250 and log p = -1e-300 it returns 5.889, where
  // log P(Y <= 5.889) is -696.6 rather than -690.8 (R 4.5.1).  Given the far
  // tail it is right.
  const bool flip = log_p && p > -M_LN2;
  const double y = flip
      ? R::qgamma(mhn::log1m_exp(p), shape, /*scale=*/1.0,
                  lower_tail ? 0 : 1, 1)
      : R::qgamma(p, shape, /*scale=*/1.0, lower_tail ? 1 : 0, log_p ? 1 : 0);
  if (y >= DBL_MIN) return std::sqrt(y) / sqrt_beta;
  // y = beta x^2 below DBL_MIN, where P(a, y) = y^a / Gamma(a + 1) exactly
  // (see pmhn_sqrt_gamma_one), so the quantile is explicit:
  //   log x = (log p + lgamma(a + 1)) / (2 a) - log sqrt(beta).
  // The square root of a subnormal y has lost most of its digits already, and
  // an exactly 0 y (reached for x below 1e-162) would give 0 for an x that is
  // an ordinary number.
  double log_p_low;
  if (lower_tail) log_p_low = log_p ? p : std::log(p);
  else            log_p_low = log_p ? mhn::log1m_exp(p) : std::log1p(-p);
  return std::exp((log_p_low + R::lgamma1p(shape)) / (2.0 * shape)
                  - std::log(sqrt_beta));
}

Rcpp::NumericVector qmhn_sqrt_gamma(const Rcpp::NumericVector& p,
                                    double alpha, double beta,
                                    bool lower_tail, bool log_p) {
  const double shape = alpha / 2.0;
  const double sqrt_beta = std::sqrt(beta);
  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = qmhn_sqrt_gamma_one(p[i], shape, sqrt_beta, lower_tail, log_p);
  }
  return out;
}

namespace {

// Three of the solvers of qmhn_truncated_normal_one, for the small lower
// tail, for Newton's method past z0 = kTnQnormMaxZ0 and for the Newton step
// that refines the normal quantile from z0 = 3 on, are kept in functions of
// their own, so that the erfcx and series evaluations they call are inlined
// into them.

// The small lower tail (see qmhn_truncated_normal_one): the root d of
// I(d) = T, for DBL_MIN <= T and h = T (|z0| + T) <= 1/4.
//
// Start from the reversion of I(d) = d - z0 d^2/2 + (z0^2 - 1) d^3/6 + ...,
//   d = T + z0 T^2/2 + (2 z0^2 + 1) T^3/6 + O(T^4),
// and refine by Halley's method on I(d) - T, whose derivatives
// I' = exp(-z0 d - d^2/2) and I'' = -(z0 + d) I' come with the series.
// The error after a step is about (d^2 (z0 + d)^2/12 + d^2/6) times the
// cube of the step, at most 0.05 step^3 here, so a relative step of 4e-6
// leaves less than 1e-17 and the loop can stop without another pass.
// Past z0 = 1e150, where z0^2 overflows, the start is written with
// h = z0 T <= 1/4 instead, as T (1 + h/2 + (2 h^2 + T^2)/6).
double tn_small_tail_offset(double T, double z0) {
  double d;
  if (z0 <= 1e150) {
    d = T * (1.0 + T * (0.5 * z0 + T * (2.0 * z0 * z0 + 1.0) / 6.0));
  } else {
    const double h = z0 * T;
    d = T * (1.0 + 0.5 * h + (2.0 * h * h + T * T) / 6.0);
  }
  for (int it = 0; it < 8; ++it) {
    double slope;
    const double J = tn_mean_kernel(d, z0, &slope);
    const double r = d * J - T;
    const double step = r / (slope + 0.5 * r * (z0 + d));
    d -= step;
    if (std::fabs(step) <= 4e-6 * d) break;
  }
  return d;
}

// For z0 > kTnQnormMaxZ0: the offset t = z - z0 of the standardised
// quantile z from the truncation point, given log_upper = log(1 - p), the
// log of the upper-tail probability.  Returns NaN where Newton's method does
// not converge.
//
// Solving for t directly, rather than recovering it as z - z0, avoids the
// mu + sigma z cancellation; assembling the answer as sigma * (z - z0) only
// moves it: at gamma = -1000, beta = 1e-4 both z and z0 are about 70710
// while their difference is 1.4e-7, so the subtraction keeps four digits of
// the answer.
//
// The defining equation, in the erfcx form that pmhn_truncated_normal_one
// uses for the same reason,
//
//     G(t) = -t (2 z0 + t)/2 + log erfcx((z0+t)/sqrt2)
//            - log erfcx(z0/sqrt2) - log(1 - p) = 0,
//
// contains no such subtraction, and Newton on it converges in a few steps.
// d/dt of the left side is -sqrt(2/pi)/erfcx((z0+t)/sqrt2), the (z0+t)
// terms cancelling, since d/du log erfcx(u) = 2u - 2/(sqrt(pi) erfcx(u)): it
// is minus the normal hazard lambda(z0 + t).  G is then concave and
// decreasing.  With L = -log(1 - p), the start is the root L / z0 of the
// linear part of G, t z0 = L, and where L exceeds 1e-3 z0^2 the root of its
// quadratic part, t (z0 + t/2) = L,
//   t = 2 L / (z0 + sqrt(z0^2 + 2 L)).
// At either start G is below 0: at the second it is the log of the erfcx
// ratio, below 0 since erfcx decreases, and at the first that log less
// L^2 / (2 z0^2).  So the start lies to the right of the root and the
// iteration moves monotonically onto it.  Up to L = 1e-3 z0^2 the first
// start exceeds the second by a relative L / (2 z0^2) or less, 5e-4, and
// takes as many Newton steps (but for L from about 0.9e-3 z0^2 with z0 just
// above 30), while it saves the square root.  Past a point between
// 1e-3 z0^2 and 2e-3 z0^2, depending on z0, it takes one step more, and
// where L is far above z0^2 it lies a factor of about sqrt(L / 2) / z0 too
// far out, and each Newton step would then only halve it.  Newton's error
// after a step is G''/(2 G') = lambda'/(2 lambda) times the square of the
// step, and for z = z0 + t >= 3 that factor is below 1/(2z) <= 1/(2t), so
// once a step is below 1e-9 t the relative error left is below 1e-18.
//
// The start is formed from L / 2, and G and its slope divided by 4.  Each
// scaling is by a power of 2, so it is exact and leaves every Newton step as
// it was, and nothing overflows up to L = DBL_MAX, where 2 L would, and so
// would t (z0 + t/2), about L, once rounded up.
double tn_upper_offset(const TruncNormal& tn, double log_upper) {
  const double z0 = tn.z0;
  const double L = -log_upper;
  const double hz = 0.5 * z0;
  // z0^2 overflows past 1.3e154 and then selects L / z0.
  double d = (L <= 1e-3 * (z0 * z0))
      ? L / z0
      : L / (hz + std::hypot(hz, std::sqrt(0.5 * L)));
  if (!(d > 0.0)) d = 0.0;
  for (int it = 0; it < 60; ++it) {
    const double ea = erfcx_tail((z0 + d) * M_SQRT1_2);
    // Past z0 = kTnLargeZ0 the two logarithms of erfcx, each about
    // -log(z0), are taken as the log of their quotient, which loses nothing
    // of their difference d / z0, and z0 + d/2 replaces 2 z0 + d, which would
    // overflow past 9e307.  G / 4 performs the operations of G in the same
    // order, each operand scaled by a power of 2, and so is a quarter of G
    // to the last bit.
    double G4;
    if (z0 <= kTnLargeZ0) {
      G4 = -(0.5 * d) * (z0 + 0.5 * d) * 0.5 + 0.25 * std::log(ea)
           - 0.25 * tn.log_erfcx0() - 0.25 * log_upper;
    } else {
      G4 = -(0.5 * d) * (hz + 0.25 * d) + 0.25 * std::log(ea / tn.erfcx0)
           - 0.25 * log_upper;
    }
    const double dG4 = -(0.5 * M_SQRT1_2) / (std::sqrt(M_PI) * ea);
    const double step = G4 / dG4;
    if (!std::isfinite(step)) break;
    d -= step;
    if (!(d > 0.0)) d = 0.0;
    if (std::fabs(step) <= 1e-9 * d) return d;
  }
  return R_NaN;
}

// For 3 <= z0 <= kTnQnormMaxZ0: one Newton step on the offset t = zq - z0
// that the normal quantile zq gave, with log_upper = log(1 - p).
//
// From z0 = 3 on t = zq - z0 subtracts two quantities of size z0 or more, and
// the error of zq, a few eps of it, is a relative error of about eps z0^2 / h
// in t (h = t (z0 + t)).  One Newton step on the erfcx form of the defining
// equation (see tn_upper_offset),
//   G(t) = -t (z0 + t/2) + log(erfcx((z0 + t)/sqrt2) / erfcx(z0/sqrt2))
//          - log(1 - p),
// leaves a relative error of about the square of that one (Newton's factor
// G''/(2 G') is below 1/(2t) here).  What remains is the rounding of G, a
// few eps times the larger of 1 and h, over its slope lambda(z0 + t) >
// z0 + t: a relative few eps (1 + h) / h in t, a few eps for h of order one
// or more and about 100 eps at the smallest h this branch meets, about 0.01
// near z0 = 3 (below it the small-tail series is used).  G is the logarithm
// of the upper tail that pmhn forms from z0 = 3 on, less log(1 - p), so the
// quantile is also the root of pmhn's own function.
double tn_refine_offset(const TruncNormal& tn, double t, double log_upper) {
  const double z0 = tn.z0;
  const double ratio = erfcx_tail((z0 + t) * M_SQRT1_2, tn.erfcx0);
  const double G = -(t * (z0 + 0.5 * t)) + std::log(ratio) - log_upper;
  const double dG = -M_SQRT2 / (kSqrtPi * (ratio * tn.erfcx0));
  const double step = G / dG;
  if (std::isfinite(step) && t - step > 0.0) return t - step;
  return t;
}

}  // namespace

// Work with the standardised lower bound z0 = -mu/sigma and with
// t = x/sigma >= 0, so that the truncation identity reads
//
//     Phibar(t + z0) = Phibar(z0) * (1 - p),
//
// and hence  t + z0 = qnorm(log Phibar(z0) + log1p(-p), upper, log).
//
// The truncation mass is taken from pnorm's own upper tail on the log scale,
// which keeps it representable however far out z0 goes (1 - Phi(z0) is
// exactly 0 once z0 passes about 8.3), and the answer is assembled as
// sigma * t rather than mu + sigma * (t + z0): those are equal, since
// sigma * z0 = -mu, but the second subtracts two nearly equal quantities once
// mu is far below zero.
double qmhn_truncated_normal_one(double p, const TruncNormal& tn,
                                 bool lower_tail, bool log_p) {
  if (std::isnan(p)) return R_IsNA(p) ? NA_REAL : R_NaN;
  if (!log_p && (p < 0.0 || p > 1.0)) return R_NaN;
  if (log_p && p > 0.0)                return R_NaN;
  // p_low = 0 gives 0 and p_low = 1 gives Inf.
  if (log_p) {
    if (p == R_NegInf) return lower_tail ? 0.0 : R_PosInf;
    if (p == 0.0)      return lower_tail ? R_PosInf : 0.0;
  } else {
    if (p == 0.0) return lower_tail ? 0.0 : R_PosInf;
    if (p == 1.0) return lower_tail ? R_PosInf : 0.0;
  }

  // Each of these is formed only where it is needed, and from the tail the
  // caller supplied, so that a small probability is never rebuilt as 1 minus
  // one near 1.
  auto log_p_low = [&]() -> double {
    if (lower_tail) return log_p ? p : std::log(p);
    return log_p ? mhn::log1m_exp(p) : std::log1p(-p);
  };
  auto log_p_up = [&]() -> double {
    if (!lower_tail) return log_p ? p : std::log(p);
    return log_p ? mhn::log1m_exp(p) : std::log1p(-p);
  };

  // Exponential(rate -gamma): the quantile of an upper-tail probability P is
  // log(P) / gamma.
  if (tn.exp_limit) return log_p_up() / tn.gamma;

  const double z0 = tn.z0;
  const double sigma = tn.sigma;

  // The small lower tail.  F = lambda(z0) I(d) (see mhn_special_cases.h), so
  // the quantile solves I(d) = T with T = p_low / lambda(z0) and no
  // subtraction anywhere.  The forms below instead recover d as z - z0 from a
  // z computed through a tail probability of order one, which loses about
  // eps (1 + |z0|) / d of it relative: they put qmhn(1e-12, 1, 1, -3) 6e-4
  // off and return 0 at p = 1e-16.  The direct solution costs a few series
  // evaluations against one qnorm.  It is used within h = d (|z0| + d) <= 1/4,
  // where the series needs about twenty terms, and there where that loss would
  // exceed about 1000 eps (2e-13), d < (1 + |z0|) / 1000.  For z0 below about
  // -15 the bound h <= 1/4 ends the region first, and between the two the
  // qnorm form loses up to about 4 eps z0^2: a few times the eps z0^2 by which
  // the quantile itself moves when gamma changes in its last place.  That is
  // below 1e-12 up to z0 = kTnQnormMaxZ0 = 30.  Past it the alternative is an
  // iteration on erfcx that costs more than this one, so there this is used
  // wherever h <= 1/4.
  const double hazard = log_p ? 0.0 : tn.hazard();
  const double T = (hazard >= DBL_MIN)
      ? (lower_tail ? p : 1.0 - p) / hazard
      : std::exp(log_p_low() - tn.log_hazard());
  if (T * (std::fabs(z0) + T) <= kTnSmall &&
      (z0 > kTnQnormMaxZ0 || T <= 1e-3 * (1.0 + std::fabs(z0)))) {
    // Below DBL_MIN, d^2 is nothing next to d, and I(d) = (1 - exp(-h)) / z0
    // with h = z0 d, so d = -log(1 - z0 T) / z0, which is T to double
    // precision unless z0 T passes 1e-17, from z0 of about 1e290.
    if (T < DBL_MIN) {
      const double log_T = log_p_low() - tn.log_hazard();
      const double hT = std::exp(log_T + std::log(z0));
      if (!(z0 > 0.0) || !(hT > 1e-17)) {
        return std::exp(log_p_low() - tn.log_hazard() + std::log(sigma));
      }
      return std::exp(std::log(-std::log1p(-hT)) - std::log(z0) +
                      std::log(sigma));
    }
    return sigma * tn_small_tail_offset(T, z0);
  }

  double t;
  if (z0 > kTnQnormMaxZ0) {
    t = tn_upper_offset(tn, log_p_up());
    if (std::isnan(t)) {
      mhn::warn("qmhn: the truncated-normal quantile did not converge "
                "for beta=%g, gamma=%g; returning NA", tn.beta, tn.gamma);
      return NA_REAL;
    }
  } else {
    // What the identity needs is Phibar(z) = Phibar(z0) (1 - p_low), the
    // upper-tail probability 1 - p_low as supplied when the caller gave the
    // upper tail.  For a probability on the linear scale, where that product
    // is at most 1/2 and a normal double, qnorm takes it as it stands: it is
    // formed to within 2 eps relative, as accurately as its logarithm below,
    // and qnorm on the linear scale costs about half as much.
    double zq;                                  // the standard normal quantile
    bool linear = false;
    if (!log_p) {
      const double upper_z = tn.mass * (lower_tail ? 1.0 - p : p);
      if (upper_z <= 0.5 && upper_z >= DBL_MIN) {
        zq = R::qnorm(upper_z, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
        linear = true;
      }
    }
    if (!linear) {
      const double log_upper_z = tn.log_mass + log_p_up();   // log Phibar(z)
      if (-log_upper_z >= DBL_MIN) {
        zq = R::qnorm(log_upper_z, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/1);
      } else {
        // Phi(z) = -log Phibar(z) is below DBL_MIN, which needs z0 below
        // about -37, and qnorm could no longer recover it from
        // log Phibar(z).  Take the same point from the lower tail,
        //   Phi(z) = Phi(z0) + p_low Phibar(z0),
        // formed on the log scale.
        const double a = tn.log_lower0();
        const double b = log_p_low() + tn.log_mass;
        const double log_lower_z =
            std::max(a, b) + std::log1p(std::exp(-std::fabs(a - b)));
        zq = R::qnorm(log_lower_z, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
      }
    }
    // The quantile x = mu + sigma zq, taken below as sigma (zq - z0).  Under
    // a strong positive tilt (fma_offset, z0 below -kTnLargeZ0) mu is far
    // above 0 and sigma zq a correction to it; zq - z0 and its product with
    // sigma then each round at the size of mu, and x comes out a unit or two
    // in its last place off, rather than the double nearest it, which is the
    // only answer once sigma is below the spacing of the doubles at mu (|z0|
    // beyond about 1e16).  There mu is taken as the pair
    // gamma / (2 beta) + r / beta, with the remainder
    // r = gamma / 2 - beta mu exact (fma; the halving of gamma is exact), and
    // x as one rounding of mu plus the rest.  Neither needs 2 beta, which
    // overflows for beta above about 9e307.
    if (tn.fma_offset && std::isfinite(tn.mu)) {
      const double mu_lo = mhn::fma(-tn.mu, tn.beta, 0.5 * tn.gamma) / tn.beta;
      return std::max(tn.mu + (mu_lo + sigma * zq), 0.0);
    }
    t = zq - z0;
    if (z0 >= 3.0 && t > 0.0) t = tn_refine_offset(tn, t, log_p_up());
  }
  return std::max(sigma * t, 0.0);
}

Rcpp::NumericVector qmhn_truncated_normal(const Rcpp::NumericVector& p,
                                          double beta, double gamma,
                                          bool lower_tail, bool log_p) {
  TruncNormal tn;
  tn.set(beta, gamma);
  const R_xlen_t n = p.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = qmhn_truncated_normal_one(p[i], tn, lower_tail, log_p);
  }
  return out;
}

}  // namespace mhn

// ---------------------------------------------------------------------------
// R-visible exports
// ---------------------------------------------------------------------------

// [[Rcpp::export(.is_sqrt_gamma, rng = false)]]
bool is_sqrt_gamma_R(double gamma, double beta = 1.0) {
  return mhn::is_sqrt_gamma(gamma, beta);
}

// [[Rcpp::export(.is_truncated_normal, rng = false)]]
bool is_truncated_normal_R(double alpha) {
  return mhn::is_truncated_normal(alpha);
}
