# Tests for .mhn_log_normalizing_const() and its internal methods.
#
# The normalizing constant is Psi[alpha/2, gamma/sqrt(beta)].
# Helper: nc(alpha, beta, gamma) = exp(.mhn_log_normalizing_const(...)).

nc <- function(alpha, beta, gamma, ...) {
  exp(mhn:::.mhn_log_normalizing_const(alpha, beta, gamma, ...))
}

# One unit in the last place of a double of the size of x.
ulp_of <- function(x) 2^(floor(log2(abs(x))) - 52)

# ============================================================
# 1. gamma = 0 special case: Psi[alpha/2, 0] = Gamma(alpha/2)
# ============================================================

test_that("gamma = 0: Psi[alpha/2, 0] = Gamma(alpha/2)", {
  alpha_vals <- c(0.25, 0.5, 1, 1.5, 2, 3, 5, 10)
  for (alpha in alpha_vals) {
    expect_equal(nc(alpha, 1, 0), gamma(alpha / 2),
                 tolerance = 1e-12,
                 label = sprintf("alpha=%.2f", alpha))
  }
})

test_that("gamma = 0 result is independent of beta", {
  for (beta in c(0.01, 0.1, 1, 10, 100)) {
    expect_equal(nc(2, beta, 0), gamma(1),
                 tolerance = 1e-12,
                 label = sprintf("beta=%.2f", beta))
  }
})

# ============================================================
# 2. alpha = 1 closed form (Sun et al. 2023, Supplementary Lemma 9(c))
#    Psi[1/2, z] = 2 sqrt(pi) exp(z^2 / 4) Phi(z / sqrt(2))
# ============================================================

test_that("alpha = 1: matches Sun et al. (2023, Supplementary Lemma 9(c))", {
  closed_form <- function(gamma, beta) {
    z <- gamma / sqrt(beta)
    2 * sqrt(pi) * exp(z^2 / 4) * pnorm(z / sqrt(2))
  }

  test_cases <- expand.grid(
    gamma = c(-10, -5, -1, 0, 1, 5, 10),
    beta = c(0.5, 1, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    g <- test_cases$gamma[i]
    b <- test_cases$beta[i]
    expect_equal(nc(1, b, g), closed_form(g, b),
                 tolerance = 1e-12,
                 label = sprintf("beta=%.1f, gamma=%.1f", b, g))
  }
})

# ============================================================
# 3. alpha = 2, gamma >= 0 closed form (Sun et al. 2023, Supplementary
#    Lemma 9(c) ratio)
# ============================================================

test_that("alpha = 2, gamma >= 0: matches the Supplementary Lemma 9(c) ratio", {
  closed_form <- function(gamma, beta) {
    z <- gamma / sqrt(beta)
    psi_half <- 2 * sqrt(pi) * exp(z^2 / 4) * pnorm(z / sqrt(2))
    ratio <- z / 2 + exp(-z^2 / 4) /
      (2 * sqrt(pi) * pnorm(z / sqrt(2)))
    psi_half * ratio
  }

  gamma_vals <- c(0.5, 1, 2, 5, 10)
  for (g in gamma_vals) {
    expect_equal(nc(2, 1, g), closed_form(g, 1),
                 tolerance = 1e-12,
                 label = sprintf("gamma=%.1f", g))
  }
})

# ============================================================
# 4. Positivity: Psi > 0 across the full parameter grid
# ============================================================

test_that("normalizing constant is positive for all valid inputs", {
  test_cases <- expand.grid(
    alpha = c(0.1, 0.5, 1, 2, 3, 5, 10),
    beta = c(0.1, 1, 10),
    gamma = c(-15, -10, -5, -1, 0, 1, 5, 10, 15)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    val <- nc(a, b, g)
    expect_true(is.finite(val) && val > 0,
                label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f",
                                a, b, g))
  }
})

# ============================================================
# 5. Monotonicity: Psi increases with gamma
#    (differentiating the series Sum Gamma(a + n/2) z^n / n! term by term
#     gives d Psi[a, z] / dz = Psi[a + 1/2, z], which is positive)
# ============================================================

test_that("normalizing constant increases with gamma", {
  alpha_vals <- c(0.5, 1, 2, 5)
  gamma_seq <- c(-10, -5, -1, 0, 1, 5, 10)
  for (a in alpha_vals) {
    vals <- sapply(gamma_seq, function(g) nc(a, 1, g))
    for (j in seq_along(vals)[-1]) {
      expect_true(vals[j] > vals[j - 1],
                  label = sprintf("alpha=%.1f, gamma=%.1f > %.1f",
                                  a, gamma_seq[j], gamma_seq[j - 1]))
    }
  }
})

# ============================================================
# 6. Recurrence identity (Sun et al. 2023, Supplementary Lemma 9(a)):
#    gamma >= 0 only
#    Psi[(alpha + 2)/2, z] = (alpha/2) Psi[alpha/2, z]
#                          + (z/2)     Psi[(alpha + 1)/2, z]
#
#    For gamma < 0 the recurrence RHS becomes a difference of nearly
#    equal positive numbers; cancellation amplifies numerical-integration
#    error, so we restrict the check to gamma >= 0.
# ============================================================

test_that("recurrence identity holds for gamma >= 0", {
  test_cases <- expand.grid(
    alpha = c(1, 2, 3, 5),
    gamma = c(0, 1, 5)
  )
  beta <- 1
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    g <- test_cases$gamma[i]
    z <- g / sqrt(beta)
    lhs <- nc(a + 2, beta, g)
    rhs <- (a / 2) * nc(a, beta, g) +
      (z / 2) * nc(a + 1, beta, g)
    expect_equal(lhs, rhs,
                 tolerance = 1e-12,
                 label = sprintf("alpha=%d, gamma=%.1f", a, g))
  }
})

# ============================================================
# 7. Accuracy of the numerical-integration path for gamma < 0
#    Uses the alpha = 1 closed form as ground truth (exact for all gamma).
# ============================================================

test_that("integration matches alpha=1 closed form for gamma < 0", {
  closed_form <- function(gamma, beta) {
    z <- gamma / sqrt(beta)
    2 * sqrt(pi) * exp(z^2 / 4) * pnorm(z / sqrt(2))
  }

  # alpha = 1 normally dispatches to closed form; here we call the
  # integration kernel directly to exercise its accuracy.
  gamma_vals <- c(-1, -5, -10, -30)
  for (g in gamma_vals) {
    log_int <- mhn:::.psi_integrate(1, 1, g)
    expect_equal(exp(log_int), closed_form(g, 1),
                 tolerance = 1e-11,
                 label = sprintf("gamma=%.1f", g))
  }
})

# ============================================================
# 8. Scaling consistency in beta
#    Different (beta, gamma) pairs that share z = gamma/sqrt(beta) must
#    yield the same Psi value.
# ============================================================

test_that("consistent across beta scalings for same z", {
  # Psi[alpha/2, z] depends only on z; .mhn_log_normalizing_const
  # accepts the full (alpha, beta, gamma) triple, so we sweep different
  # (beta, gamma) pairs that share the same z = gamma / sqrt(beta).  Forming
  # gamma = z sqrt(beta) and dividing again recovers z to an ulp or two,
  # which moves log Psi by about E[U] * 2 eps |z| <= eps z^2, 2e-13 at z = 30.
  for (alpha in c(0.5, 3, 7)) {
    for (z in c(-20, -2, 3, 30)) {
      ref <- mhn:::.mhn_log_normalizing_const(alpha, 1, z)
      for (b in c(1e-6, 0.3, 7, 1e6)) {
        got <- mhn:::.mhn_log_normalizing_const(alpha, b, z * sqrt(b))
        expect_lt(abs(got - ref), 1e-11,
                  label = sprintf("alpha=%g, z=%g, beta=%g", alpha, z, b))
      }
    }
  }
})

# ============================================================
# 9. Direct test of the series path for gamma > 0
#    Cross-checks .psi_series against the alpha = 1 closed form.
# ============================================================

test_that("series method matches closed form for gamma > 0", {
  closed_form <- function(gamma, beta) {
    z <- gamma / sqrt(beta)
    2 * sqrt(pi) * exp(z^2 / 4) * pnorm(z / sqrt(2))
  }

  gamma_vals <- c(0.5, 1, 5, 10)
  for (g in gamma_vals) {
    log_ser <- mhn:::.psi_series(1, 1, g, .Machine$double.eps^0.5)
    expect_equal(exp(log_ser), closed_form(g, 1),
                 tolerance = 1e-12,
                 label = sprintf("gamma=%.1f", g))
  }
})

# ============================================================
# 10. log Psi where Psi itself overflows
#    The closed forms above, taken in logs: log Gamma(alpha/2) at z = 0,
#    log(2 sqrt(pi)) + z^2/4 + log Phi(z/sqrt(2)) at alpha = 1, and at
#    alpha = 2 log Psi[1, z] = log(1 + z sqrt(pi) e^(z^2/4) Phi(z/sqrt(2))).
# ============================================================

test_that("log Psi matches the closed forms where Psi overflows", {
  log_nc <- function(alpha, z) mhn:::.mhn_log_normalizing_const(alpha, 1, z)
  for (a in c(400, 1e4, 1e6)) {
    ref <- lgamma(a / 2)
    expect_lt(abs(log_nc(a, 0) - ref) / ulp_of(ref), 8,
              label = sprintf("ulp error, alpha=%g, z=0", a))
  }
  for (z in c(60, 100, 1000, 1e5)) {
    ref1 <- log(2 * sqrt(pi)) + z^2 / 4 + pnorm(z / sqrt(2), log.p = TRUE)
    expect_lt(abs(log_nc(1, z) - ref1) / ulp_of(ref1), 8,
              label = sprintf("ulp error, alpha=1, z=%g", z))
    lead <- log(z * sqrt(pi)) + pnorm(z / sqrt(2), log.p = TRUE)
    ref2 <- z^2 / 4 + lead + log1p(exp(-z^2 / 4 - lead))
    expect_lt(abs(log_nc(2, z) - ref2) / ulp_of(ref2), 8,
              label = sprintf("ulp error, alpha=2, z=%g", z))
  }
})

# ============================================================
# 11. Edge cases
# ============================================================

test_that("very small alpha works", {
  val <- nc(0.01, 1, 0)
  expect_equal(val, gamma(0.005), tolerance = 1e-12)
})

test_that("very large alpha works", {
  val <- nc(100, 1, 0)
  expect_equal(val, gamma(50), tolerance = 1e-12)
})

# ============================================================
# 12. The gamma > 0 series against numerical integration on a grid
#
# The dispatcher takes these points from the Lemma 10 series, truncated at
# lengths that follow from the critical indices of lemma10_C
# (mhn_psi_series.cpp).  Their discriminant carries alpha_adj * z^2, as in the
# proof in the Supplementary Material; the statement of Lemma 10 prints
# alpha_adj * z (see Appendix B.1 of the package's theory vignette).
# stats::integrate at rel.tol = 1e-12 is the
# independent check.
# ============================================================

test_that("gamma > 0 series agrees with numerical integration on a grid", {
  psi_via_integration <- function(alpha, beta, gamma) {
    integrand <- function(x) x^(alpha - 1) * exp(-beta * x^2 + gamma * x)
    val <- stats::integrate(integrand, lower = 0, upper = Inf,
                            rel.tol = 1e-12)$value
    2 * beta^(alpha / 2) * val
  }

  grid <- expand.grid(
    alpha = c(0.5, 1.5, 3, 5, 8),
    z     = c(0.1, 0.5, 1, 2, 5, 10)
  )
  for (i in seq_len(nrow(grid))) {
    a <- grid$alpha[i]
    z <- grid$z[i]
    beta  <- 1
    gamma <- z * sqrt(beta)
    expect_equal(
      nc(a, beta, gamma),
      psi_via_integration(a, beta, gamma),
      tolerance = 1e-11,
      label = sprintf("alpha=%.2f, z=%.2f", a, z)
    )
  }
})

# ============================================================
# Reference values for the sections below
#
# Unless a test says otherwise, every log Psi reference was computed with
# Rmpfr at 160 bits or more, by one of three routes:
#   the power series sum_n Gamma(alpha/2 + n/2) z^n / n!, for |z| <= 40,
#     with enough extra bits to absorb the cancellation for z < 0;
#   the expansion 2 |z|^-alpha sum_k (-1)^k Gamma(alpha + 2k) / (k! z^(2k)),
#     for z < -40, summed to its smallest term (below 1e-45);
#   2 int exp(alpha v - e^(2v) + z e^v) dv (x = e^v) by the trapezoid rule,
#     which converges geometrically for this entire integrand at either sign
#     of z, for z > 40 and wherever a test says so; its value is unchanged
#     when the step is halved.
# The routes agree with one another at overlapping points to 1e-45 or
# better.  A few tests at extreme parameters compare with limits formed in
# double precision instead, and say so.
# The accuracy of log Psi is stated as an absolute error: it is the relative
# error of Psi, and hence of every density and probability built on it.
# ============================================================

log_psi_abs_err <- function(alpha, beta, gamma, ref) {
  abs(mhn:::.mhn_log_normalizing_const(alpha, beta, gamma) - ref)
}

# Relative difference, stated explicitly: expect_equal() switches to an
# absolute difference once the expected value is below the tolerance, and a
# variance of 1e-12 is.
rel_diff <- function(x, y) abs(x - y) / abs(y)

# ============================================================
# 13. Shapes next to alpha = 1 and alpha = 2
#
# The closed forms of Sun et al. (2023, Supplementary Lemma 9(c)) hold at
# alpha = 1 and alpha = 2 exactly.  Since d log Psi[alpha/2, z] / d alpha =
# E[log U] with U = sqrt(beta) X, a shape delta away from either point
# differs in log Psi by about delta E[log U], which is -delta (log|z| + 0.577)
# under a strong negative tilt.  The moment formulas divide Psi at alpha + 1
# by Psi at alpha, so any such offset in one of them goes straight into the
# mean, and Lemma 2c amplifies it in the variance.
# ============================================================

test_that("log Psi next to alpha = 1 and 2 matches its reference", {
  # Shapes 1.4e-8 from the special points, inside sqrt(eps) = 1.49e-8.  The
  # offset of a closed form taken at the special point instead would be
  # 2.0e-7 at z = -1e6 and 5.5e-8 at z = 100.  The general paths are within
  # 4e-14 at z = -1e6 and one unit in the last place (4.5e-13) at z = 100,
  # so 1e-11 leaves a wide margin and still catches either offset.
  refs <- list(
    list(alpha = 1 + 1.4e-8, gamma = -1e6, ref = -13.12236357890449645),
    list(alpha = 1 - 1.4e-8, gamma = -1e6, ref = -13.12236317590816082),
    list(alpha = 2 + 1.4e-8, gamma = 100,  ref = 2505.177535183682513),
    list(alpha = 2 - 1.4e-8, gamma = 100,  ref = 2505.177535074143069)
  )
  for (r in refs) {
    expect_lt(log_psi_abs_err(r$alpha, 1, r$gamma, r$ref), 1e-11,
              label = sprintf("alpha=%.9f, gamma=%g", r$alpha, r$gamma))
  }
})

test_that("log Psi is continuous through alpha = 1 and alpha = 2", {
  # The closed form at the special point and the general path beside it must
  # meet.  A shape change of 2^-40 (9.1e-13) moves log Psi by at most
  # 9.1e-13 |E[log U]|, below 2e-11 for these tilts; the tolerance allows that
  # plus a few units in the last place of log Psi.
  h <- 2^-40
  for (gamma in c(-1e6, -1e3, -5, 5, 15, 100, 1e3)) {
    lp <- function(a) mhn:::.mhn_log_normalizing_const(a, 1, gamma)
    tol <- 3e-11 + 8 * ulp_of(lp(1))
    expect_lt(abs(lp(1 + h) - lp(1)), tol, label = sprintf("alpha=1+, gamma=%g", gamma))
    expect_lt(abs(lp(1 - h) - lp(1)), tol, label = sprintf("alpha=1-, gamma=%g", gamma))
    if (gamma > 0) {
      tol2 <- 3e-11 + 8 * ulp_of(lp(2))
      expect_lt(abs(lp(2 + h) - lp(2)), tol2, label = sprintf("alpha=2+, gamma=%g", gamma))
      expect_lt(abs(lp(2 - h) - lp(2)), tol2, label = sprintf("alpha=2-, gamma=%g", gamma))
    }
  }
})

test_that("mean and variance next to alpha = 1 agree with their neighbours", {
  # 1 +- 1.4e-8 lies inside sqrt(eps) of 1 and 1 +- 1.6e-8 just outside it.
  # The true values differ by 2e-9 of relative shape, i.e. about 2e-9
  # relative in the mean and the variance (E[X] ~ alpha / |gamma| and
  # Var ~ alpha / gamma^2 here).  An offset of 2e-7 in one log Psi is a
  # relative error of 1e-7 to 2e-7 in the mean, and in the variance it gives
  # 1.0e-7 in place of 1.0e-12 at gamma = -1e6 and a 5 percent error at
  # gamma = -1e3.  The bounds (1e-8 for the mean, 1e-6 for the variance,
  # which takes the quadrature at gamma = -1e6) sit well above the true
  # difference and the accuracy of either value, and well below those
  # errors.
  for (gamma in c(-1e6, -1e3)) {
    for (s in c(1, -1)) {
      inside  <- 1 + s * 1.4e-8
      outside <- 1 + s * 1.6e-8
      expect_lt(rel_diff(mhn_var(inside, 1, gamma), mhn_var(outside, 1, gamma)),
                1e-6,
                label = sprintf("var, alpha-1=%g, gamma=%g", inside - 1, gamma))
      expect_lt(rel_diff(mhn_mean(inside, 1, gamma), mhn_mean(outside, 1, gamma)),
                1e-8,
                label = sprintf("mean, alpha-1=%g, gamma=%g", inside - 1, gamma))
    }
  }
})

test_that("mean and variance next to alpha = 2 agree with their neighbours", {
  # As above, for the alpha = 2 closed form (gamma >= 0).  The mean moves by
  # d E[X] / d alpha * 2e-9 = Cov(X, log X) * 2e-9, about 4e-13 relative at
  # gamma = 100, and the two log Psi values each carry about one unit in the
  # last place (4.5e-13).  An offset of 5.5e-8 in log Psi would show up as a
  # relative error of 5.5e-8 in the mean.  Lemma 2c, alpha / (2 beta) +
  # E[X] (gamma / (2 beta) - E[X]), leaves 0.5 out of terms of 2500, so the
  # variance carries the relative error of the mean times |gamma / (2 beta) -
  # 2 E[X]| E[X] / Var = 5.0e3: about 5e-9 for an accurate mean, 2.7e-4 for
  # the offset.
  for (s in c(1, -1)) {
    inside  <- 2 + s * 1.4e-8
    outside <- 2 + s * 1.6e-8
    expect_lt(rel_diff(mhn_mean(inside, 1, 100), mhn_mean(outside, 1, 100)),
              1e-10, label = sprintf("mean, alpha-2=%g", inside - 2))
    expect_lt(rel_diff(mhn_var(inside, 1, 100), mhn_var(outside, 1, 100)),
              1e-6, label = sprintf("var, alpha-2=%g", inside - 2))
  }
})

test_that("the density next to alpha = 2 integrates to one", {
  # dmhn combines log Psi with the kernel at the same shape, so the two must
  # belong together.  At (alpha, beta, gamma) = (2 +- 1.4e-8, 1, 100) the law
  # is close to a normal with mean 50 and variance 1/2, so [40, 60] holds all
  # but about erfc(10) = 2e-45 of the mass.  stats::integrate reaches 1e-13 on
  # this smooth integrand; 1e-10 is well above that and well below an offset
  # of 5.5e-8.
  for (alpha in c(2 + 1.4e-8, 2 - 1.4e-8)) {
    total <- stats::integrate(function(x) dmhn(x, alpha, 1, 100), 40, 60,
                              rel.tol = 1e-13)$value
    expect_equal(total, 1, tolerance = 1e-10,
                 label = sprintf("alpha-2=%g", alpha - 2))
  }
})

test_that("log Psi next to gamma = 0 carries the tilt", {
  # Psi[alpha/2, 0] = Gamma(alpha/2) holds at z = 0 exactly.  Since
  # d log Psi[alpha/2, z] / dz = E U and its second derivative is Var U, with
  # U^2 ~ Gamma(alpha/2, 1) at z = 0,
  #   log Psi[alpha/2, z] = lgamma(alpha/2) + z E U + (z^2 / 2) Var U,
  # up to a term of order z^3, 1e-28 here, with E U = Gamma((alpha+1)/2) /
  # Gamma(alpha/2) and Var U = alpha/2 - (E U)^2.  Taking Gamma(alpha/2) for
  # a tilt of 1e-9 would be 1e-9 off; the series (z > 0) and the quadrature
  # (z < 0) are within a few units in the last place of log Psi, and eight of
  # them, with 1e-15 for a log Psi near 0, allow for that.  beta = 1e6 checks
  # that the tilt is taken as z = gamma / sqrt(beta).
  for (a in c(0.3, 2.5, 50)) {
    eu <- exp(lgamma((a + 1) / 2) - lgamma(a / 2))
    vu <- a / 2 - eu^2
    ref0 <- lgamma(a / 2)
    for (z in c(-1e-9, 1e-9)) {
      for (b in c(1, 1e6)) {
        info <- sprintf("alpha=%g beta=%g z=%g", a, b, z)
        got <- mhn:::.mhn_log_normalizing_const(a, b, z * sqrt(b))
        expect_lt(abs(got - (ref0 + z * eu + z^2 / 2 * vu)),
                  1e-15 + 8 * ulp_of(ref0), label = info)
      }
    }
  }
})

# ============================================================
# 14. Positive tilt beyond z = 20
#
# For z = gamma / sqrt(beta) > 20, and for alpha > 100, the normalizing
# constant comes from the quadrature rather than from the series of
# Sun et al. (2023, Supplementary Lemma 10): the number of series terms
# grows like z^2, and the rounding of their running sum grows with that
# count, reaching a few hundred units in the last place of log Psi near
# z = 1e3.
# ============================================================

test_that("log Psi at a large positive tilt is accurate to a few ulp", {
  # The quadrature is within about one unit in the last place of log Psi at
  # these points; 8 ulp allows for rounding that differs between platforms.
  # A summation error of 81 to 248 ulp (the series at these points) fails.
  refs <- list(
    list(alpha = 1.5, z = 999,  ref = 249504.6223156720277),
    list(alpha = 0.5, z = 1000, ref = 249998.1582088242765),
    list(alpha = 300, z = 999,  ref = 251359.4734112481109),
    list(alpha = 2.5, z = 1170, ref = 342235.8229304422176)
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
  # The same z reached with beta != 1 (499.5 / sqrt(0.25) = 999 exactly): Psi
  # depends on (beta, gamma) only through z, and so does the route taken.
  expect_lt(log_psi_abs_err(1.5, 0.25, 499.5, 249504.6223156720277) /
              ulp_of(249504.6223156720277), 8)
})

test_that("the mean at a large positive tilt keeps its accuracy", {
  # E[X] = Psi[(alpha+1)/2, z] / (sqrt(beta) Psi[alpha/2, z]) (Sun et al.
  # 2023, Lemma 2a).  Reference: exp of the difference of the two Rmpfr log
  # Psi values at (2.5, 900) and (1.5, 900).  Each log Psi is about 2.0e5,
  # one unit in its last place is 2.9e-11, and the difference of two carries
  # a few of them; 5e-10 allows eight on each side.  An error of 9e-9 in the
  # ratio exceeds it.
  expect_equal(mhn_mean(1.5, 1, 900), 450.0005555562414300, tolerance = 5e-10)
})

test_that("the quadrature resolves a peak far from the origin", {
  # For z > 0 the integrand u^(alpha-1) exp(-u^2 + z u) peaks near z/2 and is
  # about one unit wide there, while a panel reaching from the peak to the
  # origin, or far past it, is tens or hundreds of units long.  An adaptive
  # rule applied to such a panel can accept an estimate with the peak
  # under-resolved when its own error estimate happens to come out small; at
  # these points that costs 3e-7 to 4e-4 in log Psi.  References: Rmpfr at
  # 200 bits, by the power series for z <= 40 and the trapezoid rule in
  # x = e^v above (see the head of this part), each cross-checked on the
  # completed square to 1e-53.  The quadrature is within one unit in the last
  # place; 8 ulp allows for rounding that differs between platforms.
  refs <- list(
    list(alpha = 71.4,  z = 189,    ref = 9251.8731917488440),
    list(alpha = 36.2,  z = 58,     ref = 961.14524237761194),
    list(alpha = 50,    z = 640,    ref = 102683.91898177358),
    list(alpha = 82.8,  z = 965,    ref = 233312.96324779528),
    list(alpha = 0.15848931924611134, z = 23.901, ref = 141.99510060856439),
    list(alpha = 1e5,   z = 20.301, ref = 495575.75945024799)
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
  # At the user level: E[X] = exp(log Psi at alpha + 1 - log Psi at alpha)
  # (Sun et al. 2023, Lemma 2a), the reference being the exp of the
  # difference of the Rmpfr values at (72.4, 189) and (71.4, 189).  Each log
  # Psi is about 9.3e3 with one unit in its last place 1.8e-12, so the ratio
  # is good to about 1e-11; an under-resolved peak puts it 3.6e-4 off.
  expect_equal(mhn_mean(71.4, 1, 189), 94.871050472947534, tolerance = 1e-10)
})

test_that("log Psi at the strongest positive tilts follows the bell at z/2", {
  # With u = z/2 + t the integrand is (z/2 + t)^(alpha-1) exp(z^2/4 - t^2), so
  # log Psi = z^2/4 + (alpha - 1) log(z/2) + log(2 sqrt(pi)) + O(alpha^2 / z^2).
  # At these tilts the last term is far below a unit in the last place, and
  # the reference, formed in double precision, is within one or two of them;
  # 8 allows for that and for the quadrature's own rounding.  In u the
  # exponent is formed from u^2 and z u, each about z^2/4, and their rounding,
  # of order eps z^2/4, passes the range of exp from z of a few times 1e9,
  # where tanh_sinh (which takes a panel for alpha < 1) stops at an infinite
  # value; past z = 1e3 the integral is taken in log u instead.  For these
  # shapes log Psi is returned as Inf from z of about 1.9e154, where
  # (z/2)^2 reaches half the largest double; from 2.7e154 (z/2)^2 itself
  # overflows.
  for (a in c(0.5, 3, 150)) {
    for (z in c(1e9, 1e12, 1e16, 1e100, 1e154)) {
      ref <- (z / 2)^2 + (a - 1) * log(z / 2) + log(2 * sqrt(pi))
      err <- log_psi_abs_err(a, 1, z, ref)
      expect_lt(err / ulp_of(ref), 8,
                label = sprintf("ulp error, alpha=%g, z=%g", a, z))
    }
    for (z in c(3e154, 1e200, 1.7e308)) {
      expect_identical(mhn:::.mhn_log_normalizing_const(a, 1, z), Inf,
                       label = sprintf("log Psi, alpha=%g, z=%g", a, z))
    }
  }
})

test_that("log Psi at a large shape is accurate to a few ulp", {
  # The series truncation length of Lemma 10(d) bounds an absolute error, so
  # it grows with the size of the terms, by about log Gamma(alpha/2) / log 2,
  # and so does the rounding of the sum: at (1e5, 10) the series is 16 units
  # in the last place off.  The quadrature, which serves every alpha > 100,
  # is within one unit at these points; 8 allows for rounding that differs
  # between platforms.  References: the power series sum_n Gamma(alpha/2 +
  # n/2) z^n / n! in Rmpfr at 200 bits, cross-checked on the completed square
  # to 1e-55.
  refs <- list(
    list(alpha = 150, z = 5,  ref = 294.07203522683054),
    list(alpha = 1e4, z = 5,  ref = 37939.314203273897),
    list(alpha = 1e5, z = 10, ref = 493233.03221165301)
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
  # Far larger shapes.  In u the log integrand is formed from terms of size
  # alpha log u, whose rounding passes the range of exp from alpha near 1e17,
  # and the integral in log u, which forms no such terms, takes over.
  # References: the trapezoid rule in x = e^v (see the head of this part) in
  # Rmpfr at 240 bits, unchanged when its step is halved; and at the two
  # largest shapes lgamma(alpha/2) in double precision, from which log Psi
  # differs by about z sqrt(alpha/2), far below a unit in its last place.
  refs <- list(
    list(alpha = 3e17,  z = 10, ref = 5.7824117572240241e18),
    list(alpha = 1e18,  z = 1,  ref = 1.9876692247373545e19),
    list(alpha = 1e20,  z = 2,  ref = 2.2179277339801905e21),
    list(alpha = 1e100, z = 1,  ref = lgamma(5e99)),
    list(alpha = 1e300, z = 5,  ref = lgamma(5e299))
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
})

# ============================================================
# 15. The Lemma 10 series at a very small shape
#
# The ratio of consecutive A terms is (a + k - 1) z^2 / ((2k - 1) 2k) with
# a = alpha / 2.  For k = 1 the factor a + k - 1 is a itself; formed as
# (a + 1) - 1 it would lose every digit of a below the rounding of 1, so the
# test checks shapes where that matters: 1e-8 (a relative loss of 6e-9) and
# 1e-16 (where (a + 1) - 1 is 0).
# ============================================================

test_that("the Psi series keeps its accuracy for tiny alpha", {
  # log Psi is about 55.5 at z = 15 and 398 at z = 40.  The series is within
  # 1.5e-14 at z = 15 and 4e-13 (seven units in the last place) at z = 40,
  # where it sums some 2300 terms, so 2e-12 is a safe bound.  The loss
  # described above costs 3.0e-9 at alpha = 1e-8 and log 2 at alpha = 1e-16.
  # The series is called directly; through the dispatcher z = 15 takes the
  # series and z = 40 the quadrature, and both must meet the same bound.
  refs <- list(
    list(alpha = 1e-8,  z = 15, ref = 55.50970480619454753),
    list(alpha = 1e-16, z = 15, ref = 55.50970480179527835),
    list(alpha = 1e-8,  z = 40, ref = 398.2710338104256747),
    list(alpha = 1e-16, z = 40, ref = 398.2710337804871770)
  )
  for (r in refs) {
    got <- mhn:::.psi_series(r$alpha, 1, r$z)
    expect_lt(abs(got - r$ref), 2e-12,
              label = sprintf("series alpha=%g, z=%g", r$alpha, r$z))
    expect_lt(log_psi_abs_err(r$alpha, 1, r$z, r$ref), 2e-12,
              label = sprintf("dispatcher alpha=%g, z=%g", r$alpha, r$z))
  }
})

# ============================================================
# 16. Negative tilt
#
# For gamma < 0 the quadrature works in v = log u.  There the integrand
# u^alpha exp(-u^2 + z u) has a single peak, at the positive root u* of
# 2 u^2 - z u - alpha = 0, and nothing singular at the origin, and the tilt
# enters only through u*.  To the right of the peak the integrand is cut off
# where it has fallen by L = -log(tol) + 40.  To the left it falls only like
# exp(-alpha |v - log u*|), and it is followed down to the closed form for the
# head near u = 0 unless what lies beyond is below exp(-L) of the integral.
# The first test's points include shapes and tilts where a much shorter range
# leaves most of the mass out (alpha = 60, z = -15, and alpha = 300,
# z = -1e-6, where u^(alpha-1) exp(-u^2 + z u) peaks at u = 2.9 and 12.2),
# and ones with a strong tilt where the peak is squeezed towards the origin.
# ============================================================

test_that("the quadrature for gamma < 0 captures the whole integral", {
  # The quadrature is within a unit or two in the last place of these
  # references (4e-15 or less below alpha = 300); 1e-11 is a safe bound and
  # far below the error of a truncation that misses mass.
  refs <- list(
    list(alpha = 3,   z = -150,  ref = -13.64614464132463916),
    list(alpha = 3,   z = -500,  ref = -17.25757793241879826),
    list(alpha = 3.7, z = -150,  ref = -16.41890361287085715),
    list(alpha = 10,  z = -3,    ref = -2.393866566840880289),
    list(alpha = 60,  z = -15,   ref = 11.41740547302475706),
    list(alpha = 300, z = -1e-6, ref = 600.0094583180807761)
  )
  for (r in refs) {
    expect_lt(log_psi_abs_err(r$alpha, 1, r$z, r$ref), 1e-11,
              label = sprintf("alpha=%g, z=%g", r$alpha, r$z))
  }
})

test_that("the quadrature for gamma < 0 resolves a peak far from the origin", {
  # A large shape carries the peak of u^(alpha-1) exp(-u^2 + z u) away from
  # the origin under a negative tilt too, to about sqrt(alpha / 2) when the
  # tilt is weak, and the integrand is a bell narrower than one unit there.
  # On a panel that holds such a bell at one end of a much longer range --
  # from the origin to the peak, or tens of units past it -- an adaptive rule
  # can accept an estimate whose two embedded rules agree by chance while
  # both are off: at these points by 4e-8 to 1e-4 in log Psi, and by log 2
  # (the whole of one side of the bell) at the three largest shapes.
  # References: the trapezoid rule in x = e^v (see the head of this part) in
  # Rmpfr at 240 bits, unchanged when its step is halved.  The quadrature is
  # within one unit in the last place; 8 allows for rounding that differs
  # between platforms.
  refs <- list(
    list(alpha = 2795.0130456459246, z = -9.1122237683838687,
         ref = 8390.7086593793811),
    list(alpha = 213.6,     z = -30,       ref = 167.95178870731820),
    list(alpha = 515.81101, z = -46.06336, ref = 638.75849954639455),
    list(alpha = 14967.461, z = -53.589,   ref = 54975.918352304128),
    list(alpha = 1.2e7,     z = -10,       ref = 87619130.877287135),
    list(alpha = 1.5e7,     z = -1,        ref = 111225356.35478511),
    list(alpha = 1e8,       z = -10,       ref = 836305972.04537773)
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
  # At the user level: E[X] = exp(log Psi at alpha + 1 - log Psi at alpha)
  # (Sun et al. 2023, Lemma 2a), the reference being the exp of the
  # difference of the Rmpfr values.  Each log Psi is about 8.4e3, with one
  # unit in its last place 1.8e-12, so the ratio is good to about 1e-11; an
  # under-resolved bell puts it 9.6e-5 off.
  expect_equal(mhn_mean(2795.0130456459246, 1, -9.1122237683838687),
               35.171393872780271, tolerance = 1e-10)
})

test_that("log Psi next to alpha = 1 under a negative tilt is accurate", {
  # In u the integrand u^(alpha-1) exp(-u^2 + z u) is singular at the origin
  # for alpha < 1, and for alpha just above 1 it bends sharply within
  # (alpha - 1) / |z| of it, which a quadrature in u resolves only to about
  # its tolerance: 6e-14 to 1.4e-10 in log Psi at these points.  In log u the
  # integrand is analytic.  References: the power series
  # sum_n Gamma(alpha/2 + n/2) z^n / n! in Rmpfr, with the extra bits the
  # cancellation needs.  log Psi is of order one here, and its absolute error
  # is the relative error of Psi; 8 eps (1.8e-15) allows for rounding that
  # differs between platforms.
  refs <- list(
    list(alpha = 0.88,       z = -0.5, ref = 0.46201587301459357),
    list(alpha = 0.9,        z = -5,   ref = -0.74757077627934831),
    list(alpha = 0.99,       z = -2,   ref = -0.26072974399522425),
    list(alpha = 1 + 1e-8,   z = -0.5, ref = 0.31145012783356629),
    list(alpha = 1 + 1.4e-8, z = -15,  ref = -2.0236025771230017),
    list(alpha = 1 + 1e-6,   z = -5,   ref = -0.98445261494084635)
  )
  for (r in refs) {
    expect_lt(log_psi_abs_err(r$alpha, 1, r$z, r$ref),
              8 * .Machine$double.eps,
              label = sprintf("alpha=%.9g, z=%g", r$alpha, r$z))
  }
})

test_that("log Psi under a strong negative tilt follows Watson's lemma", {
  # Psi[alpha/2, z] = 2 int_0^Inf u^(alpha-1) exp(-u^2 - |z| u) du has the
  # asymptotic expansion 2 |z|^-alpha sum_k (-1)^k Gamma(alpha + 2k) /
  # (k! z^(2k)) (Watson's lemma), whose first term is the Gamma(alpha) limit
  # of the law.  References: that series in Rmpfr at 300 bits, summed until
  # its terms fall below 1e-60 of the first; at |z| >= 1e20 the first term
  # alone gives log Psi to double precision.  In log u the tilt only moves
  # the peak of the integrand, here to about alpha / |z|, down to 5e-301, and
  # the quadrature stays within a unit in the last place of log Psi; 8 allows
  # for rounding that differs between platforms.
  refs <- list(
    list(alpha = 0.5, z = -1e20,  ref = -21.760338806455813),
    list(alpha = 2,   z = -1e20,  ref = -91.410256539201882),
    list(alpha = 1e4, z = -1e20,  ref = -378416.60795518622),
    list(alpha = 0.5, z = -1e70,  ref = -79.324966131306951),
    list(alpha = 2,   z = -1e70,  ref = -321.66876583860648),
    list(alpha = 1e4, z = -1e70,  ref = -1529709.1544522091),
    list(alpha = 0.5, z = -1e300, ref = -344.12225182562219),
    list(alpha = 2,   z = -1e300, ref = -1380.8579086158675),
    list(alpha = 1e4, z = -1e300, ref = -6825654.8683385141),
    list(alpha = 2,   z = -38000, ref = -20.397535701012224),
    list(alpha = 2,   z = -38700, ref = -20.434042581481972),
    list(alpha = 0.5, z = -18000, ref = -3.6335513972693203),
    list(alpha = 0.5, z = -18700, ref = -3.6526272800814503)
  )
  for (r in refs) {
    err <- log_psi_abs_err(r$alpha, 1, r$z, r$ref)
    expect_lt(err / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
  # The density follows: at gamma / sqrt(beta) = -1e70 the law is
  # Gamma(2, rate 1e70) to double precision (the factor exp(-beta x^2) differs
  # from 1 by 1e-139 at these x).  dmhn is within 3e-14 of dgamma here.
  x <- c(1e-71, 1e-70, 5e-70)
  expect_equal(dmhn(x, 2, 1, -1e70), dgamma(x, 2, rate = 1e70),
               tolerance = 1e-12)
})

test_that("log Psi under a negative tilt holds at extreme parameters", {
  # Two limits give log Psi to double precision here.  Where
  # alpha (alpha + 1) / z^2 is below 1e-17, the first term of Watson's lemma
  # above, log 2 + log Gamma(alpha) - alpha log|z|.  And for a huge shape at
  # a modest tilt, log Gamma(alpha/2): log Psi exceeds it by about
  # z Gamma((alpha+1)/2) / Gamma(alpha/2), near z sqrt(alpha/2), which at
  # these points is far below a unit in the last place of log Gamma(alpha/2).
  # The shapes run from 1e-300, where Psi is about 2 / alpha, to 1e300, and
  # the tilts out to -1.7e308, where z^2 overflows.  Both references are
  # formed in double precision, to within a unit or two in the last place; 8
  # allows for that and for the quadrature's own rounding.
  gamma_limit <- list(c(1e-300, -1e20), c(1e-300, -1.7e308), c(1e-10, -1e300),
                      c(0.5, -1.7e308), c(2, -1.7e308), c(1e4, -1e150))
  for (p in gamma_limit) {
    ref <- log(2) + lgamma(p[1]) - p[1] * log(-p[2])
    err <- log_psi_abs_err(p[1], 1, p[2], ref)
    expect_lt(err / ulp_of(ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", p[1], p[2]))
  }
  large_shape <- list(c(1e100, -1), c(1e300, -1), c(1e300, -1e100))
  for (p in large_shape) {
    ref <- lgamma(p[1] / 2)
    err <- log_psi_abs_err(p[1], 1, p[2], ref)
    expect_lt(err / ulp_of(ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", p[1], p[2]))
  }
})

# ============================================================
# 17. The quadrature at a small shape and a moderate positive tilt
#
# The dispatcher sends 0 < z <= 20 to the Lemma 10 series whenever
# alpha <= 100, but the quadrature is also called directly.  For alpha < 1
# and z^2 > 8 (1 - alpha) the integrand u^(alpha-1) exp(-u^2 + z u) is
# singular at the origin, falls to a local minimum near (1 - alpha) / z and
# rises to its peak near z / 2.  Until that peak lies sqrt(-log(tol) + 40) =
# 7.6 or more from the origin, the quadrature takes the integral in log u,
# as for z < 0.
# ============================================================

test_that("the quadrature serves a small shape at a moderate positive tilt", {
  # In u a single rule over the range from the origin to the peak has to
  # resolve the singular end and the bell together: tanh_sinh over it stops
  # with a "singular point" error for alpha <= 0.06 and z in about
  # [12.8, 14.3], and at z = 8 the quadrature in u is about 300 to 700 units
  # in the last place off.  In log u the integrand has a single peak.
  # References: the power series in Rmpfr (see the head of this part); 8
  # units in the last place allow for rounding that differs between
  # platforms.
  refs <- list(
    list(alpha = 0.05, z = 12.8, ref = 40.473695927026448),
    list(alpha = 0.01, z = 13,   ref = 41.67445076932043),
    list(alpha = 1e-8, z = 14,   ref = 48.330080546648951),
    list(alpha = 0.03, z = 14.2, ref = 49.783941267633352),
    list(alpha = 0.01, z = 8,    ref = 15.926772109029756),
    list(alpha = 0.05, z = 8,    ref = 15.980143239471145)
  )
  for (r in refs) {
    got <- mhn:::.psi_integrate(r$alpha, 1, r$z)
    expect_lt(abs(got - r$ref) / ulp_of(r$ref), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", r$alpha, r$z))
  }
})

test_that("the quadrature serves a shape below 6 whose peak is near the origin", {
  # For 1 <= alpha < 6 and a peak within 7.6 of the origin the power
  # u^(alpha - 1) and the bell share one range in u, and the derivative of
  # order ceil(alpha - 1) of the power, unbounded at the origin, cost a single
  # rule there hundreds to thousands of units in the last place (2679 at
  # alpha = 3.77, z = 6e-5); the integral is taken in log u instead.  The
  # Lemma 10 series, which the dispatcher uses at these points, is within a
  # few units in the last place of log Psi (or of 1, where |log Psi| < 1), and
  # serves as the reference; 8 such units allow for both.
  for (az in list(c(3.7699740463867784, 6.1488785252442231e-05),
                  c(2.1677922326150254, 5.5067966841161251),
                  c(4.5, 1), c(5.9, 3), c(2, 0.5))) {
    a <- az[1]
    z <- az[2]
    ser <- mhn:::.psi_series(a, 1, z)
    got <- mhn:::.psi_integrate(a, 1, z)
    expect_lt(abs(got - ser) / ulp_of(max(1, abs(ser))), 8,
              label = sprintf("ulp error, alpha=%g, z=%g", a, z))
  }
})

test_that("the Lemma 10 series declines, rather than overflows, at a strong tilt", {
  # The truncation lengths of Lemma 10, about 0.72 z^2, are held as doubles
  # and tested against the ceiling of 1e6 terms before any is converted to
  # an integer type, which a 32-bit long could not hold from |z| of about
  # 5.5e4.  Past the ceiling the series returns NaN for the dispatcher to
  # hand over to the quadrature.
  for (z in c(2e3, 6e4, 1e5, 4e9, 1e10, 1e100)) {
    expect_true(is.nan(mhn:::.psi_series(3, 1, z)),
                label = sprintf("z=%g", z))
  }
})

test_that("the Lemma 10 series at a large shape keeps only the terms that count", {
  # Its absolute truncation bound would keep about log Gamma(alpha / 2) /
  # log 2 terms past the index C of Lemma 10, thousands at alpha = 1000; it
  # stops instead once the terms left are below 2^-56 of those kept, where
  # adding them could not change the sum.  Against the quadrature, within
  # about a unit in the last place here, the series is within a few; 16
  # allows for the rounding of its running sum over some hundred terms.
  for (az in list(c(300, 1), c(1000, 5), c(1e4, 10), c(1e4, 0.1))) {
    a <- az[1]
    z <- az[2]
    ser <- mhn:::.psi_series(a, 1, z)
    got <- mhn:::.psi_integrate(a, 1, z)
    expect_lt(abs(ser - got) / ulp_of(got), 16,
              label = sprintf("ulp error, alpha=%g, z=%g", a, z))
  }
})
