# Tests for mhn_mean, mhn_var, mhn_skewness, mhn_kurtosis, mhn_mode.
#
# Cross-checks each helper -- a ratio of Fox-Wright values or the Lemma 2c
# form where those are accurate, quadrature of the kernel centred on its peak
# otherwise -- against mhn_ref() of helper-mhn-cdf.R, which integrates the
# kernel of the law alone, and against high-precision references, plus
# specialised tests for special-case identities, Lemma 4 bounds, and
# parameter validation.
#
# Tolerances.  The package computes the variance, skewness and kurtosis either
# from exact closed forms or by quadratures that are asked for a relative
# accuracy of sqrt(eps) = 1.5e-8, and 1e-8 relative is the working tolerance
# it is written to.  Where a check is looser than that, the reference is the
# limiting factor, and the comment at the check says why.

# ============================================================
# 1. Mean: cross-check against numerical integration
# ============================================================

test_that("mhn_mean matches the kernel reference", {
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 3, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    # The reference asks integrate() for a relative 1e-12, or 1e-10 where
    # rounding prevents that, and agrees with the package to 1.1e-15 on this
    # grid.
    expect_equal(mhn_mean(a, b, g), mhn_ref(a, b, g)$mean,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 2. Variance: cross-check against numerical integration
# ============================================================

test_that("mhn_var matches the kernel reference", {
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    # As for the mean: the package agrees with the reference to 2.5e-14 here.
    expect_equal(mhn_var(a, b, g), mhn_ref(a, b, g)$var,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 3. Variance upper bound: Var(X) <= 1/(2 beta) when alpha >= 1
#    (Sun et al. 2023, Lemma 4c)
# ============================================================

test_that("variance upper bound holds for alpha >= 1", {
  # The grid reaches z = gamma / sqrt(beta) from -3200 to 28000, so it covers
  # the strong positive tilts at which the Lemma 2c form magnifies any error
  # in the mean (by about z^2 / 2), and alpha up to 500, where it magnifies it
  # by about 4 alpha.  At alpha = 1 and a large z the law is a normal truncated
  # far below its mean, whose variance falls short of the bound 1/(2 beta) by
  # a relative amount of order exp(-z^2 / 4), so there the bound is tight and
  # an error of 1e-9 upwards would show.  The slack is relative, because the
  # bound itself runs from 0.1 to 500 on this grid.
  test_cases <- expand.grid(
    alpha = c(1, 1.5, 2, 3, 10, 100, 500),
    beta = c(0.001, 0.01, 1, 5),
    gamma = c(-100, -20, -5, 0, 5, 20, 50, 100, 300, 700, 900)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    v <- mhn_var(a, b, g)
    expect_true(v <= (1 + 1e-9) / (2 * b),
                label = sprintf("alpha=%g, beta=%g, gamma=%g", a, b, g))
  }
})

# ============================================================
# 4. Special-case moments: gamma = 0 (sqrt-Gamma reduction)
# ============================================================

test_that("mhn_mean for gamma=0 matches sqrt-Gamma mean", {
  # E(X) = Gamma((alpha+1)/2) / (sqrt(beta) * Gamma(alpha/2))
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 5)) {
      expected <- gamma((a + 1) / 2) / (sqrt(b) * gamma(a / 2))
      expect_equal(mhn_mean(a, b, 0), expected,
                   tolerance = 1e-10,
                   label = sprintf("alpha=%.1f, beta=%.1f", a, b))
    }
  }
})

test_that("mhn_var for gamma=0 matches sqrt-Gamma variance", {
  # Var(X) = alpha/(2 beta) - [Gamma((alpha+1)/2) / (sqrt(beta) Gamma(alpha/2))]^2
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 5)) {
      mu <- gamma((a + 1) / 2) / (sqrt(b) * gamma(a / 2))
      expected <- a / (2 * b) - mu^2
      expect_equal(mhn_var(a, b, 0), expected,
                   tolerance = 1e-10,
                   label = sprintf("alpha=%.1f, beta=%.1f", a, b))
    }
  }
})

# ============================================================
# 5. Moment recurrence consistency (Sun et al. 2023, Lemma 2b)
# ============================================================

test_that("mhn_mean satisfies E(X^2) = alpha/(2 beta) + gamma/(2 beta) E(X)", {
  # Sun et al. (2023, Lemma 2b), with E(X^2) from the kernel reference and
  # E(X) from the package: an error in mhn_mean moves the right-hand side by
  # gamma / 2 times it.  (mhn_var is formed from mhn_mean by Lemma 2c, so
  # mhn_var + mhn_mean^2 would satisfy the recurrence whatever the mean.)
  # The two sides agree to 1.3e-15 here.
  for (a in c(0.5, 2, 5)) {
    for (g in c(-2, 0, 2)) {
      ref <- mhn_ref(a, 1, g)
      ex2_from_kernel <- ref$var + ref$mean^2
      ex2_from_rec <- a / 2 + g / 2 * mhn_mean(a, 1, g)  # beta = 1
      expect_equal(ex2_from_kernel, ex2_from_rec,
                   tolerance = 1e-10,
                   label = sprintf("alpha=%.1f, gamma=%.1f", a, g))
    }
  }
})

# ============================================================
# 6. Skewness: cross-check against numerical integration
# ============================================================

test_that("mhn_skewness matches the kernel reference", {
  test_cases <- list(
    list(alpha = 2, beta = 1, gamma = 0),
    list(alpha = 2, beta = 1, gamma = 2),
    list(alpha = 3, beta = 2, gamma = -1),
    list(alpha = 0.5, beta = 1, gamma = 1),
    list(alpha = 5, beta = 1, gamma = 3)
  )
  for (tc in test_cases) {
    a <- tc$alpha; b <- tc$beta; g <- tc$gamma
    # The reference forms the third central moment about its own mean; the
    # package agrees with it to 4e-15 here, and with 400-bit references to
    # 4e-15.
    expect_equal(mhn_skewness(a, b, g), mhn_ref(a, b, g)$skewness,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 7. Excess kurtosis: cross-check against numerical integration
# ============================================================

test_that("mhn_kurtosis matches the kernel reference", {
  test_cases <- list(
    list(alpha = 2, beta = 1, gamma = 0),
    list(alpha = 2, beta = 1, gamma = 2),
    list(alpha = 3, beta = 2, gamma = -1),
    list(alpha = 5, beta = 1, gamma = 3)
  )
  for (tc in test_cases) {
    a <- tc$alpha; b <- tc$beta; g <- tc$gamma
    # Subtracting 3 from a fourth moment near 2.95 magnifies the reference's
    # error by up to about 60 here; the package agrees with it to 3e-14, and
    # 1e-10 leaves room for that magnification.
    expect_equal(mhn_kurtosis(a, b, g), mhn_ref(a, b, g)$excess_kurtosis,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 8. Parameter validation
# ============================================================

test_that("moment functions error on invalid parameters", {
  expect_error(mhn_mean(-1, 1, 0), "alpha must be positive")
  expect_error(mhn_var(1, -1, 0), "beta must be positive")
  expect_error(mhn_skewness(0, 1, 0), "alpha must be positive")
  expect_error(mhn_kurtosis(1, 0, 0), "beta must be positive")
})

# ============================================================
# 9. Mode: closed-form verification
# ============================================================

test_that("mhn_mode for alpha > 1: closed form", {
  test_cases <- expand.grid(
    alpha = c(1.5, 2, 3, 5, 10),
    beta = c(0.5, 1, 5),
    gamma = c(-5, -1, 0, 1, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    expected <- (g + sqrt(g^2 + 8 * b * (a - 1))) / (4 * b)
    expect_equal(mhn_mode(a, b, g), expected,
                 tolerance = 1e-12,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

test_that("mhn_mode for alpha = 1: max(0, gamma/(2 beta))", {
  expect_equal(mhn_mode(1, 1, 2), 1)     # gamma/(2 beta) = 1 > 0
  expect_equal(mhn_mode(1, 1, -2), 0)    # gamma/(2 beta) = -1 < 0
  expect_equal(mhn_mode(1, 1, 0), 0)     # gamma/(2 beta) = 0
  expect_equal(mhn_mode(1, 2, 3), 0.75)  # gamma/(2 beta) = 0.75
})

test_that("mhn_mode for alpha < 1: NA when no interior mode", {
  # gamma <= 0 with alpha < 1: always NA (density monotone decreasing).
  expect_true(is.na(mhn_mode(0.5, 1, -1)))
  expect_true(is.na(mhn_mode(0.5, 1, 0)))
  expect_true(is.na(mhn_mode(0.3, 2, -5)))
})

test_that("mhn_mode for alpha < 1, gamma > 0: conditional", {
  # An interior local maximum exists iff alpha > 1 - gamma^2 / (8 beta)
  # (Sun et al. 2023, Lemma 3c, which prints >=; the equality case is tested
  # below).
  # alpha=0.9, beta=1, gamma=2: 1 - 4/8 = 0.5 < 0.9, mode exists.
  mode_val <- mhn_mode(0.9, 1, 2)
  expect_true(!is.na(mode_val) && mode_val > 0)

  # alpha < 1 - gamma^2 / (8 beta): NA.
  # alpha=0.1, beta=1, gamma=0.5: 1 - 0.25/8 = 0.96875 > 0.1, no mode.
  expect_true(is.na(mhn_mode(0.1, 1, 0.5)))
})

test_that("mhn_mode for alpha < 1, gamma > 0: the Lemma 3c value", {
  # The larger root of 2 beta x^2 - gamma x + 1 - alpha = 0 (Sun et al. 2023,
  # Lemma 3c) at ordinary parameters, and a local maximum of the density there.
  for (p in list(c(0.5, 1, 5), c(0.9, 1, 2), c(0.3, 1, 4), c(0.7, 2.5, 3))) {
    alpha <- p[1]
    beta <- p[2]
    gamma <- p[3]
    m <- mhn_mode(alpha, beta, gamma)
    expect_equal(m,
                 (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta),
                 tolerance = 1e-14)
    expect_true(all(dmhn(m, alpha, beta, gamma, log = TRUE) >
                      dmhn(m * c(0.999, 1.001), alpha, beta, gamma,
                           log = TRUE)))
  }
})

test_that("mhn_mode is NA at the Lemma 3c boundary, where there is only an inflection point", {
  # At alpha = 1 - gamma^2 / (8 beta) the two stationary points
  # (gamma +- sqrt(gamma^2 + 8 beta (alpha - 1))) / (4 beta) coincide at
  # gamma / (4 beta), and d/dx log f = -(2 beta / x) (x - gamma / (4 beta))^2
  # is negative on both sides of it: the density is strictly decreasing, so
  # there is no local maximum to report.  The three triples are exact in
  # binary, so the discriminant is exactly 0.
  expect_true(is.na(mhn_mode(0.875, 1, 1)))
  expect_true(is.na(mhn_mode(0.75, 2, 2)))
  expect_true(is.na(mhn_mode(0.5, 0.25, 1)))
  # The density really is decreasing through gamma / (4 beta) = 0.25.
  x <- 0.25 + c(-1e-2, -1e-3, 0, 1e-3, 1e-2)
  expect_true(all(diff(dmhn(x, 0.875, 1, 1, log = TRUE)) < 0))
  # Just inside the region the maximum exists, next to the merged point.
  m <- mhn_mode(0.875 + 1e-4, 1, 1)
  expect_true(!is.na(m) && abs(m - 0.25) < 0.01)
})

test_that("mhn_mode does not snap alpha near 1 to the alpha = 1 formula", {
  # Near (alpha, gamma) = (1, 0) the Lemma 3b mode sqrt((alpha - 1) / (2 beta))
  # has an infinite slope in alpha, so any tolerance band around alpha = 1
  # returns 0 for a mode far from 0; and below 1 there is no mode at all when
  # gamma < 0 (Lemma 3d) or gamma = 0 (the derivative of log f is negative),
  # or when gamma^2 < 8 beta (1 - alpha) (Lemma 3d).  The expected
  # values are the Lemma 3b/3c formulas at the alpha actually stored, which
  # R evaluates to within a few ulps here.
  a <- 1 + 1e-8
  expect_equal(mhn_mode(a, 1, 0), sqrt((a - 1) / 2), tolerance = 1e-12)
  expect_equal(mhn_mode(a, 1e-6, 0), sqrt((a - 1) / 2e-6), tolerance = 1e-12)
  expect_equal(mhn_mode(a, 1, 1e-4), (1e-4 + sqrt(1e-8 + 8 * (a - 1))) / 4,
               tolerance = 1e-12)
  expect_true(dmhn(mhn_mode(a, 1e-6, 0), a, 1e-6, 0) > 0)
  a <- 1 - 1e-8
  expect_true(is.na(mhn_mode(a, 1, -1)))
  expect_true(is.na(mhn_mode(a, 1, 0)))
  expect_true(is.na(mhn_mode(a, 1, 1e-5)))  # 1e-10 < 8e-8: no local maximum
  a <- 1 - 1e-9
  expect_equal(mhn_mode(a, 1, 1e-4), (1e-4 + sqrt(1e-8 + 8 * (a - 1))) / 4,
               tolerance = 1e-12)
})

test_that("mhn_mode holds where gamma^2 or 8 beta (alpha - 1) overflows", {
  # The mode is the positive root of 2 beta x^2 - gamma x - (alpha - 1) = 0,
  # about gamma / (2 beta) or (alpha - 1) / |gamma| under a strong tilt, and
  # representable although gamma^2 overflows from |gamma| = 1.34e154, and
  # 8 beta (alpha - 1) or 2 beta can overflow on their own.  The first four
  # values are exact to double precision: (1e200 + sqrt(1e400 + 8)) / 4 and
  # 2 / (sqrt(1e400 + 8) + 1e200) differ from 5e199 and 1e-200 by a relative
  # 1e-400, and at (1e300, 1, +-1e150) the discriminant is 9e300 to within
  # 8, whose root is 3e150.  With a huge beta the mode follows from its value
  # at beta = 1 by scale equivariance, sqrt(beta) X ~ MHN(alpha, 1,
  # gamma / sqrt(beta)); 1e303, 1e308 and their square roots carry a relative
  # rounding of about 1e-16, hence the tolerance.  At -1.7e308 the mode,
  # 1/1.7e308, is subnormal and held to its spacing, 4.9e-324.  The tiny
  # values are compared as ratios to 1: expect_equal() compares absolute
  # differences once the expected value is smaller than the tolerance, and
  # would then pass for any tiny value.
  expect_equal(mhn_mode(2, 1, 1e200), 5e199, tolerance = 1e-15)
  expect_equal(mhn_mode(2, 1, -1e200) / 1e-200, 1, tolerance = 1e-15)
  expect_equal(mhn_mode(1e300, 1, 1e150), 1e150, tolerance = 1e-15)
  expect_equal(mhn_mode(1e300, 1, -1e150), 5e149, tolerance = 1e-15)
  expect_equal(mhn_mode(2, 1, -1.7e308) / (1 / 1.7e308), 1, tolerance = 1e-12)
  sb <- sqrt(1e303)
  expect_equal(mhn_mode(1e10, 1e303, 1e5 * sb) / (mhn_mode(1e10, 1, 1e5) / sb),
               1, tolerance = 1e-14)
  expect_equal(mhn_mode(1e10, 1e303, -1e5 * sb) /
                 (mhn_mode(1e10, 1, -1e5) / sb), 1, tolerance = 1e-14)
  sb <- sqrt(1e308)
  expect_equal(mhn_mode(0.5, 1e308, 1e6 * sb) / (mhn_mode(0.5, 1, 1e6) / sb),
               1, tolerance = 1e-14)
  expect_equal(mhn_mode(1, 1e308, 1e300), 1e300 / 2 / 1e308, tolerance = 1e-15)
})

# ============================================================
# 10. Density at mode is a local maximum
# ============================================================

test_that("density at mode is a local maximum", {
  test_cases <- list(
    list(alpha = 2, beta = 1, gamma = 0),
    list(alpha = 2, beta = 1, gamma = 3),
    list(alpha = 5, beta = 2, gamma = -1),
    list(alpha = 3, beta = 1, gamma = -3)
  )
  eps <- 1e-5
  for (tc in test_cases) {
    a <- tc$alpha; b <- tc$beta; g <- tc$gamma
    mode_val <- mhn_mode(a, b, g)
    d_mode <- dmhn(mode_val, a, b, g)
    d_left <- dmhn(mode_val - eps, a, b, g)
    d_right <- dmhn(mode_val + eps, a, b, g)
    expect_true(d_mode >= d_left - 1e-10,
                label = sprintf("mode >= left, a=%.1f g=%.1f", a, g))
    expect_true(d_mode >= d_right - 1e-10,
                label = sprintf("mode >= right, a=%.1f g=%.1f", a, g))
  }
})

# ============================================================
# 11. Mean >= mode for alpha > 1 (Sun et al. 2023, Lemma 4a)
# ============================================================

test_that("mean >= mode for alpha > 1 (Sun et al. 2023, Lemma 4a)", {
  test_cases <- expand.grid(
    alpha = c(1.5, 2, 3, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    expect_true(mhn_mean(a, b, g) >= mhn_mode(a, b, g) - 1e-10,
                label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

test_that("mean >= mode holds under a strong positive tilt", {
  # Under a strong tilt the mean exceeds the mode by very little: the gap is
  # about 2 / z^3 at alpha = 2 (see the exact values in section 15), and over
  # this grid it is at least 2.5e-12 of the mode (alpha = 1.5, z = 945; from a
  # 200-bit trapezoid rule on the peak-centred kernel).  A mean held to
  # rounding keeps it; one that is only 1e-9 accurate does not.  The slack of
  # 1e-13 relative covers the rounding of the two sides.
  test_cases <- expand.grid(
    alpha = c(1.5, 2, 3, 5, 10, 30),
    beta = c(0.01, 1, 100),
    z = c(100, 300, 600, 900, 945)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$z[i] * sqrt(b)
    mode <- mhn_mode(a, b, g)
    expect_true(mhn_mean(a, b, g) >= mode * (1 - 1e-13),
                label = sprintf("alpha=%g, beta=%g, z=%g", a, b, test_cases$z[i]))
  }
})

test_that("mean and variance at alpha = 1 keep to the truncated-normal mean and Lemma 4c next to z = 12", {
  # At alpha = 1 the law is a normal truncated to x > 0, with mean
  # z / 2 + lambda / sqrt(2) and variance (1 - z lambda / sqrt(2) - lambda^2)
  # / 2 at beta = 1, lambda = phi(z / sqrt(2)) / Phi(z / sqrt(2)), and for
  # z from 10.8 to 12 the mean exceeds the mode z / 2 by less than 1e-13 and
  # the variance falls short of the Lemma 4c bound 1/2 by less than 1e-12,
  # less than the error the moments are allowed (1e-8 relative, as checked
  # against the exact values below), so the bounds are checked apart.  At
  # alpha = 1 mean >= mode holds because the truncated-normal mean exceeds
  # z / 2 by lambda / sqrt(2) > 0; Lemma 4a is stated for alpha > 1.  Just
  # above alpha = 1 the same bounds hold (Sun et al. 2023, Lemmas 4a and 4c).
  z <- seq(10.8, 12, by = 0.01)
  lam <- dnorm(z / sqrt(2)) / pnorm(z / sqrt(2))
  tn_mean <- z / 2 + lam / sqrt(2)
  tn_var <- (1 - z * lam / sqrt(2) - lam^2) / 2
  m <- vapply(z, function(g) mhn_mean(1, 1, g), numeric(1))
  v <- vapply(z, function(g) mhn_var(1, 1, g), numeric(1))
  expect_lt(max(abs(m / tn_mean - 1)), 1e-8)
  expect_lt(max(abs(v / tn_var - 1)), 1e-8)
  for (a in c(1, 1 + 1e-12)) {
    m <- vapply(z, function(g) mhn_mean(a, 1, g), numeric(1))
    v <- vapply(z, function(g) mhn_var(a, 1, g), numeric(1))
    mode <- vapply(z, function(g) mhn_mode(a, 1, g), numeric(1))
    expect_true(all(m >= mode), label = sprintf("mean >= mode, alpha=%.13g", a))
    expect_true(all(v <= 0.5), label = sprintf("var <= 1/2, alpha=%.13g", a))
  }
})

test_that("the variance next to alpha = 1 keeps to Lemma 4c past z = 12", {
  # Above z = 12 the variance comes from the quadrature.  At alpha = 1 and
  # just above it is the truncated-normal variance, (1 - z lambda / sqrt(2)
  # - lambda^2) / 2 with lambda = phi(z / sqrt(2)) / Phi(z / sqrt(2)), which
  # lies within 4e-16 of the Lemma 4c bound 1/2 here (Sun et al. 2023,
  # Lemma 4c).  That is less than the quadrature's error of a few ulps, and
  # at 62 of these points, for either shape, the quadrature by itself lands
  # above 1/2, by up to 2.2e-16; for alpha >= 1 its result is capped at 1/2.
  # 1e-15 against the exact values allows for those few ulps.
  z <- seq(12.01, 14, by = 0.01)
  lam <- dnorm(z / sqrt(2)) / pnorm(z / sqrt(2))
  tn_var <- (1 - z * lam / sqrt(2) - lam^2) / 2
  for (a in c(1, 1 + 2^-52)) {
    v <- vapply(z, function(g) mhn_var(a, 1, g), numeric(1))
    info <- sprintf("alpha=%.17g", a)
    expect_true(all(v <= 0.5), label = paste("var <= 1/2,", info))
    expect_lt(max(abs(v - tn_var)), 1e-15, label = info)
  }
})

test_that("the mean at a strong negative tilt stays below the Lemma 4a bound", {
  # Under a strong negative tilt the mean, about alpha/|z| (1 - 2 (alpha +
  # 1)/z^2), lies below the upper bound (z + sqrt(z^2 + 8 alpha)) / 4 of
  # Sun et al. (2023, Lemma 4a; it holds for every alpha > 0, as the text
  # after it notes) by only a relative 2/z^2, 2e-12 at z = -1e6.
  # References: the expansion of the mean in powers of 1/z^2 about the
  # Gamma(alpha, |z|) limit, summed in Rmpfr at 256 bits to a truncation
  # error below 1e-30, and a trapezoid rule in log x at 200 bits, which
  # agree to every digit.  The package meets them to 3e-16 here; 1e-14
  # allows for that, and is far inside the gap to the bound.
  ref <- list(c(10, -100, 0.0997810060131453513),
              c(300, -1e3, 0.2998196174753803844),
              c(1e4, -1e5, 0.0999997999808001914),
              c(1e4, -1e6, 0.0099999997999800083),
              c(3e6, -1e6, 2.9999820002099970573))
  for (r in ref) {
    a <- r[1]; z <- r[2]
    info <- sprintf("alpha=%g z=%g", a, z)
    m <- mhn_mean(a, 1, z)
    expect_lt(abs(m / r[3] - 1), 1e-14, label = info)
    expect_lt(m, 2 * a / (sqrt(z^2 + 8 * a) - z), label = info)
  }
})

test_that("the mean for alpha < 6 at a strong negative tilt obeys Lemma 4a", {
  # Below alpha = 6 the mean for z < -30 is the Lemma 2a ratio, whose error,
  # up to about 5e-14 relative for |z| up to 1e14, exceeds the gap 2/z^2
  # between the mean and the upper bound of Lemma 4a once |z| passes about
  # 1e7; the ratio is then capped at the bound.  At these points the ratio
  # by itself lies 1.3e-14 to 3e-14 above the bound, and the mean lies
  # within 3e-22 of it.
  # References: the expansion of the mean in powers of 1/z^2 about the
  # Gamma(alpha, |z|) limit, summed in Rmpfr at 256 bits until a term falls
  # below 2^-250 of the sum; its terms alternate, and the remainder of the
  # series of exp(-x^2) is below its next term, so that bounds the
  # truncation error.  The bound is formed here to within 1.2 eps, and the
  # mean is allowed 4 eps above it; the package returns its own rounding of
  # the bound, and 1e-15 against the references allows for that.
  ref <- list(c(5.7, -8e11, 7.1250000000000002220e-12),
              c(4.5, -4e13, 1.1250000000000000000e-13),
              c(5.99, -9e10, 6.6555555555555557924e-11),
              c(3, -1e12, 3.0000000000000000000e-12),
              c(1.5, -5e11, 3.0000000000000000000e-12))
  for (r in ref) {
    a <- r[1]; z <- r[2]
    info <- sprintf("alpha=%g z=%g", a, z)
    m <- mhn_mean(a, 1, z)
    bound <- 2 * a / (sqrt(z^2 + 8 * a) - z)
    expect_lte(m, bound * (1 + 4 * .Machine$double.eps), label = info)
    expect_lt(abs(m / r[3] - 1), 1e-15, label = info)
  }
})

# ============================================================
# 12. Mode parameter validation
# ============================================================

test_that("mhn_mode errors on invalid parameters", {
  expect_error(mhn_mode(-1, 1, 0), "alpha must be positive")
  expect_error(mhn_mode(1, -1, 0), "beta must be positive")
})

# ============================================================
# 13. Large alpha: the kernel is integrated far enough from its peak
# ============================================================
#
# The central moments are integrals of the kernel over a window around its
# peak m.  Relative to the peak its logarithm falls by
# (alpha-1) [s/m - log(1 + s/m)] + beta s^2 over a distance s.  Near the peak
# the first term is only about (alpha-1) s^2 / (2 m^2), so on its own it needs
# |s| of about m sqrt(2 L / (alpha-1)) to fall by L nats; a window cut shorter
# than that at large alpha stops a few nats from the peak, where the tails
# still carry much of the third and fourth moments.  At alpha = 1000 such a
# window gives skewness 0.0011 and excess kurtosis -0.63, against exact values
# of 0.022 and 7.5e-7.

test_that("skewness and kurtosis at gamma = 0 match exact values for large alpha", {
  # At gamma = 0, X^2 ~ Gamma(alpha/2, rate = beta), so
  # E X^k = Gamma((alpha + k)/2) / (Gamma(alpha/2) beta^(k/2)) exactly.  The
  # central moments below were formed from these raw moments in Rmpfr at 512
  # bits -- in double the differences cancel completely at this alpha -- and
  # rounded to 17 significant digits.
  ref <- list(
    list(alpha = 100,  skew = 0.071152520231415717,  exkurt = 7.611758446649681e-05),
    list(alpha = 1000, skew = 0.022374655011123699,  exkurt = 7.5112429301003366e-07),
    list(alpha = 1e4,  skew = 0.0070715097530468498, exkurt = 7.5011249296488305e-09)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g", r$alpha)
    # Each integral is requested to a relative tolerance of sqrt(eps) = 1.5e-8,
    # and Gauss-Kronrod on these smooth bells does far better than it is asked
    # (about 1e-13 here).  1e-8 relative on the skewness keeps that margin and
    # still catches a window cut short, whose error is 3.5e-6 relative at
    # alpha = 100 and of order one at alpha = 1000.
    expect_equal(mhn_skewness(r$alpha, 1, 0), r$skew, tolerance = 1e-8,
                 info = info)
    # The excess kurtosis is c4 / c2^2 - 3 with c4 / c2^2 close to 3, and it
    # is compared in absolute terms because the value itself is as small as
    # 7.5e-9.  Rounding the ratio alone costs about 1e-15, and the package
    # meets these values to 2e-15; 1e-12 leaves room for that while being
    # only 1.3e-4 of the value at alpha = 1e4 (1.3e-6 at alpha = 1000).
    expect_lt(abs(mhn_kurtosis(r$alpha, 1, 0) - r$exkurt), 1e-12, label = info)
  }
})

test_that("skewness and kurtosis for large alpha under a moderate tilt match references", {
  # References: E X^k = Psi[(alpha+k)/2, z] / Psi[alpha/2, z] with the
  # Fox-Wright series of section 17, in Rmpfr at 400 bits with 2000 terms and
  # again at 700 bits with 3000, with the same digits.  The window must reach
  # as far from the peak under a tilt as at gamma = 0, and at (100, 1, -1) one
  # cut short leaves the excess kurtosis 0.4 percent low.  Tolerances as in
  # the test above; the package meets these to 6e-15 relative (skewness) and
  # 2e-15 absolute (excess kurtosis).
  ref <- list(
    list(alpha = 100, z = -1, skew = 0.074935509973977382, exkurt = 0.00065604476332243981),
    list(alpha = 100, z = 5, skew = 0.053456668457780056, exkurt = -0.0017594951231349827),
    list(alpha = 50, z = -8, skew = 0.16140127856820974, exkurt = 0.021530946422444557),
    list(alpha = 1000, z = 2, skew = 0.021628933575181995, exkurt = -3.1374343828150563e-05)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g z=%g", r$alpha, r$z)
    expect_equal(mhn_skewness(r$alpha, 1, r$z), r$skew, tolerance = 1e-8,
                 info = info)
    expect_lt(abs(mhn_kurtosis(r$alpha, 1, r$z) - r$exkurt), 1e-12, label = info)
  }
})

test_that("moments under a strong negative tilt match the Gamma limit for large alpha", {
  # For gamma < 0, substituting y = |gamma| x turns the kernel into
  # y^(alpha-1) e^{-y} e^{-eps y^2}, eps = beta / gamma^2, and expanding the
  # last factor gives
  #   int x^(alpha+k-1) e^{-beta x^2 + gamma x} dx
  #     = |gamma|^-(alpha+k) sum_j (-eps)^j / j! Gamma(alpha + k + 2 j).
  # The series diverges in the end, but its terms keep falling until j is
  # about 1 / (4 eps); with eps = 1e-10 each of the first terms is at most 1e-4
  # times the one before, so 60 terms summed in Rmpfr at 512 bits give the
  # values below to all digits shown.  (At alpha = 100, gamma = -1000 the same
  # sum agrees with R's integrate() at rel.tol = 1e-13 to 14 digits.)  They lie
  # within 1e-6 relative of the Gamma(alpha, |gamma|) limit -- skewness
  # 2 / sqrt(alpha), excess kurtosis 6 / alpha, variance alpha / gamma^2 --
  # which is the law this tilt makes the MHN approach.
  ref <- list(
    list(alpha = 100,  var = 9.9999993940000405e-09, skew = 0.1999999939400004,
         exkurt = 0.05999999515200042),
    list(alpha = 1000, var = 9.9999939940040102e-08, skew = 0.063245534210740334,
         exkurt = 0.0059999951952040944)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g gamma=-1e5", r$alpha)
    # mhn_var takes this point from the same quadrature as the higher moments,
    # because the Lemma 2c closed form cancels at this tilt; the relative
    # tolerances follow the sqrt(eps) request of that quadrature with the same
    # margin as above.
    expect_close(mhn_var(r$alpha, 1, -1e5), r$var, 1e-8, info = info)
    expect_equal(mhn_skewness(r$alpha, 1, -1e5), r$skew, tolerance = 1e-8,
                 info = info)
    expect_equal(mhn_kurtosis(r$alpha, 1, -1e5), r$exkurt, tolerance = 1e-7,
                 info = info)
  }
})

test_that("mean and variance at alpha = 1e6 with a strong positive tilt", {
  # At z = gamma / sqrt(beta) = 1000 both the mean and the variance are taken
  # from the central-moment quadrature, because the ratio of Fox-Wright values
  # the mean would otherwise use cancels at this tilt, and the Lemma 2c form of
  # the variance has derivative gamma / (2 beta) - 2 mu, about -1500, with
  # respect to the mean: an error of 1e-4 in the mean would move it from 1/3
  # to 1/2.  The peak is at m = 999.99967, where the Laplace approximation
  # 1 / (2 beta + (alpha - 1) / m^2) gives 1/3.  The references are the
  # trapezoid rule on the peak-centred kernel in Rmpfr at 200 bits, over
  # +-25 (43 standard deviations); the trapezoid error for such an integrand
  # falls like exp(-c / h), and steps of 0.03 and 0.02 agree in every digit
  # shown.  The same procedure reproduces the exact alpha = 1000, gamma = 0
  # values above to 17 digits.
  #
  # The mean is the peak m plus an offset of +1.1e-4 from quadrature, so
  # 1e-12 relative (1e-9 absolute) is loose for it.  The variance is a
  # quadrature asked for sqrt(eps) relative, hence 1e-8.
  expect_equal(mhn_mean(1e6, 1, 1000), 999.99977777776132, tolerance = 1e-12)
  expect_equal(mhn_var(1e6, 1, 1000), 0.33333330864194055, tolerance = 1e-8)
})

# ============================================================
# 14. Scale equivariance in beta
# ============================================================

test_that("skewness and kurtosis depend on beta only through z = gamma / sqrt(beta)", {
  # If X ~ MHN(alpha, beta, gamma) then sqrt(beta) X ~ MHN(alpha, 1, z) with
  # z = gamma / sqrt(beta) (by a change of variable in the density of Sun et
  # al. 2023, Lemma 1a), and skewness and
  # kurtosis are unchanged by rescaling.  So at fixed (alpha, z) they must not
  # move with beta at all.  The points cover a positive tilt, an interior peak
  # with alpha < 1, a mild negative tilt and a strong one; beta runs over
  # 32 orders of magnitude.  The moments are computed for beta = 1 from z, so
  # what remains is the rounding of z itself, and 1e-10 relative leaves ample
  # room for that while catching any dependence of the quadrature on beta.
  betas <- 10^c(-16, -8, 0, 8, 16)
  for (p in list(c(1.5, 5), c(0.5, 5), c(2.5, -1), c(1.5, -20))) {
    a <- p[1]; z <- p[2]
    s1 <- mhn_skewness(a, 1, z)
    k1 <- mhn_kurtosis(a, 1, z)
    for (b in betas) {
      info <- sprintf("alpha=%g z=%g beta=%g", a, z, b)
      expect_equal(mhn_skewness(a, b, z * sqrt(b)), s1, tolerance = 1e-10,
                   info = info)
      expect_equal(mhn_kurtosis(a, b, z * sqrt(b)), k1, tolerance = 1e-10,
                   info = info)
    }
  }
})

test_that("moments at extreme beta equal the beta = 1 values at the same z", {
  # The same scale identity, far outside the range above.  At beta = 1e300
  # the fourth central moment is of order 1e-600, and a quadrature in absolute
  # units underflows (excess kurtosis exactly -3, skewness exactly 0); at
  # beta = 1e-300 it overflows.  Skewness and kurtosis must equal their
  # beta = 1 values, and so must beta Var(X) and sqrt(beta) E(X); 1e-12
  # relative covers the rounding of z = gamma / sqrt(beta).
  betas <- 10^c(-300, -200, -150, 150, 200, 300)
  for (p in list(c(2, 0.7), c(3, 1.3), c(1.5, 5), c(1.5, -20), c(0.5, 0.7),
                 c(0.5, -1), c(0.05, 5))) {
    a <- p[1]; z <- p[2]
    ref <- c(mhn_mean(a, 1, z), mhn_var(a, 1, z), mhn_skewness(a, 1, z),
             mhn_kurtosis(a, 1, z))
    for (b in betas) {
      g <- z * sqrt(b)
      got <- c(mhn_mean(a, b, g) * sqrt(b), mhn_var(a, b, g) * b,
               mhn_skewness(a, b, g), mhn_kurtosis(a, b, g))
      expect_close(got, ref, 1e-12,
                   info = sprintf("alpha=%g z=%g beta=%g", a, z, b))
    }
  }
  # Under a strong negative tilt the variance comes from the quadrature (the
  # Lemma 2c form cancels), which is where absolute units hurt it most.
  for (p in list(c(1.5, -1e7), c(0.5, -1e8))) {
    a <- p[1]; z <- p[2]
    v1 <- mhn_var(a, 1, z)
    for (b in 10^c(-200, -150, 150, 200)) {
      expect_close(mhn_var(a, b, z * sqrt(b)) * b, v1, 1e-12,
                   info = sprintf("alpha=%g z=%g beta=%g", a, z, b))
    }
  }
})

test_that("skewness, kurtosis and variance reach the Gamma limit at an extreme negative tilt", {
  # As z -> -Inf the law tends to Gamma(alpha, rate = |z|) (substitute
  # y = |z| x: the factor exp(-y^2 / z^2) tends to 1), with skewness
  # 2 / sqrt(alpha), excess kurtosis 6 / alpha and variance alpha / z^2.  At
  # z = -1e70 the corrections are of order (alpha + 1) / z^2, about 1e-139.
  # The scale of the law is then 1e-70 and its fourth central moment 1e-280,
  # which, taken over a window that narrow, falls below the smallest double;
  # the quadrature measures offsets in units of its window, so no such power
  # appears.  1e-8 is the working tolerance.
  for (a in c(0.5, 2, 10)) {
    info <- sprintf("alpha=%g", a)
    expect_equal(mhn_skewness(a, 1, -1e70), 2 / sqrt(a), tolerance = 1e-8,
                 info = info)
    expect_equal(mhn_kurtosis(a, 1, -1e70), 6 / a, tolerance = 1e-8, info = info)
    expect_equal(mhn_var(a, 1, -1e70) * 1e140, a, tolerance = 1e-8, info = info)
  }
})

# ============================================================
# 15. Strong positive tilt: exact values at alpha = 2
# ============================================================
#
# For alpha = 2 and beta = 1 write J_n = int_0^Inf x^n exp(-x^2 + g x) dx.
# Integrating the derivative of x^n exp(-x^2 + g x) over (0, Inf) gives
#   2 J_(n+1) = g J_n + n J_(n-1) + [n = 0],
# and J_0 = sqrt(pi) exp(g^2 / 4) Phi(g / sqrt(2)), in agreement with the
# closed forms of Sun et al. (2023, Supplementary Lemma 9(c)).  Since
# E X^k = J_(k+1) / J_1, eliminating J gives
#   E X = g/2 + J_0 / (g J_0 + 1),   Var X = 1 + g E X / 2 - (E X)^2
# (the second is Lemma 2c with alpha = 2),
# and as J_0 grows like exp(g^2 / 4), E X = g/2 + 1/g and Var X = 1/2 - 1/g^2
# up to terms of order exp(-g^2 / 4), below 1e-4000 here.  The recurrence
# evaluated in Rmpfr at 400 bits agrees with these to every digit.
#
# At these tilts both logarithms in the Lemma 2a ratio for the mean are about
# g^2 / 4 while the mean is about g / 2, and the Lemma 2c form of the variance
# multiplies the mean's relative error by about g^2 / 2 (see var_impl in
# src/mhn_moments.cpp).  A mean 6e-9 too low at g = 900 would put the
# variance at 0.5012 there, above the Lemma 4c bound of 1/2.

test_that("mean and variance at alpha = 2 match the exact values for a large tilt", {
  for (g in c(200, 600, 900)) {
    info <- sprintf("gamma=%g", g)
    # The mean is held as the peak plus a quadrature offset of order 1/g^3,
    # so it is exact to rounding; 1e-12 relative is loose for it.
    expect_equal(mhn_mean(2, 1, g), g / 2 + 1 / g, tolerance = 1e-12, info = info)
    # The variance is a quadrature asked for sqrt(eps) relative.
    expect_equal(mhn_var(2, 1, g), 1 / 2 - 1 / g^2, tolerance = 1e-8, info = info)
  }
})

test_that("variance under a strong positive tilt matches references for alpha != 2", {
  # References: the trapezoid rule on the peak-centred kernel in Rmpfr at 200
  # bits over +-40 standard deviations; steps of 1/60 and 1/90 of a standard
  # deviation agree to 6e-17, and the same rule reproduces the exact alpha = 2
  # values of the test above to 17 digits.  (10, 0.001, 20) is z = 632.46,
  # checked in the original units.  Tolerances as in the test above.
  expect_equal(mhn_var(3, 1, 700), 0.49999591841732577, tolerance = 1e-8)
  expect_equal(mhn_mean(3, 1, 700), 350.00285713119536, tolerance = 1e-12)
  expect_equal(mhn_var(1.5, 1, 800), 0.49999921874633785, tolerance = 1e-8)
  expect_equal(mhn_mean(1.5, 1, 800), 400.00062500097658, tolerance = 1e-12)
  expect_equal(mhn_var(10, 0.001, 20), 499.97750269966257, tolerance = 1e-8)
  expect_equal(mhn_mean(10, 0.001, 20), 316.24199569714136 / sqrt(0.001),
               tolerance = 1e-12)
})

# ============================================================
# 16. Large alpha: the variance does not inherit the mean's rounding
# ============================================================
#
# At gamma = 0, X^2 ~ Gamma(alpha/2, rate = beta) (Sun et al. 2023,
# Lemma 6a), so at beta = 1 E X = Gamma((alpha+1)/2) / Gamma(alpha/2) and
# Var X = alpha/2 - (E X)^2; the references are these evaluated in Rmpfr at
# 600 bits (900 bits agree to every digit shown).  In the Lemma 2c form
# Var = alpha/2 - mu^2 the mean's relative error is multiplied by about
# 4 alpha, and the mean formed as exp(lgamma((alpha+1)/2) - lgamma(alpha/2))
# carries a rounding of about eps alpha log(alpha) from its two logarithms:
# 5e-10 at alpha = 1e6, which puts the variance at 0.2505 against 0.24999994.
# The package forms it as sqrt(pi) / B(alpha/2, 1/2) with R's lbeta instead,
# which is good to a few units in the last place at any alpha.

test_that("mean and variance at gamma = 0 match exact values for large alpha", {
  ref <- list(
    list(alpha = 1e4, mean = 70.708910373801629, var = 0.24999374968751953),
    list(alpha = 1e6, mean = 707.10660440987431, var = 0.24999993749996874),
    list(alpha = 1e8, mean = 7071.0677941878057, var = 0.249999999375)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g", r$alpha)
    # The mean is exact to rounding (7e-16 here); 1e-12 relative is loose for
    # it and would still catch the lgamma difference at alpha = 1e6.  The
    # variance comes from the Lemma 2c form at alpha = 1e4 (2e-11 adrift, the
    # mean's rounding times 4 alpha) and from the quadrature at 1e6 and 1e8;
    # 1e-8 relative is the working tolerance.
    expect_equal(mhn_mean(r$alpha, 1, 0), r$mean, tolerance = 1e-12, info = info)
    expect_equal(mhn_var(r$alpha, 1, 0), r$var, tolerance = 1e-8, info = info)
  }
  # With a tilt the same holds, and so must the bound Var <= 1/(2 beta) of
  # Lemma 4c.  Reference: the trapezoid rule of section 15.
  expect_equal(mhn_var(1e9, 1, -1), 0.2499972048525276, tolerance = 1e-8)
  expect_equal(mhn_mean(1e9, 1, -1), 22360.429770805331, tolerance = 1e-12)
})

test_that("the variance at a huge alpha is about 1/4 and within the Lemma 4c bound", {
  # At alpha = 1e16 the two logarithms in the Lemma 2a ratio are about 1.8e17,
  # and their difference, about 18, is lost to rounding.  An error estimate
  # for the Lemma 2c form that is linear in the mean's error, and formed from
  # that mean, then accepts alpha/2 - 1 = 5e15 as the variance, which Lemma 4c
  # bounds by 1/(2 beta).
  #
  # For a large alpha, Var X = 1/4 - 1/(16 alpha) + O(alpha^-2) at z = 0 (the
  # formula of the test above, expanded by Stirling's series), and
  # d Var / dz = E[(X - E X)^3] (at beta = 1 the tilt z generates the
  # cumulants), about 1/(8 sqrt(2 alpha)) here.  So Var X =
  # 1/4 + z / (8 sqrt(2 alpha)) to within about 1e-14 at these points; 1e-8
  # relative is the working tolerance.  The last check repeats the first at
  # beta = 1e-10, where alpha/2 - 1 would make beta Var X 5e99.
  for (p in list(c(1e16, 0), c(1e16, 12), c(1e20, -1), c(1e100, 0))) {
    a <- p[1]; z <- p[2]
    info <- sprintf("alpha=%g z=%g", a, z)
    expected <- 1 / 4 + z / (8 * sqrt(2 * a))
    expect_equal(mhn_var(a, 1, z), expected, tolerance = 1e-8, info = info)
  }
  v <- mhn_var(1e100, 1e-10, 0)
  expect_true(v <= 1 / (2 * 1e-10))
  expect_equal(v * 1e-10, 1 / 4, tolerance = 1e-8)
})

# ============================================================
# 17. alpha < 1: the unbounded factor x^(alpha-1) at the origin
# ============================================================
#
# References: E X^k = Psi[(alpha+k)/2, z] / Psi[alpha/2, z] (Sun et al. 2023,
# Lemma 2a, at beta = 1), with the Fox-Wright function
#   Psi[a, z] = sum_n Gamma(a + n/2) z^n / n!
# of their Supplementary Eq. (1), summed in Rmpfr at 500 bits and again at
# 700 bits to 3000 terms, with the same digits; the central moments are formed
# from the raw ones at that precision.
#
# With 0 < alpha < 1 and gamma > 0 large enough for an interior peak
# (Lemma 3c), the density falls from +Inf at the origin to a local minimum and
# rises to the peak.  When the peak lies within a few standard deviations of
# the origin, much of the mass sits closer to the origin than an offset from
# the peak can resolve -- at (0.01, 1, 3) the peak is at 1.01 and
# P(X < 2.2e-16) = 0.62 -- and a rule that treats the stretch between origin
# and peak as smooth gives skewness 1.20 against 3.47.

test_that("skewness and kurtosis for alpha < 1 with an interior peak near the origin", {
  ref <- list(
    list(alpha = 0.05, gamma = 3, skew = 1.3094906882856663, kurt = 0.72701703003463569),
    list(alpha = 0.01, gamma = 3, skew = 3.4664856698969233, kurt = 12.149158551891004),
    list(alpha = 1e-6, gamma = 3, skew = 363.19363730586815, kurt = 147843.37892832942),
    list(alpha = 0.1, gamma = 2.8, skew = 0.98232785308643111, kurt = 0.032238198502032404),
    list(alpha = 0.5, gamma = 5, skew = -0.050262735484446444, kurt = 0.023707109085474697)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g gamma=%g", r$alpha, r$gamma)
    expect_false(is.na(mhn_mode(r$alpha, 1, r$gamma)), info = info)
    # The quadratures are asked for sqrt(eps) relative; they meet the
    # references to 3e-12 or better here.
    expect_equal(mhn_skewness(r$alpha, 1, r$gamma), r$skew, tolerance = 1e-8,
                 info = info)
    expect_equal(mhn_kurtosis(r$alpha, 1, r$gamma), r$kurt, tolerance = 1e-8,
                 info = info)
  }
})

test_that("central moments for alpha < 1 without a peak include the mass at the origin", {
  # Without an interior peak the mass near the origin is just as far out of a
  # quadrature's reach.  The normalising constant accounts for it in closed
  # form, and so must every central moment, in which that mass carries the
  # factor (0 - mu)^k.  References and tolerances as in the test above.
  ref <- list(
    list(alpha = 1e-3, gamma = 1, var = 0.0013608247169128297,
         skew = 30.667755126738005, kurt = 1145.6481270518598),
    list(alpha = 1e-3, gamma = -1, var = 0.00022698000060549994,
         skew = 46.491593352402091, kurt = 2857.3624434244525),
    list(alpha = 0.02, gamma = -1, var = 0.0044650671375853325,
         skew = 10.350690700972178, kurt = 141.1533652334362),
    list(alpha = 1e-6, gamma = -1, var = 2.2717912013720304e-07,
         skew = 1470.5316858855508, kurt = 2859183.2561681429),
    list(alpha = 0.02, gamma = 2.79, var = 0.25283202661152343,
         skew = 2.7652473950957774, kurt = 7.384606396839323),
    list(alpha = 1e-3, gamma = 2.8, var = 0.017385201408661392,
         skew = 13.08714784289549, kurt = 192.0628600861925)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g gamma=%g", r$alpha, r$gamma)
    expect_true(is.na(mhn_mode(r$alpha, 1, r$gamma)), info = info)
    expect_equal(mhn_var(r$alpha, 1, r$gamma), r$var, tolerance = 1e-8, info = info)
    expect_equal(mhn_skewness(r$alpha, 1, r$gamma), r$skew, tolerance = 1e-8,
                 info = info)
    expect_equal(mhn_kurtosis(r$alpha, 1, r$gamma), r$kurt, tolerance = 1e-8,
                 info = info)
  }
})

test_that("skewness and kurtosis for alpha just below 1 match references", {
  # As alpha rises to 1 the local minimum x_min = (1 - alpha) / (2 m) falls
  # far below the peak m, and near it the factor (x/m)^(alpha-1) bends
  # sharply: its log-derivative (alpha-1)/x is -2 m at x_min against
  # (alpha-1)/m at the peak.  The kernel there weighs only about exp(-m^2) of
  # the peak, 1e-10 at z = 9.6, but a rule that takes the stretch from x_min
  # to the peak as smooth misses the values below by 1.2e-8 to 3e-5
  # relative.  References as in the tests above, at 300 and 500 bits (400 and
  # 700 at z = 9.6) with the same digits, for the alpha actually stored.  The
  # package meets them to 5e-10 relative or better (the value at
  # (0.999, 9.6) is only 3.6e-6 in size), so 1e-8 is the working tolerance
  # with room to spare.
  ref <- list(
    list(fn = mhn_skewness, alpha = 0.9999, z = 5, value = 0.0088939260080368593),
    list(fn = mhn_kurtosis, alpha = 0.9999, z = 5, value = -0.026059458442726859),
    list(fn = mhn_kurtosis, alpha = 0.99999, z = 4, value = -0.11094079832445684),
    list(fn = mhn_skewness, alpha = 0.99999, z = 3, value = 0.1868376844743049),
    list(fn = mhn_kurtosis, alpha = 0.999, z = 6.5, value = -0.00083457328827481853),
    list(fn = mhn_skewness, alpha = 0.99, z = 7, value = -0.00018206435963738415),
    list(fn = mhn_kurtosis, alpha = 0.99, z = 7, value = -2.8396531283964692e-05),
    list(fn = mhn_kurtosis, alpha = 0.999, z = 9.6, value = 3.6214377465259817e-06),
    list(fn = mhn_kurtosis, alpha = 0.99, z = 9.6, value = 3.6389998090745096e-05)
  )
  for (r in ref) {
    expect_equal(r$fn(r$alpha, 1, r$z), r$value, tolerance = 1e-8,
                 info = sprintf("alpha=%g z=%g", r$alpha, r$z))
  }
})

test_that("skewness and kurtosis near alpha = 1 at a small tilt are accurate", {
  # Just past the tilt z = sqrt(8 (1 - alpha)) at which the peak appears
  # (Sun et al. 2023, Lemma 3c) the peak m lies close to the origin, and to
  # its right the factor (1 + s/m)^(alpha-1) bends on the scale m, far inside
  # a window that reaches sqrt(45) from the peak; a rule that takes that
  # stretch as smooth misses the values below by up to 1.7e-10.  References:
  # E X^j = Psi[(alpha + j)/2, z] / Psi[alpha/2, z] (Sun et al. 2023,
  # Lemma 2a), with Psi summed as its power series in z in Rmpfr at 320 bits,
  # for the alpha actually stored.  The package meets them to 3e-15 here, and
  # 1e-13 allows for that.
  ref <- list(
    list(alpha = 0.999, z = 0.0939, skew = 0.97102442477146389,
         kurt = 0.79965218434126983),
    list(alpha = 0.9999, z = 0.0297, skew = 0.98750919173935392,
         kurt = 0.84676851582682966),
    list(alpha = 0.99999, z = 0.00939, skew = 0.99280797634400320,
         kurt = 0.86204942895990089),
    list(alpha = 0.999, z = 0.134, skew = 0.96027369614789226,
         kurt = 0.76971162520618641),
    list(alpha = 0.9999, z = 0.0424, skew = 0.98414994412398094,
         kurt = 0.83715889272985267),
    list(alpha = 0.99999, z = 0.0134, skew = 0.99175183053083016,
         kurt = 0.85900265167789058),
    list(alpha = 0.999, z = 0.268, skew = 0.92385617484555405,
         kurt = 0.67147876532284900),
    list(alpha = 0.9999, z = 0.0849, skew = 0.97285597312988592,
         kurt = 0.80515975298938292),
    list(alpha = 0.99999, z = 0.0268, skew = 0.98821727209926380,
         kurt = 0.84883647746109636)
  )
  for (r in ref) {
    info <- sprintf("alpha=%g z=%g", r$alpha, r$z)
    expect_lt(abs(mhn_skewness(r$alpha, 1, r$z) - r$skew), 1e-13, label = info)
    expect_lt(abs(mhn_kurtosis(r$alpha, 1, r$z) - r$kurt), 1e-13, label = info)
  }
})

# ============================================================
# 18. Tilts below sqrt(eps) in size
# ============================================================

test_that("mean and variance at a tilt below sqrt(eps) are not those of z = 0", {
  # The normalising constant treats |z| < sqrt(eps) = 1.5e-8 as z = 0, which
  # moves the Lemma 2a ratio by about z Var X / E X relative, and the Lemma 2c
  # variance, which magnifies that, by about 2 E X |z| relative: 1.8e-8 at
  # (2, 1e-8) and 4e-7 at (400, -1.4e-8).  References: the series of
  # section 17 in Rmpfr at 400 bits (60 terms, which at |z| < 1.5e-8 leave
  # nothing to add).  The package meets them to 2e-15; the tolerances leave
  # that ample room and still catch the substitution.
  expect_equal(mhn_var(2, 1, 1e-8), 0.21460183722996781, tolerance = 1e-10)
  expect_equal(mhn_var(100, 1, 1e-8), 0.24937189506634264, tolerance = 1e-10)
  expect_equal(mhn_var(400, 1, -1.4e-8), 0.24984355493251761, tolerance = 1e-10)
  # The mean moves by only 2.4e-9 relative at (2, 1e-8).
  expect_equal(mhn_mean(2, 1, 1e-8), 0.88622692759877641, tolerance = 1e-12)
})

# ============================================================
# 19. Strong negative tilt near alpha = 1
# ============================================================
#
# Under a negative tilt the Lemma 2c form of the variance magnifies the
# relative error of the mean by about z^2 / 2, and for alpha near 1 the
# normalising constants behind the Lemma 2a ratio can leave that ratio 1e-11
# adrift, which the form turns into 1e-6 at z = -800.  Below z = -30 the
# variance comes from the central-moment quadrature instead.  References: the
# expansion of section 13, sum_j (-1/z^2)^j / j! Gamma(alpha + k + 2 j), in
# Rmpfr at 300 bits, truncated where its terms stop falling; 700 bits agree
# to 1e-89.

test_that("variance under a strong negative tilt near alpha = 1 matches references", {
  ref <- list(
    list(alpha = 0.88, z = -200, var = 2.1993798459696343e-05),
    list(alpha = 0.9, z = -800, var = 1.4062249517980760e-06),
    list(alpha = 0.95, z = -400, var = 5.9370658646293675e-06),
    list(alpha = 1 + 1e-6, z = -700, var = 2.0407683897965300e-06)
  )
  for (r in ref) {
    # The quadrature is asked for sqrt(eps) relative and meets these to 4e-15;
    # the Lemma 2c form misses them by 5e-8 to 1e-6.
    expect_equal(mhn_var(r$alpha, 1, r$z), r$var, tolerance = 1e-8,
                 info = sprintf("alpha=%g z=%g", r$alpha, r$z))
  }
})

# ============================================================
# 20. 1 < alpha < 3: the corner of x^(alpha-1) at the origin
# ============================================================
#
# For 1 < alpha < 2 the factor x^(alpha-1) has an unbounded derivative at the
# origin, and for 2 < alpha < 3 an unbounded second derivative.  When the
# peak is close enough to the origin for the kernel there to matter, a rule
# on the peak-centred window alone loses accuracy at that corner or bisects
# towards it at length; near alpha = 1 the power factor also keeps changing
# far to the right of a peak squeezed against the origin.  References: the
# series of section 17 (z >= 0) and the expansion of section 19 (z < 0), in
# Rmpfr at 600 and 900 bits (300 and 700 for the expansion), which agree to
# more than 80 digits.

test_that("moments for 1 < alpha < 3 with the peak near the origin match references", {
  ref <- list(
    # The kernel by the origin weighs about exp((alpha - 1) - m^2) = 1e-7 of
    # the peak at (1.5, 8), and 1e-7 at (1.1, 8).
    list(fn = mhn_skewness, alpha = 1.5, z = 8, value = 0.0063056144465324808),
    list(fn = mhn_kurtosis, alpha = 1.5, z = 8, value = -0.0036568859031795040),
    list(fn = mhn_skewness, alpha = 1.1, z = 8, value = 0.0013585242144882114),
    list(fn = mhn_kurtosis, alpha = 1.1, z = 8, value = -0.00084641875798222197),
    list(fn = mhn_kurtosis, alpha = 1 + 1e-6, z = 0, value = 0.86917584050580086),
    list(fn = mhn_var, alpha = 1.001, z = -1000, value = 1.0009879821943760e-06),
    list(fn = mhn_skewness, alpha = 1.001, z = -1000, value = 1.9989887496501146),
    list(fn = mhn_kurtosis, alpha = 1.001, z = -1000, value = 5.9939100447893413),
    list(fn = mhn_var, alpha = 2.5, z = -1000, value = 2.4999475013999565e-06),
    list(fn = mhn_skewness, alpha = 2.5, z = -1000, value = 1.2648977829361328)
  )
  for (r in ref) {
    # The quadratures are asked for sqrt(eps) relative and meet these to
    # 3e-9 (the excess kurtosis at (1.1, 8), whose size is 8e-4) or better;
    # a rule on the window alone misses the first five by 1e-8 to 9e-7.
    expect_equal(r$fn(r$alpha, 1, r$z), r$value, tolerance = 1e-8,
                 info = sprintf("alpha=%g z=%g", r$alpha, r$z))
  }
})

# ============================================================
# 21. Tilts near the overflow of z^2
# ============================================================

test_that("moments past the overflow of z^2 under a positive tilt are those of the normal limit", {
  # From z = 1.34e154 on z^2 overflows, and from about 2.7e154 so do the
  # logarithms of the Fox-Wright values, about z^2/4, so the ratio of Lemma 2a
  # cannot be formed; the moments come from the quadrature about the peak,
  # whose position must then be found without squaring z.  The law is normal
  # there to far below double precision: its log density about the peak m is
  # -(x - m)^2 (2 + (alpha - 1)/m^2)/2 up to a cubic term of relative size
  # (alpha - 1)/m^3, below 1e-300 here.  So the mean is m and the variance
  # 1/(2 + (alpha - 1)/m^2), each to within about 1e-300 relative, and the
  # skewness and excess kurtosis are 0 to within 1e-150.  m is the positive
  # root of 2 x^2 - z x - (alpha - 1) = 0, written as
  # (z/4)(1 + sqrt(1 + 8 (alpha - 1)/z^2)) with the ratio formed first, so
  # that nothing overflows; at these tilts it is z/2 to double precision.
  # 1e-8 is the working tolerance, absolute for the two shape measures.
  for (a in c(1e-3, 0.5, 2, 1e3)) {
    for (z in c(1.35e154, 1e155, 1e200, 1e300, .Machine$double.xmax)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      m <- (z / 4) * (1 + sqrt(1 + (8 * (a - 1) / z) / z))
      expect_equal(mhn_mean(a, 1, z), m, tolerance = 1e-8, info = info)
      expect_equal(mhn_var(a, 1, z), 1 / (2 + ((a - 1) / m) / m),
                   tolerance = 1e-8, info = info)
      expect_lt(abs(mhn_skewness(a, 1, z)), 1e-8, label = info)
      expect_lt(abs(mhn_kurtosis(a, 1, z)), 1e-8, label = info)
    }
  }
})

test_that("moments under a negative tilt match the Gamma limit up to the largest double", {
  # For alpha <= 1 and z < 0 the kernel has no interior peak, and the
  # quadrature runs from the origin to where exp(z x - x^2) has fallen by 45
  # nats, the positive root of x^2 - z x - 45 = 0, about 45/|z|; for
  # alpha > 1 it is centred on the peak, about (alpha - 1)/|z|.  Both roots
  # must be formed without squaring z, which overflows from |z| = 1.34e154
  # on: a root of 0 would cut the Gamma tail off (at alpha = 1e-3 the mean
  # would come out 0.04 of its value), or leave no peak to centre on.  At
  # these tilts the law is Gamma(alpha, rate = |z|), the limit of the test
  # at z = -1e70 in section 14, up to relative corrections of order
  # (alpha + 1) / z^2 (2 (alpha + 1) / z^2 in the mean, 8 (alpha + 1) / z^2
  # in the excess kurtosis), below 2e-303 here, so its mean, variance,
  # skewness and excess kurtosis are alpha / |z|, alpha / z^2, 2 / sqrt(alpha)
  # and 6 / alpha.  1e-8 is the working tolerance.  The variance is compared
  # only where alpha / z^2 is a normal double; below 2.2e-308 it is held to
  # the spacing of the subnormal doubles (see "a subnormal variance near
  # z = -1e154 is accurate to the spacing of doubles"), and below 4.9e-324
  # it rounds to 0.
  zs <- c(-6.8e153, -1e154, -1.3e154, -1e155, -1e200, -1e300,
          -.Machine$double.xmax)
  for (a in c(1e-3, 0.5, 1, 2, 3, 1e3)) {
    for (z in zs) {
      info <- sprintf("alpha=%g z=%g", a, z)
      got <- c(mhn_mean(a, 1, z) * abs(z) / a,
               mhn_skewness(a, 1, z) * sqrt(a) / 2, mhn_kurtosis(a, 1, z) * a / 6)
      expect_equal(got, rep(1, 3), tolerance = 1e-8, info = info)
      v <- mhn_var(a, 1, z)
      vref <- (a / abs(z)) / abs(z)
      if (vref >= .Machine$double.xmin) {
        expect_equal(v / vref, 1, tolerance = 1e-8, info = info)
      } else {
        expect_lte(abs(v - vref), 1e-8 * vref + 2 * 4.94e-324, label = info)
      }
    }
  }
})

test_that("moments with a subnormal peak just above alpha = 1 match the Gamma limit", {
  # For alpha just above 1 and z near -1e308 the peak, about (alpha - 1)/|z|,
  # is subnormal and holds only some of its digits (17 bits at
  # (1 + 1e-10, -1e308), where it is 1e-318).  The law is Gamma(alpha,
  # rate |z|) to within (alpha + 1)/z^2, below 1e-300, so the mean,
  # skewness and excess kurtosis are alpha/|z|, 2/sqrt(alpha) and 6/alpha;
  # the variance, alpha/z^2, underflows.  1e-8 is the working tolerance.
  for (a in c(1 + 1e-10, 1 + 1e-6, 1.5)) {
    for (z in c(-1e308, -1.7e308, -.Machine$double.xmax)) {
      info <- sprintf("alpha=%.12g z=%g", a, z)
      got <- c(mhn_mean(a, 1, z) * abs(z) / a,
               mhn_skewness(a, 1, z) * sqrt(a) / 2, mhn_kurtosis(a, 1, z) * a / 6)
      expect_equal(got, rep(1, 3), tolerance = 1e-8, info = info)
    }
  }
})

test_that("a subnormal variance near z = -1e154 is accurate to the spacing of doubles", {
  # Without an interior peak the quadrature holds the variance as h^2 c2,
  # with the unit h about alpha / |z| and c2 about 1 / alpha.  Near
  # z = -1e154 h^2 falls below the smallest normal double, 2.2e-308, and for
  # small alpha below the smallest subnormal, while alpha / z^2 is still
  # representable: h^2 formed first gives 0 at (1e-8, -1e154) and
  # (1e-10, -1e152), and 3e-4 too much at (1e-6, -1.2e154).  The law is
  # Gamma(alpha, rate = |z|) to within 6 (alpha + 1) / z^2, at most 6e-304,
  # in the variance.  These variances are subnormal, held only to the
  # spacing of doubles there, 4.9e-324: that is 4.9e-8 of 1e-316, 4.9e-10 of
  # 1e-314 and 7.1e-10 of 6.9e-315, and each tolerance is about twice it.
  expect_equal(mhn_var(1e-8, 1, -1e154) * 1e308 / 1e-8, 1, tolerance = 1e-7)
  expect_equal(mhn_var(1e-10, 1, -1e152) * 1e304 / 1e-10, 1, tolerance = 1e-9)
  expect_equal(mhn_var(1e-6, 1, -1.2e154) * 1.44e308 / 1e-6, 1,
               tolerance = 1.5e-9)
})

test_that("moments where gamma / sqrt(beta) overflows or underflows at beta = 1 are those of X", {
  # The moments are computed at beta = 1 and z = gamma / sqrt(beta), then
  # rescaled.  z itself overflows once |gamma| exceeds the largest double
  # times sqrt(beta).  Under such a negative tilt the law of X is
  # Gamma(alpha, rate |gamma|) to within a relative of order (alpha + 1) / z^2,
  # below 1e-307; under such a positive one it is normal about its mode m, with
  # variance 1 / (2 beta + (alpha - 1) / m^2), to within 1 / (beta m^2),
  # below 1e-600.  At (2, 0.9, DBL_MAX) the mode is DBL_MAX / 1.8 to a
  # relative 1e-616.  And where alpha / |z|, the mean at beta = 1, or
  # alpha / z^2, its variance, falls below the normal range, the moments of
  # X, 1 / sqrt(beta) and 1 / beta times as large, can still lie inside it:
  # at (0.5, 1e-300, -1e50), z = -1e200, the variance is 5e-101, and at
  # (1e-10, 1e-100, -1.7e258), z = -1.7e308, the mean is 1e-10 / 1.7e258.
  # Both are Gamma limits, as above.  1e-8 is the working tolerance.
  expect_close(mhn_mean(2, 1e-320, -1e150), 2e-150, 1e-8)
  expect_close(mhn_var(2, 1e-320, -1e150), 2e-300, 1e-8)
  expect_equal(mhn_skewness(2, 1e-320, -1e150), 2 / sqrt(2), tolerance = 1e-8)
  expect_equal(mhn_kurtosis(2, 1e-320, -1e150), 3, tolerance = 1e-8)
  xmax <- .Machine$double.xmax
  expect_equal(mhn_mean(2, 0.9, xmax), xmax / 1.8, tolerance = 1e-8)
  expect_equal(mhn_var(2, 0.9, xmax), 1 / 1.8, tolerance = 1e-8)
  expect_equal(mhn_skewness(2, 0.9, xmax), 0)
  expect_equal(mhn_kurtosis(2, 0.9, xmax), 0)
  expect_close(mhn_var(0.5, 1e-300, -1e50), 5e-101, 1e-8)
  expect_close(mhn_mean(1e-10, 1e-100, -1.7e258), 1e-10 / 1.7e258, 1e-8)
})

test_that("moments at a tiny shape stay in range", {
  # Under a strong negative tilt the law is Gamma(alpha, rate |z|) to within
  # relative corrections of order (alpha + 1) / z^2, below 1e-19 here, so
  # the skewness and excess kurtosis are 2 / sqrt(alpha) and 6 / alpha, up
  # to 6e300.  The range of the quadrature then runs out to about 45 / |z|,
  # 4.5e101 times the scale alpha / |z| at (1e-100, -1e10), and its fourth
  # power leaves the double range; at alpha = 1e-300, z = -1e10 the scale
  # itself is subnormal.  The variance alpha / z^2 underflows to 0 at
  # (1e-300, -1e100).  1e-8 is the working tolerance.
  for (a in c(1e-100, 1e-200, 1e-300)) {
    for (z in c(-1e10, -1e100)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      expect_equal(mhn_skewness(a, 1, z) * sqrt(a) / 2, 1, tolerance = 1e-8,
                   info = info)
      expect_equal(mhn_kurtosis(a, 1, z) * a / 6, 1, tolerance = 1e-8,
                   info = info)
    }
  }
  expect_equal(mhn_var(1e-100, 1, -1e10) * 1e20 / 1e-100, 1, tolerance = 1e-8)
  expect_identical(mhn_var(1e-300, 1, -1e100), 0)
  # Under a positive tilt a tiny alpha puts almost all the mass at the
  # origin.  E[X^k] = Psi(alpha + k) / Psi(alpha) with Psi(a) =
  # int_0^inf x^(a-1) exp(-x^2 + z x) dx, and Psi(alpha) is
  # (1 + O(alpha e^(z^2/4))) / alpha, e^(z^2/4) being 7e10 at z = 10, while
  # Psi(alpha + k) = I_k (1 + O(alpha)), I_k = Psi(k).  So the skewness is
  # I_3 / (I_2^(3/2) sqrt(alpha)) and the excess kurtosis
  # I_4 / (alpha I_2^2) - 3 to within a relative 1e-100 here.  The variance,
  # about alpha I_2, is squared on the way to the kurtosis and passes below
  # the smallest double from alpha = 1e-155 on.
  I <- function(k, z) {
    integrate(function(x) x^(k - 1) * exp(-x^2 + z * x), 0, Inf,
              rel.tol = 1e-12)$value
  }
  for (z in c(2, 5, 10)) {
    i2 <- I(2, z)
    i3 <- I(3, z)
    i4 <- I(4, z)
    for (a in c(1e-120, 1e-200, 1e-300)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      expect_close(mhn_skewness(a, 1, z) * sqrt(a), i3 / i2^1.5, 1e-8,
                   info = info)
      expect_close(mhn_kurtosis(a, 1, z) * a, i4 / i2^2, 1e-8, info = info)
    }
  }
})

test_that("mhn_mean holds far below the peak (z > 12, tiny alpha)", {
  # For z > 12 the mean comes from a quadrature that holds it as the peak plus
  # an offset.  Under a shape so small that almost all the mass lies by the
  # origin, the mean is far below the peak, and that sum keeps the rounding of
  # the peak: 1.1e-4 of the mean at alpha = 1e-27, z = 12.1.  The ratio of
  # Lemma 2a (Sun et al. 2023) forms no such sum.  As alpha -> 0,
  # Psi[alpha/2, z] is 2 / alpha and Psi[(alpha+1)/2, z] is 2 A, with
  # A = int_0^Inf exp(-x^2 + z x) dx = sqrt(pi) exp(z^2/4) pnorm(z / sqrt(2)),
  # so E X is alpha A, up to a relative alpha A / (z / 2), below 1e-11 here.
  A <- function(z) sqrt(pi) * exp(z^2 / 4) * pnorm(z / sqrt(2))
  for (p in list(c(1e-27, 12.1), c(1e-28, 12.5), c(1e-29, 13), c(1e-35, 15),
                 c(1e-54, 20), c(1e-184, 40))) {
    info <- sprintf("alpha=%g z=%g", p[1], p[2])
    expect_equal(mhn_mean(p[1], 1, p[2]) / (p[1] * A(p[2])), 1,
                 tolerance = 1e-9, info = info)
  }
  # The same law at another scale: sqrt(beta) X ~ MHN(alpha, 1, gamma /
  # sqrt(beta)).
  expect_equal(mhn_mean(1e-27, 4, 24.2) * 2 / (1e-27 * A(12.1)), 1,
               tolerance = 1e-9)
})
