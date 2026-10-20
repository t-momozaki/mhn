# Tests for dmhn() -- MHN density function.
#
# Covers: integrate-to-1, non-negativity, special values, log = TRUE,
# vectorization, and NA propagation.

# ============================================================
# 1. integrate(dmhn, 0, Inf) ~= 1 over a parameter grid
# ============================================================

test_that("density integrates to 1 for diverse parameter combinations", {
  test_params <- expand.grid(
    alpha = c(0.5, 1, 1.5, 2, 3, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-5, -1, 0, 1, 5)
  )
  # 90 parameter combinations: alpha below, at and above 1, three values of
  # beta, and gamma negative, zero and positive.
  for (i in seq_len(nrow(test_params))) {
    a <- test_params$alpha[i]
    b <- test_params$beta[i]
    g <- test_params$gamma[i]
    result <- integrate(
      function(x) dmhn(x, alpha = a, beta = b, gamma = g),
      lower = 0, upper = Inf,
      rel.tol = 1e-12, subdivisions = 200
    )
    # 4e-13 off at most, as measured.
    expect_equal(result$value, 1,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 2. Non-negativity
# ============================================================

test_that("density is finite and non-negative for all x", {
  # The one infinite value allowed is the pole of x^(alpha - 1) at x = 0 for
  # alpha < 1; NaN, -Inf and Inf anywhere else fail.  The tilts reach -1000,
  # where the normalising constant and the kernel are each far from one.
  x <- c(-2, -1, -0.01, 0, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50)
  for (a in c(0.5, 1, 2)) {
    for (b in c(1e-4, 1)) {
      for (g in c(-1000, -100, -2, 0, 2)) {
        d <- dmhn(x, alpha = a, beta = b, gamma = g)
        expect_true(all((is.finite(d) & d >= 0) |
                          (x == 0 & a < 1 & d == Inf)),
                    label = sprintf("alpha=%g, beta=%g, gamma=%g", a, b, g))
      }
    }
  }
})

# ============================================================
# 3. Handling of x < 0
# ============================================================

test_that("density is 0 for x < 0", {
  expect_equal(dmhn(-1, alpha = 2, beta = 1, gamma = 0), 0)
  expect_equal(dmhn(-0.001, alpha = 2, beta = 1, gamma = 0), 0)
  expect_equal(dmhn(c(-2, -1), alpha = 2, beta = 1, gamma = 1), c(0, 0))
})

test_that("log-density is -Inf for x < 0", {
  expect_equal(dmhn(-1, alpha = 2, beta = 1, gamma = 0, log = TRUE), -Inf)
  expect_equal(dmhn(c(-2, -1), alpha = 2, beta = 1, gamma = 1, log = TRUE),
               c(-Inf, -Inf))
})

# ============================================================
# 4. Boundary x = 0 (case split by alpha)
# ============================================================

test_that("density at x=0: alpha > 1 gives 0", {
  expect_equal(dmhn(0, alpha = 2, beta = 1, gamma = 1), 0)
  expect_equal(dmhn(0, alpha = 5, beta = 2, gamma = -1), 0)
})

test_that("density at x=0: alpha < 1 gives Inf", {
  expect_equal(dmhn(0, alpha = 0.5, beta = 1, gamma = 1), Inf)
  expect_equal(dmhn(0, alpha = 0.3, beta = 2, gamma = -1), Inf)
})

test_that("density at x=0: alpha=1 gives finite positive value", {
  d <- dmhn(0, alpha = 1, beta = 1, gamma = 0)
  expect_true(is.finite(d) && d > 0)
  expect_equal(d, 2 * sqrt(1 / pi), tolerance = 1e-12)
})

# ============================================================
# 5. log = TRUE consistency
# ============================================================

test_that("log=TRUE is consistent with log(dmhn(...))", {
  x <- c(0.1, 0.5, 1, 2, 5)
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    log_d <- dmhn(x, alpha = a, beta = b, gamma = g, log = TRUE)
    d <- dmhn(x, alpha = a, beta = b, gamma = g)
    # Compare only at points where density is strictly positive.
    pos <- d > 0 & is.finite(d)
    if (any(pos)) {
      expect_equal(log_d[pos], log(d[pos]),
                   tolerance = 1e-12,
                   label = sprintf("a=%.1f, b=%.1f, g=%.1f", a, b, g))
    }
  }
})

# ============================================================
# 6. Vectorization
# ============================================================

test_that("dmhn is vectorized over x", {
  x <- c(0.5, 1, 2)
  d <- dmhn(x, alpha = 2, beta = 1, gamma = 1)
  expect_length(d, 3)
  # Each element matches the scalar call.
  for (j in seq_along(x)) {
    expect_equal(d[j], dmhn(x[j], alpha = 2, beta = 1, gamma = 1),
                 tolerance = 1e-15)
  }
})

# ============================================================
# 7. NA propagation
# ============================================================

test_that("dmhn returns NA for NA input and NaN for NaN input", {
  # is.na() is TRUE for NaN as well, so the two have to be told apart: as in
  # base R's d/p/q family, NA in gives NA out and NaN in gives NaN out, on the
  # general path and on each special case (gamma = 0, alpha = 1).
  result <- dmhn(c(1, NA, 2), alpha = 2, beta = 1, gamma = 0)
  expect_equal(length(result), 3)
  expect_true(is.na(result[2]) && !is.nan(result[2]))
  expect_false(is.na(result[1]))
  expect_false(is.na(result[3]))

  for (p in list(c(2, 1, 1), c(2, 1, 0), c(1, 1, 1))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    na  <- dmhn(NA_real_, p[1], p[2], p[3])
    nan <- dmhn(NaN, p[1], p[2], p[3])
    expect_true(is.na(na) && !is.nan(na), info = info)
    expect_true(is.nan(nan), info = info)
  }
  expect_true(is.nan(dmhn(NaN, 2, 1, 1, log = TRUE)))
})

# ============================================================
# 8. Parameter validation
# ============================================================

test_that("dmhn gives NaN with a warning for parameters outside the parameter space", {
  # As dgamma does; test-input-validation.R tests the rule in full.
  for (args in list(list(alpha = -1), list(alpha = 0), list(beta = -1),
                    list(beta = 0))) {
    expect_warning(out <- do.call(dmhn, c(list(1), args)), "NaNs produced")
    expect_identical(out, NaN)
  }
})

# ============================================================
# 9. Continuity at special-case boundaries
# ============================================================

test_that("general case is continuous near special-case boundaries", {
  x <- c(0.5, 1, 2, 3)
  # Near gamma = 0 boundary.
  d_at_zero <- dmhn(x, alpha = 2, beta = 1, gamma = 0)
  d_near_zero <- dmhn(x, alpha = 2, beta = 1, gamma = 1e-6)
  expect_equal(d_at_zero, d_near_zero, tolerance = 1e-4)

  # Near alpha = 1 boundary.
  d_at_one <- dmhn(x, alpha = 1, beta = 1, gamma = 2)
  d_near_one <- dmhn(x, alpha = 1 + 1e-6, beta = 1, gamma = 2)
  expect_equal(d_at_one, d_near_one, tolerance = 1e-4)
})


# ============================================================
# 10. Parameter vectorization
# ============================================================

test_that("dmhn accepts vector parameters and matches per-element scalar calls", {
  x <- c(0.5, 1, 2, 3)
  alpha <- c(1.5, 2, 3, 4)
  beta <- c(0.5, 1, 1, 2)
  gamma <- c(-1, 0, 1, 2)

  expected <- vapply(seq_along(x), function(i) {
    dmhn(x[i], alpha = alpha[i], beta = beta[i], gamma = gamma[i])
  }, numeric(1))

  result <- dmhn(x, alpha = alpha, beta = beta, gamma = gamma)
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("dmhn recycles shorter parameter vectors (R recycling rules)", {
  x <- 1:6
  # alpha recycles c(1,2) -> c(1,2,1,2,1,2)
  result <- dmhn(x, alpha = c(1, 2), beta = 1, gamma = 0)
  expected <- c(dmhn(1, 1), dmhn(2, 2), dmhn(3, 1),
                dmhn(4, 2), dmhn(5, 1), dmhn(6, 2))
  expect_close(result, expected, 1e-12)
})

test_that("dmhn vector parameters with x recycled (length(x)=1)", {
  # MCMC-style: scalar x, parameters as a chain
  x_obs <- 1.5
  alpha_chain <- c(1.5, 2, 3, 5)
  beta_chain <- c(1, 1, 0.5, 2)
  gamma_chain <- c(-1, 0, 2, -3)

  result <- dmhn(x_obs, alpha = alpha_chain, beta = beta_chain,
                 gamma = gamma_chain)
  expected <- vapply(seq_along(alpha_chain), function(i) {
    dmhn(x_obs, alpha_chain[i], beta_chain[i], gamma_chain[i])
  }, numeric(1))
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("dmhn scalar and length-1-vector parameters give identical results", {
  x <- seq(0.1, 5, length.out = 50)
  expect_equal(dmhn(x, alpha = 2, beta = 1, gamma = 1),
               dmhn(x, alpha = c(2), beta = c(1), gamma = c(1)),
               tolerance = 1e-12)
})

test_that("dmhn vector parameters honor log = TRUE", {
  x <- c(0.5, 1, 2)
  alpha <- c(1.5, 2, 3)
  log_d <- dmhn(x, alpha = alpha, beta = 1, gamma = 0, log = TRUE)
  d <- dmhn(x, alpha = alpha, beta = 1, gamma = 0)
  expect_equal(log_d, log(d), tolerance = 1e-12)
})

test_that("dmhn answers an invalid element of a vector parameter on its own", {
  ok <- dmhn(1, alpha = 1, beta = 1, gamma = 0)
  expect_warning(out <- dmhn(1, alpha = c(1, -1), beta = 1, gamma = 0),
                 "NaNs produced")
  expect_identical(out, c(ok, NaN))
  expect_identical(dmhn(1, alpha = c(1, NA), beta = 1, gamma = 0), c(ok, NA))
  expect_warning(out <- dmhn(1, alpha = 1, beta = c(1, 0), gamma = 0),
                 "NaNs produced")
  expect_identical(out, c(ok, NaN))
  expect_identical(dmhn(1, alpha = 1, beta = 1, gamma = c(0, NA)), c(ok, NA))
  expect_warning(out <- dmhn(1, alpha = 1, beta = 1, gamma = c(0, Inf)),
                 "NaNs produced")
  expect_identical(out, c(ok, NaN))
})

test_that("dmhn handles special-case dispatch element-wise", {
  # Mix gamma=0 (sqrt-Gamma), alpha=1 (truncated-normal), and a general triple.
  x <- c(1, 1, 1)
  alpha <- c(2, 1, 3)
  beta  <- c(1, 1, 1)
  gamma <- c(0, 2, 1)
  result <- dmhn(x, alpha = alpha, beta = beta, gamma = gamma)
  expected <- c(
    dmhn(1, alpha = 2, beta = 1, gamma = 0),
    dmhn(1, alpha = 1, beta = 1, gamma = 2),
    dmhn(1, alpha = 3, beta = 1, gamma = 1)
  )
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("dmhn returns numeric(0) for empty x", {
  expect_equal(dmhn(numeric(0), alpha = c(1, 2)), numeric(0))
})

# ============================================================
# 11. gamma = 0 closed form at the extremes of x and beta
# ============================================================

test_that("the gamma = 0 density stays exact where x^2 underflows", {
  # With t = sqrt(beta) x, f(x) = 2 sqrt(beta) t^(alpha - 1) exp(-t^2) /
  # Gamma(alpha/2) (Sun et al. 2023, Lemma 6a).  At x = 1e-170, x^2 = 1e-340
  # is below the smallest double, so a density that squares x first is
  # evaluated at 0 instead: Inf for alpha < 2, 0 for alpha > 2.  References:
  # the formula above in Rmpfr at 256 bits (exp(-t^2) is 1 to every digit
  # here).  Tolerance: a few units in the last place of log f.
  expect_equal(dmhn(1e-170), 2 / sqrt(pi), tolerance = 1e-14)
  expect_equal(dmhn(1e-170, 3, 1, 0, log = TRUE), -782.06500219978034,
               tolerance = 1e-14)
  expect_equal(dmhn(1e-170, 0.5, 1, 0, log = TRUE), 195.12485756035576,
               tolerance = 1e-14)
  expect_equal(dmhn(1e-170, 2.5, 1, 0, log = TRUE), -586.36777969649984,
               tolerance = 1e-14)
  # x = 1e-161 squares into the subnormal band, 1e-322 = 20 units of
  # 4.9e-324, where x^2 keeps barely two digits.
  expect_equal(dmhn(1e-161), 1.1283791670955126, tolerance = 1e-14)
  # What matters is t = sqrt(beta) x, not x: a large beta lifts t back into
  # range, a small one pushes it further down.
  expect_equal(dmhn(1e-170, 1, 1e100, 0, log = TRUE), 115.25003688733753,
               tolerance = 1e-14)
  expect_equal(dmhn(1e-120, 1, 1e-100, 0, log = TRUE), -115.00847241206704,
               tolerance = 1e-14)
  expect_equal(dmhn(1e-300, 0.5, 1e-100, 0, log = TRUE), 287.22826128011758,
               tolerance = 1e-14)
  # Recycled parameters go through the same per-element evaluation.
  expect_equal(dmhn(c(1e-170, 1e-161), c(1, 1), 1, 0),
               c(2 / sqrt(pi), 1.1283791670955126), tolerance = 1e-14)
  # Above alpha = 20 the density goes through R's dgamma wherever beta x^2 is
  # a normal double, and through the same closed form below that.  Here alpha
  # and beta change from one element to the next, so each element is
  # evaluated with the constants of its own triple.  References: the formula
  # above in Rmpfr at 256 bits.  Tolerance: a few units in the last place.
  got <- dmhn(rep(1e-170, 3), c(50, 50, 3), c(1, 4, 4), 0, log = TRUE)
  ref <- c(-19234.625406857955, -19199.968047829956, -779.98556065810055)
  expect_lt(max(abs(got / ref - 1)), 1e-14)
})

test_that("the gamma = 0 density agrees with its dgamma form across shapes", {
  # f(x) = 2 beta x g(beta x^2), g the Gamma(alpha/2, 1) density, is evaluated
  # in closed form up to alpha = 20 and through R's dgamma above it.  Where x^2
  # is an ordinary number, R's dgamma form is an independent reference.  Each
  # side rounds terms whose sizes add up to S = |log(2 beta x)| +
  # |(alpha/2 - 1) log(beta x^2)| + beta x^2 + |lgamma(alpha/2)|, the closed
  # form cancelling them and dgamma forming beta x^2 in its own order, so the
  # tolerance is 16 eps S: 6e-15 to 2e-11 here, and the largest difference
  # measured is 6% of it.
  for (alpha in c(0.01, 0.5, 1, 2.5, 10, 20, 50, 1000)) {
    for (beta in c(1e-6, 1, 1e6)) {
      x <- sqrt(qgamma(c(1e-6, 0.01, 0.5, 0.99, 1 - 1e-6), alpha / 2,
                       rate = beta))
      x <- x[x^2 > 1e-300]
      ref <- log(2 * beta * x) + dgamma(beta * x^2, alpha / 2, 1, log = TRUE)
      got <- dmhn(x, alpha, beta, 0, log = TRUE)
      y <- beta * x^2
      size <- abs(log(2 * beta * x)) + abs((alpha / 2 - 1) * log(y)) + y +
        abs(lgamma(alpha / 2))
      expect_lt(max(abs(got - ref) / (16 * .Machine$double.eps * size)), 1,
                label = sprintf("alpha=%g beta=%g", alpha, beta))
      # Scalar and recycled parameters give the same doubles.
      expect_identical(dmhn(x, rep(alpha, length(x)), beta, 0, log = TRUE),
                       got)
    }
  }
})

# ============================================================
# 12. Tilts past the overflow of z^2
# ============================================================

test_that("the density at a tilt past the overflow of z^2 is its normal limit", {
  # From z = gamma / sqrt(beta) = 1.34e154 on, z^2 overflows, and the peak of
  # the kernel, which the normalising constant is integrated about, must be
  # found without it.  With beta = 1 the law is normal to far below double
  # precision, with mean m = z/2 + delta, delta = 2 (alpha - 1) / (z +
  # sqrt(z^2 + 8 (alpha - 1))), and variance s2 = 1 / (2 + (alpha - 1)/m^2):
  # the cubic term of its log density is of relative size (alpha - 1)/m^3,
  # below 1e-300.  z/2 is exact, and delta, below 1e-154 here, is far inside
  # the standard deviation, so the log density at x = z/2 is
  # -log(2 pi s2)/2 = -log(pi)/2 to double precision.  At x = z/2 (1 + 1e-10)
  # it is -(x - m)^2 / (2 s2) - log(pi)/2, with x - m = z/2 1e-10 to a
  # relative 1e-16, -2.5e289 at z = 1e155, and -Inf, below -DBL_MAX, at
  # z = 1e200.
  for (a in c(0.5, 2, 10)) {
    for (z in c(1e155, 1e200, 1e300, .Machine$double.xmax)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      expect_equal(dmhn(z / 2, a, 1, z, log = TRUE), -0.5 * log(pi),
                   tolerance = 1e-12, info = info)
      expect_equal(dmhn(z / 2, a, 1, z), 1 / sqrt(pi), tolerance = 1e-12,
                   info = info)
    }
    x <- 1e155 / 2 * (1 + 1e-10)
    d <- x - 1e155 / 2
    expect_equal(dmhn(x, a, 1, 1e155, log = TRUE), -d * d - 0.5 * log(pi),
                 tolerance = 1e-12, info = sprintf("alpha=%g", a))
    expect_identical(dmhn(1e200 / 2 * (1 + 1e-10), a, 1, 1e200, log = TRUE),
                     -Inf)
  }
})

test_that("the density at a negative tilt past the overflow of z^2 is its Gamma limit", {
  # Under a negative tilt the law tends to Gamma(alpha, rate |gamma|), and at
  # these tilts its relative departure from that limit, of order
  # (alpha + 1) / z^2, is below 1e-300.  The quantiles of the limit are
  # divided by |z| in exact binary scaling only for z a power of two, so the
  # reference is taken at the doubles x actually passed: |z| x is formed in
  # double precision, a relative rounding of 1.1e-16 that moves the log
  # density by at most 1.1e-16 (alpha - 1 + |z| x), below 1e-14 here.
  for (a in c(0.5, 2, 10)) {
    for (z in c(-1e155, -1e200, -1e300)) {
      x <- qgamma(c(0.01, 0.5, 0.99), a) / abs(z)
      ref <- dgamma(abs(z) * x, a, log = TRUE) + log(abs(z))
      expect_equal(dmhn(x, a, 1, z, log = TRUE), ref, tolerance = 1e-12,
                   info = sprintf("alpha=%g z=%g", a, z))
    }
  }
})

test_that("the density where gamma / sqrt(beta) overflows is that of its limit law", {
  # z = gamma / sqrt(beta) passes the largest double once |gamma| exceeds it
  # times sqrt(beta), which needs beta < 1, while gamma and
  # log |z| = log |gamma| - log(beta) / 2 are still in range.  Under such a
  # negative tilt the law is Gamma(alpha, rate |gamma|) to within a relative
  # alpha (alpha + 1) / z^2, 0 in double precision, so the log density is
  # dgamma's at |gamma| x plus log |gamma|, the product rounded as in the
  # previous test.  Under such a positive tilt it is normal with variance
  # 1 / (2 beta) about its mode, which at beta = 1/2 is gamma itself to
  # within (alpha - 1) / gamma < 1e-307, far below the spacing of the doubles
  # there, so the log density at x = gamma is log(sqrt(beta / pi)).  The
  # normalising constant is formed from log |z| in both.
  X <- .Machine$double.xmax
  for (a in c(0.5, 2, 10)) {
    for (bg in list(c(1e-100, -X), c(1e-200, -1e250), c(0.5, -X))) {
      b <- bg[1]
      g <- bg[2]
      x <- qgamma(c(0.01, 0.5, 0.99), a) / abs(g)
      expect_equal(dmhn(x, a, b, g, log = TRUE),
                   dgamma(abs(g) * x, a, log = TRUE) + log(abs(g)),
                   tolerance = 1e-12,
                   info = sprintf("alpha=%g beta=%g gamma=%g", a, b, g))
    }
    for (g in c(1.3e308, 1.7e308, X)) {
      expect_equal(dmhn(g, a, 0.5, g, log = TRUE), 0.5 * log(0.5 / pi),
                   tolerance = 1e-12, info = sprintf("alpha=%g gamma=%g", a, g))
    }
  }
})

test_that("the density with a completed square holds where 2 beta overflows", {
  # For z above 1e3, and below the switch to the peak frame (at this beta
  # from z of about 3e4, where z (1 + |log v*|) passes 1e7), the kernel is
  # written about its centre gamma / (2 beta), and 2 beta overflows for beta
  # above about 9e307, where the centre is formed as (gamma / 2) / beta.
  # The density scales as
  # f(x | alpha, beta, gamma) = sqrt(beta) f(sqrt(beta) x | alpha, 1, z), so
  # the reference is taken at beta = 1; forming sqrt(beta) x rounds it by
  # 1.1e-16 relative, which moves the log density by about 1.1e-16 z |k|,
  # below 1e-10 here at k standard deviations from the mode.
  b <- .Machine$double.xmax
  sb <- sqrt(b)
  for (a in c(0.5, 2, 10)) {
    for (z in c(2e3, 2e4)) {
      u <- mhn_mode(a, 1, z) + c(-3, -1, 0, 1, 3) / sqrt(2)
      expect_lt(max(abs(dmhn(u / sb, a, b, z * sb, log = TRUE) -
                          (log(sb) + dmhn(u, a, 1, z, log = TRUE)))), 1e-9,
                label = sprintf("alpha=%g z=%g", a, z))
    }
  }
})

test_that("the density far below a peak keeps log(x / v*) where it leaves the normal range", {
  # Measured from the peak v*, log f(x) needs log(x / v*), and x / v*
  # underflows for a subnormal x, or at x = 1e-300 with v* above 1e8; it is
  # then taken as log x - log v*.  At (alpha, 1, 2e6) the law is measured
  # from the peak (z = 2e6), and
  #   log f(x) = (alpha - 1) log x + gamma x - x^2 - log C,
  #   log C = gamma^2 / 4 + (alpha - 1) log(gamma / 2) + log(sqrt(pi))
  #           + log(1 + (alpha - 1)(alpha - 2) / gamma^2 + ...),
  # the last term 3e-13 at alpha = 0.3 and 0 at alpha = 2.  log C is 1e12 in
  # size, so the result carries about 1e-4 of rounding; 1e-3 allows for it,
  # far below the (alpha - 1) log x term (500) that would otherwise be lost,
  # and the density for alpha < 1, then +Inf, is 0 here.
  g <- 2e6
  for (a in c(0.3, 2)) {
    log_C <- g^2 / 4 + (a - 1) * log(g / 2) + 0.5 * log(pi)
    for (x in c(1e-320, 1e-300)) {
      info <- sprintf("alpha=%g x=%g", a, x)
      ref <- (a - 1) * log(x) + g * x - x^2 - log_C
      expect_lt(abs(dmhn(x, a, 1, g, log = TRUE) - ref), 1e-3, label = info)
      expect_lt(abs(dmhn(x, c(a, a), 1, g, log = TRUE)[1] - ref), 1e-3,
                label = info)
      expect_identical(dmhn(x, a, 1, g), 0, label = info)
    }
  }
})

test_that("the density holds at the smallest positive shape", {
  # alpha / 2 underflows to 0 at alpha = 2^-1074, and with it the first term
  # of the series for the normalising constant.  As alpha -> 0, Psi[alpha/2, z]
  # is 2 / alpha up to a relative alpha, so the log density at x = 1 is
  # log(alpha) + gamma - beta.
  alpha <- 2^-1074
  for (gamma in c(-1, 0.001, 1, 20, 21)) {
    expect_equal(dmhn(1, alpha, 1, gamma, log = TRUE), log(alpha) + gamma - 1,
                 tolerance = 1e-13)
  }
  expect_lt(dmhn(1, alpha, 1, 1), 1e-300)
})
