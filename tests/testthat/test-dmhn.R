# Tests for dmhn() — MHN density function.
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
  # 90 parameter combinations (well beyond the 20+ baseline required by spec).
  for (i in seq_len(nrow(test_params))) {
    a <- test_params$alpha[i]
    b <- test_params$beta[i]
    g <- test_params$gamma[i]
    result <- integrate(
      function(x) dmhn(x, alpha = a, beta = b, gamma = g),
      lower = 0, upper = Inf,
      rel.tol = 1e-8, subdivisions = 200
    )
    expect_equal(result$value, 1,
                 tolerance = 1e-5,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 2. Non-negativity
# ============================================================

test_that("density is non-negative for all x", {
  x <- c(-2, -1, -0.01, 0, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50)
  for (a in c(0.5, 1, 2)) {
    for (g in c(-2, 0, 2)) {
      d <- dmhn(x, alpha = a, beta = 1, gamma = g)
      # Inf is >= 0
      expect_true(all(d >= 0 | is.nan(d)),
                  label = sprintf("alpha=%.1f, gamma=%.1f", a, g))
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

test_that("dmhn returns NA for NA input", {
  result <- dmhn(c(1, NA, 2), alpha = 2, beta = 1, gamma = 0)
  expect_equal(length(result), 3)
  expect_true(is.na(result[2]))
  expect_false(is.na(result[1]))
  expect_false(is.na(result[3]))
})

# ============================================================
# 8. Parameter validation
# ============================================================

test_that("dmhn errors on invalid parameters", {
  expect_error(dmhn(1, alpha = -1), "alpha must be positive")
  expect_error(dmhn(1, alpha = 0), "alpha must be positive")
  expect_error(dmhn(1, beta = -1), "beta must be positive")
  expect_error(dmhn(1, beta = 0), "beta must be positive")
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
# 10. Parameter vectorization (B5)
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
  expect_equal(result, expected, tolerance = 1e-12)
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

test_that("dmhn rejects invalid elements in vector parameters", {
  expect_error(dmhn(1, alpha = c(1, -1), beta = 1, gamma = 0),
               "alpha must be positive")
  expect_error(dmhn(1, alpha = c(1, NA), beta = 1, gamma = 0),
               "alpha must be positive")
  expect_error(dmhn(1, alpha = 1, beta = c(1, 0), gamma = 0),
               "beta must be positive")
  expect_error(dmhn(1, alpha = 1, beta = 1, gamma = c(0, NA)),
               "gamma must be a finite")
  expect_error(dmhn(1, alpha = 1, beta = 1, gamma = c(0, Inf)),
               "gamma must be a finite")
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
