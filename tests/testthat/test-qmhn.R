# Tests for qmhn() -- the MHN quantile function.

# ============================================================
# 1. Boundary values
# ============================================================

test_that("qmhn(0) == 0 and qmhn(1) == Inf", {
  for (a in c(0.5, 1, 2, 5)) for (b in c(0.5, 1, 3)) for (g in c(-2, 0, 1, 3)) {
    expect_equal(qmhn(0, a, b, g), 0,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
    expect_equal(qmhn(1, a, b, g), Inf,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("qmhn handles log.p and lower.tail at boundaries", {
  expect_equal(qmhn(-Inf, 2, 1, 1, log.p = TRUE), 0)
  expect_equal(qmhn(0, 2, 1, 1, log.p = TRUE), Inf)         # log(1) = 0
  expect_equal(qmhn(0, 2, 1, 1, lower.tail = FALSE), Inf)
  expect_equal(qmhn(1, 2, 1, 1, lower.tail = FALSE), 0)
})

test_that("qmhn returns NaN for out-of-range probabilities", {
  expect_true(is.nan(qmhn(-0.5, 2, 1, 1)))
  expect_true(is.nan(qmhn(2,    2, 1, 1)))
  expect_true(is.nan(qmhn(1.5,  2, 1, 1, log.p = TRUE)))   # log p > 0
})

# ============================================================
# 2. Round-trip: pmhn(qmhn(p)) ≈ p
# ============================================================

test_that("round-trip pmhn(qmhn(p)) ~= p", {
  p_vals <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
  params <- expand.grid(alpha = c(0.5, 1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-3, -1, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    x <- qmhn(p_vals, a, b, g)
    back <- pmhn(x, a, b, g)
    expect_equal(back, p_vals, tolerance = 1e-6,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 3. Reverse round-trip: qmhn(pmhn(x)) ≈ x
# ============================================================

test_that("reverse round-trip qmhn(pmhn(x)) ~= x", {
  x_vals <- c(0.1, 0.3, 0.5, 1, 2, 3)
  params <- expand.grid(alpha = c(1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0.5, 2))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    p <- pmhn(x_vals, a, b, g)
    back <- qmhn(p, a, b, g)
    # Allow looser absolute tolerance because we are inverting through F.
    expect_equal(back, x_vals, tolerance = 1e-5,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 4. Monotonicity in p
# ============================================================

test_that("qmhn is monotone non-decreasing in p", {
  p <- seq(0.001, 0.999, length.out = 200)
  params <- expand.grid(alpha = c(0.5, 1, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    x <- qmhn(p, a, b, g)
    expect_true(all(diff(x) >= -1e-8),
                label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 5. lower.tail / log.p consistency
# ============================================================

test_that("qmhn lower.tail equivalence", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  params <- list(c(2, 1, 1), c(1.5, 0.5, -1), c(0.5, 1, 2))
  for (par in params) {
    a <- par[1]; b <- par[2]; g <- par[3]
    lower <- qmhn(p, a, b, g)
    upper <- qmhn(1 - p, a, b, g, lower.tail = FALSE)
    expect_equal(lower, upper, tolerance = 1e-8,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("qmhn log.p equivalence", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  params <- list(c(2, 1, 1), c(1.5, 0.5, -1), c(0.5, 1, 2))
  for (par in params) {
    a <- par[1]; b <- par[2]; g <- par[3]
    linear <- qmhn(p, a, b, g)
    logp   <- qmhn(log(p), a, b, g, log.p = TRUE)
    expect_equal(linear, logp, tolerance = 1e-8,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 6. Special-case dispatch
# ============================================================

test_that("gamma = 0: qmhn == sqrt(qgamma(p, alpha/2, scale = 1/beta))", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  for (a in c(0.5, 1, 1.5, 2, 4)) for (b in c(0.5, 1, 3)) {
    ref <- sqrt(qgamma(p, shape = a/2, scale = 1/b))
    got <- qmhn(p, alpha = a, beta = b, gamma = 0)
    expect_equal(got, ref, tolerance = 1e-10,
                 label = sprintf("a=%.1f b=%.1f", a, b))
  }
})

test_that("alpha = 1: qmhn matches truncated normal inverse", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  for (b in c(0.5, 1, 2)) for (g in c(-2, 0, 1, 3)) {
    mu <- g / (2 * b); sigma <- 1 / sqrt(2 * b)
    Phi_low <- pnorm(0, mu, sigma)
    # Reference via direct inversion of the truncated-normal CDF.
    ref <- qnorm(Phi_low + p * (1 - Phi_low), mu, sigma)
    got <- qmhn(p, alpha = 1, beta = b, gamma = g)
    expect_equal(got, ref, tolerance = 1e-8,
                 label = sprintf("b=%.1f g=%+.1f", b, g))
  }
})

# ============================================================
# 7. Vectorization / recycling
# ============================================================

test_that("qmhn vectorizes p", {
  p <- c(0.1, 0.5, 0.9)
  out <- qmhn(p, 2, 1, 1)
  expect_length(out, 3)
  expect_equal(out, sapply(p, function(pp) qmhn(pp, 2, 1, 1)),
               tolerance = 1e-10)
})

test_that("qmhn recycles parameter vectors", {
  p <- c(0.3, 0.5, 0.7)
  a <- c(1, 2)
  b <- 1
  g <- c(0, 1, -1)
  n <- max(length(p), length(a), length(b), length(g))
  out <- qmhn(p, a, b, g)
  expect_length(out, n)
  ref <- vapply(seq_len(n), function(i) {
    qmhn(p[((i-1) %% length(p)) + 1],
         a[((i-1) %% length(a)) + 1],
         b[((i-1) %% length(b)) + 1],
         g[((i-1) %% length(g)) + 1])
  }, numeric(1))
  expect_equal(out, ref, tolerance = 1e-10)
})

test_that("qmhn with empty p returns numeric(0)", {
  expect_equal(qmhn(numeric(0), 2, 1, 1), numeric(0))
})

# ============================================================
# 8. NA / NaN propagation and parameter validation
# ============================================================

test_that("qmhn propagates NA in p", {
  out <- qmhn(c(0.5, NA, 0.9), 2, 1, 1)
  expect_true(is.na(out[2]) && !is.nan(out[2]))
  expect_false(is.na(out[1]))
})

test_that("qmhn propagates NaN in p", {
  out <- qmhn(c(0.5, NaN, 0.9), 2, 1, 1)
  expect_true(is.nan(out[2]))
})

test_that("qmhn errors on invalid parameters", {
  expect_error(qmhn(0.5, alpha = 0, beta = 1, gamma = 0), "alpha")
  expect_error(qmhn(0.5, alpha = -1, beta = 1, gamma = 0), "alpha")
  expect_error(qmhn(0.5, alpha = 1, beta = 0, gamma = 0), "beta")
  expect_error(qmhn(0.5, alpha = 1, beta = 1, gamma = Inf), "gamma")
  expect_error(qmhn(0.5, alpha = NA, beta = 1, gamma = 0), "alpha")
})
