# Tests for special-case detection and density dispatch.
#
# Special cases (Sun et al. 2023, Lemma 6):
#   gamma = 0        -> sqrt-Gamma         (Lemma 6a)
#   alpha = 1        -> truncated normal   (Lemma 6b)
#   alpha = 1, g = 0 -> half-normal        (Lemma 6c)

# ============================================================
# 1. Detection predicates
# ============================================================

test_that(".is_sqrt_gamma detects gamma ~ 0", {
  expect_true(mhn:::.is_sqrt_gamma(0))
  expect_true(mhn:::.is_sqrt_gamma(1e-10))
  expect_true(mhn:::.is_sqrt_gamma(-1e-10))
  expect_false(mhn:::.is_sqrt_gamma(0.01))
  expect_false(mhn:::.is_sqrt_gamma(-0.01))
})

test_that(".is_truncated_normal detects alpha ~ 1", {
  expect_true(mhn:::.is_truncated_normal(1))
  expect_true(mhn:::.is_truncated_normal(1 + 1e-10))
  expect_true(mhn:::.is_truncated_normal(1 - 1e-10))
  expect_false(mhn:::.is_truncated_normal(1.01))
  expect_false(mhn:::.is_truncated_normal(0.99))
})

# ============================================================
# 2. Half-normal distribution (alpha = 1, gamma = 0)
# ============================================================

test_that("dmhn matches half-normal density for alpha=1, gamma=0", {
  x <- seq(0.01, 5, length.out = 50)
  for (beta in c(0.5, 1, 5)) {
    # Half-normal density: 2 sqrt(beta / pi) exp(-beta x^2).
    expected <- 2 * sqrt(beta / pi) * exp(-beta * x^2)
    result <- dmhn(x, alpha = 1, beta = beta, gamma = 0)
    expect_equal(result, expected,
                 tolerance = 1e-12,
                 label = sprintf("beta=%.1f", beta))
  }
})

test_that("half-normal density at x=0", {
  for (beta in c(0.5, 1, 5)) {
    expected <- 2 * sqrt(beta / pi)
    expect_equal(dmhn(0, alpha = 1, beta = beta, gamma = 0), expected,
                 tolerance = 1e-12,
                 label = sprintf("beta=%.1f", beta))
  }
})

# ============================================================
# 3. Truncated normal distribution (alpha = 1)
# ============================================================

test_that("dmhn matches truncated normal for alpha=1", {
  x <- seq(0.01, 5, length.out = 50)
  test_cases <- expand.grid(
    beta = c(0.5, 1, 5),
    gamma = c(-5, -1, 1, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    mu <- g / (2 * b)
    sigma <- 1 / sqrt(2 * b)
    # Truncated normal density: dnorm(x, mu, sigma) / pnorm(mu / sigma).
    expected <- dnorm(x, mean = mu, sd = sigma) / pnorm(mu / sigma)
    result <- dmhn(x, alpha = 1, beta = b, gamma = g)
    expect_equal(result, expected,
                 tolerance = 1e-10,
                 label = sprintf("beta=%.1f, gamma=%.1f", b, g))
  }
})

test_that("truncated normal log-density is consistent", {
  x <- c(0.1, 0.5, 1, 2, 5)
  for (gamma in c(-5, 2)) {
    log_d <- dmhn(x, alpha = 1, beta = 1, gamma = gamma, log = TRUE)
    d <- dmhn(x, alpha = 1, beta = 1, gamma = gamma)
    expect_equal(log_d, log(d),
                 tolerance = 1e-12,
                 label = sprintf("gamma=%.1f", gamma))
  }
})

# ============================================================
# 4. sqrt-Gamma distribution (gamma = 0)
# ============================================================

test_that("dmhn matches sqrt-Gamma density for gamma=0", {
  x <- seq(0.01, 5, length.out = 50)
  test_cases <- expand.grid(
    alpha = c(0.5, 1.5, 2, 3, 5),
    beta = c(0.5, 1, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    # sqrt-Gamma density: 2 x * dgamma(x^2, alpha/2, rate = beta).
    expected <- dgamma(x^2, shape = a / 2, rate = b) * 2 * x
    result <- dmhn(x, alpha = a, beta = b, gamma = 0)
    expect_equal(result, expected,
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, beta=%.1f", a, b))
  }
})

test_that("sqrt-Gamma x=0 boundary behavior", {
  # alpha > 1: density = 0.
  expect_equal(dmhn(0, alpha = 2, beta = 1, gamma = 0), 0)
  # alpha < 1: density = Inf.
  expect_equal(dmhn(0, alpha = 0.5, beta = 1, gamma = 0), Inf)
  # alpha = 1: half-normal density at 0.
  expect_equal(dmhn(0, alpha = 1, beta = 1, gamma = 0),
               2 * sqrt(1 / pi), tolerance = 1e-12)
})

# ============================================================
# 5. Cross-validation across the three special cases
# ============================================================

test_that("alpha=1 gamma=0: half-normal = truncated normal = sqrt-Gamma", {
  x <- seq(0.01, 4, length.out = 30)
  for (beta in c(0.5, 1, 5)) {
    # Half-normal density.
    d_hn <- 2 * sqrt(beta / pi) * exp(-beta * x^2)
    # Truncated normal density (gamma = 0 implies mu = 0).
    sigma <- 1 / sqrt(2 * beta)
    d_tn <- dnorm(x, mean = 0, sd = sigma) / pnorm(0 / sigma)
    # sqrt-Gamma density.
    d_sg <- dgamma(x^2, shape = 0.5, rate = beta) * 2 * x
    # dmhn.
    d_mhn <- dmhn(x, alpha = 1, beta = beta, gamma = 0)

    expect_equal(d_mhn, d_hn, tolerance = 1e-12,
                 label = sprintf("mhn vs HN, beta=%.1f", beta))
    expect_equal(d_mhn, d_tn, tolerance = 1e-10,
                 label = sprintf("mhn vs TN, beta=%.1f", beta))
    expect_equal(d_mhn, d_sg, tolerance = 1e-10,
                 label = sprintf("mhn vs sqGamma, beta=%.1f", beta))
  }
})
