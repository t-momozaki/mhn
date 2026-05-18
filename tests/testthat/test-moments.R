# Tests for mhn_mean, mhn_var, mhn_skewness, mhn_kurtosis, mhn_mode.
#
# Cross-checks each closed-form / recurrence result against numerical
# integration of the corresponding density-weighted integral, plus
# specialised tests for special-case identities, Lemma 4 bounds, and
# parameter validation.

# ============================================================
# 1. Mean: cross-check against numerical integration
# ============================================================

test_that("mhn_mean matches numerical integration", {
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 3, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    int_mean <- integrate(
      function(x) x * dmhn(x, alpha = a, beta = b, gamma = g),
      lower = 0, upper = Inf,
      rel.tol = 1e-8, subdivisions = 200
    )$value
    expect_equal(mhn_mean(a, b, g), int_mean,
                 tolerance = 1e-5,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 2. Variance: cross-check against numerical integration
# ============================================================

test_that("mhn_var matches numerical integration", {
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 5),
    beta = c(0.5, 1, 5),
    gamma = c(-2, 0, 2)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    mu <- mhn_mean(a, b, g)
    int_var <- integrate(
      function(x) (x - mu)^2 * dmhn(x, alpha = a, beta = b, gamma = g),
      lower = 0, upper = Inf,
      rel.tol = 1e-8, subdivisions = 200
    )$value
    expect_equal(mhn_var(a, b, g), int_var,
                 tolerance = 1e-4,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 3. Variance upper bound: Var(X) <= 1/(2 beta) when alpha >= 1
#    (Sun et al. 2023, Lemma 4c)
# ============================================================

test_that("variance upper bound holds for alpha >= 1", {
  test_cases <- expand.grid(
    alpha = c(1, 1.5, 2, 5, 10),
    beta = c(0.5, 1, 5),
    gamma = c(-5, -1, 0, 1, 5)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    b <- test_cases$beta[i]
    g <- test_cases$gamma[i]
    v <- mhn_var(a, b, g)
    expect_true(v <= 1 / (2 * b) + 1e-10,
                label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
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

test_that("recurrence E(X^2) = alpha/(2 beta) + gamma/(2 beta) E(X)", {
  for (a in c(0.5, 2, 5)) {
    for (g in c(-2, 0, 2)) {
      mu <- mhn_mean(a, 1, g)
      v <- mhn_var(a, 1, g)
      ex2_from_var <- v + mu^2
      ex2_from_rec <- a / 2 + g / 2 * mu  # beta = 1
      expect_equal(ex2_from_var, ex2_from_rec,
                   tolerance = 1e-8,
                   label = sprintf("alpha=%.1f, gamma=%.1f", a, g))
    }
  }
})

# ============================================================
# 6. Skewness: cross-check against numerical integration
# ============================================================

test_that("mhn_skewness matches numerical integration", {
  test_cases <- list(
    list(alpha = 2, beta = 1, gamma = 0),
    list(alpha = 2, beta = 1, gamma = 2),
    list(alpha = 3, beta = 2, gamma = -1),
    list(alpha = 0.5, beta = 1, gamma = 1),
    list(alpha = 5, beta = 1, gamma = 3)
  )
  for (tc in test_cases) {
    a <- tc$alpha; b <- tc$beta; g <- tc$gamma
    mu <- mhn_mean(a, b, g)
    sigma <- sqrt(mhn_var(a, b, g))
    int_skew <- integrate(
      function(x) ((x - mu) / sigma)^3 * dmhn(x, a, b, g),
      lower = 0, upper = Inf,
      rel.tol = 1e-8, subdivisions = 200
    )$value
    expect_equal(mhn_skewness(a, b, g), int_skew,
                 tolerance = 1e-3,
                 label = sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f", a, b, g))
  }
})

# ============================================================
# 7. Excess kurtosis: cross-check against numerical integration
# ============================================================

test_that("mhn_kurtosis matches numerical integration", {
  test_cases <- list(
    list(alpha = 2, beta = 1, gamma = 0),
    list(alpha = 2, beta = 1, gamma = 2),
    list(alpha = 3, beta = 2, gamma = -1),
    list(alpha = 5, beta = 1, gamma = 3)
  )
  for (tc in test_cases) {
    a <- tc$alpha; b <- tc$beta; g <- tc$gamma
    mu <- mhn_mean(a, b, g)
    sigma <- sqrt(mhn_var(a, b, g))
    int_kurt <- integrate(
      function(x) ((x - mu) / sigma)^4 * dmhn(x, a, b, g),
      lower = 0, upper = Inf,
      rel.tol = 1e-8, subdivisions = 200
    )$value - 3
    expect_equal(mhn_kurtosis(a, b, g), int_kurt,
                 tolerance = 1e-2,
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
  # Interior mode exists iff alpha >= 1 - gamma^2 / (8 beta).
  # alpha=0.9, beta=1, gamma=2: 1 - 4/8 = 0.5 <= 0.9, mode exists.
  mode_val <- mhn_mode(0.9, 1, 2)
  expect_true(!is.na(mode_val) && mode_val > 0)

  # alpha < 1 - gamma^2 / (8 beta): NA.
  # alpha=0.1, beta=1, gamma=0.5: 1 - 0.25/8 = 0.96875 > 0.1, no mode.
  expect_true(is.na(mhn_mode(0.1, 1, 0.5)))
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

# ============================================================
# 12. Mode parameter validation
# ============================================================

test_that("mhn_mode errors on invalid parameters", {
  expect_error(mhn_mode(-1, 1, 0), "alpha must be positive")
  expect_error(mhn_mode(1, -1, 0), "beta must be positive")
})
