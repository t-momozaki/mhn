# Tests for .mhn_log_normalizing_const() and its internal methods.
#
# The normalizing constant is Psi[alpha/2, gamma/sqrt(beta)].
# Helper: nc(alpha, beta, gamma) = exp(.mhn_log_normalizing_const(...)).

nc <- function(alpha, beta, gamma, ...) {
  exp(mhn:::.mhn_log_normalizing_const(alpha, beta, gamma, ...))
}

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
# 2. alpha = 1 closed form (Sun et al. 2023, Lemma 9c)
#    Psi[1/2, z] = 2 sqrt(pi) exp(z^2 / 4) Phi(z / sqrt(2))
# ============================================================

test_that("alpha = 1: matches Sun et al. 2023 Lemma 9c closed form", {
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
                 tolerance = 1e-8,
                 label = sprintf("beta=%.1f, gamma=%.1f", b, g))
  }
})

# ============================================================
# 3. alpha = 2, gamma >= 0 closed form (Sun et al. 2023, Lemma 9c ratio)
# ============================================================

test_that("alpha = 2, gamma >= 0: matches Sun et al. 2023 Lemma 9c ratio", {
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
                 tolerance = 1e-8,
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
#    (in the series Sum Gamma(a + n/2) z^n / n!, increasing z leaves
#     even-power terms unchanged and increases the odd-power terms)
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
# 6. Recurrence identity (Sun et al. 2023, Lemma 9a): gamma >= 0 only
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
                 tolerance = 1e-6,
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
                 tolerance = 1e-6,
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
  # (beta, gamma) pairs that share the same z = gamma / sqrt(beta).
  alpha <- 3
  z_target <- -2
  beta_vals <- c(0.25, 1, 4, 16)
  vals <- numeric(length(beta_vals))
  for (j in seq_along(beta_vals)) {
    b <- beta_vals[j]
    g <- z_target * sqrt(b)  # gamma = z * sqrt(beta)
    vals[j] <- nc(alpha, b, g)
  }
  for (j in seq_along(vals)[-1]) {
    expect_equal(vals[j], vals[1],
                 tolerance = 1e-4,
                 label = sprintf("beta=%.2f vs beta=%.2f",
                                 beta_vals[j], beta_vals[1]))
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
                 tolerance = 1e-6,
                 label = sprintf("gamma=%.1f", g))
  }
})

# ============================================================
# 10. Log-output consistency
# ============================================================

test_that("log output is consistent with exp", {
  test_cases <- expand.grid(
    alpha = c(0.5, 1, 2, 5),
    gamma = c(-5, 0, 3)
  )
  for (i in seq_len(nrow(test_cases))) {
    a <- test_cases$alpha[i]
    g <- test_cases$gamma[i]
    log_val <- mhn:::.mhn_log_normalizing_const(a, 1, g)
    expect_equal(log_val, log(nc(a, 1, g)),
                 tolerance = 1e-10,
                 label = sprintf("alpha=%.1f, gamma=%.1f", a, g))
  }
})

# ============================================================
# 11. Edge cases
# ============================================================

test_that("very small alpha works", {
  val <- nc(0.01, 1, 0)
  expect_equal(val, gamma(0.005), tolerance = 1e-10)
})

test_that("very large alpha works", {
  val <- nc(100, 1, 0)
  expect_equal(val, gamma(50), tolerance = 1e-6)
})

test_that("very small beta works", {
  val <- nc(2, 0.001, 0)
  expect_equal(val, gamma(1), tolerance = 1e-10)
})

test_that("very large beta works", {
  val <- nc(2, 1000, 0)
  expect_equal(val, gamma(1), tolerance = 1e-10)
})

# ============================================================
# 11. Regression: gamma > 0 series matches numerical integration
#     across a parameter grid.
#
# This catches a bug where the C1/C2 critical-index discriminant
# in .lemma10_C() incorrectly used `alpha_adj * z` instead of
# `alpha_adj * z^2` (the original paper main text has a typesetting
# error; the supplementary proof uses z^2; see the "Errata in
# Sun et al. (2023)" subsection of the package's theory vignette).
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
      tolerance = 1e-7,
      label = sprintf("alpha=%.2f, z=%.2f", a, z)
    )
  }
})
