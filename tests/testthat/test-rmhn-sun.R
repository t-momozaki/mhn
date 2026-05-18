# Tests for Algorithms 1 (gamma > 0, alpha > 1) and 3 (gamma <= 0,
# alpha > 0) of Sun et al. (2023).  Exercises the internal
# `.rmhn_sun_algo1_cpp` and `.rmhn_sun_algo3_cpp` exports directly,
# ahead of the full rmhn() dispatcher.
#
# Acceptance probability lower bounds are asserted only where the
# theoretical bound is proven:
#   - Algo 1, alpha >= 4: >= 0.8 (Sun et al. 2023, Theorem 2e)
#   - Algo 3, alpha >  1: >= 1 / sqrt(2) (Sun et al. 2023, Theorem 4c)
# For 1 < alpha < 4 (Algo 1) and alpha <= 1 (Algo 3) the paper does not
# prove a uniform bound; we only check KS in those regimes.

mhn_ecdf <- function(alpha, beta, gamma) {
  function(qs) {
    sapply(qs, function(qi)
      stats::integrate(function(u) dmhn(u, alpha, beta, gamma),
                       0, qi, rel.tol = 1e-6)$value)
  }
}

# =====================================================================
# Algorithm 1: gamma > 0, alpha > 1
# =====================================================================

run_sun_algo1 <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- mhn:::.rmhn_sun_algo1_cpp(n, alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    proposal      = attr(x, "sun_proposal"),
    retries       = attr(x, "sun_retries"),
    acceptance    = n / (n + attr(x, "sun_retries")),
    ks            = ks.test(x, mhn_ecdf(alpha, beta, gamma)),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = mhn_mean(alpha, beta, gamma),
    theory_var    = mhn_var(alpha, beta, gamma)
  )
}

test_that("Algo 1: KS, mean/var, sample shape across representative points", {
  skip_on_cran()
  cases <- list(c(3, 1, 2), c(5, 1, 5), c(10, 1, 10), c(2, 1, 0.5))
  for (params in cases) {
    r <- run_sun_algo1(params[1], params[2], params[3])
    info <- sprintf("alpha=%g, gamma=%g", params[1], params[3])
    expect_gt(r$ks$p.value, 0.001)
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se, label = info)
    expect_true(r$proposal %in% c("normal", "sqrt_gamma"))
  }
})

test_that("Algo 1: alpha >= 4 yields acceptance >= 0.8 (Theorem 2e)", {
  skip_on_cran()
  for (params in list(c(5, 1, 5), c(10, 1, 10))) {
    r <- run_sun_algo1(params[1], params[2], params[3])
    expect_gt(r$acceptance, 0.8)
  }
})

test_that("Algo 1: rejects alpha <= 1 and gamma <= 0", {
  expect_error(mhn:::.rmhn_sun_algo1_cpp(10, 1.0, 1, 1),
               regexp = "alpha > 1")
  expect_error(mhn:::.rmhn_sun_algo1_cpp(10, 0.5, 1, 1),
               regexp = "alpha > 1")
  expect_error(mhn:::.rmhn_sun_algo1_cpp(10, 2.0, 1, 0),
               regexp = "gamma > 0")
  expect_error(mhn:::.rmhn_sun_algo1_cpp(10, 2.0, 1, -1),
               regexp = "gamma > 0")
})

test_that("Algo 1 dump exposes consistent fields", {
  d <- mhn:::.dump_sun_algo1_cpp(5, 1, 5)
  expect_true(d$chosen %in% c("normal", "sqrt_gamma"))
  expect_gt(d$mu_opt, 0)
  expect_gt(d$delta_opt, 0)
  expect_gt(d$sigma, 0)
  expect_true(is.finite(d$log_K1_minus_K2))
})

# =====================================================================
# Algorithm 3: gamma <= 0, alpha > 0
# =====================================================================

run_sun_algo3 <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- mhn:::.rmhn_sun_algo3_cpp(n, alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    m             = attr(x, "sun_m"),
    r             = attr(x, "sun_r"),
    used_inflex   = attr(x, "sun_used_inflex"),
    retries       = attr(x, "sun_retries"),
    acceptance    = n / (n + attr(x, "sun_retries")),
    ks            = ks.test(x, mhn_ecdf(alpha, beta, gamma)),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = mhn_mean(alpha, beta, gamma),
    theory_var    = mhn_var(alpha, beta, gamma)
  )
}

test_that("Algo 3 alpha > 1: KS + acceptance >= 1/sqrt(2) (Theorem 4c)", {
  skip_on_cran()
  for (params in list(c(3, 1, -2), c(5, 1, -5), c(1.5, 1, -1))) {
    r <- run_sun_algo3(params[1], params[2], params[3])
    info <- sprintf("alpha=%g, gamma=%g", params[1], params[3])
    expect_gt(r$ks$p.value, 0.001)
    expect_gt(r$acceptance, 1 / sqrt(2))
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se, label = info)
  }
})

test_that("Algo 3 alpha <= 1: KS only (no proven uniform acceptance bound)", {
  skip_on_cran()
  for (params in list(c(0.5, 1, -1), c(0.7, 1, -2))) {
    r <- run_sun_algo3(params[1], params[2], params[3])
    info <- sprintf("alpha=%g, gamma=%g", params[1], params[3])
    expect_gt(r$ks$p.value, 0.001)
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se, label = info)
  }
})

test_that("Algo 3 gamma = 0: acceptance is 1", {
  skip_on_cran()
  r <- run_sun_algo3(3, 1, 0)
  expect_gt(r$ks$p.value, 0.001)
  # Sun et al. (2023) remark following Theorem 4: A_neg(m, alpha, beta,
  # gamma) becomes unity when gamma = 0, irrespective of m > 0.
  # Implementation gives 1.0 (no rejections) but allow tiny numerical slack.
  expect_gt(r$acceptance, 0.999)
})

test_that("Algo 3: rejects gamma > 0", {
  expect_error(mhn:::.rmhn_sun_algo3_cpp(10, 2, 1, 1),
               regexp = "gamma <= 0")
  expect_error(mhn:::.rmhn_sun_algo3_cpp(10, 2, 1, 5),
               regexp = "gamma <= 0")
})

test_that("Algo 3 dump exposes consistent fields and ordering invariants", {
  d <- mhn:::.dump_sun_algo3_cpp(3, 1, -2)
  expect_gt(d$m, 0)
  expect_gt(d$m_init, 0)
  expect_gt(d$r, 0); expect_lt(d$r, 1)   # r in (0, 1) by formula
  expect_gt(d$shape, 0)
  expect_gt(d$rate, 0)
  expect_lt(abs(d$m_betam_gam - d$rate), 1e-12)  # both equal m*(beta*m+|gamma|)
  expect_true(is.logical(d$used_inflex_heuristic))
})

# =====================================================================
# Newton-step refinement: uses l(m) = log(A_neg(m)) from Sun et al.
# (2023, Supplementary Section 2.17) as the target.  Verify
# l(m_recommend) >= l(m_init), since Newton on l targets the maximum
# of l (= maximum acceptance probability).  The implementation lives in
# src/mhn_sun.cpp::newton_refine_m().
# =====================================================================
test_that("Algo 3 Newton refinement does not lower l(m)", {
  # Paper l(m) drops the constant log Psi term (omitted in our C++ helper).
  l_no_psi <- function(m, alpha, beta, gamma_abs) {
    A <- beta * m + gamma_abs
    B <- 2 * beta * m + gamma_abs
    C <- beta * m^2 + m * gamma_abs
    sh <- alpha * A / B
    sh * log(C) + log(B) - log(2) -
      (alpha / 2) * log(beta * m^2) - lgamma(sh) - log(A)
  }
  for (params in list(c(3, -2), c(5, -5), c(100, -10), c(10000, -100))) {
    alpha <- params[1]; gamma <- params[2]
    d <- mhn:::.dump_sun_algo3_cpp(alpha, 1, gamma)
    li <- l_no_psi(d$m_init, alpha, 1, abs(gamma))
    lm <- l_no_psi(d$m,      alpha, 1, abs(gamma))
    # Newton on l should not decrease l (within tiny numerical slack).
    # On failure, the test name plus alpha/gamma identify the case.
    expect_gte(lm, li - 1e-9)
  }
})

# Regression guard for the alpha > 100 collapse (acc dropped to 0.001
# under an earlier ad-hoc Newton target that diverged for extreme
# alpha; the current l(m)-based Newton holds the Sun et al. (2023,
# Theorem 4c) bound).
test_that("Algo 3 acceptance maintains Theorem 4c bound for extreme alpha", {
  skip_on_cran()
  for (params in list(c(1000, -100), c(1000, -10), c(10000, -100), c(10000, -10))) {
    alpha <- params[1]; gamma <- params[2]
    set.seed(42)
    x <- mhn:::.rmhn_sun_algo3_cpp(2000, alpha, 1, gamma)
    acc <- 2000 / (2000 + attr(x, "sun_retries"))
    # Sun et al. (2023, Theorem 4c) (alpha > 1, gamma <= 0):
    # acc >= 1 / sqrt(2) ~ 0.707.  Use 0.6 as a safety margin to
    # avoid spurious noise failures.
    expect_gt(acc, 0.6)
    expect_equal(sum(is.na(x)), 0L)
  }
})
