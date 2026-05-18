# Tests for the public rmhn() dispatcher.
#
# Coverage:
#   - general goodness-of-fit (KS, mean, var) across regimes
#   - method-agreement (RTDR vs Sun) where both apply
#   - region (c)/(d) boundary continuity
#   - density-vs-histogram chi-square consistency
#   - vectorisation reproducibility (cache-on-change <-> scalar)
#   - edge cases (n=0, n<0, NA propagation, Inf gamma -> NA)
#   - method argument validation
#   - Robert (1995) truncated-normal kernel cross-check
#
# Lower-level kernels (RTDR per region, Sun Algo 1/3) are covered by
# test-rmhn-rtdr-regions.R and test-rmhn-sun.R; this file focuses on the
# integration surface.

mhn_cdf <- function(alpha, beta, gamma) {
  function(qs) {
    vapply(qs, function(qi) {
      if (qi <= 0) return(0)
      stats::integrate(function(u) dmhn(u, alpha, beta, gamma),
                       0, qi, rel.tol = 1e-6)$value
    }, numeric(1))
  }
}

# =====================================================================
# General goodness-of-fit (method-agnostic, default "auto")
# =====================================================================

run_default <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- rmhn(n, alpha, beta, gamma)
  list(
    samples     = x,
    ks          = ks.test(x, mhn_cdf(alpha, beta, gamma)),
    sample_mean = mean(x),
    sample_var  = var(x),
    theory_mean = mhn_mean(alpha, beta, gamma),
    theory_var  = mhn_var(alpha, beta, gamma)
  )
}

regimes <- list(
  list(alpha = 1.5, beta = 1, gamma = -2,  label = "region a (alpha=1.5, gamma=-2)"),
  list(alpha = 5,   beta = 1, gamma =  0.5, label = "region a (alpha=5, gamma=0.5)"),
  list(alpha = 0.7, beta = 1, gamma =  0,   label = "region b/c (alpha=0.7, gamma=0)"),
  list(alpha = 0.3, beta = 1, gamma = -5,   label = "region c (alpha=0.3, gamma=-5)"),
  list(alpha = 0.3, beta = 1, gamma =  5,   label = "region d (alpha=0.3, gamma=5)")
)

for (r in regimes) {
  test_that(paste("default rmhn KS+moments:", r$label), {
    skip_on_cran()
    out <- run_default(r$alpha, r$beta, r$gamma, n = 10000)
    expect_gt(out$ks$p.value, 0.001)
    # mean: 5 sigma tolerance
    se_mean <- sqrt(out$theory_var / 10000)
    expect_lt(abs(out$sample_mean - out$theory_mean), 5 * se_mean)
    # variance: 10 percent relative tolerance
    expect_lt(abs(out$sample_var - out$theory_var) / out$theory_var, 0.10)
  })
}

# =====================================================================
# method=rtdr vs method=sun agreement (Sun-applicable regimes)
# =====================================================================

method_agreement_cases <- list(
  list(alpha = 2, beta = 1, gamma = 1,  label = "Sun Algo 1 (alpha=2, gamma=1)"),
  list(alpha = 5, beta = 1, gamma = 5,  label = "Sun Algo 1 (alpha=5, gamma=5)"),
  list(alpha = 3, beta = 1, gamma = -2, label = "Sun Algo 3 (alpha=3, gamma=-2)"),
  list(alpha = 5, beta = 1, gamma = -5, label = "Sun Algo 3 (alpha=5, gamma=-5)")
)

for (r in method_agreement_cases) {
  test_that(paste("RTDR vs Sun two-sample KS:", r$label), {
    skip_on_cran()
    set.seed(1); xr <- rmhn(10000, r$alpha, r$beta, r$gamma, method = "rtdr")
    set.seed(2); xs <- rmhn(10000, r$alpha, r$beta, r$gamma, method = "sun")
    expect_gt(suppressWarnings(ks.test(xr, xs))$p.value, 0.001)
  })
}

# =====================================================================
# region (c)/(d) boundary continuity
# =====================================================================

test_that("region (c)/(d) boundary: sample mean varies smoothly", {
  skip_on_cran()
  # Region (c)/(d) boundary at gamma_star = 2(1 - sqrt(1 - 2*alpha)).
  # The grid is +/-10% around gamma_star (wider than spec's +/-5%) to
  # keep the signal-to-noise ratio comfortably above 1: at +/-5% the
  # adjacent-mean change (~3e-3) is the same magnitude as the standard
  # error at N=5000 (~6e-3), so monotonicity gets dominated by noise.
  alpha <- 0.3
  gamma_star <- 2 * (1 - sqrt(1 - 2 * alpha))
  gammas <- gamma_star * c(0.90, 0.95, 1.0, 1.05, 1.10)
  N <- 10000
  means <- vapply(gammas, function(g) {
    set.seed(1); mean(rmhn(N, alpha = alpha, gamma = g))
  }, numeric(1))
  expect_gt(stats::cor(gammas, means, method = "spearman"), 0.95)
})

# =====================================================================
# density-vs-histogram chi-square consistency
# =====================================================================

test_that("chi-square: histogram of rmhn matches dmhn", {
  skip_on_cran()
  alpha <- 2; beta <- 1; gamma <- 0.5
  N <- 50000
  set.seed(1)
  x <- rmhn(N, alpha, beta, gamma)
  q_hi <- stats::uniroot(function(q)
    stats::integrate(function(u) dmhn(u, alpha, beta, gamma),
                     0, q, rel.tol = 1e-6)$value - 0.999,
    c(0.01, 50))$root
  bins <- seq(0, q_hi, length.out = 51)
  obs <- as.integer(table(cut(x[x < q_hi], breaks = bins, include.lowest = TRUE)))
  expected <- vapply(seq_len(50), function(k) {
    stats::integrate(function(u) dmhn(u, alpha, beta, gamma),
                     bins[k], bins[k + 1], rel.tol = 1e-6)$value
  }, numeric(1))
  expected <- expected / sum(expected) * sum(obs)
  # Pool low-expectation bins (chi-square requires E >= ~5 per bin).
  keep <- expected >= 5
  cs <- sum((obs[keep] - expected[keep])^2 / expected[keep])
  df <- sum(keep) - 1
  expect_gt(stats::pchisq(cs, df, lower.tail = FALSE), 0.001)
})

# =====================================================================
# vectorization reproducibility
# =====================================================================

test_that("vectorized rmhn equals scalar repetition under same seed", {
  alphas <- c(2, 3, 5, 2, 1)
  set.seed(1)
  v <- rmhn(5, alpha = alphas)
  set.seed(1)
  s <- vapply(alphas, function(a) rmhn(1, alpha = a), numeric(1))
  expect_equal(v, s)
})

test_that("length-1 vector equals scalar default", {
  set.seed(1); a <- rmhn(100, alpha = 2)
  set.seed(1); b <- rmhn(100, alpha = c(2))
  expect_identical(a, b)
})

# =====================================================================
# edge cases
# =====================================================================

test_that("rmhn(0, ...) returns numeric(0)", {
  expect_identical(rmhn(0, alpha = 2), numeric(0))
  expect_identical(rmhn(0, alpha = c(1, 2, 3)), numeric(0))
})

test_that("rmhn(-1, ...) errors with 'non-negative'", {
  expect_error(rmhn(-1, alpha = 2), "non-negative")
})

test_that("NA in alpha propagates to NA outputs", {
  set.seed(1)
  x <- rmhn(5, alpha = c(1, NA, 2, NA, 3))
  expect_equal(is.na(x), c(FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_true(all(x[!is.na(x)] > 0))
})

test_that("NA in beta and gamma propagates to NA outputs", {
  set.seed(1)
  x_b <- rmhn(4, alpha = 2, beta = c(1, NA, 1, NA))
  expect_equal(is.na(x_b), c(FALSE, TRUE, FALSE, TRUE))

  x_g <- rmhn(4, alpha = 2, gamma = c(0.5, NA, 0.5, NA))
  expect_equal(is.na(x_g), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("non-finite gamma (Inf, -Inf, NaN) propagates to NA", {
  x <- rmhn(5, alpha = 2,
            gamma = c(0.5, Inf, 0.5, -Inf, NaN))
  expect_equal(is.na(x), c(FALSE, TRUE, FALSE, TRUE, TRUE))
})

test_that("invalid alpha or beta still throws", {
  expect_error(rmhn(5, alpha = -1),  "alpha must be positive")
  expect_error(rmhn(5, alpha = 0),   "alpha must be positive")
  expect_error(rmhn(5, beta = 0),    "beta must be positive")
  expect_error(rmhn(5, beta = -0.5), "beta must be positive")
})

# =====================================================================
# method argument validation
# =====================================================================

test_that("method=sun rejects scalar (alpha<1, gamma>0)", {
  expect_error(rmhn(10, alpha = 0.5, gamma = 1, method = "sun"),
               "alpha<1 and gamma>0")
})

test_that("method=sun rejects vector containing (alpha<1, gamma>0)", {
  expect_error(
    rmhn(5, alpha = c(2, 0.3, 2, 2, 2), gamma = 1, method = "sun"),
    "alpha<1 and gamma>0"
  )
})

test_that("invalid method string errors via match.arg", {
  expect_error(rmhn(10, method = "invalid"))
})

test_that("method=sun NA elements skip prescan", {
  set.seed(1)
  expect_silent(
    rmhn(5, alpha = c(2, NA, 2, NA, 3), gamma = 1, method = "sun")
  )
})

# =====================================================================
# Special-case dispatch produces reasonable output
# =====================================================================

test_that("Special-case half-normal (alpha=1, gamma=0) matches HN(1/sqrt(2))", {
  skip_on_cran()
  set.seed(1)
  x <- rmhn(20000, alpha = 1, beta = 1, gamma = 0)
  # HN(sigma=1/sqrt(2)) has mean sigma * sqrt(2/pi) = 1/sqrt(pi)
  expect_lt(abs(mean(x) - 1 / sqrt(pi)), 0.02)
})

test_that("Special-case sqrt-Gamma (alpha=3, gamma=0) matches sqrt(Gamma(1.5, 1))", {
  skip_on_cran()
  set.seed(1)
  x <- rmhn(20000, alpha = 3, beta = 1, gamma = 0)
  # KS against the reference distribution
  ref_cdf <- function(q) {
    vapply(q, function(qi) {
      if (qi <= 0) return(0)
      stats::pgamma(qi^2, shape = 1.5, rate = 1)
    }, numeric(1))
  }
  expect_gt(ks.test(x, ref_cdf)$p.value, 0.001)
})

test_that("Special-case truncated normal (alpha=1, gamma=1) matches TN", {
  skip_on_cran()
  set.seed(1)
  x <- rmhn(20000, alpha = 1, beta = 1, gamma = 1)
  mu <- 1 / 2; sigma <- 1 / sqrt(2)
  ref_cdf <- function(q) {
    p0 <- stats::pnorm(0, mu, sigma)
    vapply(q, function(qi) {
      if (qi <= 0) return(0)
      (stats::pnorm(qi, mu, sigma) - p0) / (1 - p0)
    }, numeric(1))
  }
  expect_gt(ks.test(x, ref_cdf)$p.value, 0.001)
})

# =====================================================================
# Robert (1995) truncated-normal kernel cross-check
# =====================================================================

test_that("rtnorm_robert_zero matches truncated normal CDF: mu > 0", {
  skip_on_cran()
  set.seed(1)
  x <- mhn:::.rtnorm_robert_cpp(10000, mu = 1, sigma = 0.5)
  cdf <- function(q) {
    p0 <- stats::pnorm(0, 1, 0.5)
    vapply(q, function(qi) {
      if (qi <= 0) return(0)
      (stats::pnorm(qi, 1, 0.5) - p0) / (1 - p0)
    }, numeric(1))
  }
  expect_gt(ks.test(x, cdf)$p.value, 0.001)
})

test_that("rtnorm_robert_zero matches truncated normal CDF: mu < 0", {
  skip_on_cran()
  set.seed(1)
  x <- mhn:::.rtnorm_robert_cpp(10000, mu = -2, sigma = 1)
  cdf <- function(q) {
    p0 <- stats::pnorm(0, -2, 1)
    vapply(q, function(qi) {
      if (qi <= 0) return(0)
      (stats::pnorm(qi, -2, 1) - p0) / (1 - p0)
    }, numeric(1))
  }
  expect_gt(ks.test(x, cdf)$p.value, 0.001)
})

test_that("rtnorm_robert_zero edge cases", {
  expect_identical(mhn:::.rtnorm_robert_cpp(0, mu = 0, sigma = 1), numeric(0))
  expect_error(mhn:::.rtnorm_robert_cpp(-1, mu = 0, sigma = 1), "non-negative")
  expect_error(mhn:::.rtnorm_robert_cpp(5, mu = 0, sigma = 0), "sigma must be positive")
  expect_error(mhn:::.rtnorm_robert_cpp(5, mu = Inf, sigma = 1), "mu must be finite")
})

# =====================================================================
# auto-path n-dependent dispatch verification
# Strategy: rebuild_cache is RNG-deterministic, so under the same seed
# `method = "auto"` and the method that auto picks must produce
# bit-identical outputs.
# =====================================================================

test_that("auto -> Sun A1 for gamma>0 & alpha>1 (any n)", {
  for (n_val in c(1L, 100L, 1000L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 5, gamma = 2,  method = "auto")
    set.seed(7L); s <- rmhn(n_val, alpha = 5, gamma = 2,  method = "sun")
    expect_identical(a, s)
  }
})

test_that("auto -> RTDR for gamma>0 & alpha<1 (Sun A2 unimplemented)", {
  for (n_val in c(1L, 100L, 1000L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 0.3, gamma = 5, method = "auto")
    set.seed(7L); r <- rmhn(n_val, alpha = 0.3, gamma = 5, method = "rtdr")
    expect_identical(a, r)
  }
})

test_that("auto -> Sun A3 for gamma<0 with small n (Gibbs)", {
  # samples_per_setup < 25 -> Sun A3
  for (n_val in c(1L, 5L, 24L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 5, gamma = -10, method = "auto")
    set.seed(7L); s <- rmhn(n_val, alpha = 5, gamma = -10, method = "sun")
    expect_identical(a, s)
  }
})

test_that("auto -> RTDR for gamma<0 with large n (batch)", {
  # samples_per_setup >= 25 -> RTDR
  for (n_val in c(25L, 100L, 1000L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 5, gamma = -10, method = "auto")
    set.seed(7L); r <- rmhn(n_val, alpha = 5, gamma = -10, method = "rtdr")
    expect_identical(a, r)
  }
})

test_that("auto special-case intercepts unaffected by n", {
  # alpha = 1 -> trunc-normal special case (intercepted before dispatch)
  for (n_val in c(1L, 100L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 1, gamma = 2, method = "auto")
    set.seed(7L); r <- rmhn(n_val, alpha = 1, gamma = 2, method = "rtdr")
    expect_identical(a, r)
  }
  # gamma = 0 -> sqrt-Gamma special case
  for (n_val in c(1L, 100L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 5, gamma = 0, method = "auto")
    set.seed(7L); r <- rmhn(n_val, alpha = 5, gamma = 0, method = "rtdr")
    expect_identical(a, r)
  }
})

test_that("auto vectorised params: per-sample-amortisation drives dispatch", {
  # max(na, nb, ng) >= n -> samples_per_setup = 1 -> Sun A3 for gamma<0
  set.seed(7L)
  a <- rmhn(5L, alpha = c(5, 5, 5, 5, 5), gamma = -10, method = "auto")
  set.seed(7L)
  s <- rmhn(5L, alpha = c(5, 5, 5, 5, 5), gamma = -10, method = "sun")
  expect_identical(a, s)
})
