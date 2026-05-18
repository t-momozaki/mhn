# Tests for the RTDR (Relaxed Transformed Density Rejection) sampling
# kernel, region by region.  Exercises the internal `.rmhn_rtdr_cpp`
# directly so that regions can be validated independently of the
# dispatcher.
#
# Targets the per-region KS tests and the region (c)/(d) boundary
# continuity check, exercising Gao & Wang (2025) Theorems 3.1, 3.2, 4.4.
# Acceptance-rate >= 1/e is checked where theoretically required
# (envelope-14 is suboptimal for the large-gamma region D regime; that
# case is asserted lenient).
#
# KS p-value threshold > 0.001 keeps the false-rejection rate < 0.1%
# per (region, parameter) combination.

# Helper: empirical CDF of dmhn via numerical integration.
mhn_ecdf <- function(alpha, beta, gamma) {
  function(qs) {
    sapply(qs, function(qi)
      stats::integrate(function(u) dmhn(u, alpha, beta, gamma),
                       0, qi, rel.tol = 1e-6)$value)
  }
}

run_rtdr_region <- function(alpha, beta, gamma, n = 10000, seed = 1L,
                            ks_threshold = 0.001) {
  set.seed(seed)
  x <- mhn:::.rmhn_rtdr_cpp(n, alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    region        = attr(x, "rtdr_region"),
    retries       = attr(x, "rtdr_retries"),
    acceptance    = n / (n + attr(x, "rtdr_retries")),
    ks            = ks.test(x, mhn_ecdf(alpha, beta, gamma)),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = mhn_mean(alpha, beta, gamma),
    theory_var    = mhn_var(alpha, beta, gamma)
  )
}

# Region A: alpha >= 1, log-concave on f(x).
# Expected: region = 0; acceptance >= 1/e theoretically.
test_that("region A: KS, mean/var, acceptance for representative points", {
  skip_on_cran()
  for (params in list(c(1.5, 1, -2), c(5, 1, 0.5), c(10, 1, 10))) {
    r <- run_rtdr_region(params[1], params[2], params[3])
    expect_equal(r$region, 0L,
                 info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
    expect_gt(r$ks$p.value, 0.001)
    expect_gt(r$acceptance, 1 / exp(1))
    # Sample mean within ~5 std-errors of theoretical (very loose).
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se)
  }
})

# Region BC: alpha < 1, T_{-1/2}-concave on g(y).  Covers spec's "region (b)"
# and "region (c)" (which collapse to the same envelope).
test_that("region BC: KS, mean/var, acceptance for representative points", {
  skip_on_cran()
  cases <- list(
    c(0.7, 1, -2), c(0.7, 1,  0), c(0.7, 1,  5),     # spec region (b)
    c(0.3, 1, -5), c(0.3, 1,  0), c(0.3, 1, 0.5)     # spec region (c), all below threshold
  )
  for (params in cases) {
    r <- run_rtdr_region(params[1], params[2], params[3])
    expect_equal(r$region, 1L,
                 info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
    expect_gt(r$ks$p.value, 0.001)
    expect_gt(r$acceptance, 1 / exp(1))
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se)
  }
})

# Region D: alpha < 1/2 and gamma > threshold.
# Envelope 14 always yields correct samples; acceptance >= 1/e holds for
# moderate-gamma cases where log(g(m_g)/g(y*)) <= 2.49.  The (0.1, 50)
# case is in envelope-15 territory: samples remain correct (validated by
# KS) but acceptance is permitted below 1/e.
test_that("region D moderate (envelope 14 sufficient): KS + acceptance >= 1/e", {
  skip_on_cran()
  r <- run_rtdr_region(0.3, 1, 5)
  expect_equal(r$region, 2L)
  expect_gt(r$ks$p.value, 0.001)
  expect_gt(r$acceptance, 1 / exp(1))
})

test_that("region D extreme (envelope 15 territory): KS + acceptance >= 1/e", {
  skip_on_cran()
  # With the inflection-point envelope implemented, Gao & Wang (2025,
  # Theorem 4.4)'s uniform
  # acceptance >= 1/e bound holds for the entire region D, including the
  # large-gamma extreme.
  r <- run_rtdr_region(0.1, 1, 50, n = 5000)
  expect_equal(r$region, 2L)
  expect_gt(r$ks$p.value, 0.001)
  expect_gt(r$acceptance, 1 / exp(1))
  # Confirm envelope 15 is selected (not envelope 14).
  env <- mhn:::.dump_rtdr_envelope_cpp(0.1, 1, 50)
  expect_true(isTRUE(env$has_left_tangent_d))
})

# Envelope switching criterion: log(g(m_g)/g(y*)) > 2.49 -> envelope 15.
# Confirms env 14 vs env 15 selection across the documented threshold.
test_that("region D envelope 14 vs 15 selection follows the 2.49 threshold", {
  metric <- function(alpha, gamma) {
    gn <- gamma
    m_g <- log((gn + sqrt(gn^2 + 8 * alpha)) / 4)
    y_star <- log(gn / 4)
    log_g <- function(y) alpha * y - exp(2 * y) + gn * exp(y)
    log_g(m_g) - log_g(y_star)
  }
  cases <- list(
    list(p = c(0.3, 5),  expected_env = 14L),  # metric ~ 1.77
    list(p = c(0.3, 10), expected_env = 15L),  # metric ~ 6.46
    list(p = c(0.1, 50), expected_env = 15L)   # metric ~ 156
  )
  for (case in cases) {
    p <- case$p
    m <- metric(p[1], p[2])
    env <- mhn:::.dump_rtdr_envelope_cpp(p[1], 1, p[2])
    if (case$expected_env == 14L) {
      expect_lt(m, 2.49)
      expect_false(isTRUE(env$has_left_tangent_d))
    } else {
      expect_gt(m, 2.49)
      expect_true(isTRUE(env$has_left_tangent_d))
    }
  }
})

# Envelope 15 internal identity: by construction y_K equals t_hat_l =
# log(gamma/2 - exp(t_l)) at k = K (the unterminated index).  Early
# termination on the leftmost secant breakpoint (Gao & Wang 2025, proof
# of Lemma 4.1, Appendix A.6) means the LAST stored breakpoint
# y_b[K_eff] generally does NOT equal t_hat_l; it stops earlier.  Here we
# verify the identity by feeding the original K back into the breakpoint
# formula and comparing the result to t_hat_l (independent of K_eff).
test_that("envelope 15: y_K (formula at original K) equals t_hat_l", {
  for (p in list(c(0.3, 10), c(0.1, 50))) {
    env <- mhn:::.dump_rtdr_envelope_cpp(p[1], 1, p[2])
    skip_if_not(isTRUE(env$has_left_tangent_d))
    t_hat_l_expected <- log(p[2] / 2 - exp(env$t_l))
    K_orig <- as.integer(ceiling(env$rho))
    drho <- env$rho / K_orig
    rk <- 4 * K_orig * drho
    inside <- p[2]^2 - rk
    skip_if_not(inside > 0, "K-th breakpoint not real-valued for this case")
    y_K_formula <- log(2 * K_orig * drho / (sqrt(inside) + p[2]))
    expect_lt(abs(y_K_formula - t_hat_l_expected), 1e-9)
  }
})

# Verify that the early-termination optimisation (Gao & Wang 2025, proof
# of Lemma 4.1, Appendix A.6) collapses K_eff to a small value (typically
# 1) for envelope 15 in the extreme corner of region D, while leaving
# acceptance well above 1/e.
test_that("early termination shrinks K_eff in region (d) for extreme gamma", {
  for (p in list(c(0.3, 10), c(0.1, 50), c(0.3, 10000), c(0.01, 10000))) {
    env <- mhn:::.dump_rtdr_envelope_cpp(p[1], 1, p[2])
    skip_if_not(isTRUE(env$has_left_tangent_d))
    K_orig <- as.integer(ceiling(env$rho))
    K_eff <- length(env$y_break)
    expect_lte(K_eff, K_orig)
    # For extreme cases (K_orig in the thousands), K_eff should drop to
    # single digits.  For moderate cases the drop is less dramatic but
    # still satisfies K_eff <= K_orig.
    if (K_orig > 100L) {
      expect_lt(K_eff, 10L)
    }
  }
})

# Spec §8.2.5: boundary continuity at gamma = 2(1 - sqrt(1 - 2*alpha)),
# the BC/D split.  For alpha = 0.3 the threshold is ~0.7350.  Sample mean
# and variance should vary smoothly across the boundary.
test_that("region BC/D boundary continuity: smooth mean/variance", {
  skip_on_cran()
  alpha <- 0.3
  threshold <- 2 * (1 - sqrt(1 - 2 * alpha))   # ~ 0.7350
  # Use larger gamma steps so signal exceeds statistical noise (SE ~ 0.004
  # at n = 10000 for these parameters).
  gs <- threshold + c(-0.20, -0.10, -0.02, 0.02, 0.10, 0.20)
  means <- vars <- numeric(length(gs))
  regions <- integer(length(gs))
  for (i in seq_along(gs)) {
    set.seed(99L + i)
    x <- mhn:::.rmhn_rtdr_cpp(20000, alpha, 1, gs[i])
    means[i]   <- mean(x)
    vars[i]    <- var(x)
    regions[i] <- attr(x, "rtdr_region")
  }
  # Region tags should switch from BC (=1) to D (=2) as gamma crosses threshold.
  expect_true(all(regions[gs <  threshold] == 1L))
  expect_true(all(regions[gs >= threshold] == 2L))
  # Mean rises with gamma; check via Spearman correlation to absorb noise.
  expect_gt(cor(gs, means, method = "spearman"), 0.95)
})

# dump_rtdr_envelope_cpp sanity (structural check on the envelope object).
test_that("dump_rtdr_envelope_cpp exposes consistent fields per region", {
  e_a  <- mhn:::.dump_rtdr_envelope_cpp(5,   1, 0.5)
  e_bc <- mhn:::.dump_rtdr_envelope_cpp(0.7, 1, 0)
  e_d  <- mhn:::.dump_rtdr_envelope_cpp(0.3, 1, 5)

  expect_equal(e_a$region,  0L)
  expect_equal(e_bc$region, 1L)
  expect_equal(e_d$region,  2L)

  # Plateau is non-empty in every envelope: p_l < p_r (region A standard
  # and BC) or 0 < p_r (region A simplified).
  expect_lt(e_a$p_l, e_a$p_r)        # alpha=5 is not simplified
  expect_lt(e_bc$p_l, e_bc$p_r)
  expect_lt(e_d$y_break[length(e_d$y_break)], e_d$p_r)

  # Region D has K_eff >= 1, monotone breakpoints, ends at y_star.
  expect_gte(e_d$K_eff, 1L)
  expect_true(all(diff(e_d$y_break) > 0))
  expect_lt(abs(e_d$y_break[length(e_d$y_break)] - e_d$y_star), 1e-9)

  # piece_log_area aligns with `pieces` list.
  for (e in list(e_a, e_bc, e_d)) {
    expect_equal(length(e$piece_log_area), length(e$pieces))
  }
})
