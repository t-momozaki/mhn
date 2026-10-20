# Tests for Algorithms 1 (gamma > 0, alpha > 1) and 3 (gamma <= 0,
# alpha > 0) of Sun et al. (2023).  Exercises the internal
# `.rmhn_sun_algo1_cpp` and `.rmhn_sun_algo3_cpp` exports directly,
# ahead of the full rmhn() dispatcher.
#
# Acceptance probability lower bounds are asserted only where Sun et al.
# prove a bound:
#   - Algo 1, alpha >= 4: >= 0.8 (Sun et al. 2023, Theorem 2e), where the
#     better of its two proposals is taken; that holds up to alpha = 1e18,
#     above which the Normal proposal serves throughout and accepts about
#     1/sqrt(2) under a small tilt (see mhn_sun.h).
#   - Algo 3, alpha >  1: >= 1 / sqrt(2) (Sun et al. 2023, Theorem 4c,
#     stated for alpha > 1; the supplementary proof treats alpha >= 3
#     only).  Theorem 4c covers the matching points m_opt and m_star, not
#     the one-step Newton m used here, so for Algo 3 the check is
#     empirical.
# For 1 < alpha < 4 (Algo 1) and alpha <= 1 (Algo 3) the paper does not
# prove a uniform bound; we only check KS in those regimes.  The references,
# mhn_ref() and mhn_cdf(), are in helper-mhn-cdf.R and are formed from the
# kernel of the law alone.

# =====================================================================
# Algorithm 1: gamma > 0, alpha > 1
# =====================================================================

run_sun_algo1 <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- mhn:::.rmhn_sun_algo1_cpp(n, alpha, beta, gamma)
  ref <- mhn_ref(alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    proposal      = attr(x, "sun_proposal"),
    retries       = attr(x, "sun_retries"),
    acceptance    = n / (n + attr(x, "sun_retries")),
    ref           = ref,
    ks            = ks.test(x, ref$cdf),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = ref$mean,
    theory_var    = ref$var
  )
}

test_that("Algo 1: KS, mean/var, sample shape across representative points", {
  skip_on_cran()
  cases <- list(c(3, 1, 2), c(5, 1, 5), c(10, 1, 10), c(2, 1, 0.5),
                c(3, 4, 4), c(5, 0.25, 2.5))
  for (params in cases) {
    r <- run_sun_algo1(params[1], params[2], params[3])
    info <- sprintf("alpha=%g, gamma=%g", params[1], params[3])
    expect_gt(r$ks$p.value, 0.001)
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se, label = info)
    expect_var_near(r$sample_var, r$ref, r$n, info = info)
  }
})

test_that("Algo 1: alpha >= 4 yields acceptance >= 0.8 (Theorem 2e)", {
  # Up to alpha = 1e18, where Algorithm 1 takes the better of its two
  # proposals; z / sqrt(alpha) = 1 is about where the acceptance is least
  # (0.815 at alpha = 1e6, 0.825 at 1e17 with z / sqrt(alpha) = 0.95).  4e4
  # draws put the acceptance measured within about 0.002 of its value.
  n <- 4e4
  for (params in list(c(5, 1, 5), c(10, 1, 10), c(1e6, 1, 1e3),
                      c(1e17, 1, 3e8))) {
    set.seed(1)
    x <- mhn:::.rmhn_sun_algo1_cpp(n, params[1], params[2], params[3])
    expect_gt(n / (n + attr(x, "sun_retries")), 0.8,
              label = sprintf("acceptance at alpha = %g", params[1]))
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

test_that("Algo 1 setup matches Theorem 1 evaluated from its definition", {
  # With h(x) = x^(alpha-1) exp(-beta x^2 + gamma x), the Normal proposal
  # g1 = N(mu_opt, 1/(2 beta)) and the sqrt-Gamma proposal
  # g2(x) = 2 delta^(alpha/2) x^(alpha-1) exp(-delta x^2) / Gamma(alpha/2),
  # the envelope constants of Sun et al. (2023, Theorem 1) are
  # K_i = sup_x h(x) / g_i(x), leaving out the normalising factor
  # 2 beta^(alpha/2) / Psi of the density, which is common to both.
  # Sun et al. print f_sqrtGam without the factor 2, so that it integrates to
  # 1/2; taken literally it would put the reference low by log 2.  The optimal
  # constants are taken from their closed forms and the suprema found by
  # optimize(), so the comparison checks the selector's algebra rather than
  # restating it.  Each log ratio is concave, its maximiser is located to
  # 1e-10, and a maximum moves by the square of that; the logarithms
  # themselves are of size 100 at most here, so the reference is good to
  # 1e-13, and 1e-8 leaves a wide margin.  Algorithm 1 takes the Normal
  # proposal when K2 > K1, that is when log K1 - log K2 < 0.  Both proposals
  # occur among these points.
  pts <- list(c(5, 1, 5), c(2, 1, 0.5), c(20, 100, 30), c(1.2, 4, 0.1),
              c(50, 1, 7.1))
  for (p in pts) {
    a <- p[1]; b <- p[2]; g <- p[3]
    info <- sprintf("alpha=%g beta=%g gamma=%g", a, b, g)
    mu <- (g + sqrt(g^2 + 8 * (a - 1) * b)) / (4 * b)
    delta <- 8 * a * b^2 / (sqrt(g^2 + 8 * a * b) + g)^2
    log_h <- function(x) (a - 1) * log(x) - b * x^2 + g * x
    log_g1 <- function(x) 0.5 * log(b / pi) - b * (x - mu)^2
    log_g2 <- function(x) log(2) + (a / 2) * log(delta) + (a - 1) * log(x) -
      delta * x^2 - lgamma(a / 2)
    upper <- 100 * (mu + 1 / sqrt(b))
    sup_log <- function(f) {
      optimize(f, c(0, upper), maximum = TRUE, tol = 1e-10)$objective
    }
    ref <- sup_log(function(x) log_h(x) - log_g1(x)) -
      sup_log(function(x) log_h(x) - log_g2(x))
    d <- mhn:::.dump_sun_algo1_cpp(a, b, g)
    expect_equal(d$log_K1_minus_K2, ref, tolerance = 1e-8, info = info)
    expect_identical(d$chosen, if (ref < 0) "normal" else "sqrt_gamma",
                     info = info)
    expect_equal(d$mu_opt, mu, tolerance = 1e-12, info = info)
    expect_equal(d$delta_opt, delta, tolerance = 1e-12, info = info)
  }
})

test_that("Algo 1 constants are finite and scale with beta at extreme beta", {
  # delta_opt is a rate on the X^2 axis, so at fixed (alpha, z = gamma /
  # sqrt(beta)) it is beta times its beta = 1 value, and mu_opt is 1 / sqrt(beta)
  # times its own; the proposal choice and log K1 - log K2 do not change.
  # beta = 1e160 is past the point where 8 alpha beta^2 overflows a double, and
  # 1e-170 past the point where it underflows to zero, although delta_opt
  # itself is of order beta at both.  delta_opt and mu_opt are a few roundings
  # from their beta = 1 values, hence 1e-14.  log K1 - log K2 (about 0.28) is
  # a sum of several terms of size log(beta), about 390, that cancel; each
  # carries a rounding of about 390 eps = 1e-13, hence 1e-11.
  a <- 1.5; z <- 0.5
  d1 <- mhn:::.dump_sun_algo1_cpp(a, 1, z)
  for (b in c(1e160, 1e-170)) {
    d <- mhn:::.dump_sun_algo1_cpp(a, b, z * sqrt(b))
    info <- sprintf("beta=%g", b)
    expect_true(is.finite(d$delta_opt) && d$delta_opt > 0, info = info)
    expect_equal(d$delta_opt / b, d1$delta_opt, tolerance = 1e-14, info = info)
    expect_equal(d$mu_opt * sqrt(b), d1$mu_opt, tolerance = 1e-14, info = info)
    expect_equal(d$log_K1_minus_K2, d1$log_K1_minus_K2, tolerance = 1e-11,
                 info = info)
    expect_identical(d$chosen, d1$chosen, info = info)
  }
})

test_that("rmhn(method = 'sun') draws correctly at beta = 1e160", {
  # At beta = 1e160 Algorithm 1 must still run with the sqrt-Gamma proposal it
  # chooses at beta = 1.  sqrt(beta) X then has the beta = 1 law, so the
  # rescaled draws are tested against the reference law at beta = 1, at the
  # 0.001 level used by the other KS checks in this file.
  a <- 1.5; b <- 1e160; z <- 0.5
  set.seed(1)
  x <- expect_silent(rmhn(2000, a, b, z * sqrt(b), method = "sun"))
  expect_true(all(is.finite(x)) && all(x > 0))
  skip_on_cran()
  ks <- ks.test(sqrt(b) * x, mhn_cdf(a, 1, z))
  expect_gt(ks$p.value, 0.001)
})

test_that("Algo 1 samples the law where its setup leaves the double range", {
  # At (2, 1e300, 1e155) gamma^2 passes the largest double, and the setup
  # written with beta and gamma cannot be formed, although the law is an
  # ordinary one: that of X = Y / sqrt(beta) with Y ~ MHN(2, 1, 1e5).  It is
  # normal to far below the resolution of the test (standardised skewness
  # about 2 (alpha-1) sd^3 / m^3 = 6e-15), so sqrt(beta) X is compared with
  # N(m, 1 / ((alpha-1)/m^2 + 2)), m the mode of Y, at the 0.001 level used
  # throughout this file.
  set.seed(15)
  x <- expect_silent(rmhn(2000, 2, 1e300, 1e155))
  expect_false(anyNA(x))
  # At gamma / sqrt(beta) = 1e300, gamma^2 overflows in every form of the
  # setup, and at (2, 0.81, 1.7e308) gamma / sqrt(beta) itself does.  At
  # (5e307, 1.7e308, 1.7e154) and (1e308, 1e308, 1e154) alpha is above
  # DBL_MAX / 8, so that 8 alpha beta overflows for any beta >= 1, while the
  # mode is of order one (0.38 and 0.71).  In each case the standard
  # deviation of the law is below 1e-150 times its mode m, the positive root
  # of 2 beta x^2 - gamma x - (alpha - 1) = 0, so every draw must be m
  # rounded to double precision.  m is written here as q + sqrt(q^2 + s^2),
  # with q = gamma / (4 beta) and s^2 = (alpha - 1) / (2 beta), and the root
  # taken as max(q, s) sqrt(1 + (min(q, s) / max(q, s))^2), so that nothing
  # overflows; at the first two points it is gamma / (2 beta) to a relative
  # 1e-600.  It is formed here and in the setup with a few roundings each,
  # hence 4 eps.
  for (p in list(c(2, 1, 1e300), c(2, 0.81, 1.7e308),
                 c(5e307, 1.7e308, 1.7e154), c(1e308, 1e308, 1e154))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    set.seed(16)
    y <- expect_silent(rmhn(50, p[1], p[2], p[3]))
    q <- (p[3] / 4) / p[2]
    s <- sqrt(((p[1] - 1) / 2) / p[2])
    m <- q + max(q, s) * sqrt(1 + (min(q, s) / max(q, s))^2)
    expect_true(all(abs(y / m - 1) <= 4 * .Machine$double.eps), info = info)
  }
  # At beta = 1e-310 and gamma = 1 the mode, about gamma / (2 beta) = 5e309,
  # is beyond the largest double; no draw can represent the law, and the call
  # stops rather than returning NaN or Inf.
  expect_error(rmhn(3, 2, 1e-310, 1), "largest double precision number")
  m <- mhn_mode(2, 1, 1e5)
  s <- 1 / sqrt(1 / m^2 + 2)
  expect_gt(ks.test(sqrt(1e300) * x, "pnorm", m, s)$p.value, 0.001)
})

test_that("Algo 1 places its Normal draws at the mode to double precision", {
  # At (2, 7, gamma) with gamma / sqrt(beta) = 1e15 the doubles near the
  # mode are 0.117 standard deviations apart, and the double nearest the mode
  # lies 0.0134, 0.050 standard deviations, below it.  A Normal proposal
  # centred at that double would move the whole law by as much, which 1e5
  # draws resolve at 16 standard errors.  At (3e28, 1, 1) that offset is
  # 0.0012 standard deviations, but a sqrt-Gamma proposal, the one Theorem 1
  # of Sun et al. (2023) prefers there, is placed by R::rgamma from a rounded
  # sqrt(alpha/2 - 1/2), and moves the law by 0.026 standard deviations,
  # 8 standard errors here.  Each mode is the positive root of
  # 2 beta x^2 - gamma x - (alpha - 1) = 0, computed with Rmpfr at 300 bits
  # from these doubles and written as the double hi plus the remainder lo;
  # the law's mean exceeds its mode by about (alpha - 1) sd^4 / m^3, at most
  # 2e-15 standard deviations here.  The draws are measured from hi, a
  # subtraction that is exact for draws within a factor of two of it, and
  # their mean must lie within 5 standard errors of lo.
  refs <- list(
    list(p = c(2, 7, 2645751311064590.5), hi = 188982236504613.59,
         lo = 0.01339285714285752, sd = 0.2672612419124244),
    list(p = c(3e28, 1, 1), hi = 122474487139159.16,
         lo = 0.00059574180608350149, sd = 0.5000000000000006)
  )
  n <- 1e5
  for (r in refs) {
    set.seed(17)
    x <- rmhn(n, r$p[1], r$p[2], r$p[3])
    expect_lt(abs(mean(x - r$hi) - r$lo), 5 * r$sd / sqrt(n),
              label = sprintf("offset of the mean at alpha=%g beta=%g gamma=%g",
                              r$p[1], r$p[2], r$p[3]))
  }
  # 2^510 X, for X of the first law, is MHN(2, 7 / 2^1020, gamma / 2^510).
  # With beta below 1e-300 the setup is built for 2^k X, 2^k a power of two
  # within a factor of two of sqrt(beta).  Scaling by a power of two is exact
  # and commutes with every rounding in the setup and in the draws, so under
  # one seed the draws there must be 2^510 times those above, placed at the
  # mode as they are.  (A setup rescaled by the rounded sqrt(beta), for
  # comparison, moves this law by about 0.1 standard deviations.)
  set.seed(17)
  x1 <- rmhn(1000, 2, 7, refs[[1]]$p[3])
  set.seed(17)
  x2 <- rmhn(1000, 2, 7 * 2^-1020, refs[[1]]$p[3] * 2^-510)
  expect_identical(x2, x1 * 2^510)
})

test_that("Algo 1 draws the mode where the law is narrower than the doubles", {
  # At these shapes the standard deviation, at most 1 / sqrt(2 beta), is
  # below 1e-49 times the mode m = q + sqrt(q^2 + (alpha - 1) / (2 beta)),
  # q = gamma / (4 beta), so every draw must be m rounded to double
  # precision.  A sqrt-Gamma proposal, which Theorem 1 of Sun et al. (2023)
  # prefers at the first two points, takes a single value there, off the
  # centre of its acceptance test, and every proposal would be rejected.  At
  # (1e100, 1e200, 1e100) alpha beta is 1e300, at the bound up to which the
  # setup is formed with beta and gamma directly, and at alpha = 1.7e308
  # 8 alpha overflows.  m is formed here and in the setup with a few roundings
  # each, hence 4 eps.
  for (p in list(c(1e100, 1, 1), c(1e100, 1e200, 1e100), c(1e300, 1, 1),
                 c(1.7e308, 1, 1))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    set.seed(23)
    y <- expect_silent(rmhn(50, p[1], p[2], p[3]))
    q <- (p[3] / 4) / p[2]
    m <- q + sqrt(q^2 + ((p[1] - 1) / 2) / p[2])
    expect_true(all(abs(y / m - 1) <= 4 * .Machine$double.eps), info = info)
  }
})

test_that("Algo 1 proposal selector matches a 600-bit evaluation", {
  # log K1 - log K2 of Sun et al. (2023, Theorem 1) at the optimal constants,
  # with Psi cancelled.  The references were computed with Rmpfr at 600 bits
  # from the literal expression
  #   log 2 + log(pi)/2 + (a-1) log(sqrt(b) (a-1) / (2 b mu - g)) - (a-1)
  #   + b mu^2 - (a/2) log b - lgamma(a/2) - g^2 / (4 (b - delta))
  #   + (a/2) log delta,
  # with mu_opt and delta_opt also in 600 bits, the inputs being the doubles
  # below.  In double precision that expression is a sum of terms of size
  # a log a that cancel to order one, and from alpha of about 1e15 its
  # rounding decides the sign that picks the proposal.  Above alpha = 1e4 the
  # package evaluates a cancellation-free form, whose terms are of order one
  # and which agrees with the references to a few roundings; hence 1e-10.  Up
  # to 1e4 it keeps the terms of size a log a.  For alpha <= 100 these are at
  # most about 1e3, and at these points their rounding is at most of order
  # 1e-12 relative, well inside 1e-10.  At alpha = 1e4 they reach about 1e5
  # (lgamma(5000) = 3.8e4, (a-1) log(a-1) = 9.2e4), each carrying a rounding
  # of a few 1e-11 once the libm's own error in lgamma and log is counted, so
  # their sum can be off by about 1e-10 absolute, 6e-10 of the value 0.17
  # there; hence 2e-9 for alpha in (100, 1e4].
  refs <- list(
    list(c(2, 1, 0.5),         0.27302576245340127),
    list(c(4, 1, 2),           0.0033392718951550915),
    list(c(1.5, 1, 0.05),      0.44735741203739132),
    list(c(4, 1, 20),         -1.6678335252224028),
    list(c(100, 0.01, 3000),  -7.3148872526483224),
    list(c(1 + 1e-7, 1, 1),    0.096573584483401531),
    list(c(3, 4, 4),          -0.049171847379161381),
    list(c(1e4, 1, 50),        0.17070887726826159),
    list(c(1e6, 2, 1e3),       0.099107151432189886),
    list(c(1e15, 1, 1e8),     -0.61585005983923435),
    list(c(1e16, 1, 1e10),    -3.9122229454547988),
    list(c(1e18, 1, 1),        0.34657358992641929),
    list(c(1e20, 1, 1e14),    -8.5171932114162363),
    list(c(1e25, 1, 1e12),     0.23500181462286779)
  )
  for (r in refs) {
    p <- r[[1]]
    got <- mhn:::.dump_sun_algo1_cpp(p[1], p[2], p[3])$log_K1_minus_K2
    tol <- if (p[1] > 100 && p[1] <= 1e4) 2e-9 else 1e-10
    expect_equal(got, r[[2]], tolerance = tol,
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
})

test_that("Algo 1 samples the right law at huge alpha", {
  # At these shapes the law is normal to far below the resolution of the
  # test: its standardised skewness is about 2 (alpha-1) sd^3 / m^3, at most
  # 3e-9 here.  So the draws are compared with N(m, 1 / ((alpha-1)/m^2 +
  # 2 beta)), m the mode, by a Kolmogorov-Smirnov test at the 0.001 level used
  # throughout this file.  Both acceptance tests of Algorithm 1, as printed,
  # are differences of terms of size (alpha - 1)|X/m - 1| or gamma X that
  # cancel to order one, and so is the comparison that picks the proposal;
  # evaluated that way their rounding exceeds one at these shapes, and the
  # accepted law has a standard deviation between 1.04 and 50 times too large
  # at these four points.  The draws lie on the grid of doubles at the mode,
  # whose spacing is between 1e-8 and 9e-5 standard deviations here: it makes
  # ties, about which ks.test warns, but moves the statistic by at most the
  # normal density times that spacing, below 4e-5 against a critical value of
  # 0.0195.
  pts <- list(c(1e15, 1, 1e8), c(1e16, 1, 1e10), c(1e18, 1, 1),
              c(1e20, 1, 1e12))
  draws <- lapply(pts, function(p) {
    set.seed(3)
    mhn:::.rmhn_sun_algo1_cpp(10000, p[1], p[2], p[3])
  })
  for (x in draws) expect_false(anyNA(x))
  for (k in seq_along(pts)) {
    p <- pts[[k]]
    m <- mhn_mode(p[1], p[2], p[3])
    s <- 1 / sqrt((p[1] - 1) / m^2 + 2 * p[2])
    ks <- suppressWarnings(ks.test(draws[[k]], "pnorm", m, s))
    expect_gt(ks$p.value, 0.001,
              label = sprintf("KS p at alpha=%g beta=%g gamma=%g",
                              p[1], p[2], p[3]))
  }
})

# =====================================================================
# Algorithm 3: gamma <= 0, alpha > 0
#
# Theorem 4 of Sun et al. (2023) is stated for alpha > 1, but its part
# (a), a ratio of integrals, and l, l', l'' hold for every alpha > 0; the
# tests below that use the acceptance probability of Theorem 4a or the
# Newton step on l at alpha <= 1 rely on that.
# =====================================================================

run_sun_algo3 <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- mhn:::.rmhn_sun_algo3_cpp(n, alpha, beta, gamma)
  ref <- mhn_ref(alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    m             = attr(x, "sun_m"),
    r             = attr(x, "sun_r"),
    used_inflex   = attr(x, "sun_used_inflex"),
    retries       = attr(x, "sun_retries"),
    acceptance    = n / (n + attr(x, "sun_retries")),
    ref           = ref,
    ks            = ks.test(x, ref$cdf),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = ref$mean,
    theory_var    = ref$var
  )
}

test_that("Algo 3 alpha > 1: KS + acceptance >= 1/sqrt(2) (Theorem 4c)", {
  skip_on_cran()
  # The 1/sqrt(2) is Sun et al.'s (2023) bound (Theorem 4c, stated for
  # alpha > 1; the supplementary proof treats alpha >= 3 only), proved at
  # m_opt and m_star rather than at the Newton m used here, so the check is
  # empirical, at alpha = 1.5 and 1.05 as elsewhere.
  for (params in list(c(3, 1, -2), c(5, 1, -5), c(1.5, 1, -1),
                      c(3, 4, -4), c(5, 0.01, -0.5), c(1.05, 100, -10))) {
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
  # The proposal is then the law itself, and no proposal is rejected.
  expect_identical(as.numeric(r$retries), 0)
})

test_that("Algo 3: rejects gamma > 0", {
  expect_error(mhn:::.rmhn_sun_algo3_cpp(10, 2, 1, 1),
               regexp = "gamma <= 0")
  expect_error(mhn:::.rmhn_sun_algo3_cpp(10, 2, 1, 5),
               regexp = "gamma <= 0")
})

test_that("Algo 3 dump exposes consistent fields and ordering invariants", {
  # r, shape and rate are checked against the Gamma proposal of Sun et al.
  # (2023, Algorithm 3, steps 2-3; the shape is alpha (beta m + |gamma|) /
  # (2 beta m + |gamma|) as in their Theorem 4a, step 2 printing
  # (beta m + |gamma| alpha) / (2 beta m + |gamma|), a misprint), recomputed
  # here in R from the dumped m;
  # comparing one dumped field with another could not fail, because the setup
  # assigns them from a single expression.  The starting point m_init is
  # checked against its definition in Sun et al. (2023, Section 4.2), with
  # the closed form rescaled to beta:
  #   alpha <= 1.1:  alpha^2 / ((1 + alpha) sqrt(beta));
  #   alpha > 1.1:   1.5 lambda x_mode + (1 - 1.5 lambda) x_inflex, with
  #                  lambda = f(x_left) / (f(x_left) + f(x_inflex)) and
  #                  x_left = 2 x_mode - x_inflex, or x_inflex itself where
  #                  x_left <= 0, as at (1.5, 1, -20).
  # Here x_inflex, the inflection point right of the mode, where
  # (log f)'^2 + (log f)'' = 0, is found by uniroot to a relative 1e-14; the
  # setup's bisection stops at a relative width of 1e-12, hence 1e-10.
  for (par in list(c(3, 1, -2), c(5, 0.25, -5), c(1.05, 4, -0.5),
                   c(100, 1, -10), c(1.5, 1, -20))) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    d <- mhn:::.dump_sun_algo3_cpp(alpha, beta, gamma)
    g_abs <- abs(gamma)
    A <- beta * d$m + g_abs
    B <- 2 * beta * d$m + g_abs
    info <- sprintf("alpha=%g beta=%g gamma=%g", alpha, beta, gamma)

    expect_gt(d$m, 0)
    expect_gt(d$r, 0); expect_lt(d$r, 1)
    expect_equal(d$r, A / B, tolerance = 1e-12, info = info)
    expect_equal(d$shape, alpha * A / B, tolerance = 1e-12, info = info)
    expect_equal(d$rate, d$m * (beta * d$m + g_abs), tolerance = 1e-12,
                 info = info)
    expect_identical(d$used_inflex_heuristic, alpha > 1.1, info = info)
    if (alpha <= 1.1) {
      m_init <- alpha^2 / ((1 + alpha) * sqrt(beta))
    } else {
      x_mode <- 2 * (alpha - 1) / (sqrt(gamma^2 + 8 * beta * (alpha - 1)) -
                                     gamma)
      curv <- function(x) {
        ((alpha - 1) / x - 2 * beta * x + gamma)^2 - (alpha - 1) / x^2 -
          2 * beta
      }
      x_inflex <- uniroot(curv, c(x_mode * (1 + 1e-8), 2 * x_mode + 1),
                          extendInt = "upX", tol = 1e-14 * x_mode)$root
      x_left <- 2 * x_mode - x_inflex
      m_init <- if (x_left > 0) {
        log_f <- function(x) (alpha - 1) * log(x) - beta * x^2 + gamma * x
        lambda <- 1 / (1 + exp(log_f(x_inflex) - log_f(x_left)))
        1.5 * lambda * x_mode + (1 - 1.5 * lambda) * x_inflex
      } else {
        x_inflex
      }
    }
    expect_equal(d$m_init, m_init, tolerance = 1e-10, info = info)
  }
})

test_that("Algo 3 setup is scale-equivariant in beta", {
  # sqrt(beta) X ~ MHN(alpha, 1, z) with z = gamma / sqrt(beta), and the
  # acceptance probability of Algorithm 3 (Sun et al. 2023, Theorem 4a)
  # depends on (beta, gamma, m) only through z and sqrt(beta) m.  So at a
  # fixed (alpha, z) the matching point times sqrt(beta), and r, the Gamma
  # shape and the rate, which are functions of (alpha, z, sqrt(beta) m), must
  # not depend on beta.  The alpha > 1.1 cases go through the inflection-point
  # heuristic, whose bisection has to measure its bracket and its stopping
  # width relative to the density's own scale; (0.7, -3) goes through the
  # alpha <= 1.1 closed form.  Beyond beta of about 1e154 (and below 1e-154)
  # the terms of l'(m) and l''(m) leave the double range unless they are
  # formed at beta = 1.  The only difference left is the rounding of
  # z sqrt(beta) / sqrt(beta), which can move the bisection's last step, a
  # relative 1e-12 before the Newton step shrinks it; hence 1e-10.
  for (p in list(c(1.5, -1), c(5, -20), c(20, -20), c(3, -0.5), c(0.7, -3))) {
    a <- p[1]; z <- p[2]
    d1 <- mhn:::.dump_sun_algo3_cpp(a, 1, z)
    for (b in 10^c(-200, -30, 30, 100, 200)) {
      d <- mhn:::.dump_sun_algo3_cpp(a, b, z * sqrt(b))
      info <- sprintf("alpha=%g z=%g beta=%g", a, z, b)
      expect_equal(d$m * sqrt(b), d1$m, tolerance = 1e-10, info = info)
      expect_equal(d$m_init * sqrt(b), d1$m_init, tolerance = 1e-10,
                   info = info)
      expect_equal(c(d$r, d$shape, d$rate), c(d1$r, d1$shape, d1$rate),
                   tolerance = 1e-10, info = info)
    }
  }
})

test_that("Algo 3 keeps its acceptance at beta = 1e28", {
  # The law at (alpha, beta, z sqrt(beta)) is the beta = 1 law scaled by
  # 1 / sqrt(beta), so the sampler has to perform as it does at beta = 1: no
  # draw may exhaust the retry budget of 1000 proposals, which would come back
  # as NaN with a warning.  A matching point that has lost its scale leaves
  # an acceptance of 6e-4 at (5, z = -20) and 2e-9 at (20, z = -20).  The
  # rescaled draws are then compared with the beta = 1 law.
  b <- 1e28
  for (p in list(c(5, -20), c(20, -20))) {
    a <- p[1]; z <- p[2]
    set.seed(11)
    x <- expect_silent(rmhn(2000, a, b, z * sqrt(b), method = "sun"))
    expect_false(anyNA(x))
  }
  skip_on_cran()
  set.seed(12)
  x <- rmhn(2000, 20, b, -20 * sqrt(b), method = "sun")
  expect_gt(ks.test(sqrt(b) * x, mhn_cdf(20, 1, -20))$p.value, 0.001)
})

test_that("Algo 3 samples the Gamma limit at a huge tilt", {
  # With t = |gamma| / sqrt(beta) the kernel is x^(alpha-1) exp(-|gamma| x)
  # times exp(-beta x^2), and beta x^2 = (|gamma| x)^2 / t^2.  Here t is from
  # 1e157 to 1e200, so for |gamma| x up to 1e5 that factor is 1 to within
  # 1e-300, and the law is Gamma(alpha, rate |gamma|) to double precision.
  # Algorithm 3's setup must then work in a unit in which neither t^2 nor
  # 1 / (sqrt(beta) m)^2 leaves the double range.  The cases include a
  # subnormal beta, both sides of alpha = 1.1 (the two m_init formulas), and
  # the default path, which takes Algorithm 3 past a tilt of 1e150 whatever
  # the number of draws per setup, up to the shape 2^59.  No draw may exhaust
  # the retry budget, which would come back as NaN with a warning, and
  # |gamma| x must follow Gamma(alpha, 1) at the 0.001 level used throughout
  # this file.
  cases <- list(c(20, 1e-300, -1e10), c(0.5, 1e-300, -1e10),
                c(3, 1, -1e200), c(1.5, 1e-310, -100))
  methods <- c("auto", "sun", "sun", "sun")
  set.seed(13)
  expect_false(anyNA(expect_silent(rmhn(20, 0.5, 1e-300, -1e10))))
  draws <- vector("list", length(cases))
  for (k in seq_along(cases)) {
    p <- cases[[k]]
    set.seed(14)
    draws[[k]] <- expect_silent(rmhn(2000, p[1], p[2], p[3],
                                     method = methods[k]))
    expect_false(anyNA(draws[[k]]),
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
  for (k in seq_along(cases)) {
    p <- cases[[k]]
    expect_gt(ks.test(draws[[k]] * abs(p[3]), "pgamma", p[1])$p.value, 0.001,
              label = sprintf("KS p at alpha=%g beta=%g gamma=%g",
                              p[1], p[2], p[3]))
  }
})

test_that("Algo 3 samples the Gamma limit where the tilt leaves the double range", {
  # |gamma| / sqrt(beta) overflows to Inf at the first three points, and the
  # default path takes Algorithm 3 at the fourth, as it does past a tilt of
  # 1e150.  The law is Gamma(alpha, rate |gamma|) to double precision (see the
  # test above), and the setup has to be formed in units of alpha / |gamma|
  # without the tilt.
  # At the last two points the tilt is finite, 1.7e308, but a matching point
  # of order one, as Sun et al.'s closed form for alpha <= 1.1 gives, would
  # put the rate of the Gamma proposal near the largest double, where most of
  # its draws fall below DBL_MIN and are drawn again from their conditional
  # law there (see gamma_precision_floor in mhn_sun.h).  The Gibbs-style call
  # rmhn(1, ...), which the default path answers with Algorithm 3, is
  # included.  No draw may exhaust the retry budget, and |gamma| x must follow
  # Gamma(alpha, 1) at the 0.001 level used throughout this file.
  cases <- list(c(2, 1e-300, -1e300), c(0.5, 1e-300, -1e300),
                c(2, 5e-324, -1e200), c(20, 1e-300, -1e300),
                c(1.05, 1, -1.7e308), c(0.5, 1, -1.7e308))
  methods <- c("sun", "sun", "sun", "auto", "sun", "sun")
  draws <- vector("list", length(cases) + 1)
  for (k in seq_along(cases)) {
    p <- cases[[k]]
    set.seed(18)
    draws[[k]] <- expect_silent(rmhn(2000, p[1], p[2], p[3],
                                     method = methods[k]))
    expect_false(anyNA(draws[[k]]),
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
  set.seed(19)
  draws[[length(cases) + 1]] <- expect_silent(
    vapply(1:2000, function(i) rmhn(1, 1.05, 1, -1.7e308), numeric(1)))
  cases <- c(cases, list(c(1.05, 1, -1.7e308)))
  for (k in seq_along(cases)) {
    p <- cases[[k]]
    expect_gt(ks.test(draws[[k]] * abs(p[3]), "pgamma", p[1])$p.value, 0.001,
              label = sprintf("KS p at alpha=%g beta=%g gamma=%g (case %d)",
                              p[1], p[2], p[3], k))
  }
})

test_that("Algo 3 returns the zeros of a tiny shape at a huge tilt", {
  # Gamma(1e-8, rate 1e308) puts all but 3.6e-7 of its mass below 2^-1075,
  # where a double rounds to 0: P(X > x) is 1 - (1e308 x)^1e-8 to first order.
  # So the draws must be zeros, not NaN from proposals that are all rejected,
  # with method = "sun" and on the default path, which takes Algorithm 3 at
  # this tilt (-gamma / sqrt(beta) > 1e150) for any n.
  for (n in c(200, 50)) {
    m <- if (n == 200) "sun" else "auto"
    set.seed(20)
    x <- expect_silent(rmhn(n, 1e-8, 1, -1e308, method = m))
    expect_identical(x, numeric(n), label = paste("method", m))
  }
})

test_that("Algo 3 samples the right law at huge alpha", {
  # The acceptance exponent of Algorithm 3, rate T - beta X^2 - |gamma| X, is
  # a sum of three terms of size rate, about alpha, that cancel to order one;
  # evaluated that way its rounding reaches order one near alpha = 1e16, and
  # the accepted law is distorted.  At these shapes the law is normal to
  # far below the resolution of the test (standardised skewness at most
  # 2 (alpha-1) sd^3 / m^3 < 1e-7), so the draws are compared with
  # N(m, 1 / ((alpha-1)/m^2 + 2 beta)), m the mode, at the 0.001 level.  The
  # acceptance must also stay above the 0.9 asserted by the extreme-alpha
  # guard below (0.99 is asserted here).  The draws lie on the grid of
  # doubles at the mode, whose spacing is between 3e-8 and 2.4e-7 standard
  # deviations here (the largest for the third of the draws at the last point
  # that lie above 2^29): it makes ties, about which ks.test warns, but
  # moves the statistic by at most the normal
  # density times that spacing, 0.4 x 2.4e-7, below 1e-7, against a critical
  # value of 0.0195.  The last point is the largest shape at which
  # Algorithm 3 is used.
  for (p in list(c(1e16, 1, -1), c(1e17, 1, -1e10), c(5e17, 2, -3),
                 c(2^59, 1, -1))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    set.seed(8)
    x <- mhn:::.rmhn_sun_algo3_cpp(10000, p[1], p[2], p[3])
    expect_false(anyNA(x), info = info)
    # The acceptance measured is 0.9994 to 1 here; 0.99 is some 40 of its
    # standard errors below.
    expect_gt(10000 / (10000 + attr(x, "sun_retries")), 0.99, label = info)
    m <- mhn_mode(p[1], p[2], p[3])
    s <- 1 / sqrt((p[1] - 1) / m^2 + 2 * p[2])
    ks <- suppressWarnings(ks.test(x, "pnorm", m, s))
    expect_gt(ks$p.value, 0.001, label = info)
  }
})

test_that("Algo 3 is not used above a shape of 2^59", {
  # Algorithm 3 proposes a Gamma variate of shape up to alpha, and R's rgamma
  # builds one of shape a from sqrt(a - 1/2) rounded to double precision and
  # multiplies it by the scale 1 / rate, rounded too; together these move
  # every draw alike by up to about 3 x 2^-53 sqrt(a) standard deviations:
  # 2.5e-7 at a = 2^59, but a whole standard deviation near a = 2^103, where
  # the law is still wider than the spacing of the doubles.  Above 2^59,
  # method = "sun" stops with that reason, and the default method draws from
  # RTDR, which places each draw at the mode computed to about twice double
  # precision, for a single draw as for many.
  for (p in list(c(2^59 * (1 + 2^-52), 1, -1), c(1e20, 1, -1e10),
                 c(1e250, 1, -1), c(1e160, 1, -1e300))) {
    expect_error(rmhn(3, p[1], p[2], p[3], method = "sun"),
                 "Algorithm 3 is not used",
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
  expect_silent(rmhn(3, 2^59, 1, -1, method = "sun"))
  # At (2^101, 1/2, -2^50) the mode m solves m^2 + 2^50 m - (2^101 - 1) = 0,
  # so m = 2^50 - 1/(3 * 2^50) to first order, and the law is normal about it
  # with standard deviation 1/sqrt(3) to far below the resolution of the test
  # (standardised skewness below 1e-15).  x - 2^50 is exact, and the mean of the
  # draws must lie within five standard errors of 2^50.  Drawn from the Gamma
  # proposal, it lies 0.34 standard deviations below.
  for (n in c(1L, 20000L)) {
    set.seed(9)
    x <- if (n == 1L) vapply(1:20000, function(i) rmhn(1, 2^101, 0.5, -2^50), 0)
         else rmhn(n, 2^101, 0.5, -2^50)
    expect_lt(abs(mean(x - 2^50)), 5 / sqrt(3 * 20000),
              label = sprintf("mean offset, %d draw(s) per call", n))
    expect_lt(abs(sd(x) * sqrt(3) - 1), 0.03)
  }
})

test_that("Algo 3 samples a shape below 1e-150", {
  # For a shape this small the law puts all but a fraction of order 1e3 alpha
  # of its mass below 2^-1075, where a double rounds to 0, whatever the tilt,
  # so every draw must be 0.  The Newton step of Sun et al. (2023,
  # Supplementary Section 2.17) needs the digamma and trigamma functions at a
  # shape between alpha / 2 and alpha, where they overflow (trigamma below
  # 7.46e-155, digamma below 5.6e-309); and the closed-form starting point
  # alpha^2 / (1 + alpha) makes the rate of the Gamma proposal subnormal or 0.
  # Neither may surface as an error or as NaN draws, on the default path,
  # which takes Algorithm 3 for fewer than 25 draws at alpha < 0.1, or with
  # method = "sun".  The matching point is raised until the rate
  # u (u + |gamma| / sqrt(beta)) is normal, and from a tilt of about 1.8e16
  # that point, about 4.5e-308 / |gamma|, lies below half the smallest
  # positive double and rounds to 0; the tilts from -1e20 to -1e149 cover
  # that band, up to 1e150, past which Algorithm 3 changes its unit of length.
  for (a in c(1e-153, 1e-200, 1e-300, 5e-324)) {
    for (g in c(-2e-8, -2, -1e10, -1e20, -1e100, -1e149, -1e200)) {
      for (m in c("auto", "sun")) {
        set.seed(10)
        x <- expect_silent(rmhn(5, a, 1, g, method = m))
        expect_identical(x, numeric(5),
                         label = sprintf("alpha=%g gamma=%g %s", a, g, m))
      }
    }
  }
  # The same band reached through a tiny beta: -gamma / sqrt(beta) = 1e100.
  set.seed(10)
  x <- expect_silent(rmhn(5, 1e-200, 1e-300, -1e-50))
  expect_identical(x, numeric(5))
  # For 25 draws and more the default method takes RTDR at alpha < 0.1,
  # which builds its envelope on the log axis about the mode of
  # log(sqrt(beta) X), the log of about alpha sqrt(beta) / |gamma| here.  At
  # these points that number underflows to 0, and the default method must
  # fall back on Algorithm 3 rather than stop.
  for (p in list(c(5e-324, 1, -2), c(1e-315, 1, -1e10), c(1e-310, 4, -2e15),
                 c(1e-300, 1, -1e100))) {
    set.seed(11)
    x <- expect_silent(rmhn(200, p[1], p[2], p[3]))
    expect_identical(x, numeric(200),
                     label = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
})

test_that("Algo 3 draws the far left tail at a small shape", {
  # At alpha = 0.01 the Gamma proposal of Algorithm 3 has a shape near 0.005,
  # and R's rgamma then returns its unit-scale variate as exactly 0 about
  # 2.4% of the time; X = m T^r would turn those into zeros, although the law
  # puts only 5.8e-4 of its mass below the smallest positive double.  The mass
  # below 1e-200 is P(X < 1e-200) = x^alpha / (alpha Z) at x = 1e-200, with
  # Z = int_0^Inf t^(alpha-1) exp(-t^2 + gamma t) dt
  #   = (1/2) sum_k Gamma((alpha + k)/2) gamma^k / k!,
  # the exponential being 1 to within 1e-205 below x.  In 256-bit Rmpfr
  # arithmetic this is 0.0100286976449 at (0.01, 1, -1e-5).  The fraction
  # of 2e5 draws below 1e-200 must lie within 5 standard errors of it.
  n <- 2e5
  ref <- 0.0100286976449
  set.seed(5)
  x <- rmhn(n, 0.01, 1, -1e-5, method = "sun")
  expect_lt(abs(mean(x < 1e-200) - ref), 5 * sqrt(ref * (1 - ref) / n))
  # Exact zeros: X underflows to 0 below 2^-1075, which carries 5.8e-4 of
  # the mass (same computation); 5 standard errors is 2.7e-4.
  p0 <- 5.82333915701e-4
  expect_lt(abs(mean(x == 0) - p0), 5 * sqrt(p0 * (1 - p0) / n))
  # The same law scaled by 1e150 (beta = 1e-300, gamma = -1e-155).  Its
  # matching point m is 1e150 times larger, so m T^r is representable where
  # T^r is not.  The fraction below 1e-50 is the reference above, and the mass
  # below 2^-1075 is p0 (1e-150)^0.01 = 1.84e-5: 3.7 zeros are expected
  # in 2e5 draws, and more than 20 has Poisson probability 5e-10.
  set.seed(6)
  y <- rmhn(n, 0.01, 1e-300, -1e-155, method = "sun")
  expect_lt(abs(mean(y < 1e-50) - ref), 5 * sqrt(ref * (1 - ref) / n))
  expect_lte(sum(y == 0), 20)
})

# =====================================================================
# Newton-step refinement: uses l(m) = log(A_neg(m)) from Sun et al.
# (2023, Supplementary Section 2.17) as the target.  Verify
# l(m_recommend) >= l(m_init), since Newton on l targets the maximum
# of l (= maximum acceptance probability).  The implementation lives in
# src/mhn_sun.cpp::newton_refine_m().
# =====================================================================
test_that("Algo 3 Newton refinement does not lower l(m)", {
  # l(m) of Sun et al. (2023, Supplementary Section 2.17) without its log Psi
  # term, which does not depend on m and so cannot change the comparison
  # below; the C++ Newton step uses only l' and l''.
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

test_that("Algo 3 takes one Newton step on l(m) from m_init", {
  skip_if_not_installed("Rmpfr")
  # Sun et al. (2023, Section 4.2) recommend m_recommend, one Newton-Raphson
  # step on l(m) = log A_neg(m) from m_init (Supplementary Section 2.17).  The
  # step is recomputed here from l itself, differentiated numerically in
  # 300-bit arithmetic with h = 1e-30, so that it shares neither of the closed
  # forms of l' and l'' with the C++ code (Sun et al. print l'' with a sign
  # missing).  The two agree to a unit in the last place; 1e-12 leaves room
  # for the digamma and lgamma of other platforms.
  l_mp <- function(m, alpha, beta, gamma_abs) {
    A <- beta * m + gamma_abs
    B <- 2 * beta * m + gamma_abs
    C <- beta * m^2 + m * gamma_abs
    sh <- alpha * A / B
    sh * log(C) + log(B) - log(Rmpfr::mpfr(2, 300)) -
      (alpha / 2) * log(beta * m^2) - lgamma(sh) - log(A)
  }
  for (params in list(c(3, -2), c(5, -5), c(100, -10), c(10000, -100),
                      c(1.5, -0.5), c(0.5, -3))) {
    alpha <- params[1]; ga <- abs(params[2])
    info <- sprintf("alpha=%g, gamma=%g", alpha, params[2])
    d <- mhn:::.dump_sun_algo3_cpp(alpha, 1, params[2])
    m0 <- Rmpfr::mpfr(d$m_init, 300)
    h <- Rmpfr::mpfr(1e-30, 300)
    lp <- (l_mp(m0 + h, alpha, 1, ga) - l_mp(m0 - h, alpha, 1, ga)) / (2 * h)
    lpp <- (l_mp(m0 + h, alpha, 1, ga) - 2 * l_mp(m0, alpha, 1, ga) +
              l_mp(m0 - h, alpha, 1, ga)) / h^2
    expected <- Rmpfr::asNumeric(m0 - lp / lpp)
    expect_false(d$m == d$m_init, info = info)
    expect_lt(abs(d$m / expected - 1), 1e-12, label = info)
  }
})

# Guard for extreme alpha (1000 and 10000): with the one-step l(m)-based
# Newton matching point the acceptance stays above the 1/sqrt(2) level of
# Sun et al. (2023, Theorem 4c), which the theorem proves at m_opt and
# m_star rather than at this m.
test_that("Algo 3 keeps a high acceptance at large alpha", {
  for (params in list(c(1000, -100), c(1000, -10), c(10000, -100), c(10000, -10))) {
    alpha <- params[1]; gamma <- params[2]
    set.seed(42)
    x <- mhn:::.rmhn_sun_algo3_cpp(2000, alpha, 1, gamma)
    acc <- 2000 / (2000 + attr(x, "sun_retries"))
    # Sun et al. (2023, Theorem 4c) bound the acceptance by 1/sqrt(2) at the
    # matching points m_opt and m_star, not at the one-step Newton m used
    # here, so the bound below is empirical: the acceptance measured is 0.94
    # to 0.98 at these points, and 0.9 is some five of its standard errors
    # (about 0.005 at 2000 draws) below the least of them.
    expect_gt(acc, 0.9, label = sprintf("alpha=%g, gamma=%g", alpha, gamma))
    expect_equal(sum(is.na(x)), 0L)
  }
})
