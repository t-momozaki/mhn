# Tests for special-case detection and density dispatch.
#
# Special cases (Sun et al. 2023, Lemma 6):
#   gamma = 0        -> sqrt-Gamma         (Lemma 6a)
#   alpha = 1        -> truncated normal   (Lemma 6b)
#   alpha = 1, g = 0 -> half-normal        (Lemma 6c)

# ============================================================
# 1. Detection predicates
# ============================================================
#
# These are the tests rmhn routes on: it draws from the truncated normal for
# a shape within sqrt(eps) of 1 and from the square root of a Gamma variate
# for a tilt gamma / sqrt(beta) within sqrt(eps) of 0.  dmhn, pmhn and qmhn
# take those laws only at alpha = 1 and gamma = 0 exactly (section 9).

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

test_that("the special-case thresholds sit at sqrt(eps), on gamma / sqrt(beta)", {
  # |alpha - 1| < e and |gamma| < e sqrt(beta), e = sqrt(.Machine$double.eps).
  # The alpha points are 1 +- 0.99 e and 1 +- 1.01 e: 1 + e (1 - 1e-12) would
  # round to 1 + e itself.
  e <- sqrt(.Machine$double.eps)
  expect_true(mhn:::.is_sqrt_gamma(e * (1 - 1e-12)))
  expect_false(mhn:::.is_sqrt_gamma(e))
  for (b in c(1e-6, 1, 1e6)) {
    info <- sprintf("beta=%g", b)
    expect_true(mhn:::.is_sqrt_gamma(0.99 * e * sqrt(b), b), info = info)
    expect_true(mhn:::.is_sqrt_gamma(-0.99 * e * sqrt(b), b), info = info)
    expect_false(mhn:::.is_sqrt_gamma(1.01 * e * sqrt(b), b), info = info)
    expect_false(mhn:::.is_sqrt_gamma(-1.01 * e * sqrt(b), b), info = info)
  }
  for (s in c(-1, 1)) {
    expect_true(mhn:::.is_truncated_normal(1 + s * 0.99 * e))
    expect_false(mhn:::.is_truncated_normal(1 + s * 1.01 * e))
  }
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
    # Each element against its own size: exp(-25 beta) at x = 5, the
    # densities span 5.4 decades at beta = 0.5, 10.9 at 1 and 54 at 5.
    expect_close(result, expected, 1e-12, info = sprintf("beta=%.1f", beta))
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
    # Each side carries a rounding of about eps times the size of the terms
    # of log f, below 3e-14 here; 1e-12 is some 40 times that.
    expected <- dnorm(x, mean = mu, sd = sigma) / pnorm(mu / sigma)
    result <- dmhn(x, alpha = 1, beta = b, gamma = g)
    expect_close(result, expected, 1e-12,
                 info = sprintf("beta=%.1f, gamma=%.1f", b, g))
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
    # sqrt-Gamma density: 2 x * dgamma(x^2, alpha/2, rate = beta).  The
    # rounding of each side is below 4e-14 here, and 1e-12 some 30 times it.
    expected <- dgamma(x^2, shape = a / 2, rate = b) * 2 * x
    result <- dmhn(x, alpha = a, beta = b, gamma = 0)
    expect_close(result, expected, 1e-12,
                 info = sprintf("alpha=%.1f, beta=%.1f", a, b))
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

test_that("dmhn at alpha = 1, gamma = 0 agrees with three reference formulas", {
  # The half-normal, truncated-normal and sqrt-Gamma densities.  dmhn takes
  # the sqrt-Gamma route here (gamma = 0 is checked first).  Each side
  # carries a rounding of about eps times the size of the terms of log f,
  # beta x^2 up to 80 here; 1e-12 is some 50 times that.
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

    expect_close(d_mhn, d_hn, 1e-12, info = sprintf("mhn vs HN, beta=%.1f", beta))
    expect_close(d_mhn, d_tn, 1e-12, info = sprintf("mhn vs TN, beta=%.1f", beta))
    expect_close(d_mhn, d_sg, 1e-12,
                 info = sprintf("mhn vs sqGamma, beta=%.1f", beta))
  }
})

# ============================================================
# 6. sqrt-Gamma CDF and quantile where beta x^2 underflows
# ============================================================

# Largest elementwise relative error.  The references in sections 6 and 7 are
# compared through this rather than expect_equal(): testthat applies its
# tolerance in absolute terms once the reference is smaller than the
# tolerance, which would pass 0 for 1e-200, and it averages over a vector.
rel_err <- function(x, ref) max(abs(x / ref - 1))

test_that("pmhn with gamma = 0 stays exact where beta q^2 underflows", {
  # F(q) = P(alpha/2, t^2) with t = sqrt(beta) q (Lemma 6a).  Squaring
  # q = 1e-170 gives 1e-340, below the smallest double, so a CDF that forms
  # q^2 first sees q = 0.  Once t^2 < DBL_MIN, exp(-t^2) = 1 and P(a, t^2) =
  # t^(2a) / Gamma(a + 1) hold exactly in double precision, and F is an
  # ordinary number: at alpha = 0.01 about 2 percent of the mass lies below
  # 1e-170.  References: the incomplete-gamma series in Rmpfr at 512 bits.
  # Tolerance: a few units in the last place of F, log F and their complements.
  expect_lt(rel_err(pmhn(1e-170, 0.01, 1, 0), 0.020009880717994675), 1e-14)
  expect_lt(rel_err(pmhn(1e-170, 0.01, 1, 0, lower.tail = FALSE),
                    0.9799901192820053), 1e-14)
  expect_lt(rel_err(pmhn(1e-170, 0.01, 1, 0, log.p = TRUE),
                    -3.9115290915239691), 1e-14)
  expect_lt(rel_err(pmhn(1e-170, 0.01, 1, 0, lower.tail = FALSE, log.p = TRUE),
                    -0.020212789733647523), 1e-14)
  # Half-normal: F(q) = (2 / sqrt(pi)) q to every digit at this q.
  expect_lt(rel_err(pmhn(1e-170), 1.1283791670955121e-170), 1e-14)
  # The threshold is on t, not on q: here t = 1e-13 is ordinary while q^2
  # underflows.
  expect_lt(rel_err(pmhn(1e-163, 0.9, 1e300, 0, log.p = TRUE),
                    -26.818824997290029), 1e-14)
})

test_that("qmhn with gamma = 0 recovers quantiles whose square underflows", {
  # Where beta x^2 < DBL_MIN the quantile is explicit,
  #   x = (p Gamma(a + 1))^(1 / (2a)) / sqrt(beta),  a = alpha / 2,
  # and the references are that formula in Rmpfr at 256 bits (the 512-bit
  # series CDF at the first one returns 0.01 to 5e-16).  At alpha = 0.01 the
  # exponent 1/(2a) = 100 multiplies the relative rounding of p, or of log p,
  # by 100, hence the tolerance.  A reference below 1e-300 is written as
  # m * 10^e: R's parser can misread a 17-digit literal that small (by 2e-7
  # where long double is double, as on arm64 macOS), while 10^e is taken by
  # pow() to the last unit.
  expect_lt(rel_err(qmhn(0.01, 0.01, 1, 0), 7.508445370342333e-201), 1e-12)
  expect_lt(rel_err(qmhn(1e-3, 0.01, 1, 0), 7.5084453703423695 * 10^-301),
            1e-12)
  expect_lt(rel_err(qmhn(1e-300, 1.7, 1, 0), 3.2743496750577213e-177), 1e-12)
  expect_lt(rel_err(qmhn(1e-170), 8.8622692545275776e-171), 1e-12)
  # 1e-160 squares into the subnormal band, where x^2 keeps only a few digits.
  expect_lt(rel_err(qmhn(1e-160), 8.8622692545275778e-161), 1e-12)
  # Every tail flag leads to the same point.
  expect_lt(rel_err(qmhn(log(0.01), 0.01, 1, 0, log.p = TRUE),
                    7.508445370342333e-201), 1e-12)
  expect_lt(rel_err(qmhn(0.99, 0.01, 1, 0, lower.tail = FALSE),
                    7.508445370342333e-201), 1e-12)
  expect_lt(rel_err(qmhn(log(0.99), 0.01, 1, 0, lower.tail = FALSE,
                         log.p = TRUE),
                    7.508445370342333e-201), 1e-12)
})

# ============================================================
# 7. Truncated normal in its small lower tail (alpha = 1)
# ============================================================

test_that("the alpha = 1 CDF keeps its relative accuracy in the small lower tail", {
  # With z0 = -mu/sigma and d = q/sigma, F(q) = lambda(z0) int_0^d
  # exp(-z0 s - s^2/2) ds, lambda the normal hazard: a product of positive
  # factors.  Written as a difference of two log tail probabilities, F loses
  # about eps L / F of itself, L = |log Phibar(z0)| (4.1 at gamma = -3): 3e-3
  # at F = 3.5e-13 and all of it at 3.5e-16.
  # References: (erfc(z0/sqrt2) - erfc(zq/sqrt2)) / erfc(z0/sqrt2) in
  # Rmpfr at 2048 bits.  The results are within a few units in the last
  # place (0 to 4.4e-16); the tolerance, 1e-13, leaves room for the
  # mathematical libraries of other platforms.
  q <- c(1e-8, 1e-13, 1e-16)
  ref <- list(
    "-3" = c(3.508800749147516e-08, 3.5088008017790013e-13,
             3.508800801779527e-16),
    "1"  = c(5.7795636563504452e-09, 5.7795636274529168e-14,
             5.7795636274526271e-17),
    "3"  = c(6.0490300440336692e-10, 6.0490299532991264e-15,
             6.0490299532982202e-18))
  for (g in names(ref)) {
    gamma <- as.numeric(g)
    expect_lt(rel_err(pmhn(q, 1, 1, gamma), ref[[g]]), 1e-13,
              label = sprintf("gamma=%s", g))
    expect_lt(rel_err(pmhn(q, 1, 1, gamma, log.p = TRUE), log(ref[[g]])),
              1e-13, label = sprintf("gamma=%s log.p", g))
    # Recycled parameters take the same per-element function.
    expect_lt(rel_err(pmhn(q, c(1, 1, 1), 1, gamma), ref[[g]]), 1e-13,
              label = sprintf("gamma=%s recycled", g))
  }
  # z0 >= 3, where the upper tail is taken in its erfcx form.
  expect_lt(rel_err(pmhn(1e-17, 1, 1, -10), 1.0192700054863682e-16), 1e-13)
  expect_lt(rel_err(pmhn(1e-17, 1, 1, -60), 6.003329639875621e-16), 1e-13)
  expect_lt(rel_err(pmhn(1e-20, 1, 1, -30, log.p = TRUE),
                    -42.648294468439282), 1e-14)
})

test_that("the alpha = 1 quantile keeps its relative accuracy in the small lower tail", {
  # The quantile solves lambda(z0) I(d) = p for d directly; recovering it as
  # z - z0 from a z computed through a tail probability of order one loses
  # about eps (1 + |z0|) / d of it.  References: Newton's method on the
  # 2048-bit Rmpfr CDF above, at p = exp(log p) with log p the double given.
  # Tolerance: the returned quantiles are accurate to a few 1e-15; the margin
  # covers the rounding of log(1e-300) and of the p passed on each line.
  expect_lt(rel_err(qmhn(1e-12, 1, 1, -3), 2.8499765489487667e-13), 1e-12)
  expect_lt(rel_err(qmhn(1e-16, 1, 1, -3), 2.849976548947536e-17), 1e-12)
  expect_lt(rel_err(qmhn(log(1e-300), 1, 1, -3, log.p = TRUE),
                    2.8499765489476131 * 10^-301), 1e-12)
  expect_lt(rel_err(qmhn(1e-20, 1, 1, 1), 1.7302344337036988e-20), 1e-12)
  expect_lt(rel_err(qmhn(1e-20, 1, 1, -10), 9.8109430731538723e-22), 1e-12)
  expect_lt(rel_err(qmhn(1e-20, 1, 1, -60), 1.6657422796805094e-22), 1e-12)
  expect_lt(rel_err(qmhn(c(1e-12, 1e-16), 1, c(1, 1), -3),
                    c(2.8499765489487667e-13, 2.849976548947536e-17)), 1e-12)
  # A lower tail that exists only on the log scale: F(1.797) = exp(-800) at
  # gamma = 60, where exp(-800) itself underflows.
  expect_lt(rel_err(qmhn(-800, 1, 1, 60, log.p = TRUE), 1.797261814312614),
            1e-12)
})

test_that("the alpha = 1 quantile for z0 >= 3 inverts the upper tail", {
  # At beta = 1, gamma = -10 the truncation point is z0 = 7.07, and a p that is
  # not small solves log Phibar(z0 + t) - log Phibar(z0) = log(1 - p): up to
  # z0 = 30 through qnorm of Phibar(z0) (1 - p), refined by one Newton step
  # on the erfcx form of that equation, and beyond by Newton's method on the
  # erfcx form alone; its slope is minus the normal hazard.  Either leaves a
  # few eps of t here.
  # References: Newton's method on the 2048-bit Rmpfr CDF.
  ref <- c(0.034875973022088669, 0.067564378135163486, 0.22118919889674787)
  expect_lt(rel_err(qmhn(c(0.3, 0.5, 0.9), 1, 1, -10), ref), 1e-14)
  expect_lt(rel_err(qmhn(c(0.3, 0.5, 0.9), c(1, 1, 1), 1, -10), ref), 1e-14)
  expect_lt(rel_err(qmhn(1e-30, 1, 1, -10, lower.tail = FALSE),
                    4.6660321637506863), 1e-14)
})

test_that("the alpha = 1 quantile past z0 = 30 is accurate at ordinary p", {
  # Past z0 = 30 a p outside the small lower tail is solved by Newton's
  # method on the erfcx form of log Phibar(z0 + t) - log Phibar(z0) =
  # log(1 - p) alone, from the root L / z0 of its linear part,
  # L = -log(1 - p), up to L = 1e-3 z0^2 and from the root of its quadratic
  # part beyond.  At beta = 1 the truncation point is z0 = -gamma / sqrt(2):
  # 100, 3000 and 35 here, the second past z0 = 1000, where the two
  # logarithms of erfcx are taken as the log of their quotient.  At z0 = 100
  # p = 0.001, where h = T (z0 + T) = 1e-3 <= 1/4 (T = p / lambda(z0)), is
  # on the small-tail series; the other p there and at 3000 start from the
  # linear root, the upper tail 1e-30 from the quadratic one.  At z0 = 35,
  # where 1e-3 z0^2 is 1.2, p = 0.5 starts from the linear root and 0.9,
  # 0.999 and the upper tails 1e-100 and 1e-300 from the quadratic one.
  # Newton leaves a few eps of t.
  # References: Newton's method on the 1024-bit Rmpfr CDF.
  g <- -100 * sqrt(2)
  p <- c(0.001, 0.3, 0.5, 0.9, 0.999)
  ref <- c(7.0738981023591194e-06, 0.0025217756253144295,
           0.0049006309491189939, 0.016278234326346473,
           0.048823470792455088)
  expect_lt(max(rel_err(qmhn(p, 1, 1, g), ref)), 1e-14)
  expect_lt(max(rel_err(qmhn(1 - p, 1, 1, g, lower.tail = FALSE), ref)),
            1e-14)
  expect_lt(max(rel_err(qmhn(log(p), 1, 1, g, log.p = TRUE), ref)), 1e-14)
  expect_lt(rel_err(qmhn(1e-30, 1, 1, g, lower.tail = FALSE),
                    0.48672839559302161), 1e-14)
  g <- -3000 * sqrt(2)
  expect_lt(max(rel_err(qmhn(c(0.3, 0.9), 1, 1, g),
                        c(8.4069079505944862e-05, 0.00054272438144298491))),
            1e-14)
  g <- -35 * sqrt(2)
  ref <- c(0.013988337287603720, 0.046437864961701644, 0.13905407176194118)
  expect_lt(max(rel_err(qmhn(c(0.5, 0.9, 0.999), 1, 1, g), ref)), 1e-14)
  expect_lt(max(rel_err(qmhn(c(0.5, 0.1, 0.001), 1, 1, g,
                             lower.tail = FALSE), ref)), 1e-14)
  ref <- c(4.2788230137251553, 11.346948831236769)
  expect_lt(max(rel_err(qmhn(c(1e-100, 1e-300), 1, 1, g, lower.tail = FALSE),
                        ref)), 1e-14)
  expect_lt(max(rel_err(qmhn(log(c(1e-100, 1e-300)), 1, 1, g,
                             lower.tail = FALSE, log.p = TRUE), ref)), 1e-14)
})

test_that("the alpha = 1 CDF and quantile are accurate just past z0 = 3", {
  # At beta = 1, gamma = -4.5 the truncation point is z0 = 3.18.  Outside
  # h = d (z0 + d) <= 1/4 the upper tail is taken as log Phibar(z) -
  # log Phibar(z0) in erfcx form, and the quantile through qnorm of
  # Phibar(z0) (1 - p), refined by one Newton step on that erfcx form.
  # References: erfc in Rmpfr at 2048 bits, and
  # Newton's method on that CDF.  Tolerance: log P(X > q) is good to a few
  # units in its last place, and P(X > q) and F carry eps |log P(X > q)| on
  # top of that, 3e-15 at q = 2.
  q <- c(0.1, 0.5, 2)
  log_upper <- c(-0.49765309758535869, -2.6766805139584777,
                 -13.580310005103888)
  upper <- c(0.60795579966052615, 0.068791126754051621,
             1.2651624752445207e-06)
  lower <- c(0.3920442003394739, 0.93120887324594837, 0.99999873483752477)
  for (beta in list(1, c(1, 1, 1))) {
    expect_lt(rel_err(pmhn(q, 1, beta, -4.5, lower.tail = FALSE,
                           log.p = TRUE), log_upper), 2e-14)
    expect_lt(rel_err(pmhn(q, 1, beta, -4.5, lower.tail = FALSE), upper),
              2e-14)
    expect_lt(rel_err(pmhn(q, 1, beta, -4.5), lower), 2e-14)
  }
  expect_lt(rel_err(qmhn(c(0.3, 0.5, 0.9), 1, 1, -4.5),
                    c(0.072050600143254082, 0.13828538820421732,
                      0.43510743547570618)), 2e-14)
})

test_that("the alpha = 1 CDF and quantile invert each other across the tails", {
  # A round trip over thirteen decades of q at tilts on both sides of the
  # truncation point (z0 = -gamma / sqrt(2) from -21 to 21), in both tails.
  # Each direction is accurate to a few 1e-15 relative.  The quantile's
  # relative error is the absolute error of log p divided by
  # |d log p / d log q|, and on this grid both are of the same order (in the
  # upper tail at small q, log p is about -F and accurate relative to F), so
  # the round trip stays near 1e-15; the tolerance leaves a wide margin.
  q <- 10^seq(-12, 0, by = 1)
  for (gamma in c(-30, -3, 0.5, 3, 30)) {
    for (lt in c(TRUE, FALSE)) {
      p <- pmhn(q, 1, 1, gamma, lower.tail = lt, log.p = TRUE)
      expect_lt(rel_err(qmhn(p, 1, 1, gamma, lower.tail = lt, log.p = TRUE),
                        q), 1e-12,
                label = sprintf("gamma=%g lower.tail=%s", gamma, lt))
    }
  }
})

# ============================================================
# 8. One implementation per special case
# ============================================================

test_that("pmhn special cases give the same double for scalar and recycled parameters", {
  # Each special case is evaluated by one per-element function that both the
  # all-scalar path and the recycling loop call, so the two agree exactly under
  # every tail flag.  The points cover F rounding to 1 (log F about -1.6e-28),
  # F below the double range, a strongly tilted truncated normal (z0 = 7e4 and
  # 7e6, taken in erfcx form), and beta q^2 below DBL_MIN.
  pts <- list(c(8, 2, 1, 0), c(1e-40, 10, 1, 0), c(1e-3, 1, 1, 60),
              c(1e-3, 1, 1e-4, -1000), c(1e-7, 1, 1e-6, -1e4),
              c(6, 1, 1, 1), c(1e-170, 0.01, 1, 0))
  for (pt in pts) {
    for (lt in c(TRUE, FALSE)) {
      for (lp in c(TRUE, FALSE)) {
        scalar <- pmhn(pt[1], pt[2], pt[3], pt[4], lower.tail = lt, log.p = lp)
        recycled <- pmhn(pt[1], c(pt[2], pt[2]), pt[3], pt[4],
                         lower.tail = lt, log.p = lp)
        expect_identical(recycled, c(scalar, scalar),
                         info = sprintf("(%g, %g, %g, %g) lt=%s lp=%s",
                                        pt[1], pt[2], pt[3], pt[4], lt, lp))
      }
    }
  }
})

test_that("pmhn special cases with recycled parameters match independent references", {
  # log F at alpha = 2, gamma = 0 is log P(1, 64) = log1p(-exp(-64)), from
  # R's pgamma; F itself rounds to 1.
  expect_lt(rel_err(pmhn(8, c(2, 2), 1, 0, log.p = TRUE),
                    pgamma(64, 1, log.p = TRUE)), 1e-14)
  # log P(5, 1e-80) = 5 log(1e-80) - lgamma(6) + O(1e-80), from Rmpfr; F
  # itself is below the double range.
  expect_lt(rel_err(pmhn(1e-40, c(10, 10), 1, 0, log.p = TRUE),
                    -925.8215289404003), 1e-14)
  # Truncated normal, Rmpfr at 2048 bits.
  expect_lt(rel_err(pmhn(1e-3, c(1, 1), 1, 60, log.p = TRUE),
                    -907.44997056474961), 1e-14)
  # z0 = 70711 and 7.1e6: log Phibar from erfcx by its continued fraction in
  # Rmpfr at 512 bits.  The tolerance allows for the rounding of z0 and of
  # sigma, to which F is proportional here.
  expect_lt(rel_err(pmhn(1e-3, c(1, 1), 1e-4, -1000, lower.tail = FALSE,
                         log.p = TRUE), -1.0000000003), 1e-13)
  expect_lt(rel_err(pmhn(1e-5, c(1, 1), 1e-4, -1000),
                    0.0099501662528219473), 1e-12)
  expect_lt(rel_err(pmhn(c(1e-7, 1e-6), c(1, 1), 1e-6, -1e4),
                    c(0.00099950016662502818, 0.0099501662508321449)), 1e-12)
})

test_that("qmhn special cases give the same double for scalar and recycled parameters", {
  # The recycling loop hands each special-case element to the same
  # per-element quantile the all-scalar path loops over, with the constants
  # its slot holds for the triple.
  p <- c(1e-300, 1e-20, 1e-12, 0.01, 0.3, 0.5, 0.9, 1 - 1e-12)
  pts <- list(c(1, 1, 1), c(1, 1, -3), c(1, 1, -10), c(1, 1, 60),
              c(1, 1e-4, -1000), c(0.01, 1, 0), c(2.5, 3, 0), c(1, 1, 0))
  for (pt in pts) {
    for (lt in c(TRUE, FALSE)) {
      for (lp in c(TRUE, FALSE)) {
        pp <- if (lp) log(p) else p
        scalar <- qmhn(pp, pt[1], pt[2], pt[3], lower.tail = lt, log.p = lp)
        recycled <- qmhn(pp, pt[1], c(pt[2], pt[2]), pt[3],
                         lower.tail = lt, log.p = lp)
        expect_identical(recycled, scalar,
                         info = sprintf("(%g, %g, %g) lower.tail=%s log.p=%s",
                                        pt[1], pt[2], pt[3], lt, lp))
      }
    }
  }
})

# ============================================================
# 9. The closed forms serve alpha = 1 and gamma = 0 exactly
# ============================================================
#
# dmhn, pmhn and qmhn take the truncated normal and the square root of a
# Gamma law only at alpha = 1 and gamma = 0 exactly.  A shape or a tilt
# however close to those values is a different law, and is evaluated as the
# general case.  Near them the law moves smoothly.  With U = sqrt(beta) X
# and Psi[alpha/2, z] = 2 int_0^Inf u^(alpha-1) exp(-u^2 + z u) du,
# differentiating under the integral gives
#
#   d log Psi / d alpha = E log U,   d^2 log Psi / d alpha^2 = Var log U,
#   d log Psi / d z = E U,           d^2 log Psi / d z^2 = Var U,
#
# so that, with the moments taken at alpha = 1 or at z = 0,
#
#   log f(x | 1 + delta) - log f(x | 1)
#       = delta (log(sqrt(beta) x) - E log U) - (delta^2 / 2) Var log U,
#   log f(x | z) - log f(x | 0) = z (sqrt(beta) x - E U) - (z^2 / 2) Var U,
#
# up to terms of order delta^3 or z^3, below 1e-23 at the delta and z of
# 1e-8 used here.  Taken as the special case instead, the difference would
# be 0, which is off by about 7e-6 next to x = 1e-300 under the shape, where
# log(sqrt(beta) x) = -691, and by 5e-8 at sqrt(beta) x = 6 under the tilt.

# E log U and Var log U under alpha = 1, beta = 1 and the tilt z, from
# stats::integrate of the truncated normal kernel exp(-u^2 + z u), split at
# 1 so that the logarithmic singularity at 0 sits at an end point.  The two
# agree with a 200-bit Gauss-Legendre sum to 6e-15 at z = -1, and they enter
# multiplied by delta = 1e-8.
log_u_moments <- function(z) {
  k <- function(u) exp(-u^2 + z * u)
  int <- function(f) {
    stats::integrate(f, 0, 1, rel.tol = 1e-13)$value +
      stats::integrate(f, 1, Inf, rel.tol = 1e-13)$value
  }
  Z <- int(k)
  m <- int(function(u) log(u) * k(u)) / Z
  c(mean = m, var = int(function(u) (log(u) - m)^2 * k(u)) / Z)
}

test_that("dmhn and pmhn near alpha = 1 are the law requested", {
  # Each log density and log probability is formed to a few units of eps of
  # its size, at most 700 here, so to about 2e-13; 1e-11 allows for that.
  d <- 1e-8
  mom <- log_u_moments(-1)
  x <- c(1e-300, 1e-5, 0.5, 3)
  got <- dmhn(x, 1 + d, 1, -1, log = TRUE) - dmhn(x, 1, 1, -1, log = TRUE)
  ref <- d * (log(x) - mom[["mean"]]) - d^2 / 2 * mom[["var"]]
  expect_lt(max(abs(got - ref)), 1e-11)
  # Next to the origin F(q) = 2 beta^(alpha/2) q^alpha / (alpha Psi) to a
  # relative q, so log F moves as the log density does, less log(1 + delta).
  q <- c(1e-300, 1e-200)
  got <- pmhn(q, 1 + d, 1, -1, log.p = TRUE) - pmhn(q, 1, 1, -1, log.p = TRUE)
  ref <- d * (log(q) - mom[["mean"]]) - log1p(d) - d^2 / 2 * mom[["var"]]
  expect_lt(max(abs(got - ref)), 1e-11)
  # Below 1 the same holds, with the moments under z = 2.
  mom <- log_u_moments(2)
  got <- dmhn(x, 1 - d, 1, 2, log = TRUE) - dmhn(x, 1, 1, 2, log = TRUE)
  ref <- -d * (log(x) - mom[["mean"]]) - d^2 / 2 * mom[["var"]]
  expect_lt(max(abs(got - ref)), 1e-11)
  # The recycled-parameter path makes the same choice.
  expect_identical(dmhn(x, c(1 - d, 1 - d), 1, 2, log = TRUE),
                   dmhn(x, 1 - d, 1, 2, log = TRUE))
  expect_identical(pmhn(q, c(1 + d, 1 + d), 1, -1, log.p = TRUE),
                   pmhn(q, 1 + d, 1, -1, log.p = TRUE))
})

test_that("dmhn and pmhn near gamma = 0 are the law requested", {
  # alpha = 2.5, beta = 4 and gamma = 2e-8, so z = 1e-8.  Under gamma = 0,
  # U^2 is Gamma(alpha/2, 1): E U = Gamma((alpha+1)/2) / Gamma(alpha/2),
  # Var U = alpha/2 - (E U)^2, and past t the upper tail has
  # E(U | U > t) = Gamma((alpha+1)/2, t^2) / Gamma(alpha/2, t^2), the
  # upper incomplete Gamma functions, whose logarithms pgamma gives.
  # log P(X > q) moves by z (E(U | U > t) - E U), t = sqrt(beta) q, and
  # next to the origin log F by -z E U.  Tolerances as in the test above.
  a <- 2.5; b <- 4; z <- 1e-8; g <- z * sqrt(b)
  eu <- exp(lgamma((a + 1) / 2) - lgamma(a / 2))
  vu <- a / 2 - eu^2
  x <- c(1e-300, 0.25, 1.5, 3)
  got <- dmhn(x, a, b, g, log = TRUE) - dmhn(x, a, b, 0, log = TRUE)
  ref <- z * (sqrt(b) * x - eu) - z^2 / 2 * vu
  expect_lt(max(abs(got - ref)), 1e-11)
  t <- c(1, 3, 6)
  cond <- exp(lgamma((a + 1) / 2) +
                stats::pgamma(t^2, (a + 1) / 2, lower.tail = FALSE,
                              log.p = TRUE) -
                lgamma(a / 2) -
                stats::pgamma(t^2, a / 2, lower.tail = FALSE, log.p = TRUE))
  got <- pmhn(t / sqrt(b), a, b, g, lower.tail = FALSE, log.p = TRUE) -
    pmhn(t / sqrt(b), a, b, 0, lower.tail = FALSE, log.p = TRUE)
  expect_lt(max(abs(got - z * (cond - eu))), 1e-11)
  got <- pmhn(1e-300, a, b, g, log.p = TRUE) -
    pmhn(1e-300, a, b, 0, log.p = TRUE)
  expect_lt(abs(got - (-z * eu - z^2 / 2 * vu)), 1e-11)
  expect_identical(pmhn(t / sqrt(b), c(a, a, a), b, g, lower.tail = FALSE,
                        log.p = TRUE),
                   pmhn(t / sqrt(b), a, b, g, lower.tail = FALSE,
                        log.p = TRUE))
})

test_that("qmhn near alpha = 1 inverts the law requested", {
  # At p = 1e-300 and 1e-200 the quantile lies next to the origin, where
  # log F(q) = log 2 + alpha log q - log alpha - log Psi[alpha/2, z] at
  # beta = 1, to a relative q.  So log q = (log p - log 2 + log alpha +
  # log Psi) / alpha, with log Psi at alpha = 1 + delta from the expansion
  # above and at alpha = 1 in closed form: Psi[1/2, z] =
  # 2 sqrt(pi) exp(z^2/4) Phi(z / sqrt2).  The root is found to a relative
  # 2^-46 = 1.4e-14, and 1e-12 allows for it; the truncated normal's
  # quantile would be 7e-6 too large at p = 1e-300.
  d <- 1e-8
  z <- -1
  mom <- log_u_moments(z)
  log_psi1 <- log(2 * sqrt(pi)) + z^2 / 4 + pnorm(z / sqrt(2), log.p = TRUE)
  log_psi <- log_psi1 + d * mom[["mean"]] + d^2 / 2 * mom[["var"]]
  p <- c(1e-300, 1e-200)
  ref <- exp((log(p) - log(2) + log1p(d) + log_psi) / (1 + d))
  expect_lt(max(abs(qmhn(p, 1 + d, 1, z) / ref - 1)), 1e-12)
  expect_lt(max(abs(qmhn(p, 1, 1, z) / exp(log(p) - log(2) + log_psi1) - 1)),
            1e-12)
  expect_identical(qmhn(p, c(1 + d, 1 + d), 1, z), qmhn(p, 1 + d, 1, z))
})

# Across 1 +- sqrt(eps) in the shape and +- sqrt(eps) in the tilt, where rmhn
# switches between the special case and the general law, d, p, q and the
# moments are the general law on both sides.  Between 0.99 and 1.01
# sqrt(eps) the law moves by 3e-10 of the shape or the tilt times the
# sensitivity above, |log(sqrt(beta) x) - E log U| or |sqrt(beta) x - E U|,
# below 8 on these grids, so by at most 2.4e-9 relative; 1e-8 allows for
# that and for the rounding of the general forms.  Taking the special case
# on one side would leave a jump of 6e-8 to 1.1e-7 on these grids, and a
# general path that lost mass next to the special values, or a closed form
# that underflowed at a tiny x, a far larger one.

rel_jump <- function(a, b) max(abs(a / b - 1))

test_that("d, p, q and the moments are continuous across the alpha = 1 threshold", {
  e <- sqrt(.Machine$double.eps)
  for (bg in list(c(1, -1e4), c(1e-4, -50), c(1e4, -1e4), c(1, 5))) {
    b <- bg[1]; g <- bg[2]
    # The body and the far upper tail of the law, where pmhn is close to 1.
    x <- qmhn(c(1e-3, 0.5, 0.99, 1 - 1e-12), 1, b, g)
    p <- c(1e-3, 0.5, 0.99)
    for (s in c(-1, 1)) {
      a_in <- 1 + s * 0.99 * e
      a_out <- 1 + s * 1.01 * e
      info <- sprintf("beta=%g gamma=%g alpha=1%+.2ge", b, g, s)
      expect_lt(rel_jump(dmhn(x, a_in, b, g), dmhn(x, a_out, b, g)), 1e-8,
                label = info)
      expect_lt(rel_jump(pmhn(x, a_in, b, g), pmhn(x, a_out, b, g)), 1e-8,
                label = info)
      expect_lt(rel_jump(pmhn(x, a_in, b, g, lower.tail = FALSE),
                         pmhn(x, a_out, b, g, lower.tail = FALSE)), 1e-8,
                label = info)
      expect_lt(rel_jump(qmhn(p, a_in, b, g), qmhn(p, a_out, b, g)), 1e-8,
                label = info)
      expect_lt(rel_jump(mhn_mean(a_in, b, g), mhn_mean(a_out, b, g)), 1e-8,
                label = info)
      expect_lt(rel_jump(mhn_var(a_in, b, g), mhn_var(a_out, b, g)), 1e-8,
                label = info)
    }
  }
})

test_that("d, p, q and the moments are continuous across the gamma = 0 threshold", {
  e <- sqrt(.Machine$double.eps)
  for (ab in list(c(0.01, 1), c(0.01, 1e-6), c(2.5, 1), c(2.5, 1e6))) {
    a <- ab[1]; b <- ab[2]
    x <- qmhn(c(1e-3, 0.5, 0.99, 1 - 1e-12), a, b, 0)
    # At alpha = 0.01 F(x) is still about 0.02 at sqrt(beta) x = 1e-170,
    # where beta x^2 is far below the double range.
    if (a < 1) x <- c(1e-170 / sqrt(b), x)
    p <- c(1e-3, 0.5, 0.99)
    for (s in c(-1, 1)) {
      g_in <- s * 0.99 * e * sqrt(b)
      g_out <- s * 1.01 * e * sqrt(b)
      info <- sprintf("alpha=%g beta=%g gamma=%+.2ge sqrt(beta)", a, b, s)
      expect_lt(rel_jump(dmhn(x, a, b, g_in), dmhn(x, a, b, g_out)), 1e-8,
                label = info)
      expect_lt(rel_jump(pmhn(x, a, b, g_in), pmhn(x, a, b, g_out)), 1e-8,
                label = info)
      expect_lt(rel_jump(pmhn(x, a, b, g_in, lower.tail = FALSE),
                         pmhn(x, a, b, g_out, lower.tail = FALSE)), 1e-8,
                label = info)
      expect_lt(rel_jump(qmhn(p, a, b, g_in), qmhn(p, a, b, g_out)), 1e-8,
                label = info)
      expect_lt(rel_jump(mhn_mean(a, b, g_in), mhn_mean(a, b, g_out)), 1e-8,
                label = info)
      expect_lt(rel_jump(mhn_var(a, b, g_in), mhn_var(a, b, g_out)), 1e-8,
                label = info)
    }
  }
})

# ============================================================
# The truncated normal far from its truncation point
# ============================================================

test_that("the alpha = 1 law under a strong negative tilt is the Exponential it tends to", {
  # With z0 = -gamma / sqrt(2 beta) the truncation point, the density is
  # phi(z0 + d) / (sigma Phibar(z0)), d = x / sigma, sigma = 1 / sqrt(2 beta),
  # and by the asymptotic series of the Mills ratio,
  #   log f(x) = log|gamma| + gamma x - beta x^2
  #              - log(1 - 1/z0^2 + 3/z0^4 - 15/z0^6 + ...),
  # whose omitted terms are below 105 / z0^8, 2.6e-37 at z0 = 7e4.  Taken as
  # dnorm - log Phibar(z0), each of those two logarithms is about -z0^2 / 2
  # and their difference carries eps z0^2 / 2 of rounding: 5.5e-7 at
  # z0 = 7e4.  1e-12 absolute is within a few units of the result's own
  # rounding.
  for (g in c(-1e4, -1e5, -1e20, -1e155)) {
    z0 <- -g / sqrt(2)
    x <- c(0.01, 1, 5) / -g
    ref <- log(-g) + g * x - x^2 - log1p(-1 / z0^2 + 3 / z0^4 - 15 / z0^6)
    info <- sprintf("gamma=%g", g)
    expect_lt(max(abs(dmhn(x, 1, 1, g, log = TRUE) - ref)), 1e-12,
              label = info)
    expect_lt(max(abs(dmhn(x, c(1, 1, 1), 1, g, log = TRUE) - ref)), 1e-12,
              label = info)
  }
  # Past z0 of about 1e8 the law is Exponential(rate -gamma) to within a
  # relative 1 / z0^2 < 1e-16 at every quantile p <= 0.9 and every q below
  # 10 / |gamma|, so the quantiles, the tails and the density are qexp's,
  # pexp's and dexp's to double precision: 1e-14 relative is a few units in
  # the last place, and the quantile 1e-310 at gamma = -1e300 is a subnormal
  # double, resolved to 5e-324.  The quantile's start from the reversion of
  # the series of its lower tail forms z0^2 below z0 = 1e150, and z0 T with
  # T the target over the hazard beyond.
  close <- function(x, ref) {
    all(abs(x - ref) <= 1e-14 * abs(ref) + 1e-323)
  }
  p <- c(1e-10, 0.1, 0.5, 0.9)
  for (g in c(-1e10, -1e100, -1e155, -1e300)) {
    info <- sprintf("gamma=%g", g)
    expect_true(close(qmhn(p, 1, 1, g), qexp(p, -g)), label = info)
    expect_true(close(qmhn(p, 1, 1, g, lower.tail = FALSE),
                      qexp(p, -g, lower.tail = FALSE)), label = info)
    q <- c(1e-3, 1, 10) / -g
    expect_true(close(pmhn(q, 1, 1, g), pexp(q, -g)), label = info)
    expect_true(close(pmhn(q, 1, 1, g, lower.tail = FALSE, log.p = TRUE),
                      pexp(q, -g, lower.tail = FALSE, log.p = TRUE)),
                label = info)
    expect_true(close(dmhn(q, 1, 1, g, log = TRUE), dexp(q, -g, log = TRUE)),
                label = info)
  }
  # Where gamma / (2 beta) overflows while z0 does not, and where
  # z0 = -gamma / sqrt(2 beta) itself overflows (gamma / sqrt(beta) = -Inf),
  # the law is still that Exponential: (1, 1e-300, -1e10) has z0 = 7e159,
  # and (1, 1e-100, -1e300) z0 = Inf.
  for (bg in list(c(1e-300, -1e10), c(1e-100, -1e300))) {
    b <- bg[1]
    g <- bg[2]
    info <- sprintf("beta=%g gamma=%g", b, g)
    q <- c(0.1, 1, 10) / -g
    expect_true(close(pmhn(q, 1, b, g), pexp(q, -g)), label = info)
    expect_true(close(pmhn(q, c(1, 1, 1), b, g), pexp(q, -g)), label = info)
    expect_true(close(pmhn(q, 1, b, g, lower.tail = FALSE, log.p = TRUE),
                      pexp(q, -g, lower.tail = FALSE, log.p = TRUE)),
                label = info)
    expect_true(close(qmhn(p, 1, b, g), qexp(p, -g)), label = info)
    expect_true(close(dmhn(q, 1, b, g, log = TRUE), dexp(q, -g, log = TRUE)),
                label = info)
  }
})

test_that("the alpha = 1 law under a strong positive tilt is measured from its centre", {
  skip_if_not_installed("Rmpfr")
  # Under a strong positive tilt, z0 = -gamma / sqrt(2 beta) far below -3,
  # the law is N(mu, sigma^2) with mu = gamma / (2 beta) and
  # sigma = 1 / sqrt(2 beta), to within the truncated mass Phi(z0), which
  # underflows here.  Unless beta is a power of 2, mu rounds to half a
  # spacing of the doubles at mu, a fraction eps |z0| / 2 of sigma: 1e-6 of
  # it at z0 = 7e9.  The reference takes (x - mu) / sigma in Rmpfr at 256 bits
  # from the doubles passed.  1e-12 is a few units of the rounding of the
  # results.
  b <- 3
  for (z in c(1e4, 1e10, 1e14)) {
    g <- z * sqrt(b)
    info <- sprintf("z=%g", z)
    mu <- Rmpfr::mpfr(g, 256) / (2 * Rmpfr::mpfr(b, 256))
    s <- 1 / sqrt(2 * b)
    x <- Rmpfr::asNumeric(mu) + c(-3, -0.5, 0, 0.5, 3) * s
    t <- Rmpfr::asNumeric((Rmpfr::mpfr(x, 256) - mu) / s)
    expect_lt(max(abs(pmhn(x, 1, b, g) - pnorm(t))), 1e-12, label = info)
    expect_lt(max(abs(pmhn(x, 1, b, g, lower.tail = FALSE) -
                        pnorm(t, lower.tail = FALSE))), 1e-12, label = info)
    expect_lt(max(abs(pmhn(x, c(1, 1, 1, 1, 1), b, g) - pnorm(t))), 1e-12,
              label = info)
    expect_lt(max(abs(dmhn(x, 1, b, g, log = TRUE) -
                        (dnorm(t, log = TRUE) - log(s)))), 1e-12, label = info)
    expect_lt(max(abs(dmhn(x, c(1, 1, 1, 1, 1), b, g, log = TRUE) -
                        (dnorm(t, log = TRUE) - log(s)))), 1e-12, label = info)
  }
})

test_that("the alpha = 1 quantile under a strong positive tilt is the double nearest it", {
  skip_if_not_installed("Rmpfr")
  # As in the previous test the law is N(mu, sigma^2), so its quantile is
  # mu + sigma qnorm(p), rounded here from Rmpfr at 256 bits.  Formed as
  # sigma (qnorm(p) - z0), it would round twice at the size of mu and land a
  # unit or two in its last place away.  At z = 1e4 to 1e14 sigma spans
  # 1e12 to 100 spacings of the doubles at mu, and the quantile is within
  # one spacing of the reference (qnorm itself is accurate to a few units in
  # its last place).  At the three laws below sigma is 1e-54, 1e-79 and
  # 1e-242 of that spacing, so every quantile strictly inside (0, 1) rounds
  # to the double nearest mu, gamma / (2 beta) as R divides it.
  p <- c(1e-3, 0.1, 0.5, 0.9)
  b <- 3
  for (z in c(1e4, 1e10, 1e14)) {
    g <- z * sqrt(b)
    info <- sprintf("z=%g", z)
    mu <- Rmpfr::mpfr(g, 256) / (2 * Rmpfr::mpfr(b, 256))
    s <- 1 / sqrt(Rmpfr::mpfr(2 * b, 256))
    spacing <- 2^(floor(log2(Rmpfr::asNumeric(mu))) - 52)
    expect_lte(max(abs(qmhn(p, 1, b, g) -
                         Rmpfr::asNumeric(mu + s * qnorm(p)))),
               spacing, label = info)
    expect_lte(max(abs(qmhn(p, 1, b, g, lower.tail = FALSE) -
                         Rmpfr::asNumeric(mu + s * qnorm(p, lower.tail = FALSE)))),
               spacing, label = info)
  }
  for (bg in list(c(1e-100, 1e20), c(1e10, 1e100),
                  c(1e100, .Machine$double.xmax))) {
    info <- sprintf("beta=%g gamma=%g", bg[1], bg[2])
    mu <- bg[2] / (2 * bg[1])
    expect_identical(qmhn(p, 1, bg[1], bg[2]), rep(mu, 4), label = info)
    expect_identical(qmhn(p, 1, bg[1], bg[2], lower.tail = FALSE),
                     rep(mu, 4), label = info)
  }
})

test_that("the alpha = 1 quantile from z0 = 3 to 30 is accurate to a few eps", {
  skip_if_not_installed("Rmpfr")
  # From z0 = 3 to 30 the quantile away from the small lower tail is taken
  # through qnorm of Phibar(z0) (1 - p), whose t = z - z0 carries a relative
  # error of about eps z0^2 / h (h = t (z0 + t)), and one Newton step on the
  # erfcx form of the upper tail then leaves a few eps (1 + h) / h of t:
  # about 50 eps at the smallest h on this grid (0.02, at z0 = 3.5 and
  # p = 0.02) and a few eps elsewhere.  The error is measured as
  # (P(q) - p) / (q f(q)), the relative error of q to first order, with P
  # the tail asked for, in Rmpfr at 400 bits.  3e-14 allows for that bound.
  err <- function(q, b, g, p, lower, log_p) {
    B <- Rmpfr::mpfr(b, 400)
    s2b <- sqrt(2 * B)
    z0 <- -Rmpfr::mpfr(g, 400) / s2b
    zq <- Rmpfr::mpfr(q, 400) * s2b + z0
    mass <- Rmpfr::pnorm(z0, lower.tail = FALSE)
    up <- Rmpfr::pnorm(zq, lower.tail = FALSE) / mass
    P <- if (lower) 1 - up else up
    target <- if (log_p) exp(Rmpfr::mpfr(p, 400)) else Rmpfr::mpfr(p, 400)
    dens <- Rmpfr::dnorm(zq) * s2b / mass
    e <- (P - target) / (dens * Rmpfr::mpfr(q, 400))
    Rmpfr::asNumeric(abs(if (lower) e else -e))
  }
  p <- c(0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1 - 1e-6)
  for (z0 in c(3.5, 4, 5, 7, 10, 15, 20, 25, 29, 30)) {
    for (b in c(1e-3, 1, 1e3)) {
      g <- -z0 * sqrt(2 * b)
      info <- sprintf("z0=%g beta=%g", z0, b)
      e_lo <- mapply(err, qmhn(p, 1, b, g), p,
                     MoreArgs = list(b = b, g = g, lower = TRUE, log_p = FALSE))
      e_up <- mapply(err, qmhn(p, 1, b, g, lower.tail = FALSE), p,
                     MoreArgs = list(b = b, g = g, lower = FALSE,
                                     log_p = FALSE))
      e_lg <- mapply(err, qmhn(log(p), 1, b, g, log.p = TRUE), log(p),
                     MoreArgs = list(b = b, g = g, lower = TRUE, log_p = TRUE))
      expect_lt(max(e_lo, e_up, e_lg), 3e-14, label = info)
    }
  }
})

test_that("the alpha = 1 law holds its offset where 2 beta overflows", {
  skip_if_not_installed("Rmpfr")
  # At beta = 1e308, 2 beta is beyond the largest double.  Under a strong
  # positive tilt the offset (x - mu) / sigma is taken from beta x - gamma / 2
  # rounded once, and never from a rounded mu = gamma / (2 beta), which is
  # eps |z0| / 2 standard deviations off: 5e-5 at z0 = -1e12 and 1.1 at
  # -1e16.  At the double nearest mu and its two neighbours the law is
  # N(mu, 1 / (2 beta)) to within the truncated mass Phi(z0), which
  # underflows; the reference takes the offset in Rmpfr at 256 bits.  1e-12
  # is a few units of the rounding of the results, and the quantile of F(x)
  # is x itself, which the quantile returns as the double nearest it.
  b <- 1e308
  for (z0 in c(-1e12, -1e16)) {
    g <- -z0 * sqrt(2) * sqrt(b)
    info <- sprintf("z0=%g", z0)
    mu <- (0.5 * g) / b
    x <- c(mu * (1 - 2^-52), mu, mu * (1 + 2^-52))
    B <- Rmpfr::mpfr(b, 256)
    t <- Rmpfr::asNumeric((Rmpfr::mpfr(x, 256) - Rmpfr::mpfr(g, 256) / (2 * B)) *
                            sqrt(2 * B))
    log_sigma <- -0.5 * log(2) - 0.5 * log(b)
    expect_lt(max(abs(pmhn(x, 1, b, g) - pnorm(t))), 1e-12, label = info)
    expect_lt(max(abs(pmhn(x, 1, b, g, lower.tail = FALSE) -
                        pnorm(t, lower.tail = FALSE))), 1e-12, label = info)
    expect_lt(max(abs(dmhn(x, 1, b, g, log = TRUE) -
                        (dnorm(t, log = TRUE) - log_sigma))), 1e-12,
              label = info)
    expect_lt(max(abs(dmhn(x, c(1, 1, 1), b, g, log = TRUE) -
                        (dnorm(t, log = TRUE) - log_sigma))), 1e-12,
              label = info)
    expect_identical(qmhn(pnorm(t), 1, b, g), x, label = info)
    expect_identical(qmhn(pnorm(t, lower.tail = FALSE), 1, b, g,
                          lower.tail = FALSE), x, label = info)
  }
})

test_that("the alpha = 1 lower tail stays finite next to 0 under a strong positive tilt", {
  skip_if_not_installed("Rmpfr")
  # With z0 = -gamma / sqrt(2 beta) far below 0 the lower tail at d = q / sigma
  # standard deviations above 0 is (Phi(z0 + d) - Phi(z0)) / Phibar(z0),
  # whose logarithm is about -z0^2 / 2 + log(expm1(h)) - log |z0|,
  # h = d |z0|.  z0 + d is resolved only to eps |z0|, which is of the order
  # of d itself at h of order eps z0^2, and there the difference of the two
  # normal tails would be 0; the quantity is finite and ordinary.  The
  # reference forms log Phi(z) = -z^2 / 2 + log(erfcx(-z / sqrt2) / 2) from
  # the asymptotic series of erfcx in Rmpfr, and log F from the exact
  # d = q sqrt(2 beta).  The branch forms log F from terms at most about four
  # times its size, each good to a few eps, so 1e-14 relative allows for it.
  log_erfcx <- function(u) {
    s <- 1
    term <- 1
    w <- 1 / (2 * u^2)
    for (k in 1:30) {
      term <- -term * (2 * k - 1) * w
      s <- s + term
    }
    log(s) - log(u) - 0.5 * log(Rmpfr::Const("pi", Rmpfr::getPrec(u)[1]))
  }
  log_Phi_neg <- function(z) -z^2 / 2 + log_erfcx(-z / sqrt(2)) - log(2)
  for (z in c(1e4, 1e8, 1e16, 1e100)) {
    for (b in c(1e-100, 1, 1e300)) {
      g <- z * sqrt(b)
      if (!is.finite(g)) next
      sigma <- 1 / sqrt(2 * b)
      z0 <- -g / sqrt(2 * b)
      pr <- 2 * ceiling(log2(z)) + 200
      B <- Rmpfr::mpfr(b, pr)
      Z0 <- -Rmpfr::mpfr(g, pr) / sqrt(2 * B)
      for (q in c(c(0.3, 3, 30) / abs(z0), c(0.1, 0.49) * abs(z0)) * sigma) {
        info <- sprintf("z=%g beta=%g q=%g", z, b, q)
        ZQ <- Rmpfr::mpfr(q, pr) * sqrt(2 * B) + Z0
        r <- log_Phi_neg(ZQ) - log_Phi_neg(Z0)
        ref <- Rmpfr::asNumeric(log_Phi_neg(Z0) + r + log1p(-exp(-r)))
        expect_lt(abs(pmhn(q, 1, b, g, log.p = TRUE) / ref - 1), 1e-14,
                  label = info)
        expect_lt(abs(pmhn(q, c(1, 1), b, g, log.p = TRUE)[1] / ref - 1),
                  1e-14, label = info)
      }
    }
  }
  # At z = 1e8, beta = 1e-100 and h = 0.3, and at the far end of the range.
  expect_true(is.finite(pmhn(3e41, 1, 1e-100, 1e-42, log.p = TRUE)))
  # There log F is -z0^2 / 2 = -2.5e299 to a relative 1e-296, the rest of it
  # being of the size of log z0.
  expect_lt(abs(pmhn(1e-300, 1, 1e300, 1e300, log.p = TRUE) / -2.5e299 - 1),
            1e-14)
})

test_that("the alpha = 1 quantile far into the upper tail of a far truncation converges", {
  # For z0 > 30 the upper-tail quantile solves
  #   -t (z0 + t/2) + log(erfcx((z0 + t)/sqrt2) / erfcx(z0/sqrt2)) = log p
  # by Newton's method, from the root of its quadratic part wherever -log p
  # exceeds 1e-3 z0^2, which lies just to the right of the root however far
  # out log p is (below that, from the root of its linear part).  The start
  # is formed from L / 2 and the equation divided by 4, so that nothing
  # overflows up to L = -log p = DBL_MAX, 1.8e308.  The
  # distribution function, taken from the same erfcx form to a few eps of
  # log p, returns log p at the quantile; a relative error e in q moves it by
  # about 2 e this far out, so 1e-14 bounds e by a few eps.  At log p = -L
  # with L far above z0^2 the quantile is sigma sqrt(2 L) to a relative
  # z0 / sqrt(2 L).
  for (z0 in c(35, 1e3, 1e6)) {
    g <- -z0 * sqrt(2)
    for (L in c(1e10, 1e40, 1e100, 1e300, 1.7e308, .Machine$double.xmax)) {
      info <- sprintf("z0=%g log p=%g", z0, -L)
      q <- qmhn(-L, 1, 1, g, lower.tail = FALSE, log.p = TRUE)
      expect_lt(abs(pmhn(q, 1, 1, g, lower.tail = FALSE, log.p = TRUE) / -L -
                      1), 1e-14, label = info)
      if (L >= 1e40) {
        expect_lt(abs(q / sqrt(L) - 1), 1e-8, label = info)
      }
    }
  }
})
