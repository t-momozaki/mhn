# Tests for the RTDR (Relaxed Transformed Density Rejection) sampling
# kernel, region by region.  Exercises the internal `.rmhn_rtdr_cpp`
# directly so that regions can be validated independently of the
# dispatcher.
#
# (a)-(d) are the vignette's parameter regions; A, BC and D the RTDR setups
# ((b) and (c) share BC).
#
# Targets the per-setup KS tests and the continuity check at the boundary
# of regions (c) and (d), that is of setups BC and D, exercising Gao & Wang
# (2025) Theorems 3.1, 3.2, 4.4.  Acceptance > 1/e is asserted in every
# setup.  Gao & Wang prove that bound for contact points in the band of
# their Eq. (7) and, in region D, for the full secant chain of their
# Theorem 4.4, which the package truncates (see the early-termination
# tests below); in region D, the large-gamma corner included, the check is
# therefore empirical.
#
# KS p-value threshold > 0.001 keeps the false-rejection rate < 0.1%
# per (region, parameter) combination.  The references, mhn_ref() and
# mhn_cdf(), are in helper-mhn-cdf.R and are formed from the kernel of the
# law alone.

run_rtdr_region <- function(alpha, beta, gamma, n = 10000, seed = 1L,
                            ks_threshold = 0.001) {
  set.seed(seed)
  x <- mhn:::.rmhn_rtdr_cpp(n, alpha, beta, gamma)
  ref <- mhn_ref(alpha, beta, gamma)
  list(
    samples       = x,
    n             = n,
    region        = attr(x, "rtdr_region"),
    retries       = attr(x, "rtdr_retries"),
    acceptance    = n / (n + attr(x, "rtdr_retries")),
    ref           = ref,
    ks            = ks.test(x, ref$cdf),
    sample_mean   = mean(x),
    sample_var    = var(x),
    theory_mean   = ref$mean,
    theory_var    = ref$var
  )
}

# Region A: alpha > 1 (alpha = 1 is sampled in region BC), log-concave on f(x).
# Expected: region = 0; acceptance >= 1/e when the contact points lie in the
# band of Gao & Wang's Eq. (7).
test_that("region A: KS, mean/var, acceptance for representative points", {
  skip_on_cran()
  for (params in list(c(1.5, 1, -2), c(5, 1, 0.5), c(10, 1, 10),
                      c(1.5, 4, -4), c(5, 0.25, 0.25))) {
    r <- run_rtdr_region(params[1], params[2], params[3])
    expect_equal(r$region, 0L,
                 info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
    expect_gt(r$ks$p.value, 0.001)
    expect_gt(r$acceptance, 1 / exp(1))
    # Sample mean and variance within five of their standard errors.
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se)
    expect_var_near(r$sample_var, r$ref, r$n,
                    info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
  }
})

# Region BC: alpha in [1/2, 1], or alpha < 1/2 at or below the D threshold,
# on g(y), serving regions (b) and (c).  For gamma <= 0 the density is
# log-concave and a T_0 (log-tangent) hat is used; for gamma > 0 it is only
# T_{-1/2}-concave (Gao & Wang 2025, Lemma 3.3 and Theorem 3.2) and a
# T_{-1/2} (inverse-square) hat is used.  gamma = 2 is included because a
# log-tangent hat would not dominate g there, and would bias the sample
# upward.
test_that("region BC: KS, mean/var, acceptance for representative points", {
  skip_on_cran()
  cases <- list(
    c(0.7, 1, -2), c(0.7, 1, 0),                     # region (b), gamma <= 0: T_0
    c(0.7, 1,  2), c(0.7, 1, 5),                     # region (b), gamma > 0: T_{-1/2}
    c(0.3, 1, -5), c(0.3, 1, 0), c(0.3, 1, 0.5),     # region (c), all below the BC/D threshold
    c(0.7, 4, -4), c(0.7, 0.25, 1), c(0.3, 9, 1.5)   # the same regions at beta != 1
  )
  for (params in cases) {
    r <- run_rtdr_region(params[1], params[2], params[3])
    expect_equal(r$region, 1L,
                 info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
    expect_gt(r$ks$p.value, 0.001)
    expect_gt(r$acceptance, 1 / exp(1))
    se <- sqrt(r$theory_var / r$n)
    expect_lt(abs(r$sample_mean - r$theory_mean), 5 * se)
    expect_var_near(r$sample_var, r$ref, r$n,
                    info = sprintf("alpha=%g, gamma=%g", params[1], params[3]))
  }
})

# Region D: alpha < 1/2 and gamma > threshold.
# Gao & Wang (2025, Eqs. (14)-(15)) use envelope (14) when
# log(g(m_g)/g(y*)) <= 2.49 and envelope (15) otherwise; their Theorem 4.4
# bounds acceptance below by 1/e for both, with the full secant chain.  The
# package truncates that chain, so here the bound is checked empirically,
# and the (0.1, 50) case, which uses (15), is held to it as well.
test_that("region D moderate (envelope (14)): KS + acceptance >= 1/e", {
  skip_on_cran()
  r <- run_rtdr_region(0.3, 1, 5)
  expect_equal(r$region, 2L)
  expect_gt(r$ks$p.value, 0.001)
  expect_gt(r$acceptance, 1 / exp(1))
})

test_that("region D extreme (envelope (15)): KS + acceptance >= 1/e", {
  skip_on_cran()
  # Gao & Wang (2025, Theorem 4.4) bound the acceptance below by 1/e over the
  # whole of region D, the large-gamma extreme included, for the full secant
  # chain; the package truncates the chain, so the check is empirical.
  r <- run_rtdr_region(0.1, 1, 50, n = 5000)
  expect_equal(r$region, 2L)
  expect_gt(r$ks$p.value, 0.001)
  expect_gt(r$acceptance, 1 / exp(1))
  # Confirm envelope (15) is selected, not envelope (14).
  env <- mhn:::.dump_rtdr_envelope_cpp(0.1, 1, 50)
  expect_true(isTRUE(env$has_left_tangent_d))
})

# Envelope switching criterion: log(g(m_g)/g(y*)) > 2.49 -> envelope (15).
# Confirms envelope (14) vs (15) selection across the documented threshold.
test_that("region D envelope (14) vs (15) selection follows its threshold", {
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

# Envelope (15) internal identity: by construction y_K at k = K (the
# unterminated index) equals t_dot_l = log(gamma/2 - exp(t_l)), the dual
# point of the contact point t_l (Gao & Wang 2025, Lemma 4.5, printed
# \dot t_l; its closed form is in the paragraph before their Eq. (14) and
# in the definition of \hat y* after Eq. (15)).
# Early termination on the leftmost secant breakpoint (Gao & Wang 2025,
# Appendix B code; Section 4, paragraph before Lemma 4.5) means the LAST
# stored breakpoint y_b[K_eff] generally does NOT equal t_dot_l; it stops
# earlier.  Here we verify the identity by feeding the original K back
# into the breakpoint formula and comparing the result to t_dot_l
# (independent of K_eff).
test_that("envelope (15): y_K (formula at original K) equals t_dot_l", {
  for (p in list(c(0.3, 10), c(0.1, 50))) {
    env <- mhn:::.dump_rtdr_envelope_cpp(p[1], 1, p[2])
    expect_true(isTRUE(env$has_left_tangent_d))
    t_dot_l_expected <- log(p[2] / 2 - exp(env$t_l))
    K_orig <- as.integer(ceiling(env$rho))
    drho <- env$rho / K_orig
    rk <- 4 * K_orig * drho
    inside <- p[2]^2 - rk
    expect_gt(inside, 0)
    y_K_formula <- log(2 * K_orig * drho / (sqrt(inside) + p[2]))
    expect_lt(abs(y_K_formula - t_dot_l_expected), 1e-9)
  }
})

# Verify that the early-termination optimisation (Gao & Wang 2025,
# Appendix B code; Section 4, paragraph before Lemma 4.5) collapses K_eff
# to a small value (typically 1) for envelope (15) in the extreme corner of
# region D, while leaving acceptance well above 1/e.
test_that("early termination shrinks K_eff in region D for extreme gamma", {
  for (p in list(c(0.3, 10), c(0.1, 50), c(0.3, 10000), c(0.01, 10000))) {
    env <- mhn:::.dump_rtdr_envelope_cpp(p[1], 1, p[2])
    expect_true(isTRUE(env$has_left_tangent_d))
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

# The BC/D split lies at gamma = 2 (1 - sqrt(1 - 2 alpha)), about 0.7350 for
# alpha = 0.3.  Close to it on either side the region is the one expected,
# and 2e5 draws keep the reference law (chi-square on 100 bins) and its mean
# (five standard errors).
test_that("region BC/D boundary: each side keeps its law", {
  skip_on_cran()
  alpha <- 0.3
  threshold <- 2 * (1 - sqrt(1 - 2 * alpha))
  for (g in threshold + c(-0.02, 0.02)) {
    ref <- mhn_ref(alpha, 1, g)
    set.seed(99L)
    x <- mhn:::.rmhn_rtdr_cpp(2e5, alpha, 1, g)
    info <- sprintf("gamma = %.4f", g)
    expect_identical(attr(x, "rtdr_region"), if (g < threshold) 1L else 2L,
                     info = info)
    expect_gt(chisq_ref(x, ref), 1e-4, label = info)
    expect_lt(abs(mean(x) - ref$mean) / sqrt(ref$var / 2e5), 5, label = info)
  }
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

# ---------------------------------------------------------------------
# Large tilts and shapes: the grid of doubles at the mode.
# ---------------------------------------------------------------------
#
# Once the law's standard deviation is a modest multiple of the spacing of the
# doubles at its mode, its draws can only take the values of those doubles.
# With log f(m + d) - log f(m) = (alpha - 1) [log1p(d/m) - d/m] - beta d^2,
# the law of X near its mode m is a normal of variance
# 1/(2 beta + (alpha - 1)/m^2) -- 1/(2 beta) on the tilt axis, about
# 1/(4 beta) on the shape axis -- up to a log-density term
# (alpha - 1) d^3 / (3 m^3), which at the points below is under 1e-37 on the
# tilt axis and under 1e-14 over +-5 standard deviations on the shape axis.
# The tests compare draws with that normal rounded to the nearest double: a
# draw exact to double precision is distributed exactly so.

# Chi-square test of draws x against N(mode, sd^2) rounded to the nearest
# double.  The exact mode is mode_hi + mode_lo, with mode_hi a double and
# mode_lo below one spacing; every draw must lie in the binade of mode_hi, so
# that x - mode_hi is exact and no digit of the sub-spacing offset is lost.
# Neighbouring cells are pooled until each expects at least five draws.
rounded_normal_chisq_p <- function(x, mode_hi, mode_lo, sd) {
  spacing <- 2^(floor(log2(mode_hi)) - 52)
  i <- (x - mode_hi) / spacing
  stopifnot(all(i == round(i)), all(floor(log2(x)) == floor(log2(mode_hi))))
  cells <- seq(min(i), max(i))
  z <- cells * spacing - mode_lo
  upper <- c(pnorm((z[-length(z)] + spacing / 2) / sd), 1)
  lower <- c(0, pnorm((z[-1] - spacing / 2) / sd))
  expected <- length(x) * (upper - lower)
  observed <- tabulate(i - min(i) + 1, nbins = length(cells))
  o <- e <- numeric(0)
  acc_o <- acc_e <- 0
  for (j in seq_along(cells)) {
    acc_o <- acc_o + observed[j]
    acc_e <- acc_e + expected[j]
    if (acc_e >= 5) {
      o <- c(o, acc_o); e <- c(e, acc_e); acc_o <- acc_e <- 0
    }
  }
  o[length(o)] <- o[length(o)] + acc_o
  e[length(e)] <- e[length(e)] + acc_e
  pchisq(sum((o - e)^2 / e), df = length(e) - 1, lower.tail = FALSE)
}

# The mode at beta = 1 and a large tilt, gamma/2 + w with
# w = 2 (alpha - 1) / (S + gamma) and S = sqrt(gamma^2 + 8 (alpha - 1)): an
# even gamma below 2^54 makes gamma/2 a double, and w, a sum of positive terms,
# is formed to full relative precision.
tilt_mode <- function(alpha, gamma) {
  list(hi = gamma / 2,
       lo = 2 * (alpha - 1) / (sqrt(gamma^2 + 8 * (alpha - 1)) + gamma),
       sd = 1 / sqrt(2 + (alpha - 1) / (gamma / 2)^2))
}

test_that("for alpha < 1 large-tilt draws fall on the grid of x, not of log x", {
  # For alpha < 1 the envelope is built on the log axis y = log(sqrt(beta) x).
  # A draw formed as exp(y) from an absolute log-axis double can take only the
  # values exp(k ulp(y)), a grid about x ulp(log x) wide: 0.018 at a tilt of
  # 1e13 -- eighteen spacings of the doubles at the mode and 2.5% of the
  # standard deviation -- and 0.18 at 1e14.  Formed from the offset from the
  # mode, a draw can reach every double near the mode.
  #
  # The standardised offsets are tested against pnorm; x - m is exact because
  # x and m lie within a factor of two of each other.  At 1e14 the grid of x
  # is itself 1.1% of a standard deviation, which moves the Kolmogorov-Smirnov
  # statistic of an exact sampler by about 0.002, against a 0.1% critical value
  # of 0.014 at n = 2e4.  The sd tolerance is six standard errors.
  n <- 2e4
  for (alpha in c(0.3, 0.7)) {                   # region D, region BC
    for (tilt in c(1e13, 1e14)) {
      info <- sprintf("alpha = %g, gamma = %g", alpha, tilt)
      set.seed(1)
      x <- rmhn(n, alpha, 1, tilt)
      m <- mhn_mode(alpha, 1, tilt)
      spacing <- 2^(floor(log2(m)) - 52)        # ulp of every draw
      expect_identical(min(diff(sort(unique(x)))), spacing, info = info)
      z <- (x - m) * sqrt(2)
      expect_gt(suppressWarnings(ks.test(z, "pnorm"))$p.value, 0.001,
                label = info)
      expect_lt(abs(mean(z)) * sqrt(n), 5, label = info)
      expect_lt(abs(sd(z) - 1), 0.03, label = info)
    }
  }
  # The same law at beta != 1.
  set.seed(2)
  x <- rmhn(n, 0.3, 1e-8, 1e9)
  z <- (x - mhn_mode(0.3, 1e-8, 1e9)) * sqrt(2e-8)
  expect_gt(suppressWarnings(ks.test(z, "pnorm"))$p.value, 0.001)
})

test_that("for alpha < 1 draws are the law rounded to double up to the grid limit", {
  # An envelope placed on the absolute log axis cannot resolve the mode once
  # one ulp of log x exceeds the law's width there, about 1.4 / gamma (gamma
  # ~ 3e14): its draws then take three or four values, or the construction
  # fails.  Built on the offset from the mode, the envelope has no such
  # floor, and a draw is the continuous draw rounded once to the nearest
  # double.  The only limit left is the grid of x: at 5e14 its spacing is
  # 1/32 (4.4% of a standard deviation), at 4e15 it is 1/4 and at 8e15 it is
  # 1/2 (0.71 of one).
  n <- 2e4
  for (alpha in c(0.3, 0.7)) {                   # region D, region BC
    for (tilt in c(5e14, 4e15, 8e15)) {
      set.seed(3)
      x <- rmhn(n, alpha, 1, tilt)
      m <- tilt_mode(alpha, tilt)
      expect_gt(rounded_normal_chisq_p(x, m$hi, m$lo, m$sd), 0.001,
                label = sprintf("alpha = %g, gamma = %g", alpha, tilt))
    }
  }
})

# Reference modes for the two tests below: the root
# (gamma + sqrt(gamma^2 + 8 beta (alpha - 1))) / (4 beta), or
# 2 (alpha - 1) / (sqrt(gamma^2 + 8 beta (alpha - 1)) - gamma) for gamma < 0,
# evaluated with Rmpfr at 300 bits from the double inputs and split into the
# nearest double `hi` and the remainder `lo`, which is given to six digits.
# hi is written as an integer plus a binary fraction, both exact, because R
# does not read every 17-digit decimal to the nearest double (it reads
# 912866666666666.75 as 912866666666666.875).  The standard deviation
# 1/sqrt(2 beta + (alpha - 1)/m^2) is formed from hi.
exact_modes <- list(
  list(alpha = 0.3, beta = 0.3, gamma = 5.4772e14,
       hi = 912866666666666 + 0.75, lo = -0.0495505),
  list(alpha = 0.7, beta = 0.3, gamma = 5.4772e14,
       hi = 912866666666666 + 0.75, lo = -0.0495505),
  list(alpha = 2, beta = 0.3, gamma = 5.4772e14,
       hi = 912866666666666 + 0.75, lo = -0.0495505),
  list(alpha = 0.3, beta = 2, gamma = 1.4142e15,
       hi = 353550000000000, lo = -4.94979e-16),
  list(alpha = 0.7, beta = 2, gamma = 1.4142e15,
       hi = 353550000000000, lo = -2.12134e-16),
  list(alpha = 2, beta = 2, gamma = 1.4142e15,
       hi = 353550000000000, lo = 7.07114e-16),
  list(alpha = 3e13 + 1, beta = 1, gamma = 1e15,
       hi = 5e14, lo = 0.03),
  list(alpha = 1e29, beta = 1, gamma = -1,
       hi = 223606797749978 + 0.71875, lo = -0.00868711),
  list(alpha = 5e28, beta = 1, gamma = -1,
       hi = 158113883008418 + 0.71875, lo = -0.00892274),
  list(alpha = 1e30, beta = 1, gamma = 1,
       hi = 707106781186547 + 0.75, lo = 0.0314311),
  list(alpha = 8e30, beta = 1, gamma = -1,
       hi = 1999999999999999 + 0.75, lo = 0.0198846),
  list(alpha = 0.7, beta = 0.3, gamma = 5.2034e15,
       hi = 8672333333333334, lo = -0.345726)
)

# Draws of rmhn(n, alpha, beta, gamma, method) against the law rounded to
# double around the reference mode: the chi-square above, and the sample mean
# within five standard errors of the mode (the mean of X lies within
# 1/(8 m) standard deviations of the mode here, far below one standard error).
check_rounded_law <- function(ref, method, n = 2e4, seed = 6) {
  info <- sprintf("alpha = %g, beta = %g, gamma = %g, method = %s",
                  ref$alpha, ref$beta, ref$gamma, method)
  sd <- 1 / sqrt(2 * ref$beta + (ref$alpha - 1) / ref$hi^2)
  set.seed(seed)
  x <- rmhn(n, ref$alpha, ref$beta, ref$gamma, method = method)
  z <- ((x - ref$hi) - ref$lo) / sd
  expect_lt(abs(mean(z)) * sqrt(n), 5, label = info)
  expect_gt(rounded_normal_chisq_p(x, ref$hi, ref$lo, sd), 0.001,
            label = info)
}

test_that("draws at beta != 1 are placed at the mode in the caller's coordinate", {
  # The envelope is built for sqrt(beta) X, whose mode is about gamma_norm/2
  # with gamma_norm = gamma / sqrt(beta).  The doubles gamma_norm and
  # sqrt(beta) each carry a rounding of up to half a unit in the last place,
  # so a draw returned as (sqrt(beta) X) / sqrt(beta) sits up to about
  # 2e-16 gamma_norm standard deviations from where it belongs -- 0.2 at
  # gamma_norm = 1e15, where the doubles at the mode are about a tenth of a
  # standard deviation apart.  The draws are therefore placed at the mode
  # computed from the caller's own beta and gamma.  The first six points have
  # gamma_norm = 1e15 at beta = 0.3 and 2, for regions D, BC and A.  The last
  # has gamma_norm = 9.5e15 at beta = 0.3, past the grid limit of beta = 1 but
  # not of beta = 0.3 (see the next test but one): the doubles at its mode are
  # a unit apart, 0.77 of its standard deviation.
  for (ref in exact_modes[c(1:6, 12)]) {
    check_rounded_law(ref, if (ref$alpha > 1) "rtdr" else "auto")
  }
})

test_that("region A draws are the law rounded to double around the exact mode", {
  # For alpha > 1 this sampler draws the offset of X from the mode and adds it
  # to the mode carried to about twice double precision, so a draw is the
  # exact draw rounded once, whether or not the mode is itself a double and
  # however few doubles the law spans.  At the first point the mode is 0.03,
  # 0.042 standard deviations, above the double 5e14; the shape-axis points
  # have standard deviation 1/2 against spacings of 1/32 (5e28), 1/32 (1e29),
  # 1/8 (1e30) and 1/4 (8e30).  On the tilt axis at alpha = 2 the doubles at
  # the mode are 1/16 apart at a tilt of 1e15 and 1/2, 0.71 standard
  # deviations, at 6e15.
  for (ref in exact_modes[7:11]) check_rounded_law(ref, "rtdr")
  for (tilt in c(1e15, 6e15)) {
    m <- tilt_mode(2, tilt)
    set.seed(5)
    x <- rmhn(2e4, 2, 1, tilt, method = "rtdr")
    expect_gt(rounded_normal_chisq_p(x, m$hi, m$lo, m$sd), 0.001,
              label = sprintf("alpha = 2, gamma = %g", tilt))
  }
})

test_that("past the grid limit rmhn says why it cannot sample", {
  # Every region stops once the law's standard deviation is no larger than the
  # spacing of the doubles at its mode, in the caller's coordinate.  At
  # beta = 1 that is from a tilt of 2^53 = 9.0e15, where the mode exceeds 2^52
  # and the doubles there are a whole unit apart against a standard deviation
  # of 1/sqrt(2), and on the shape axis at a small tilt at the latest from
  # 2^105 = 4.1e31, where they are a unit apart against 1/2.  (A strong
  # negative tilt widens the law relative to its mode, and moves the stop to
  # between 2^104 and 2^106.)
  for (alpha in c(0.3, 0.7)) {
    for (tilt in c(1e16, 1e20)) {
      expect_error(rmhn(5, alpha, 1, tilt),
                   "narrower than the spacing of the double precision numbers")
      expect_error(rmhn(5, alpha, 1e-6, tilt * 1e-3),
                   "narrower than the spacing of the double precision numbers")
    }
  }
  for (args in list(c(2, 1e16), c(2, 1e20), c(5e31, -1), c(1e33, 1))) {
    expect_error(rmhn(5, args[1], 1, args[2], method = "rtdr"),
                 "narrower than the spacing of the double precision numbers",
                 info = sprintf("alpha = %g, gamma = %g", args[1], args[2]))
  }
  # The limit is decided in the caller's coordinate, not in the beta = 1
  # coordinate of the envelope.  At beta = 1 a tilt of 9.5e15 puts the mode at
  # 4.75e15, where the doubles are a unit apart, and is refused; at
  # beta = 0.3, gamma = 5.2034e15 has the same gamma / sqrt(beta), but the
  # mode, 8.67e15, still lies below 2^53, where the doubles are a unit apart
  # against a standard deviation of 1/sqrt(0.6) = 1.29, and is sampled.
  expect_error(rmhn(5, 0.7, 1, 9.5e15),
               "narrower than the spacing of the double precision numbers")
  expect_true(all(is.finite(rmhn(5, 0.7, 0.3, 5.2034e15))))
  # Under a strong negative tilt the law is close to Gamma(alpha, -gamma), of
  # relative width 1/sqrt(alpha - 1).  With gamma = -(alpha - 1) / 1.999 the
  # mode lies just below 2, where the doubles are 2^-52 apart, and the
  # standard deviation, 1.999 / sqrt(alpha - 1), is 1.08 spacings at
  # alpha = 7e31, past 2^105, and 0.995 at 1.01 * 2^106.  The default method
  # draws there from RTDR, which must sample the first, every draw within
  # 8 standard deviations of the mode, and refuse the second.
  a <- 7e31
  g <- -(a - 1) / 1.999
  set.seed(3)
  x <- rmhn(200, a, 1, g)
  expect_lt(max(abs(x - mhn_mode(a, 1, g))), 8 * 1.999 / sqrt(a - 1))
  a <- 1.01 * 2^106
  expect_error(rmhn(5, a, 1, -(a - 1) / 1.999),
               "narrower than the spacing of the double precision numbers")
  # The normal limit it offers has variance 1 / (-(log f)''(m)), which is
  # 1/(2 beta + (alpha - 1)/m^2): 1/(2 beta) only on the tilt axis, about
  # 1/(4 beta) on the shape axis.
  msg <- tryCatch(rmhn(5, 1e33, 1, 1, method = "rtdr"),
                  error = conditionMessage)
  expect_match(msg, "1/(2*beta + (alpha - 1)/m^2)", fixed = TRUE)
  expect_false(grepl("defect", msg, fixed = TRUE))
})

test_that("a raw peak that overflows is refused, not reported as a defect", {
  # Past |gamma| / sqrt(beta) of about 1.3e154 the square of the tilt
  # overflows.  With a positive tilt the raw log-density at the mode is then
  # Inf - Inf = NaN, and the law is far narrower than one spacing of the
  # doubles at its mode.  With a negative tilt the law is Gamma(alpha, rate
  # |gamma|) to double precision, of relative width 1/sqrt(alpha), which is
  # narrower than the spacing of the doubles wherever its mode lies once alpha
  # passes 2^106 = 8.1e31, as in the last three cases.  That is what the
  # message must say.
  for (args in list(c(0.3, 1, 1e160), c(0.7, 1, 1e200), c(2, 1, 1e300),
                    c(0.3, 1e-310, 1), c(1e33, 1, -1e160),
                    c(1e40, 1, -1e160), c(1e300, 1, -1e300))) {
    msg <- tryCatch({
      rmhn(3, args[1], args[2], args[3], method = "rtdr")
      "no error"
    }, error = conditionMessage)
    info <- paste(args, collapse = ", ")
    expect_match(msg, "narrower than the spacing", info = info)
    expect_false(grepl("defect", msg, fixed = TRUE), info = info)
  }
})

# ---------------------------------------------------------------------
# Strong negative tilts: the gamma limit.
# ---------------------------------------------------------------------
#
# With gamma / sqrt(beta) far below zero the draws lie near alpha / |gamma|.
# At every point below where the draws are compared with pgamma, the factor
# exp(-beta x^2) of the density differs from 1 by less than 1e-35 over the
# range they occupy, so the law is Gamma(alpha, rate |gamma|) to double
# precision and |gamma| X is compared with pgamma(shape = alpha).
gamma_limit_ks_p <- function(x, alpha, gamma) {
  ks.test(abs(gamma) * x, "pgamma", shape = alpha)$p.value
}

test_that("at a strong negative tilt region A draws follow the gamma limit", {
  # rmhn(2e4, 2, 1e-56, -1) is Gamma(2, 1), and rmhn's default method sends
  # it to RTDR, as it does for gamma < 0 with 1 < alpha < 1e3 at any n up to
  # a tilt -gamma / sqrt(beta) of 1e150, here 1e28.  The envelope is built
  # for sqrt(beta) X, whose mode is 1e-28 here.
  set.seed(7)
  x <- rmhn(2e4, 2, 1e-56, -1)
  expect_false(anyNA(x))
  expect_gt(gamma_limit_ks_p(x, 2, -1), 0.001)
  # With method = "rtdr", modes of 1e-28 to 1e-18: alpha below e, where the
  # envelope is a plateau from x = 0 and a right tail, alpha above it, where a
  # left tangent is added, and a large shape.
  for (args in list(c(1.2, -1e27), c(2, -1e28), c(3, -1e28), c(5, -1e28),
                    c(1e12, -1e30))) {
    info <- sprintf("alpha = %g, gamma = %g", args[1], args[2])
    set.seed(8)
    x <- rmhn(2e4, args[1], 1, args[2], method = "rtdr")
    expect_false(anyNA(x), info = info)
    expect_gt(gamma_limit_ks_p(x, args[1], args[2]), 0.001, label = info)
  }
})

test_that("at a strong negative tilt region A keeps its acceptance bound", {
  # Gao & Wang (2025, Eq. (7)) bound the rejection constant by e once the
  # log-density at each contact point lies between 0.46 and 2.49 below its
  # value at the mode, so the contact search has to reach that band however
  # narrow the law is.  The law's standard deviation, about
  # sqrt(alpha) / |gamma|, is 3e-14 and 1e-15 at these two points.  The
  # expected acceptance is above 0.8, so a sample of 1e4 draws cannot fall to
  # 1/e by chance.
  for (args in list(c(1e3, -1e15), c(1e6, -1e18))) {
    info <- sprintf("alpha = %g, gamma = %g", args[1], args[2])
    set.seed(9)
    x <- mhn:::.rmhn_rtdr_cpp(1e4, args[1], 1, args[2])
    expect_false(anyNA(x), info = info)
    expect_gt(1e4 / (1e4 + attr(x, "rtdr_retries")), 1 / exp(1),
              label = info)
  }
})

test_that("a negative tilt whose square overflows is sampled from its gamma limit", {
  # Past |gamma| / sqrt(beta) of about 1.3e154 the square of the tilt
  # overflows, but the law, Gamma(alpha, rate |gamma|) to double precision,
  # has a relative width of about 1/sqrt(alpha), which stays above the
  # spacing of the doubles for alpha below 2^104, and the sampler finds the
  # mode without squaring the tilt.  Region BC (alpha <= 1), region A on both
  # sides of alpha = e, and two large shapes, the second far enough out that
  # region A measures its envelope from the mode, where (alpha - 1)/m^2 is
  # about gamma^2/alpha = 1e388 and the derivative has to be formed without
  # it.
  for (args in list(c(0.5, -1e160), c(2, -1e160), c(5, -1e200),
                    c(1e10, -1e200), c(1e12, -1e200))) {
    info <- sprintf("alpha = %g, gamma = %g", args[1], args[2])
    set.seed(10)
    x <- rmhn(2e4, args[1], 1, args[2], method = "rtdr")
    expect_false(anyNA(x), info = info)
    expect_gt(gamma_limit_ks_p(x, args[1], args[2]), 0.001, label = info)
  }
  # At alpha = 1e31 the standard deviation, m / sqrt(alpha - 1) = 3.2e-145,
  # is two spacings of the doubles at the mode m = 1e-129, so the draws are
  # compared with the normal limit rounded to double, as above.  Over +-5
  # standard deviations its log-density differs from the gamma law's by the
  # cubic term (alpha - 1) d^3 / (3 m^3), under 2e-14, and the mean of the
  # gamma law lies 1/sqrt(alpha) standard deviations from its mode.  The mode
  # is from Rmpfr at 400 bits, in hexadecimal, as is the tilt: the decimal
  # -1e160 is read as the nearest double, -0x1.6c2d4256ffcc3p+531, by R on
  # Linux, but two spacings away by R on macOS (arm64), and a mode two
  # spacings off is a standard deviation off.
  check_rounded_law(list(alpha = 1e31, beta = 1,
                         gamma = -0x1.6c2d4256ffcc3p+531,
                         hi = 0x1.62e6d72d6dfbp-429,
                         lo = 0x1.8eb1cbd7ece6cp-484),
                    "rtdr")
})

test_that("region A samples a tilt next to the largest double at a large shape", {
  # Once (alpha - 1)/m^2 overflows, region A forms the derivative of its
  # log-density about the mode m from factors that stay in range.  Left of
  # the mode (s < 0) it is -(s/m) (alpha - 1)/(m + s) - 2 s, whose factor
  # (alpha - 1)/(m + s) is about |gamma_norm| (1 - s/m).  That factor exceeds
  # |gamma_norm| by the relative |s|/m, of order 1/sqrt(alpha - 1) at the
  # left tangent, and would overflow within that of the largest double,
  # leaving the tangent's slope infinite and the triple refused.  Grouped as
  # (alpha - 1)/2 / m, about |gamma_norm|/2, times s/(m + s), it does not.
  # The law is Gamma(alpha, rate |gamma|) to double precision: its relative
  # width, 1/sqrt(alpha), is at least 1e-7, some 4.5e8 spacings of the
  # doubles, so the draws have no ties for the KS test.
  for (a in c(1e11, 1e12, 1e14)) {
    for (g in c(-.Machine$double.xmax, -(1 - 1e-7) * .Machine$double.xmax)) {
      info <- sprintf("alpha = %g, gamma = %g", a, g)
      set.seed(12)
      x <- rmhn(2e4, a, 1, g, method = "rtdr")
      expect_false(anyNA(x), info = info)
      expect_gt(gamma_limit_ks_p(x, a, g), 0.001, label = info)
    }
  }
})

test_that("region A is built from the mode where the terms of log f(m) cancel", {
  # At a tilt near -alpha/e the mode m is near e, where the three terms of
  # the raw log-density at the mode, (alpha - 1) log m, -m^2 and
  # gamma_norm m, are of order alpha and cancel to a small log f(m).  That
  # value still carries their rounding error, about alpha DBL_EPSILON (22 at
  # alpha = 1e17, against the contact drop of 1), so the envelope must be
  # built relative to the mode, and its contact points searched from next to
  # it.  Built from the raw values, the draws at 1e17 fail a KS test against
  # Sun et al.'s (2023) Algorithm 3 with p = 0, and at 3e17 the setup stops
  # as a defect.
  dump <- function(a, g) mhn:::.dump_rtdr_envelope_cpp(a, 1, g)
  for (a in c(1e17, 3e17)) {
    g <- -a / exp(1)
    info <- sprintf("alpha = %g", a)
    env <- dump(a, g)
    expect_true(env$centred, info = info)
    expect_true(env$quick_start, info = info)
    set.seed(30)
    x <- rmhn(50, a, 1, g, method = "rtdr")
    expect_false(anyNA(x), info = info)
  }
  # Same law, 2e4 draws from each sampler; the mean is held to the
  # reference within five standard errors.
  a <- 1e17
  g <- -a / exp(1)
  set.seed(31)
  x <- rmhn(2e4, a, 1, g, method = "rtdr")
  set.seed(32)
  y <- rmhn(2e4, a, 1, g, method = "sun")
  expect_gt(suppressWarnings(ks.test(x, y))$p.value, 0.001)
  expect_lt(abs(mean(x) - mhn_ref(a, 1, g)$mean) / (sd(y) / sqrt(2e4)), 5)
})

test_that("a tilt out of the envelope's range is refused, not called a defect", {
  # The envelope is built for sqrt(beta) X, whose law lies near
  # alpha / |gamma_norm|, gamma_norm = gamma / sqrt(beta).  gamma_norm is -Inf
  # in the first three cases; in the fourth the left tangent's slope, a few
  # times |gamma_norm|, is past the largest double.  The law of X is
  # Gamma(alpha, rate -gamma) to double precision in each, and the message
  # must say so.
  for (args in list(c(2, 1e-300, -1e300), c(0.5, 1e-300, -1e300),
                    c(1e-10, 1e-300, -1e300), c(3, 1, -1.7e308))) {
    msg <- tryCatch({
      rmhn(30, args[1], args[2], args[3], method = "rtdr")
      "no error"
    }, error = conditionMessage)
    info <- paste(args, collapse = ", ")
    expect_match(msg, "gamma distribution with shape alpha and rate -gamma",
                 fixed = TRUE, info = info)
    expect_false(grepl("defect", msg, fixed = TRUE), info = info)
  }
})

test_that("method = 'rtdr' samples a shape whose log-axis mode is lost", {
  # For alpha <= 1 RTDR works on the log axis, and under a negative tilt the
  # mode of log(sqrt(beta) X) lies at the log of about alpha / |gamma_norm|.
  # Below a shape of about 2.5e-324 |gamma_norm| that number rounds to 0, so
  # that no envelope can be placed.  All but a fraction of at most about
  # 1e3 alpha of the mass of X then lies below 2^-1075 and rounds to 0, so
  # every draw is 0 except with a probability below 1e-294 in the first four
  # cells and about 1e-16 in the last, and Sun et al.'s (2023) Algorithm 3,
  # which samples the law exactly, gives them, as it does on the default
  # path.  gamma_norm is -1e3, -1e10, -1e50 and -1e30 in the first four, and
  # the scale alpha / |gamma_norm| is 5e-327, 1e-330, 1e-350 and 1e-330.  In
  # the last gamma_norm = -1e300 / 1e-150 overflows to -Inf; below a shape of
  # 1e-15 that is sampled the same way (at alpha = 1e-10 it is refused, see
  # the test above).
  cells <- list(c(5e-324, 1e-300, -1e-147), c(1e-320, 1e-150, -1e-65),
                c(1e-300, 1e-300, -1e-100), c(1e-300, 1, -1e30),
                c(1e-20, 1e-300, -1e300))
  for (args in cells) {
    for (n in c(5, 200)) {
      info <- sprintf("n=%d alpha=%g beta=%g gamma=%g", n, args[1], args[2],
                      args[3])
      set.seed(3)
      x <- expect_silent(rmhn(n, args[1], args[2], args[3], method = "rtdr"))
      expect_identical(x, numeric(n), info = info)
    }
  }
})

test_that("a tiny shape at a negative tilt is sampled, not called a defect", {
  # For alpha < 1 and gamma <= 0 region BC works on y = log(sqrt(beta) x),
  # whose log-density alpha y - e^{2y} + gamma_norm e^y falls off to the left
  # only like alpha y.  The left contact search starts about sqrt(2 / alpha)
  # left of the mode, 1.4e20 at alpha = 1e-40, and at alpha = 1e-310 that
  # offset overflows; the tangent there meets the plateau at
  # p_l = m_g + (gamma_norm u - u^2) / alpha, u = e^{m_g}, which must be
  # formed without terms of the size of the offset.  At (1e-40, 1, -1)
  # u^2 / alpha is 1e-40, so p_l - m_g is -1, and the doubles near
  # m_g = -92 are 1.4e-14 apart: hence the tolerance.
  env <- mhn:::.dump_rtdr_envelope_cpp(1e-40, 1, -1)
  expect_equal(env$p_l - env$mode, -1, tolerance = 1e-12)
  # X exceeds the smallest positive double, 4.9e-324, with probability about
  # 745 alpha, at most 7.5e-38 here, so all 500 draws are 0 except with a
  # probability below 1e-34.  rmhn's default method sends gamma < 0,
  # alpha < 0.1 and n >= 25 to RTDR.
  for (alpha in c(1e-40, 1e-200, 1e-310)) {
    set.seed(11)
    x <- tryCatch(rmhn(500, alpha, 1, -1), error = conditionMessage)
    expect_identical(x, numeric(500), info = sprintf("alpha = %g", alpha))
  }
})

test_that("region BC keeps a draw whose value at beta = 1 is subnormal", {
  # Region BC samples log(sqrt(beta) X) and returns X.  At
  # (1e-3, 1e-300, -100), z = gamma / sqrt(beta) = -1e152, seven in ten draws
  # of sqrt(beta) X lie below the smallest normal double, where they keep
  # only some of their digits or round to 0, although X itself, 1e150 times
  # larger, is representable.  The law of X is Gamma(1e-3, rate 100) to
  # within a relative 1e-300, and X rounds to 0 exactly when it lies below
  # half the smallest subnormal double, 2^-1075, with probability
  # p0 = pgamma(100 2^-1075, 1e-3) = 0.477, and it is subnormal with
  # probability p1 - p0, p1 = pgamma(100 DBL_MIN, 1e-3) = 0.495.  So the
  # numbers of zeros and of draws below DBL_MIN must be binomial, here within
  # 5 standard deviations, and the normal draws must follow the law
  # conditioned on X >= DBL_MIN.  (The subnormal draws crowd onto the
  # smallest few subnormal doubles, and a KS test would see ties there.)
  set.seed(13)
  n <- 5000
  x <- rmhn(n, 1e-3, 1e-300, -100, method = "rtdr")
  xmin <- .Machine$double.xmin
  p0 <- pgamma((100 * xmin) * 2^-53, 1e-3)
  p1 <- pgamma(100 * xmin, 1e-3)
  expect_lt(abs(sum(x == 0) - n * p0), 5 * sqrt(n * p0 * (1 - p0)))
  expect_lt(abs(sum(x < xmin) - n * p1), 5 * sqrt(n * p1 * (1 - p1)))
  u <- (pgamma(100 * x[x >= xmin], 1e-3) - p1) / (1 - p1)
  expect_gt(ks.test(u, "punif")$p.value, 0.001)
})

test_that("region A never reports the shape-axis grid limit as a defect", {
  # On the shape axis the law of X has standard deviation about 1/2 and mode
  # about sqrt(alpha/2).  For alpha in [2^103, 2^105) the doubles at the mode
  # are 1/2 apart, so the standard deviation and the spacing agree to within a
  # relative 1/(8 m), below double precision, and whether the sampler refuses
  # depends on the last bit of the computed standard deviation.  Either
  # outcome is right -- a draw there is the exact draw rounded once -- but a
  # message calling this a defect is not: it is what an envelope built in x
  # itself, whose plateau edges round onto the mode here, would report.  If
  # it samples, every draw lies within 4 (8 standard deviations) of
  # sqrt(alpha/2), which is within 1/4 of the mode and computed to within 1.
  set.seed(4)
  for (alpha in c(1.05e31, 2e31, 3.5e31, 4.05e31)) {
    for (gamma in c(-1, 1)) {
      info <- sprintf("alpha = %g, gamma = %g", alpha, gamma)
      x <- tryCatch(rmhn(200, alpha, 1, gamma, method = "rtdr"),
                    error = conditionMessage)
      if (is.character(x)) {
        expect_match(x, "cannot sample at alpha", info = info)
        expect_false(grepl("defect", x, fixed = TRUE), info = info)
      } else {
        expect_true(all(is.finite(x)), info = info)
        expect_lt(max(abs(x - sqrt(alpha / 2))), 4, label = info)
      }
    }
  }
})

test_that("the tangent pieces dominate the target where their mass is", {
  # Region D's envelope (15) tangent covers (y_K, p_l], with y_K the last
  # secant breakpoint, and region A's left tangent covers [0, p_l].  Their
  # mass lies within a few 1/slope of p_l, but y_K (about -log(gamma) once the
  # secant chain stops at k = 1) and 0 lie far to the left.  A piece anchored
  # at its left end evaluates log h as base + slope * (y - a), two terms of
  # size slope * (p_l - a) -- about 6e13 at a tilt of 1e12 in region D, and at
  # a shape of 1e27 in region A -- that cancel to the O(1) value that matters,
  # which leaves errors of order 1e-2 in log h and lets the hat dip below the
  # target.  Anchored at p_l, where it meets the plateau, the same line has no
  # cancellation.  The right tails are anchored at p_r; built from the offset
  # to the mode, p_r keeps full relative precision, where a p_r rounded to the
  # grid of x at the mode would drop the tail by |slope| times that rounding
  # (1e-2 in log h at alpha = 9.68377e27, gamma = 170.165).  Each piece is
  # evaluated here as its type prescribes, on 2001 points spanning 20/|slope|
  # away from its anchor, against the target measured from the mode; 1e-9
  # allows for rounding in the target only.
  log1p_minus_u <- function(u) {
    ifelse(abs(u) < 1e-4, -u * u * (0.5 - u / 3 + u * u / 4), log1p(u) - u)
  }
  min_gaps <- function(alpha, gamma, log_target) {
    env <- mhn:::.dump_rtdr_envelope_cpp(alpha, 1, gamma)
    expect_true(env$centred)
    types <- vapply(env$pieces, function(p) p$type, integer(1))
    plateau <- which(types == PIECE_PLATEAU)
    gap <- function(piece, s) {
      anchor <- if (piece$type == PIECE_EXP_BOUNDED_RIGHT) piece$b else piece$a
      log_h <- piece$base_log_dens + piece$slope * (s - anchor)
      d <- s - (env$mode - env$origin)     # offset from the mode, exact
      min(log_h - log_target(d, env))
    }
    left <- env$pieces[[plateau - 1L]]
    right <- env$pieces[[plateau + 1L]]
    steps <- seq(0, 20, length.out = 2001)
    c(left = gap(left, left$b - steps / left$slope),
      right = gap(right, right$a - steps / right$slope))
  }
  for (tilt in c(1e10, 1e12, 1e13)) {
    gaps <- min_gaps(0.3, tilt, function(d, env) {
      E <- expm1(d)
      -(env$u_mode^2 * E^2 + 0.3 * (E - d))
    })
    expect_gt(min(gaps), -1e-9, label = sprintf("region D, gamma = %g", tilt))
  }
  for (args in list(c(1e27, 1), c(1e29, 1), c(9.68377e27, 170.165),
                    c(2, 1e13))) {
    alpha <- args[1]
    gaps <- min_gaps(alpha, args[2], function(d, env) {
      (alpha - 1) * log1p_minus_u(d / env$mode) - d^2
    })
    expect_gt(min(gaps), -1e-9,
              label = sprintf("region A, alpha = %g, gamma = %g",
                              alpha, args[2]))
  }
})

test_that("every piece of a region D envelope lies above the target", {
  # Gao & Wang (2025, Theorem 4.4, Eqs. (14) and (15)): the hat dominates the
  # law of y = log(sqrt(beta) X), log g(y) = alpha y - e^{2y} + z e^y up to a
  # constant, on the whole axis, the secants between breakpoints included.
  # Each piece is evaluated as its type prescribes on 401 points of its
  # support, reaching 20 / |slope| into an unbounded tail, with no random
  # numbers.  A secant slope of 0.9 times the right one, whose law is off by
  # D = 0.0034 and passes the sampling tests, puts the hat below g here.
  piece_log_h <- function(p, y) {
    if (p$type == PIECE_PLATEAU) return(rep(p$base_log_dens, length(y)))
    anchor <- if (p$type == PIECE_EXP_BOUNDED_RIGHT) p$b else p$a
    p$base_log_dens + p$slope * (y - anchor)
  }
  piece_grid <- function(p) {
    if (p$type %in% c(PIECE_EXP_LEFT, PIECE_EXP_RIGHT)) {
      p$a - seq(0, 20, length.out = 401) / p$slope
    } else {
      seq(p$a, p$b, length.out = 401)
    }
  }
  # Envelope (14) with five and three breakpoints, envelope (15) with two
  # and five.
  for (cs in list(c(0.3, 5), c(0.45, 4), c(0.3, 20), c(0.05, 8))) {
    alpha <- cs[1]
    z <- cs[2]
    info <- sprintf("alpha = %g, gamma = %g", alpha, z)
    env <- mhn:::.dump_rtdr_envelope_cpp(alpha, 1, z)
    expect_identical(env$region, 2L, info = info)
    expect_false(env$centred, info = info)
    expect_gte(env$K_eff, 2L)
    types <- vapply(env$pieces, function(p) p$type, integer(1))
    expect_true(all(types %in% c(PIECE_PLATEAU, PIECE_EXP_RIGHT,
                                 PIECE_EXP_LEFT, PIECE_EXP_BOUNDED,
                                 PIECE_SECANT)), info = info)
    for (p in env$pieces) {
      y <- piece_grid(p)
      log_g <- alpha * y - exp(2 * y) + z * exp(y)
      expect_gt(min(piece_log_h(p, y) - log_g), -1e-9 * max(1, abs(log_g)),
                label = sprintf("%s, piece type %d on [%g, %g]", info,
                                p$type, min(y), max(y)))
    }
  }
})

test_that("every piece of a region BC envelope lies above the target", {
  # Region BC draws y = log(sqrt(beta) X), log g(y) = alpha y - e^{2y} + z e^y
  # up to a constant, from tangents of T(g) at two contact points and a
  # plateau between them.  For z <= 0 g is log-concave and T = log (T_0); for
  # z > 0 it is only T_{-1/2}-concave (Gao & Wang 2025, Theorem 3.2, whose
  # beta is -z), and the tails are h(y) = exp(base) / (s (y - t) - 1)^2 with
  # s = L'(t) / 2.  A tangent to a concave T(g) lies above it, so each piece
  # dominates g on its support; a tail whose slope is off does not.  Each
  # piece is evaluated as its type prescribes on its support, reaching 20
  # (T_0) or 40 (T_{-1/2}) units of 1/|slope| into a tail, with no random
  # numbers.  A right T_{-1/2} slope of 1.3 times the right one puts the hat
  # 0.037 below g at (0.7, 1, 2).  Past the switch to the centred path
  # (gamma = 1e8 here) the pieces are measured from the mode, and so is g.
  piece_log_h <- function(p, y) {
    switch(as.character(p$type),
      "0" = rep(p$base_log_dens, length(y)),
      "1" = , "2" = p$base_log_dens + p$slope * (y - p$a),
      "5" = , "6" = p$base_log_dens - 2 * log(abs(p$slope * (y - p$a) - 1)),
      stop("unexpected piece type ", p$type))
  }
  piece_grid <- function(p) {
    switch(as.character(p$type),
      "0" = seq(p$a, p$b, length.out = 401),
      "1" = , "2" = p$a - seq(0, 20, length.out = 2001) / p$slope,
      "5" = (p$b - p$aux) - seq(0, 40, length.out = 4001) / p$slope,
      "6" = (p$b + p$aux) - seq(0, 40, length.out = 4001) / p$slope)
  }
  cases <- list(c(0.7, -2), c(0.3, -5), c(0.9, 0), c(0.05, -1),
                c(0.7, 2), c(0.7, 5), c(0.3, 0.5), c(0.6, 10), c(0.99, 3),
                c(0.7, 1e8))
  for (cs in cases) {
    alpha <- cs[1]
    z <- cs[2]
    info <- sprintf("alpha = %g, gamma = %g", alpha, z)
    env <- mhn:::.dump_rtdr_envelope_cpp(alpha, 1, z)
    expect_identical(env$region, 1L, info = info)
    expect_identical(env$centred, z > 1e7, info = info)
    types <- vapply(env$pieces, function(p) p$type, integer(1))
    expect_identical(types, if (z > 0) {
      c(PIECE_TNEGHALF_LEFT, PIECE_PLATEAU, PIECE_TNEGHALF_RIGHT)
    } else {
      c(PIECE_EXP_LEFT, PIECE_PLATEAU, PIECE_EXP_RIGHT)
    }, info = info)
    log_target <- if (env$centred) {
      function(s) {
        E <- expm1(s)
        -(env$u_mode^2 * E^2 + alpha * (E - s))
      }
    } else {
      function(y) alpha * y - exp(2 * y) + z * exp(y)
    }
    for (p in env$pieces) {
      y <- piece_grid(p)
      log_g <- log_target(y)
      expect_gt(min(piece_log_h(p, y) - log_g), -1e-9 * max(1, abs(log_g)),
                label = sprintf("%s, piece type %d on [%g, %g]", info,
                                p$type, min(y), max(y)))
    }
  }
})

test_that("a left-tail draw past the doubles is accepted with its g / h", {
  # On the log axis y = log(sqrt(beta) X) the leftmost piece of regions BC
  # and D is h(y) = h(a) exp(s (y - a)), and a draw from it, a + log(u) / s,
  # passes the largest double, and is -Inf, where |log u| > s DBL_MAX: for a
  # shape near 2^-1022 (s = alpha in region D) at u below e^-4.  There
  # log g(y) - log h(y) is (alpha - s) y + phi(e^y) + const, with
  # phi(t) = gamma t - t^2 negligible for y below -745, so the draws past
  # -DBL_MAX are accepted with probability (s / alpha) g(y0) / h(y0) on
  # average, y0 = -DBL_MAX (exactly g / h for s = alpha), and all of them
  # are 0.  Taken as -Inf <= -Inf they would all be accepted, and 0 would
  # come out too often.  The acceptance rate of the sampler is
  # int g / int h, with
  # int g = (1 + w) / alpha on the log axis, w = alpha K0(gamma), since
  # t^alpha = 1 to double precision there, and 1 / (1 + w) from pmhn() at the
  # smallest positive double.  (Region BC, which a tiny shape reaches for
  # gamma <= 0, has s = alpha to double precision wherever a draw can pass
  # the doubles, and all its draws are 0 there.)
  a1 <- 2^-1022
  # Each case has a seed of its own: the first two laws are nearly the same.
  for (cs in list(c(a1, 0.5, 23), c(1.25 * a1, 0.5, 31), c(a1, 0.9, 37))) {
    a <- cs[1]
    z <- uniroot(function(z) pmhn(5e-324, a, 1, z) - cs[2], c(50, 56),
                 tol = 1e-10)$root
    info <- sprintf("alpha = %g, z = %g", a, z)
    env <- mhn:::.dump_rtdr_envelope_cpp(a, 1, z)
    log_H <- max(env$piece_log_area) +
      log(sum(exp(env$piece_log_area - max(env$piece_log_area))))
    q0 <- pmhn(5e-324, a, 1, z)
    accept <- exp(-log(a) - log(q0) - log_H)
    n <- 1e6
    set.seed(cs[3])
    x <- mhn:::.rmhn_rtdr_cpp(n, a, 1, z)
    tries <- n + attr(x, "rtdr_retries")
    expect_true(abs(n / tries - accept) < 5 * sqrt(accept * (1 - accept) / tries),
                info = paste(info, "| acceptance", n / tries, "against", accept))
    expect_true(abs(mean(x == 0) - q0) < 5 * sqrt(q0 * (1 - q0) / n),
                info = paste(info, "| zeros", mean(x == 0), "against", q0))
  }
})
