# Tests for the public rmhn() dispatcher.
#
# Coverage:
#   - general goodness-of-fit (KS, mean, var) across regimes
#   - method-agreement (RTDR vs Sun) where both apply
#   - the law on either side of the boundary of regions (c) and (d), that
#     is of the RTDR setups BC and D
#   - density-vs-histogram chi-square consistency
#   - vectorisation reproducibility (cache-on-change <-> scalar)
#   - edge cases (n=0, n<0, an n that is not a count, elements past n,
#     NA propagation, a non-finite parameter -> NaN)
#   - method argument validation
#   - Robert (1995) truncated-normal kernel cross-check
#
# Lower-level kernels (RTDR per region, Sun Algo 1/3) are covered by
# test-rmhn-rtdr-regions.R and test-rmhn-sun.R; this file focuses on the
# integration surface.  The references, mhn_ref() and mhn_cdf(), are in
# helper-mhn-cdf.R and are formed from the kernel of the law alone.
#
# (a)-(d) are the vignette's parameter regions; A, BC and D the RTDR setups
# ((b) and (c) share BC).

# =====================================================================
# General goodness-of-fit (method-agnostic, default "auto")
# =====================================================================

run_default <- function(alpha, beta, gamma, n = 10000, seed = 1L) {
  set.seed(seed)
  x <- rmhn(n, alpha, beta, gamma)
  ref <- mhn_ref(alpha, beta, gamma)
  list(
    samples     = x,
    ref         = ref,
    ks          = ks.test(x, ref$cdf),
    sample_mean = mean(x),
    sample_var  = var(x),
    theory_mean = ref$mean,
    theory_var  = ref$var
  )
}

# Each label names the sampler that the default method uses for 10000 draws.
regimes <- list(
  list(alpha = 1.5, beta = 1, gamma = -2,  label = "RTDR region A (alpha=1.5, gamma=-2)"),
  list(alpha = 5,   beta = 1, gamma = 10,   label = "Sun Algorithm 1 (alpha=5, gamma=10)"),
  list(alpha = 0.7, beta = 1, gamma =  0,   label = "sqrt-Gamma closed form (alpha=0.7, gamma=0)"),
  list(alpha = 0.7, beta = 1, gamma =  2,   label = "RTDR region b (alpha=0.7, gamma=2)"),
  list(alpha = 0.3, beta = 1, gamma = -5,   label = "RTDR region c (alpha=0.3, gamma=-5)"),
  list(alpha = 0.3, beta = 1, gamma =  5,   label = "RTDR region d (alpha=0.3, gamma=5)"),
  # Non-unit beta: a setup that drops a factor of sqrt(beta) still samples the
  # right law at beta = 1.
  list(alpha = 1.5, beta = 4,    gamma = -2,  label = "RTDR region A, beta=4"),
  list(alpha = 5,   beta = 0.25, gamma = 5,   label = "Sun Algorithm 1, beta=0.25"),
  list(alpha = 0.7, beta = 9,    gamma = 0,   label = "sqrt-Gamma closed form, beta=9"),
  list(alpha = 0.3, beta = 0.01, gamma = 5,   label = "RTDR region d, beta=0.01")
)

for (r in regimes) {
  test_that(paste("default rmhn KS+moments:", r$label), {
    skip_on_cran()
    out <- run_default(r$alpha, r$beta, r$gamma, n = 10000)
    expect_gt(out$ks$p.value, 0.001)
    # mean: 5 sigma tolerance
    se_mean <- sqrt(out$theory_var / 10000)
    expect_lt(abs(out$sample_mean - out$theory_mean), 5 * se_mean)
    # variance: five standard errors of a sample variance
    expect_var_near(out$sample_var, out$ref, 10000, info = r$label)
  })
}

# =====================================================================
# Each sampler on its own, by the method that reaches it at any n
# =====================================================================

sampler_points <- list(
  list(p = c(1.5, 1, -2),  method = "rtdr", label = "RTDR region A"),
  list(p = c(0.7, 1, 2),   method = "rtdr", label = "RTDR region BC, gamma > 0"),
  list(p = c(0.3, 1, -5),  method = "rtdr", label = "RTDR region BC, gamma < 0"),
  list(p = c(0.3, 1, 5),   method = "rtdr", label = "RTDR region D, envelope (14)"),
  list(p = c(0.05, 1, 8),  method = "rtdr", label = "RTDR region D, envelope (15)"),
  list(p = c(5, 1, 10),    method = "sun",  label = "Sun Algorithm 1, Normal proposal"),
  list(p = c(5, 1, 1),     method = "sun",  label = "Sun Algorithm 1, sqrt-Gamma proposal"),
  list(p = c(2, 1, -3),    method = "sun",  label = "Sun Algorithm 3, alpha > 1"),
  list(p = c(0.5, 1, -10), method = "sun",  label = "Sun Algorithm 3, alpha < 1"),
  list(p = c(1, 1, -3),    method = "auto", label = "truncated normal (Robert 1995)"),
  list(p = c(0.7, 1, 0),   method = "auto", label = "sqrt-Gamma closed form")
)

test_that("each sampler keeps its law, on CRAN too", {
  # 2000 draws a point under a fixed seed against the reference CDF; a
  # correct sampler falls below p = 1e-4 for one seed in 1e4.  These run on
  # every platform CRAN checks, where the larger tests are skipped.
  for (i in seq_along(sampler_points)) {
    sp <- sampler_points[[i]]
    set.seed(100 + i)
    x <- rmhn(2000, sp$p[1], sp$p[2], sp$p[3], method = sp$method)
    expect_gt(ks.test(x, mhn_cdf(sp$p[1], sp$p[2], sp$p[3]))$p.value, 1e-4,
              label = sp$label)
  }
})

test_that("each sampler keeps its law and moments at 2e5 draws", {
  skip_on_cran()
  # A chi-square test on 100 bins of near-equal reference probability, and
  # the mean and the variance to five of their standard errors.  2e5 draws
  # tell apart laws a few thousandths apart in distribution function, which
  # 1e4 draws do not.
  for (i in seq_along(sampler_points)) {
    sp <- sampler_points[[i]]
    ref <- mhn_ref(sp$p[1], sp$p[2], sp$p[3])
    set.seed(200 + i)
    x <- rmhn(2e5, sp$p[1], sp$p[2], sp$p[3], method = sp$method)
    expect_gt(chisq_ref(x, ref), 1e-4, label = sp$label)
    expect_lt(abs(mean(x) - ref$mean) / sqrt(ref$var / 2e5), 5,
              label = sp$label)
    expect_var_near(var(x), ref, 2e5, info = sp$label)
  }
})

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
    set.seed(1); xr <- rmhn(10000, r$alpha, r$beta, r$gamma, method = "rtdr")
    set.seed(2); xs <- rmhn(10000, r$alpha, r$beta, r$gamma, method = "sun")
    expect_gt(suppressWarnings(ks.test(xr, xs))$p.value, 0.001)
  })
}

# =====================================================================
# region (c)/(d) boundary continuity
# =====================================================================

test_that("region (c)/(d) boundary: the law holds on either side", {
  # The boundary of regions (c) and (d), where RTDR passes from its setup BC
  # to its setup D, lies at
  # gamma_star = 2 (1 - sqrt(1 - 2 alpha)).  On each side, 5% from it, 2e5
  # draws are held to the reference law by a chi-square test on 100 bins,
  # and their mean to five standard errors.
  alpha <- 0.3
  gamma_star <- 2 * (1 - sqrt(1 - 2 * alpha))
  for (g in gamma_star * c(0.95, 1.05)) {
    ref <- mhn_ref(alpha, 1, g)
    set.seed(1)
    x <- rmhn(2e5, alpha = alpha, gamma = g, method = "rtdr")
    info <- sprintf("gamma = %.4f", g)
    expect_gt(chisq_ref(x, ref), 1e-4, label = info)
    expect_lt(abs(mean(x) - ref$mean) / sqrt(ref$var / 2e5), 5, label = info)
  }
})

# =====================================================================
# density-vs-histogram chi-square consistency
# =====================================================================

test_that("chi-square: histogram of rmhn matches the kernel reference", {
  # 100 nearly equiprobable bins whose probabilities come from the kernel
  # reference of helper-mhn-cdf.R, which does not use dmhn or pmhn, so every
  # bin expects 500 draws and none needs to be dropped or pooled.
  set.seed(1)
  x <- rmhn(50000, 2, 1, 0.5)
  expect_gt(chisq_ref(x, mhn_ref(2, 1, 0.5)), 0.001)
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

test_that("vectorised rmhn equals the element-wise loop on the general samplers", {
  # The test above only reaches the closed forms (beta = 1, gamma = 0).  This
  # set reaches Sun Algorithms 1 and 3, RTDR in its three setups, and the
  # truncated-normal and sqrt-Gamma closed forms, with non-unit beta,
  # consecutive repeated triples (a reused setup) and method = "sun".  The
  # setups consume no random numbers, so under one seed the vectorised call
  # and the loop of single draws must agree bit for bit.  n equals the
  # parameter length: with shorter recycled vectors, auto may choose a
  # different sampler wherever its choice depends on the draws per setup
  # (gamma < 0 with alpha < 1, or with alpha >= 1e3 where RTDR's setup walks
  # its contact points in, and gamma > 0 with alpha > 1 at a tilt up to
  # 2 sqrt(alpha)), by design, since more draws share each setup.
  al <- c(2, 2, 3, 5, 1, 1.5, 0.3, 0.3, 0.7, 10, 4)
  be <- c(1, 1, 4, 0.25, 1, 9, 1, 1, 1, 0.04, 2)
  for (s in c(-1, 1)) {
    ga <- s * c(2, 2, 2, 0.5, 1, 3, 2, 2, 5, 1, 0)
    for (m in c("auto", "rtdr", if (s < 0) "sun")) {
      set.seed(1)
      v <- rmhn(length(al), al, be, ga, method = m)
      set.seed(1)
      w <- mapply(function(a, b, g) rmhn(1, a, b, g, method = m), al, be, ga)
      expect_identical(v, w, label = paste("method", m, "sign", s))
    }
  }
})

test_that("recycled parameters give the draws of the same parameters written out", {
  # Recycled vectors repeat their triples with period lcm of their lengths, so
  # one setup per position within the period can serve all the draws at that
  # position; written out in full, the same triples alternate and a setup is
  # rebuilt at every element.  The setups consume no random numbers, so both
  # calls must return the same draws.  Only choices that do not depend on the
  # number of draws per setup are compared: method = "rtdr" and "sun", and
  # "auto" with gamma > 0 where alpha < 1 or the tilt gamma / sqrt(beta)
  # exceeds 2 sqrt(alpha).  Lengths 2 and 3 give period 6.  The second case
  # repeats a triple inside a period, the third has the same triple at two
  # consecutive residues, the fourth mixes a closed form (alpha = 1) with a
  # rejection sampler, and in the last the NA element must stay NaN in place.
  n <- 60
  cases <- list(
    list(a = c(1.5, 3),       b = c(1, 4, 0.5), g = -2,       m = "rtdr"),
    list(a = c(2.5, 3.5, 2.5), b = 2,           g = c(-1, -3), m = "sun"),
    list(a = c(2, 2, 3),      b = 1,            g = -1,       m = "sun"),
    list(a = c(1, 1.5),       b = c(1, 2, 3),   g = 6,        m = "auto"),
    list(a = c(0.3, 0.7),     b = 1,            g = c(2, 5, 1), m = "auto"),
    list(a = c(2, 5),         b = c(1, 0.25),   g = 6,        m = "auto"),
    list(a = c(0.5, 2),       b = 1,            g = c(-1, NA, 3), m = "rtdr")
  )
  draw <- function(...) {
    # An NA parameter draws NaN with a warning, which is not what is compared.
    withCallingHandlers(rmhn(...), warning = function(w) {
      if (identical(conditionMessage(w), "NAs produced")) {
        invokeRestart("muffleWarning")
      }
    })
  }
  for (cs in cases) {
    set.seed(9)
    v <- draw(n, cs$a, cs$b, cs$g, method = cs$m)
    set.seed(9)
    w <- draw(n, rep_len(cs$a, n), rep_len(cs$b, n), rep_len(cs$g, n),
              method = cs$m)
    expect_identical(v, w, label = paste("method", cs$m))
  }
})

test_that("an exact method string and its partial match draw the same", {
  # An exact method string skips match.arg; a partial match, the full choice
  # vector, and invalid values still go through it, so partial matches are
  # accepted and invalid values raise match.arg's error.
  set.seed(2); a <- rmhn(20, 3, 1, -2, method = "sun")
  set.seed(2); b <- rmhn(20, 3, 1, -2, method = "su")
  expect_identical(a, b)
  set.seed(2); a <- rmhn(20, 0.5, 1, 1)
  set.seed(2); b <- rmhn(20, 0.5, 1, 1, method = c("auto", "rtdr", "sun"))
  expect_identical(a, b)
  # The error is match.arg's; its call says so in any language, where its
  # message is translated.
  for (m in list(NA_character_, c("sun", "rtdr"), "AUTO")) {
    e <- expect_error(rmhn(3, method = m))
    expect_identical(conditionCall(e), quote(match.arg(method)),
                     label = deparse(m))
  }
})

test_that("the sqrt-Gamma closed form draws the far left tail at a small shape", {
  # X = sqrt(T) with T ~ Gamma(alpha/2, rate beta).  At alpha = 0.01, R's
  # rgamma returns T as exactly 0 about 2.4% of the time, and sqrt(0) = 0,
  # although the law puts only 5.8e-4 of its mass below the smallest
  # positive double.  The exact masses are, to first order in the incomplete
  # gamma function (the next term is smaller by x^2 = 1e-400),
  #   P(X < x) = x^alpha / Gamma(1 + alpha/2),
  # 0.0100286967623 at x = 1e-200 and 5.82333864448e-4 at x = 2^-1075, below
  # which X rounds to 0 (256-bit Rmpfr).  The fractions of 2e5 draws must lie
  # within 5 standard errors of them.
  n <- 2e5
  set.seed(4)
  x <- rmhn(n, 0.01, 1, 0)
  ref <- 0.0100286967623
  expect_lt(abs(mean(x < 1e-200) - ref), 5 * sqrt(ref * (1 - ref) / n))
  p0 <- 5.82333864448e-4
  expect_lt(abs(mean(x == 0) - p0), 5 * sqrt(p0 * (1 - p0) / n))
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

test_that("an n that is not a count is refused with a true message", {
  # rmhn() takes a single n through as.double(), which gives NA for NULL and
  # for a string that is not a number (the latter after a coercion warning,
  # as in rnorm).  test-input-validation.R tests the rest of the rule.
  msg <- "n must be a finite, non-negative number"
  for (n in list(NA, NaN, Inf, 2^31, NULL, "a")) {
    expect_error(suppressWarnings(rmhn(n, alpha = 2)), msg,
                 label = deparse(n))
  }
  expect_error(rmhn(-Inf, alpha = 2), "n must be non-negative")
  expect_identical(rmhn(0.5, alpha = 2), numeric(0))
  expect_identical(rmhn(numeric(0), alpha = 2), numeric(0))
})

test_that("parameter elements past n are ignored, as rnorm ignores them", {
  # Recycled to length n, the draws use only the first n elements of each
  # vector, and only those can draw NaN.
  set.seed(1)
  x <- rmhn(1, alpha = c(2, -1))
  set.seed(1)
  expect_identical(x, rmhn(1, alpha = 2))
  set.seed(1)
  y <- rmhn(2, alpha = 2, beta = c(1, 3, -1), gamma = c(0.5, -2, 0))
  set.seed(1)
  expect_identical(y, rmhn(2, alpha = 2, beta = c(1, 3), gamma = c(0.5, -2)))
  expect_warning(x <- rmhn(2, alpha = c(2, -1)), "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE))
  expect_warning(x <- rmhn(3, beta = c(1, 3, -1)), "NAs produced")
  expect_identical(is.nan(x), c(FALSE, FALSE, TRUE))
})

test_that("NA in alpha draws NaN in place, with a warning", {
  set.seed(1)
  expect_warning(x <- rmhn(5, alpha = c(1, NA, 2, NA, 3)), "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_true(all(x[!is.na(x)] > 0))
})

test_that("NA in beta and gamma draws NaN in place, with a warning", {
  set.seed(1)
  expect_warning(x_b <- rmhn(4, alpha = 2, beta = c(1, NA, 1, NA)),
                 "NAs produced")
  expect_identical(is.nan(x_b), c(FALSE, TRUE, FALSE, TRUE))

  expect_warning(x_g <- rmhn(4, alpha = 2, gamma = c(0.5, NA, 0.5, NA)),
                 "NAs produced")
  expect_identical(is.nan(x_g), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("a non-finite alpha, beta, or gamma draws NaN", {
  expect_warning(x <- rmhn(5, alpha = 2, gamma = c(0.5, Inf, 0.5, -Inf, NaN)),
                 "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, FALSE, TRUE, TRUE))

  # The guard covers all three parameters, not gamma alone: an infinite alpha
  # or beta would otherwise reach the sampler and come back as something the
  # caller could not tell from a draw.
  expect_warning(x_a <- rmhn(5, alpha = c(2, Inf, 2, NaN, -Inf)),
                 "NAs produced")
  expect_identical(is.nan(x_a), c(FALSE, TRUE, FALSE, TRUE, TRUE))
  expect_warning(x_b <- rmhn(5, alpha = 2, beta = c(1, Inf, 1, NaN, -Inf)),
                 "NAs produced")
  expect_identical(is.nan(x_b), c(FALSE, TRUE, FALSE, TRUE, TRUE))
})

test_that("a zero or negative alpha or beta draws NaN, as in rgamma", {
  for (args in list(list(alpha = -1), list(alpha = 0), list(beta = 0),
                    list(beta = -0.5))) {
    expect_warning(x <- do.call(rmhn, c(list(5), args)), "NAs produced")
    expect_identical(x, rep(NaN, 5))
  }
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
  e <- expect_error(rmhn(10, method = "invalid"))
  expect_identical(conditionCall(e), quote(match.arg(method)))
})

test_that("method=sun NA elements skip prescan", {
  # The NA elements draw NaN with a warning; they do not reach the check that
  # refuses alpha < 1 with gamma > 0.
  set.seed(1)
  expect_warning(
    x <- rmhn(5, alpha = c(2, NA, 2, NA, 3), gamma = 1, method = "sun"),
    "NAs produced"
  )
  expect_identical(is.nan(x), c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

# =====================================================================
# Special-case dispatch produces reasonable output
# =====================================================================

test_that("Special-case half-normal (alpha=1, gamma=0) is |N(0, 1/(2 beta))|", {
  # The half-normal is drawn as |sigma Z|, sigma = 1 / sqrt(2 beta), with Z
  # from R's normal generator: the law's own definition.  Under one seed the
  # draws are therefore those of abs(rnorm()), at any beta; a wrong scale
  # changes every one of them.
  for (b in c(0.5, 1, 3, 1e-300, 1e300)) {
    set.seed(1)
    x <- rmhn(2000, alpha = 1, beta = b, gamma = 0)
    set.seed(1)
    expect_identical(x, abs(rnorm(2000, 0, 1 / sqrt(2 * b))),
                     info = sprintf("beta = %g", b))
  }
})

test_that("alpha = 1 with gamma / (2 beta) beyond the doubles is Exp(-gamma)", {
  # Under a negative tilt with gamma / (2 beta) past the largest double, the
  # truncated normal is exponential with rate -gamma to double precision and
  # is drawn as E / (-gamma), with E from R's exponential generator.  Under
  # one seed the draws are therefore those of rexp(n) / -gamma.
  for (cs in list(c(1e-300, -1e300), c(1e-200, -1e250))) {
    set.seed(1)
    x <- rmhn(2000, alpha = 1, beta = cs[1], gamma = cs[2])
    set.seed(1)
    expect_identical(x, rexp(2000) / -cs[2],
                     info = sprintf("beta = %g, gamma = %g", cs[1], cs[2]))
  }
})

test_that("Special-case sqrt-Gamma (alpha=3, gamma=0) matches sqrt(Gamma(1.5, 1))", {
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

test_that("the closed forms keep their law at extreme beta and tilt", {
  # Each law is written through a variable that has a fixed law:
  #   alpha = 1: sqrt(2 beta) X ~ N(mu, 1) truncated to (0, Inf), with
  #     mu = gamma / sqrt(2 beta) (Sun et al. 2023, Lemma 6b), |N(0, 1)| at
  #     gamma = 0;
  #   gamma = 0: beta X^2 ~ Gamma(alpha / 2, 1) (Sun et al. 2023, Lemma 6a);
  #   alpha = 1 with gamma / sqrt(beta) = -1e160, or with gamma / (2 beta)
  #     beyond the largest double: |gamma| X ~ Exp(1) to within a relative
  #     (sqrt(beta) / gamma)^2, the normal truncated far to the right of its
  #     mean.
  # beta = 1.7e308 makes 2 beta overflow, beta = 5e-324 and 1e-310 make
  # 1 / beta overflow, and at beta = 1e308 the Gamma variate of rate beta
  # falls below DBL_MIN more often than not, where it is drawn again from its
  # conditional law.  The variables are formed without leaving the double
  # range, and compared with their laws at the 0.001 level.
  hn <- function(q) 2 * pnorm(q) - 1
  mu <- 1e155 / (sqrt(2) * sqrt(1.7e308))
  tn <- function(q) (pnorm(q - mu) - pnorm(-mu)) / pnorm(mu)
  cases <- list(
    list(p = c(1,   1.7e308, 0),     f = function(x) sqrt(2) * sqrt(1.7e308) * x, cdf = hn),
    list(p = c(1,   1.7e308, 1e155), f = function(x) sqrt(2) * sqrt(1.7e308) * x, cdf = tn),
    list(p = c(3,   5e-324,  0),     f = function(x) (sqrt(5e-324) * x)^2,
         cdf = function(q) pgamma(q, 1.5)),
    list(p = c(0.5, 1e-310,  0),     f = function(x) (sqrt(1e-310) * x)^2,
         cdf = function(q) pgamma(q, 0.25)),
    list(p = c(1.5, 1e308,   0),     f = function(x) (sqrt(1e308) * x)^2,
         cdf = function(q) pgamma(q, 0.75)),
    list(p = c(1,   1,       -1e160), f = function(x) 1e160 * x,  cdf = pexp),
    list(p = c(1,   1e-10,   -1e300), f = function(x) 1e300 * x,  cdf = pexp)
  )
  draws <- lapply(cases, function(cs) {
    set.seed(21)
    rmhn(2000, cs$p[1], cs$p[2], cs$p[3])
  })
  for (k in seq_along(cases)) {
    expect_true(all(is.finite(draws[[k]]) & draws[[k]] > 0),
                info = paste(cases[[k]]$p, collapse = ", "))
  }
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    expect_gt(ks.test(cs$f(draws[[k]]), cs$cdf)$p.value, 0.001,
              label = paste("KS p at", paste(cs$p, collapse = ", ")))
  }
})

test_that("the closed forms stop where the law lies beyond the largest double", {
  # At alpha = 1 the law is N(gamma / (2 beta), 1 / (2 beta)) truncated to
  # x > 0; at (1, 1e-310, 1) and (1, 1e-300, 1e10) its mean is 5e309, with a
  # standard deviation below 3.2e161.  At gamma = 0, beta X^2 is
  # Gamma(alpha / 2, 1), so that X is about sqrt(alpha / (2 beta)) = 3e311 at
  # (1e300, 5e-324, 0).  alpha within sqrt(eps) of 1, and a tilt below
  # sqrt(eps) sqrt(beta), are drawn in closed form too.  No double can
  # represent these laws, and the call must stop with the reason, as the
  # general samplers do at (2, 1e-310, 1) (see test-rmhn-sun.R), rather than
  # return Inf.  Just inside the range, at (1, 1e-300, 1) and
  # (1e300, 1e-300, 0), the draws, about 5e299 and 7e299, must be finite.
  beyond <- list(c(1, 1e-310, 1), c(1 + 1e-10, 1e-310, 1), c(1, 1e-300, 1e10),
                 c(1e300, 5e-324, 0), c(1e300, 5e-324, -1e-310))
  for (p in beyond) {
    expect_error(rmhn(3, p[1], p[2], p[3]),
                 "beyond the largest double precision number",
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
  set.seed(22)
  for (p in list(c(1, 1e-300, 1), c(1e300, 1e-300, 0))) {
    x <- rmhn(3, p[1], p[2], p[3])
    expect_true(all(is.finite(x) & x > 1e299),
                info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
})

test_that("the truncated normal is centred at its exact mean at a steep tilt", {
  # At (1, 3, 6 * 2^48 + 0.25), a tilt gamma / sqrt(beta) of 9.7e14, the law
  # is N(mu, 1/6) truncated far below its mean, with mu = gamma / 6 =
  # 2^48 + 1/24 exactly.  The doubles near mu are 1/16 apart, so mu itself
  # rounds to 2^48 + 1/16, a shift of 0.051 standard deviations that would
  # move every draw alike.  A draw is the exact draw rounded once, and the
  # mean of such draws is mu to far below one standard error, since the
  # standard deviation spans 6.5 spacings.  x - 2^48 is exact, so the sample
  # mean of the offsets is compared with 1/24, within five standard errors
  # (0.0065 for 1e5 draws).  alpha = 1 + 1e-10 is drawn from the same law.
  n <- 1e5
  for (a in c(1, 1 + 1e-10)) {
    set.seed(31)
    x <- rmhn(n, a, 3, 6 * 2^48 + 0.25)
    expect_lt(abs(mean(x - 2^48) - 1 / 24), 5 / sqrt(6 * n),
              label = sprintf("offset of the mean at alpha = %.10g", a))
    expect_lt(abs(sd(x) * sqrt(6) - 1), 0.02)
  }
})

test_that("a negligible tilt at a huge shape is centred at the exact mode", {
  # With gamma = 0 and beta = 1/2, X^2 / 2 is Gamma(alpha / 2, 1).  R's
  # rgamma builds a Gamma variate of shape a from sqrt(a - 1/2) rounded to
  # double precision, which moves every draw alike by up to 2^-52 sqrt(a)
  # standard deviations, 0.18 at a = 2^99.  At alpha = 4 M^2 with
  # M = 2^49 + 253 * 2^23, a = 2 M^2 and that rounding is 0.0625 (half a unit
  # in the last place of sqrt(a), found by a 200-bit scan), which would shift
  # X by 0.088, 0.125 standard deviations.  The mode is sqrt(alpha - 1),
  # 2M - 2^-51, the mean lies within 1/(8M) of it, and the standard
  # deviation is 1/sqrt(2) to a relative 1e-30; x - 2M is exact.  The sample
  # mean of 1e5 draws must lie within five standard errors (0.011) of 2M
  # under every method.
  M <- 2^49 + 253 * 2^23
  n <- 1e5
  for (m in c("auto", "sun", "rtdr")) {
    set.seed(32)
    x <- rmhn(n, 4 * M^2, 0.5, 0, method = m)
    expect_lt(abs(mean(x - 2 * M)), 5 / sqrt(2 * n), label = paste("method", m))
    expect_lt(abs(sd(x) * sqrt(2) - 1), 0.02, label = paste("method", m))
  }
})

# =====================================================================
# Robert (1995) truncated-normal kernel cross-check
# =====================================================================

test_that("rtnorm_robert_zero matches truncated normal CDF: mu > 0", {
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

test_that("rmhn at alpha = 1 with mu < 0 keeps the truncated normal's law", {
  # There rmhn draws by Robert (1995), Proposition 2.3: an exponential
  # proposal of rate alpha* accepted with probability
  # exp(-(z - alpha*)^2 / 2), where the truncation point a = -mu / sigma
  # sets how often a draw is rejected (a = 0.35 and 2 here).  With 2e5 draws
  # the KS test tells apart a law whose acceptance exponent is 1/4 instead of
  # 1/2 (D = 0.027 and p = 8e-127 at a = 0.35).
  for (g in c(-0.5, -2 * sqrt(2))) {
    mu <- g / 2
    s <- 1 / sqrt(2)
    log_q0 <- pnorm(0, mu, s, lower.tail = FALSE, log.p = TRUE)
    cdf <- function(q) {
      -expm1(pnorm(q, mu, s, lower.tail = FALSE, log.p = TRUE) - log_q0)
    }
    set.seed(1)
    x <- rmhn(2e5, alpha = 1, beta = 1, gamma = g)
    # R's uniform generator takes 2^32 values, so 2e5 draws carry a few
    # ties (3 here), which move the KS statistic by at most 3 / 2e5; only
    # the warning about them is muffled.
    ks <- withCallingHandlers(ks.test(x, cdf), warning = function(w) {
      if (grepl("ties", conditionMessage(w))) invokeRestart("muffleWarning")
    })
    expect_gt(ks$p.value, 1e-4, label = sprintf("KS p at gamma = %g", g))
  }
})

test_that("rtnorm_robert_zero edge cases", {
  expect_identical(mhn:::.rtnorm_robert_cpp(0, mu = 0, sigma = 1), numeric(0))
  expect_error(mhn:::.rtnorm_robert_cpp(-1, mu = 0, sigma = 1), "non-negative")
  expect_error(mhn:::.rtnorm_robert_cpp(5, mu = 0, sigma = 0), "sigma must be positive")
  expect_error(mhn:::.rtnorm_robert_cpp(5, mu = Inf, sigma = 1), "mu must be finite")
  # A finite mu whose standardised bound -mu / sigma overflows would leave the
  # exponential proposal without a rate, and the direct entry point refuses
  # it.  Within rmhn that bound, -gamma / sqrt(2 beta), is finite wherever
  # mu = gamma / (2 beta) is.
  expect_error(mhn:::.rtnorm_robert_cpp(5, mu = -1e300, sigma = 1e-10),
               "mu / sigma must be finite")
})

# =====================================================================
# auto-path n-dependent dispatch verification
# Strategy: rebuild_cache is RNG-deterministic, so under the same seed
# `method = "auto"` and the method that auto picks must produce
# bit-identical outputs.
# =====================================================================

# Whether rmhn(n, a, b, g) with the default method draws what method m draws
# under the same seed.
same_as <- function(n, a, b, g, m) {
  set.seed(7L); x <- rmhn(n, a, b, g, method = "auto")
  set.seed(7L); y <- rmhn(n, a, b, g, method = m)
  identical(x, y)
}

test_that("auto -> Sun Algorithm 1 for gamma>0 & alpha>1 with few draws per setup or a strong tilt", {
  # Below the threshold of the shape band and tilt half (see the threshold
  # tests below), and at any n once gamma / sqrt(beta) > 2 sqrt(alpha).
  for (n_val in c(1L, 4L)) expect_true(same_as(n_val, 5, 1, 2, "sun"))
  for (n_val in c(1L, 24L)) expect_true(same_as(n_val, 50, 1, 2, "sun"))
  expect_true(same_as(24L, 4, 1, 4, "sun"))                 # tilt 2 sqrt(alpha), below 25
  for (n_val in c(1L, 100L, 1000L)) expect_true(same_as(n_val, 5, 1, 10, "sun"))
  expect_true(same_as(100L, 4, 1, 4 * (1 + 1e-9), "sun"))   # just past 2 sqrt(alpha)
  expect_true(same_as(100L, 5, 100, 50, "sun"))             # tilt 5 > 2 sqrt(5)
  expect_true(same_as(100L, 2^60, 1, 1, "sun"))             # above kRgammaShapeMax
})

test_that("auto -> RTDR for gamma>0 & alpha>1 with a mild tilt and enough draws per setup", {
  for (n_val in c(5L, 100L, 1000L)) expect_true(same_as(n_val, 5, 1, 2, "rtdr"))
  for (n_val in c(25L, 1000L)) expect_true(same_as(n_val, 50, 1, 2, "rtdr"))
  expect_true(same_as(25L, 4, 1, 4, "rtdr"))                # tilt exactly 2 sqrt(alpha)
  expect_true(same_as(100L, 5, 0.01, 0.2, "rtdr"))          # tilt 2 at beta = 0.01
})

test_that("auto -> RTDR for gamma>0 & alpha<1 (Sun Algorithm 2 unimplemented)", {
  for (n_val in c(1L, 100L, 1000L)) {
    set.seed(7L); a <- rmhn(n_val, alpha = 0.3, gamma = 5, method = "auto")
    set.seed(7L); r <- rmhn(n_val, alpha = 0.3, gamma = 5, method = "rtdr")
    expect_identical(a, r)
  }
})

test_that("auto -> Sun Algorithm 3 for gamma<0 & alpha<1 with few draws per setup", {
  # Below 5 draws per setup for 0.1 <= alpha < 1, and below 25 for
  # alpha < 0.1.
  for (n_val in c(1L, 4L)) expect_true(same_as(n_val, 0.5, 1, -10, "sun"))
  for (n_val in c(1L, 10L, 24L)) expect_true(same_as(n_val, 0.01, 1, -10, "sun"))
})

test_that("auto thresholds for gamma>0 & alpha>1 follow the shape band and the tilt half", {
  # RTDR from T draws per setup, T = 5 / 25 (alpha < 10), 25 / 50
  # (10 <= alpha < 1e3), 50 / 100 (1e3 <= alpha < 1e8) and 100 / 1000
  # (alpha >= 1e8), for gamma / sqrt(beta) <= sqrt(alpha) / up to
  # 2 sqrt(alpha); the tilt sqrt(alpha) itself belongs to the lower half.
  # Each band is tested inside and on either side of its edges.  At beta = 4
  # gamma / sqrt(beta) is formed exactly, so the tilts sqrt(alpha) and
  # 2 sqrt(alpha) are met exactly too.
  below <- 1 - 1e-12
  bands <- list(c(5, 5, 25), c(10 * below, 5, 25),
                c(10, 25, 50), c(50, 25, 50), c(1e3 * below, 25, 50),
                c(1e3, 50, 100), c(1e4, 50, 100), c(1e8 * below, 50, 100),
                c(1e8, 100, 1000), c(1e10, 100, 1000))
  for (bd in bands) {
    a <- bd[1]
    for (h in 1:2) {
      t <- bd[1 + h]
      for (zr in if (h == 1) c(0.5, 1) else c(1 + 1e-9, 1.5, 2)) {
        for (b in c(1, 4)) {
          g <- zr * sqrt(a) * sqrt(b)
          info <- sprintf("alpha=%.17g z/sqrt(alpha)=%g beta=%g", a, zr, b)
          expect_true(same_as(as.integer(t) - 1L, a, b, g, "sun"), label = info)
          expect_true(same_as(as.integer(t), a, b, g, "rtdr"), label = info)
        }
      }
    }
  }
})

test_that("auto -> RTDR for gamma<0 from 5 draws per setup (25 for alpha < 0.1; any n for 1 < alpha < 1e3)", {
  for (n_val in c(5L, 100L)) expect_true(same_as(n_val, 0.5, 1, -10, "rtdr"))
  for (n_val in c(25L, 1000L)) expect_true(same_as(n_val, 0.01, 1, -10, "rtdr"))
  for (a in c(1.5, 5, 10, 100, 999)) for (n_val in c(1L, 5L, 25L, 1000L))
    expect_true(same_as(n_val, a, 1, -10, "rtdr"),
                label = sprintf("alpha=%g n=%d", a, n_val))
})

test_that("auto for gamma<0 & alpha >= 1e3: Sun Algorithm 3 below 10 / 25 / 50 draws per setup", {
  # From alpha = 1e3 RTDR's region-A setup walks its contact points in from
  # m/2 and m + max(1, m/2), and Algorithm 3 wins small batches: RTDR from
  # 10 draws per setup (alpha < 1e5), 25 (alpha < 1e8) or 50, and at any n
  # below alpha = 1e3.  Each band is tested inside and on either side of its
  # edges.
  below <- 1 - 1e-12
  bands <- list(c(1e3 * below, 1), c(1e3, 10), c(1e4, 10), c(1e5 * below, 10),
                c(1e5, 25), c(1e7, 25), c(1e8 * below, 25), c(1e8, 50),
                c(1e12, 50))
  for (p in bands) {
    t <- as.integer(p[2])
    for (z in c(-2, -1e4)) {
      info <- sprintf("alpha=%.17g z=%g", p[1], z)
      if (t > 1L) {
        expect_true(same_as(1L, p[1], 1, z, "sun"), label = info)
        expect_true(same_as(t - 1L, p[1], 1, z, "sun"), label = info)
      }
      expect_true(same_as(t, p[1], 1, z, "rtdr"), label = info)
      expect_true(same_as(t, p[1], 4, 2 * z, "rtdr"), label = info)
    }
  }
  # Recycled shapes keep one setup per residue, each serving n / 2 draws
  # here: 9 at n = 18 (Algorithm 3) and 10 at n = 20 (RTDR).
  set.seed(3); v <- rmhn(18, c(1e4, 2e4), 1, -2)
  set.seed(3); s <- rmhn(18, c(1e4, 2e4), 1, -2, method = "sun")
  expect_identical(v, s)
  set.seed(3); v <- rmhn(20, c(1e4, 2e4), 1, -2)
  set.seed(3); r <- rmhn(20, c(1e4, 2e4), 1, -2, method = "rtdr")
  expect_identical(v, r)
})

test_that("auto takes RTDR's quick setup start into account", {
  # RTDR's region-A setup starts its contact points next to the mode on its
  # centred path (from a shape of about 2.4e12 at moderate tilts, and of
  # about 2.3e13 at gamma / sqrt(beta) = -alpha/e, where the terms of the
  # log-density at the mode cancel) and at gamma / sqrt(beta) <
  # -1e10 sqrt(alpha); there RTDR serves gamma < 0 at any n, and gamma > 0
  # from 5 / 25 draws per setup.  The choice must follow the setup's own
  # decision, which the envelope dump reports.
  dump <- function(a, b, g) mhn:::.dump_rtdr_envelope_cpp(a, b, g)
  cases <- list(c(1e4, 1, -1e12 * (1 + 1e-9)), c(1e4, 1, -1e12 * (1 - 1e-9)),
                c(1e4, 100, -1e13 * (1 + 1e-9)), c(1e13, 1, -2), c(1e13, 1, -1e12),
                c(1e13, 1, -1e16), c(1e12, 1, -1e8), c(2.8e12, 1, -2),
                c(1e13, 0.01, -0.2), c(1e14, 1, -100), c(1e17, 1, -1e17 / exp(1)),
                c(1e17, 4, -2e17 / exp(1)))
  for (p in cases) {
    q <- dump(p[1], p[2], p[3])$quick_start
    info <- sprintf("alpha=%g beta=%g gamma=%g quick=%s", p[1], p[2], p[3], q)
    expect_true(same_as(1L, p[1], p[2], p[3], if (q) "rtdr" else "sun"), label = info)
  }
  expect_true(dump(1e4, 1, -1e12 * (1 + 1e-9))$quick_start)
  expect_false(dump(1e4, 1, -1e12 * (1 - 1e-9))$quick_start)
  expect_true(dump(1e13, 1, -2)$centred)
  expect_true(dump(1e17, 1, -1e17 / exp(1))$quick_start)
  expect_true(same_as(49L, 1e17, 1, -1e17 / exp(1), "rtdr"))
  # At (1e13, -1e12) the setup walks in: the mode is 10, and the magnitudes
  # of the three terms of the log-density there, (alpha - 1) log m, m^2 and
  # |gamma| m, sum to 3.3e13, below the switch at 0.01 / DBL_EPSILON =
  # 4.5e13.  The band's 50 draws apply.
  expect_false(dump(1e13, 1, -1e12)$quick_start)
  expect_true(same_as(49L, 1e13, 1, -1e12, "sun"))
  expect_true(same_as(50L, 1e13, 1, -1e12, "rtdr"))
  # gamma > 0: 5 / 25 on the centred path, against 100 / 1000 below it
  for (p in list(c(1e13, 0.5, 5), c(1e13, 1.5, 25), c(1e14, 1, 5))) {
    g <- p[2] * sqrt(p[1])
    info <- sprintf("alpha=%g z/sqrt(alpha)=%g", p[1], p[2])
    expect_true(dump(p[1], 1, g)$centred, label = info)
    expect_true(same_as(as.integer(p[3]) - 1L, p[1], 1, g, "sun"), label = info)
    expect_true(same_as(as.integer(p[3]), p[1], 1, g, "rtdr"), label = info)
  }
  expect_false(dump(1e12, 1, 0.5e6)$centred)
  expect_true(same_as(99L, 1e12, 1, 0.5e6, "sun"))
  expect_true(same_as(100L, 1e12, 1, 0.5e6, "rtdr"))
  # One setup per draw (S = 1), as in a vectorised Gibbs update: RTDR for
  # gamma < 0 on the centred path, Algorithm 3 where the setup walks in, and
  # Algorithm 1 for gamma > 0 even on the centred path.
  vec_same_as <- function(al, g, m) {
    set.seed(7L); x <- rmhn(length(al), al, 1, g)
    set.seed(7L); y <- rmhn(length(al), al, 1, g, method = m)
    identical(x, y)
  }
  expect_true(vec_same_as(1e13 * (1 + (0:2) * 1e-12), -2, "rtdr"))
  expect_true(vec_same_as(1e4 * (1 + (0:2) * 1e-12), -2, "sun"))
  expect_true(vec_same_as(1e13 * (1 + (0:2) * 1e-12), 0.5 * sqrt(1e13), "sun"))
})

test_that("auto follows the setup's start on either side of a flip of its centred path", {
  # The setup's switch to the centred path is a comparison of a computed
  # double, so it flips between two adjacent shapes.  Found here by bisection
  # on the doubles from the setup's own report, at a fixed tilt: the default
  # method must agree with the setup on both of them.
  flip <- function(z, lo, hi) {
    centred <- function(a) mhn:::.dump_rtdr_envelope_cpp(a, 1, z)$centred
    stopifnot(!centred(lo), centred(hi))
    repeat {
      mid <- lo + (hi - lo) / 2
      if (mid <= lo || mid >= hi) break
      if (centred(mid)) hi <- mid else lo <- mid
    }
    c(lo, hi)
  }
  dump <- function(a, b, g) mhn:::.dump_rtdr_envelope_cpp(a, b, g)
  # gamma < 0: RTDR at any n on the centred side, from 50 draws below it
  ab <- flip(-2, 2.8e12, 1e13)
  expect_identical(ab[2], ab[1] + 2^-11)                # adjacent doubles
  for (b in c(1, 4)) {
    g <- -2 * sqrt(b)
    expect_false(dump(ab[1], b, g)$quick_start, label = sprintf("beta=%g", b))
    expect_true(dump(ab[2], b, g)$quick_start, label = sprintf("beta=%g", b))
    expect_true(same_as(1L, ab[1], b, g, "sun"))
    expect_true(same_as(49L, ab[1], b, g, "sun"))
    expect_true(same_as(50L, ab[1], b, g, "rtdr"))
    expect_true(same_as(1L, ab[2], b, g, "rtdr"))
  }
  # gamma > 0 at a tilt below sqrt(alpha): RTDR from 5 draws on the centred
  # side, from 100 below it
  z <- 1e6
  ab <- flip(z, 2.8e12, 1e13)
  expect_identical(ab[2], ab[1] + 2^-11)
  expect_true(same_as(5L, ab[1], 1, z, "sun"))
  expect_true(same_as(99L, ab[1], 1, z, "sun"))
  expect_true(same_as(100L, ab[1], 1, z, "rtdr"))
  expect_true(same_as(4L, ab[2], 1, z, "sun"))
  expect_true(same_as(5L, ab[2], 1, z, "rtdr"))
})

test_that("auto -> Sun Algorithm 3 for gamma<0 past a tilt of 1e150 at any n", {
  # Past a tilt of 1e150 Algorithm 3 draws the Gamma(alpha, rate |gamma|)
  # limit of the law and accepts it at once, so the default method takes it
  # there whatever the number of draws.  The cutoff is on the tilt, so the
  # last case, with a small gamma at a smaller beta, is routed the same way;
  # the fourth has alpha < 0.1, where the batch cutoff is 25.  The density
  # of Y = |gamma| X is that of Gamma(alpha, 1) times
  # exp(-beta Y^2 / gamma^2), which is 1 to within 1e-290 for Y < 1000,
  # beyond which Gamma(alpha, 1) has no mass to speak of; Y is compared with
  # Gamma(alpha, 1) at the 0.001 level.
  cases <- list(c(2, 1, -1e160), c(0.5, 1, -1e160), c(5, 1e10, -1e200),
                c(0.05, 1, -1e160), c(2, 1e-300, -1.1))
  for (p in cases) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    for (n_val in c(1L, 25L, 1000L)) {
      set.seed(7L); a <- rmhn(n_val, p[1], p[2], p[3], method = "auto")
      set.seed(7L); s <- rmhn(n_val, p[1], p[2], p[3], method = "sun")
      expect_identical(a, s, info = info)
    }
  }
  for (p in cases) {
    set.seed(8L)
    x <- rmhn(2000, p[1], p[2], p[3])
    expect_gt(ks.test(abs(p[3]) * x, "pgamma", p[1])$p.value, 0.001,
              label = sprintf("KS p at alpha=%g beta=%g gamma=%g",
                              p[1], p[2], p[3]))
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
  # max(na, nb, ng) >= n -> samples_per_setup = 1 -> Sun Algorithm 3
  # (gamma < 0, alpha < 1) and Sun Algorithm 1 (gamma > 0, alpha > 1, mild
  # tilt), while the scalar call of the same size counts 5 draws per setup
  # and takes RTDR.
  set.seed(7L); a <- rmhn(5L, alpha = rep(0.5, 5), gamma = -10, method = "auto")
  set.seed(7L); s <- rmhn(5L, alpha = rep(0.5, 5), gamma = -10, method = "sun")
  expect_identical(a, s)
  expect_true(same_as(5L, 0.5, 1, -10, "rtdr"))
  set.seed(7L); a <- rmhn(5L, alpha = rep(5, 5), gamma = 2, method = "auto")
  set.seed(7L); s <- rmhn(5L, alpha = rep(5, 5), gamma = 2, method = "sun")
  expect_identical(a, s)
})

test_that("recycled parameters count their draws per kept setup", {
  # With the setups kept one per residue, each serves n / period draws; past
  # 256 residues one setup is kept and rebuilt as the triples change, and the
  # call counts one draw per setup.
  set.seed(3); v <- rmhn(60, alpha = c(0.5, 0.6), gamma = -10)
  set.seed(3); r <- rmhn(60, alpha = c(0.5, 0.6), gamma = -10, method = "rtdr")
  expect_identical(v, r)
  al <- seq(0.3, 0.9, length.out = 300)
  set.seed(3); v <- rmhn(3000, alpha = al, gamma = -10)
  set.seed(3); s <- rmhn(3000, alpha = al, gamma = -10, method = "sun")
  expect_identical(v, s)
  # The period is the lcm of the lengths, 6 here, not the longest length:
  # 24 draws count 4 per setup (Algorithm 3), 30 count 5 (RTDR).
  set.seed(3); v <- rmhn(24, alpha = c(0.5, 0.6), beta = c(1, 2, 3), gamma = -10)
  set.seed(3); s <- rmhn(24, alpha = c(0.5, 0.6), beta = c(1, 2, 3), gamma = -10,
                         method = "sun")
  expect_identical(v, s)
  set.seed(3); v <- rmhn(30, alpha = c(0.5, 0.6), beta = c(1, 2, 3), gamma = -10)
  set.seed(3); r <- rmhn(30, alpha = c(0.5, 0.6), beta = c(1, 2, 3), gamma = -10,
                         method = "rtdr")
  expect_identical(v, r)
})

test_that("the samplers are scale-equivariant in beta", {
  # If X ~ MHN(alpha, beta, gamma) then sqrt(beta) X ~ MHN(alpha, 1, gamma/sqrt(beta)):
  # the family depends on the tilt only through Delta = gamma/sqrt(beta)
  # (a change of variables in the density of Sun et al. 2023, Lemma 1a).  A
  # setup that has lost a factor of sqrt(beta) still samples the right law at
  # beta = 1 and drifts away from it elsewhere.  Comparing the two sides with a
  # two-sample Kolmogorov-Smirnov test needs no reference implementation at
  # all.
  cases <- list(
    c(1.5, 4,    -2), c(5,   0.25, 0.5), c(0.7, 9,    0),
    c(0.3, 0.01, 5),  c(3,   4,    -4),  c(2,   100,  50),
    c(0.45, 0.5, 1),  c(10,  0.04, -1)
  )
  for (par in cases) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    delta <- gamma / sqrt(beta)
    for (method in c("auto", "rtdr")) {
      set.seed(101)
      x <- rmhn(6000, alpha, beta, gamma, method = method) * sqrt(beta)
      set.seed(202)
      y <- rmhn(6000, alpha, 1, delta, method = method)
      info <- sprintf("alpha=%g beta=%g gamma=%g method=%s",
                      alpha, beta, gamma, method)
      expect_gt(suppressWarnings(stats::ks.test(x, y)$p.value), 0.001,
                label = info)
    }
  }
})
