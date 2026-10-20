# Numerical-stability tests: parameter regions where a direct formula loses
# precision or overflows, and the behaviour of the samplers, their warnings
# and their error messages at extreme parameters.  The grids elsewhere in the
# suite stop well short of most of these regions.

test_that("mhn_mode stays exact for strongly negative gamma", {
  # Sun et al. (2023) Lemma 3b as printed evaluates
  #   (gamma + sqrt(gamma^2 + 8 beta (alpha-1))) / (4 beta),
  # whose numerator cancels for gamma < 0.  Gao & Wang (2025, Eq. (9)) also give
  # the conjugate form, which is what the package uses.
  exact <- function(alpha, beta, gamma) {
    2 * (alpha - 1) / (sqrt(gamma^2 + 8 * beta * (alpha - 1)) - gamma)
  }
  for (alpha in c(1.1, 2, 10)) {
    for (beta in c(0.1, 1, 10)) {
      for (gamma in -10^c(4, 6, 8, 10)) {
        got <- mhn_mode(alpha, beta, gamma)
        expect_gt(got, 0)
        expect_equal(got, exact(alpha, beta, gamma), tolerance = 1e-13)
      }
    }
  }
  # The specific point at which the literal form collapses to exactly zero.
  expect_equal(mhn_mode(1.1, 1, -1e8), 1e-9, tolerance = 1e-12)
})

test_that("dmhn and pmhn remain finite and positive at alpha = 2 with large gamma", {
  # Psi[1, z] = 1 + sqrt(pi) z exp(z^2/4) Phi(z/sqrt(2)); the middle factor
  # overflows a double for z > 53.3, which, formed as written, drives the
  # normalising constant to +Inf and the whole density to zero.
  for (gamma in c(53.3, 60, 100, 1000)) {
    mode_x <- (gamma + sqrt(gamma^2 + 8)) / 4    # alpha = 2, beta = 1
    d <- dmhn(mode_x, alpha = 2, beta = 1, gamma = gamma)
    expect_true(is.finite(d))
    expect_gt(d, 0)
    # The density at the mode of an MHN(2, 1, gamma) tends to a constant as the
    # tilt grows, because the law approaches a normal of variance 1/2.
    expect_equal(d, 1 / sqrt(pi), tolerance = 1e-3)

    p <- pmhn(mode_x, alpha = 2, beta = 1, gamma = gamma)
    expect_true(is.finite(p))
    expect_gt(p, 0)
    expect_lt(p, 1)
  }
  # No collapse across that overflow boundary.  One unit above the mode
  # is sqrt(2) standard deviations once the law has become normal with variance
  # 1/2, so every value must sit at Phi(sqrt(2)) and stay there; a one-sided
  # bound on the change would admit the collapse to zero this block guards.
  gs <- c(50, 53, 53.3, 55, 60)
  ms <- (gs + sqrt(gs^2 + 8)) / 4
  v <- vapply(seq_along(gs), function(i) pmhn(ms[i] + 1, 2, 1, gs[i]), 0)
  expect_lt(max(abs(v - stats::pnorm(sqrt(2)))), 1e-3)
  expect_lt(max(abs(diff(v))), 1e-4)
})

test_that("Sun Algorithm 1 keeps delta_opt inside (0, beta)", {
  # The first part of Sun et al.'s (2023) Theorem 1 (printed "a.") requires
  # delta in (0, beta).  Evaluated literally the optimum delta_opt of the
  # second part (printed "a." too) is a difference of two quantities of size
  # gamma^2 and leaves the interval for large tilt; a negative value would
  # make the proposal a Gamma with negative scale, and every draw NaN.
  for (alpha in c(1.5, 2, 10, 100)) {
    for (beta in c(0.01, 1, 100)) {
      for (gamma in 10^seq(2, 12, by = 0.5)) {
        s <- mhn:::.dump_sun_algo1_cpp(alpha, beta, gamma)
        expect_gt(s$delta_opt, 0)
        expect_lt(s$delta_opt, beta)
        expect_gt(s$mu_opt, 0)
      }
    }
  }
})

test_that("rmhn draws are finite at extreme positive tilt", {
  # Triples at which the literal Theorem 1 formula for delta_opt loses every
  # digit: its sign then depends on how the compiler evaluates the
  # subtraction, and a negative value makes the proposal a Gamma with negative
  # scale and every draw NaN, with no error or warning.
  cases <- list(c(1.5, 0.01, 3162.3), c(2, 0.01, 5623.4),
                c(1.5, 1, 56234.1),   c(10, 0.01, 31622.8))
  for (p in cases) {
    set.seed(42)
    x <- rmhn(200, alpha = p[1], beta = p[2], gamma = p[3])
    expect_false(anyNA(x))
    expect_true(all(is.finite(x)))
    expect_true(all(x > 0))
    # At this tilt the law concentrates tightly at its mode, so the exact mean
    # sits there too.  The sample mean is judged in its own standard errors,
    # since any relative band would be thousands of them wide.
    mode_x <- (p[3] + sqrt(p[3]^2 + 8 * p[2] * (p[1] - 1))) / (4 * p[2])
    expect_equal(mhn_mean(p[1], p[2], p[3]), mode_x, tolerance = 1e-6)
    z <- (mean(x) - mhn_mean(p[1], p[2], p[3])) /
         sqrt(mhn_var(p[1], p[2], p[3]) / length(x))
    expect_lt(abs(z), 5)
  }
})

test_that("pmhn stays bounded in cost and range for a strong tilt", {
  # With K the Lemma 10(d) length, the CDF series runs to index 2K + 1,
  # which grows like 1.44 z^2: at alpha = 2.5, beta = 1 it is about 9.2e5
  # terms at |gamma| = 800, 1.2e6 at 900, 1.4e8 at 1e4 and 1.4e10 at 1e5.
  # Past a ceiling of 1e6 terms the series returns its NaN sentinel without
  # summing, and pmhn() takes quadrature instead, which answers in a bounded
  # number of evaluations.  The
  # ceiling is asserted through that sentinel rather than the wall clock, which
  # under valgrind or on a loaded machine varies by an order of magnitude.
  #
  # For gamma > 0 every term is positive and there is no cancellation guard,
  # so at this shape and q the sentinel can come only from the ceiling:
  # without it the series would run at gamma = 900 (stopping early, since for
  # a positive tilt it also stops once the terms fall below 2^-54 of Psi) and
  # return a number.  At gamma = 800, under the ceiling, it does return one,
  # and it matches the quadrature.  gamma = 1e4 is refused on z^2 alone,
  # before any truncation length is formed.  The CDF at the mode is 0.5 to
  # within the skewness of the near-normal law there.
  for (gamma in c(900, 1e4)) {
    q <- mhn_mode(2.5, 1, gamma)
    expect_true(is.nan(mhn:::.pmhn_force(q, 2.5, 1, gamma, "series")),
                info = sprintf("gamma=%g", gamma))
    expect_equal(pmhn(q, 2.5, 1, gamma), 0.5, tolerance = 1e-3)
  }
  q <- mhn_mode(2.5, 1, 800)
  series <- mhn:::.pmhn_force(q, 2.5, 1, 800, "series")
  expect_false(is.nan(series))
  expect_equal(series, mhn:::.pmhn_force(q, 2.5, 1, 800, "integrate"),
               tolerance = 1e-8)

  # For gamma < 0 the cancellation guard also returns the sentinel at these
  # tilts, so there it shows only that the series declines, not which exit
  # fired: at gamma = -800, under the ceiling, it is the guard.
  for (gamma in c(-800, -2000, -10000, -1e5)) {
    v <- pmhn(0.5, alpha = 2.5, beta = 1, gamma = gamma)
    expect_true(is.finite(v))
    expect_gte(v, 0)
    expect_lte(v, 1)
    expect_true(is.nan(mhn:::.pmhn_force(0.5, 2.5, 1, gamma, "series")))
  }
  # For gamma < 0 the alternating series is cancellation-prone and returns the
  # NaN sentinel once it cannot meet tolerance, which is how it hands over to
  # quadrature.  The contract is therefore one-sided: wherever the series does
  # produce a number, it must agree with quadrature.
  for (gamma in c(-1, -5, -20, -50, -200, -600)) {
    series <- mhn:::.pmhn_force(0.3, 2.5, 1, gamma, "series")
    quad   <- mhn:::.pmhn_force(0.3, 2.5, 1, gamma, "integrate")
    expect_true(is.finite(quad))
    if (!is.nan(series)) expect_equal(series, quad, tolerance = 1e-10)
    # Either way the public entry point must return a valid probability.
    p <- pmhn(0.3, 2.5, 1, gamma)
    expect_true(is.finite(p) && p >= 0 && p <= 1)
  }
})

test_that("a sampler that exhausts its retry budget reports once, not per draw", {
  # A rejection sampler that rejects max_retries proposals for one draw
  # returns it as NaN, which no successful draw is; the dispatcher counts those
  # and raises a single warning naming the count.  rmhn() leaves the budget at
  # 1000, which no supported law exhausts, so the internal entry point is
  # given a smaller one.  With a budget of 0 every draw of a rejection sampler
  # fails; the NA parameter in the middle draws NaN too, as rnorm() does, but
  # that is not a failure of the sampler: it is reported by its own warning
  # and must not be counted.
  set.seed(1)
  w <- testthat::capture_warnings(
    x <- mhn:::.rmhn_cpp(3L, c(2, NA, 2), 1, 1, "auto", max_retries = 0L))
  expect_length(w, 2L)
  expect_identical(w[1], "NAs produced")
  expect_match(w[2], "exhausted its retry budget for 2 of 3 draws",
               fixed = TRUE)
  expect_true(all(is.nan(x)))
  # The same for each sampler: Algorithm 1 above, RTDR and Algorithm 3 here.
  for (args in list(list(0.5, 1, 2, "auto"), list(2, 1, -1, "sun"),
                    list(2, 1, -1, "rtdr"))) {
    w <- testthat::capture_warnings(
      x <- mhn:::.rmhn_cpp(4L, args[[1]], args[[2]], args[[3]], args[[4]],
                           max_retries = 0L))
    expect_length(w, 1L)
    expect_match(w, "for 4 of 4 draws", fixed = TRUE)
    expect_true(all(is.nan(x)))
  }
  # With a budget of one proposal, a draw fails exactly when its proposal is
  # rejected, about once in five times for RTDR at (0.5, 1, 2), whose
  # acceptance rate is 0.81: the single warning must carry the count of NaN
  # draws, whatever it is.
  set.seed(2)
  w <- testthat::capture_warnings(
    x <- mhn:::.rmhn_cpp(200L, 0.5, 1, 2, "rtdr", max_retries = 1L))
  k <- sum(is.nan(x))
  expect_gt(k, 0L)
  expect_lt(k, 200L)
  expect_length(w, 1L)
  expect_match(w, sprintf("for %d of 200 draws", k), fixed = TRUE)
  # The closed forms have no budget to exhaust: the half-normal at
  # (1, 1, 0) draws finite values even with a budget of 0.
  expect_silent(x <- mhn:::.rmhn_cpp(3L, 1, 1, 0, "auto", max_retries = 0L))
  expect_true(all(is.finite(x)))
  expect_error(mhn:::.rmhn_cpp(1L, 2, 1, 1, "auto", max_retries = -1L),
               "must be non-negative")
  # A working sampler with the default budget stays silent.
  set.seed(1)
  expect_silent(rmhn(500, alpha = 2, beta = 1, gamma = 1))
})

test_that("region D envelope pieces stay finite and usable at high tilt", {
  # An envelope piece integrates to exp(base) * expm1(slope * width) / slope,
  # and its inverse CDF needs log(1 + u (exp(slope * width) - 1)).  Formed by
  # exponentiating slope * width first, both overflow past about 709: the area
  # becomes +Inf, which turns the whole piece-selection cumulative sum into
  # NaN, and the inverse CDF returns +Inf for every u, so the widest piece
  # proposes nothing inside its own support.  The package forms both on the
  # log scale (log_abs_expm1 and log1p_u_expm1 in mhn_stable.h).
  for (alpha in c(0.01, 0.05, 0.3, 0.45)) {
    for (gamma in c(50, 100, 1000, 1e4)) {
      env <- mhn:::.dump_rtdr_envelope_cpp(alpha, 1, gamma)
      expect_true(all(is.finite(env$piece_log_area)),
                  info = sprintf("alpha=%g gamma=%g", alpha, gamma))
      # A piece whose inverse CDF overflows proposes +Inf for every uniform
      # draw, so it never yields an accepted sample.  Drawing through the
      # kernel exercises every piece and would surface that as a non-finite
      # value or an exhausted retry budget.
      set.seed(3)
      x <- mhn:::.rmhn_rtdr_cpp(300, alpha, 1, gamma)
      expect_true(all(is.finite(as.numeric(x))),
                  info = sprintf("alpha=%g gamma=%g", alpha, gamma))
      expect_true(all(as.numeric(x) > 0))
    }
  }
})

test_that("rmhn is unbiased in region D at high tilt", {
  skip_on_cran()
  # rmhn() under its default method = "auto" sends alpha < 1 with gamma > 0 to
  # RTDR, region D for alpha = 0.3 at these tilts.  An envelope piece whose
  # area overflows (see the test above) leaves the draws biased, which the
  # Kolmogorov-Smirnov test and the five-standard-error bound on the mean
  # below both detect.
  for (gamma in c(50, 100, 200, 1000)) {
    set.seed(11)
    x <- rmhn(4000, alpha = 0.3, beta = 1, gamma = gamma)
    expect_false(anyNA(x))
    ref <- mhn_ref(0.3, 1, gamma)
    expect_gt(stats::ks.test(x, ref$cdf)$p.value, 0.001)
    z <- (mean(x) - ref$mean) / sqrt(ref$var / length(x))
    expect_lt(abs(z), 5)
  }
})

test_that("Sun Algorithm 1 meets the Theorem 2e acceptance bound", {
  # Theorem 2e of Sun et al. (2023) guarantees acceptance >= 0.8 for alpha >= 4
  # and any gamma > 0.  Which proposal is used depends on comparing K_2 with the
  # envelope constant K_1 of the first part of Theorem 1 (printed "a."),
  # whose power has base
  # sqrt(beta) (alpha-1) / (2 beta mu - gamma).  Putting (alpha-1) under the
  # square root instead understates log K_1 by (alpha-1)/2 * log(alpha-1),
  # selects the Normal proposal over most of the small-tilt half of the region,
  # and drops acceptance to about 0.70.
  for (alpha in c(4, 10)) {
    for (gamma in c(0.05, 0.5, 2, 20)) {
      set.seed(7)
      x <- mhn:::.rmhn_sun_algo1_cpp(5000, alpha, 1, gamma)
      acceptance <- 5000 / (5000 + attr(x, "sun_retries"))
      expect_gte(acceptance, 0.8)
    }
  }
})

test_that("the Psi quadrature resolves the peak at strongly negative gamma", {
  # For gamma < 0 the normalising constant is evaluated by quadrature of
  # u^(alpha-1) exp(-u^2 - |z| u).  That integrand is close to a
  # Gamma(alpha, rate |z|) kernel, whose mass lies within a few multiples of
  # (alpha + 1)/|z| of the origin.  An upper limit near |z|/2 -- the peak of
  # exp(-u^2 + |z| u), the opposite sign -- would leave a fixed-order rule
  # sampling only the flat tail, and log Psi would collapse to -Inf from |z|
  # about 600, taking dmhn, pmhn and the moments with it.
  reference <- function(alpha, gamma) {
    cc <- abs(gamma)
    I <- stats::integrate(function(s) exp((alpha - 1) * log(s) - s^2 / cc^2 - s),
                          0, 60, rel.tol = 1e-12, subdivisions = 4000)$value
    log(2) - alpha * log(cc) + log(I)
  }
  for (alpha in c(0.3, 0.8, 1.5, 2.5, 20)) {
    for (gamma in c(-100, -600, -1000, -1e4, -1e5)) {
      got <- mhn:::.mhn_log_normalizing_const(alpha, 1, gamma)
      expect_true(is.finite(got), info = sprintf("alpha=%g gamma=%g", alpha, gamma))
      expect_equal(got, reference(alpha, gamma), tolerance = 1e-8,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
  # Downstream consumers must come back finite and sane.
  expect_true(is.finite(dmhn(0.005, 2.5, 1, -1000)))
  expect_gt(dmhn(0.005, 2.5, 1, -1000), 0)
  expect_true(is.finite(mhn_mean(2.5, 1, -1000)))
  expect_gt(mhn_mean(2.5, 1, -1000), 0)
  expect_gt(mhn_var(2.5, 1, -1000), 0)
  # At this tilt the law is Gamma(alpha, rate |gamma|) up to the factor
  # exp(-beta x^2), which moves F by a relative beta E[X^2] at most,
  # alpha (alpha + 1) / gamma^2 = 8.8e-8 here.
  expect_equal(pmhn(1e-4, 2.5, 1, -1e4), pgamma(1, 2.5), tolerance = 1e-6)
})

test_that("the special-case tests are made on the scale-free tilt", {
  # The shape depends on gamma only through Delta = gamma / sqrt(beta)
  # (sqrt(beta) X ~ MHN(alpha, 1, Delta), from the density in Sun et al. 2023,
  # Lemma 1a), so "is the tilt negligible?" cannot be a test on |gamma| alone.
  # At beta = 4e-18 a gamma of 1e-8 is Delta = 5 -- an ordinary tilt -- yet a
  # test on |gamma| would send it to the gamma = 0 closed form of rmhn, which
  # draws from that law for a Delta within sqrt(eps) of 0.  The other
  # exported functions take gamma = 0 only where it is 0 exactly, and the
  # mean below checks that they see the tilt at such a beta.
  alpha <- 2; gamma <- 1e-8; beta <- (gamma / 5)^2      # Delta = 5
  expect_false(mhn:::.is_sqrt_gamma(gamma, beta))
  expect_true(mhn:::.is_sqrt_gamma(gamma, 1e40))        # genuinely negligible

  # Reference in the scaled variable u = x * sqrt(beta), where the tilt is 5.
  s <- 1 / sqrt(beta)
  lk <- function(u) (alpha - 1) * log(u) - u * u + 5 * u
  pk <- max(lk(seq(0.01, 20, length.out = 2e4)))
  Z <- stats::integrate(function(u) exp(lk(u) - pk), 0, 50, rel.tol = 1e-12)$value
  M <- stats::integrate(function(u) u * exp(lk(u) - pk), 0, 50, rel.tol = 1e-12)$value
  expect_equal(mhn_mean(alpha, beta, gamma), s * M / Z, tolerance = 1e-8)

  # And beta = 1, where gamma and Delta coincide, is unchanged.
  expect_true(mhn:::.is_sqrt_gamma(1e-10))
  expect_false(mhn:::.is_sqrt_gamma(0.01))
})

test_that("the alpha = 1 sampler stays inside its support at extreme tilt", {
  # Robert (1995) draws y = a + e with a = -mu/sigma and returns mu + sigma*y.
  # That subtracts two nearly equal quantities once mu is far below zero, but
  # the mu cancels analytically: mu + sigma(a + e) = sigma*e.  The literal
  # form rounds each draw to the spacing of the doubles near |mu|, about
  # 7e-9 at gamma = -1e8 against a mean of 1e-8: of 2e5 such draws some 40%
  # are exactly zero, outside (0, Inf), and there are only about a dozen
  # distinct values.
  for (gamma in c(-1e4, -1e6, -1e8, -1e10)) {
    set.seed(31)
    x <- rmhn(2e4, alpha = 1, beta = 1, gamma = gamma)
    expect_true(all(x > 0), info = sprintf("gamma=%g", gamma))
    expect_gt(length(unique(x)), 1e4)
    # As the tilt grows the law tends to Exponential(rate = |gamma|).
    expect_close(mean(x), 1 / abs(gamma), 0.05,
                 info = sprintf("gamma=%g", gamma))
  }
})

test_that("Psi at alpha = 1 avoids the z^2/4 cancellation", {
  # log Psi[1/2, z] = log(2) + log(pi)/2 + z^2/4 + log Phi(z/sqrt(2)).  For
  # z << 0 the last two terms cancel to every significant digit; at z = -1e8
  # this form is off by 0.73 in a quantity of size 17.7.  Writing
  # Phi(z/sqrt2) = exp(-z^2/4) erfcx(|z|/2) / 2 removes the large terms.
  reference <- function(gamma) {
    cc <- abs(gamma)
    I <- stats::integrate(function(s) exp(-s^2 / cc^2 - s), 0, Inf,
                          rel.tol = 1e-13)$value
    log(2) - log(cc) + log(I)
  }
  for (gamma in c(-4, -10, -100, -1e4, -1e6, -1e8, -1e10)) {
    expect_equal(mhn:::.mhn_log_normalizing_const(1, 1, gamma),
                 reference(gamma), tolerance = 1e-10,
                 info = sprintf("gamma=%g", gamma))
  }
  # No seam where the two evaluations meet.  A jump would show up as one
  # second difference out of line with its neighbours, so compare their spread
  # against their size rather than against an absolute figure -- the curve has
  # real curvature here, and second differences over a step of 0.01 are
  # legitimately of order 1e-6.  Measure that spread about the median: a jump
  # puts +J into one second difference and -J into the next, so a max/min ratio
  # turns negative and passes for every jump larger than the differences
  # themselves, which are about 3.8e-6 here.
  near <- vapply(seq(-4.04, -3.96, by = 0.01),
                 function(g) mhn:::.mhn_log_normalizing_const(1, 1, g), 0)
  second <- diff(diff(near))
  expect_true(all(second > 0))
  expect_lt(max(abs(second - stats::median(second))) / stats::median(second),
            0.05)
})

test_that("skewness and kurtosis survive a large tilt", {
  # Formed from raw moments, the third and fourth central moments cancel as
  # soon as sigma is small next to mu.  At (3, 1, 1000) the mean is about 500
  # and the standard deviation about 0.7, so mu^4 is 6e10 against a fourth
  # central moment below 1: even from a mean exact to rounding, that expansion
  # gives a skewness of exactly 0 there and an excess kurtosis of the wrong
  # sign at gamma = 200.
  for (gamma in c(20, 100, 200, 1000)) {
    k <- mhn_kurtosis(3, 1, gamma)
    s <- mhn_skewness(3, 1, gamma)
    expect_true(is.finite(k)); expect_true(is.finite(s))
    # Excess kurtosis is bounded below by -2 for every distribution, and at
    # alpha = 3 the MHN is right-skewed (as it is for alpha >= 1; for alpha < 1
    # a positive tilt can make it left-skewed, e.g. -0.30 at (0.1, 1, 5)), so
    # both have a known sign here.
    expect_gt(k, -2)
    expect_gt(s, 0)
    # A strong positive tilt drives the law to a normal, so both go to 0.
    expect_lt(abs(k), 0.01)
    expect_lt(abs(s), 0.01)
  }
  # Pearson's inequality, kurtosis >= skewness^2 - 2, holds for every
  # distribution.  Near a normal law the bound is close to -2, so it catches
  # only gross failures; the values themselves are checked against exact ones
  # in test-moments.R.  Skewness and kurtosis depend on beta only through
  # z = gamma / sqrt(beta) (Sun et al. 2023, Lemmas 1a and 2a), so each point
  # is also taken with beta scaled by 1e-160 and 1e140 and gamma by the square
  # root of that, which leaves z unchanged: a quadrature in the units of x
  # would leave the double range at those scales.  1e-10 relative covers the
  # rounding of z.
  for (p in list(c(0.5, 1, -5), c(2, 1, -30), c(3, 1, 0), c(5, 1, 2),
                 c(10, 0.1, 5), c(100, 1, -1))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3])
    k <- mhn_kurtosis(p[1], p[2], p[3]); s <- mhn_skewness(p[1], p[2], p[3])
    expect_gte(k, s^2 - 2 - 1e-6, label = info)
    for (e in c(-160, 140)) {
      b <- p[2] * 10^e
      g <- p[3] * 10^(e / 2)
      expect_equal(mhn_kurtosis(p[1], b, g), k, tolerance = 1e-10,
                   info = paste(info, "scaled by", 10^e))
      expect_equal(mhn_skewness(p[1], b, g), s, tolerance = 1e-10,
                   info = paste(info, "scaled by", 10^e))
    }
  }
})

test_that("qmhn reaches the quantile for small alpha", {
  # For alpha < 1 the distribution function behaves like x^alpha near zero, so
  # the quantile of a small p sits near p^(1/alpha) -- 2.3e-40 at alpha = 0.1,
  # p = 1e-4 -- far below any lower bracket set from the scale of the law.
  # The search has to keep halving that end until F there falls below p; a
  # bracket halved only 30 times from sqrt(eps) stops at 1.4e-17, and a solver
  # that returned that endpoint as the root would give pmhn 0.019 instead of
  # the 1e-4 asked for.
  for (alpha in c(0.5, 0.2, 0.1, 0.05)) {
    for (p in c(1e-2, 1e-4, 1e-8)) {
      q <- qmhn(p, alpha, 1, 1)
      expect_true(is.finite(q))
      expect_gt(q, 0)
      expect_close(pmhn(q, alpha, 1, 1), p, 1e-12,
                   info = sprintf("alpha=%g p=%g", alpha, p))
    }
  }
  # Far enough down the true quantile is below the smallest representable
  # double; the search must stop cleanly there rather than raising.
  expect_no_error(qmhn(1e-8, 0.01, 1, 1))
})

test_that("a large positive tilt, in gamma or in 1 / beta, reaches the normal limit", {
  # The Lemma 10(d) truncation length K of the series for Psi grows like z^2,
  # and at these tilts it passes the ceiling of 1e6 that psi_series puts on K.
  # Psi depends on the tilt only through z = gamma / sqrt(beta), and the
  # integrand is well behaved for either sign, so above z = 20 the quadrature
  # used for gamma < 0 serves instead.  dmhn, pmhn and the mean must then be
  # those of the normal limit.
  for (gamma in c(3e3, 1e4, 1e5)) {
    expect_no_error(mhn:::.mhn_log_normalizing_const(2.5, 1, gamma))
    mode_x <- (gamma + sqrt(gamma^2 + 8 * 1.5)) / 4
    d <- dmhn(mode_x, 2.5, 1, gamma)
    expect_true(is.finite(d)); expect_gt(d, 0)
    # A strong tilt drives the law to a normal of variance 1/2, whose density
    # at the mode is 1/sqrt(pi).
    expect_equal(d, 1 / sqrt(pi), tolerance = 1e-3)
    expect_equal(pmhn(mode_x, 2.5, 1, gamma), 0.5, tolerance = 1e-3)
    expect_equal(mhn_mean(2.5, 1, gamma), mode_x, tolerance = 1e-5)
  }
  # The same tilt, z = 2000, reached through a small beta rather than a large
  # gamma.  The normal limit then has variance 1/(2 beta), so the density at
  # the mode is sqrt(beta / pi).
  beta <- 1e-4
  mode_x <- (20 + sqrt(20^2 + 8 * beta * 1.5)) / (4 * beta)
  expect_equal(dmhn(mode_x, 2.5, beta, 20), sqrt(beta / pi), tolerance = 1e-3)
  expect_equal(pmhn(mode_x, 2.5, beta, 20), 0.5, tolerance = 1e-3)
  expect_equal(mhn_mean(2.5, beta, 20), mode_x, tolerance = 1e-5)
})

test_that("the mean stays exact at extreme tilt", {
  # E[X] is a ratio of two Fox-Wright values whose logarithms are both about
  # z^2/4 while their difference is only about log z, so subtracting them
  # costs an absolute error of order (z^2/4) times the unit roundoff, which is
  # a relative error in the ratio: 31% by gamma = 1e8.  For z = gamma /
  # sqrt(beta) above 12 the mean comes instead from the same quadrature the
  # central moments use, taken in the variable centred on the peak; every
  # gamma below is on that side.
  for (gamma in c(1e3, 1e5, 1e7, 1e8, 1e9, 1e10)) {
    got <- mhn_mean(2.5, 1, gamma)
    # A strong tilt drives the law to a normal centred on the mode, and the
    # mean approaches it to a relative O(1/gamma^2).
    mode_x <- (gamma + sqrt(gamma^2 + 8 * 1.5)) / 4
    expect_equal(got, mode_x, tolerance = 1e-8,
                 info = sprintf("gamma=%g", gamma))
  }
  # In the ordinary range the mean is the Lemma 2a ratio: at gamma = 0 it is
  # Gamma(1.75) / Gamma(1.25), and at gamma = -1 it is taken from the kernel
  # reference of helper-mhn-cdf.R.
  expect_equal(mhn_mean(2.5, 1, 0), gamma(1.75) / gamma(1.25),
               tolerance = 1e-12)
  expect_equal(mhn_mean(2.5, 1, -1), mhn_ref(2.5, 1, -1)$mean,
               tolerance = 1e-10)
})

test_that("the distribution function stays monotone past the support", {
  # A quadrature range trimmed below the anchor but not above it leaves the
  # upper panel running out to x, so once x is far past the support the rule
  # sees nothing there and F collapses to the value of the first panel alone:
  # at (2.5, 1, -30), 0.301 for every q beyond about 3600.
  for (p in list(c(2.5, 1, -30), c(1.5, 1, -50), c(5, 1, -10), c(0.7, 1, -20),
                 c(2, 0.1, -30), c(3, 1, 5))) {
    q <- c(0.01, 0.1, 1, 10, 100, 1e3, 3.6e3, 3.65e3, 1e4, 1e6, 1e8)
    v <- pmhn(q, p[1], p[2], p[3])
    expect_false(is.unsorted(v),
                 info = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
    expect_true(all(v >= 0 & v <= 1))
    expect_equal(v[length(v)], 1, tolerance = 1e-8)
  }
})

test_that("the moments match the Gamma limit that a strong negative tilt implies", {
  # As gamma -> -Inf the quadratic term stops mattering and the law tends to
  # Gamma(alpha, |gamma|), whose skewness is 2/sqrt(alpha) and whose excess
  # kurtosis is 6/alpha.  That is an independent check on the whole moment
  # path, and it covers both shapes the quadrature must handle: no interior
  # peak (alpha <= 1), where the raw-moment expansion gives -58485 for a
  # kurtosis of 6, and a peak squeezed towards the origin (alpha > 1), where a
  # window taken from the curvature alone would be far too wide.
  #
  # The limit is approached at a known rate.  Y = |gamma| X has density
  # proportional to y^(alpha-1) e^(-y) e^(-eps y^2), eps = beta / gamma^2, and
  # to first order in eps the factor e^(-eps y^2) lowers the n-th cumulant of
  # Y from alpha (n-1)! by eps alpha (alpha+1) (n+1)!.  So the skewness is
  # 2/sqrt(alpha) (1 - 3 (alpha+1) eps) and the excess kurtosis
  # 6/alpha (1 - 8 (alpha+1) eps), up to terms in eps^2 that stay below 2e-9
  # relative on this grid.  1e-8 leaves room for those and still resolves the
  # first-order term, 3.2e-5 of the kurtosis at alpha = 3, gamma = -1000, some
  # three thousand times over.
  for (alpha in c(0.2, 0.5, 1, 1.01, 1.2, 1.5, 3)) {
    for (gamma in c(-1000, -3000)) {
      eps <- 1 / gamma^2
      info <- sprintf("alpha=%g gamma=%g", alpha, gamma)
      expect_equal(mhn_kurtosis(alpha, 1, gamma),
                   6 / alpha * (1 - 8 * (alpha + 1) * eps),
                   tolerance = 1e-8, info = info)
      expect_equal(mhn_skewness(alpha, 1, gamma),
                   2 / sqrt(alpha) * (1 - 3 * (alpha + 1) * eps),
                   tolerance = 1e-8, info = info)
    }
  }
})

test_that("mhn_mean under a strong negative tilt is the Lemma 2a ratio", {
  # The cancellation in the ratio of two Fox-Wright values is one-sided: it
  # afflicts gamma > 0, where both logarithms grow like z^2/4.  For gamma < 0
  # they stay modest and the ratio is accurate well within the working
  # tolerance, so it is what the mean returns, and the two must agree to
  # rounding.
  for (gamma in c(-940, -960, -1000, -1e4)) {
    ratio <- exp(mhn:::.mhn_log_normalizing_const(2.2, 1, gamma) -
                 mhn:::.mhn_log_normalizing_const(1.2, 1, gamma))
    expect_equal(mhn_mean(1.2, 1, gamma), ratio, tolerance = 1e-12,
                 info = sprintf("gamma=%g", gamma))
  }
})

test_that("the alpha = 1 quantile survives a strongly negative tilt", {
  # Formed as 1 - Phi(-mu/sigma), the truncation mass rounds to exactly 0 once
  # -mu/sigma passes about 8.3 -- reached at gamma = -12 with beta = 1 -- and
  # qmhn would then return Inf for every p.  Taking it from pnorm's own upper
  # tail on the log scale, and assembling the answer as sigma * t rather than
  # mu + sigma * (t + z0), keeps both the mass and the subtraction in range.
  for (gamma in c(-2, -10, -30, -50, -200, -1000)) {
    for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
      q <- qmhn(p, 1, 1, gamma)
      expect_true(is.finite(q), info = sprintf("gamma=%g p=%g", gamma, p))
      expect_gt(q, 0)
      expect_equal(pmhn(q, 1, 1, gamma), p, tolerance = 1e-7,
                   info = sprintf("gamma=%g p=%g", gamma, p))
    }
  }
  # As gamma -> -Inf the law tends to Exponential(|gamma|) at beta = 1.
  expect_equal(qmhn(0.5, 1, 1, -1000), stats::qexp(0.5, 1000), tolerance = 1e-4)
})

test_that("the upper tail is evaluated, not derived from one minus the lower", {
  # Deriving P(X > q) as 1 - F loses it entirely once F rounds to 1, which
  # happens as soon as the survival probability falls below about 1e-16: at
  # alpha = 2.5, beta = 1, gamma = 1 that happens from q = 6.6, and at q = 8
  # 1 - F is exactly 0, and -Inf on the log scale, for a true value near
  # 5e-25.
  for (tri in list(c(2.5, 1, 1), c(1.5, 1, -2), c(0.7, 1, 3), c(3, 0.5, 0))) {
    a <- tri[1]; b <- tri[2]; g <- tri[3]
    for (q in c(3, 5, 8, 12, 20, 30)) {
      lq <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
      expect_true(is.finite(lq), info = sprintf("a=%g b=%g g=%g q=%g", a, b, g, q))
      expect_lt(lq, 0)
    }
    # The two tails still sum to one where both are representable.
    for (q in c(0.5, 1, 2)) {
      expect_equal(pmhn(q, a, b, g) + pmhn(q, a, b, g, lower.tail = FALSE), 1,
                   tolerance = 1e-12)
    }
  }

  # And qmhn inverts that tail rather than inverting 1 - p, which rounds to 1
  # below p = 1e-16 and would give Inf.
  for (tri in list(c(2.5, 1, 1), c(2, 1, 0), c(1, 1, 2), c(1, 1, -10))) {
    for (p in c(1e-8, 1e-16, 1e-20, 1e-30)) {
      q <- qmhn(p, tri[1], tri[2], tri[3], lower.tail = FALSE)
      expect_true(is.finite(q), info = sprintf("a=%g g=%g p=%g", tri[1], tri[3], p))
      expect_close(pmhn(q, tri[1], tri[2], tri[3], lower.tail = FALSE), p,
                   1e-11, info = sprintf("a=%g g=%g p=%g", tri[1], tri[3], p))
    }
    # log.p, where the point of the flag is that p need not be representable.
    q <- qmhn(-100, tri[1], tri[2], tri[3], lower.tail = FALSE, log.p = TRUE)
    expect_equal(pmhn(q, tri[1], tri[2], tri[3], lower.tail = FALSE, log.p = TRUE),
                 -100, tolerance = 1e-12)
  }

  # The gamma = 0 case must agree with base R to rounding, since both reduce to
  # the same qgamma call.
  for (p in c(1e-16, 1e-20, 1e-30, 1e-100)) {
    expect_equal(qmhn(p, 2, 1, 0, lower.tail = FALSE),
                 sqrt(stats::qgamma(p, 1, 1, lower.tail = FALSE)),
                 tolerance = 1e-12)
  }
})

# P(X > q) for a small shape, from the kernel on the log axis, x = exp(u):
# int exp(alpha u - beta e^{2u} + gamma e^u) du, whose part below U = -40 is
# exp(alpha U) / alpha to double precision, since e^{2U} and |gamma| e^U are
# below rounding there.  It calls no function of the package.
tail_reference <- function(alpha, beta, gamma, q) {
  kernel <- function(u) exp(alpha * u - beta * exp(2 * u) + gamma * exp(u))
  U <- -40
  tail_int <- stats::integrate(kernel, U, 40, rel.tol = 1e-13,
                               subdivisions = 10000L)$value
  upper <- stats::integrate(kernel, log(q), 40, rel.tol = 1e-13,
                            subdivisions = 10000L)$value
  upper / (exp(alpha * U) / alpha + tail_int)
}

test_that("the normalising constant keeps the mass below the smallest double", {
  # For gamma < 0 the normalising constant is a quadrature in log x that
  # stops at x of about 4.5e-9 / |z| and takes everything below that in
  # closed form.  A fraction 1/e of the mass of x^(alpha-1) on (0, 1) sits
  # below x = exp(-1/alpha).  At alpha = 1e-6 that is exp(-1e6), which no
  # abscissa in x can reach: a quadrature in x that stops at the smallest
  # doubles saturates at Psi = 1.3e3 against a true value of 2.0e6, and
  # would scale every density, probability and moment at small alpha by the
  # shortfall.
  #
  # Reference: substituting x = exp(u) turns the integral into
  # int exp(alpha u - beta e^{2u} + gamma e^u) du, whose left half is exactly
  # exp(alpha U)/alpha once e^{2U} and |gamma| e^U fall below rounding.
  psi_reference <- function(alpha, beta, gamma) {
    kernel <- function(u) exp(alpha * u - beta * exp(2 * u) + gamma * exp(u))
    U <- -40
    tail_int <- stats::integrate(kernel, U, 40, rel.tol = 1e-13,
                                 subdivisions = 10000L)$value
    log(2) + (alpha / 2) * log(beta) + log(exp(alpha * U) / alpha + tail_int)
  }

  for (alpha in c(1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1)) {
    for (gamma in c(-10, -1, -0.1)) {
      expect_equal(.mhn_log_normalizing_const(alpha, 1, gamma),
                   psi_reference(alpha, 1, gamma),
                   tolerance = 1e-8,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }

  # Tail probabilities are the user-visible consequence, and a missing head
  # would make them wrong by the same factor, about 1.5e3 at alpha = 1e-6.
  for (alpha in c(1e-6, 1e-4, 1e-3, 0.01)) {
    for (gamma in c(-10, -1)) {
      for (q in c(0.01, 0.5)) {
        expect_close(pmhn(q, alpha, 1, gamma, lower.tail = FALSE),
                     tail_reference(alpha, 1, gamma, q),
                     1e-11,
                     info = sprintf("alpha=%g gamma=%g q=%g", alpha, gamma, q))
      }
    }
  }
})

test_that("rmhn covers the whole support when the mode is tiny or alpha is small", {
  # rmhn() sends these (gamma < 0 with 1 < alpha < 1e3, at any n) to RTDR
  # region A, whose left contact search starts at half the mode.  The modes
  # here run from 1e-9 to 1e-6, so even an absolute floor of 1e-6 on that
  # start would put it at or to the right of the mode, and the search, which
  # tests the acceptance band before it clamps into range, would return a
  # contact point with the wrong sign of slope.
  for (par in list(c(4, 1, -1e7), c(4, 1, -3e6), c(3, 1, -1e7),
                   c(4, 1, -1e8), c(1.1, 1, -1e8))) {
    set.seed(20)
    draws <- rmhn(4000, par[1], par[2], par[3])
    expect_true(all(is.finite(draws)) && all(draws > 0))
    # At these tilts the law is Gamma(alpha, rate |gamma|) to a relative
    # O(alpha^2 / gamma^2), below 1e-11, so the sample mean is judged in
    # standard errors of that law: mean alpha / |gamma|, variance
    # alpha / gamma^2.
    z <- (mean(draws) - par[1] / abs(par[3])) /
      sqrt(par[1] / par[3]^2 / length(draws))
    expect_lt(abs(z), 5, label = sprintf("alpha=%g gamma=%g", par[1], par[3]))
  }

  # The right contact search runs on the log axis, where the density carries
  # exp(2y) and overflows past y = 354.  The start Gao & Wang (2025,
  # Section 2) recommend for their Eq. (8) iteration is
  # m + sqrt(-2 delta / L''(m)), which for a flat mode lands well beyond
  # that -- 447 above the mode at alpha = 1e-5.  The slope there is NaN, a
  # sign guard written as a plain
  # comparison lets it through because NaN fails every comparison, and the
  # plateau piece would be dropped without any error, leaving a sampler whose
  # support excluded the mode and the entire right tail.
  for (par in list(c(1e-4, 1, 0.01), c(3e-4, 1, 0.01), c(1e-5, 1, -1),
                   c(1e-6, 1, 0.001))) {
    set.seed(21)
    draws <- rmhn(100000, par[1], par[2], par[3])
    expect_true(all(is.finite(draws)) && all(draws >= 0))
    for (q in c(0.005, 0.05)) {
      expected <- tail_reference(par[1], par[2], par[3], q)
      # Binomial three-sigma band, floored so that a probability of a few
      # draws in 1e5 does not make the test flaky.
      band <- max(3 * sqrt(expected * (1 - expected) / 100000), 3e-5)
      expect_lt(abs(mean(draws > q) - expected), band)
    }
  }
})

test_that("the CDF quadrature keeps the mass below the smallest double", {
  # The same unreachable head as in the normalising constant: for alpha < 1
  # with a negative tilt the lower-tail quadrature runs from 0, and a fraction
  # 1/e of the mass of x^(alpha-1) on (0, 1) lies below x = exp(-1/alpha).  A
  # quadrature missing that head, divided by a normalising constant that has
  # it, gives F(1e-3) at (1e-6, 1, -1000) as 6.3e-4 against a true value of 1.
  #
  # Reference: with beta x^2 negligible against |gamma| x the law is
  # Gamma(alpha, rate = |gamma|), which base R evaluates exactly.
  for (alpha in c(1e-6, 1e-4, 1e-2)) {
    for (gamma in c(-1000, -100)) {
      for (q in c(1e-4, 1e-3)) {
        expect_equal(pmhn(q, alpha, 1, gamma),
                     stats::pgamma(q, alpha, rate = -gamma),
                     tolerance = 1e-5,
                     info = sprintf("alpha=%g gamma=%g q=%g", alpha, gamma, q))
      }
    }
  }

  # The upper tail on the same path.  Its range must be set from the slope of
  # the kernel, not from the Gaussian reach (894 at beta = 1e-4, where the
  # kernel dies within 0.08): set from the latter it leaves the whole mass in
  # a sliver the rule never samples, and the tail comes out as 0 for every q.
  for (q in c(1e-6, 1e-4, 1e-3, 1e-2)) {
    expect_close(pmhn(q, 0.05, 1e-4, -1000, lower.tail = FALSE),
                 stats::pgamma(q, 0.05, rate = 1000, lower.tail = FALSE),
                 1e-5, info = sprintf("q=%g", q))
  }
  # qmhn inverts that tail; inverting one that is 0 everywhere would return its
  # own initial bracket -- the same number for every p, with no warning.
  for (p in c(1e-2, 1e-3, 1e-4, 1e-6)) {
    expect_equal(qmhn(p, 0.05, 1e-4, -1000, lower.tail = FALSE),
                 stats::qgamma(p, 0.05, rate = 1000, lower.tail = FALSE),
                 tolerance = 1e-5, info = sprintf("p=%g", p))
  }
})

test_that("the variance stays accurate on both sides of the Lemma 2c cancellation", {
  # Var = alpha/(2 beta) + mu gamma/(2 beta) - mu^2.  Under a strong negative
  # tilt the first two terms are each about alpha/(2 beta) while the variance
  # is about alpha/gamma^2, so the closed form amplifies the rounding of its
  # own terms by about gamma^2/(2 beta).  Measured over the shapes below, at
  # gamma = -1e7 it is off by 2 to 34 percent, depending on alpha (17 percent
  # high at alpha = 1), and at -1e8 it is 8 to 80 times too large, or goes
  # negative and is clamped to 0.  The law then tends to Gamma(alpha, |gamma|),
  # whose variance is alpha/gamma^2.
  for (alpha in c(1e-4, 0.01, 0.3, 1, 3, 10)) {
    for (gamma in c(-1e4, -1e6, -1e8)) {
      expect_close(mhn_var(alpha, 1, gamma), alpha / gamma^2,
                   1e-3,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
  # Under a strong positive tilt the last two terms are each about
  # gamma^2/(4 beta^2) while the variance is about 1/(2 beta), and the form
  # multiplies the relative error of the mean by about gamma^2/(2 beta).  A
  # mean formed as the ratio of two Fox-Wright values of size exp(z^2/4),
  # z = gamma / sqrt(beta), carries a relative error of up to about eps z^2/4
  # (half an ulp in each logarithm), so the variance the form would give from
  # it is off by up to about eps z^4/8; measured, 2e-9 at gamma = 100, 1e-7 at
  # 300 and 1e-5 at 900.
  # Reference: the second central moment by stats::integrate in the offset t
  # from the mode m, with the log kernel written as
  # (alpha-1) log1p(t/m) + t (gamma - 2m) - t^2, in which nothing cancels, over
  # +-12 standard-deviation-sized units (the mass beyond is below 1e-60).  The
  # quadrature is asked for 1e-13; 1e-8 is the working tolerance of the
  # moment helpers.
  central_var <- function(alpha, gamma) {
    m <- mhn_mode(alpha, 1, gamma)
    k <- function(t) exp((alpha - 1) * log1p(t / m) + t * (gamma - 2 * m) - t^2)
    br <- seq(-12, 12, by = 3)
    panel <- function(f) sum(vapply(seq_len(length(br) - 1), function(i)
      stats::integrate(f, br[i], br[i + 1], rel.tol = 1e-13)$value, 0))
    z0 <- panel(k)
    d <- panel(function(t) t * k(t)) / z0
    panel(function(t) (t - d)^2 * k(t)) / z0
  }
  for (alpha in c(1e-4, 0.01, 0.3, 1, 3, 10)) {
    for (gamma in c(100, 300, 900)) {
      expect_equal(mhn_var(alpha, 1, gamma), central_var(alpha, gamma),
                   tolerance = 1e-8,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
  # Where the form loses little to cancellation -- |gamma| <= 10 at beta = 1,
  # where its terms are at most about 80 times the variance -- it is accurate,
  # and mhn_var must agree with it to rounding.
  for (alpha in c(1e-4, 0.01, 0.3, 1, 3, 10)) {
    for (gamma in c(-10, -1, 0.5, 10)) {
      mu <- mhn_mean(alpha, 1, gamma)
      expect_equal(mhn_var(alpha, 1, gamma),
                   alpha / 2 + mu * (gamma / 2 - mu),
                   tolerance = 1e-12,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
})

test_that("pmhn agrees element-wise and vectorised under every tail flag", {
  # Each special case has one per-element implementation, called both by the
  # all-scalar path and by the loop that recycles vector parameters, so the
  # two must return the same double under every tail flag.  The grid is the
  # truncated normal at strong tilts of either sign.  At q = 3.0184056,
  # gamma = 18 the lower tail is 1.35e-17: derived as 1 minus the upper tail it
  # would be exactly 0.
  for (gamma in c(18, 30, 10, -10)) {
    for (q in c(0.5, 1, 2, 3.0184056)) {
      for (lt in c(TRUE, FALSE)) {
        for (lp in c(TRUE, FALSE)) {
          scalar <- pmhn(q, 1, 1, gamma, lower.tail = lt, log.p = lp)
          vector <- pmhn(q, alpha = c(1, 1), beta = 1, gamma = gamma,
                         lower.tail = lt, log.p = lp)[1]
          expect_identical(vector, scalar,
                           info = sprintf("g=%g q=%g lower.tail=%s log.p=%s",
                                          gamma, q, lt, lp))
        }
      }
    }
  }
  expect_gt(pmhn(3.0184056, 1, 1, 18), 1e-17)
})

test_that("the RTDR region boundary at alpha exactly 1 is sampled, not rejected", {
  # Region A needs an interior mode in x, and at alpha = 1 with a non-positive
  # tilt there is none: the mode formula returns the origin.  An envelope
  # built around it would degenerate to a plateau spanning [0, Inf] and
  # reject every proposal, so alpha = 1 belongs to region BC, whose log-axis
  # mode is positive for every alpha > 0.  The public API is shielded by the
  # truncated-normal interception, so this is reached through the kernel.
  for (alpha in c(1 - 1.49e-8, 1, 1 + 1.49e-8)) {
    for (gamma in c(-2, 0, 2)) {
      set.seed(4)
      draws <- mhn:::.rmhn_rtdr_cpp(4000, alpha, 1, gamma)
      expect_false(anyNA(draws),
                   info = sprintf("alpha-1=%g gamma=%g", alpha - 1, gamma))
      expect_true(all(draws > 0))
      expect_equal(mean(draws), mhn_mean(alpha, 1, gamma), tolerance = 0.05,
                   info = sprintf("alpha-1=%g gamma=%g", alpha - 1, gamma))
    }
  }
  # alpha = 1/2, the other branch point of classify_region.
  for (gamma in c(-2, 0.5, 5)) {
    set.seed(5)
    draws <- mhn:::.rmhn_rtdr_cpp(4000, 0.5, 1, gamma)
    expect_false(anyNA(draws))
    expect_equal(mean(draws), mhn_mean(0.5, 1, gamma), tolerance = 0.06,
                 info = sprintf("gamma=%g", gamma))
  }
})

test_that("method = \"sun\" accepts what the dispatcher answers in closed form", {
  # The pre-scan that refuses method = "sun" for alpha < 1 with gamma > 0 must
  # test the same scale-free predicates the dispatcher routes on.  Tested on
  # the bare floats it would reject parameters the sun path handles perfectly
  # well -- alpha within sqrt(eps) of 1, or a negligible tilt.
  for (par in list(c(1 - 1e-9, 1, 2), c(1, 1, 2), c(0.5, 1, 1e-12),
                   c(0.5, 1e18, 1))) {
    set.seed(6)
    expect_silent(draws <- rmhn(5, par[1], par[2], par[3], method = "sun"))
    expect_true(all(is.finite(draws)) && all(draws > 0),
                info = sprintf("alpha=%g beta=%g gamma=%g",
                               par[1], par[2], par[3]))
  }
  # The genuinely unsupported combination must still be refused.
  expect_error(rmhn(5, 0.5, 1, 2, method = "sun"), "alpha<1 and gamma>0")
})

test_that("the log density does not cancel against the normalising constant", {
  # log Psi and the kernel exponent -beta x^2 + gamma x are each about
  # gamma^2/(4 beta) near the mode, while the log density they combine into is
  # of order one.  Assembled directly, the result would carry an absolute
  # error of about eps gamma^2/(4 beta): 3.7e-8 at |z| = 3e4, and 4.6 percent
  # at 3e7.  For z above 1e3 both sides have that term removed analytically.
  #
  # Reference: completing the square gives
  #   log f(x) = log 2 + (alpha/2) log beta - P + (alpha-1) log x
  #              - beta (x - gamma/(2 beta))^2,
  # with P = log Psi - z^2/4 = log(2 int u^(alpha-1) exp(-(u - z/2)^2) du),
  # in which nothing of size z^2 appears.
  shifted_log_density <- function(alpha, beta, gamma, x) {
    z <- gamma / sqrt(beta)
    P <- log(2 * stats::integrate(
      function(u) exp((alpha - 1) * log(u) - (u - z / 2)^2),
      max(0, z / 2 - 40), z / 2 + 40,
      rel.tol = 1e-12, subdivisions = 4000L)$value)
    log(2) + (alpha / 2) * log(beta) - P + (alpha - 1) * log(x) -
      beta * (x - gamma / (2 * beta))^2
  }

  for (par in list(c(3, 1, 3e3), c(3, 1, 3e4), c(1.5, 1, 5e3), c(10, 1, 1e4),
                   c(3, 1e-6, 3e3), c(3, 1e-6, 3e4))) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    mode_x <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta)
    ref <- shifted_log_density(alpha, beta, gamma, mode_x)
    expect_equal(dmhn(mode_x, alpha, beta, gamma, log = TRUE), ref,
                 tolerance = 1e-10,
                 info = sprintf("alpha=%g beta=%g gamma=%g", alpha, beta, gamma))
  }

  # The scalar fast path and the recycling loop share one implementation of
  # the general log density, so they must return the same doubles.
  for (par in list(c(3, 1, 3e4), c(3, 1e-6, 3e3), c(2.5, 1, 5))) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    xs <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta) *
      c(0.5, 0.9, 1, 1.1, 2)
    expect_identical(dmhn(xs, alpha = c(alpha, alpha), beta = beta,
                          gamma = gamma, log = TRUE),
                     dmhn(xs, alpha, beta, gamma, log = TRUE))
  }
})

test_that("the CDF cancellation guard fires before the tolerance is spent", {
  # For gamma < 0 the CDF series alternates, and its two sums, of the even and
  # of the odd terms, nearly cancel.  The guard compares their relative
  # difference with the rounding carried by the terms, weighted by the size of
  # each term (see mhn_cdf_series.cpp), and returns the NaN sentinel when the
  # result could miss the working tolerance.  A guard that charged only one
  # unit of roundoff to each sum would admit results many times worse than
  # that tolerance.
  #
  # The assertion is on the guard's contract rather than on a reference value:
  # wherever the series is allowed through, it must agree with the integration
  # path, which the cdf_series_accuracy cross-check in test-pmhn.R validates
  # against an arbitrary-precision reference across the whole grid.
  for (alpha in c(0.3, 1.5, 3, 10)) {
    for (gamma in c(-18, -20, -22, -25, -50, -100)) {
      for (qf in c(0.5, 1, 2)) {
        q <- qf * mhn_mean(alpha, 1, gamma)
        series <- mhn:::.pmhn_force(q, alpha, 1, gamma, "series")
        if (is.nan(series)) next          # guard fired; nothing to check
        quad <- mhn:::.pmhn_force(q, alpha, 1, gamma, "integrate")
        expect_equal(series, quad, tolerance = 1e-8,
                     info = sprintf("alpha=%g gamma=%g q=%g", alpha, gamma, q))
      }
    }
  }
})

test_that("the two tails of pmhn sum to one", {
  # The cheapest check there is on the survival function.
  #
  # An upper-tail range anchored at max(peak, q) fails for q at or below the
  # mode: the anchor is then the mode, where the slope is zero by definition,
  # so a reach derived from the slope collapses to the Gaussian width -- 894
  # at beta = 1e-4, where the kernel dies within 0.08 -- and the panel from the
  # mode outward sees nothing.  pmhn(q, 3, 1e-4, -1000, lower.tail = FALSE)
  # would then report only the mass between q and the mode, the two tails
  # would sum to F(mode) = 0.323 for every q below it, and qmhn would return
  # the denormal floor for every p.  A range spanning sixty decades, at
  # alpha < 1 with a tiny q, fails the same way.  The package takes the
  # integral on the log axis, with panels that double in width away from the
  # peak, which handles both.
  for (alpha in c(0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 3, 10)) {
    for (beta in c(1e-4, 0.01, 1, 100)) {
      for (gamma in c(-1000, -100, -10, -1, 0, 1, 10, 100)) {
        mu <- tryCatch(mhn_mean(alpha, beta, gamma), error = function(e) NA_real_)
        if (!is.finite(mu) || mu <= 0) next
        for (f in c(1e-60, 1e-6, 0.1, 0.5, 1, 2, 10)) {
          q <- f * mu
          lo <- pmhn(q, alpha, beta, gamma)
          up <- pmhn(q, alpha, beta, gamma, lower.tail = FALSE)
          expect_lt(abs(lo + up - 1), 1e-8,
                    label = sprintf("alpha=%g beta=%g gamma=%g q=%.4g",
                                    alpha, beta, gamma, q))
        }
      }
    }
  }
})

test_that("the survival function matches its Gamma limit where beta is negligible", {
  # With beta x^2 small against |gamma| x the law is Gamma(alpha, rate=|gamma|),
  # which base R evaluates in both tails without cancelling.  This pins the
  # upper tail against something outside the package.
  #
  # The approximation is only as good as beta * x / |gamma| is small, and x runs
  # out to the 0.001 upper quantile, so the triples below are chosen to keep
  # that ratio under 1e-7 there.  (0.05, 0.01, -1) does not: its ratio reaches
  # 0.027 and the Gamma law is then 3 percent away from the MHN one, which says
  # nothing about the package.
  for (par in list(c(3, 1e-4, -1000), c(0.05, 1e-8, -1), c(1, 1e-4, -1000),
                   c(0.5, 1e-6, -100), c(10, 1e-4, -1000))) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    # Guard the guard: the reference must actually be the right law here.
    expect_lt(beta * stats::qgamma(1e-3, alpha, rate = -gamma,
                                   lower.tail = FALSE) / -gamma, 1e-7)
    for (p in c(0.99, 0.5, 0.1, 1e-3)) {
      q <- stats::qgamma(p, alpha, rate = -gamma, lower.tail = FALSE)
      expect_close(pmhn(q, alpha, beta, gamma, lower.tail = FALSE), p,
                   1e-5,
                   info = sprintf("alpha=%g beta=%g gamma=%g p=%g",
                                  alpha, beta, gamma, p))
      expect_close(qmhn(p, alpha, beta, gamma, lower.tail = FALSE), q,
                   1e-5,
                   info = sprintf("alpha=%g beta=%g gamma=%g p=%g",
                                  alpha, beta, gamma, p))
    }
  }
})

test_that("the density vanishes at infinity on every path", {
  # Every base R density returns 0 there.  Evaluated as written, the general
  # expression gives NaN, since (alpha-1) log(x) and gamma x are +Inf while
  # -beta x^2 is -Inf, and the sqrt-Gamma form gives NaN for the same reason.
  for (par in list(c(2, 1, 1), c(2, 1, 0), c(1, 1, 2), c(0.3, 1, -1),
                   c(0.3, 1, 5), c(10, 4, -20))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", par[1], par[2], par[3])
    expect_identical(dmhn(Inf, par[1], par[2], par[3]), 0, info = info)
    expect_identical(dmhn(Inf, par[1], par[2], par[3], log = TRUE), -Inf,
                     info = info)
    expect_identical(dmhn(-Inf, par[1], par[2], par[3]), 0, info = info)
  }
  # And in a vector, alongside the other boundary values.
  expect_equal(dmhn(c(1, Inf, -Inf), alpha = c(2, 2, 2), beta = 1, gamma = 1),
               c(dmhn(1, 2, 1, 1), 0, 0))
})

test_that("a steep tilt samples its Gaussian limit, and reports nothing", {
  # At alpha = 0.7 with beta = 1e-8 and gamma = 1e8 the standardised tilt
  # gamma / sqrt(beta) is 1e12, and region BC's contact drop, read off as a
  # difference of two values of log g whose own size is
  # gamma_norm^2 / 4 = 2.5e23, would be a small multiple of ulp(2.5e23) =
  # 3.4e7 against a target drop of log 4.  The envelope would then be decided
  # by the last bit of exp(), and differently on different platforms: a
  # NaN-bounded envelope and 50 NaN draws with a retry-budget warning on one,
  # an error from the p_l >= p_r guard on another.  The drop is measured from
  # the mode instead.
  #
  # There is nothing exotic about the density here: at this tilt it is
  # Normal(mode, 1 / (2 beta)) to within a relative 1e-24 in the variance, so
  # the draws are checked against that.
  mu <- mhn_mode(0.7, 1e-8, 1e8)      # 5e15
  sigma <- 1 / sqrt(2 * 1e-8)         # 7071.07
  set.seed(2)
  draws <- expect_silent(rmhn(50, 0.7, 1e-8, 1e8))
  expect_length(draws, 50L)
  expect_true(all(is.finite(draws)))
  expect_lt(abs(mean(draws) - mu), 5 * sigma / sqrt(50))
  expect_gt(sd(draws) / sigma, 0.5)
  expect_lt(sd(draws) / sigma, 2)

  # A working sampler stays quiet.  No draw here exhausts its retry budget, so
  # any warning would be spurious.
  set.seed(2)
  expect_silent(rmhn(50, 2, 1, 1))
  set.seed(2)
  expect_silent(rmhn(50, 0.3, 1, 100))
})

test_that("the region BC envelope stays well posed across the tilt range", {
  skip_on_cran()

  # The same invariant as a sweep.  Read off the uncentred log g, the envelope
  # would be decided by rounding noise from gamma / sqrt(beta) ~ 1e8 upward,
  # so the standardised tilt is swept from 1e6 to 1e13, requiring at each
  # point an ordered finite envelope with finite piece areas and draws
  # consistent with the density's Gaussian limit.  None of these depends on
  # libm rounding.  The sweep stops at 1e13 to bound its run time; the sampler
  # itself refuses only once the law's standard deviation is no larger than
  # the spacing of the doubles at its mode, in the caller's coordinate -- at
  # beta = 1 from gamma / sqrt(beta) = 2^53 = 9.0e15, where the mode passes
  # 2^52 and the doubles there lie a unit apart, farther than the standard
  # deviation 1/sqrt(2), and at beta = 1e-8 from about 7.4e15 (see the refusal
  # test below).
  for (log_tilt in seq(6, 13, by = 0.5)) {
    for (beta in c(1e-8, 1)) {
      gamma <- 10^log_tilt * sqrt(beta)
      info <- sprintf("alpha = 0.7, beta = %g, gamma = %g", beta, gamma)
      env <- mhn:::.dump_rtdr_envelope_cpp(0.7, beta, gamma)
      expect_identical(env$region, 1L, info = info)
      expect_true(is.finite(env$p_l) && is.finite(env$p_r), info = info)
      expect_lt(env$p_l, env$p_r)
      expect_true(all(is.finite(env$piece_log_area)), info = info)

      set.seed(4)
      draws <- rmhn(2000, 0.7, beta, gamma)
      expect_true(all(is.finite(draws)), info = info)
      mu <- mhn_mode(0.7, beta, gamma)
      sigma <- 1 / sqrt(2 * beta)
      expect_lt(abs(mean(draws) - mu) / (sigma / sqrt(2000)), 5)
      expect_lt(abs(sd(draws) / sigma - 1), 0.15)
    }
  }
})

test_that("the region A envelope stays well posed across the tilt range", {
  skip_on_cran()

  # Region BC's invariant, applied to region A.  Region A works on f(x) in x
  # rather than on the log axis, and log f(x) = (alpha-1) log x - x^2 +
  # gamma_norm x is gamma_norm^2/4 - (x - gamma_norm/2)^2 up to the power term,
  # so every value of it carries gamma_norm^2/4 -- 2.5e19 at gamma_norm = 1e10,
  # one ulp of which is 4096 -- while the contact drop, the two intersections,
  # the piece areas and the accept test are all O(1) differences of two such
  # values.  Region A therefore measures its ordinate from the mode.
  #
  # Two of the assertions below, on the width of the plateau and on the
  # acceptance rate, catch an unconverged contact search: started a fixed
  # fraction of the mode away from it, the search would need log2(m/2) Newton
  # steps against a budget of 30, so past gamma_norm ~ 2e9 it would return a
  # contact whose drop is hundreds instead of one.  The draws would stay
  # correctly distributed -- any tangent of a concave log-density dominates --
  # so neither the KS test nor the moments see it.  Only the width of the
  # plateau and the acceptance rate do, and the acceptance rate is the quantity
  # Gao & Wang (2025, Theorem 3.1) bound below by 1/e.
  #
  # The default method sends gamma > 0 with alpha > 1 to Sun's Algorithm 1 at
  # tilts above 2 sqrt(alpha), which covers this sweep, so region A has to be
  # asked for by name here.  The sweep stops at 1e13 to bound its run time;
  # region A returns the exact draw rounded once to double precision and
  # refuses only where the law is narrower than the spacing of the doubles at
  # its mode (at beta = 1 from gamma / sqrt(beta) = 2^53 = 9.0e15; see the
  # refusal test below).
  for (log_tilt in seq(6, 13, by = 0.5)) {
    for (beta in c(1e-8, 1)) {
      for (alpha in c(1.5, 2, 10)) {
        gamma <- 10^log_tilt * sqrt(beta)
        info <- sprintf("alpha = %g, beta = %g, gamma = %g", alpha, beta, gamma)
        env <- mhn:::.dump_rtdr_envelope_cpp(alpha, beta, gamma)
        expect_identical(env$region, 0L, info = info)
        expect_true(is.finite(env$p_l) && is.finite(env$p_r), info = info)
        expect_gt(env$p_l, 0)
        expect_lt(env$p_l, env$p_r)
        expect_lt(env$t_l, env$mode)
        expect_gt(env$t_r, env$mode)
        expect_true(all(is.finite(env$piece_log_area)), info = info)

        # The envelope is built in the beta = 1 coordinate, where the limiting
        # standard deviation is 1/sqrt(2).  An unconverged contact search shows
        # up here and nowhere else: from the fixed-fraction start this ratio
        # would be 1646 at gamma_norm = 1e13.
        expect_lt((env$p_r - env$p_l) / 2 * sqrt(2), 20)

        set.seed(4)
        draws <- rmhn(2000, alpha, beta, gamma, method = "rtdr")
        expect_true(all(is.finite(draws)), info = info)
        mu <- mhn_mode(alpha, beta, gamma)
        sigma <- 1 / sqrt(2 * beta)
        expect_lt(abs(mean(draws) - mu) / (sigma / sqrt(2000)), 5)
        expect_lt(abs(sd(draws) / sigma - 1), 0.15)
        expect_lt(suppressWarnings(
          stats::ks.test((draws - mu) / sigma, "pnorm")$statistic), 0.05)

        # Acceptance stays above the 1/e floor of Gao & Wang (2025,
        # Theorem 3.1).  This is the invariant that the wrong contact point
        # violates: from the fixed-fraction start it would fall to 0.0076 at
        # gamma_norm = 1e12.
        set.seed(4)
        d2 <- mhn:::.rmhn_rtdr_cpp(2000L, alpha, beta, gamma)
        expect_gt(2000 / (2000 + attr(d2, "rtdr_retries")), 1 / exp(1))
      }
    }
  }
})

test_that("the region D envelope stays well posed across the tilt range", {
  skip_on_cran()

  # The same invariant for region D, which is on the *default* path: rmhn's
  # auto dispatch sends every alpha < 1 with gamma > 0 to RTDR, and
  # classify_region sends alpha < 1/2 above the gamma_d threshold to region D.
  # So every failure this block asserts against would be reachable from a
  # plain rmhn(n, 0.3, 1, 1e8).
  #
  # Region D is the inflection-point envelope on log g(y) for the log-axis
  # density g(y) = exp(alpha y - e^{2y} + gamma_norm e^y), with log-linear
  # pieces throughout: a secant chain over the convex half, a tangent over the
  # concave half, a plateau, and a right tail.
  # It needs the region BC mode-centred ordinate and three safeguards of its
  # own: a contact search started where it converges, a dual point and rho
  # formed without differences of gamma_norm-sized values, and the count
  # ceil(rho) kept in double, since narrowing it into an int is undefined
  # behaviour past gamma_norm ~ 1.8e9 and saturates one way on arm64 and the
  # other on x86-64.  Without them, at gamma_norm = 1e10, the plateau's right
  # edge would lie to the left of its left edge, the piece would be dropped
  # without any error, and the mode would lie outside the sampler's support.
  #
  # The sweep stops at 1e12 to bound its run time; the sampler refuses only
  # where the law is narrower than the spacing of the doubles at its mode (at
  # beta = 1 from gamma / sqrt(beta) = 2^53 = 9.0e15; see the refusal test
  # below).
  for (log_tilt in seq(6, 12, by = 0.5)) {
    for (beta in c(1e-8, 1)) {
      for (alpha in c(0.05, 0.3, 0.49)) {
        gamma <- 10^log_tilt * sqrt(beta)
        info <- sprintf("alpha = %g, beta = %g, gamma = %g", alpha, beta, gamma)
        env <- mhn:::.dump_rtdr_envelope_cpp(alpha, beta, gamma)
        expect_identical(env$region, 2L, info = info)
        expect_gte(env$K_eff, 1L)
        expect_true(is.finite(env$p_l) && is.finite(env$p_r), info = info)
        expect_lt(env$p_l, env$p_r)
        expect_true(all(is.finite(env$y_break)), info = info)
        expect_true(all(is.finite(env$log_dens_break)), info = info)
        expect_true(all(is.finite(env$alpha_k)), info = info)
        expect_true(all(is.finite(env$piece_log_area)), info = info)

        # rho is v (gamma_norm - v) with v <= gamma_norm/2, so it is positive
        # and below gamma_norm^2/4; formed as a difference of gamma_norm-sized
        # values it would be 41 times its true value at gamma_norm = 1e10.
        gn <- gamma / sqrt(beta)
        expect_gt(env$rho, 0)
        expect_lt(env$rho, gn * gn / 4)

        # Exactly one plateau piece, spanning the mode.  Dropping it silently
        # would leave the mode outside the support.
        types <- vapply(env$pieces, function(p) p$type, integer(1))
        expect_identical(sum(types == 0L), 1L, info = info)
        expect_lt(env$p_l, env$mode)
        expect_gt(env$p_r, env$mode)

        # Plateau half-width against the log-axis standard deviation,
        # 1/sqrt(-L''(m_g)) with -L''(m_g) = u_mode (4 u_mode - gamma_norm).
        u_m <- exp(env$mode)
        sd_y <- 1 / sqrt(u_m * (4 * u_m - gn))
        expect_lt((env$p_r - env$p_l) / 2 / sd_y, 20)

        set.seed(4)
        draws <- rmhn(2000, alpha, beta, gamma)
        expect_true(all(is.finite(draws)), info = info)
        mu <- mhn_mode(alpha, beta, gamma)
        sigma <- 1 / sqrt(2 * beta)
        expect_lt(abs(mean(draws) - mu) / (sigma / sqrt(2000)), 5)
        expect_lt(abs(sd(draws) / sigma - 1), 0.15)
        expect_lt(suppressWarnings(
          stats::ks.test((draws - mu) / sigma, "pnorm")$statistic), 0.05)

        set.seed(4)
        d2 <- mhn:::.rmhn_rtdr_cpp(2000L, alpha, beta, gamma)
        expect_gt(2000 / (2000 + attr(d2, "rtdr_retries")), 1 / exp(1))
      }
    }
  }
})

test_that("a tilt past the representable range is refused in the user's terms", {
  # Every region refuses to sample once the law's standard deviation is no
  # larger than the spacing of the doubles at its mode, in the caller's
  # coordinate: at beta = 1 from a tilt of 2^53 = 9.0e15, and on the shape
  # axis at the latest from alpha = 2^105 = 4.1e31.  Past that point draws
  # could take only a few distinct values, and finite draws of the wrong
  # spread would be worse than an error because nothing would say so.  The
  # message must not tell the user they have found a defect: they have found
  # the edge of double precision, and the normal limit is the answer there.
  expect_error(rmhn(10, 2, 1, 1e16, method = "rtdr"),
               "narrower than the spacing of the double")
  expect_error(rmhn(10, 0.3, 1, 1e16), "narrower than the spacing of the double")
  expect_error(rmhn(10, 1e33, 1, 1, method = "rtdr"),
               "narrower than the spacing of the double")
  # Inside that limit it samples, with no NaN (test-rmhn-rtdr-regions.R checks
  # the law of the draws up to 8e15).
  set.seed(1)
  expect_true(all(is.finite(rmhn(200, 2, 1, 1e14, method = "rtdr"))))
  set.seed(1)
  expect_true(all(is.finite(rmhn(200, 0.3, 1, 1e12))))
})

test_that("a non-finite shape or rate is taken as outside the parameter space", {
  # alpha > 0 is true of Inf, so a test of positivity alone would pass an
  # infinite shape on to the evaluators, which have no meaning for it, and the
  # user would see whatever failed first deep inside them -- an exception from
  # a container or a quadrature, whose message may be as bare as the word
  # "vector", or a NaN without a warning.  The d/p/q functions give NaN with
  # a warning, as for any parameter outside the space, and the moment
  # functions stop with an error that names the argument.
  for (f in list(dmhn, pmhn)) {
    expect_warning(out <- f(1, Inf, 1, 1), "NaNs produced")
    expect_identical(out, NaN)
    expect_warning(out <- f(1, 1, Inf, 1), "NaNs produced")
    expect_identical(out, NaN)
    expect_identical(f(1, NaN, 1, 1), NaN)
  }
  expect_warning(out <- qmhn(0.5, Inf, 1, 1), "NaNs produced")
  expect_identical(out, NaN)
  for (f in list(mhn_mean, mhn_var, mhn_skewness, mhn_kurtosis, mhn_mode)) {
    expect_error(f(Inf, 1, 1), "alpha must be positive and finite")
    expect_error(f(1, Inf, 1), "beta must be positive and finite")
  }
  # The vector path too, where the offending element is not the first.
  expect_warning(out <- dmhn(1, alpha = c(2, Inf), beta = 1, gamma = 1),
                 "NaNs produced")
  expect_identical(out, c(dmhn(1, 2, 1, 1), NaN))

  # rmhn draws NaN there, with a warning, as rnorm() does.
  set.seed(1)
  expect_warning(x <- rmhn(4, alpha = c(2, Inf, 2, NaN), beta = 1, gamma = 1),
                 "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("the distribution function approaches its normal limit at a strong tilt", {
  # log Psi and the log of the integral over [0, x] are each about
  # gamma^2/(4 beta) -- 2.5e13 at gamma = 1e7 -- while the log probability they
  # combine into is of order one.  Assembled directly the answer carries the
  # rounding of both, and how that rounding falls depends on the platform's
  # quadrature and libm: it can cost four digits, or enough for F at the mode
  # to be clamped to 1 rather than 0.5.  For z above 1e3 both sides have the
  # term removed analytically, and the integrand is written in the variable
  # centred on its own anchor.
  #
  # As the tilt grows the law tends to a normal of variance 1/(2 beta) centred
  # on the mode, so the CDF must approach Phi at every multiple of that scale.
  for (alpha in c(1.5, 2.5, 10)) {
    for (beta in c(1, 100)) {
      for (gamma in c(1e4, 1e5, 1e7)) {
        mode_x <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta)
        sdev <- 1 / sqrt(2 * beta)
        for (k in c(-2, -1, -0.5, 0, 0.5, 1, 2)) {
          expect_equal(pmhn(mode_x + k * sdev, alpha, beta, gamma),
                       stats::pnorm(k), tolerance = 1e-4,
                       info = sprintf("alpha=%g beta=%g gamma=%g k=%g",
                                      alpha, beta, gamma, k))
        }
      }
    }
  }
})
