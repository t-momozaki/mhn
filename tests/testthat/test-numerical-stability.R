# Regression tests for closed forms that lose precision in double arithmetic.
#
# Each block below pins a parameter region where a textbook algebraic form
# silently returns a wrong answer -- zero, NaN, or a value outside the domain
# the source theorem requires -- while the mathematically identical form used
# by the package stays exact.  The grids elsewhere in the suite stop well short
# of these regions, which is why they went unnoticed.

test_that("mhn_mode stays exact for strongly negative gamma", {
  # Sun et al. (2023) Lemma 3b as printed evaluates
  #   (gamma + sqrt(gamma^2 + 8 beta (alpha-1))) / (4 beta),
  # whose numerator cancels for gamma < 0.  Gao & Wang (2025, Eq. 9) also give
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
  # overflows a double for z > 53.3, which drove the normalising constant to
  # +Inf and the whole density to zero.
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
  # No collapse across the former overflow boundary.  One unit above the mode
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
  # Sun et al. (2023) Theorem 1a requires delta in (0, beta).  Evaluated
  # literally the optimum is a difference of two quantities of size gamma^2
  # and leaves the interval for large tilt; a negative value made the proposal
  # a Gamma with negative scale, so every draw came back NaN.
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
  # The (alpha, beta, gamma) triples at which the old delta_opt turned
  # negative.  Each returned an all-NaN vector, with no error or warning.
  cases <- list(c(1.5, 0.01, 3162.3), c(2, 0.01, 5623.4),
                c(1.5, 1, 56234.1),   c(10, 0.01, 31622.8))
  for (p in cases) {
    set.seed(42)
    x <- rmhn(200, alpha = p[1], beta = p[2], gamma = p[3])
    expect_false(anyNA(x))
    expect_true(all(is.finite(x)))
    expect_true(all(x > 0))
    # At this tilt the law concentrates tightly at its mode, so the exact mean
    # sits there too.  Judge the sample mean on the scale of its own standard
    # error: 5% of a mean of 1e5 is more than 1e4 standard errors wide, a band
    # no sampler could fail.
    mode_x <- (p[3] + sqrt(p[3]^2 + 8 * p[2] * (p[1] - 1))) / (4 * p[2])
    expect_equal(mhn_mean(p[1], p[2], p[3]), mode_x, tolerance = 1e-6)
    z <- (mean(x) - mhn_mean(p[1], p[2], p[3])) /
         sqrt(mhn_var(p[1], p[2], p[3]) / length(x))
    expect_lt(abs(z), 5)
  }
})

test_that("pmhn stays bounded in cost and range for strongly negative gamma", {
  # The Lemma 10(d) truncation length grows like z^2, so the CDF series asked
  # for ~5e7 terms at gamma = -1e4.  Past a ceiling the series now defers to
  # quadrature, which answers in a bounded number of evaluations.
  for (gamma in c(-800, -2000, -10000, -1e5)) {
    elapsed <- system.time(v <- pmhn(0.5, alpha = 2.5, beta = 1, gamma = gamma))[["elapsed"]]
    expect_true(is.finite(v))
    expect_gte(v, 0)
    expect_lte(v, 1)
    expect_lt(elapsed, 5)
  }
  # For gamma < 0 the alternating series is cancellation-prone and returns the
  # NaN sentinel once it cannot meet tolerance, which is how it hands over to
  # quadrature.  The contract is therefore one-sided: wherever the series does
  # produce a number, it must agree with quadrature.
  for (gamma in c(-1, -5, -20, -50, -200, -600)) {
    series <- mhn:::.pmhn_force(0.3, 2.5, 1, gamma, "series")
    quad   <- mhn:::.pmhn_force(0.3, 2.5, 1, gamma, "integrate")
    expect_true(is.finite(quad))
    if (!is.nan(series)) expect_equal(series, quad, tolerance = 1e-6)
    # Either way the public entry point must return a valid probability.
    p <- pmhn(0.3, 2.5, 1, gamma)
    expect_true(is.finite(p) && p >= 0 && p <= 1)
  }
})

test_that("a sampler that exhausts its retry budget reports once, not per draw", {
  # Successful draws are never NaN, so NaN is the failure signal; the wrapper
  # counts them and raises a single warning rather than one per element.
  set.seed(1)
  expect_silent(rmhn(500, alpha = 2, beta = 1, gamma = 1))
})

test_that("region D envelope pieces stay finite and usable at high tilt", {
  # An envelope piece integrates to exp(base) * expm1(slope * width) / slope,
  # and its inverse CDF needs log(1 + u (exp(slope * width) - 1)).  Both were
  # formed by exponentiating slope * width first, which overflows past about
  # 709.  The area became +Inf, which turned the whole piece-selection
  # cumulative sum into NaN; the inverse CDF returned +Inf for every u, so the
  # widest piece proposed nothing inside its own support.
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
  # alpha < 1 with gamma > 0 routes to RTDR under every method.  Before the
  # envelope-overflow fix, gamma = 100 gave a Kolmogorov-Smirnov p-value of
  # about 1e-121 and a mean 20 standard errors high.
  for (gamma in c(50, 100, 200, 1000)) {
    set.seed(11)
    x <- rmhn(4000, alpha = 0.3, beta = 1, gamma = gamma)
    expect_false(anyNA(x))
    expect_gt(stats::ks.test(x, function(q) pmhn(q, 0.3, 1, gamma))$p.value, 0.001)
    z <- (mean(x) - mhn_mean(0.3, 1, gamma)) /
         sqrt(mhn_var(0.3, 1, gamma) / length(x))
    expect_lt(abs(z), 5)
  }
})

test_that("Sun Algorithm 1 meets the Theorem 2e acceptance bound", {
  skip_on_cran()
  # Theorem 2e of Sun et al. (2023) guarantees acceptance >= 0.8 for alpha >= 4
  # and any gamma > 0.  The envelope constant K_1 had been transcribed from the
  # main-text statement of Theorem 1a, whose power has base sqrt(beta (alpha-1));
  # a direct maximisation of f_MHN / f_Normal gives sqrt(beta) (alpha-1), which
  # is also what the paper's own proof and Theorem 1c carry.  With the printed
  # form, log K_1 is short by (alpha-1)/2 * log(alpha-1), so the comparison
  # against K_2 picked the Normal proposal over most of the small-tilt half of
  # the region and acceptance fell to about 0.70.
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
  # u^(alpha-1) exp(-u^2 - |z| u).  That integrand's mass sits in a window of
  # width about 1/|z| around u = (alpha-1)/|z|, but the upper limit had been set
  # near |z|/2 -- the peak of exp(-u^2 + |z| u), the opposite sign.  A
  # fixed-order rule then sampled only the flat tail, and log Psi collapsed to
  # -Inf from |z| about 600, taking dmhn, pmhn and the moments with it.
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
  expect_equal(pmhn(1e-4, 2.5, 1, -1e4), 0.150855, tolerance = 1e-4)
})

test_that("the special-case tests are made on the scale-free tilt", {
  # The family depends on gamma only through Delta = gamma / sqrt(beta)
  # (Sun et al. 2023, Theorem 1c), so "is the tilt negligible?" cannot be a
  # test on |gamma| alone.  At beta = 4e-18 a gamma of 1e-8 is Delta = 5 -- an
  # ordinary tilt -- yet it was routed to the gamma = 0 closed form by both the
  # special-case detector and the Psi dispatcher, and every exported function
  # returned the untilted answer.
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
  # the mu cancels analytically: mu + sigma(a + e) = sigma*e.  Returning the
  # literal form gave negative draws -- outside (0, Inf) -- and collapsed
  # 2e5 draws onto 13 distinct doubles.
  for (gamma in c(-1e4, -1e6, -1e8, -1e10)) {
    set.seed(31)
    x <- rmhn(2e4, alpha = 1, beta = 1, gamma = gamma)
    expect_true(all(x > 0), info = sprintf("gamma=%g", gamma))
    expect_gt(length(unique(x)), 1e4)
    # As the tilt grows the law tends to Exponential(rate = |gamma|).
    expect_equal(mean(x), 1 / abs(gamma), tolerance = 0.05)
  }
})

test_that("Psi at alpha = 1 avoids the z^2/4 cancellation", {
  # log Psi[1/2, z] = log(2) + log(pi)/2 + z^2/4 + log Phi(z/sqrt(2)).  For
  # z << 0 the last two terms cancel to every significant digit; at z = -1e8
  # the absolute error reached 0.73 in a quantity of size 17.7.  Writing
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
  # central moment below 1: mhn_kurtosis returned -701.6, for a quantity that
  # cannot be less than -2, and mhn_skewness returned the wrong sign at 200.
  for (gamma in c(20, 100, 200, 1000)) {
    k <- mhn_kurtosis(3, 1, gamma)
    s <- mhn_skewness(3, 1, gamma)
    expect_true(is.finite(k)); expect_true(is.finite(s))
    # Excess kurtosis is bounded below by -2 for every distribution, and the
    # MHN is right-skewed, so both have a known sign here.
    expect_gt(k, -2)
    expect_gt(s, 0)
    # A strong positive tilt drives the law to a normal, so both go to 0.
    expect_lt(abs(k), 0.01)
    expect_lt(abs(s), 0.01)
  }
  # The Pearson inequality kurtosis >= skewness^2 - 2 must hold everywhere.
  for (p in list(c(0.5,1,-5), c(2,1,-30), c(3,1,0), c(5,1,2), c(10,0.1,5), c(100,1,-1))) {
    k <- mhn_kurtosis(p[1], p[2], p[3]); s <- mhn_skewness(p[1], p[2], p[3])
    expect_gte(k, s^2 - 2 - 1e-6, label = sprintf("alpha=%g beta=%g gamma=%g", p[1], p[2], p[3]))
  }
})

test_that("qmhn reaches the quantile for small alpha", {
  # For alpha < 1 the distribution function behaves like x^alpha near zero, so
  # the quantile of a small p sits near p^(1/alpha) -- 2.3e-40 at alpha = 0.1,
  # p = 1e-4.  The lower bracket started at sqrt(eps) and could be halved only
  # 30 times, reaching 1.4e-17, and the solver then returned that endpoint as
  # if it were the root.  pmhn of the answer came back as 0.019 instead of the
  # 1e-4 asked for.
  for (alpha in c(0.5, 0.2, 0.1, 0.05)) {
    for (p in c(1e-2, 1e-4, 1e-8)) {
      q <- qmhn(p, alpha, 1, 1)
      expect_true(is.finite(q))
      expect_gt(q, 0)
      expect_equal(pmhn(q, alpha, 1, 1), p, tolerance = 1e-5,
                   info = sprintf("alpha=%g p=%g", alpha, p))
    }
  }
  # Far enough down the true quantile is below the smallest representable
  # double; the search must stop cleanly there rather than raising.
  expect_no_error(qmhn(1e-8, 0.01, 1, 1))
})

test_that("a large positive tilt no longer raises from the Psi series", {
  # The Lemma 10 truncation length grows like z^2, and past a ceiling the
  # series raised an error with no alternative path -- so dmhn, pmhn, qmhn and
  # the moment helpers simply failed on ordinary parameter sets such as
  # (2.5, 1e-4, 20).  Psi depends on the tilt only through z = gamma/sqrt(beta)
  # and the integrand is well behaved for either sign, so the quadrature that
  # already served gamma < 0 answers these too.
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
  # The parameter set the audit named, where beta rather than gamma is extreme.
  expect_no_error(dmhn(1, 2.5, 1e-4, 20))
  expect_true(is.finite(dmhn(1, 2.5, 1e-4, 20)))
})

test_that("the mean stays exact at extreme tilt", {
  # E[X] is a ratio of two Fox-Wright values whose logarithms are both about
  # z^2/4 while their difference is only about log z, so subtracting them costs
  # an absolute error of (z^2/4) times the unit roundoff: 31% in the answer by
  # gamma = 1e8.  Past that point the mean comes from the same quadrature the
  # central moments use, taken in the variable centred on the peak.
  for (gamma in c(1e3, 1e5, 1e7, 1e8, 1e9, 1e10)) {
    got <- mhn_mean(2.5, 1, gamma)
    # A strong tilt drives the law to a normal centred on the mode, and the
    # mean approaches it to a relative O(1/gamma^2).
    mode_x <- (gamma + sqrt(gamma^2 + 8 * 1.5)) / 4
    expect_equal(got, mode_x, tolerance = 1e-8,
                 info = sprintf("gamma=%g", gamma))
  }
  # The ordinary range keeps the closed form and is unchanged.
  expect_equal(mhn_mean(2.5, 1, 0), 1.0139673601, tolerance = 1e-9)
  expect_equal(mhn_mean(2.5, 1, -1), 0.819124503718, tolerance = 1e-9)
})

test_that("the distribution function stays monotone past the support", {
  # Trimming the quadrature range below the anchor but not above it left the
  # upper panel running out to x, so once x left the support the rule saw
  # nothing there and F collapsed to the value of the first panel alone:
  # pmhn(q, 2.5, 1, -30) held at 1 until q = 3600, then dropped to 0.301 for
  # every larger q.
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
  # path, and it is the region the earlier repair missed: with no interior peak
  # the code fell back to the raw-moment expansion, which gave -58485 for a
  # kurtosis of 6, and where a peak did exist the integration window was taken
  # from the curvature and so was far too wide once the tilt squeezed the peak
  # towards the origin.
  for (alpha in c(0.2, 0.5, 1, 1.01, 1.2, 1.5, 3)) {
    for (gamma in c(-1000, -3000)) {
      expect_equal(mhn_kurtosis(alpha, 1, gamma), 6 / alpha, tolerance = 0.02,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
      expect_equal(mhn_skewness(alpha, 1, gamma), 2 / sqrt(alpha), tolerance = 0.02,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
})

test_that("mhn_mean keeps the closed form where it is exact", {
  # The cancellation in the ratio of two Fox-Wright values is one-sided: it
  # afflicts gamma > 0, where both logarithms grow like z^2/4.  Switching to
  # quadrature for gamma < 0 as well replaced an exact answer with a truncation
  # error of 1.4e-3.
  for (gamma in c(-940, -960, -1000, -1e4)) {
    ratio <- exp(mhn:::.mhn_log_normalizing_const(2.2, 1, gamma) -
                 mhn:::.mhn_log_normalizing_const(1.2, 1, gamma))
    expect_equal(mhn_mean(1.2, 1, gamma), ratio, tolerance = 1e-12,
                 info = sprintf("gamma=%g", gamma))
  }
})

test_that("the alpha = 1 quantile survives a strongly negative tilt", {
  # The truncation mass was formed as 1 - Phi(-mu/sigma), which underflows to
  # exactly 0 once -mu/sigma passes about 8.3 -- reached at gamma = -12 with
  # beta = 1 -- after which qmhn returned Inf for every p.  Taking it from
  # pnorm's own upper tail on the log scale, and assembling the answer as
  # sigma * t rather than mu + sigma * (t + z0), keeps both the mass and the
  # subtraction in range.
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
  # alpha = 2.5, beta = 1, gamma = 1 that is q = 8, where pmhn returned exactly
  # 0 and -Inf on the log scale for a true value near 5e-25.
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

  # And qmhn inverts that tail rather than inverting 1 - p, which is what made
  # it return Inf below p = 1e-16.
  for (tri in list(c(2.5, 1, 1), c(2, 1, 0), c(1, 1, 2), c(1, 1, -10))) {
    for (p in c(1e-8, 1e-16, 1e-20, 1e-30)) {
      q <- qmhn(p, tri[1], tri[2], tri[3], lower.tail = FALSE)
      expect_true(is.finite(q), info = sprintf("a=%g g=%g p=%g", tri[1], tri[3], p))
      expect_equal(pmhn(q, tri[1], tri[2], tri[3], lower.tail = FALSE), p,
                   tolerance = 1e-6)
    }
    # log.p, where the point of the flag is that p need not be representable.
    q <- qmhn(-100, tri[1], tri[2], tri[3], lower.tail = FALSE, log.p = TRUE)
    expect_equal(pmhn(q, tri[1], tri[2], tri[3], lower.tail = FALSE, log.p = TRUE),
                 -100, tolerance = 1e-6)
  }

  # The gamma = 0 case must agree with base R exactly, since both reduce to the
  # same qgamma call.
  for (p in c(1e-16, 1e-20, 1e-30, 1e-100)) {
    expect_equal(qmhn(p, 2, 1, 0, lower.tail = FALSE),
                 sqrt(stats::qgamma(p, 1, 1, lower.tail = FALSE)),
                 tolerance = 1e-12)
  }
})

test_that("the normalising constant keeps the mass below the smallest double", {
  # For gamma < 0 the normalising constant is a quadrature in x, and half the
  # mass of x^(alpha-1) sits below x = exp(-1/alpha).  At alpha = 1e-6 that is
  # exp(-1e6), so no abscissa can reach it: the quadrature saturated at
  # Psi = 1.3e3 against a true value of 2.0e6, and every density, probability
  # and moment at small alpha was scaled by the shortfall.
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

  # Tail probabilities are the user-visible consequence, and they were wrong by
  # the same factor -- 1545 times too large at alpha = 1e-6.
  tail_reference <- function(alpha, beta, gamma, q) {
    kernel <- function(u) exp(alpha * u - beta * exp(2 * u) + gamma * exp(u))
    U <- -40
    tail_int <- stats::integrate(kernel, U, 40, rel.tol = 1e-13,
                                 subdivisions = 10000L)$value
    upper <- stats::integrate(kernel, log(q), 40, rel.tol = 1e-13,
                              subdivisions = 10000L)$value
    upper / (exp(alpha * U) / alpha + tail_int)
  }

  for (alpha in c(1e-6, 1e-4, 1e-3, 0.01)) {
    for (gamma in c(-10, -1)) {
      for (q in c(0.01, 0.5)) {
        expect_equal(pmhn(q, alpha, 1, gamma, lower.tail = FALSE),
                     tail_reference(alpha, 1, gamma, q),
                     tolerance = 1e-6,
                     info = sprintf("alpha=%g gamma=%g q=%g", alpha, gamma, q))
      }
    }
  }
})

test_that("rmhn covers the whole support when the mode is tiny or alpha is small", {
  skip_on_cran()

  # Region A starts its left contact search at half the mode.  An absolute
  # floor of 1e-6 was applied as well, which put the start to the right of the
  # mode whenever the mode was smaller than 2e-6 and the search then returned a
  # point with the wrong sign of slope.  These all threw on the default path.
  for (par in list(c(4, 1, -1e7), c(4, 1, -3e6), c(3, 1, -1e7),
                   c(4, 1, -1e8), c(1.1, 1, -1e8))) {
    set.seed(20)
    draws <- rmhn(4000, par[1], par[2], par[3])
    expect_true(all(is.finite(draws)) && all(draws > 0))
    expect_equal(mean(draws), mhn_mean(par[1], par[2], par[3]),
                 tolerance = 0.05,
                 info = sprintf("alpha=%g gamma=%g", par[1], par[3]))
  }

  # The right contact search runs on the log axis, where the density carries
  # exp(2y) and overflows past y = 354.  The recommended start of Gao & Wang
  # (2025), Eq. (8) is m + sqrt(-2 delta / L''(m)), which for a flat mode lands
  # well beyond that -- 447 above the mode at alpha = 1e-5.  The slope came
  # back as NaN, the sign guard let it through because NaN fails every
  # comparison, and the plateau piece was dropped without any error, leaving a
  # sampler whose support excluded the mode and the entire right tail.
  for (par in list(c(1e-4, 1, 0.01), c(3e-4, 1, 0.01), c(1e-5, 1, -1),
                   c(1e-6, 1, 0.001))) {
    set.seed(21)
    draws <- rmhn(100000, par[1], par[2], par[3])
    expect_true(all(is.finite(draws)) && all(draws >= 0))
    for (q in c(0.005, 0.05)) {
      expected <- pmhn(q, par[1], par[2], par[3], lower.tail = FALSE)
      # Binomial three-sigma band, floored so that a probability of a few
      # draws in 1e5 does not make the test flaky.
      band <- max(3 * sqrt(expected * (1 - expected) / 100000), 3e-5)
      expect_lt(abs(mean(draws > q) - expected), band)
    }
  }
})

test_that("the CDF quadrature keeps the mass below the smallest double", {
  # The same unreachable head that the normalising constant had: for alpha < 1
  # with a negative tilt the lower-tail quadrature runs from 0, and half the
  # mass of x^(alpha-1) lies below x = exp(-1/alpha).  It was masked while the
  # normalising constant carried the same shortfall, and appeared as soon as
  # that was fixed -- F(1e-3) at (1e-6, 1, -1000) came back as 6.3e-4 against a
  # true value of 1.
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

  # The upper tail on the same path, which underflowed to 0 for every q because
  # its range was set from the Gaussian width (894 at beta = 1e-4) rather than
  # from the slope, leaving the whole mass in a sliver the rule never sampled.
  for (q in c(1e-6, 1e-4, 1e-3, 1e-2)) {
    expect_equal(pmhn(q, 0.05, 1e-4, -1000, lower.tail = FALSE),
                 stats::pgamma(q, 0.05, rate = 1000, lower.tail = FALSE),
                 tolerance = 1e-5, info = sprintf("q=%g", q))
  }
  # qmhn inverted that, so it returned its own initial bracket -- the same
  # number for every p, with no warning.
  for (p in c(1e-2, 1e-3, 1e-4, 1e-6)) {
    expect_equal(qmhn(p, 0.05, 1e-4, -1000, lower.tail = FALSE),
                 stats::qgamma(p, 0.05, rate = 1000, lower.tail = FALSE),
                 tolerance = 1e-5, info = sprintf("p=%g", p))
  }
})

test_that("the variance stays accurate on both sides of the Lemma 2c cancellation", {
  # alpha/(2 beta) and mu gamma/(2 beta) are each about alpha/2 while the
  # variance is about alpha/gamma^2, so the closed form amplifies the rounding
  # of its own terms by gamma^2/(2 beta alpha): it was 17 percent high at
  # gamma = -1e7 and clamped to exactly 0 by -1e8.  Under a strong tilt the law
  # tends to Gamma(alpha, |gamma|), whose variance is alpha/gamma^2.
  for (alpha in c(1e-4, 0.01, 0.3, 1, 3, 10)) {
    for (gamma in c(-1e4, -1e6, -1e8)) {
      expect_equal(mhn_var(alpha, 1, gamma), alpha / gamma^2,
                   tolerance = 1e-3,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
  # Away from the cancellation the closed form is the accurate one and must
  # still be what is used, to the last bit.
  for (alpha in c(1e-4, 0.01, 0.3, 1, 3, 10)) {
    for (gamma in c(-10, -1, 0.5, 100)) {
      mu <- mhn_mean(alpha, 1, gamma)
      expect_equal(mhn_var(alpha, 1, gamma),
                   alpha / 2 + mu * (gamma / 2 - mu),
                   tolerance = 1e-12,
                   info = sprintf("alpha=%g gamma=%g", alpha, gamma))
    }
  }
})

test_that("pmhn agrees element-wise and vectorised under every tail flag", {
  # pmhn dispatches the special cases twice, and the two implementations
  # disagreed: the vector-parameter path derived the truncated-normal lower tail
  # as 1 - upper, which is exactly 0 once the lower tail falls below the
  # rounding of 1.  pmhn(3.0184, 1, 1, 18) was 1.35e-17 scalar and 0 vectorised.
  #
  # The comparison is made on whichever tail is the small one.  A probability
  # that rounds to 1 carries no information about its own distance from 1, so
  # its logarithm is meaningless to compare at any tolerance -- at q = 2 the
  # lower tail is 1 - 2.7e-11 and log(F) is limited to about 4e-6 relative by
  # the double that holds F.  It is the small tail that a wrong branch destroys,
  # and the small tail that has to survive.
  for (gamma in c(18, 30, 10, -10)) {
    for (q in c(0.5, 1, 2, 3.0184056)) {
      for (lt in c(TRUE, FALSE)) {
        for (lp in c(TRUE, FALSE)) {
          scalar <- pmhn(q, 1, 1, gamma, lower.tail = lt, log.p = lp)
          vector <- pmhn(q, alpha = c(1, 1), beta = 1, gamma = gamma,
                         lower.tail = lt, log.p = lp)[1]
          info <- sprintf("g=%g q=%g lower.tail=%s log.p=%s", gamma, q, lt, lp)
          if ((lp && scalar > -0.693) || (!lp && scalar > 0.5)) {
            expect_lt(abs(vector - scalar), 1e-8, label = info)
          } else {
            expect_equal(vector, scalar, tolerance = 1e-10, info = info)
          }
        }
      }
    }
  }
})

test_that("the RTDR region boundary at alpha exactly 1 is sampled, not rejected", {
  skip_on_cran()

  # Region A needs an interior mode in x, and at alpha = 1 there is none: the
  # mode is 0 for a non-positive tilt, log f(0) is -Inf, and the envelope
  # degenerated to a plateau of height -Inf spanning [0, Inf].  Every proposal
  # was rejected and the draws came back NaN.  The public API is shielded by the
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
  # The pre-scan tested the bare floats while the dispatcher routes on the
  # scale-free predicates, so it rejected parameters the sun path handles
  # perfectly well -- alpha within sqrt(eps) of 1, or a negligible tilt.
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
  # of order one.  Assembled directly, the result carried an absolute error of
  # about eps gamma^2/(4 beta): 3.7e-8 at |z| = 3e4, and 4.6 percent at 3e7.
  # Both sides now have that term removed analytically.
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

  # The scalar fast path and the recycling loop share one implementation, so
  # they cannot disagree the way pmhn's two truncated-normal paths once did.
  for (par in list(c(3, 1, 3e4), c(3, 1e-6, 3e3), c(2.5, 1, 5))) {
    alpha <- par[1]; beta <- par[2]; gamma <- par[3]
    xs <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta) *
      c(0.5, 0.9, 1, 1.1, 2)
    expect_equal(dmhn(xs, alpha = c(alpha, alpha), beta = beta, gamma = gamma,
                      log = TRUE),
                 dmhn(xs, alpha, beta, gamma, log = TRUE),
                 tolerance = 1e-14)
  }
})

test_that("the CDF cancellation guard fires before the tolerance is spent", {
  # The guard charged one unit of roundoff to each accumulator, but the two sums
  # run over i_max terms and their errors accumulate as a random walk.  Without
  # the sqrt(i_max) it admitted results some twenty times worse than the
  # tolerance asked for -- at (10, 1, -100) the dispatcher returned 0.542194283
  # where an independent quadrature gives 0.542194428.
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
  # The cheapest check there is on the survival function, and the one that would
  # have caught the defect below the moment it was introduced.
  #
  # log_ccdf_integrate anchored its range at max(peak, q).  For q at or below
  # the mode that anchor IS the mode, where the slope is zero by definition, so
  # a reach derived from the slope collapsed to the Gaussian width -- 894 at
  # beta = 1e-4 where the kernel dies within 0.08 -- and the panel from the mode
  # outward returned nothing.  pmhn(q, 3, 1e-4, -1000, lower.tail = FALSE) then
  # reported only the mass between q and the mode: the two tails summed to
  # F(mode) = 0.323 for every q below it, and qmhn inverted that, returning the
  # denormal floor for every p.
  #
  # A second failure of the same shape sat at alpha < 1 with a tiny q, where the
  # range spanned sixty decades. Both are gone now that the integral is taken on
  # the log axis with panels that double in width away from the peak.
  for (alpha in c(0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.5, 3, 10)) {
    for (beta in c(1e-4, 0.01, 1, 100)) {
      for (gamma in c(-1000, -100, -10, -1, 0, 1, 10, 100)) {
        mu <- tryCatch(mhn_mean(alpha, beta, gamma), error = function(e) NA_real_)
        if (!is.finite(mu) || mu <= 0) next
        for (f in c(1e-60, 1e-6, 0.1, 0.5, 1, 2, 10)) {
          q <- f * mu
          lo <- pmhn(q, alpha, beta, gamma)
          up <- pmhn(q, alpha, beta, gamma, lower.tail = FALSE)
          expect_lt(abs(lo + up - 1), 1e-6,
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
      expect_equal(pmhn(q, alpha, beta, gamma, lower.tail = FALSE), p,
                   tolerance = 1e-5,
                   info = sprintf("alpha=%g beta=%g gamma=%g p=%g",
                                  alpha, beta, gamma, p))
      expect_equal(qmhn(p, alpha, beta, gamma, lower.tail = FALSE), q,
                   tolerance = 1e-5,
                   info = sprintf("alpha=%g beta=%g gamma=%g p=%g",
                                  alpha, beta, gamma, p))
    }
  }
})

test_that("the density vanishes at infinity on every path", {
  # Every base R density returns 0 there.  The general expression gave NaN,
  # since (alpha-1) log(x) and gamma x are +Inf while -beta x^2 is -Inf; the
  # sqrt-Gamma branch gave NaN for the same reason.  Only the truncated-normal
  # branch was right, because it delegates to dnorm.
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

test_that("an exhausted retry budget is reported once, and only then", {
  skip_on_cran()

  # rmhn's only user-visible failure signal.  The counter behind it was guarded
  # by `ISNAN(x) && !Rcpp::NumericVector::is_na(x)`, and for doubles that second
  # test is the first one negated, so the condition was never true and the
  # sampler returned NaN draws in silence.  Nothing in the suite exercised the
  # warning, because nothing in the suite reached a parameter set that exhausts
  # the budget.
  #
  # This one does: at alpha = 0.7 with beta = 1e-8 and gamma = 1e8 the
  # standardised tilt is 1e12, and no envelope proposes an acceptable draw
  # inside the retry limit.  The contract is that the call still returns, one
  # vector of NaN, with exactly one warning naming how many draws failed.
  set.seed(2)
  warnings_seen <- character(0)
  draws <- withCallingHandlers(
    rmhn(50, 0.7, 1e-8, 1e8),
    warning = function(cnd) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_length(warnings_seen, 1L)
  expect_match(warnings_seen[1], "exhausted its retry budget for 50 of 50")
  expect_length(draws, 50L)
  expect_true(all(is.nan(draws)))

  # And it stays quiet when the sampler is working, which is the other half of
  # "reported once": a per-draw warning would fire fifty times here.
  set.seed(2)
  expect_silent(rmhn(50, 2, 1, 1))
  set.seed(2)
  expect_silent(rmhn(50, 0.3, 1, 100))
})

test_that("a non-finite shape or rate is refused, in the caller's vocabulary", {
  # alpha > 0 is true of Inf, so an infinite shape passed validation and reached
  # the series, where a truncation length computed from it overflowed the
  # size_t handed to std::vector.  What the user saw was that container's own
  # exception, whose entire message is the word "vector".
  for (f in list(dmhn, pmhn)) {
    expect_error(f(1, Inf, 1, 1), "alpha must be positive and finite")
    expect_error(f(1, 1, Inf, 1), "beta must be positive and finite")
    expect_error(f(1, NaN, 1, 1), "alpha must be positive")
  }
  expect_error(qmhn(0.5, Inf, 1, 1), "alpha must be positive and finite")
  for (f in list(mhn_mean, mhn_var, mhn_skewness, mhn_kurtosis, mhn_mode)) {
    expect_error(f(Inf, 1, 1), "alpha must be positive and finite")
    expect_error(f(1, Inf, 1), "beta must be positive and finite")
  }
  # The vector path too, where the offending element is not the first.
  expect_error(dmhn(1, alpha = c(2, Inf), beta = 1, gamma = 1),
               "alpha must be positive and finite")

  # rmhn keeps the contract its documentation states: an NA draw, not an error.
  set.seed(1)
  expect_silent(x <- rmhn(4, alpha = c(2, Inf, 2, NaN), beta = 1, gamma = 1))
  expect_equal(is.na(x), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("the distribution function approaches its normal limit at a strong tilt", {
  # win-builder found this on all three Windows R versions while the local check
  # passed: pmhn at the mode returned 1 rather than 0.5.
  #
  # log Psi and the log of the integral over [0, x] are each about
  # gamma^2/(4 beta) -- 2.5e13 at gamma = 1e7 -- while the log probability they
  # combine into is of order one.  Assembled directly the answer carries the
  # rounding of both, and how that rounding falls depends on the platform's
  # quadrature: here it cost four digits, there enough to be clamped to 1.
  # Both sides now have the term removed analytically, and the integrand is
  # written in the variable centred on its own anchor.
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
