# The closed forms of mhn_stable.h at the two ends of the double range: where
# sqrt(8 beta (alpha - 1)) passes the largest double, and where gamma^2 and
# 8 beta (alpha - 1) fall below the smallest normal one.
#
# Several blocks use the scale identity of the family: sqrt(beta) X is
# MHN(alpha, 1, gamma / sqrt(beta)).  At beta = 1 none of these quantities
# leaves the range.  Tiny and huge values are compared as ratios:
# expect_equal() compares absolute differences once the expected value is
# below its tolerance.

test_that("mhn_mode holds where sqrt(8 beta (alpha - 1)) itself overflows", {
  # With alpha - 1 = beta = b in double precision and gamma = s b, the mode
  # is the positive root of 2 x^2 - s x - 1 = 0, whatever b is.
  b <- 1e308
  expect_identical(b - 1, b)
  expect_equal(mhn_mode(b, b, -b), 0.5, tolerance = 1e-14)
  expect_equal(mhn_mode(b, b, b), 1, tolerance = 1e-14)
  expect_equal(mhn_mode(b, b, 1e300), (1e-8 + sqrt(1e-16 + 8)) / 4,
               tolerance = 1e-14)
  # alpha - 1 = 2 beta = gamma: x^2 - x - 1 = 0
  expect_equal(mhn_mode(b, b / 2, b), (1 + sqrt(5)) / 2, tolerance = 1e-14)
  # alpha - 1 = 6.8 beta, gamma = -6.8 beta: 2 x^2 + 6.8 x - 6.8 = 0
  expect_equal(mhn_mode(1.7e308, 2.5e307, -1.7e308),
               (-6.8 + sqrt(6.8^2 + 8 * 6.8)) / 4, tolerance = 1e-14)
})

test_that("mhn_mode holds where gamma^2 and 8 beta (alpha - 1) underflow", {
  # Lemmas 3b and 3c of Sun et al. (2023) at beta = 1, in the form that only
  # adds for a negative tilt.
  unit_mode <- function(alpha, z) {
    d <- sqrt(z^2 + 8 * (alpha - 1))
    if (z >= 0) (z + d) / 4 else 2 * (alpha - 1) / (d - z)
  }
  cases <- list(c(2, 1e-322, 3e-161), c(2, 1e-320, -1e-160),
                c(1 + 2^-52, 1e-320, 1e-170), c(1 + 2^-52, 1e-320, -1e-170),
                c(1 + 2^-52, 1e-300, 4e-158),
                c(0.5, 5e-324, 5e-162), c(1 - 1e-12, 1e-310, 3e-161),
                c(1 - 2^-53, 1e-320, 1e-167))
  for (p in cases) {
    s <- sqrt(p[2])
    expect_equal(mhn_mode(p[1], p[2], p[3]) * s, unit_mode(p[1], p[3] / s),
                 tolerance = 1e-12,
                 info = sprintf("(%.17g, %g, %g)", p[1], p[2], p[3]))
  }
})

test_that("dmhn and pmhn hold with alpha and beta near the largest double", {
  # alpha - 1 = beta = gamma = b: the mode is 1 and the law is narrower than
  # the spacing of the doubles there.  log f(x) - log f(1) is
  # b (log x - x^2 + x), and log f(1), about 354, is below its rounding.
  b <- 1e308
  x <- c(0.5, 0.9, 1.1)
  expect_equal(dmhn(x, b, b, b, log = TRUE) / (b * (log(x) - x^2 + x)),
               rep(1, 3), tolerance = 1e-10)
  expect_identical(pmhn(c(0.9, 1.1), b, b, b), c(0, 1))
  # alpha - 1 = 2 beta = gamma: the mode m is the golden ratio, and
  # log f(x) - log f(m) is b (log(x / m) - (x^2 - m^2) / 2 + (x - m)).
  m <- (1 + sqrt(5)) / 2
  expect_identical(pmhn(c(0.9, 1.1) * m, b, b / 2, b), c(0, 1))
  expect_equal(dmhn(0.9 * m, b, b / 2, b, log = TRUE) /
                 (b * (log(0.9) + 0.095 * m^2 - 0.1 * m)),
               1, tolerance = 1e-10)
})

test_that("dmhn, pmhn and qmhn at a subnormal beta match the law at beta = 1", {
  # The largest difference between a quantity at (alpha, beta, gamma) and the
  # same quantity at (alpha, 1, gamma / sqrt(beta)), at the quantiles of the
  # second law: differences of logs for the density and the two tails,
  # relative differences for the two quantiles.
  scale_gap <- function(alpha, beta, gamma, probs) {
    s <- sqrt(beta)
    z <- gamma / s
    # a quantile that underflows at beta = 1 is 0 at both scales
    rel <- function(a, b) ifelse(b == 0, abs(a), abs(a / b - 1))
    q1 <- qmhn(probs, alpha, 1, z)
    x <- q1 / s
    c(density = max(abs(dmhn(x, alpha, beta, gamma, log = TRUE) - log(s) -
                          dmhn(q1, alpha, 1, z, log = TRUE))),
      lower = max(abs(pmhn(x, alpha, beta, gamma, log.p = TRUE) -
                        pmhn(q1, alpha, 1, z, log.p = TRUE))),
      upper = max(abs(pmhn(x, alpha, beta, gamma, lower.tail = FALSE,
                           log.p = TRUE) -
                        pmhn(q1, alpha, 1, z, lower.tail = FALSE,
                             log.p = TRUE))),
      quantile = max(rel(qmhn(probs, alpha, beta, gamma) * s, q1)),
      upper_quantile = max(rel(
        qmhn(probs, alpha, beta, gamma, lower.tail = FALSE) * s,
        qmhn(probs, alpha, 1, z, lower.tail = FALSE))))
  }
  probs <- c(1e-6, 0.05, 0.5, 0.95, 1 - 1e-9)
  cases <- list(
    c(2, 1e-322, 3e-161),            # the upper tail is centred on the peak
    c(0.5, 5e-324, 5e-162),          # a Lemma 3c peak
    c(3e4, 1e-323, 1e-161),          # measured from the peak
    c(270253.9, 6.85e-318, -5.42e-171)
  )
  for (p in cases) {
    gap <- scale_gap(p[1], p[2], p[3], probs)
    expect_lt(max(gap), 1e-7)
  }
  # A shape so small that the quantiles below 0.9 underflow at beta = 1.
  gap <- scale_gap(3.045e-4, 3.46e-323, -6.7e-174, c(0.9, 0.95, 0.99))
  expect_lt(max(gap), 1e-7)
  # Next to the origin the lower tail is a closed form up to the point where
  # beta x^2 + |gamma| x reaches 4.5e-9; at a subnormal beta that point must
  # still be found.  beta x^2 is 0.003 and 0.06 here.
  gap <- scale_gap(150, 1e-320, 1e-170, c(1e-300, 1e-200))
  expect_lt(max(gap[c("lower", "quantile")]), 1e-7)
})

test_that("pmhn's upper tail at a subnormal beta matches quadrature", {
  # For alpha = 2 the kernel at beta = 1 is t exp(-t^2 + z t), and its tail
  # integrals are taken by integrate(), independently of the package.
  alpha <- 2; beta <- 1e-322; gamma <- 3e-161
  s <- sqrt(beta)
  z <- gamma / s
  kernel <- function(t) t * exp(-t^2 + z * t)
  total <- integrate(kernel, 0, Inf, rel.tol = 1e-12)$value
  for (x in c(1, 2.5, 4)) {
    upper <- integrate(kernel, x, Inf, rel.tol = 1e-12)$value / total
    expect_equal(pmhn(x / s, alpha, beta, gamma, lower.tail = FALSE) / upper,
                 1, tolerance = 1e-8)
  }
})
