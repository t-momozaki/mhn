# Tests for pmhn() -- the MHN cumulative distribution function.

# ============================================================
# 1. Boundary values
# ============================================================

test_that("pmhn(0) == 0", {
  for (a in c(0.5, 1, 2, 5)) for (b in c(0.5, 1, 3)) for (g in c(-2, 0, 1, 3)) {
    expect_equal(pmhn(0, a, b, g), 0,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("pmhn(Inf) == 1", {
  for (a in c(0.5, 1, 2)) for (b in c(0.5, 1)) for (g in c(-2, 0, 2)) {
    expect_equal(pmhn(Inf, a, b, g), 1,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("pmhn(q < 0) == 0", {
  expect_equal(pmhn(-1, 2, 1, 1), 0)
  expect_equal(pmhn(c(-Inf, -2, -0.1), 2, 1, 1), c(0, 0, 0))
})

test_that("pmhn handles log.p at boundaries", {
  expect_equal(pmhn(0, 2, 1, 1, log.p = TRUE), -Inf)
  expect_equal(pmhn(Inf, 2, 1, 1, log.p = TRUE), 0)
  expect_equal(pmhn(0, 2, 1, 1, lower.tail = FALSE), 1)
  expect_equal(pmhn(Inf, 2, 1, 1, lower.tail = FALSE), 0)
})

# ============================================================
# 2. Monotonicity
# ============================================================

test_that("pmhn is monotone non-decreasing in q", {
  q <- seq(0, 8, length.out = 200)
  params <- expand.grid(alpha = c(0.5, 1, 1.5, 2, 4),
                        beta  = c(0.5, 1, 3),
                        gamma = c(-3, -1, 0, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    info <- sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g)
    F_q <- pmhn(q, a, b, g)
    # Over the whole grid, F = 1 included: above the peak of the density of
    # log X the lower tail is one minus the directly computed upper tail, so
    # it cannot wobble where F approaches 1.  The count of points strictly
    # inside (1e-6, 1 - 1e-6) keeps a CDF collapsed to 0 or 1 from passing.
    expect_false(is.unsorted(F_q), label = info)
    expect_gt(sum(F_q > 1e-6 & F_q < 1 - 1e-6), 10, label = info)
    # The upper tail is non-increasing, down to its smallest values.
    expect_false(is.unsorted(rev(pmhn(q, a, b, g, lower.tail = FALSE))),
                 label = info)
  }
})

# ============================================================
# 3. Tail complement: pmhn(q) + pmhn(q, lower.tail = FALSE) == 1
# ============================================================

test_that("lower.tail complement", {
  q <- c(0.1, 0.5, 1, 2, 4)
  params <- expand.grid(alpha = c(0.5, 1, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    lower <- pmhn(q, a, b, g)
    upper <- pmhn(q, a, b, g, lower.tail = FALSE)
    expect_equal(lower + upper, rep(1, length(q)), tolerance = 1e-10,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 4. log.p consistency: pmhn(log.p=TRUE) == log(pmhn())
# ============================================================

test_that("log.p == log of linear value", {
  q <- c(0.1, 0.5, 1, 2, 4)
  params <- expand.grid(alpha = c(0.5, 1, 2),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    linear <- pmhn(q, a, b, g)
    logp   <- pmhn(q, a, b, g, log.p = TRUE)
    # Compare where linear > 0
    idx <- linear > 0
    # Element by element.  Where linear is close to 1, log(linear) carries
    # linear's own error absolutely, a few units of the spacing of the
    # doubles just below 1 (1.1e-16), which logp does not.  The comparison
    # allows four such units; on this grid 29 of the 180 elements need some
    # of it, at most 5.2e-17.
    expect_close(logp[idx], log(linear[idx]), 1e-12, abs_tol = 4e-16,
                 info = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
    # Upper-tail log.  The direct path (lower.tail=FALSE, log.p=TRUE) is
    # often more accurate than log1p(-linear) because the latter loses
    # precision when linear is close to 1; allow a moderate tolerance.
    upper_logp <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
    idx2 <- linear < 1 - 1e-10
    expect_close(upper_logp[idx2], log1p(-linear[idx2]), 1e-6,
                 info = sprintf("upper a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 5. pmhn against the kernel reference
#    mhn_ref() of helper-mhn-cdf.R integrates the kernel alone, to a relative
#    1e-12 (1e-10 where rounding prevents that), and calls no function of
#    the package.
# ============================================================

test_that("pmhn matches the kernel reference", {
  params <- expand.grid(alpha = c(0.5, 1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-3, -1, 1, 3))
  q_vals <- c(0.1, 0.5, 1, 2, 4)
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    expect_close(pmhn(q_vals, a, b, g), mhn_cdf(a, b, g)(q_vals), 1e-9,
                 info = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 6. Special-case dispatch
# ============================================================

test_that("gamma = 0: pmhn == pgamma(q^2, alpha/2, scale = 1/beta)", {
  q <- c(0.01, 0.5, 1, 2, 5)
  for (a in c(0.5, 1, 1.5, 2, 4)) for (b in c(0.5, 1, 3)) {
    ref <- pgamma(q^2, shape = a/2, scale = 1/b)
    got <- pmhn(q, alpha = a, beta = b, gamma = 0)
    expect_equal(got, ref, tolerance = 1e-12,
                 label = sprintf("a=%.1f b=%.1f", a, b))
  }
})

test_that("alpha = 1: pmhn == truncated normal CDF", {
  q <- c(0.01, 0.5, 1, 2, 5)
  for (b in c(0.5, 1, 2)) for (g in c(-2, 0, 1, 3)) {
    mu <- g / (2 * b); sigma <- 1 / sqrt(2 * b)
    ref <- (pnorm(q, mu, sigma) - pnorm(0, mu, sigma)) /
           pnorm(0, mu, sigma, lower.tail = FALSE)
    got <- pmhn(q, alpha = 1, beta = b, gamma = g)
    expect_equal(got, ref, tolerance = 1e-10,
                 label = sprintf("b=%.1f g=%+.1f", b, g))
  }
})

# ============================================================
# 7. Vectorization / recycling
# ============================================================

test_that("pmhn vectorizes q", {
  q <- c(0.1, 0.5, 1, 2)
  out <- pmhn(q, 2, 1, 1)
  expect_length(out, 4)
  expect_equal(out, sapply(q, function(x) pmhn(x, 2, 1, 1)), tolerance = 1e-14)
})

test_that("pmhn recycles parameter vectors", {
  q <- c(0.5, 1.0, 1.5)
  a <- c(1, 2)
  b <- 1
  g <- c(0, 1, -1)
  n <- max(length(q), length(a), length(b), length(g))
  out <- pmhn(q, a, b, g)
  expect_length(out, n)
  # Reference: element-by-element scalar calls
  ref <- vapply(seq_len(n), function(i) {
    pmhn(q[((i-1) %% length(q)) + 1],
         a[((i-1) %% length(a)) + 1],
         b[((i-1) %% length(b)) + 1],
         g[((i-1) %% length(g)) + 1])
  }, numeric(1))
  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("pmhn with empty q returns numeric(0)", {
  expect_equal(pmhn(numeric(0), 2, 1, 1), numeric(0))
})

# ============================================================
# 8. NA / NaN propagation
# ============================================================

test_that("pmhn propagates NA in q", {
  out <- pmhn(c(0.5, NA, 1), 2, 1, 1)
  expect_true(is.na(out[2]))
  expect_false(is.na(out[1]))
  expect_false(is.na(out[3]))
})

test_that("pmhn gives NaN with a warning for parameters outside the parameter space", {
  # As pgamma does; test-input-validation.R tests the rule in full.
  for (args in list(list(alpha = 0), list(alpha = -1), list(beta = 0),
                    list(gamma = Inf))) {
    expect_warning(out <- do.call(pmhn, c(list(1), args)), "NaNs produced")
    expect_identical(out, NaN)
  }
  expect_identical(pmhn(1, alpha = NA, beta = 1, gamma = 0), NA_real_)
})

# ============================================================
# 9. Series-vs-integration cross-check across the accuracy grid
# ============================================================
#
# cdf_series_accuracy_reference.csv, beside this file, holds 396 points
# (alpha, beta, gamma, q) and the reference F(q) at each, `val_ref_quad`: an
# arbitrary-precision quadrature (Rmpfr, 200 bits) on the log axis.  It shares
# no code with the package: substituting x = exp(u) lets it compute the
# normalising constant from scratch alongside the numerator, so there is no
# Psi, no series, and no package call anywhere in it.  The values were
# computed on 2026-09-28 and are fixed; the test asserts that the dispatch
# (`pmhn(q, alpha, beta, gamma)`) reproduces them at every grid point.

test_that("pmhn matches an independent reference across the cdf_series_accuracy grid", {
  df <- utils::read.csv(test_path("cdf_series_accuracy_reference.csv"))
  expect_identical(nrow(df), 396L)

  # Skip rows where F sits at the double-precision floor -- below 1e-10, or
  # within 1e-10 of 1 -- since there a relative comparison measures rounding
  # rather than accuracy.  The two tails are tested directly, and to much
  # greater depth, in test-numerical-stability.R.
  skip_row <- is.na(df$val_ref_quad) |
                df$val_ref_quad < 1e-10 |
                (1 - df$val_ref_quad) < 1e-10
  df <- df[!skip_row, ]
  expect_gt(nrow(df), 100)

  # 1e-8 is the working tolerance the package is written to, and the
  # dispatcher meets it across this grid with several digits to spare.  A
  # series summed past its cancellation limit, or a quadrature that misses the
  # mass, is off by 0.1 to 1 or returns exactly zero, far outside this bound.
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    val <- pmhn(r$q, alpha = r$alpha, beta = r$beta, gamma = r$gamma)
    expect_lt(
      abs(val - r$val_ref_quad) / max(abs(r$val_ref_quad), 1e-300),
      1e-8,
      label = sprintf("pmhn(q=%g, alpha=%g, beta=%g, gamma=%g)",
                      r$q, r$alpha, r$beta, r$gamma)
    )
  }
})

# ============================================================
# 10. Scale equivariance in beta
# ============================================================

test_that("pmhn depends on beta only through sqrt(beta) q and z = gamma / sqrt(beta)", {
  # If X ~ MHN(alpha, beta, gamma) then sqrt(beta) X ~ MHN(alpha, 1, z) with
  # z = gamma / sqrt(beta) (by a change of variable in the density of
  # Sun et al. 2023, Lemma 1a), so
  # pmhn(y / sqrt(beta), alpha, beta, z sqrt(beta)) must equal
  # pmhn(y, alpha, 1, z) for every beta.  The quadrature ranges scale exactly
  # with 1 / sqrt(beta), so the two sides differ only by rounding, and 1e-10 in
  # absolute terms leaves ample room for that while catching any dependence of
  # the quadrature's refinement on the width of its range.
  #
  # At (20, -10) the lower tail comes from the Gauss-Kronrod quadrature at all
  # but the first few points, where the Lemma 1b series answers; (0.5, 5) has
  # an interior peak with alpha < 1 and is answered by the series throughout.
  # The upper tail is always quadrature, on the log axis.
  betas <- 10^c(-16, -8, 0, 8, 16)
  grid <- function(a, z) seq(0.02, 3, length.out = 40) * qmhn(0.999, a, 1, z)
  rescaled <- function(a, z, b, y, ...) pmhn(y / sqrt(b), a, b, z * sqrt(b), ...)
  for (p in list(c(20, -10), c(0.5, 5))) {
    a <- p[1]; z <- p[2]
    y <- grid(a, z)
    lower1 <- pmhn(y, a, 1, z)
    upper1 <- pmhn(y, a, 1, z, lower.tail = FALSE)
    for (b in betas) {
      info <- sprintf("alpha=%g z=%g beta=%g", a, z, b)
      expect_lt(max(abs(rescaled(a, z, b, y) - lower1)), 1e-10, label = info)
      expect_lt(max(abs(rescaled(a, z, b, y, lower.tail = FALSE) - upper1)),
                1e-10, label = info)
    }
  }

  # At (1.5, -20) and (3, -50) pmhn hands the lower part of the grid to the
  # series.  For gamma < 0 its terms alternate in sign, and the rounding they
  # leave -- up to about 5e-10, inside the 1e-8 working tolerance -- differs
  # from one scale to another.  There the upper tail is compared through pmhn
  # and the lower-tail quadrature directly, over the whole grid.
  for (p in list(c(1.5, -20), c(3, -50))) {
    a <- p[1]; z <- p[2]
    y <- grid(a, z)
    quad <- function(b) {
      vapply(y / sqrt(b), function(q) mhn:::.pmhn_force(q, a, b, z * sqrt(b), "integrate"),
             numeric(1))
    }
    quad1 <- quad(1)
    upper1 <- pmhn(y, a, 1, z, lower.tail = FALSE)
    for (b in betas) {
      info <- sprintf("alpha=%g z=%g beta=%g", a, z, b)
      expect_lt(max(abs(quad(b) - quad1)), 1e-10, label = info)
      expect_lt(max(abs(rescaled(a, z, b, y, lower.tail = FALSE) - upper1)),
                1e-10, label = info)
    }
  }
})

test_that("pmhn is non-decreasing far into the upper tail at a tiny beta", {
  # (1.5, 1e-14, -2e-6) is (alpha, z) = (1.5, -20) with the axis stretched by
  # 1e7.  These q lie at 0.74, 1, 2 and 5 on the beta = 1 axis, where the
  # upper tail is 1e-6, 3.6e-9, 4.7e-19 and 3.9e-54, and the lower-tail
  # quadrature answers each of them.  The distribution function must be
  # non-decreasing there and equal to its beta = 1 values.
  y <- c(0.735194438, 1, 2, 5)
  v <- pmhn(1e7 * y, 1.5, 1e-14, -2e-6)
  expect_false(is.unsorted(v))
  expect_lt(max(abs(v - pmhn(y, 1.5, 1, -20))), 1e-10)
})

# ============================================================
# 11. Strong tilts: the quadrature path against exact limits
# ============================================================
#
# For gamma < 0 the MHN(alpha, beta, gamma) law is Gamma(alpha, rate = -gamma)
# tilted by exp(-beta x^2), and under that Gamma law E(beta X^2) =
# alpha (alpha + 1) beta / gamma^2, so the two agree up to a relative
# correction of that order, and base R's pgamma and qgamma are then exact
# references from outside the package.  Every triple below that is compared
# with them has alpha (alpha + 1) beta / gamma^2 <= 1e-8.

test_that("pmhn matches the Gamma limit when the mode lies below sqrt(eps)", {
  # The modes (alpha - 1) / |gamma| run from 1e-9 to 1e-8 here, below
  # sqrt(.Machine$double.eps) = 1.5e-8, and the q from 1.5e-10 to 5e-7.  The
  # quadrature has to normalise the kernel at the point it expands it about,
  # however close to the origin that point is.  The tolerance is the working
  # tolerance sqrt(eps) of the quadrature, rounded up.
  p <- c(0.01, 0.5, 0.99)
  for (par in list(c(2, -1e8), c(2, -1e9), c(5, -1e9), c(1.5, -1e8),
                   c(1.01, -1e7), c(3, -3e8))) {
    a <- par[1]; g <- par[2]
    q <- stats::qgamma(p, a, rate = -g)
    info <- sprintf("alpha=%g gamma=%g", a, g)
    expect_equal(pmhn(q, a, 1, g), p, tolerance = 2e-8, info = info)
    expect_equal(pmhn(q, a, 1, g, lower.tail = FALSE), 1 - p,
                 tolerance = 2e-8, info = info)
  }
  # q below sqrt(eps) while the mode is above it: (3, 1, -1e8) has its mode at
  # 2e-8, and F(2e-9) is pgamma(0.2, 3) = 1.15e-3.
  expect_equal(pmhn(2e-9, 3, 1, -1e8), stats::pgamma(0.2, 3), tolerance = 2e-8)
})

test_that("pmhn is accurate in relative terms deep in the left tail", {
  # At (3, 1, -1000) the mode is 2e-3 and these q lie far below sqrt(eps).
  # References: the ratio of int_0^q t^2 exp(-t^2 - 1000 t) dt to the same
  # integral over (0, Inf), each by Romberg integration in Rmpfr at 320 bits
  # in the variable u = 1000 t.  They lie 1.2e-5 above pgamma(q, 3, 1000),
  # the O(beta / gamma^2) of the Gamma limit.  An absolute check such as
  # lower + upper = 1 cannot see an error in a probability of 1e-19, so the
  # comparison is relative, at the working tolerance.  (expect_equal would
  # compare these in absolute terms, since they are below its tolerance.)
  ref <- c(1.6666854165921688e-19, 1.6666741665066692e-16)
  expect_lt(max(abs(pmhn(c(1e-9, 1e-8), 3, 1, -1000) / ref - 1)), 2e-8)
})

test_that("pmhn depends on beta only through sqrt(beta) q and z at strong tilts", {
  # The identity of section 10, on triples that the quadrature answers:
  # z = -1000 (past the length ceiling of the series), z = 2000 (the shifted
  # form for a large positive tilt) and z = -30 (series and quadrature both in
  # play).  At beta = 1e24 the mode of (2.5, beta, -1000 sqrt(beta)) is
  # 1.5e-15 and that of (10, beta, -30 sqrt(beta)) is 3e-13.  1e-8 is the
  # working tolerance: the series, where it answers, may use up to about half
  # of it, and which side of its cancellation guard a point falls on can
  # differ from one scale to another.
  betas <- c(1e-10, 1e10, 1e24)
  for (p in list(c(2.5, -1000), c(2, 2000), c(10, -30))) {
    a <- p[1]; z <- p[2]
    y <- qmhn(c(1e-6, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1 - 1e-6), a, 1, z)
    lower1 <- pmhn(y, a, 1, z)
    upper1 <- pmhn(y, a, 1, z, lower.tail = FALSE)
    for (b in betas) {
      info <- sprintf("alpha=%g z=%g beta=%g", a, z, b)
      lower <- pmhn(y / sqrt(b), a, b, z * sqrt(b))
      upper <- pmhn(y / sqrt(b), a, b, z * sqrt(b), lower.tail = FALSE)
      expect_lt(max(abs(lower - lower1)), 1e-8, label = info)
      expect_lt(max(abs(upper - upper1) / upper1), 1e-8, label = info)
    }
  }
  # (2, 1e26, 1e17) is z = 1e4, where the law is symmetric about its mode up
  # to O((alpha - 1) / z^3), so F(mode) = 1/2 to 1e-12 (integrate() in the
  # variable centred on the mode gives 0.4999999999992).  1e-8 is the working
  # tolerance.
  expect_equal(pmhn(mhn_mode(2, 1e26, 1e17), 2, 1e26, 1e17), 0.5,
               tolerance = 1e-8)
})

test_that("pmhn keeps the mass above a mode squeezed by a strong negative tilt", {
  # At (2.5, 1, -3e4) the mode is 5e-5 and the law is Gamma(2.5, 3e4) to a
  # relative 1e-8, so F is 1 to double precision from q = 0.01 upward: the
  # upper tail at q = 0.01 is below exp(-290).  The part of the range above
  # the mode must be sized by how fast the kernel falls there, over about
  # L m / (alpha - 1), and not by the Gaussian width sqrt(L / beta), which is
  # thousands of times wider and puts every quadrature node where the kernel
  # has already vanished.
  expect_equal(pmhn(c(0.01, 1, 6, 10, 1000), 2.5, 1, -3e4), rep(1, 5),
               tolerance = 1e-12)
  expect_equal(pmhn(1000, 3, 1e-4, -1000), 1, tolerance = 1e-12)
  q <- stats::qgamma(c(0.1, 0.5, 0.9), 2.5, rate = 3e4)
  expect_equal(pmhn(q, 2.5, 1, -3e4), c(0.1, 0.5, 0.9), tolerance = 2e-8)
  for (a in c(1.01, 1.5, 5, 20)) {
    for (g in c(-2.2e4, -1e5, -1e7)) {
      expect_equal(pmhn(100, a, 1, g), 1, tolerance = 1e-12,
                   info = sprintf("alpha=%g gamma=%g", a, g))
    }
  }
})

test_that("pmhn reaches one past the support when alpha < 1 and gamma < 0", {
  # Without an interior peak (Sun et al. 2023, Lemma 3d) the kernel falls from
  # the origin, and once beta t^2 - gamma t exceeds L the rest of the range
  # holds less than exp(-L) of the mass.  For q >= 1e3 at these tilts the
  # upper tail is below exp(-2e4), so F is 1 to double precision.  A range
  # running out to q itself leaves tanh_sinh sampling a region where the
  # integrand has vanished, and F then falls back below one.
  expect_equal(pmhn(c(4356, 5003, 1e8), 0.05, 1, -1000), rep(1, 3),
               tolerance = 1e-12)
  expect_equal(pmhn(c(42210, 2.357e10), 0.8, 1, -1e4), c(1, 1),
               tolerance = 1e-12)
  expect_equal(pmhn(2.47e8, 0.7, 1, -20), 1, tolerance = 1e-12)
  q <- 10^seq(3, 12, length.out = 200)
  for (a in c(0.05, 0.5, 0.8, 0.95)) {
    for (g in c(-1e3, -1e4, -1e5)) {
      v <- pmhn(q, a, 1, g)
      expect_lt(max(abs(v - 1)), 1e-12,
                label = sprintf("alpha=%g gamma=%g", a, g))
    }
  }
})

test_that("pmhn keeps the mass next to x below a Lemma 3c peak at a strong tilt", {
  # For alpha < 1 and a large positive tilt the kernel below its peak
  # (Sun et al. 2023, Lemma 3c) rises from a local minimum near the origin to
  # x across a range gamma / (2 beta) wide, while the mass lies within a few
  # units of x.  A single tanh_sinh range over the whole of it samples that
  # sliver too sparsely, so the range below x must be trimmed to where the
  # kernel has fallen below exp(-58) of its value at x.  The upper tail is
  # computed separately on a logarithmic axis, so the two tails must sum to
  # one; 1e-8 is the working tolerance.  The tilts stop short of about
  # 7.3e5, beyond which both tails are measured from the peak of the density
  # of log X (mhn_cdf_peak.h) and sum to one by construction.
  for (a in c(0.2, 0.7)) {
    for (z in c(6606.934, 87096.359, 138038.426, 147295, 7e5)) {
      q <- mhn_mode(a, 1, z) - c(2, 1.5, 0.906, 0.5, 0.1)
      lo <- pmhn(q, a, 1, z)
      up <- pmhn(q, a, 1, z, lower.tail = FALSE)
      expect_lt(max(abs(lo + up - 1)), 1e-8,
                label = sprintf("alpha=%g z=%g", a, z))
    }
  }
})

test_that("the lower tail is monotone and complementary to the upper near F = 1", {
  # For gamma < 0 the Lemma 1b series alternates, and near F = 1 its rounding
  # is of order 1e-8.  Above the split point F is one minus the directly
  # computed upper tail, so the two tails sum to one to rounding and F rises
  # wherever the upper tail falls.  References: the Lemma 1b series summed in
  # Rmpfr at 320 bits (400 terms; 600 give the same digits).
  x <- seq(2.9, 3.0, by = 0.001)
  lo <- pmhn(x, 5, 1, -5)
  up <- pmhn(x, 5, 1, -5, lower.tail = FALSE)
  expect_false(is.unsorted(lo))
  expect_lt(max(abs(lo + up - 1)), 1e-14)
  expect_equal(pmhn(c(2.96, 2.972, 2.973), 5, 1, -5),
               c(0.99999986036732003, 0.99999987589631752,
                 0.99999987711153869),
               tolerance = 1e-13)
  # The quantile near p = 1 depends on that: the density at 2.973 is about
  # 1.2e-6, so an error of 1e-8 in F moves the quantile by about 1e-2.
  p <- 1 - pmhn(2.973, 5, 1, -5, lower.tail = FALSE)
  expect_equal(qmhn(p, 5, 1, -5), 2.973, tolerance = 1e-8)
})

test_that("the log upper tail keeps its relative accuracy where it is near 0", {
  # At or below the median log P(X > q) = log(1 - F) is about -F, and it is
  # formed from log F as log1p(-F), which keeps the relative accuracy of F.
  # The logarithm of the upper-tail quadrature over [q, Inf), a value close
  # to 1 there, would carry the absolute error of that value into
  # log P(X > q), which at F = 1e-10 is off by up to 9 percent at these
  # points.  They cover the series (alpha < 1, and |z| <= 20 with
  # alpha <= 100), the quadrature (alpha = 150 and 7730), and the quadrature
  # with gamma^2/(4 beta) removed (z = 2000), at F from 1e-10 to 0.3.
  # References: log(1 - F) from the integral of the kernel on the axis
  # log(sqrt(beta) x), by Gauss-Legendre panels in 256-bit arithmetic.  1e-8
  # lies below the working tolerance, sqrt(eps) = 1.5e-8.  On the linear
  # scale P(X > q) stays the quadrature, whose absolute error is below 1e-10.
  cases <- list(
    list(2, 1, 1, 1.65e-5, -9.9717808618990636e-11),
    list(3, 1, 10, 0.92, -1.0244979760653198e-10),
    list(0.3, 1, -5, 6.2e-35, -1.0009139689368066e-10),
    list(50, 1, -15, 0.927, -1.0023144086084133e-10),
    list(150, 2, 30, c(8.2, 9.3),
         c(-9.1288056346165749e-11, -8.7461579183415104e-5)),
    list(2, 1, 2000, 995.5, -9.7855195607955021e-11),
    list(7730, 1.8e-4, 1.7576, c(7392, 7511, 7655),
         c(-9.5761133596843752e-11, -1.0099374942669675e-4,
           -0.35975912393828986)))
  for (cs in cases) {
    a <- cs[[1]]; b <- cs[[2]]; g <- cs[[3]]; q <- cs[[4]]
    info <- sprintf("alpha=%g beta=%g gamma=%g", a, b, g)
    got <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
    expect_lt(max(abs(got / cs[[5]] - 1)), 1e-8, label = info)
    expect_lt(max(abs(pmhn(q, a, b, g, lower.tail = FALSE) - exp(cs[[5]]))),
              1e-10, label = info)
    # The recycled-parameter path takes the same per-element route.
    expect_identical(pmhn(q, rep(a, length(q)), b, g, lower.tail = FALSE,
                          log.p = TRUE), got, label = info)
  }
})

test_that("the log upper tail below the median is accurate on either route", {
  # Below the median log P(X > q) is the logarithm of the upper-tail
  # quadrature wherever the rounding of the terms that quadrature is
  # assembled from, eps times about |log Psi| + (alpha / 2) |log beta| +
  # max(alpha, z) (1 + |log v*|), with gamma^2 / (4 beta) removed from
  # log Psi for z > 1000, is below 1e-9 of log P(X > q) itself, and
  # log1p(-F) elsewhere.  At (2e4, 1e300, 5e151) those terms reach 7e6, and
  # the logarithm of the quadrature is off by 2.6e-8 of log P(X > q) at
  # F = 0.05; F itself is accurate there to about 1e-9 relative (see
  # Range in ?pmhn), which log1p(-F) keeps, and 3e-9 allows for that.  At
  # (2.5, 1, 1.5) and (0.3, 1, -300) the terms are of order one and the
  # logarithm of the quadrature is within 3e-14 at F = 0.01 to 0.3, where
  # it is used; 1e-12 allows for it.  References as in the test above, in
  # 160-bit arithmetic.
  cases <- list(
    list(2e4, 1e300, 5e151, c(1.124055e-148, 1.129986e-148),
         c(-0.051291561297172605, -0.35669413691881596), 3e-9),
    list(2.5, 1, 1.5, c(0.550397, 1.08975),
         c(-0.051293319185341725, -0.35667571199364523), 1e-12),
    list(0.3, 1, -300, c(5.00734e-10, 4.24216e-05),
         c(-0.01005033795650298, -0.35667502247568539), 1e-12))
  for (cs in cases) {
    a <- cs[[1]]; b <- cs[[2]]; g <- cs[[3]]; q <- cs[[4]]
    info <- sprintf("alpha=%g beta=%g gamma=%g", a, b, g)
    got <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
    expect_lt(max(abs(got / cs[[5]] - 1)), cs[[6]], label = info)
    expect_identical(pmhn(q, rep(a, length(q)), b, g, lower.tail = FALSE,
                          log.p = TRUE), got, label = info)
  }
})

test_that("the log upper tail keeps its accuracy where the series is costly", {
  # Below the median, where the logarithm of the upper-tail quadrature cannot
  # resolve F to 1e-9 of itself, log P(X > q) is log1p(-F).  For alpha < 1
  # past |z| = 20 the Lemma 1b series costs milliseconds, and F is taken
  # there from the lower-tail quadrature wherever that agrees with the upper
  # tail to within the error bound of the latter.  At (0.3, 1, 300) and
  # (0.01, 1, 500), where log Psi is about z^2 / 4, this is so below F of
  # about 0.06 and 0.17, and at (0.3, 1, -300) below about 5e-6; the points
  # are at F = 1e-3 and 0.015, 1e-4 and 0.04, and 1e-6 and 2e-6.  The
  # quadrature gives F to about 1e-11 of itself there, which log1p(-F)
  # keeps; 1e-10 allows for it.  References as in the test above, in
  # 256-bit arithmetic.
  cases <- list(
    list(0.3, 1, 300, c(147.8125, 148.4632),
         c(-0.0010003801618320024, -0.01511532580718665)),
    list(0.01, 1, 500, c(247.3683, 248.7601),
         c(-0.00010002308009331027, -0.040822875216645363)),
    list(0.3, 1, -300, c(2.324199e-23, 2.342646e-22),
         c(-1.0000004404356619e-06, -2.0000019332766965e-06)))
  for (cs in cases) {
    a <- cs[[1]]; b <- cs[[2]]; g <- cs[[3]]; q <- cs[[4]]
    info <- sprintf("alpha=%g beta=%g gamma=%g", a, b, g)
    got <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
    expect_lt(max(abs(got / cs[[5]] - 1)), 1e-10, label = info)
    expect_identical(pmhn(q, rep(a, length(q)), b, g, lower.tail = FALSE,
                          log.p = TRUE), got, label = info)
  }
})

test_that("pmhn keeps the lower tail at a very large shape", {
  # At alpha = 1e8 the mass is about one unit wide next to the mode
  # m = 7071, and the law is measured from the peak of the density of log X
  # (mhn_cdf_peak.h), so no term of the size of alpha log m = 9e8 is formed.
  # References: the integral of exp((alpha-1) (log1p(s/m) - s/m) - s^2), the
  # kernel measured from its peak m = (1 + sqrt(1 + 8 (alpha - 1))) / 4 at
  # beta = gamma = 1, over s from -15 to the exact offset of each x from m,
  # against [-15, 15], by Romberg integration in Rmpfr at 160 bits
  # (integrateR, rel.tol = 1e-25).  Both tails agree with them to about
  # 1e-12, and 1e-10 allows for that.
  a <- 1e8
  ref <- c(0.022748223028665679, 0.27424738441459101, 0.49999059733921408,
           0.72573408068790435, 0.97724414154069583)
  m <- (1 + sqrt(1 + 8 * (a - 1))) / 4
  x <- m + c(-1, -0.3, 0, 0.3, 1)
  expect_lt(max(abs(pmhn(x, a, 1, 1) / ref - 1)), 1e-10)
  expect_lt(max(abs(pmhn(x, a, 1, 1, lower.tail = FALSE) / (1 - ref) - 1)),
            1e-10)
})

test_that("the log lower tail stays finite where F rounds to 1 or underflows", {
  # log F is taken on the log scale throughout, not as log of F: above the
  # peak of the density of log X as log1p(-Q) of the upper tail Q, which is
  # -Q when Q is below 1e-16 and F itself rounds to 1, and below it from the
  # log-scale series or quadrature, which stay finite where F is below the
  # smallest double.  References: Q by stats::integrate of the kernel
  # x^(alpha-1) exp(-beta x^2 + gamma x) over [q, Inf) and [0, Inf) at
  # rel.tol = 1e-13; at q = 1e-200, F(q) = q^2 / (2 Z) to a relative 1e-199,
  # with Z = 25.297364396949772 (see test-qmhn.R).  The quadrature tolerance
  # of the upper tail is 1.5e-8 relative, and 1e-7 allows for it; the last
  # value comes from the closed form taken next to the origin and is good to
  # rounding.  The first two references are far below any tolerance, so the
  # comparison is made in relative terms explicitly (expect_equal would fall
  # back to an absolute one and accept 0).
  cases <- list(list(q = 10, a = 2.5, g = 1, ref = -9.69801078473816e-40),
                list(q = 6, a = 5, g = -5, ref = -5.2109227399623e-25),
                list(q = 1e-200, a = 2, g = 3,
                     ref = 2 * log(1e-200) - log(2 * 25.297364396949772)))
  for (cs in cases) {
    info <- sprintf("q=%g alpha=%g gamma=%g", cs$q, cs$a, cs$g)
    tol <- if (cs$q < 1) 1e-12 else 1e-7
    got <- pmhn(cs$q, cs$a, 1, cs$g, log.p = TRUE)
    expect_lt(abs(got / cs$ref - 1), tol, label = info)
    # The recycled-parameter path takes the same per-element route.
    got_vec <- pmhn(cs$q, c(cs$a, cs$a), 1, cs$g, log.p = TRUE)
    expect_identical(got_vec, rep(got, 2), label = info)
  }
})

test_that("pmhn gives NA on the log scale wherever it gives NA on the linear", {
  # Where the distribution function cannot be evaluated pmhn returns NA, and
  # log.p = TRUE must say the same: log F is limited to at most 0 to absorb
  # rounding, and a limit taken as min(0, log F) would turn a NaN into
  # log F = 0, that is F = 1.  The triples include a shape of 1e100 and
  # extreme scales, whose normalising constant is out of reach, and points
  # below the mode, where log F is far below 0.
  cases <- list(list(q = 3.535534e49 * c(1, 2, 4), a = 1e100, b = 1, g = -10),
                list(q = 3.535534e49 * c(1, 2, 4), a = 1e100, b = 1, g = 1e4),
                list(q = 2.5e9, a = 1e10, b = 1, g = 1e10),
                list(q = 0.25, a = 1e10, b = 1e300, g = 1e300))
  for (cs in cases) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", cs$a, cs$b, cs$g)
    F <- pmhn(cs$q, cs$a, cs$b, cs$g)
    log_F <- pmhn(cs$q, cs$a, cs$b, cs$g, log.p = TRUE)
    expect_identical(is.na(log_F), is.na(F), label = info)
    # The recycled-parameter path shares the per-element evaluation.
    log_F_vec <- pmhn(cs$q, c(cs$a, cs$a), cs$b, cs$g, log.p = TRUE)
    expect_identical(is.na(log_F_vec), is.na(rep_len(F, length(log_F_vec))),
                     label = info)
  }
})

test_that("pmhn at a huge shape and tilt is its normal limit", {
  # At alpha = 1e10, beta = 1, gamma = 1e10 the law is normal to far better
  # than 1e-12: its mode, (gamma + sqrt(gamma^2 + 8 (alpha - 1))) / 4, is
  # 5e9 + 1 - 3e-10 to within 1e-19, the curvature of the log density there
  # is 2 + (alpha - 1) / m^2, and its skewness is of order 1e-19.  The
  # points q are formed in double precision, 9.5e-7 apart near m, and
  # q - (5e9 + 1) is exact.  Assembled from log Psi and the log kernel, each
  # about alpha log m = 2.2e11, the tails would carry 5e-5 of rounding;
  # measured from the peak they are good to about 1e-15.  1e-12 is far
  # inside the working tolerance.
  m <- 5e9 + 1
  sd <- 1 / sqrt(2 + 1e10 / m^2)
  q <- m + c(-3, -1.4142, 0, 1, 3) * sd
  t <- ((q - m) + 3e-10) / sd
  expect_lt(max(abs(pmhn(q, 1e10, 1, 1e10) - pnorm(t))), 1e-12)
  expect_lt(max(abs(pmhn(q, 1e10, 1, 1e10, lower.tail = FALSE) -
                      pnorm(t, lower.tail = FALSE))), 1e-12)
})

test_that("d, p and q are NA with a warning where the law cannot be formed", {
  # At alpha = 1e300, the subnormal beta = 1e-323 and gamma = -1e-160 the
  # peak of the law, sqrt(alpha / (2 beta)) about 2e311, lies beyond the
  # largest double, where the frame the large shapes are measured in cannot
  # be set up.  Rather than a value it cannot vouch for, each function
  # returns NA (not the NaN of an invalid argument), with a warning.
  expect_warning(F <- pmhn(c(1, 1e300), 1e300, 1e-323, -1e-160),
                 "could not be evaluated")
  expect_true(all(is.na(F) & !is.nan(F)))
  expect_warning(F <- pmhn(c(1, 1e300), c(1e300, 1e300), 1e-323, -1e-160,
                           lower.tail = FALSE, log.p = TRUE),
                 "could not be evaluated")
  expect_true(all(is.na(F) & !is.nan(F)))
  expect_warning(d <- dmhn(c(1, 1e300), 1e300, 1e-323, -1e-160),
                 "could not be evaluated")
  expect_true(all(is.na(d) & !is.nan(d)))
  expect_warning(d <- dmhn(c(1, 1e300), c(1e300, 1e300), 1e-323, -1e-160),
                 "could not be evaluated")
  expect_true(all(is.na(d) & !is.nan(d)))
  # An NA or NaN passed in is passed through silently, as elsewhere, and so
  # are 0 and Inf, where the density is 0 whatever the law.
  expect_silent(F <- pmhn(c(NA, NaN), 1e300, 1e-323, -1e-160))
  expect_identical(is.nan(F), c(FALSE, TRUE))
  expect_silent(F <- pmhn(c(NA, NaN), c(1e300, 1e300), 1e-323, -1e-160))
  expect_identical(is.nan(F), c(FALSE, TRUE))
  expect_silent(d <- dmhn(c(NA, NaN, 0, Inf), c(1e300, 1e300), 1e-323,
                          -1e-160))
  expect_identical(d, c(NA, NaN, 0, 0))
})

test_that("the series guard declines series values beyond the working tolerance", {
  # The rounding of each series term is proportional to the magnitudes of the
  # logarithms it is formed from, and at a large alpha lgamma((alpha + i)/2)
  # dominates them.  At alpha = 2000 the series is off by 5e-8 and 1.2e-7
  # relative at these points, 3 to 8 times the working tolerance.
  # References: the Lemma 1b series summed in Rmpfr at 320 bits (200 terms;
  # 300 give the same digits).
  expect_equal(pmhn(31.5666, 2000, 1, -0.2036507), 0.49996919763182685,
               tolerance = 1e-8)
  expect_equal(pmhn(30.41166598652667, 2000, 1, -0.2036507),
               0.010000000869548123, tolerance = 1e-8)
})

test_that("dmhn, pmhn and qmhn keep the attributes of their first argument", {
  # As in dnorm, pnorm and qnorm, the names, dim and dimnames of x, q or p
  # carry over when the result has the same length.  A classed first argument
  # (a Date, say) is used through as.numeric and gives a plain vector.
  # Unlike pnorm, attributes of a longer parameter vector are not carried
  # over; the last line checks that documented departure.
  m <- matrix(c(0.5, 1, 1.5, 2), 2, dimnames = list(c("a", "b"), NULL))
  for (f in list(function(x) dmhn(x, 2, 1, 1), function(x) pmhn(x, 2, 1, 1),
                 function(x) pmhn(x, 2, 1, -1, lower.tail = FALSE),
                 function(x) qmhn(x / 3, 2, 1, 1))) {
    out <- f(m)
    expect_equal(dim(out), c(2L, 2L))
    expect_equal(dimnames(out), dimnames(m))
    expect_equal(as.vector(out), f(as.vector(m)))
    expect_named(f(c(u = 0.5, v = 1)), c("u", "v"))
  }
  expect_null(attributes(dmhn(as.Date(c("1970-01-04", "1970-01-06")), 2, 1, 1)))
  expect_null(names(pmhn(0.5, alpha = c(a = 2, b = 3), beta = 1, gamma = 1)))
})

# ============================================================
# 12. Tilts past the overflow of z^2
# ============================================================

test_that("pmhn follows the Gamma limit at a negative tilt past the overflow of z^2", {
  # From |z| = 1.34e154 on z^2 overflows, and the anchors of the quadrature
  # -- the peak of the density of log X, the root where exp(gamma x - beta
  # x^2) has fallen by L, and the head of the range -- must be found without
  # it.  The law is Gamma(alpha, rate |gamma|) to within relative corrections
  # of order (alpha + 1) / z^2, below 1e-300, so both tails at q are those of
  # pgamma at |z| q; forming |z| q rounds it by 1.1e-16 relative, which moves
  # either tail by less than 1e-14 here.  1e-8 is the working tolerance.
  for (a in c(0.5, 2, 10)) {
    for (z in c(-1e155, -1e200, -1e300)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      q <- qgamma(c(0.001, 0.1, 0.5, 0.9, 0.999), a) / abs(z)
      y <- abs(z) * q
      expect_equal(pmhn(q, a, 1, z), pgamma(y, a), tolerance = 1e-8,
                   info = info)
      expect_equal(pmhn(q, a, 1, z, lower.tail = FALSE),
                   pgamma(y, a, lower.tail = FALSE), tolerance = 1e-8,
                   info = info)
      expect_equal(pmhn(q, a, 1, z, log.p = TRUE),
                   pgamma(y, a, log.p = TRUE), tolerance = 1e-8, info = info)
    }
  }
})

test_that("pmhn at a positive tilt past the overflow of z^2 is its normal limit", {
  # With beta = 1 the law is normal to far below double precision, with mean
  # m = z/2 + delta, delta = 2 (alpha - 1) / (z + sqrt(z^2 + 8 (alpha - 1))),
  # below 1e-154 here, and variance 1 / (2 + (alpha - 1)/m^2) = 1/2 to a
  # relative 1e-300.  It is far narrower than the spacing of the doubles at
  # its mode, and both tails and the normalising constant are measured from
  # the peak in exact offsets (mhn_cdf_peak.h), here in closed form.  z/2 is
  # exact and lies within 1e-154 of m, so each tail there is 1/2.
  # x = z/2 (1 +- 1e-10) lies d = z/2 1e-10 from z/2, exactly as R forms the
  # difference, and t = d sqrt(2) standard deviations from m: the near tail
  # is Phi(-t) to within a relative 1e-150, which underflows, and its
  # logarithm is R's pnorm(-t, log.p = TRUE), -2.5e289 at z = 1e155; the far
  # tail is 1.  1e-8 is the working tolerance; the tails at z/2 are taken to
  # about 1e-15.
  for (a in c(0.5, 2, 10)) {
    for (z in c(1e155, 1e200, 1e300, .Machine$double.xmax)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      expect_equal(pmhn(z / 2, a, 1, z), 0.5, tolerance = 1e-8, info = info)
      expect_equal(pmhn(z / 2, a, 1, z, lower.tail = FALSE), 0.5,
                   tolerance = 1e-8, info = info)
      x <- z / 2 * (1 + c(-1e-10, 1e-10))
      expect_equal(pmhn(x, a, 1, z), c(0, 1), tolerance = 1e-12, info = info)
      expect_equal(pmhn(x, a, 1, z, lower.tail = FALSE), c(1, 0),
                   tolerance = 1e-12, info = info)
    }
    z <- 1e155
    x <- z / 2 * (1 - 1e-10)
    t <- (z / 2 - x) * sqrt(2)
    expect_equal(pmhn(x, a, 1, z, log.p = TRUE), pnorm(-t, log.p = TRUE),
                 tolerance = 1e-8, info = sprintf("alpha=%g", a))
    x <- z / 2 * (1 + 1e-10)
    t <- (x - z / 2) * sqrt(2)
    expect_equal(pmhn(x, a, 1, z, lower.tail = FALSE, log.p = TRUE),
                 pnorm(-t, log.p = TRUE), tolerance = 1e-8,
                 info = sprintf("alpha=%g", a))
  }
})

test_that("the upper tail below a peak narrower than the doubles returns one", {
  # At each of these points the law is narrower than the spacing of the
  # doubles at its peak, and q lies far below it (by 1e84 standard deviations
  # in the first), so P(X > q) is 1 to double precision.  Integrated on the
  # log axis about a peak rounded to the doubles, the slope there, a
  # rounding residue of about eps z^2 / 2, or the curvature, overflows when
  # squared, and a panel loop stepping by a reach formed from that square
  # need never end.  (0.49e86, 2, 1e-100, 1e36) is the first at z = 1e86;
  # (0.99 sqrt(5e169), 1e170, 1, 1) has its peak near sqrt((alpha - 1)/2) and
  # a relative width of about 1/sqrt(2 alpha); at (1, 1e-300, 1e-300, 1e-10)
  # the mass lies near 5e289.  The median at (2, 1, 1e86) is its mode, 5e85,
  # to a relative 1e-86, and the root finder returns it to within a spacing
  # of the doubles.  The tails are taken to about 1e-15 (see the previous
  # test), far inside the working tolerance of 1e-8.
  expect_equal(pmhn(4e84, 2, 1, 1e85, lower.tail = FALSE), 1, tolerance = 1e-12)
  expect_equal(pmhn(0.49e86, 2, 1e-100, 1e36, lower.tail = FALSE), 1,
               tolerance = 1e-12)
  expect_equal(pmhn(0.99 * sqrt(5e169), 1e170, 1, 1, lower.tail = FALSE), 1,
               tolerance = 1e-12)
  expect_equal(pmhn(1, 1e-300, 1e-300, 1e-10, lower.tail = FALSE), 1,
               tolerance = 1e-12)
  expect_equal(qmhn(0.5, 2, 1, 1e86, lower.tail = FALSE), 5e85,
               tolerance = 1e-12)
})

test_that("pmhn depends on beta only through sqrt(beta) q and z where alpha beta overflows", {
  # At beta = 1e308, 8 beta (alpha - 1) and gamma^2 both pass the largest
  # double while z = gamma / sqrt(beta) is 10 or -10, an ordinary tilt; the
  # peak and every quadrature range must be formed without them.  The
  # identity of section 10 then ties the two scales.  1e-8 is the working
  # tolerance.  sqrt(1e308) carries a relative rounding of about 1e-16, which
  # moves q / sqrt(beta) by that much, the tails by below 1e-12; the terms
  # (alpha/2) log(beta) and log Psi that the normalising constant is formed
  # from are about 3.5e5 at alpha = 1000, and their rounding costs below 1e-10.
  b <- 1e308
  sb <- sqrt(b)
  for (p in list(c(2, -10), c(1000, -10), c(2, 10), c(1000, 10))) {
    a <- p[1]; z <- p[2]
    info <- sprintf("alpha=%g z=%g", a, z)
    y <- qmhn(c(0.01, 0.3, 0.5, 0.7, 0.99), a, 1, z)
    expect_lt(max(abs(pmhn(y / sb, a, b, z * sb) - pmhn(y, a, 1, z))), 1e-8,
              label = info)
    expect_lt(max(abs(pmhn(y / sb, a, b, z * sb, lower.tail = FALSE) -
                        pmhn(y, a, 1, z, lower.tail = FALSE))), 1e-8,
              label = info)
    expect_equal(qmhn(c(0.01, 0.5, 0.99), a, b, z * sb) * sb,
                 y[c(1, 3, 5)], tolerance = 1e-8, info = info)
  }
})

# ============================================================
# 13. Laws narrower than the spacing of the doubles
# ============================================================

test_that("a narrow tail where the curvature of the log kernel is negative is exponential", {
  # For alpha < 1 the curvature of the log kernel, 2 beta + (alpha - 1)/x^2,
  # is negative below x = sqrt((1 - alpha) / (2 beta)).  At these tilts the
  # kernel is far narrower than the spacing of the doubles at each x, each
  # tail is taken in closed form, and the slope there carries the kernel
  # down so steeply that the curvature term is negligible.  Under the tilt
  # gamma = -1e20 the law is Gamma(0.5, rate 1e20) to within relative
  # corrections of order (alpha + 1) / z^2 = 1.5e-40, so P(X <= x) is 1 at
  # these x and log P(X > x) is pgamma's at 1e20 x (the term beta x^2 it
  # leaves out is at most 0.25, a relative 5e-21).  Under the tilt 1e20 the
  # law is normal about 5e19 with standard deviation 1/sqrt(2): P(X <= x) is
  # 0 and its logarithm is -(5e19 - x)^2 + O(log z), -2.5e39 to a relative
  # 2e-20.  The same tail at beta = 1e-300 is the one at beta = 1 with x
  # scaled by sqrt(beta) (section 10): at sqrt(beta) x = 0.5, z = -1e17 its
  # logarithm is -5e16 - 0.5 log(5e16) - log(sqrt(pi)), -5e16 to a relative
  # 4e-16.  At the subnormal beta = 1e-320, where half the curvature is
  # beta itself, the upper tail at sqrt(beta) x = 1e10 with z = 1 has
  # logarithm -1e20 + 1e10 + O(log 1e10), -1e20 to a relative 1e-10.  1e-8
  # is the working tolerance.
  x <- c(0.01, 0.1, 0.5)
  expect_identical(pmhn(x, 0.5, 1, -1e20), c(1, 1, 1))
  expect_equal(pmhn(x, 0.5, 1, -1e20, lower.tail = FALSE, log.p = TRUE),
               pgamma(1e20 * x, 0.5, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-8)
  x <- c(0.01, 0.1, 0.4)
  expect_identical(pmhn(x, 0.5, 1, 1e20), c(0, 0, 0))
  expect_equal(pmhn(x, 0.5, 1, 1e20, log.p = TRUE), rep(-2.5e39, 3),
               tolerance = 1e-8)
  sb <- sqrt(1e-300)
  expect_equal(pmhn(0.5 / sb, 0.5, 1e-300, -1e17 * sb, lower.tail = FALSE,
                    log.p = TRUE), -5e16, tolerance = 1e-8)
  sb <- sqrt(1e-320)
  expect_equal(pmhn(1e10 / sb, 2, 1e-320, sb, lower.tail = FALSE,
                    log.p = TRUE), -1e20, tolerance = 1e-8)
  # (1e-300, 1e100, -1e250) has z = -1e200 and is Gamma(1e-300, rate 1e250)
  # to within a relative 1e-400, and at x = 1e-51 the curvature is
  # 1e100 - 0.5e102 < 0.  P(X <= 1e-51) = P(G <= 1e199) = 1, and every
  # quantile, (p Gamma(1 + alpha))^(1 / alpha) / 1e250, underflows to 0.
  expect_identical(pmhn(1e-51, 1e-300, 1e100, -1e250), 1)
  expect_identical(qmhn(c(0.1, 0.5, 0.9), 1e-300, 1e100, -1e250), c(0, 0, 0))
})

test_that("a law narrower than the doubles at a beta near the bottom of the range", {
  # With z = gamma / sqrt(beta) = 1e18 or 1e100 the law is normal about its
  # mode m with standard deviation 1/sqrt(2 beta), a relative sqrt(2) / z of
  # m or less, far below the spacing of the doubles there.  m (1 - 1e-10)
  # and m (1 + 1e-10) lie 7e7 standard deviations or more away, so the
  # tails there are 0 and 1; at m itself the two tails must add to 1.  The
  # reach of the kernel from its peak, sqrt(L / beta), would overflow if
  # formed from L / beta, from beta = 3e-307 down, and so would pi / beta,
  # from 1.7e-308 down; both ranges are exercised here, down to the
  # subnormal 1e-320.  The two tails at m, one taken directly and the other
  # as its complement, add to 1 to rounding, hence the tolerance.
  for (b in c(1e-306, 1e-307, 1e-308, 1e-320)) {
    for (z in c(1e18, 1e100)) {
      info <- sprintf("beta=%g z=%g", b, z)
      g <- z * sqrt(b)
      m <- mhn_mode(2, b, g)
      x <- m * c(1 - 1e-10, 1 + 1e-10)
      expect_equal(pmhn(x, 2, b, g), c(0, 1), tolerance = 1e-12, info = info)
      expect_equal(pmhn(x, 2, b, g, lower.tail = FALSE), c(1, 0),
                   tolerance = 1e-12, info = info)
      expect_equal(pmhn(m, 2, b, g) + pmhn(m, 2, b, g, lower.tail = FALSE), 1,
                   tolerance = 1e-12, info = info)
    }
  }
})

test_that("a law narrower than the doubles at a huge shape keeps its step in place", {
  # At these shapes the law's standard deviation is a relative 3.2e-17 or
  # less of its mode m (1 / sqrt(alpha) for the Gamma-like law under the
  # tilt -1e100, 1 / (2 m) with m = sqrt(alpha / 2) otherwise), so
  # m (1 - 1e-10) and m (1 + 1e-10) lie 3e6 standard deviations or more below
  # and above it, where the tails are 0 and 1.  The log kernel and log Psi
  # are each about 1e35 or more in size here, and assembled from them the
  # tails would carry a rounding of 1e19 nats or more; they are measured from
  # the peak instead (mhn_cdf_peak.h).  At (1e200, 1e300, -1e300) the mode is
  # 1e-100 and (alpha - 1) / m^2 there overflows; the tail at the mode
  # depends on where the law's peak lies within the spacing of the doubles
  # there, and must be a number, the two tails adding to 1.
  for (p in list(c(1e33, 1, -1e100), c(1e35, 1, -1), c(1e40, 1, 10))) {
    info <- paste(p, collapse = " ")
    m <- mhn_mode(p[1], p[2], p[3])
    x <- m * c(1 - 1e-10, 1 + 1e-10)
    expect_identical(pmhn(x, p[1], p[2], p[3]), c(0, 1), info = info)
    expect_identical(pmhn(x, p[1], p[2], p[3], lower.tail = FALSE), c(1, 0),
                     info = info)
  }
  m <- mhn_mode(1e200, 1e300, -1e300)
  f <- pmhn(m, 1e200, 1e300, -1e300)
  expect_true(is.finite(f) && f >= 0 && f <= 1)
  expect_equal(f + pmhn(m, 1e200, 1e300, -1e300, lower.tail = FALSE), 1,
               tolerance = 1e-12)
})

test_that("a law narrower than the doubles is placed by where its peak lies between them", {
  skip_if_not_installed("Rmpfr")
  # At a shape of 1e35 the law is normal about its peak v* to within
  # relative corrections of order 1 / sqrt(alpha) = 3e-18, with standard
  # deviation sd = v* / sqrt(alpha + 2 beta v*^2), a relative 2.2e-18 of v*
  # at gamma = 0, while the doubles next to v* are a relative 1.1e-16 to
  # 2.2e-16 apart.  At each of the three doubles around v*, then,
  # F = pnorm((x - v*) / sd) and log f = dnorm((x - v*) / sd, log = TRUE) -
  # log(sd), and at these beta the middle one lies 20.6 and 1.9 standard
  # deviations from v*.  The closed form P(alpha / 2, beta x^2) of gamma = 0
  # cannot give them: beta x^2 is rounded to a relative 3.3e-16, 75 standard
  # deviations of the Gamma law of beta X^2.  Under the tilt -1e300 at
  # (1e200, 1e300) the law is narrower still, sd = 1e-200 about
  # v* = 1e-100 (1 - 2e-100), and the two doubles around v* lie 6.6e83
  # standard deviations below and 6.1e83 above it, where F is 0 and 1; there
  # the offset of v* from a double v next to it is formed relative to v, as
  # Q'(v) / (2 v) = 2 beta - gamma / (2 v) passes the largest double.  v*
  # and the offsets are formed in Rmpfr, at 1024 bits for the second law (sd
  # is a relative 1e-100 of v*).
  around <- function(v) v + c(-1, 0, 1) * 2^(floor(log2(v)) - 52)
  for (ab in list(c(1e35, 1e100), c(1e35, 1e300))) {
    a <- ab[1]
    b <- ab[2]
    info <- sprintf("alpha=%g beta=%g", a, b)
    V <- sqrt(Rmpfr::mpfr(a, 256) / (2 * Rmpfr::mpfr(b, 256)))
    x <- around(Rmpfr::asNumeric(V))
    sd <- Rmpfr::asNumeric(V) / sqrt(2 * a)
    zx <- Rmpfr::asNumeric((Rmpfr::mpfr(x, 256) - V) / sd)
    expect_lt(min(abs(zx)), 25, label = info)
    expect_lt(max(abs(pmhn(x, a, b, 0) - pnorm(zx))), 1e-9, label = info)
    expect_lt(max(abs(pmhn(x, a, b, 0, lower.tail = FALSE) -
                        pnorm(zx, lower.tail = FALSE))), 1e-9, label = info)
    expect_lt(max(abs(dmhn(x, a, b, 0, log = TRUE) -
                        (dnorm(zx, log = TRUE) - log(sd)))), 1e-9,
              label = info)
  }
  A <- Rmpfr::mpfr(1e200, 1024)
  B <- Rmpfr::mpfr(1e300, 1024)
  G <- Rmpfr::mpfr(-1e300, 1024)
  V <- 2 * A / (sqrt(G^2 + 8 * A * B) - G)
  x <- around(Rmpfr::asNumeric(V))
  zx <- Rmpfr::asNumeric((Rmpfr::mpfr(x, 1024) - V) /
                           (V / sqrt(A + 2 * B * V^2)))
  side <- as.numeric(zx > 0)
  expect_true(all(abs(zx) > 1e80) && any(zx > 0) && any(zx < 0))
  expect_identical(pmhn(x, 1e200, 1e300, -1e300), side)
  expect_identical(pmhn(x, 1e200, 1e300, -1e300, lower.tail = FALSE),
                   1 - side)
  expect_true(all(dmhn(x, 1e200, 1e300, -1e300, log = TRUE) < -1e100))
})

test_that("pmhn and qmhn follow the limit laws where gamma / sqrt(beta) overflows", {
  # z = gamma / sqrt(beta) passes the largest double once |gamma| exceeds it
  # times sqrt(beta), which needs beta < 1, while gamma and
  # log |z| = log |gamma| - log(beta) / 2 are still in range.  Under such a
  # negative tilt the law is Gamma(alpha, rate |gamma|) to within a relative
  # alpha (alpha + 1) / z^2, 0 in double precision, and the tails and
  # quantiles are pgamma's and qgamma's, the product |gamma| x rounded as in
  # section 12.  Under such a positive tilt the law is normal with variance
  # 1 / (2 beta) about its mode, gamma / (2 beta) to within (alpha - 1) /
  # gamma.  For (2, 1e-100, DBL_MAX) that is beyond the largest double, so
  # every finite x lies below the law.  At beta = 1/2 it is gamma itself to
  # within 1e-307, and gamma (1 -+ 1e-10) lie 1e298 standard deviations away.
  # 1e-8 is the working tolerance; the tails near the mode are taken to
  # about 1e-15 (section 12).
  X <- .Machine$double.xmax
  p <- c(0.001, 0.1, 0.5, 0.9, 0.999)
  for (a in c(0.5, 2, 10)) {
    for (bg in list(c(1e-100, -X), c(1e-200, -1e250))) {
      b <- bg[1]
      g <- bg[2]
      info <- sprintf("alpha=%g beta=%g gamma=%g", a, b, g)
      x <- qgamma(p, a) / abs(g)
      y <- abs(g) * x
      expect_equal(pmhn(x, a, b, g), pgamma(y, a), tolerance = 1e-8,
                   info = info)
      expect_equal(pmhn(x, a, b, g, lower.tail = FALSE),
                   pgamma(y, a, lower.tail = FALSE), tolerance = 1e-8,
                   info = info)
      x <- 100 * x
      expect_equal(pmhn(x, a, b, g, lower.tail = FALSE, log.p = TRUE),
                   pgamma(abs(g) * x, a, lower.tail = FALSE, log.p = TRUE),
                   tolerance = 1e-8, info = info)
      expect_equal(qmhn(p, a, b, g) * abs(g), qgamma(p, a), tolerance = 1e-8,
                   info = info)
    }
    # The upper tail there is 1, the complement of a lower tail that is 0 in
    # double precision.
    x <- c(1, 1e300, X)
    expect_identical(pmhn(x, a, 1e-100, X), c(0, 0, 0))
    expect_equal(pmhn(x, a, 1e-100, X, lower.tail = FALSE), c(1, 1, 1),
                 tolerance = 1e-10)
    for (g in c(1.3e308, X)) {
      info <- sprintf("alpha=%g gamma=%g", a, g)
      x <- g * c(1 - 1e-10, 1 + 1e-10)
      expect_equal(pmhn(x, a, 0.5, g), c(0, 1), tolerance = 1e-12,
                   info = info)
      expect_equal(pmhn(x, a, 0.5, g, lower.tail = FALSE), c(1, 0),
                   tolerance = 1e-12, info = info)
      expect_equal(pmhn(g, a, 0.5, g) + pmhn(g, a, 0.5, g, lower.tail = FALSE),
                   1, tolerance = 1e-12, info = info)
    }
  }
})

# ============================================================
# 14. Large tilts and shapes, measured from the peak
# ============================================================

# The law of MHN(alpha, beta, gamma), computed independently of the package.
# The peak v of x g(x), g the kernel, is the positive root of
# 2 beta x^2 - gamma x - alpha = 0; it is found in Rmpfr at 256 bits from the
# doubles passed, and offsets r = x - v are formed there and then rounded.
# The log kernel relative to that peak,
#   K(r) = log g(v + r) - log g(v) = (alpha - 1) (log1p(rho) - rho) - rho
#          - beta r^2,   rho = r / v,
# (by 2 beta v - gamma = alpha / v) is integrated by stats::integrate in
# units of sd = 1 / sqrt(2 beta + alpha / v^2), to a relative 1e-13, over
# 60 sd on either side of the point asked for, beyond which it has fallen
# by far more than the double range.  log1p(u) - u is summed from its series
# below |u| = 1e-3, where the difference would cancel.  Returns the lower
# tail, the upper tail and the log density at x, and v and sd.
peak_law_reference <- function(a, b, g) {
  A <- Rmpfr::mpfr(a, 256)
  B <- Rmpfr::mpfr(b, 256)
  G <- Rmpfr::mpfr(g, 256)
  D <- sqrt(G^2 + 8 * B * A)
  V <- if (g >= 0) (G + D) / (4 * B) else 2 * A / (D - G)
  v <- Rmpfr::asNumeric(V)
  sd <- 1 / sqrt(2 * b + (a / v) / v)
  lmu <- function(u) {
    out <- log1p(u) - u
    s <- abs(u) < 1e-3
    w <- u[s]
    out[s] <- -w^2 / 2 + w^3 / 3 - w^4 / 4 + w^5 / 5 - w^6 / 6 + w^7 / 7
    out
  }
  K <- function(y) {
    r <- y * sd
    rho <- r / v
    (a - 1) * lmu(rho) - rho - b * r^2
  }
  lo <- -v / sd
  I <- function(l, u) {
    if (u <= l) return(0)
    stats::integrate(function(y) exp(K(y)), l, u, rel.tol = 1e-13,
                     abs.tol = 0, subdivisions = 2000L)$value
  }
  tot <- I(max(lo, -60), 0) + I(0, 60)
  y_of <- function(x) Rmpfr::asNumeric((Rmpfr::mpfr(x, 256) - V) / sd)
  left <- function(y) I(max(lo, y - 60), y)
  right <- function(y) I(y, y + 60)
  list(
    lower = function(x) vapply(x, function(xx) {
      y <- y_of(xx)
      if (y <= 0) left(y) / tot else 1 - right(y) / tot
    }, 0),
    upper = function(x) vapply(x, function(xx) {
      y <- y_of(xx)
      if (y >= 0) right(y) / tot else 1 - left(y) / tot
    }, 0),
    log_density = function(x) vapply(x, function(xx) {
      K(y_of(xx)) - log(tot * sd)
    }, 0),
    v = v, sd = sd)
}

test_that("pmhn and dmhn keep their accuracy as the law narrows toward the spacing of the doubles", {
  skip_if_not_installed("Rmpfr")
  # Under a tilt z = gamma / sqrt(beta) the law's standard deviation is a
  # relative sqrt(2) / z of its peak, about z / 2 in units of 1 / sqrt(beta),
  # so the doubles near the peak are z 2^-53 / sqrt(2) standard deviations
  # apart: 1e-8 of one at z = 1e8, 0.08 at z = 1e15.  An abscissa, or a
  # centre gamma / (2 beta), rounded to those doubles moves each tail by that
  # fraction of the density, which beta = 3 makes certain (gamma / (2 beta)
  # is then inexact).  Measured from the peak in exact offsets
  # (mhn_cdf_peak.h), the tails and the log density agree with the reference
  # to about 1e-15 at each z here; 1e-10 allows for the reference's own
  # accuracy and is far inside the working tolerance, 1.5e-8.
  for (a in c(0.5, 2, 10)) {
    for (b in c(1, 3)) {
      for (z in c(1e8, 1e10, 1e12, 1e14, 1e15)) {
        g <- z * sqrt(b)
        info <- sprintf("alpha=%g beta=%g z=%g", a, b, z)
        ref <- peak_law_reference(a, b, g)
        x <- ref$v + c(-3, -1, -0.3, 0, 0.3, 1, 3) * ref$sd
        expect_lt(max(abs(pmhn(x, a, b, g) - ref$lower(x))), 1e-10,
                  label = info)
        expect_lt(max(abs(pmhn(x, a, b, g, lower.tail = FALSE) -
                            ref$upper(x))), 1e-10, label = info)
        expect_lt(max(abs(dmhn(x, a, b, g, log = TRUE) -
                            ref$log_density(x))), 1e-10, label = info)
      }
    }
  }
  # Far below such a peak -- 7e8 standard deviations at z = 1e9 -- the upper
  # tail is 1 to double precision.
  for (z in c(1e8, 1e9, 1e14, 3e15)) {
    expect_identical(pmhn(1, 2, 1, z, lower.tail = FALSE), 1,
                     label = sprintf("z=%g", z))
  }
})

test_that("pmhn and dmhn at a large shape are measured from the peak", {
  skip_if_not_installed("Rmpfr")
  # Assembled from log Psi and the log kernel, which are each about
  # alpha log(v) in size, log f and log F would carry a rounding of about
  # eps alpha log(v): 1e-7 at alpha = 1e8 and 1e-2 at 1e15.  Measured from
  # the peak (mhn_cdf_peak.h) nothing of that size is formed, and the tails
  # and the log density agree with the reference to about 1e-12, the
  # reference's own accuracy in log f.  The standard deviation is a relative
  # 1 / sqrt(2 alpha) or less of v, 2e-8 at 1e15, still wider than the
  # spacing of the doubles.  1e-10 is far inside the working tolerance.
  # gamma = 0 is measured from the peak too: the closed form
  # P(alpha / 2, beta x^2) carries the rounding of beta x^2, up to 3 eps
  # sqrt(alpha / 2) standard deviations of the Gamma law of beta X^2, which
  # within three standard deviations moves F by up to 0.4 times that, 6e-9,
  # and log f by up to three times it, 4.5e-8, at alpha = 1e15.
  for (a in c(1e8, 1e10, 1e12, 1e15)) {
    for (bz in list(c(1, -10), c(1, 1), c(1, 10), c(3, -1e10),
                    c(1e-10, 5), c(1, 0), c(1e300, 0))) {
      b <- bz[1]
      g <- bz[2] * sqrt(b)
      info <- sprintf("alpha=%g beta=%g z=%g", a, b, bz[2])
      ref <- peak_law_reference(a, b, g)
      x <- ref$v + c(-3, -1, 0, 1, 3) * ref$sd
      expect_lt(max(abs(pmhn(x, a, b, g) - ref$lower(x))), 1e-10,
                label = info)
      expect_lt(max(abs(pmhn(x, a, b, g, lower.tail = FALSE) -
                          ref$upper(x))), 1e-10, label = info)
      expect_lt(max(abs(dmhn(x, a, b, g, log = TRUE) -
                          ref$log_density(x))), 1e-10, label = info)
    }
  }
})

test_that("the peak frame takes over from the ordinary assembly without a step", {
  # The density and the tails are assembled from log Psi and the log kernel
  # while the larger m of alpha and z = gamma / sqrt(beta) is at most 1e6
  # and m (1 + |log v*|), the size of the largest terms, at most 1e7, where
  # their rounding is at most about 2e-9; past either they are measured from
  # the peak.  At beta = 1 the shape switches at 1e6: there and at the next
  # double above it, which moves the law by a relative 2.2e-16 (its mode by
  # 2e-13 standard deviations), the two agree to that.
  up <- 1 + .Machine$double.eps
  for (z in c(-10, 1, 10)) {
    x <- mhn_mode(1e6, 1, z) + c(-1, 0, 1)
    expect_lt(max(abs(pmhn(x, 1e6, 1, z) - pmhn(x, 1e6 * up, 1, z))), 2e-9,
              label = sprintf("z=%g", z))
    expect_lt(max(abs(dmhn(x, 1e6, 1, z, log = TRUE) -
                        dmhn(x, 1e6 * up, 1, z, log = TRUE))), 2e-9,
              label = sprintf("z=%g", z))
  }
  # The tilt switches where z (1 + log v*) passes 1e7, v* about z / 2 at
  # beta = 1: near z = 7.25e5.  A relative 1e-6 below and above it (clear of
  # any rounding in locating the switch) moves the mode by 0.36, about half
  # a standard deviation, so the two sides are compared with x moved along
  # by the same amount: the law is normal about the mode to within a
  # relative (alpha - 1) / z^2 there, 2e-11, and within 2e-9 each side then
  # agrees with the other.
  crit <- function(z) z * (1 + log((z + sqrt(z^2 + 8 * 2)) / 4))
  zs <- stats::uniroot(function(z) crit(z) - 1e7, c(5e5, 1e6),
                       tol = 1e-3)$root
  for (a in c(0.5, 2, 10)) {
    z_lo <- zs * (1 - 1e-6)
    z_hi <- zs * (1 + 1e-6)
    x_lo <- mhn_mode(a, 1, z_lo) + c(-2, -0.3, 0, 0.3, 2)
    x_hi <- x_lo + (mhn_mode(a, 1, z_hi) - mhn_mode(a, 1, z_lo))
    expect_lt(max(abs(pmhn(x_lo, a, 1, z_lo) - pmhn(x_hi, a, 1, z_hi))), 4e-9,
              label = sprintf("alpha=%g", a))
    expect_lt(max(abs(pmhn(x_lo, a, 1, z_lo, lower.tail = FALSE) -
                        pmhn(x_hi, a, 1, z_hi, lower.tail = FALSE))), 4e-9,
              label = sprintf("alpha=%g", a))
    expect_lt(max(abs(dmhn(x_lo, a, 1, z_lo, log = TRUE) -
                        dmhn(x_hi, a, 1, z_hi, log = TRUE))), 4e-9,
              label = sprintf("alpha=%g", a))
  }
})

test_that("the peak frame is taken by the size of the terms, not by alpha or z alone", {
  skip_if_not_installed("Rmpfr")
  # At beta = 1e-300 or 1e300 |log v*| is about 350, and the ordinary
  # assembly adds (alpha / 2) log beta and alpha log x, each some 3.5e8 at
  # alpha = 1e6, which cancel: assembled so, the tails would be 5e-8 to 9e-8
  # off at a shape just below 1e6 and 1e-8 to 2e-8 at a tilt just below 1e6.
  # Wherever m (1 + |log v*|) passes 1e7 the law is measured from the peak,
  # and agrees with the reference to about 1e-14; 1e-10 allows for the
  # reference.  Just below that size the ordinary assembly is used, whose
  # rounding is then about 2e-9 (eps 1e7 times a small factor); 1e-8 is
  # within the working tolerance, 1.5e-8.
  pts <- list(c(999999, 1e-300, -3), c(999999, 1e300, 1), c(5e5, 1e300, 1),
              c(2e5, 1e300, 30), c(10, 1e-300, 9.9e5), c(2, 1e300, 9.9e5))
  for (p in pts) {
    a <- p[1]
    b <- p[2]
    g <- p[3] * sqrt(b)
    info <- sprintf("alpha=%g beta=%g z=%g", a, b, p[3])
    ref <- peak_law_reference(a, b, g)
    x <- ref$v + c(-3, -1, 0, 1, 3) * ref$sd
    expect_lt(max(abs(pmhn(x, a, b, g) - ref$lower(x))), 1e-10, label = info)
    expect_lt(max(abs(pmhn(x, a, b, g, lower.tail = FALSE) - ref$upper(x))),
              1e-10, label = info)
    expect_lt(max(abs(dmhn(x, a, b, g, log = TRUE) - ref$log_density(x))),
              1e-10, label = info)
    expect_lt(max(abs(qmhn(ref$lower(x), a, b, g) - x)) / ref$sd, 1e-9,
              label = info)
  }
  for (p in list(c(2.7e4, 1e-300, -3), c(2.7e4, 1e300, 1))) {
    a <- p[1]
    b <- p[2]
    g <- p[3] * sqrt(b)
    info <- sprintf("alpha=%g beta=%g z=%g", a, b, p[3])
    ref <- peak_law_reference(a, b, g)
    x <- ref$v + c(-3, -1, 0, 1, 3) * ref$sd
    expect_lt(max(abs(pmhn(x, a, b, g) - ref$lower(x))), 1e-8, label = info)
    expect_lt(max(abs(pmhn(x, a, b, g, lower.tail = FALSE) - ref$upper(x))),
              1e-8, label = info)
    expect_lt(max(abs(dmhn(x, a, b, g, log = TRUE) - ref$log_density(x))),
              1e-8, label = info)
  }
  # gamma = 0 stays on the closed form P(alpha / 2, t^2), t = sqrt(beta) q,
  # below a shape of 1e6 at every beta: written in t, it loses nothing to the
  # size of beta, and gives the same double as at beta = 1.
  q <- sqrt(5e5 / 2e-300) * c(0.999, 1, 1.001)
  expect_identical(pmhn(q, 5e5, 1e-300, 0), pmhn(sqrt(1e-300) * q, 5e5, 1, 0))
})

test_that("the peak frame keeps log(x / v*) where it leaves the normal range", {
  # log(1 + rho) - rho is formed from log(x / v*) far below the peak, and
  # x / v* underflows for a subnormal x, or at x = 1e-300 with v* above
  # 1e8; its logarithm is then taken as log x - log v*.  At (alpha, 1, 2e6)
  # the law is measured from the peak (z = 2e6), v* is about 1e6, and below
  # the peak
  #   F(x) = x^alpha / (alpha C) (1 + O(gamma x)),
  #   C = int x^(alpha-1) exp(-x^2 + gamma x) dx
  #     = exp(gamma^2 / 4) (gamma / 2)^(alpha-1) sqrt(pi)
  #       (1 + (alpha - 1)(alpha - 2) / gamma^2 + ...),
  # the correction being 3e-13 at alpha = 0.3.  log C is 1e12 in size, so
  # the result carries about 1e-4 of rounding; 1e-3 allows for it, far below
  # the (alpha - 1) log x term (500) that a wrong log(x / v*) would lose.
  g <- 2e6
  log_C <- function(a) g^2 / 4 + (a - 1) * log(g / 2) + 0.5 * log(pi)
  for (a in c(0.3, 2)) {
    for (x in c(1e-320, 1e-300)) {
      info <- sprintf("alpha=%g x=%g", a, x)
      expect_lt(abs(pmhn(x, a, 1, g, log.p = TRUE) -
                      (a * log(x) - log(a) - log_C(a))), 1e-3, label = info)
    }
  }
  # log F stays finite at such an x, so the quantile's search towards 0 does
  # not stop there: where the quantile lies below the smallest positive
  # double it is returned as 0.
  expect_identical(qmhn(-1e100, 2, 1, 1e7, log.p = TRUE), 0)
  expect_identical(qmhn(-1e10, 3e6, 1, 0, log.p = TRUE), 0)
})

test_that("a peak placed a few standard deviations from a double is resolved", {
  skip_if_not_installed("Rmpfr")
  # At a shape far beyond 1e31 the law is narrower than the spacing of the
  # doubles near its peak v*, and its position is carried as a double plus
  # an offset found from the residual of 2 beta x^2 - gamma x - alpha at that
  # double.  The residual is summed so that its rounding is of order eps^3
  # v*, which leaves the offset resolved to well within 1e-11 standard
  # deviations up to a shape of about 1e70.  Here x is a double, alpha is
  # 2 beta x^2 rounded, and gamma is chosen so that v* lies 1.5 standard
  # deviations above x; the law is normal about v* to within relative
  # corrections of 1 / sqrt(alpha), so F(x) is pnorm of the offset and
  # log f(x) that of the normal density, formed here in Rmpfr.  1e-12
  # covers the rounding of the results; the quantile of F(x) is x itself.
  for (b in c(3, 0.7, 1)) {
    for (la in c(43, 50, 60)) {
      x <- sqrt(10^la / (2 * b)) * 1.03
      a <- 2 * b * x^2
      X <- Rmpfr::mpfr(x, 600)
      A <- Rmpfr::mpfr(a, 600)
      B <- Rmpfr::mpfr(b, 600)
      V <- X + 1.5 * X / sqrt(2 * A)
      g <- Rmpfr::asNumeric((2 * B * V^2 - A) / V)
      G <- Rmpfr::mpfr(g, 600)
      V <- (G + sqrt(G^2 + 8 * A * B)) / (4 * B)
      sd <- V / sqrt(A + 2 * B * V^2)
      t <- Rmpfr::asNumeric((X - V) / sd)
      info <- sprintf("alpha=%g beta=%g", a, b)
      expect_no_warning(F <- pmhn(x, a, b, g))
      expect_lt(abs(F - pnorm(t)), 1e-12, label = info)
      expect_lt(abs(pmhn(x, a, b, g, lower.tail = FALSE) -
                      pnorm(t, lower.tail = FALSE)), 1e-12, label = info)
      expect_no_warning(ld <- dmhn(x, a, b, g, log = TRUE))
      expect_lt(abs(ld - Rmpfr::asNumeric(-t^2 / 2 - log(sd) -
                                            0.5 * log(2 * Rmpfr::Const("pi", 600)))),
                1e-12 * max(1, abs(ld)), label = info)
      expect_identical(qmhn(pnorm(t), a, b, g), x, label = info)
    }
  }
})

test_that("a peak that falls on a double is not read as a far tail", {
  # With alpha = beta = gamma (the double nearest 1e300) the peak of the
  # density of log X, the root of 2 beta x^2 - gamma x - alpha = 0, is 1
  # exactly, and the law is normal about it with standard deviation
  # 1 / sqrt(3 beta), some 1e-134 of the spacing of the doubles there.  The
  # offset of the peak from the double next to it is found to within an
  # error far above that width but tiny against the spacing, so at x = 1 it
  # cannot say on which side x lies, and d and p are NA with a warning
  # there, never a value from a far tail.  At the neighbouring doubles,
  # 1.9e134 standard deviations or more away, log f and log F are the normal
  # ones, -t^2 / 2 to a relative 1e-12 (the terms beyond are some 1e-266 of
  # it).
  a <- 1e300
  expect_warning(F <- pmhn(1, a, a, a), "could not be evaluated")
  expect_true(is.na(F) && !is.nan(F))
  expect_warning(d <- dmhn(1, a, a, a, log = TRUE), "could not be evaluated")
  expect_true(is.na(d) && !is.nan(d))
  expect_warning(d <- dmhn(1, c(a, a), a, a, log = TRUE),
                 "could not be evaluated")
  expect_true(all(is.na(d)))
  sd <- 1 / sqrt(3 * a)
  x_lo <- 1 - 2^-53
  x_hi <- 1 + 2^-52
  expect_lt(abs(pmhn(x_lo, a, a, a, log.p = TRUE) /
                  (-((1 - x_lo) / sd)^2 / 2) - 1), 1e-12)
  expect_identical(pmhn(x_hi, a, a, a), 1)
  expect_lt(abs(dmhn(x_hi, a, a, a, log = TRUE) /
                  (-((x_hi - 1) / sd)^2 / 2) - 1), 1e-12)
})

test_that("the series length is tested before it could overflow an integer type", {
  # The Lemma 10(d) length K grows like 0.72 z^2, which a 32-bit long
  # cannot hold from |z| of about 5.5e4 and a 64-bit one from about 3.6e9.
  # The lengths are held as doubles and tested against the ceiling before
  # any is converted, so the series declines (NaN) and pmhn and qmhn take the
  # quadrature.  At (3, 1, -1e10) the law is Gamma(3, rate 1e10) to within a
  # relative 12 / z^2, so pmhn(2e-10) is pgamma(2, 3) and the median is
  # qgamma(0.5, 3) / 1e10; 1e-8 is the working tolerance.
  for (g in c(-6e4, -1e5, -4e9, -1e10)) {
    expect_true(is.nan(mhn:::.pmhn_force(1e-5, 3, 1, g, "series")),
                label = sprintf("gamma=%g", g))
  }
  expect_equal(pmhn(2e-10, 3, 1, -1e10), pgamma(2, 3), tolerance = 1e-8)
  expect_equal(qmhn(0.5, 3, 1, -1e10) * 1e10, qgamma(0.5, 3),
               tolerance = 1e-8)
})

test_that("pmhn past |z| = 20 or alpha = 100 keeps the working tolerance", {
  # Outside |z| <= 20 with alpha <= 100, z = gamma / sqrt(beta), the lower
  # tail below the peak of the density of log X of a shape of 1 or more is
  # taken by quadrature rather than the Lemma 1b series, whose length grows
  # like z^2; below 1 the series is still tried first (see the next test).
  # The points cover a positive and a negative tilt, a large shape, a lower
  # tail of 1e-292, and two shapes below 1, where the quadrature is also
  # checked directly.  At (0.3, 1, 21) two of them lie past the kernel's peak
  # by 1e-13 and 1e-11 of its abscissa, where the quadrature past the peak
  # spans a range many orders of magnitude narrower than the law.
  # References: the integral of the law of log X,
  # measured from its peak (formed in 300-bit arithmetic), by adaptive
  # quadrature to about 1e-13; the largest error of the quadrature seen at
  # these points is 1.4e-9, at (0.7, 1, -300).  1e-8 lies below the working
  # tolerance, sqrt(eps) = 1.5e-8.
  cases <- list(
    list(1.5, 300, c(124.2, 145, 149.5),
         c(8.1969701490442054e-292, 7.5556005238390443e-13,
           0.23901712487545718)),
    list(0.3, 21, c(10.466560169308019, 10.46656016941164,
                    10.361894567613904),
         c(0.5000583339376663, 0.50005833399603361, 0.44131605435931531)),
    list(20, -300, c(0.056975953634869611, 0.063306615149855125),
         c(0.27186331389403467, 0.43943103039953496)),
    list(0.7, -300, 6.4675e-05, 0.069127764986760321),
    list(3000, 5, c(35.994193064939147, 39.993547849932384),
         c(5.5848427660284779e-16, 0.49836523199711608)))
  for (cs in cases) {
    info <- sprintf("alpha=%g gamma=%g", cs[[1]], cs[[2]])
    expect_lt(max(abs(pmhn(cs[[3]], cs[[1]], 1, cs[[2]]) / cs[[4]] - 1)),
              1e-8, label = info)
    expect_lt(max(abs(pmhn(cs[[3]], cs[[1]], 1, cs[[2]], log.p = TRUE) -
                        log(cs[[4]]))), 1e-8, label = info)
    if (cs[[1]] < 1) {
      v <- vapply(cs[[3]], function(x) {
        mhn:::.pmhn_force(x, cs[[1]], 1, cs[[2]], "integrate")
      }, 0)
      expect_lt(max(abs(v / cs[[4]] - 1)), 1e-8, label = info)
    }
  }
})

test_that("pmhn takes a shape below 1 from the series at any tilt", {
  # For alpha < 1 the lower tail below the peak of the density of log X is
  # tried with the Lemma 1b series first at every tilt up to z = 1e3, and
  # not only for |z| <= 20.  There the quadrature, tanh-sinh against the
  # power t^(alpha - 1) near the origin, misses by more than the working
  # tolerance at scattered points, such as these three, where it is off by
  # 4.6e-8, 3.5e-7 and 4.7e-8 in F; the series is within 1e-13 at each.
  # References: log F from the power series in gamma of Lemma 1b, summed
  # in MPFR at up to 3128 bits, which covers the cancellation of its
  # alternating terms, and agrees with the integral of the law of log X
  # measured from its peak to 1e-15.  1e-8 lies below the working
  # tolerance, sqrt(eps) = 1.5e-8.
  cases <- list(c(0.502, -79.7, 0.00155, -0.96958701794933377),
                c(0.19869, -86.961, 1.8361e-05, -1.1947669672593078),
                c(0.031666, 24.7791, 11.1264, -3.1698820385943565))
  for (cs in cases) {
    info <- sprintf("alpha=%g gamma=%g", cs[1], cs[2])
    expect_lt(abs(pmhn(cs[3], cs[1], 1, cs[2], log.p = TRUE) - cs[4]), 1e-8,
              label = info)
    expect_lt(abs(pmhn(cs[3], cs[1], 1, cs[2]) / exp(cs[4]) - 1), 1e-8,
              label = info)
  }
})

test_that("pmhn's upper tail holds for a tiny shape under a positive tilt", {
  # As alpha -> 0 at beta = gamma = 1, P(X > q) is alpha times the integral of
  # t^-1 exp(-t^2 + t) over (q, Inf), up to a relative alpha log(1 / q).  The
  # integral is split so that integrate() never meets the 1 / t singularity:
  # t^-1 expm1(-t^2 + t) over (q, 1), plus log(1 / q), plus the tail from 1.
  tail_integral <- function(q) {
    integrate(function(t) expm1(-t^2 + t) / t, q, 1, rel.tol = 1e-12)$value +
      log(1 / q) +
      integrate(function(t) exp(-t^2 + t) / t, 1, Inf, rel.tol = 1e-12)$value
  }
  for (alpha in c(1e-12, 1e-16, 1e-20, 1e-100)) {
    for (q in c(1e-30, 1e-200)) {
      expect_equal(pmhn(q, alpha, 1, 1, lower.tail = FALSE) /
                     (alpha * tail_integral(q)), 1, tolerance = 1e-8)
    }
  }
})

test_that("the lower log tail keeps its relative accuracy near 0", {
  # Below the peak of a tiny shape under a small tilt F is 1 less alpha
  # times a moderate factor, and log F from the lower tail, assembled from
  # terms of the size of log Psi that cancel, loses its relative accuracy:
  # in 0.1.1 it was 0 at (1e-100, 1, 5), q = 2.5, and 1.1e-4 off at
  # (1e-12, 1, 5), q = 2.5.  It is now log(1 - Q) from the upper tail Q
  # where that is the more accurate.  The reference is Q by quadrature in R,
  # with 2 / Gamma(alpha / 2) = alpha / Gamma(1 + alpha / 2) taken out of
  # the normalising constant.  The error measured is up to 1.4e-13 where the
  # upper tail is taken, and 7.3e-12 at (1e-6, 1, 5), q = 2, where the lower
  # tail stays, its error being kept below 1e-9 of log F by design.  For
  # alpha <= 1e-8 the upper tail is taken at every point here, and its
  # relative error in log F is below twice the error bound of log Q, at most
  # 6.9e-12 at these points, which 1e-11 allows for; at alpha = 1e-6 the
  # lower tail stays at some of them, and 1e-10 is kept.
  upper_ref <- function(a, z, u) {
    num <- integrate(function(t) exp((a - 1) * log(t) - t^2 + z * t),
                     u, Inf, rel.tol = 1e-13)$value
    rest <- integrate(function(t) exp((a - 1) * log(t) - t^2 + z * t) *
                        -expm1(-z * t), 0, Inf, rel.tol = 1e-13)$value
    cc <- a / gamma(1 + a / 2)
    num * cc / (1 + cc * rest)
  }
  for (a in c(2^-1022, 1e-300, 1e-100, 1e-20, 1e-12, 1e-10, 1e-8, 1e-6)) {
    for (z in c(0.5, 5)) {
      u <- c(0.1, 0.4) * z
      ref <- log1p(-vapply(u, upper_ref, numeric(1), a = a, z = z))
      for (b in c(1, 1e-300)) {
        expect_close(pmhn(u / sqrt(b), a, b, z * sqrt(b), log.p = TRUE), ref,
                     if (a <= 1e-8) 1e-11 else 1e-10,
                     info = sprintf("alpha=%g z=%g beta=%g", a, z, b))
      }
    }
  }

  # Under a strong negative tilt at an extreme beta the lower tail's terms
  # include (alpha - 1) log x, of the size of |log beta| / 2, and its error
  # bound counts them, so that at (3e-6, 1e-300, -1e-146), q = 2.7e140, the
  # upper tail serves instead.  0.1.1 was 6.2e-11 off there, at
  # beta = 1 and 1e-300 (it stops with an error at 1e300), and the upper
  # tail now gives 7.5e-15.  The reference is independent of the package:
  # with s = |z| t, F = 1 - N / D, N = int_{|z| u}^Inf s^(a-1)
  # exp(-s - s^2 / z^2) ds, on the log axis, and D = Gamma(a) -
  # int_0^Inf s^(a-1) exp(-s) (1 - exp(-s^2 / z^2)) ds.  At beta = 1 the
  # lower tail may serve, within the design's 1e-9.
  ref_neg <- function(a, z, u) {
    lo <- log(abs(z) * u)
    f <- function(v) exp(a * v - exp(v) - exp(2 * v) / z^2)
    br <- sort(unique(c(lo, -5, 0, 2, 3, 4, 5)))
    br <- br[br >= lo]
    N <- sum(mapply(function(l, h) integrate(f, l, h, rel.tol = 1e-13)$value,
                    br[-length(br)], br[-1])) +
      integrate(f, 5, Inf, rel.tol = 1e-13)$value
    R <- integrate(function(s) exp((a - 1) * log(s) - s) * -expm1(-(s / z)^2),
                   0, Inf, rel.tol = 1e-13)$value
    log1p(-N / (gamma(a) - R))
  }
  a <- 3e-6; z <- -1e4
  u <- a / abs(z) * c(0.01, 0.1, 0.333, 0.9)
  ref <- vapply(u, ref_neg, numeric(1), a = a, z = z)
  for (b in c(1e-300, 1, 1e300)) {
    expect_close(pmhn(u / sqrt(b), a, b, z * sqrt(b), log.p = TRUE), ref,
                 if (b == 1) 1e-10 else 1e-12,
                 info = sprintf("alpha=%g z=%g beta=%g", a, z, b))
  }

  # At a large tilt and an extreme beta the lower tail is the more accurate
  # of the two, and stays: at the peak of z = 1e4 its scale-free values agree
  # with those at beta = 1 to 2.1e-12, where the upper tail, whose abscissae
  # on the log axis are rounded, leaves up to 6.7e-10.
  for (a in c(1e-4, 0.01, 0.3, 0.5, 0.9)) {
    ref <- pmhn(5e3, a, 1, 1e4, log.p = TRUE)
    for (b in c(1e-300, 1e300)) {
      expect_close(pmhn(5e3 / sqrt(b), a, b, 1e4 * sqrt(b), log.p = TRUE),
                   ref, 2e-11, info = sprintf("alpha=%g beta=%g", a, b))
    }
  }

  # So it does where Q is close to one half, between the median and the peak
  # of the density of log X of a large shape, where (alpha - 1) log x is the
  # larger part of both error bounds: there log F is that of the lower tail,
  # the log of F on the linear scale.  Taken from the upper tail at these
  # points it moved by up to 3.9e-9, of the order of the error of either
  # tail there (about 2e-9, see kUpperLogRelTol in mhn_pmhn_state.h).
  for (p in list(c(9e5, 1, -1e3), c(9e5, 1, 0.5), c(9e5, 1, 1e3),
                 c(1e4, 1e300, 0.5), c(1e4, 1e300, 10), c(1e4, 1e-300, -1e3))) {
    a <- p[1]; b <- p[2]; g <- p[3] * sqrt(b)
    d <- if (a == 9e5) c(0, 1e-7, 2e-7) else c(0, 1e-6, 1e-5)
    q <- (g + sqrt(g^2 + 8 * b * a)) / (4 * b) * (1 - d)
    expect_close(pmhn(q, a, b, g, log.p = TRUE), log(pmhn(q, a, b, g)), 1e-14,
                 info = sprintf("alpha=%g beta=%g z=%g", a, b, p[3]))
  }
})

test_that("a peak of the density of log X below the smallest normal double", {
  # The upper-tail quadrature measures its integrand from the peak of the
  # density of log X, with the slope there taken from the peak's
  # stationarity, alpha / v.  That holds only at the exact peak, and a peak
  # below 2.2e-308 is resolved only to the spacing of the subnormal doubles:
  # in 0.1.1 log P(X > q) was 3.5e-7 off at the first point, and the lower
  # tail -Inf.  The reference is independent of the package: with
  # s = |z| t, P(X > q) = N / D, N = int_{|z| u}^Inf s^(a-1)
  # exp(-s - s^2 / z^2) ds on the log axis, D = Gamma(a) -
  # int_0^Inf s^(a-1) exp(-s) (1 - exp(-s^2 / z^2)) ds, u = sqrt(beta) q.
  upper_neg <- function(a, z, u) {
    lo <- log(abs(z) * u)
    f <- function(v) exp(a * v - exp(v) - exp(2 * v) / z^2)
    br <- sort(unique(c(lo, -5, 0, 2, 3, 4, 5)))
    br <- br[br >= lo]
    N <- sum(mapply(function(l, h) integrate(f, l, h, rel.tol = 1e-13)$value,
                    br[-length(br)], br[-1])) +
      integrate(f, 5, Inf, rel.tol = 1e-13)$value
    R <- integrate(function(s) exp((a - 1) * log(s) - s) * -expm1(-(s / z)^2),
                   0, Inf, rel.tol = 1e-13)$value
    N / (gamma(a) - R)
  }
  pts <- list(c(1.41e-217, 7.65e208, -37.47, 1.482e-323),
              c(4.81e-204, 2.77e225, -1972, 3.458e-323),
              c(2.85e-187, 1.48e256, -208.8, 2.346e-319))
  for (p in pts) {
    a <- p[1]; b <- p[2]; z <- p[3]; q <- p[4]
    Q <- upper_neg(a, z, q * sqrt(b))
    info <- sprintf("alpha=%g beta=%g z=%g q=%g", a, b, z, q)
    expect_close(pmhn(q, a, b, z * sqrt(b), lower.tail = FALSE, log.p = TRUE),
                 log(Q), 1e-12, info = info)
    expect_close(pmhn(q, a, b, z * sqrt(b), log.p = TRUE), log1p(-Q), 1e-12,
                 info = info)
  }
})
