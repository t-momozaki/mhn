# Tests for qmhn() -- the MHN quantile function.

# ============================================================
# 1. Boundary values
# ============================================================

test_that("qmhn(0) == 0 and qmhn(1) == Inf", {
  for (a in c(0.5, 1, 2, 5)) for (b in c(0.5, 1, 3)) for (g in c(-2, 0, 1, 3)) {
    expect_equal(qmhn(0, a, b, g), 0,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
    expect_equal(qmhn(1, a, b, g), Inf,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("qmhn handles log.p and lower.tail at boundaries", {
  expect_equal(qmhn(-Inf, 2, 1, 1, log.p = TRUE), 0)
  expect_equal(qmhn(0, 2, 1, 1, log.p = TRUE), Inf)         # log(1) = 0
  expect_equal(qmhn(0, 2, 1, 1, lower.tail = FALSE), Inf)
  expect_equal(qmhn(1, 2, 1, 1, lower.tail = FALSE), 0)
})

test_that("qmhn returns NaN for out-of-range probabilities", {
  # With the warning qnorm gives for the same inputs (see section 8).
  expect_warning(q <- qmhn(-0.5, 2, 1, 1), "NaNs produced")
  expect_true(is.nan(q))
  expect_warning(q <- qmhn(2, 2, 1, 1), "NaNs produced")
  expect_true(is.nan(q))
  expect_warning(q <- qmhn(1.5, 2, 1, 1, log.p = TRUE), "NaNs produced")  # log p > 0
  expect_true(is.nan(q))
})

# ============================================================
# 2. Round-trip: pmhn(qmhn(p)) ~= p
# ============================================================

test_that("round-trip pmhn(qmhn(p)) ~= p", {
  p_vals <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
  params <- expand.grid(alpha = c(0.5, 1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-3, -1, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    x <- qmhn(p_vals, a, b, g)
    back <- pmhn(x, a, b, g)
    expect_equal(back, p_vals, tolerance = 1e-12,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 3. Reverse round-trip: qmhn(pmhn(x)) ~= x
# ============================================================

test_that("reverse round-trip qmhn(pmhn(x)) ~= x", {
  x_vals <- c(0.1, 0.3, 0.5, 1, 2, 3)
  params <- expand.grid(alpha = c(1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0.5, 2))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    p <- pmhn(x_vals, a, b, g)
    back <- qmhn(p, a, b, g)
    # A looser relative tolerance: the inversion goes through F, whose
    # rounding near 1 moves the quantile by about eps / f(x), 2.7e-8 of x at
    # x = 3 for (1.5, 2, -2).
    expect_equal(back, x_vals, tolerance = 1e-5,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 4. Monotonicity in p
# ============================================================

test_that("qmhn is monotone non-decreasing in p", {
  p <- seq(0.001, 0.999, length.out = 200)
  params <- expand.grid(alpha = c(0.5, 1, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-2, 0, 1, 3))
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    x <- qmhn(p, a, b, g)
    expect_false(is.unsorted(x),
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 5. lower.tail / log.p consistency
# ============================================================

test_that("qmhn lower.tail equivalence", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  params <- list(c(2, 1, 1), c(1.5, 0.5, -1), c(0.5, 1, 2))
  for (par in params) {
    a <- par[1]; b <- par[2]; g <- par[3]
    lower <- qmhn(p, a, b, g)
    upper <- qmhn(1 - p, a, b, g, lower.tail = FALSE)
    expect_equal(lower, upper, tolerance = 1e-8,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

test_that("qmhn log.p equivalence", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  params <- list(c(2, 1, 1), c(1.5, 0.5, -1), c(0.5, 1, 2))
  for (par in params) {
    a <- par[1]; b <- par[2]; g <- par[3]
    linear <- qmhn(p, a, b, g)
    logp   <- qmhn(log(p), a, b, g, log.p = TRUE)
    expect_equal(linear, logp, tolerance = 1e-8,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 6. Special-case dispatch
# ============================================================

test_that("gamma = 0: qmhn == sqrt(qgamma(p, alpha/2, scale = 1/beta))", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  for (a in c(0.5, 1, 1.5, 2, 4)) for (b in c(0.5, 1, 3)) {
    ref <- sqrt(qgamma(p, shape = a/2, scale = 1/b))
    got <- qmhn(p, alpha = a, beta = b, gamma = 0)
    expect_equal(got, ref, tolerance = 1e-10,
                 label = sprintf("a=%.1f b=%.1f", a, b))
  }
})

test_that("alpha = 1: qmhn matches truncated normal inverse", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  for (b in c(0.5, 1, 2)) for (g in c(-2, 0, 1, 3)) {
    mu <- g / (2 * b); sigma <- 1 / sqrt(2 * b)
    Phi_low <- pnorm(0, mu, sigma)
    # Reference via direct inversion of the truncated-normal CDF.
    ref <- qnorm(Phi_low + p * (1 - Phi_low), mu, sigma)
    got <- qmhn(p, alpha = 1, beta = b, gamma = g)
    expect_equal(got, ref, tolerance = 1e-8,
                 label = sprintf("b=%.1f g=%+.1f", b, g))
  }
})

# ============================================================
# 7. Vectorization / recycling
# ============================================================

test_that("qmhn vectorizes p", {
  p <- c(0.1, 0.5, 0.9)
  out <- qmhn(p, 2, 1, 1)
  expect_length(out, 3)
  expect_equal(out, sapply(p, function(pp) qmhn(pp, 2, 1, 1)),
               tolerance = 1e-10)
})

test_that("qmhn recycles parameter vectors", {
  p <- c(0.3, 0.5, 0.7)
  a <- c(1, 2)
  b <- 1
  g <- c(0, 1, -1)
  n <- max(length(p), length(a), length(b), length(g))
  out <- qmhn(p, a, b, g)
  expect_length(out, n)
  ref <- vapply(seq_len(n), function(i) {
    qmhn(p[((i-1) %% length(p)) + 1],
         a[((i-1) %% length(a)) + 1],
         b[((i-1) %% length(b)) + 1],
         g[((i-1) %% length(g)) + 1])
  }, numeric(1))
  expect_equal(out, ref, tolerance = 1e-10)
})

test_that("qmhn with empty p returns numeric(0)", {
  expect_equal(qmhn(numeric(0), 2, 1, 1), numeric(0))
})

# ============================================================
# 8. NA / NaN propagation and parameter validation
# ============================================================

test_that("qmhn propagates NA in p", {
  out <- qmhn(c(0.5, NA, 0.9), 2, 1, 1)
  expect_true(is.na(out[2]) && !is.nan(out[2]))
  expect_false(is.na(out[1]))
})

test_that("qmhn propagates NaN in p", {
  out <- qmhn(c(0.5, NaN, 0.9), 2, 1, 1)
  expect_true(is.nan(out[2]))
})

test_that("qmhn gives NaN with a warning for parameters outside the parameter space", {
  # As qgamma does; test-input-validation.R tests the rule in full.
  for (args in list(list(alpha = 0), list(alpha = -1), list(beta = 0),
                    list(gamma = Inf))) {
    expect_warning(out <- do.call(qmhn, c(list(0.5), args)), "NaNs produced")
    expect_identical(out, NaN)
  }
  expect_identical(qmhn(0.5, alpha = NA, beta = 1, gamma = 0), NA_real_)
})

test_that("qmhn warns when it produces NaN from a probability, as qnorm does", {
  # Probabilities outside [0, 1], and positive log probabilities, have no
  # quantile.  A NaN or NA input is passed through silently, again as in base R.
  expect_warning(out <- qmhn(c(-0.1, 1.1, 0.5), 2, 1, 1), "NaNs produced")
  expect_true(all(is.nan(out[1:2])))
  expect_warning(qmhn(0.3, 2, 1, 1, log.p = TRUE), "NaNs produced")
  expect_silent(out <- qmhn(c(NaN, NA, 0.5), 2, 1, 1))
  expect_true(is.nan(out[1]) && is.na(out[2]) && !is.nan(out[2]))
})

test_that("qmhn follows a quantile far out in the upper tail", {
  # A log upper-tail probability of -1e300 puts the quantile near 1e150,
  # 2^498 times the scale of the initial bracket, which is doubled until it
  # holds the root.  There log P(X > x) is -x^2 + x + O(log x), so the root
  # finder's relative 2^-46 in x moves it by about 3e-14 of itself; 1e-12
  # allows for that.
  q <- qmhn(-1e300, 2, 1, 1, lower.tail = FALSE, log.p = TRUE)
  expect_equal(q, 1e150, tolerance = 1e-12)
  expect_equal(pmhn(q, 2, 1, 1, lower.tail = FALSE, log.p = TRUE), -1e300,
               tolerance = 1e-12)
})

test_that("qmhn reports a quantile it cannot bracket as NA, with a warning", {
  # At alpha = 1e300, the subnormal beta = 1e-323 and gamma = -1e-160 the
  # peak of the law, sqrt(alpha / (2 beta)) about 2e311 (gamma^2 is nothing
  # next to 8 alpha beta), lies beyond the largest double.  Its
  # distribution function is then not formed (NA, with a warning), and the
  # quantile is a missing value, which the package keeps distinct from the
  # NaN it returns for an invalid probability.
  expect_warning(out <- qmhn(0.5, 1e300, 1e-323, -1e-160),
                 "could not bracket")
  expect_true(is.na(out) && !is.nan(out))
})

# ============================================================
# 9. Strong tilts and extreme scales
# ============================================================

test_that("qmhn inverts the Gamma limit when the mode lies below sqrt(eps)", {
  # The law differs from Gamma(alpha, rate = -gamma) by a relative correction
  # of order alpha (alpha + 1) beta / gamma^2 (see test-pmhn.R, section 11),
  # which is at most 4e-16 for the triples here, so qgamma is an exact
  # reference to double precision.  The modes here are 5e-9, 1e-9 and 7e-9.
  # The tolerance allows the working tolerance of the distribution function,
  # 1.5e-8, divided by x f(x) / F(x), which is at least 0.3 at these p.
  # The quantiles are of order 1e-8, so the comparison is made in relative
  # terms explicitly (expect_equal would fall back to an absolute one).
  p <- c(0.3, 0.5, 0.7, 0.9)
  for (par in list(c(1.5, -1e8), c(2, -1e9), c(3, -3e8))) {
    a <- par[1]; g <- par[2]
    expect_silent(q <- qmhn(p, a, 1, g))
    expect_lt(max(abs(q / stats::qgamma(p, a, rate = -g) - 1)), 1e-7,
              label = sprintf("alpha=%g gamma=%g", a, g))
  }
})

test_that("qmhn depends on beta only through its scale 1 / sqrt(beta)", {
  # If X ~ MHN(alpha, 1, z) then X / sqrt(beta) ~ MHN(alpha, beta,
  # z sqrt(beta)) (by a change of variable in the density of Sun et al.
  # 2023, Lemma 1a), so the quantiles must scale exactly.  The bracket
  # therefore has to be set on the scale of the distribution: an absolute
  # left end of sqrt(eps) is above every quantile once the scale falls below
  # about 1e-16 (beta near 1e32).  Both sides solve
  # to a relative 1.4e-14, and the distribution function is scale-free to
  # rounding, so 1e-10 leaves ample room.  The comparison is relative
  # element by element, since the quantiles span twenty decades.  The upper
  # tail is solved on log P(X > x), which for p near one is about -F(x) and
  # is known only to an absolute 1e-14 or so; there it fixes x to a relative
  # 1e-14 / (alpha F(x)) at best, so that end is compared for the lower tail
  # only.
  for (a in c(0.3, 3, 20)) {
    for (z in c(-1, 5)) {
      for (lt in c(TRUE, FALSE)) {
        p <- c(1e-6, 0.3, 0.5, 0.9, if (lt) 1 - 1e-6)
        ref <- qmhn(p, a, 1, z, lower.tail = lt)
        for (b in c(1e40, 1e160)) {
          got <- qmhn(p, a, b, z * sqrt(b), lower.tail = lt) * sqrt(b)
          expect_lt(max(abs(got / ref - 1)), 1e-10,
                    label = sprintf("alpha=%g z=%g beta=%g lower=%s",
                                    a, z, b, lt))
        }
      }
    }
  }
  # A very large scale: the quantile is 3.8e-251 while the bracket starts
  # near 1e142, some 1300 halvings above it.  F is proportional to x^0.01
  # there, so its tolerance of 1e-6 corresponds to a relative 1e-4 in x.
  q <- qmhn(1e-4, 0.01, 1e-300, -1e-150)
  expect_equal(pmhn(q, 0.01, 1e-300, -1e-150), 1e-4, tolerance = 1e-6)
})

test_that("qmhn resolves quantiles of tiny lower-tail probabilities", {
  # Near the origin F(x) behaves like x^alpha / (alpha Z), so for alpha > 1 the
  # quantile of p = 1e-100 lies some 50 decades below the bulk.  The bracket
  # must close in on it rather than leave the root finder a span of many
  # decades that it cannot reduce to tolerance in its iteration budget.  The
  # root is solved to a relative 1.4e-14 in x, so F(q) matches p to a relative
  # alpha times that, 6e-13 at alpha = 40; 1e-8 is the working tolerance.
  for (a in c(1.2, 2, 3, 5, 10, 40)) {
    for (g in c(-3, 3)) {
      for (p in 10^-c(70, 100, 200, 300)) {
        info <- sprintf("alpha=%g gamma=%g p=%g", a, g, p)
        q <- qmhn(p, a, 1, g)
        expect_lt(abs(pmhn(q, a, 1, g) / p - 1), 1e-8, label = info)
        q_log <- qmhn(log(p), a, 1, g, log.p = TRUE)
        expect_lt(abs(q_log / q - 1), 1e-12, label = info)
      }
    }
  }
  # An independent anchor: at (2, 1, 3), F(x) = x^2 / (2 Z) + O(x^3) with
  # Z = int_0^Inf t exp(-t^2 + 3 t) dt = 25.297364396949772 (Romberg
  # integration in Rmpfr at 320 bits), so q(1e-100) = sqrt(2 Z 1e-100) to a
  # relative 1e-49.
  expect_lt(abs(qmhn(1e-100, 2, 1, 3) / 7.1129971737587195e-50 - 1), 1e-12)
})

test_that("qmhn solves lower-tail log probabilities below log(DBL_MIN)", {
  # exp(p) loses significant bits below log(.Machine$double.xmin) = -708.4 and
  # is 0 below -745.1, while the quantile can still be a normal double.  At
  # (2, 1, 3), F(x) = x^2 / (2 Z) (1 + O(x)) with Z = 25.297364396949772 as
  # above (also 1/2 + (3/2) exp(9/4) sqrt(pi) pnorm(3 / sqrt(2)) in closed
  # form), so log F(x) = p at x = exp((p + log(2 Z)) / 2) to a relative 1e-150
  # at these p.  The root is solved to a relative 1.4e-14.  The first two p
  # are solved on F, the others on log F, so this also checks that the two
  # forms join up.  pmhn on the log scale must then give lp back: a relative
  # 1.4e-14 in q moves log F by alpha times that, 2.8e-14, a relative 4e-17
  # of these lp, so 1e-12 leaves room for the rounding of the series.
  lp <- c(-700, -708, -709, -745, -746, -1000, -1400)
  z_2 <- 2 * 25.297364396949772
  q <- qmhn(lp, 2, 1, 3, log.p = TRUE)
  expect_lt(max(abs(q / exp((lp + log(z_2)) / 2) - 1)), 1e-12)
  expect_lt(max(abs(pmhn(q, 2, 1, 3, log.p = TRUE) / lp - 1)), 1e-12)
  # The same through the recycled-parameter path.
  q_vec <- qmhn(lp, c(2, 2), 1, 3, log.p = TRUE)
  expect_close(q_vec, q, 1e-14)
})

test_that("qmhn resolves quantiles below the smallest normal double", {
  # For alpha < 1 F(x) behaves like x^alpha near the origin, so a moderate p can
  # put the quantile below .Machine$double.xmin = 2.2e-308, where the doubles
  # are evenly spaced 4.9e-324 apart.  There the quantile is determined only
  # to that spacing: a relative 6.6e-13 at 7.5e-312 and 1.2e-8 at 4.3e-316,
  # so F(q) matches p to alpha times half of that, 3e-9 in the second case.
  # 1e-7 allows for that and for the working tolerance of the distribution
  # function.
  expect_silent(q <- qmhn(1e-250, 0.8, 1, 3))
  expect_true(q > 0 && q < .Machine$double.xmin)
  expect_lt(abs(pmhn(q, 0.8, 1, 3) / 1e-250 - 1), 1e-7)
  expect_silent(q <- qmhn(-362.5, 0.5, 1, -2, log.p = TRUE))
  expect_true(q > 0 && q < .Machine$double.xmin)
  expect_lt(abs(pmhn(q, 0.5, 1, -2, log.p = TRUE) + 362.5), 1e-7)
  expect_silent(q <- qmhn(c(1e-250, 0.5), c(0.8, 0.8), 1, -3))
  expect_lt(abs(pmhn(q[1], 0.8, 1, -3) / 1e-250 - 1), 1e-7)
})

test_that("qmhn reports a quantile below the smallest positive double as 0", {
  # F(x) behaves like x^alpha near the origin, so at alpha = 0.1 the quantile
  # of p = 1e-300 lies near 1e-3000, and at alpha = 0.5 that of log p = -800
  # near exp(-1600): F at 4.9e-324, the smallest positive double, already
  # exceeds p.  Such a quantile is 0 in double precision, as qgamma and the
  # special cases here (gamma = 0, last two lines, where the quantile is
  # sqrt(qgamma(p, alpha / 2) / beta)) report it.
  expect_identical(qmhn(1e-300, 0.1, 1, -2), 0)
  expect_identical(qmhn(c(1e-300, 1e-300), 0.1, c(1, 1), -2), c(0, 0))
  expect_identical(qmhn(-800, 0.5, 1, -2, log.p = TRUE), 0)
  expect_identical(qmhn(1e-300, 0.1, 1, 0), sqrt(stats::qgamma(1e-300, 0.05)))
  expect_identical(qmhn(-800, 0.5, 1, 0, log.p = TRUE), 0)
})

test_that("qmhn solves a survival probability near 1 on the lower tail", {
  # Near the lower end of the law log P(X > x) is about -F(x), and it carries
  # an absolute error of eps times the size of its terms (see the next test),
  # so solving log P(X > x) = log p there leaves F with a relative error of
  # that divided by F: all of F = 1e-20.  Once 1 - p is below 1e-7 qmhn solves
  # F(x) = 1 - p instead, with 1 - p formed as -expm1(log p) on the log scale
  # and by a subtraction that is exact for p >= 1/2 on the linear one.  Where
  # F comes from a quadrature, as at alpha = 412, the equation is solved
  # directly, and the answer is the lower-tail quantile of that probability,
  # to the bit.  Where it comes from the Lemma 1b series, as at the other two
  # shapes, the root of the survival equation (for 1 - p, or for a 1 - p the
  # upper tail still resolves) is refined by Newton steps on log F.  The two
  # are then roots of the same equation, found to a relative 1.4e-14 by the
  # root finder and to about 7e-15 by the Newton steps, and they agree to
  # within 1.6e-14 here; 1e-12 allows for rounding that differs between
  # platforms.
  shapes <- list(c(412, 0.0122, -42.6 * sqrt(0.0122)), c(2.5, 1, 1),
                 c(0.5, 1, 1))
  lp <- c(-1e-20, -1e-14, -1e-10)
  p <- 1 - c(1e-14, 1e-10)
  for (s in shapes) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", s[1], s[2], s[3])
    q <- qmhn(lp, s[1], s[2], s[3], lower.tail = FALSE, log.p = TRUE)
    q_lin <- qmhn(p, s[1], s[2], s[3], lower.tail = FALSE)
    if (s[1] > 100) {
      expect_identical(q, qmhn(-expm1(lp), s[1], s[2], s[3]), label = info)
      expect_identical(q_lin, qmhn(1 - p, s[1], s[2], s[3]), label = info)
    } else {
      expect_lt(max(abs(q / qmhn(-expm1(lp), s[1], s[2], s[3]) - 1)), 1e-12,
                label = info)
      expect_lt(max(abs(q_lin / qmhn(1 - p, s[1], s[2], s[3]) - 1)), 1e-12,
                label = info)
    }
    # The recycled-parameter path makes the same choice.
    expect_identical(qmhn(lp, c(s[1], s[1], s[1]), s[2], s[3],
                          lower.tail = FALSE, log.p = TRUE), q, label = info)
    # And the quantile is the right one: F(q) = 1 - p.  The root is solved
    # to a relative 1.4e-14 in q, which moves F by x f(x) / F(x) times that,
    # at most about 1e3 times here, so 1e-9 leaves ample room.
    expect_lt(max(abs(pmhn(q, s[1], s[2], s[3]) / -expm1(lp) - 1)), 1e-9,
              label = info)
    # At a survival of 1 - 1e-6 the two equations give the same quantile.
    # Wherever the survival equation would leave F with a relative error
    # above 1e-9, the quantile is taken from the lower tail instead (see the
    # next test); below that, x f(x) / F(x), about alpha next to the origin,
    # divides the error of F in the quantile.  1e-8 allows for either.
    q_up <- qmhn(-1e-6, s[1], s[2], s[3], lower.tail = FALSE, log.p = TRUE)
    q_lo <- qmhn(-expm1(-1e-6), s[1], s[2], s[3])
    expect_lt(abs(q_up / q_lo - 1), 1e-8, label = info)
  }
})

test_that("an upper-tail p close to 1 keeps the relative accuracy of 1 - p", {
  # log P(X > x) is assembled from terms as large as |log Psi| (with
  # gamma^2 / (4 beta) removed from it for z > 1000), (alpha / 2) |log beta|
  # and max(alpha, z) (1 + |log v*|), v* the peak of the density of log X,
  # and carries eps times their size or so.  Near the lower end of the law
  # it is about -F(x), so solving log P(X > x) = log p for p = 1 - d leaves
  # d with a relative error of that divided by d: at d = 1e-6 the quantile
  # would move by 9e-5 of the law's standard deviation at (9e5, 1, -500),
  # 7e-5 at (2e4, 1e-300, 5e-151), and 3e-6 and 9e-7 under the strong
  # positive tilts.  Where that error would pass 1e-9 of d, the quantile is
  # instead the root of F(x) = d, found directly where F comes from a
  # quadrature, and by polishing the survival root with Newton steps on
  # log F where it comes from the Lemma 1b series.  For alpha < 1 past
  # |z| = 20, where the series costs milliseconds, the polish takes log F
  # from the lower-tail quadrature, checked against the upper tail: at
  # (0.3, 1, 300) and (0.01, 1, 500) it serves d up to about 0.06 and 0.17.
  # References: the root of F(x) = 1 - (1 - d), the d that 1 - d leaves, by
  # Newton's method on the distribution function and the density from the
  # integral of the kernel on the axis log(sqrt(beta) x), Gauss-Legendre
  # panels in 160- or 256-bit arithmetic.  The root is found to a relative
  # 1.4e-14, which is below 2e-11 of a standard deviation at these points,
  # and F itself is accurate to about 1e-9 relative at the largest terms
  # (2e4, 1e-300, 5e-151), which moves the quantile by about 2e-10 of one;
  # 1e-9 allows for both.
  cases <- list(
    list(9e5, 1, -500, 1e-6, 555.2205961114197),
    list(9e5, 1, -500, 1e-4, 555.68739941599961),
    list(2e4, 1e-300, 5e-151, 1e-6, 9.7754697842415302e+151),
    list(0.5, 1, 900, 1e-6, 446.63826379501631),
    list(0.3, 1, 300, 1e-6, 146.63646149141815),
    list(0.3, 1, 300, 1e-4, 147.36790412220671),
    list(0.3, 1, 300, 0.015, 148.46316895618682),
    list(0.01, 1, 500, 1e-6, 246.63682797989424),
    list(0.01, 1, 500, 0.04, 248.76009306246561),
    list(50, 1, 20, 1.01e-7, 8.6975389402568908))
  for (cs in cases) {
    a <- cs[[1]]; b <- cs[[2]]; g <- cs[[3]]; d <- cs[[4]]
    info <- sprintf("alpha=%g beta=%g gamma=%g d=%g", a, b, g, d)
    sd <- sqrt(mhn_var(a, b, g))
    q <- qmhn(1 - d, a, b, g, lower.tail = FALSE)
    expect_lt(abs(q - cs[[5]]) / sd, 1e-9, label = info)
    # On the log scale the same d is -expm1(log1p(-d)), to a relative eps.
    q_log <- qmhn(log1p(-d), a, b, g, lower.tail = FALSE, log.p = TRUE)
    expect_lt(abs(q_log - cs[[5]]) / sd, 1e-9, label = info)
    # The recycled-parameter path makes the same choice.
    expect_identical(qmhn(1 - d, c(a, a), b, g, lower.tail = FALSE),
                     c(q, q), label = info)
  }
})

test_that("an upper-tail p within 1e-7 of 1 is refined on the lower tail", {
  # Within 1e-7 of 1 the quantile is the root of F(x) = d, d = 1 - p.  Where
  # F comes from the Lemma 1b series (alpha < 1 with |z| up to about 830,
  # beyond which the series is refused, or |z| <= 20 with alpha <= 100), the
  # survival equation is solved first, for d or for a d the upper tail still
  # resolves, and its root is refined by Newton steps on log F, each one
  # lower-tail evaluation and one density; a lower-tail solve would evaluate
  # F some 30 times.  The last case, at |z| = 1993, takes F from its
  # quadrature and is solved on it directly.  References: the root of
  # F(x) = 1 - (1 - d), the d that 1 - d leaves, by Newton's method on the
  # distribution function and the density from the integral of the kernel on
  # the axis log(sqrt(beta) x), Gauss-Legendre panels in 256-bit arithmetic.
  # qmhn meets them to within 1.3e-14 of the quantile.  The Newton steps end
  # below 7e-15 in log x, and 1e-12 allows for the rounding of F, which they
  # divide by x f(x) / F(x), about alpha next to the origin.
  cases <- list(
    list(0.3, 1, 300, 1e-10, 145.49948361136327),
    list(0.01, 1, 500, 1e-9, 245.75691297026273),
    list(0.3, 1, -300, 1e-12, 2.3240280815899407e-43),
    list(50, 1, 20, 1e-12, 7.5521768810053898),
    list(0.5, 1, -5, 1e-15, 1.487944096945259e-31),
    list(3, 1, 10, 1e-9, 1.1413351119516588),
    list(0.16048163993817749, 1, -1993.3449059212628, 1e-8,
         4.4985309173216283e-54))
  for (cs in cases) {
    a <- cs[[1]]; b <- cs[[2]]; g <- cs[[3]]; d <- cs[[4]]
    info <- sprintf("alpha=%.17g beta=%g gamma=%.17g d=%g", a, b, g, d)
    q <- qmhn(1 - d, a, b, g, lower.tail = FALSE)
    expect_lt(abs(q / cs[[5]] - 1), 1e-12, label = info)
    # The same d on the log scale, where -expm1(log1p(-d)) returns it to a
    # relative eps.
    d_eff <- 1 - (1 - d)
    q_log <- qmhn(log1p(-d_eff), a, b, g, lower.tail = FALSE, log.p = TRUE)
    expect_lt(abs(q_log / cs[[5]] - 1), 1e-12, label = info)
  }
  # For alpha = 0.01 F is about c x^0.01 next to the origin, and the
  # quantile of d = 1e-8 lies far below the smallest double, as the
  # lower-tail quantile of the same d reports.
  expect_identical(qmhn(1 - 1e-8, 0.01, 1, -300, lower.tail = FALSE), 0)
  expect_identical(qmhn(1 - (1 - 1e-8), 0.01, 1, -300), 0)
})

test_that("an upper-tail quantile at or below the smallest double is found", {
  # For alpha < 1 under a negative tilt F is c x^alpha next to the origin,
  # c = 2 beta^(alpha/2) / (alpha Psi), to within a relative |gamma| x, so
  # the quantile of a small d = 1 - p is (d alpha Psi / 2)^(1/alpha) /
  # sqrt(beta), and it can lie at or below the smallest positive double,
  # 4.9e-324.  The root of the survival equation, or a Newton step on log F
  # from it, then falls below that double; the steps stop there, and return
  # 0 where F there is still at least d, as the lower-tail quantile of d
  # does.
  #
  # At (0.0677, 1, -267.08) log F(4.9e-324) is -49.98 by that formula,
  # above log(1e-35) = -80.6, so the quantile is 0.
  expect_identical(qmhn(log1p(-1e-35), 0.0677, 1, -267.08,
                        lower.tail = FALSE, log.p = TRUE), 0)
  # At (0.3, 1, -300) the formula gives 2.3241994614667581e-53 and
  # 2.3241994614667469e-103 for d = 1e-15 and 1e-30, with log Psi from its
  # asymptotic expansion in 1 / z^2 summed in 320-bit arithmetic, and
  # quantiles below exp(-773) for d = 1e-100 and 1e-300, which are 0.  The
  # quantile moves by the relative error of F divided by alpha; 1e-12
  # allows for rounding of F up to 3e-13.
  q <- qmhn(log1p(-c(1e-15, 1e-30, 1e-100, 1e-300)), 0.3, 1, -300,
            lower.tail = FALSE, log.p = TRUE)
  expect_lt(abs(q[1] / 2.3241994614667581e-53 - 1), 1e-12)
  expect_lt(abs(q[2] / 2.3241994614667469e-103 - 1), 1e-12)
  expect_identical(q[3:4], c(0, 0))
  # At (0.0219, 212.77, -5436.19) F is 1.0060e-7 at 4.9e-324 and 1.0214e-7
  # at the next double, 9.9e-324, so the quantile of 1.01e-7 lies between
  # them (at about 5.9e-324 by the power law); either is the answer to the
  # spacing of the doubles there.
  q <- qmhn(1 - 1.01e-7, 0.021912835747956237, 212.7699360703501,
            -5436.1925609256077, lower.tail = FALSE)
  expect_true(q %in% c(5e-324, 1e-323))
})

test_that("an upper-tail quantile at the smallest double costs less than a CDF", {
  skip_on_cran()
  # Where the Newton steps above would give up instead, the lower-tail
  # equation would be solved from scratch, halving its bracket hundreds of
  # times at one sum of the series each: 4 to 22 s per value at these
  # triples.  The steps take a few milliseconds, less than one evaluation of
  # the CDF at the median of the same triple.  The time is judged against
  # those evaluations, timed alongside, rather than against a number of
  # seconds, so that a slow or loaded machine stretches both; the bound of 20
  # of them is about 140 times what the steps take and a tenth of what the
  # bisection would.
  triples <- list(c(0.0677, 1, -267.08), c(0.545, 1, -559.68),
                  c(0.021912835747956237, 212.7699360703501,
                    -5436.1925609256077))
  medians <- vapply(triples, function(t) qmhn(0.5, t[1], t[2], t[3]),
                    numeric(1))
  time_of <- function(f) min(vapply(1:3, function(i)
    system.time(f())[["elapsed"]], numeric(1)))
  t_q <- time_of(function() {
    qmhn(log1p(-1e-35), 0.0677, 1, -267.08, lower.tail = FALSE, log.p = TRUE)
    qmhn(log1p(-1.19e-294), 0.545, 1, -559.68, lower.tail = FALSE,
         log.p = TRUE)
    qmhn(1 - 1.01e-7, 0.021912835747956237, 212.7699360703501,
         -5436.1925609256077, lower.tail = FALSE)
  })
  t_p <- time_of(function() {
    for (i in seq_along(triples)) {
      pmhn(medians[i], triples[[i]][1], triples[[i]][2], triples[[i]][3])
    }
  })
  expect_lt(t_q, 20 * t_p)
})

test_that("qmhn inverts the law at a huge shape and tilt", {
  # At alpha = 1e10, beta = 1, gamma = 1e10 the law is normal with mode
  # 5e9 + 1 and standard deviation 1/sqrt(2 + 1e10 / 2.5e19), 0.7071067,
  # to within a skewness of order (alpha - 1) / m^3, 8e-20.  The doubles
  # near the mode are 9.5e-7 apart, and the root is found to within one of
  # those spacings.
  lp <- c(-1, -1e-3)
  q <- qmhn(lp, 1e10, 1, 1e10, log.p = TRUE)
  q_normal <- 5e9 + 1 + qnorm(lp, log.p = TRUE) / sqrt(2 + 1e10 / 2.5e19)
  expect_lt(max(abs(q - q_normal)), 2e-6)
})

test_that("qmhn inverts the limit laws at tilts past the overflow of z^2", {
  # Under a negative tilt beyond |z| = 1.34e154 the law is Gamma(alpha, rate
  # |gamma|) to within a relative 1e-300 (see the matching pmhn test), so its
  # quantiles are those of qgamma divided by |z|.  Under a positive one it is
  # normal about z/2 + (alpha - 1)/z with standard deviation 1/sqrt(2), far
  # below the spacing of the doubles there, and every quantile strictly
  # inside (0, 1) lies within 3.1 of z/2, so that z/2 is the double nearest
  # to it.  The distribution function is then a step at z/2, which the
  # root finder follows to adjacent doubles, returning the nearer.  1e-8 is
  # the working tolerance.
  p <- c(0.001, 0.1, 0.5, 0.9, 0.999)
  for (a in c(0.5, 2, 10)) {
    for (z in c(-1e155, -1e200, -1e300)) {
      expect_equal(qmhn(p, a, 1, z) * abs(z), qgamma(p, a), tolerance = 1e-8,
                   info = sprintf("alpha=%g z=%g", a, z))
    }
    for (z in c(1e155, 1e200, 1e300)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      expect_identical(qmhn(p, a, 1, z), rep(z / 2, 5), label = info)
      expect_identical(qmhn(p, a, 1, z, lower.tail = FALSE), rep(z / 2, 5),
                       label = info)
    }
  }
})

test_that("qmhn resolves a law narrower than the doubles to the nearest double", {
  # Under a tilt z = gamma / sqrt(beta) the law is normal about its peak with
  # standard deviation 1 / sqrt(2 beta) to within relative corrections of
  # order (alpha - 1) / z^2, so its quantiles are
  # m + qnorm(p) / sqrt(2 beta), m = gamma / (2 beta) + (alpha - 1) / gamma.
  # Where that standard deviation is below the spacing of the doubles, the
  # distribution function is a step from 0 to 1 between two adjacent
  # doubles, on which interpolation stalls and the log tails reach the
  # largest doubles; the quantile is then the double nearest to it.  At
  # z = 1e150 the spacing near the mode, 5e149, is 1e134 and every quantile
  # is 5e149; at (2, 1e-300, 1e-50) the mode is 5e249.
  for (p in c(0.1, 0.5, 0.9)) {
    expect_identical(qmhn(p, 2, 1, 1e150, lower.tail = FALSE), 5e149)
    expect_identical(qmhn(p, 2, 1e-300, 1e100 * sqrt(1e-300),
                          lower.tail = FALSE),
                     1e100 * sqrt(1e-300) / 2e-300)
  }
  # Near the largest double the bracket is doubled only up to it, and the
  # midpoint of a bracket is formed without overflow.  The modes here are
  # DBL_MAX / 2 and 0.7 DBL_MAX / 0.8, within a relative 1e-307 (the spacing
  # there is 2e292); qnorm(0.9) / sqrt(0.8) is nothing next to it.
  X <- .Machine$double.xmax
  expect_identical(qmhn(0.5, 0.5, 1, X), X / 2)
  expect_identical(qmhn(0.5, 10, 1, X), X / 2)
  expect_identical(qmhn(0.9, 2, 0.4, 0.7 * X, lower.tail = FALSE),
                   0.7 * X / 0.8)
  # Where the whole law lies beyond the largest double, so does every
  # quantile strictly inside (0, 1): (2, 1e-100, DBL_MAX) has its mode at
  # 9e407.
  expect_identical(qmhn(c(0.1, 0.9), 2, 1e-100, X), c(Inf, Inf))
})

test_that("qmhn keeps its accuracy as the law narrows toward the spacing of the doubles", {
  # As in the previous test the quantiles at beta = 1 are
  # m + qnorm(p) / sqrt(2), m = z / 2 + (alpha - 1) / z, to within relative
  # corrections of order (alpha - 1) / z^2, far below the spacing of the
  # doubles near m, which is 2^-53 of it: 9e-10 at z = 1e7, 0.0625 at
  # z = 1e15, where it is 0.09 of the standard deviation.  qmhn finds the
  # root to within 2^-46 of the standard deviation or, where that is below
  # the spacing, to adjacent doubles, and returns the nearer; the reference,
  # formed in double precision, is within a spacing of the true value, so the
  # two are within two spacings.  A lower-tail p near 1, where F is resolved
  # only to eps, which would move the quantile by eps / f (180 spacings at
  # p = 1 - 1e-10 and z = 1e7), is solved on the upper tail 1 - p.
  p <- c(1e-10, 0.1, 0.5, 0.9, 1 - 1e-10)
  for (a in c(0.5, 2, 10)) {
    for (z in c(1e7, 1e10, 1e15)) {
      info <- sprintf("alpha=%g z=%g", a, z)
      m <- z / 2 + (a - 1) / z
      spacing <- 2^(floor(log2(z / 2)) - 52)
      expect_lte(max(abs(qmhn(p, a, 1, z) - (m + qnorm(p) / sqrt(2)))),
                 2 * spacing, label = info)
      expect_lte(max(abs(qmhn(p, a, 1, z, lower.tail = FALSE) -
                           (m + qnorm(p, lower.tail = FALSE) / sqrt(2)))),
                 2 * spacing, label = info)
    }
  }
})

test_that("qmhn at gamma = 0 and a huge shape returns the double nearest the quantile", {
  skip_if_not_installed("Rmpfr")
  # At gamma = 0 the law is normal about v* = sqrt(alpha / (2 beta)) with
  # standard deviation sd = v* / sqrt(2 alpha), to within relative
  # corrections of order 1 / sqrt(alpha), so its quantiles are
  # v* + sd qnorm(p) to within about 3e-13 sd at alpha = 1e25.  sd is 1e3
  # spacings of the doubles near v* there, and far below one spacing at
  # 1e35 and 1e100.  qmhn finds the root to within 2^-46 sd or to adjacent
  # doubles and returns the nearer, so it lies within half a spacing of the
  # reference, formed in Rmpfr; 1e-6 of a spacing allows for the rest.
  p <- c(1e-10, 0.1, 0.5, 0.9)
  for (a in c(1e25, 1e35, 1e100)) {
    for (b in c(1e-300, 1, 1e300)) {
      info <- sprintf("alpha=%g beta=%g", a, b)
      V <- sqrt(Rmpfr::mpfr(a, 512) / (2 * Rmpfr::mpfr(b, 512)))
      sd <- V / sqrt(Rmpfr::mpfr(2 * a, 512))
      spacing <- 2^(floor(log2(Rmpfr::asNumeric(V))) - 52)
      err <- function(q, z) {
        max(Rmpfr::asNumeric(abs(Rmpfr::mpfr(q, 512) - (V + sd * z))))
      }
      expect_lte(err(qmhn(p, a, b, 0), qnorm(p)), 0.500001 * spacing,
                 label = info)
      expect_lte(err(qmhn(p, a, b, 0, lower.tail = FALSE),
                     qnorm(p, lower.tail = FALSE)), 0.500001 * spacing,
                 label = info)
    }
  }
})

test_that("qmhn solves a lower-tail p near 1 on the upper tail", {
  # F near the upper end of the law is resolved only to the rounding of 1,
  # and the root of F(x) = p would move by eps / f(x): a relative 1e-6 at
  # p = 1 - 1e-15 and alpha = 1e6 + 1.  Close to 1 (within 1e-7 of 1 at
  # least) a lower-tail p is solved as the upper tail 1 - p, exact for
  # p >= 1/2 (and -expm1(p) on the log scale), on the directly computed
  # log P(X > x).  At gamma = 0 the
  # quantile is sqrt(y / beta) with y the Gamma(alpha / 2) quantile of that
  # upper tail, from qgamma; the root is found to a relative 2^-46, so 1e-13
  # allows for it and for qgamma.
  a <- 1e6 + 1
  u <- 2^-40
  ref <- sqrt(stats::qgamma(u, a / 2, lower.tail = FALSE))
  expect_lt(abs(qmhn(1 - u, a, 1, 0) / ref - 1), 1e-13)
  expect_lt(abs(qmhn(c(1 - u, 1 - u), c(a, a), 1, 0)[2] / ref - 1), 1e-13)
  ref_log <- sqrt(stats::qgamma(-expm1(-1e-15), a / 2, lower.tail = FALSE))
  expect_lt(abs(qmhn(-1e-15, a, 1, 0, log.p = TRUE) / ref_log - 1), 1e-13)
  # In the general case the answer is the upper-tail quantile of 1 - p, to
  # the bit, and P(X > q) is 1 - p to the root's tolerance times
  # q f(q) / P(X > q), about 100 here, so 1e-10 leaves ample room.
  for (s in list(c(2.5, 1, 1), c(0.3, 1, -5), c(412, 1, 3))) {
    info <- sprintf("alpha=%g beta=%g gamma=%g", s[1], s[2], s[3])
    q <- qmhn(1 - u, s[1], s[2], s[3])
    expect_identical(q, qmhn(u, s[1], s[2], s[3], lower.tail = FALSE),
                     label = info)
    expect_identical(qmhn(log1p(-u), s[1], s[2], s[3], log.p = TRUE),
                     qmhn(-expm1(log1p(-u)), s[1], s[2], s[3],
                          lower.tail = FALSE), label = info)
    expect_lt(abs(pmhn(q, s[1], s[2], s[3], lower.tail = FALSE) / u - 1),
              1e-10, label = info)
  }
})

test_that("qmhn past |z| = 20 or alpha = 100 keeps the working tolerance", {
  # There the lower tail is taken by quadrature for a shape of 1 or more and
  # from the series below 1 (see the matching tests in test-pmhn.R).
  # References: the root of the reference distribution
  # function there, the integral of the law of log X measured from its peak.
  # A relative error e in F moves the median by about e / (2 x f(x)) of
  # itself, below 2 e at each of these points, so 1e-8 is the working
  # tolerance.
  a <- c(1.5, 0.3, 20, 0.7, 3000, 2.5)
  g <- c(300, 21, -300, -300, 5, -30)
  ref <- c(150.00166667283986, 10.466456608368206, 0.065529311613451879,
           0.0013580473630783359, 39.995629466795208, 0.072073331888236017)
  for (i in seq_along(a)) {
    expect_lt(abs(qmhn(0.5, a[i], 1, g[i]) / ref[i] - 1), 1e-8,
              label = sprintf("alpha=%g gamma=%g", a[i], g[i]))
  }
})

test_that("qmhn of a law narrower than the doubles holds at a subnormal beta", {
  # Scaling by a power of two is exact, so the quantiles at beta = 4^-k and
  # gamma = z 2^-k are 2^k times those at beta = 1, to the last bit.  At
  # k = 532 the standard deviation, about 2^532, has a square beyond the
  # largest double.  A comparison to the last bit, so not run on CRAN.
  skip_on_cran()
  p <- c(1e-10, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  for (case in list(c(2, 1e20), c(1e40, 1))) {
    alpha <- case[1]
    z <- case[2]
    at_one <- qmhn(p, alpha, 1, z)
    for (k in c(500, 532)) {
      expect_identical(qmhn(p, alpha, 4^-k, z * 2^-k) / 2^k, at_one)
    }
  }
})

test_that("qmhn at gamma = 0 resolves a far upper tail on the lower log scale", {
  # At gamma = 0, X^2 ~ Gamma(alpha/2, beta) (Sun et al. 2023, Lemma 6a), and
  # the quantile is taken from qgamma.  R's qgamma loses the far upper tail
  # when it is given as a lower-tail log probability close to 0: for alpha =
  # 1e-10 and log p = log(1 - 1e-300) the quantile came out 10.92, where it is
  # 25.70, and at alpha = 10 0.2 to 0.3 percent off.  The upper tail it is
  # given now is checked against pgamma, and the quantile against the one the
  # upper-tail form returns.
  for (a in c(1e-10, 0.1, 5, 10, 1000)) {
    for (Q in c(1e-300, 1e-100, 1e-50, 1e-20)) {
      info <- sprintf("alpha = %g, Q = %g", a, Q)
      x <- qmhn(log1p(-Q), a, 1, 0, log.p = TRUE)
      expect_equal(pgamma(x^2, a / 2, lower.tail = FALSE, log.p = TRUE),
                   log(Q), tolerance = 1e-10, info = info)
      expect_equal(x, qmhn(Q, a, 1, 0, lower.tail = FALSE), tolerance = 1e-12,
                   info = info)
      expect_equal(qmhn(log1p(-Q), a, 4, 0, log.p = TRUE), x / 2,
                   tolerance = 1e-12, info = info)
    }
  }
})

test_that("qmhn at gamma = 0 resolves a far lower tail on the upper log scale", {
  # The counterpart of the block above: R's qgamma loses the far lower tail
  # when it is given as an upper-tail log probability close to 0.  At
  # alpha = 500 and log p = -1e-300 the quantile came out 2.4267, where it is
  # 2.4561 and F there 2.8e-303 rather than 1e-300, and at alpha = 300 and
  # F = 1e-100 it was 8.4 percent off.  The lower tail it is given now is
  # checked against pgamma, and the quantile against the one the lower-tail
  # form returns.
  for (a in c(5, 300, 1000, 1e4)) {
    for (F in c(1e-300, 1e-100, 1e-50, 1e-20)) {
      info <- sprintf("alpha = %g, F = %g", a, F)
      x <- qmhn(log1p(-F), a, 1, 0, lower.tail = FALSE, log.p = TRUE)
      expect_equal(pgamma(x^2, a / 2, log.p = TRUE), log(F),
                   tolerance = 1e-10, info = info)
      expect_equal(x, qmhn(F, a, 1, 0), tolerance = 1e-12, info = info)
      expect_equal(qmhn(log1p(-F), a, 4, 0, lower.tail = FALSE, log.p = TRUE),
                   x / 2, tolerance = 1e-12, info = info)
    }
  }
})

test_that("a lower-tail p close to 1 is solved on the upper tail", {
  # F near 1 is resolved only to the rounding of 1 and, for a small shape,
  # to the error of the lower tail, a relative error of that divided by
  # 1 - p: in 0.1.1 the quantile of log1p(-1e-6) at (1e-8, 1, 0.5) was
  # 1.3e-7 of itself off.  Wherever that could pass 1e-9 of 1 - p and the
  # upper tail resolves it better, the quantile is the root of
  # P(X > x) = 1 - p, so that it is the one the upper-tail call returns.
  for (p in list(c(1e-8, 0.5), c(1e-6, 0.5), c(1e-10, 5), c(1e-6, -10),
                 c(2, 0.5), c(0.5, 1))) {
    a <- p[1]; z <- p[2]
    d <- c(1e-7, 1e-6, 3e-6)
    info <- sprintf("alpha=%g z=%g", a, z)
    expect_close(qmhn(log1p(-d), a, 1, z, log.p = TRUE),
                 qmhn(log(d), a, 1, z, lower.tail = FALSE, log.p = TRUE),
                 1e-13, info = info)
    p_lin <- 1 - d
    expect_close(qmhn(p_lin, a, 1, z), qmhn(1 - p_lin, a, 1, z,
                                             lower.tail = FALSE),
                 1e-13, info = info)
  }
  # At alpha = 2 the upper tail has a closed form, in which nothing cancels:
  # P(X > x) = (exp(-x^2 + z x) / 2 + (z / 2) e^{m^2} sqrt(pi)
  # Phi(-sqrt(2) (x - m))) / Z, m = z / 2.  The quantile of 1 - p = 1e-7 at
  # (2, 1, 0.5) left it 2.6e-9 away from 1 - p in 0.1.1, and 1e-14 now.
  upper2 <- function(x, z) {
    m <- z / 2
    num <- exp(-x^2 + z * x) / 2 +
      (z / 2) * exp(m^2) * sqrt(pi) * pnorm(-sqrt(2) * (x - m))
    num / (1 / 2 + (z / 2) * exp(m^2) * sqrt(pi) * pnorm(sqrt(2) * m))
  }
  for (z in c(0.5, 5)) {
    for (d in c(1e-7, 1e-6, 3e-6)) {
      q <- qmhn(log1p(-d), 2, 1, z, log.p = TRUE)
      expect_lt(abs(upper2(q, z) / d - 1), 1e-12,
                label = sprintf("P(X > q) / (1 - p) - 1 at z=%g d=%g", z, d))
    }
  }
})
