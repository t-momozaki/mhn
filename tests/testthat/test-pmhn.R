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
    F_q <- pmhn(q, a, b, g)
    # Restrict monotonicity to the meaningful body of the CDF; near F = 1
    # the series and integration fallback differ in their last digits
    # (~1e-10), which is expected floating-point noise.
    idx <- F_q < 1 - 1e-9
    expect_true(all(diff(F_q[idx]) >= -1e-10),
                label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
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
    expect_equal(logp[idx], log(linear[idx]), tolerance = 1e-12,
                 label = sprintf("a=%.1f b=%.1f g=%+.1f", a, b, g))
    # Upper-tail log.  The direct path (lower.tail=FALSE, log.p=TRUE) is
    # often more accurate than log1p(-linear) because the latter loses
    # precision when linear is close to 1; allow a moderate tolerance.
    upper_logp <- pmhn(q, a, b, g, lower.tail = FALSE, log.p = TRUE)
    idx2 <- linear < 1 - 1e-10
    expect_equal(upper_logp[idx2], log1p(-linear[idx2]), tolerance = 1e-6,
                 label = sprintf("upper a=%.1f b=%.1f g=%+.1f", a, b, g))
  }
})

# ============================================================
# 5. pmhn ~= integrate(dmhn, 0, q)
# ============================================================

test_that("pmhn matches integrate(dmhn, 0, q)", {
  params <- expand.grid(alpha = c(0.5, 1.5, 2, 3),
                        beta  = c(0.5, 1, 2),
                        gamma = c(-3, -1, 1, 3))
  q_vals <- c(0.1, 0.5, 1, 2, 4)
  for (i in seq_len(nrow(params))) {
    a <- params$alpha[i]; b <- params$beta[i]; g <- params$gamma[i]
    for (q in q_vals) {
      ref <- tryCatch(
        integrate(function(x) dmhn(x, a, b, g), 0, q,
                  rel.tol = 1e-10, subdivisions = 300)$value,
        error = function(e) NA
      )
      if (is.na(ref)) next
      got <- pmhn(q, a, b, g)
      expect_equal(got, ref, tolerance = 1e-6,
                   label = sprintf("a=%.1f b=%.1f g=%+.1f q=%.2f", a, b, g, q))
    }
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

test_that("pmhn errors on invalid parameters", {
  expect_error(pmhn(1, alpha = 0, beta = 1, gamma = 0), "alpha")
  expect_error(pmhn(1, alpha = -1, beta = 1, gamma = 0), "alpha")
  expect_error(pmhn(1, alpha = 1, beta = 0, gamma = 0), "beta")
  expect_error(pmhn(1, alpha = 1, beta = 1, gamma = Inf), "gamma")
  expect_error(pmhn(1, alpha = NA, beta = 1, gamma = 0), "alpha")
})

# ============================================================
# 9. Series-vs-integration cross-check across the accuracy grid
# ============================================================
#
# Loads the most recent `cdf_series_accuracy_<YYYYMMDD>.csv` shipped
# under inst/audits/results/ and asserts that the dispatch
# (`pmhn(q, alpha, beta, gamma)`) agrees with the corresponding
# Boost.Math `gauss_kronrod` / `tanh_sinh` integration result
# (`val_integrate`) across every grid point.  The integration baseline
# is the same independent path the dispatcher itself falls back to
# when the Lemma 1b series triggers its cancellation guard, so
# matching it is the strongest in-package regression we can run.
#
# The mpfr reference in the same CSV (`val_ref`) is intentionally not
# used here: for the most extreme (alpha = 10, |gamma| >= 15) corner
# the mpfr-Psi via Lemma 10 series itself loses precision, occasionally
# producing |F_mpfr| > 1.  Those rows would cause spurious failures
# even though the dispatch result is mathematically correct.

test_that("pmhn matches the Boost.Math integration across the cdf_series_accuracy grid", {
  csv_dir <- system.file("audits", "results", package = "mhn")
  skip_if(csv_dir == "",
          "cdf_series_accuracy CSV is not present in this installation")
  csv_files <- list.files(
    csv_dir,
    pattern = "^cdf_series_accuracy_\\d{8}\\.csv$",
    full.names = TRUE
  )
  skip_if(length(csv_files) == 0L,
          "cdf_series_accuracy CSV is not present in this installation")

  csv_path <- tail(sort(csv_files), 1L)
  df <- utils::read.csv(csv_path)

  # Skip rows where F is well past the bulk of the distribution
  # (|F| < 1e-10 or |1 - F| < 1e-10).  In those deep tails the
  # truncated-normal closed form and the Boost.Math integration disagree
  # not because either is wrong but because each tail probability sits
  # at the double-precision floor where rounding noise dominates the
  # comparison.
  skip_row <- is.na(df$val_integrate) |
                abs(df$val_integrate) < 1e-10 |
                abs(1 - df$val_integrate) < 1e-10
  df <- df[!skip_row, ]

  # Tolerance 1e-6 reflects the practical agreement between the Sun et
  # al. (2023, Lemma 1b) series and the Boost.Math integration of the
  # unnormalised density at double precision -- for the high-alpha,
  # moderate-|gamma| corner the two independent paths intrinsically
  # differ by O(1e-7).  The pre-audit silent bugs were 1e-1 to 1 in
  # magnitude, so 1e-6 still catches every meaningful regression.
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    val <- pmhn(r$q, alpha = r$alpha, beta = r$beta, gamma = r$gamma)
    expect_lt(
      abs(val - r$val_integrate) / max(abs(r$val_integrate), 1e-300),
      1e-6,
      label = sprintf("pmhn(q=%g, alpha=%g, beta=%g, gamma=%g) vs val_integrate",
                      r$q, r$alpha, r$beta, r$gamma)
    )
  }
})
