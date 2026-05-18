#!/usr/bin/env Rscript
# cdf_series_accuracy.R --- Accuracy audit of the pmhn() Lemma 1b series.
#
# Role: this is an *accuracy audit*, not a timing benchmark.  It
# verifies that the dispatcher in src/mhn_pmhn_state.h produces
# arbitrary-precision-correct CDF values across a wide (alpha, gamma, q)
# grid.  The companion file series_breakdown_audit.R verifies the
# underlying theoretical breakdown behaviour of the Lemma 10 Psi
# series and the Lemma 1b CDF series; this one verifies the pmhn()
# output itself.
#
# Background: pmhn() uses the Sun et al. (2023, Lemma 1b) CDF series
# as its primary path.  For gamma < 0 the series has alternating-sign
# terms whose log-sum-exp accumulator can suffer catastrophic
# cancellation as |gamma / sqrt(beta)| grows: lgamma rounding errors
# accumulate to a magnitude comparable to Psi, silently corrupting the
# CDF.  The implementation in src/mhn_cdf_series.cpp guards against
# this at runtime by tracking the cancellation depth and returning
# NaN whenever the relative error would exceed the user's tolerance,
# in which case the dispatcher in src/mhn_pmhn_state.h falls back to
# a Boost.Math integration of the unnormalised density.
#
# Purpose: at each point of a grid that spans the suspected breakdown
# region for gamma < 0, call mhn:::.pmhn_force(..., method = "series")
# and .pmhn_force(..., method = "integrate") in isolation (bypassing
# the dispatcher), then compare each value against an arbitrary-
# precision reference computed from the same Lemma 1b series in
# `Rmpfr` at 200 bits of precision.  The dispatcher's default output
# is also recorded.  The audit verifies that the runtime cancellation
# guard fires correctly across the grid and that the dispatcher's
# output matches the high-precision reference within tolerance.
#
# Outputs:
#   mhn/inst/audits/results/cdf_series_accuracy_<YYYYMMDD>.csv
#   mhn/inst/audits/results/cdf_series_accuracy_diagnostics_<YYYYMMDD>.csv
#
# Invocation (from repository root):
#   Rscript mhn/inst/audits/cdf_series_accuracy.R
#
# Environment variables (optional):
#   MHN_BENCH_QUICK=1     reduce the grid to a smoke test (~30 s)
#   MHN_BENCH_OUTDIR=DIR  override the result directory
#   MHN_BENCH_PREC=N      override the mpfr precision (default 200 bits)

suppressPackageStartupMessages({
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("cdf_series_accuracy.R requires the 'mhn' package. ",
         "Run R CMD INSTALL mhn first.")
  if (!requireNamespace("Rmpfr", quietly = TRUE))
    stop("cdf_series_accuracy.R requires the 'Rmpfr' package. ",
         "Install via install.packages('Rmpfr').")
})

QUICK <- isTRUE(nchar(Sys.getenv("MHN_BENCH_QUICK")) > 0L)
PREC_BITS <- {
  v <- suppressWarnings(as.integer(Sys.getenv("MHN_BENCH_PREC", unset = "")))
  if (is.na(v) || v < 64L) 200L else v
}
# Hard ceiling on the mpfr reference series.  In practice the series
# converges to mpfr machine zero (relative magnitude < 2^-PREC_BITS) in
# fewer than 4 * z^2 terms for the worst grid point (|z| = 100), well
# below this limit.
MAX_TERMS <- 30000L

ALPHAS <- if (QUICK) c(0.5, 2) else c(0.3, 0.7, 1, 1.5, 3, 10)
GAMMAS <- if (QUICK) {
  c(-25, -20, -10, -2, 2, 10, 20, 25)
} else {
  c(-100, -50, -25, -22, -20, -18, -15, -10, -5, -2, -0.5,
    0.5, 2, 5, 10, 15, 18, 20, 22, 25, 50, 100)
}
# Three q values per (alpha, beta, gamma) at fractions of the closed-
# form upper bound on E[X] from Sun et al. (2023) Lemma 4a.  Avoid
# qmhn() for q because qmhn calls pmhn internally and would propagate
# any pmhn bug into the choice of evaluation point.
Q_FACTORS <- c(1 / 3, 1, 3)
BETA <- 1.0

# Sun et al. (2023) Lemma 4a upper bound on E[X], proven for alpha > 1.
# Used here purely as a scale to pick informative q points; for
# alpha <= 1 we extrapolate the same formula since the audit only
# needs a reasonable positive scale, not a rigorous bound.  Strictly
# positive for every (alpha > 0, beta > 0, gamma in R).
approx_mean <- function(alpha, beta, gamma) {
  (gamma + sqrt(gamma^2 + 8 * alpha * beta)) / (4 * beta)
}

OUTDIR <- Sys.getenv("MHN_BENCH_OUTDIR", unset = "")
if (!nzchar(OUTDIR)) OUTDIR <- file.path("mhn", "inst", "audits", "results")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("cdf_series_accuracy_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR,
                        sprintf("cdf_series_accuracy_diagnostics_%s.csv", TODAY))

cat(sprintf("[cdf_series_accuracy] mode=%s prec=%d bits alpha=%d gamma=%d q=%d points=%d\n",
            if (QUICK) "QUICK" else "FULL", PREC_BITS,
            length(ALPHAS), length(GAMMAS), length(Q_FACTORS),
            length(ALPHAS) * length(GAMMAS) * length(Q_FACTORS)))
cat(sprintf("[cdf_series_accuracy] result_csv=%s\n", RESULT_CSV))
cat(sprintf("[cdf_series_accuracy] diag_csv=%s\n", DIAG_CSV))

# -----------------------------------------------------------------------
# High-precision reference: Sun et al. (2023) Lemma 1b numerator
# accumulated under Rmpfr arithmetic at `prec_bits` of precision, then
# divided by the package's well-conditioned log Psi.
#
# Lemma 1b (rearranged so each iteration shares one gamma evaluation):
#   F(x | alpha, beta, gamma)
#     = (1 / Psi) sum_{i = 0..inf}  Gamma((alpha + i)/2) (z^i / i!)
#                                   * P((alpha + i)/2, beta x^2)
# with z = gamma / sqrt(beta) and P the regularised lower incomplete
# gamma.  The numerator sum can be evaluated accurately in mpfr because
# the P_i factor exponentially suppresses high-i terms (so cancellation
# depth is bounded by log10(peak |T_i|) - log10(F * Psi), which stays
# moderate for the grid we test).
#
# Psi[alpha/2, z] itself is **not** accumulated from the Lemma 1b
# series here.  For gamma < 0 that sum has hundreds of decimal digits
# of cancellation between alternating-sign terms (|z| = 100 needs
# ~ 333 digits, more than reasonable mpfr precision can deliver
# efficiently), so we instead pull log Psi from the package's own
# dispatcher, which uses the Sun et al. (2023, Lemma 11) integral for
# gamma < 0 and the Sun et al. (2023, Lemma 10) positive-term series
# for gamma > 0 -- both well-
# conditioned.  This mirrors what the C++ implementation in
# src/mhn_cdf_series.cpp does internally; it is the numerator sum, not
# Psi, that the audit really exercises.
#
# Gamma((alpha + i)/2) is recurred in place via the shift-by-2 identity
# Gamma((alpha + i + 2)/2) = ((alpha + i)/2) * Gamma((alpha + i)/2);
# we maintain two interleaved sequences for even and odd i so each
# iteration is O(1) gamma updates (no fresh Gamma call beyond the two
# bootstrap values).
#
# The regularised P((alpha + i)/2, beta x^2) is evaluated in *double
# precision* via base R's pgamma(), then promoted to mpfr.  This
# preserves at least ~15 digits per P_i, far more than the 1e-8
# accuracy threshold we test against; the full mpfr precision is
# reserved for the outer numerator sum.  Rmpfr::pgamma() at 200 bits
# is several orders of magnitude slower than base pgamma() and was
# unfit for the full grid.
# -----------------------------------------------------------------------
log_cdf_series_mpfr <- function(alpha, beta, gamma, x,
                                prec_bits = PREC_BITS,
                                max_terms = MAX_TERMS,
                                rel_tol = 2^-(prec_bits - 16L)) {
  log_psi_d <- mhn:::.mhn_log_normalizing_const(alpha, beta, gamma, -1)
  log_psi   <- Rmpfr::mpfr(log_psi_d, prec_bits)

  ma <- Rmpfr::mpfr(alpha, prec_bits)
  mb <- Rmpfr::mpfr(beta,  prec_bits)
  mg <- Rmpfr::mpfr(gamma, prec_bits)
  mx <- Rmpfr::mpfr(x,     prec_bits)
  z  <- mg / sqrt(mb)
  y_mpfr   <- mb * mx^2
  y_double <- as.double(y_mpfr)

  zero <- Rmpfr::mpfr(0, prec_bits)
  one  <- Rmpfr::mpfr(1, prec_bits)
  half <- one / 2

  # Bootstrap the two parities of Gamma((alpha + n)/2).
  g_even <- gamma(ma / 2)            # n = 0, 2, 4, ...
  g_odd  <- gamma(ma / 2 + half)     # n = 1, 3, 5, ...

  zfact  <- one                      # z^n / n!
  S_num  <- zero                     # sum_{i} coeff_i * P_i  (numerator)
  max_term_abs <- zero               # for cancellation-aware termination

  reached_max <- TRUE
  iter <- 0L
  peak_index <- 0L
  # The full numerator term |coeff_i * P_i| typically peaks at modest i
  # (where Gamma * z^i / i! still grows but the P_i damping has not yet
  # taken hold), then decays super-geometrically.  Wait 200 iterations
  # past the running peak before letting the rel_tol cap end the loop;
  # this is robust to a transient lull before the dominant term.
  PEAK_BUFFER <- 200L
  for (i in seq.int(0, max_terms)) {
    is_even <- (i %% 2L == 0L)
    g_n     <- if (is_even) g_even else g_odd
    arg_i   <- (ma + i) / 2
    arg_i_d <- (alpha + i) / 2       # same value at double precision
    P_i     <- Rmpfr::mpfr(stats::pgamma(y_double, shape = arg_i_d),
                           prec_bits)

    coeff      <- g_n * zfact
    full_term  <- coeff * P_i
    S_num      <- S_num + full_term

    term_abs <- abs(full_term)
    if (term_abs > max_term_abs) {
      max_term_abs <- term_abs
      peak_index   <- i
    }

    # Recurrence for next iteration of the same parity.
    if (is_even) g_even <- g_even * arg_i else g_odd <- g_odd * arg_i
    zfact <- zfact * z / (i + 1L)

    if (i >= peak_index + PEAK_BUFFER &&
        term_abs < max_term_abs * rel_tol) {
      reached_max <- FALSE
      iter <- i
      break
    }
  }
  if (reached_max) iter <- max_terms

  # F = S_num / Psi.  Compute in log space for safe range; bounce
  # through log(|S_num|) so the formula does not blow up on subnormal
  # mpfr intermediates.
  if (as.numeric(S_num) == 0) {
    F_val <- 0.0
  } else {
    log_F <- log(abs(S_num)) - log_psi
    sgn   <- if (S_num < 0) -1.0 else 1.0
    F_val <- sgn * as.double(exp(log_F))
  }
  attr(F_val, "iterations") <- iter
  attr(F_val, "reached_max") <- reached_max
  F_val
}

# -----------------------------------------------------------------------
# Grid + main loop
# -----------------------------------------------------------------------
grid <- expand.grid(alpha = ALPHAS, gamma = GAMMAS, q_factor = Q_FACTORS,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
N <- nrow(grid)

rows <- vector("list", N)
t0 <- Sys.time()
progress_every <- max(1L, N %/% 20L)

for (k in seq_len(N)) {
  alpha    <- grid$alpha[k]
  gamma    <- grid$gamma[k]
  q_factor <- grid$q_factor[k]
  abs_z    <- abs(gamma) / sqrt(BETA)

  # Pick q from a closed-form scale (Lemma 4a upper bound on E[X]) so
  # that the choice of evaluation point does not depend on qmhn().
  q_value <- approx_mean(alpha, BETA, gamma) * q_factor

  val_series <- tryCatch(
    mhn:::.pmhn_force(q_value, alpha, BETA, gamma, method = "series"),
    error = function(e) NA_real_
  )
  val_integrate <- tryCatch(
    mhn:::.pmhn_force(q_value, alpha, BETA, gamma, method = "integrate"),
    error = function(e) NA_real_
  )
  val_dispatch <- tryCatch(
    mhn::pmhn(q_value, alpha = alpha, beta = BETA, gamma = gamma),
    error = function(e) NA_real_
  )
  ref <- tryCatch(
    log_cdf_series_mpfr(alpha, BETA, gamma, q_value),
    error = function(e) {
      structure(NA_real_, error = conditionMessage(e),
                iterations = NA_integer_, reached_max = NA)
    }
  )
  ref_value      <- as.numeric(ref)
  ref_iterations <- attr(ref, "iterations", exact = TRUE)
  if (is.null(ref_iterations)) ref_iterations <- NA_integer_
  ref_reached    <- attr(ref, "reached_max", exact = TRUE)
  if (is.null(ref_reached)) ref_reached <- NA

  series_returned_nan <- is.na(val_series) || is.nan(val_series)

  rel <- function(v) {
    if (is.na(ref_value)) return(NA_real_)
    if (is.na(v))         return(NA_real_)
    denom <- max(abs(ref_value), 1e-300)
    abs(v - ref_value) / denom
  }

  rows[[k]] <- data.frame(
    alpha = alpha, beta = BETA, gamma = gamma,
    q_factor = q_factor, q = q_value, abs_z = abs_z,
    val_series = val_series, val_integrate = val_integrate,
    val_dispatch = val_dispatch, val_ref = ref_value,
    series_returned_nan = series_returned_nan,
    rel_err_series = rel(val_series),
    rel_err_integrate = rel(val_integrate),
    rel_err_dispatch = rel(val_dispatch),
    ref_iterations = ref_iterations,
    ref_reached_max = ref_reached,
    stringsAsFactors = FALSE
  )

  if (k %% progress_every == 0L || k == N) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- if (k < N) elapsed / k * (N - k) else 0
    cat(sprintf("[cdf_series_accuracy] %5d / %5d  elapsed=%6.1fs  eta=%6.1fs\n",
                k, N, elapsed, eta))
  }
}

results <- do.call(rbind, rows)
write.csv(results, RESULT_CSV, row.names = FALSE)

# -----------------------------------------------------------------------
# Headline diagnostics
# -----------------------------------------------------------------------
silent_bug <- !results$series_returned_nan &
              !is.na(results$rel_err_series) &
              results$rel_err_series > 1e-8
breakdown_z <- if (any(results$series_returned_nan, na.rm = TRUE)) {
  min(results$abs_z[results$series_returned_nan], na.rm = TRUE)
} else NA_real_
silent_z <- if (any(silent_bug, na.rm = TRUE)) {
  min(results$abs_z[silent_bug], na.rm = TRUE)
} else NA_real_

cat(sprintf("\n[cdf_series_accuracy] summary:\n"))
cat(sprintf("  grid points                : %d\n", N))
cat(sprintf("  series returned NaN        : %d\n",
            sum(results$series_returned_nan)))
cat(sprintf("  silent_bug_count (>1e-8)   : %d  [smallest |z| = %s]\n",
            sum(silent_bug, na.rm = TRUE),
            if (is.na(silent_z)) "n/a" else sprintf("%.3f", silent_z)))
cat(sprintf("  series breakdown |z| (NaN) : %s\n",
            if (is.na(breakdown_z)) "n/a" else sprintf("%.3f", breakdown_z)))
cat(sprintf("  max rel_err_integrate      : %.3e\n",
            max(results$rel_err_integrate, na.rm = TRUE)))
cat(sprintf("  max rel_err_dispatch       : %.3e\n",
            max(results$rel_err_dispatch, na.rm = TRUE)))

# -----------------------------------------------------------------------
# Diagnostics CSV (single-row environment record)
# -----------------------------------------------------------------------
total_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
diag <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  rmpfr_version = as.character(utils::packageVersion("Rmpfr")),
  mhn_version = as.character(utils::packageVersion("mhn")),
  mode = if (QUICK) "QUICK" else "FULL",
  prec_bits = PREC_BITS,
  max_terms = MAX_TERMS,
  grid_points = N,
  total_secs = round(total_secs, 2),
  series_nan_count = sum(results$series_returned_nan),
  silent_bug_count = sum(silent_bug, na.rm = TRUE),
  series_nan_min_z = breakdown_z,
  silent_bug_min_z = silent_z,
  max_rel_err_series_when_finite = max(
    results$rel_err_series[!results$series_returned_nan], na.rm = TRUE),
  max_rel_err_integrate = max(results$rel_err_integrate, na.rm = TRUE),
  max_rel_err_dispatch = max(results$rel_err_dispatch, na.rm = TRUE),
  ref_max_iterations = max(results$ref_iterations, na.rm = TRUE),
  ref_reached_max_count = sum(results$ref_reached_max == TRUE, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(diag, DIAG_CSV, row.names = FALSE)

cat(sprintf("\n[cdf_series_accuracy] wrote %s\n", RESULT_CSV))
cat(sprintf("[cdf_series_accuracy] wrote %s\n", DIAG_CSV))
cat(sprintf("[cdf_series_accuracy] total time = %.1fs\n", total_secs))
