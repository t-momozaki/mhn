#!/usr/bin/env Rscript
# series_breakdown_audit.R --- Independent verification of three claims
# that the package's pmhn/Psi dispatch implicitly depends on.
#
# Background: pmhn()'s dispatch implicitly depends on three claims
# about the Sun et al. (2023, Lemma 10/11) Psi series and the Lemma 1b
# CDF series.  This script verifies all three at double precision
# against Rmpfr 500-bit ground truth.
#
# Paper refs: Sun et al. (2023), Supplementary Section 1
#
# Claim A (Sun et al. 2023, Supplementary Section 1 qualitative
#          observation; paraphrased for readability -- the original
#          wording is slightly ungrammatical, e.g. "the each terms of
#          the A(k) ... that involving 'lgamma'"):
#   In the case when gamma < 0, the series approximation procedure
#   appears to be inefficient.  From empirical experiments Sun et al.
#   observed that the accumulated errors contributed by the
#   computation of each term of A(k) and B(k) (computation error that
#   involves the 'lgamma' function implemented in base R) appear to be
#   significant compared to the functional value of the Fox-Wright
#   function.
#   ==> for gamma < 0, double-precision Psi-via-Lemma-10-series breaks
#       down at some |z|.
#
# Claim B (package's own extension of the qualitative observation in
#          Claim A above):
#   The same lgamma-rounding issue applies to the Lemma 1b CDF series
#   for gamma < 0.  This is *not* stated in Sun et al.; it is the
#   reason mhn_cdf_series.cpp installs a runtime cancellation guard
#   that returns NaN when alternating-sign cancellation would exceed
#   the user's tolerance.
#
# Claim C (theoretical control):
#   For gamma > 0, every term in both series is non-negative, so
#   sign-separated log-sum-exp accumulation should be numerically
#   stable at double precision regardless of |z|.
#
# All three are verified independently by comparing double-precision
# implementations against the same algorithm at Rmpfr 500-bit precision
# (gold standard).  No C++ code is touched; the audit is entirely in R.
#
# Output:
#   mhn/inst/audits/results/series_breakdown_audit_<YYYYMMDD>.csv
#   mhn/inst/audits/results/series_breakdown_audit_diagnostics_<YYYYMMDD>.csv
#
# Invocation (from repository root):
#   Rscript mhn/inst/audits/series_breakdown_audit.R

suppressPackageStartupMessages({
  if (!requireNamespace("Rmpfr", quietly = TRUE))
    stop("series_breakdown_audit.R requires the 'Rmpfr' package")
})

OUTDIR <- Sys.getenv("MHN_BENCH_OUTDIR", unset = "")
if (!nzchar(OUTDIR)) OUTDIR <- file.path("mhn", "inst", "audits", "results")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("series_breakdown_audit_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR,
                        sprintf("series_breakdown_audit_diagnostics_%s.csv", TODAY))

PREC_BITS <- 500L
MAX_TERMS_MPFR <- 80000L
MAX_TERMS_DOUBLE <- 80000L
PEAK_BUFFER <- 200L
Q_LEMMA10 <- 0.5
TOL <- sqrt(.Machine$double.eps)

# -----------------------------------------------------------------------
# Sun et al. (2023, Lemma 10) truncation helpers (R ports of the C++
# anonymous-namespace functions in mhn/src/mhn_psi_series.cpp).  Inputs
# always use abs(z) so the formulas are well-defined for either sign of
# gamma.  The Errata correction (alpha_adj * z^2, not alpha_adj * z) is
# applied; see the "Errata in Sun et al. (2023)" subsection of the
# package's theory vignette for the derivation.
# -----------------------------------------------------------------------
lemma10_C <- function(alpha_adj, z2, q, c1, c2) {
  denom <- 8 * q
  t <- c1 - z2
  disc <- t^2 - denom * (c2 - alpha_adj * z2)
  if (disc > 0) {
    cval <- (-t + sqrt(disc)) / denom
    return(max(as.integer(ceiling(cval)), 1L))
  }
  return(1L)
}

lemma10_K <- function(a, abs_z, z2, C, q, log_eps_quarter, is_A) {
  if (is_A) {
    log_T_C <- lgamma(a + C) + C * log(z2) - lgamma(2 * C + 1)
  } else {
    log_T_C <- lgamma(a + 0.5 + C) + (2 * C + 1) * log(abs_z) -
                 lgamma(2 * C + 2)
  }
  extra <- 0L
  if (log_T_C > log_eps_quarter) {
    extra <- as.integer(ceiling((log_T_C - log_eps_quarter) / (-log(q))))
    if (extra < 0L) extra <- 0L
  }
  return(C + extra)
}

compute_lemma10_K <- function(alpha, abs_z, tol = TOL) {
  z2 <- abs_z^2
  log_eps_q4 <- log(tol / 4)
  C1 <- lemma10_C(alpha,       z2, Q_LEMMA10, 6 * Q_LEMMA10,  4 * Q_LEMMA10)
  C2 <- lemma10_C(alpha + 1.0, z2, Q_LEMMA10, 10 * Q_LEMMA10, 12 * Q_LEMMA10)
  K1 <- lemma10_K(alpha / 2,         abs_z, z2, C1, Q_LEMMA10, log_eps_q4, TRUE)
  K2 <- lemma10_K((alpha + 1.0) / 2, abs_z, z2, C2, Q_LEMMA10, log_eps_q4, FALSE)
  list(K1 = K1, K2 = K2, K = max(K1, K2))
}

log_sum_exp <- function(xs) {
  xs <- xs[is.finite(xs)]
  if (length(xs) == 0L) return(-Inf)
  m <- max(xs)
  m + log(sum(exp(xs - m)))
}

# =======================================================================
# Experiment 1: Psi via Lemma 10 series
# =======================================================================
#
# Lemma 10: Psi[alpha/2, z] = sum_{k=0..inf} A(k) + sum_{k=0..inf} B(k)
#   A(k) = Gamma(alpha/2 + k) * z^(2k) / (2k)!           -- always >= 0
#   B(k) = Gamma((alpha+1)/2 + k) * z^(2k+1) / (2k+1)!  -- sign(B) = sign(z)
#
# For gamma < 0 we accumulate A in log_pos and |B| in log_neg, then
# combine via log_diff_exp.  Catastrophic cancellation manifests as
# log_S_neg >= log_S_pos and we return NA to mark it.

log_psi_series_double <- function(alpha, beta, gamma, tol = TOL) {
  z <- gamma / sqrt(beta)
  abs_z <- abs(z)
  if (abs_z == 0) return(lgamma(alpha / 2))
  log_abs_z <- log(abs_z)
  kk <- compute_lemma10_K(alpha, abs_z, tol)

  log_pos <- numeric(0L)
  log_neg <- numeric(0L)

  for (k in seq.int(0, kk$K)) {
    log_A <- lgamma(alpha / 2 + k) + 2 * k * log_abs_z - lgamma(2 * k + 1)
    log_pos <- c(log_pos, log_A)
    log_abs_B <- lgamma((alpha + 1) / 2 + k) + (2 * k + 1) * log_abs_z -
                   lgamma(2 * k + 2)
    if (gamma >= 0) log_pos <- c(log_pos, log_abs_B)
    else            log_neg <- c(log_neg, log_abs_B)
  }

  log_S_pos <- log_sum_exp(log_pos)
  if (length(log_neg) == 0L) return(log_S_pos)
  log_S_neg <- log_sum_exp(log_neg)
  if (!is.finite(log_S_pos)) return(NA_real_)
  if (log_S_neg >= log_S_pos) return(NA_real_)
  log_S_pos + log1p(-exp(log_S_neg - log_S_pos))
}

log_psi_series_mpfr <- function(alpha, beta, gamma, prec_bits = PREC_BITS,
                                max_terms = MAX_TERMS_MPFR) {
  ma <- Rmpfr::mpfr(alpha, prec_bits)
  mb <- Rmpfr::mpfr(beta,  prec_bits)
  mg <- Rmpfr::mpfr(gamma, prec_bits)
  z  <- mg / sqrt(mb)
  zero <- Rmpfr::mpfr(0, prec_bits)
  one  <- Rmpfr::mpfr(1, prec_bits)
  half <- one / 2

  # Shift-by-2 recurrence for Gamma((alpha + n) / 2).
  g_even <- gamma(ma / 2)
  g_odd  <- gamma(ma / 2 + half)

  zfact <- one
  S     <- zero
  max_term_abs <- zero
  peak_index <- 0L
  rel_tol <- Rmpfr::mpfr(2, prec_bits)^-(prec_bits - 16L)

  for (i in seq.int(0, max_terms)) {
    is_even <- (i %% 2L == 0L)
    g_n     <- if (is_even) g_even else g_odd
    arg_i   <- (ma + i) / 2

    term <- g_n * zfact
    S    <- S + term

    term_abs <- abs(term)
    if (term_abs > max_term_abs) {
      max_term_abs <- term_abs
      peak_index   <- i
    }

    if (is_even) g_even <- g_even * arg_i else g_odd <- g_odd * arg_i
    zfact <- zfact * z / (i + 1L)

    if (i >= peak_index + PEAK_BUFFER &&
        term_abs < max_term_abs * rel_tol) break
  }

  if (S <= zero) return(NA_real_)  # mpfr couldn't resolve sign positively
  as.double(log(S))
}

# =======================================================================
# Experiment 2: CDF via Lemma 1b series (with Lemma 10 K truncation)
# =======================================================================
#
# F(x | alpha, beta, gamma) = (1 / Psi) * sum_{i=0..} T_i
#   T_i = Gamma((alpha + i)/2) * z^i / i! * P((alpha + i)/2, beta * x^2)
#
# Since |P_i| <= 1, |T_{2k}| <= A(k) and |T_{2k+1}| <= |B(k)|, so the
# Lemma 10 K for the Psi series also bounds the truncation error of the
# CDF series.

log_cdf_series_double <- function(alpha, beta, gamma, x, log_psi,
                                  tol = TOL) {
  if (!(x > 0)) return(-Inf)
  z <- gamma / sqrt(beta)
  abs_z <- abs(z)
  if (abs_z == 0) {
    log_p0 <- pgamma(beta * x^2, alpha / 2, lower.tail = TRUE, log.p = TRUE)
    return(lgamma(alpha / 2) + log_p0 - log_psi)
  }
  log_abs_z <- log(abs_z)
  y <- beta * x^2
  kk <- compute_lemma10_K(alpha, abs_z, tol)
  i_max <- 2L * kk$K + 1L

  log_pos <- numeric(0L)
  log_neg <- numeric(0L)

  for (i in seq.int(0, i_max)) {
    s_i <- (alpha + i) / 2
    log_P <- pgamma(y, s_i, lower.tail = TRUE, log.p = TRUE)
    if (!is.finite(log_P) || log_P == -Inf) next
    log_T <- i * log_abs_z + lgamma(s_i) + log_P - lgamma(i + 1)
    if (!is.finite(log_T)) next
    if (gamma < 0 && (i %% 2L == 1L)) log_neg <- c(log_neg, log_T)
    else                              log_pos <- c(log_pos, log_T)
  }

  if (length(log_pos) == 0L) return(-Inf)
  log_S_pos <- log_sum_exp(log_pos)
  if (length(log_neg) == 0L) return(log_S_pos - log_psi)
  log_S_neg <- log_sum_exp(log_neg)
  if (log_S_neg >= log_S_pos) return(NA_real_)
  log_S_pos + log1p(-exp(log_S_neg - log_S_pos)) - log_psi
}

log_cdf_series_mpfr <- function(alpha, beta, gamma, x, log_psi_mpfr,
                                prec_bits = PREC_BITS,
                                max_terms = MAX_TERMS_MPFR) {
  ma <- Rmpfr::mpfr(alpha, prec_bits)
  mb <- Rmpfr::mpfr(beta,  prec_bits)
  mg <- Rmpfr::mpfr(gamma, prec_bits)
  mx <- Rmpfr::mpfr(x,     prec_bits)
  z  <- mg / sqrt(mb)
  y_mpfr <- mb * mx^2
  y_double <- as.double(y_mpfr)

  zero <- Rmpfr::mpfr(0, prec_bits)
  one  <- Rmpfr::mpfr(1, prec_bits)
  half <- one / 2

  g_even <- gamma(ma / 2)
  g_odd  <- gamma(ma / 2 + half)
  zfact  <- one

  S_num <- zero
  max_term_abs <- zero
  peak_index <- 0L
  rel_tol <- Rmpfr::mpfr(2, prec_bits)^-(prec_bits - 16L)

  for (i in seq.int(0, max_terms)) {
    is_even <- (i %% 2L == 0L)
    g_n     <- if (is_even) g_even else g_odd
    arg_i   <- (ma + i) / 2
    arg_i_d <- (alpha + i) / 2
    P_i     <- Rmpfr::mpfr(pgamma(y_double, shape = arg_i_d), prec_bits)

    full_term <- g_n * zfact * P_i
    S_num     <- S_num + full_term

    term_abs <- abs(full_term)
    if (term_abs > max_term_abs) {
      max_term_abs <- term_abs
      peak_index   <- i
    }

    if (is_even) g_even <- g_even * arg_i else g_odd <- g_odd * arg_i
    zfact <- zfact * z / (i + 1L)

    if (i >= peak_index + PEAK_BUFFER &&
        term_abs < max_term_abs * rel_tol) break
  }

  if (S_num <= zero) return(NA_real_)
  as.double(log(S_num) - log_psi_mpfr)
}

# =======================================================================
# Grid + main loop
# =======================================================================
ALPHAS <- c(0.5, 1.5, 3, 10)
GAMMAS <- c(-150, -100, -75, -50, -30, -25, -22, -20, -18, -15,
            -10, -8, -5, -2, -1,
             1,  2,  5,  8, 10, 15, 18, 20, 22, 25, 30, 50, 75, 100, 150)
Q_FACTORS <- c(1 / 3, 1, 3)
BETA <- 1.0

# Sun et al. (2023) Lemma 4a upper bound on E[X], proven for alpha > 1.
# Used here purely as a scale to pick informative q points; for
# alpha <= 1 we extrapolate the same formula since the audit only
# needs a reasonable positive scale, not a rigorous bound.
approx_mean <- function(alpha, beta, gamma) {
  (gamma + sqrt(gamma^2 + 8 * alpha * beta)) / (4 * beta)
}

cat(sprintf("[series_breakdown_audit] prec_bits=%d alphas=%d gammas=%d q_factors=%d\n",
            PREC_BITS, length(ALPHAS), length(GAMMAS), length(Q_FACTORS)))
cat(sprintf("[series_breakdown_audit] result_csv=%s\n", RESULT_CSV))

rel <- function(v_d, v_m) {
  if (is.na(v_d) || is.na(v_m)) return(NA_real_)
  abs(v_d - v_m) / max(abs(v_m), 1e-300)
}

rows <- list()
t0 <- Sys.time()

# ---- Experiment 1 + 3 (Psi): one row per (alpha, gamma) ----
N_psi <- length(ALPHAS) * length(GAMMAS)
k_psi <- 0L
for (alpha in ALPHAS) for (gamma in GAMMAS) {
  k_psi <- k_psi + 1L
  abs_z <- abs(gamma) / sqrt(BETA)

  log_psi_d <- tryCatch(
    log_psi_series_double(alpha, BETA, gamma),
    error = function(e) NA_real_)
  log_psi_m <- tryCatch(
    log_psi_series_mpfr(alpha, BETA, gamma),
    error = function(e) NA_real_)

  rows[[length(rows) + 1L]] <- data.frame(
    experiment = "psi", alpha = alpha, gamma = gamma, sign_gamma = sign(gamma),
    abs_z = abs_z, q_factor = NA_real_, q = NA_real_,
    log_val_double = log_psi_d, log_val_mpfr = log_psi_m,
    double_returned_nan = is.na(log_psi_d),
    rel_err = rel(log_psi_d, log_psi_m),
    stringsAsFactors = FALSE
  )

  if (k_psi %% max(1L, N_psi %/% 10L) == 0L) {
    cat(sprintf("[psi] %3d/%3d  elapsed=%6.1fs\n", k_psi, N_psi,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

# ---- Experiment 2 (CDF): one row per (alpha, gamma, q_factor) ----
N_cdf <- length(ALPHAS) * length(GAMMAS) * length(Q_FACTORS)
k_cdf <- 0L
for (alpha in ALPHAS) for (gamma in GAMMAS) {
  log_psi_d <- tryCatch(log_psi_series_double(alpha, BETA, gamma),
                        error = function(e) NA_real_)
  log_psi_m <- tryCatch(log_psi_series_mpfr(alpha, BETA, gamma),
                        error = function(e) NA_real_)

  for (qf in Q_FACTORS) {
    k_cdf <- k_cdf + 1L
    q <- approx_mean(alpha, BETA, gamma) * qf
    abs_z <- abs(gamma) / sqrt(BETA)

    log_cdf_d <- tryCatch(
      log_cdf_series_double(alpha, BETA, gamma, q,
                            log_psi = if (is.na(log_psi_d)) log_psi_m else log_psi_d),
      error = function(e) NA_real_)
    log_cdf_m <- tryCatch(
      log_cdf_series_mpfr(alpha, BETA, gamma, q,
                          log_psi_mpfr = Rmpfr::mpfr(log_psi_m, PREC_BITS)),
      error = function(e) NA_real_)

    rows[[length(rows) + 1L]] <- data.frame(
      experiment = "cdf", alpha = alpha, gamma = gamma, sign_gamma = sign(gamma),
      abs_z = abs_z, q_factor = qf, q = q,
      log_val_double = log_cdf_d, log_val_mpfr = log_cdf_m,
      double_returned_nan = is.na(log_cdf_d),
      rel_err = rel(log_cdf_d, log_cdf_m),
      stringsAsFactors = FALSE
    )
  }

  if (k_cdf %% max(1L, N_cdf %/% 10L) < length(Q_FACTORS)) {
    cat(sprintf("[cdf] %3d/%3d  elapsed=%6.1fs\n", k_cdf, N_cdf,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

results <- do.call(rbind, rows)
write.csv(results, RESULT_CSV, row.names = FALSE)

# -----------------------------------------------------------------------
# Per-experiment summary
# -----------------------------------------------------------------------
report_breakdown <- function(df, label) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  rows total                  : %d\n", nrow(df)))
  cat(sprintf("  double returned NaN         : %d\n", sum(df$double_returned_nan)))

  # Empirical breakdown |z|: smallest |z| where rel_err > 1e-8 (and not NaN)
  bad <- !df$double_returned_nan & !is.na(df$rel_err) & df$rel_err > 1e-8
  if (any(bad)) {
    cat(sprintf("  silent breakdown count      : %d  (rel_err > 1e-8 and no NaN)\n",
                sum(bad)))
    cat(sprintf("  smallest |z| with breakdown : %.3f\n",
                min(df$abs_z[bad])))
  } else {
    cat("  silent breakdown count      : 0\n")
  }

  nan_rows <- df$double_returned_nan
  if (any(nan_rows)) {
    cat(sprintf("  smallest |z| where NaN fired: %.3f\n",
                min(df$abs_z[nan_rows])))
  }

  cat(sprintf("  max rel_err (finite double) : %.3e\n",
              max(df$rel_err[!df$double_returned_nan], na.rm = TRUE)))
}

cat("\n[series_breakdown_audit] summary:")

# Experiment 1: Psi for gamma < 0 (Claim A)
report_breakdown(
  results[results$experiment == "psi" & results$gamma < 0, ],
  "Claim A: Psi via Lemma 10 series, gamma < 0")

# Experiment 3a: Psi for gamma > 0 (Claim C control)
report_breakdown(
  results[results$experiment == "psi" & results$gamma > 0, ],
  "Claim C control: Psi via Lemma 10 series, gamma > 0")

# Experiment 2: CDF for gamma < 0 (Claim B)
report_breakdown(
  results[results$experiment == "cdf" & results$gamma < 0, ],
  "Claim B: CDF via Lemma 1b series, gamma < 0")

# Experiment 3b: CDF for gamma > 0 (Claim C control)
report_breakdown(
  results[results$experiment == "cdf" & results$gamma > 0, ],
  "Claim C control: CDF via Lemma 1b series, gamma > 0")

total_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

diag <- data.frame(
  timestamp     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version     = paste(R.version$major, R.version$minor, sep = "."),
  rmpfr_version = as.character(utils::packageVersion("Rmpfr")),
  prec_bits     = PREC_BITS,
  max_terms_mpfr = MAX_TERMS_MPFR,
  peak_buffer    = PEAK_BUFFER,
  rows           = nrow(results),
  total_secs     = round(total_secs, 2),
  stringsAsFactors = FALSE
)
write.csv(diag, DIAG_CSV, row.names = FALSE)

cat(sprintf("\n[series_breakdown_audit] wrote %s\n", RESULT_CSV))
cat(sprintf("[series_breakdown_audit] wrote %s\n", DIAG_CSV))
cat(sprintf("[series_breakdown_audit] total time = %.1fs\n", total_secs))
