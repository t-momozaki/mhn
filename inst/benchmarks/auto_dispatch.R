#!/usr/bin/env Rscript
# auto_dispatch.R --- Benchmark for finalising the rmhn() auto path.
#
# Finalises the rmhn(method = "auto") dispatch by benchmarking RTDR
# vs Sun across a parameter grid and a range of (n_per_call, n_calls)
# patterns.
#
# Grid:    alpha in {0.3, 0.5, 0.7, 0.8, 1, 1.5, 3, 5, 10, 100}  (now spans
#          the non-log-concave alpha < 1 region, the headline corner)
#          gamma in {-100, -10, -2, -0.5, 0.5, 2, 10, 100}      => 80 points
# Patterns (n_per_call m spans setup-dominated to per-proposal-dominated):
#   A : n_per_call=1     x n_calls=10000  (Gibbs-style, setup dominates)
#   B..B100 : m = 5,10,25,50,100          (the crossover sweep)
#   C : n_per_call=10000 x n_calls=1      (batch, per-proposal dominates)
#
# alpha < 1 needs NO new sampler: Sun's Algorithm 3 is valid for all alpha > 0
# when gamma <= 0, so alpha<1, gamma<=0 is a fair rtdr-vs-A3 comparison; and
# alpha<1, gamma>0 is the corner where Sun would need the deliberately-omitted
# Algorithm 2, so sun is recorded as sun_incompatible and only rtdr's
# acceptance (the 1/e floor that anchors the default) is measured there.
#
# This is a FAIR, within-C++ comparison: rtdr and sun are both this package's
# own C++ kernels timed on one machine, so the confound that makes an absolute
# comparison against Gao & Wang's published R timings meaningless does not
# apply here.  All timings and acceptance rates below are comparable.
#
# What it now measures (beyond raw per-pattern times):
#   - Acceptance rate for rtdr and sun (machine-independent efficiency).
#   - A setup / per-proposal cost decomposition per kernel: with the thunk
#     running k calls of m draws each, median_us = k*(T_setup + m*T_prop), so
#     the per-call time median_us/k = T_setup + m*T_prop.  T_setup is read at
#     m=1 (pattern A) and T_prop at m=10^4 (pattern C).
#   - The EMPIRICAL Gibbs<->batch crossover: the smallest m at which the
#     rtdr-vs-sun winner flips, read directly from the sweep (no assumed 25).
#   - The SHIPPED rmhn(method="auto") pick (mirroring src/mhn_rmhn.cpp) at the
#     Gibbs (m=1) and batch (m=10^4) shapes, and whether it matches the
#     measured winner -- flagging cells where the dispatch rule is suboptimal.
#
# Outputs:
#   auto_dispatch_<YYYYMMDD>.csv                (per pattern x method raw times)
#   auto_dispatch_diagnostics_<YYYYMMDD>.csv    (per point acceptance + kernel diag)
#   auto_dispatch_recommendation_<YYYYMMDD>.csv (per point: decomposition,
#       crossover, measured winner, shipped auto pick, and mismatch flags)
#
# Invocation (from repository root):
#   Rscript mhn/inst/benchmarks/auto_dispatch.R
#
# Optional environment variables:
#   MHN_BENCH_QUICK=1 -> reduce grid and iterations for smoke testing
#   MHN_BENCH_ITER=N  -> override bench::mark iterations (default 5)
#   MHN_BENCH_OUTDIR  -> override output directory

suppressPackageStartupMessages({
  if (!requireNamespace("bench", quietly = TRUE))
    stop("auto_dispatch.R requires the 'bench' package. Install via install.packages('bench').")
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("auto_dispatch.R requires the 'mhn' package. Run R CMD INSTALL mhn first.")
})

QUICK <- isTRUE(nchar(Sys.getenv("MHN_BENCH_QUICK")) > 0L)
# Iterations are set to 50 because iter=5 was noisy enough to randomly
# invert sun/rtdr decisions, and iter=20 still showed run-to-run flips
# on 5-12 of 48 borderline cells (especially alpha >= 10, gamma in
# {2, 10}).  iter=50 keeps relative IQR below 3% on this hardware so
# that the 10% decision threshold remains meaningful.
BENCH_ITER <- {
  v <- suppressWarnings(as.integer(Sys.getenv("MHN_BENCH_ITER", unset = "")))
  if (is.na(v) || v < 1L) (if (QUICK) 2L else 50L) else v
}

ALPHAS <- if (QUICK) c(0.5, 1, 3) else c(0.3, 0.5, 0.7, 0.8, 1, 1.5, 3, 5, 10, 100)
GAMMAS <- if (QUICK) c(-10, -0.5, 2) else c(-100, -10, -2, -0.5, 0.5, 2, 10, 100)

PATTERNS <- list(
  A    = list(n_per_call = 1L,     n_calls = if (QUICK) 1000L else 10000L),
  B    = list(n_per_call = 5L,     n_calls = if (QUICK) 200L  else 2000L),
  B10  = list(n_per_call = 10L,    n_calls = if (QUICK) 100L  else 1000L),
  B25  = list(n_per_call = 25L,    n_calls = if (QUICK) 40L   else 400L),
  B50  = list(n_per_call = 50L,    n_calls = if (QUICK) 20L   else 200L),
  B100 = list(n_per_call = 100L,   n_calls = if (QUICK) 10L   else 100L),
  C    = list(n_per_call = if (QUICK) 1000L else 10000L, n_calls = 1L)
)

DECISION_THRESHOLD <- 0.10  # require a 10% wall-clock advantage to flip from RTDR to Sun

OUTDIR <- Sys.getenv("MHN_BENCH_OUTDIR", unset = "")
if (!nzchar(OUTDIR)) OUTDIR <- file.path("mhn", "inst", "benchmarks", "results")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("auto_dispatch_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR, sprintf("auto_dispatch_diagnostics_%s.csv", TODAY))
# Single-row environment / provenance record.  A separate file because
# DIAG_CSV above already holds the per-point acceptance diagnostics.
PROV_CSV   <- file.path(OUTDIR, sprintf("auto_dispatch_provenance_%s.csv", TODAY))
# Per-point recommendation: cost decomposition, empirical crossover, measured
# winner per call shape, and the shipped auto pick with mismatch flags.
REC_CSV    <- file.path(OUTDIR, sprintf("auto_dispatch_recommendation_%s.csv", TODAY))

cat(sprintf("[auto_dispatch] mode=%s iterations=%d alpha=%d gamma=%d patterns=%d points=%d\n",
            if (QUICK) "QUICK" else "FULL", BENCH_ITER,
            length(ALPHAS), length(GAMMAS), length(PATTERNS),
            length(ALPHAS) * length(GAMMAS)))
cat(sprintf("[auto_dispatch] result_csv=%s\n", RESULT_CSV))
cat(sprintf("[auto_dispatch] diag_csv=%s\n", DIAG_CSV))

is_special_case <- function(alpha, gamma) {
  EPS <- sqrt(.Machine$double.eps)
  if (abs(alpha - 1.0) < EPS && abs(gamma) < EPS) return("half_normal")
  if (abs(gamma) < EPS) return("sqrt_gamma")
  if (abs(alpha - 1.0) < EPS) return("trunc_normal")
  NA_character_
}

sun_applicable <- function(alpha, gamma) !(alpha < 1.0 && gamma > 0.0)

sun_algo_for <- function(alpha, gamma) {
  if (gamma > 0.0 && alpha > 1.0) "algo1"
  else if (gamma <= 0.0)          "algo3"
  else                             NA_character_
}

make_thunk <- function(alpha, beta, gamma, method, n_per_call, n_calls) {
  force(alpha); force(beta); force(gamma); force(method)
  force(n_per_call); force(n_calls)
  function() {
    for (i in seq_len(n_calls)) {
      mhn::rmhn(n_per_call, alpha, beta, gamma, method = method)
    }
    invisible(NULL)
  }
}

bench_one <- function(alpha, gamma, method, pattern_name, pattern, seed = 1L) {
  set.seed(seed)
  thunk <- make_thunk(alpha, 1.0, gamma, method, pattern$n_per_call, pattern$n_calls)
  res <- tryCatch(
    # NB: do not pass time_unit = "us" here. That option makes
    # as.numeric(res$median) return microseconds (and collapses
    # as.numeric(res$time) to 0); combined with the `* 1e6` below it
    # inflated median_us/iqr_us by 1e6. Left as seconds, so the
    # `* 1e6` conversion yields correct microseconds.
    bench::mark(thunk(), iterations = BENCH_ITER, check = FALSE,
                filter_gc = FALSE),
    error = function(e) structure(list(error_msg = conditionMessage(e)),
                                  class = "bench_failure")
  )
  if (inherits(res, "bench_failure")) {
    return(data.frame(
      alpha = alpha, gamma = gamma, pattern = pattern_name, method = method,
      median_us = NA_real_, iqr_us = NA_real_, n_gc = NA_integer_,
      error = res$error_msg, stringsAsFactors = FALSE
    ))
  }
  median_us <- as.numeric(res$median[[1]]) * 1e6
  times_us <- as.numeric(res$time[[1]]) * 1e6
  iqr_us <- if (length(times_us) >= 2L) IQR(times_us) else NA_real_
  n_gc <- as.integer(res$n_gc[[1]])
  data.frame(
    alpha = alpha, gamma = gamma, pattern = pattern_name, method = method,
    median_us = median_us, iqr_us = iqr_us, n_gc = n_gc,
    error = NA_character_, stringsAsFactors = FALSE
  )
}

diagnostics_one <- function(alpha, gamma, n_diag = 10000L, seed = 1L) {
  beta <- 1.0
  spec <- is_special_case(alpha, gamma)
  out <- list(
    alpha = alpha, gamma = gamma,
    rtdr_acc = NA_real_, sun_acc = NA_real_,
    sun_algo = NA_character_, sun_algo1_proposal = NA_character_,
    sun_algo3_used_inflex = NA, sun_algo3_m_init = NA_real_, sun_algo3_m = NA_real_,
    newton_failure = FALSE, keff_fallback = FALSE,
    rtdr_warning = NA_character_, sun_warning = NA_character_,
    rtdr_error = NA_character_, sun_error = NA_character_,
    notes = ifelse(is.na(spec), "", spec)
  )

  set.seed(seed)
  rtdr_res <- withCallingHandlers(
    tryCatch(
      mhn:::.rmhn_rtdr_cpp(n_diag, alpha, beta, gamma),
      error = function(e) { out$rtdr_error <<- conditionMessage(e); NULL }
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      # Keep ONLY the first warning to avoid ballooning the CSV when a
      # particular setup repeatedly emits the same warning.
      if (is.na(out$rtdr_warning)) out$rtdr_warning <<- substr(msg, 1L, 200L)
      if (grepl("K_eff", msg, fixed = TRUE) || grepl("region.*fallback", msg, ignore.case = TRUE)) {
        out$keff_fallback <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  if (!is.null(rtdr_res)) {
    retries <- attr(rtdr_res, "rtdr_retries")
    if (!is.null(retries)) out$rtdr_acc <- n_diag / (n_diag + as.numeric(retries))
  } else if (!is.na(out$rtdr_error) && grepl("Newton", out$rtdr_error, ignore.case = TRUE)) {
    out$newton_failure <- TRUE
  }

  algo <- sun_algo_for(alpha, gamma)
  out$sun_algo <- if (is.na(algo)) NA_character_ else algo
  if (!is.na(algo) && sun_applicable(alpha, gamma)) {
    if (algo == "algo1") {
      set.seed(seed)
      sun_res <- withCallingHandlers(
        tryCatch(
          mhn:::.rmhn_sun_algo1_cpp(n_diag, alpha, beta, gamma),
          error = function(e) { out$sun_error <<- conditionMessage(e); NULL }
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (is.na(out$sun_warning)) out$sun_warning <<- substr(msg, 1L, 200L)
          invokeRestart("muffleWarning")
        }
      )
      if (!is.null(sun_res)) {
        retries <- attr(sun_res, "sun_retries")
        if (!is.null(retries)) out$sun_acc <- n_diag / (n_diag + as.numeric(retries))
      }
      dump <- tryCatch(mhn:::.dump_sun_algo1_cpp(alpha, beta, gamma), error = function(e) NULL)
      if (!is.null(dump) && !is.null(dump$proposal)) {
        out$sun_algo1_proposal <- as.character(dump$proposal)
      }
    } else if (algo == "algo3") {
      set.seed(seed)
      sun_res <- withCallingHandlers(
        tryCatch(
          mhn:::.rmhn_sun_algo3_cpp(n_diag, alpha, beta, gamma),
          error = function(e) { out$sun_error <<- conditionMessage(e); NULL }
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (is.na(out$sun_warning)) out$sun_warning <<- substr(msg, 1L, 200L)
          invokeRestart("muffleWarning")
        }
      )
      if (!is.null(sun_res)) {
        retries <- attr(sun_res, "sun_retries")
        if (!is.null(retries)) out$sun_acc <- n_diag / (n_diag + as.numeric(retries))
      }
      dump <- tryCatch(mhn:::.dump_sun_algo3_cpp(alpha, beta, gamma), error = function(e) NULL)
      if (!is.null(dump)) {
        if (!is.null(dump$used_inflex_heuristic))
          out$sun_algo3_used_inflex <- as.logical(dump$used_inflex_heuristic)
        if (!is.null(dump$m_init)) out$sun_algo3_m_init <- as.numeric(dump$m_init)
        if (!is.null(dump$m))      out$sun_algo3_m      <- as.numeric(dump$m)
      }
    }
  }
  out
}

decide_method <- function(rows, alpha, gamma) {
  spec <- is_special_case(alpha, gamma)
  if (!is.na(spec)) return(list(decision = "special", reason = spec))
  if (!sun_applicable(alpha, gamma)) {
    return(list(decision = "rtdr", reason = "sun_incompatible"))
  }
  pick <- function(patt) {
    r <- rows$median_us[rows$pattern == patt & rows$method == "rtdr"]
    s <- rows$median_us[rows$pattern == patt & rows$method == "sun"]
    if (length(r) != 1L || length(s) != 1L || !is.finite(r) || !is.finite(s)) {
      return(NA_character_)
    }
    diff_rel <- (r - s) / max(r, s)
    if (diff_rel >= DECISION_THRESHOLD) "sun"
    else if (-diff_rel >= DECISION_THRESHOLD) "rtdr"
    else "tie"
  }
  pa <- pick("A"); pb <- pick("B")
  if (is.na(pa) || is.na(pb)) {
    return(list(decision = "rtdr", reason = "missing_data"))
  }
  if (pa == "sun" && pb == "sun") return(list(decision = "sun", reason = "AB_consistent_sun"))
  if (pa == "rtdr" && pb == "rtdr") return(list(decision = "rtdr", reason = "AB_consistent_rtdr"))
  list(decision = "rtdr", reason = sprintf("inconsistent(A=%s,B=%s)_default_rtdr", pa, pb))
}

# -----------------------------------------------------------------------
# Cost decomposition, empirical crossover, and shipped-auto comparison
# -----------------------------------------------------------------------
# Per-call time for a pattern = median_us / n_calls = T_setup + m*T_prop.
percall_us <- function(rows, patt, method) {
  v <- rows$median_us[rows$pattern == patt & rows$method == method]
  k <- PATTERNS[[patt]]$n_calls
  if (length(v) == 1L && is.finite(v)) v / k else NA_real_
}
# T_setup from m=1 (pattern A), T_prop from m=10^4 (pattern C).
decomp <- function(rows, method) {
  yA <- percall_us(rows, "A", method); yC <- percall_us(rows, "C", method)
  mA <- PATTERNS$A$n_per_call; mC <- PATTERNS$C$n_per_call
  prop <- if (is.finite(yA) && is.finite(yC) && mC > mA) (yC - yA) / (mC - mA) else NA_real_
  setup <- if (is.finite(yA) && is.finite(prop)) yA - prop * mA else yA
  list(setup_us = setup, prop_us = prop)
}
# Patterns ordered by n_per_call, for scanning the crossover.
PATT_ORDER <- names(PATTERNS)[order(vapply(PATTERNS, function(p) p$n_per_call, numeric(1)))]
winner_at <- function(rows, patt) {
  r <- percall_us(rows, patt, "rtdr"); s <- percall_us(rows, patt, "sun")
  if (!is.finite(r) || !is.finite(s)) return(NA_character_)
  if (r < s) "rtdr" else "sun"
}
# Smallest m at which the rtdr-vs-sun winner flips from its m=1 value; NA if
# one kernel wins throughout (no crossover) or sun is unavailable.
crossover_m <- function(rows) {
  w  <- vapply(PATT_ORDER, function(p) winner_at(rows, p), character(1))
  ms <- vapply(PATT_ORDER, function(p) as.numeric(PATTERNS[[p]]$n_per_call), numeric(1))
  ok <- !is.na(w); w <- w[ok]; ms <- ms[ok]
  if (length(w) < 2L) return(NA_real_)
  flip <- which(w != w[1])
  if (!length(flip)) return(NA_real_)
  ms[flip[1]]
}
# The shipped rmhn(method="auto") pick, mirroring src/mhn_rmhn.cpp.  For a
# scalar-parameter rmhn(m, ...) call, samples_per_setup == m.
shipped_auto <- function(alpha, gamma, m) {
  if (!is.na(is_special_case(alpha, gamma))) return("special")
  if (gamma > 0.0) return(if (alpha > 1.0) "sun" else "rtdr")
  # gamma <= 0 (non-special): large batch -> RTDR, except alpha>=10 where
  # Sun A3's per-proposal cost falls below RTDR's and wins in the batch
  # regime too.  The "large batch" cutoff is 25, raised to 100 for
  # alpha<0.1 where the crossover moves right (see src/mhn_rmhn.cpp).
  n_switch <- if (alpha < 0.1) 100L else 25L
  if (m >= n_switch && alpha < 10.0) "rtdr" else "sun"
}
agree <- function(pick, measured) {
  if (pick == "special" || is.na(measured)) return(NA)
  pick == measured
}

recommend_one <- function(rows, diag_row, alpha, gamma) {
  spec <- is_special_case(alpha, gamma)
  region <- if (!is.na(spec)) "special"
            else if (alpha < 1 && gamma > 0) "corner_alpha<1_gamma>0"
            else "general"
  rda <- decomp(rows, "rtdr"); sda <- decomp(rows, "sun")
  wg <- winner_at(rows, "A"); wb <- winner_at(rows, "C")
  ag <- shipped_auto(alpha, gamma, 1L)
  ab <- shipped_auto(alpha, gamma, as.integer(PATTERNS$C$n_per_call))
  data.frame(
    alpha = alpha, gamma = gamma, region = region, sun_algo = diag_row$sun_algo,
    rtdr_acc = diag_row$rtdr_acc, sun_acc = diag_row$sun_acc,
    rtdr_setup_us = rda$setup_us, rtdr_prop_us = rda$prop_us,
    sun_setup_us  = sda$setup_us, sun_prop_us  = sda$prop_us,
    winner_gibbs = wg, winner_batch = wb, crossover_m = crossover_m(rows),
    auto_pick_gibbs = ag, auto_pick_batch = ab,
    agree_gibbs = agree(ag, wg), agree_batch = agree(ab, wb),
    stringsAsFactors = FALSE)
}

bench_rows <- list(); diag_rows <- list(); rec_rows <- list()
t_start <- Sys.time()
total_points <- length(ALPHAS) * length(GAMMAS); point_idx <- 0L

for (alpha in ALPHAS) {
  for (gamma in GAMMAS) {
    point_idx <- point_idx + 1L
    cat(sprintf("[%3d/%3d] alpha=%g gamma=%g  ", point_idx, total_points, alpha, gamma))
    flush.console()

    diag_row <- diagnostics_one(alpha, gamma)

    rows <- list()
    for (pname in names(PATTERNS)) {
      patt <- PATTERNS[[pname]]
      rows[[length(rows) + 1L]] <- bench_one(alpha, gamma, "rtdr", pname, patt)
      if (sun_applicable(alpha, gamma)) {
        rows[[length(rows) + 1L]] <- bench_one(alpha, gamma, "sun", pname, patt)
      } else {
        rows[[length(rows) + 1L]] <- data.frame(
          alpha = alpha, gamma = gamma, pattern = pname, method = "sun",
          median_us = NA_real_, iqr_us = NA_real_, n_gc = NA_integer_,
          error = "sun_incompatible", stringsAsFactors = FALSE
        )
      }
    }
    rows_df <- do.call(rbind, rows)

    dec <- decide_method(rows_df, alpha, gamma)
    rows_df$decision <- dec$decision
    rows_df$decision_reason <- dec$reason
    rows_df$acceptance <- ifelse(rows_df$method == "rtdr", diag_row$rtdr_acc,
                          ifelse(rows_df$method == "sun",  diag_row$sun_acc, NA_real_))
    rows_df$notes <- diag_row$notes

    bench_rows[[length(bench_rows) + 1L]] <- rows_df
    diag_rows[[length(diag_rows) + 1L]]   <- as.data.frame(diag_row, stringsAsFactors = FALSE)
    rec_rows[[length(rec_rows) + 1L]]     <- recommend_one(rows_df, diag_row, alpha, gamma)

    cat(sprintf("decision=%s (%s) rtdr_acc=%s sun_acc=%s newton_fail=%s keff_fb=%s\n",
                dec$decision, dec$reason,
                formatC(diag_row$rtdr_acc, digits = 3, format = "f"),
                formatC(diag_row$sun_acc,  digits = 3, format = "f"),
                diag_row$newton_failure, diag_row$keff_fallback))
  }
}

bench_df <- do.call(rbind, bench_rows)
diag_df  <- do.call(rbind, diag_rows)
rec_df   <- do.call(rbind, rec_rows)
write.csv(bench_df, RESULT_CSV, row.names = FALSE)
write.csv(diag_df,  DIAG_CSV,   row.names = FALSE)
write.csv(rec_df,   REC_CSV,    row.names = FALSE)

t_end <- Sys.time()
elapsed_min <- as.numeric(difftime(t_end, t_start, units = "mins"))

cat("\n========== SUMMARY ==========\n")
cat(sprintf("Elapsed: %.1f min   (start=%s end=%s)\n",
            elapsed_min, format(t_start, "%H:%M:%S"), format(t_end, "%H:%M:%S")))

cat("\n[Decision matrix per (alpha, gamma)]\n")
dec_table <- unique(bench_df[, c("alpha", "gamma", "decision", "decision_reason")])
dec_wide  <- reshape(dec_table[, c("alpha", "gamma", "decision")],
                     idvar = "alpha", timevar = "gamma", direction = "wide")
print(dec_wide, row.names = FALSE)

cat("\n[Shipped auto vs measured winner]  (general cells only)\n")
gen <- rec_df[rec_df$region == "general", ]
mm_g <- gen[!is.na(gen$agree_gibbs) & !gen$agree_gibbs, ]
mm_b <- gen[!is.na(gen$agree_batch) & !gen$agree_batch, ]
cat(sprintf("  Gibbs (m=1)  : shipped auto matches measured winner in %d / %d cells\n",
            sum(gen$agree_gibbs, na.rm = TRUE), sum(!is.na(gen$agree_gibbs))))
cat(sprintf("  Batch (m=1e4): shipped auto matches measured winner in %d / %d cells\n",
            sum(gen$agree_batch, na.rm = TRUE), sum(!is.na(gen$agree_batch))))
if (nrow(mm_g) || nrow(mm_b)) {
  cat("  MISMATCH cells (dispatcher may be suboptimal):\n")
  for (i in seq_len(nrow(mm_g)))
    cat(sprintf("    a=%-5g g=%-6g Gibbs: auto=%s but measured winner=%s\n",
                mm_g$alpha[i], mm_g$gamma[i], mm_g$auto_pick_gibbs[i], mm_g$winner_gibbs[i]))
  for (i in seq_len(nrow(mm_b)))
    cat(sprintf("    a=%-5g g=%-6g batch: auto=%s but measured winner=%s\n",
                mm_b$alpha[i], mm_b$gamma[i], mm_b$auto_pick_batch[i], mm_b$winner_batch[i]))
} else {
  cat("  No mismatches: the shipped dispatch rule agrees with the measured optimum.\n")
}
cat("\n[Empirical crossover m (rtdr<->sun winner flips)]  (general, gamma<=0)\n")
xo <- gen[gen$gamma <= 0 & !is.na(gen$crossover_m), c("alpha","gamma","winner_gibbs","crossover_m","winner_batch")]
if (nrow(xo)) print(xo, row.names = FALSE) else cat("  (no crossover cells; one kernel wins throughout)\n")

cat("\n[Carry-over verification]\n")
cat(sprintf("  Newton failures              : %d / %d points\n",
            sum(diag_df$newton_failure, na.rm = TRUE), nrow(diag_df)))
cat(sprintf("  K_eff fallbacks fired        : %d / %d points\n",
            sum(diag_df$keff_fallback, na.rm = TRUE), nrow(diag_df)))
algo1_idx <- diag_df$sun_algo == "algo1" & diag_df$alpha %in% c(1.5, 3) & !is.na(diag_df$sun_acc)
cat(sprintf("  Algo1 alpha in {1.5,3} acc>=0.5 : %d / %d cases\n",
            sum(diag_df$sun_acc[algo1_idx] >= 0.5, na.rm = TRUE), sum(algo1_idx)))
algo3_idx <- diag_df$sun_algo == "algo3" & diag_df$alpha == 1 & !is.na(diag_df$sun_acc)
cat(sprintf("  Algo3 alpha=1 acc>=0.5          : %d / %d cases\n",
            sum(diag_df$sun_acc[algo3_idx] >= 0.5, na.rm = TRUE), sum(algo3_idx)))
inflex_idx <- diag_df$sun_algo == "algo3" & diag_df$alpha > 1.1 & !is.na(diag_df$sun_algo3_used_inflex)
cat(sprintf("  Algo3 used_inflex (alpha>1.1)   : TRUE=%d FALSE=%d (of %d)\n",
            sum(diag_df$sun_algo3_used_inflex[inflex_idx], na.rm = TRUE),
            sum(!diag_df$sun_algo3_used_inflex[inflex_idx], na.rm = TRUE),
            sum(inflex_idx)))

# -----------------------------------------------------------------------
# Provenance CSV (single-row environment record)
# -----------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
si <- sessionInfo()
prov <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  platform = R.version$platform,
  os = si$running %||% R.version$os,
  mhn_version = as.character(utils::packageVersion("mhn")),
  bench_version = as.character(utils::packageVersion("bench")),
  mode = if (QUICK) "QUICK" else "FULL",
  iterations = BENCH_ITER,
  grid_points = nrow(diag_df),
  elapsed_min = round(elapsed_min, 2),
  stringsAsFactors = FALSE
)
write.csv(prov, PROV_CSV, row.names = FALSE)

cat(sprintf("\nResults written to:\n  %s\n  %s\n  %s\n  %s\n",
            RESULT_CSV, DIAG_CSV, REC_CSV, PROV_CSV))
cat("\nNext step: inspect the decision matrix above and update the auto path\n")
cat("in mhn/src/mhn_rmhn.cpp / R/rmhn.R.  The auto-vs-forced equivalence\n")
cat("regression block in tests/testthat/test-rmhn.R guards the dispatch rule.\n")
