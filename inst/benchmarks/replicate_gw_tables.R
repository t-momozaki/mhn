#!/usr/bin/env Rscript
# replicate_gw_tables.R --- Repeat the timing experiment behind Tables 1 and
# 2 of Gao & Wang (2025): their design, on this machine.
#
# Source: Gao & Wang (2025), Tables 1 and 2.
#
# Original paper grid (Gao & Wang notation):
#   lambda in {0.01, 0.3, 0.8, 1.1, 10, 100, 10000}
#   beta   in {-10000, -100, -2, -0.5, 10, 100, 10000}
#
# Conversion to Sun notation used by this package (Sun's f propto
# x^(alpha-1) exp(-beta x^2 + gamma x) vs Gao & Wang's f propto
# x^(lambda-1) exp(-alpha_GW x^2 - beta_GW x), so):
#   alpha (Sun) = lambda     (GW)
#   gamma (Sun) = -beta      (GW)         -- SIGN FLIP
#   beta  (Sun) = 1                        -- paper sets GW alpha = 1
#
# So the reproduction grid is:
#   alpha in {0.01, 0.3, 0.8, 1.1, 10, 100, 10000}
#   gamma in {10000, 100, 2, 0.5, -10, -100, -10000}
#
# Patterns (Gao & Wang 2025, Section 5):
#   T1 : 1 variate per call x 10000 calls  (Table 1, "as in Gibbs sampler")
#   T2 : 5 variates per call x 2000 calls  (Table 2, "varying parameter case")
#
# Methods:
#   A1 = RTDR (Gao & Wang's proposed generator) -> mhn::rmhn(method = "rtdr")
#   A2 = Sun, Kong, Pal generator               -> mhn::rmhn(method = "sun")
#
# A2 limitation: Sun et al. (2023) Algorithm 2 (the only Sun
# algorithm that covers alpha < 1 & gamma > 0; Algorithm 1 requires
# alpha > 1) is not implemented in this package: Sun et al. (2023,
# Section 4.1.2) note that its running time becomes significant when alpha
# or gamma^2 / beta is large (the truncation index of its mixture weights,
# their Lemma 8, grows like gamma^2 / beta), which is the very regime the
# rest of this script measures.
# For those cells (alpha in {0.01, 0.3, 0.8} x gamma in {10000, 100, 2,
# 0.5}, 12 cells) the A2 column is recorded as NA with reason
# "sun_algo2_not_implemented".
#
# Output: replicate_gw_tables_<YYYYMMDD>.csv (and a _diagnostics_ file that
# records the software and machine of the run, from provenance.R)
# Columns: alpha, gamma, gw_lambda, gw_beta, pattern, method,
#          time_sec_median, time_sec_iqr, paper_value, error
#
# Invocation:
#   Rscript <file>
# where <file> is
# system.file("benchmarks", "replicate_gw_tables.R", package = "mhn").
# The results are written to MHN_GW_OUTDIR when it names a directory, and not
# at all otherwise.
# Next to them it records the software and machine of the run
# (provenance.R), for which it runs `R CMD config`, the compiler with
# `--version`, and `sysctl` on macOS or `lscpu` on Linux.
#
# Optional environment variables:
#   MHN_GW_QUICK=1  -> reduce grid (3x3) and iterations for smoke test
#   MHN_GW_ITER=N   -> override bench::mark iterations (default 50; 2 with
#                      MHN_GW_QUICK=1)
#   MHN_GW_OUTDIR   -> override output directory
#   MHN_GIT_COMMIT  -> git commit of the installed package sources, recorded
#                      in the diagnostics CSV

suppressPackageStartupMessages({
  if (!requireNamespace("bench", quietly = TRUE))
    stop("replicate_gw_tables.R requires 'bench'. install.packages('bench').")
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("replicate_gw_tables.R requires the 'mhn' package; install it first.")
})

QUICK <- isTRUE(nchar(Sys.getenv("MHN_GW_QUICK")) > 0L)
# 50 iterations by default: at 5 the bench::mark IQR can approach the
# median, and at 20 the rtdr/sun ranking of near-tie cells still changes
# from one run to the next.  Even at 50 the relative IQR of a timing is
# often several per cent.
BENCH_ITER <- {
  v <- suppressWarnings(as.integer(Sys.getenv("MHN_GW_ITER", unset = "")))
  if (is.na(v) || v < 1L) (if (QUICK) 2L else 50L) else v
}

# The provenance helper ships beside this script.
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  here <- if (length(f)) dirname(gsub("~+~", " ", f[1L], fixed = TRUE)) else "."
  cand <- c(file.path(here, "provenance.R"),
            system.file("benchmarks", "provenance.R", package = "mhn"))
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (!length(cand))
    stop("cannot find provenance.R beside this script or in the installed mhn")
  source(cand[1L], local = FALSE)
})
RUN <- provenance_start()

# Full Gao & Wang grid (Sun notation).
GW_LAMBDAS <- if (QUICK) c(0.3, 1.1, 10) else c(0.01, 0.3, 0.8, 1.1, 10, 100, 10000)
GW_BETAS   <- if (QUICK) c(-100, -0.5, 100) else c(-10000, -100, -2, -0.5, 10, 100, 10000)
ALPHAS <- GW_LAMBDAS                       # Sun alpha = GW lambda
GAMMAS <- -GW_BETAS                        # Sun gamma = -GW beta (SIGN FLIP)

PATTERNS <- list(
  T1 = list(n_per_call = 1L, n_calls = if (QUICK) 1000L else 10000L,
            label = "Table 1 (1 var x 10000 call)"),
  T2 = list(n_per_call = 5L, n_calls = if (QUICK) 200L  else 2000L,
            label = "Table 2 (5 var x 2000 call)")
)

OUTDIR <- Sys.getenv("MHN_GW_OUTDIR", unset = "")
# The results are written only where MHN_GW_OUTDIR names a directory; without
# it the script prints its summary and writes nothing.
SAVE <- nzchar(OUTDIR)
if (SAVE) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, path) {
  if (SAVE) utils::write.csv(x, path, row.names = FALSE)
}
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("replicate_gw_tables_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR,
                        sprintf("replicate_gw_tables_diagnostics_%s.csv", TODAY))

cat(sprintf("[replicate_gw] mode=%s iterations=%d lambda=%d beta=%d patterns=%d points=%d\n",
            if (QUICK) "QUICK" else "FULL", BENCH_ITER,
            length(GW_LAMBDAS), length(GW_BETAS), length(PATTERNS),
            length(ALPHAS) * length(GAMMAS)))
cat(if (SAVE) sprintf("[replicate_gw] result_csv=%s\n", RESULT_CSV) else
      "[replicate_gw] results not saved (set MHN_GW_OUTDIR to save them)\n")

# ----- Reference values from Gao & Wang (2025) Tables 1 and 2 ---------------
# Indexed by (gw_lambda, gw_beta). NA where not in paper.
# Source: Gao & Wang (2025), Tables 1 and 2.
PAPER_T1 <- list(
  # rows = lambda, cols = beta (GW notation)
  "0.01"    = c(`-10000` = 0.46, `-100` = 0.36, `-2` = 0.50, `-0.5` = 0.44,
                `10` = 0.36, `100` = 0.35, `10000` = 0.34),
  "0.3"     = c(`-10000` = 0.39, `-100` = 0.36, `-2` = 0.46, `-0.5` = 0.31,
                `10` = 0.28, `100` = 0.27, `10000` = 0.27),
  "0.8"     = c(`-10000` = 0.31, `-100` = 0.30, `-2` = 0.28, `-0.5` = 0.28,
                `10` = 0.25, `100` = 0.25, `10000` = 0.27),
  "1.1"     = c(`-10000` = 0.32, `-100` = 0.25, `-2` = 0.27, `-0.5` = 0.24,
                `10` = 0.23, `100` = 0.27, `10000` = 0.25),
  "10"      = c(`-10000` = 0.29, `-100` = 0.26, `-2` = 0.28, `-0.5` = 0.28,
                `10` = 0.25, `100` = 0.26, `10000` = 0.25),
  "100"     = c(`-10000` = 0.25, `-100` = 0.28, `-2` = 0.27, `-0.5` = 0.27,
                `10` = 0.26, `100` = 0.25, `10000` = 0.27),
  "10000"   = c(`-10000` = 0.25, `-100` = 0.28, `-2` = 0.27, `-0.5` = 0.28,
                `10` = 0.26, `100` = 0.26, `10000` = 0.35)
)
PAPER_T1_A2 <- list(
  "0.01"    = c(`-10000` = 12990.09, `-100` = 137.27, `-2` = 0.57, `-0.5` = 0.47,
                `10` = 0.31, `100` = 0.28, `10000` = 0.30),
  "0.3"     = c(`-10000` = 461.28, `-100` = 5.20, `-2` = 0.53, `-0.5` = 0.56,
                `10` = 0.28, `100` = 0.31, `10000` = 0.30),
  "0.8"     = c(`-10000` = 174.34, `-100` = 2.30, `-2` = 0.54, `-0.5` = 0.60,
                `10` = 0.30, `100` = 0.29, `10000` = 0.29),
  "1.1"     = c(`-10000` = 0.45, `-100` = 0.45, `-2` = 0.45, `-0.5` = 0.45,
                `10` = 0.30, `100` = 0.31, `10000` = 0.30),
  "10"      = c(`-10000` = 0.46, `-100` = 0.44, `-2` = 0.42, `-0.5` = 0.43,
                `10` = 0.31, `100` = 0.32, `10000` = 0.37),
  "100"     = c(`-10000` = 0.47, `-100` = 0.45, `-2` = 0.42, `-0.5` = 0.42,
                `10` = 0.33, `100` = 0.33, `10000` = 0.36),
  "10000"   = c(`-10000` = 0.45, `-100` = 0.45, `-2` = 0.44, `-0.5` = 0.44,
                `10` = 0.33, `100` = 0.33, `10000` = 0.31)
)
PAPER_T2 <- list(
  "0.01"    = c(`-10000` = 0.16, `-100` = 0.13, `-2` = 0.19, `-0.5` = 0.15,
                `10` = 0.09, `100` = 0.10, `10000` = 0.09),
  "0.3"     = c(`-10000` = 0.13, `-100` = 0.13, `-2` = 0.16, `-0.5` = 0.14,
                `10` = 0.08, `100` = 0.08, `10000` = 0.09),
  "0.8"     = c(`-10000` = 0.07, `-100` = 0.09, `-2` = 0.10, `-0.5` = 0.11,
                `10` = 0.08, `100` = 0.08, `10000` = 0.08),
  "1.1"     = c(`-10000` = 0.08, `-100` = 0.10, `-2` = 0.07, `-0.5` = 0.08,
                `10` = 0.06, `100` = 0.07, `10000` = 0.06),
  "10"      = c(`-10000` = 0.08, `-100` = 0.10, `-2` = 0.09, `-0.5` = 0.10,
                `10` = 0.09, `100` = 0.10, `10000` = 0.11),
  "100"     = c(`-10000` = 0.07, `-100` = 0.08, `-2` = 0.08, `-0.5` = 0.08,
                `10` = 0.11, `100` = 0.08, `10000` = 0.08),
  "10000"   = c(`-10000` = 0.08, `-100` = 0.08, `-2` = 0.09, `-0.5` = 0.08,
                `10` = 0.08, `100` = 0.09, `10000` = 0.07)
)
PAPER_T2_A2 <- list(
  "0.01"    = c(`-10000` = 13487.79, `-100` = 137.20, `-2` = 0.29, `-0.5` = 0.25,
                `10` = 0.06, `100` = 0.06, `10000` = 0.06),
  "0.3"     = c(`-10000` = 454.01, `-100` = 4.95, `-2` = 0.28, `-0.5` = 0.30,
                `10` = 0.07, `100` = 0.05, `10000` = 0.07),
  "0.8"     = c(`-10000` = 172.46, `-100` = 1.98, `-2` = 0.26, `-0.5` = 0.30,
                `10` = 0.06, `100` = 0.06, `10000` = 0.06),
  "1.1"     = c(`-10000` = 0.09, `-100` = 0.08, `-2` = 0.11, `-0.5` = 0.10,
                `10` = 0.08, `100` = 0.08, `10000` = 0.08),
  "10"      = c(`-10000` = 0.09, `-100` = 0.08, `-2` = 0.09, `-0.5` = 0.08,
                `10` = 0.06, `100` = 0.06, `10000` = 0.08),
  "100"     = c(`-10000` = 0.10, `-100` = 0.09, `-2` = 0.09, `-0.5` = 0.09,
                `10` = 0.08, `100` = 0.06, `10000` = 0.06),
  "10000"   = c(`-10000` = 0.09, `-100` = 0.10, `-2` = 0.08, `-0.5` = 0.09,
                `10` = 0.06, `100` = 0.07, `10000` = 0.07)
)

paper_value <- function(gw_lambda, gw_beta, pattern, method) {
  tab <- if (pattern == "T1" && method == "rtdr") PAPER_T1
         else if (pattern == "T1" && method == "sun") PAPER_T1_A2
         else if (pattern == "T2" && method == "rtdr") PAPER_T2
         else if (pattern == "T2" && method == "sun") PAPER_T2_A2
         else NULL
  if (is.null(tab)) return(NA_real_)
  row <- tab[[as.character(gw_lambda)]]
  if (is.null(row)) return(NA_real_)
  v <- row[[as.character(gw_beta)]]
  if (is.null(v)) return(NA_real_) else as.numeric(v)
}

sun_unavailable <- function(alpha, gamma) alpha < 1.0 && gamma > 0.0

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
  thunk <- make_thunk(alpha, 1.0, gamma, method,
                      pattern$n_per_call, pattern$n_calls)
  res <- tryCatch(
    bench::mark(thunk(), iterations = BENCH_ITER, check = FALSE,
                filter_gc = FALSE),
    error = function(e) structure(list(error_msg = conditionMessage(e)),
                                  class = "bench_failure")
  )
  if (inherits(res, "bench_failure")) {
    return(list(median = NA_real_, iqr = NA_real_, error = res$error_msg))
  }
  # bench_time as.numeric returns SECONDS.
  med <- as.numeric(res$median[[1]])
  times <- as.numeric(res$time[[1]])
  iqr <- if (length(times) >= 2L) IQR(times) else NA_real_
  list(median = med, iqr = iqr, error = NA_character_)
}

rows <- list()
t_start <- Sys.time()
total_pts <- length(ALPHAS) * length(GAMMAS) * length(PATTERNS)
idx <- 0L

for (i_alpha in seq_along(ALPHAS)) {
  alpha <- ALPHAS[i_alpha]
  gw_lambda <- GW_LAMBDAS[i_alpha]
  for (i_gamma in seq_along(GAMMAS)) {
    gamma <- GAMMAS[i_gamma]
    gw_beta <- GW_BETAS[i_gamma]
    for (pname in names(PATTERNS)) {
      idx <- idx + 1L
      patt <- PATTERNS[[pname]]
      cat(sprintf("[%3d/%3d] lambda=%g beta=%g (alpha=%g gamma=%g) %s ",
                  idx, total_pts, gw_lambda, gw_beta, alpha, gamma, pname))
      flush.console()

      r_rtdr <- bench_one(alpha, gamma, "rtdr", pname, patt)
      cat(sprintf("rtdr=%.3fs ", r_rtdr$median))
      flush.console()

      if (sun_unavailable(alpha, gamma)) {
        r_sun <- list(median = NA_real_, iqr = NA_real_,
                      error = "sun_algo2_not_implemented")
        cat("sun=NA\n")
      } else {
        r_sun <- bench_one(alpha, gamma, "sun", pname, patt)
        cat(sprintf("sun=%.3fs\n", r_sun$median))
      }

      rows[[length(rows) + 1L]] <- data.frame(
        alpha = alpha, gamma = gamma, gw_lambda = gw_lambda, gw_beta = gw_beta,
        pattern = pname, method = "rtdr",
        time_sec_median = r_rtdr$median, time_sec_iqr = r_rtdr$iqr,
        paper_value = paper_value(gw_lambda, gw_beta, pname, "rtdr"),
        error = r_rtdr$error, stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        alpha = alpha, gamma = gamma, gw_lambda = gw_lambda, gw_beta = gw_beta,
        pattern = pname, method = "sun",
        time_sec_median = r_sun$median, time_sec_iqr = r_sun$iqr,
        paper_value = paper_value(gw_lambda, gw_beta, pname, "sun"),
        error = r_sun$error, stringsAsFactors = FALSE
      )
    }
  }
}

df <- do.call(rbind, rows)
save_csv(df, RESULT_CSV)

t_end <- Sys.time()
elapsed_min <- as.numeric(difftime(t_end, t_start, units = "mins"))

cat("\n========== SUMMARY ==========\n")
cat(sprintf("Elapsed: %.1f min   (start=%s end=%s)\n",
            elapsed_min, format(t_start, "%H:%M:%S"), format(t_end, "%H:%M:%S")))

# The package's times are those of its C++ kernels on this machine, and the
# paper's those of Gao & Wang's R code on their hardware, so the two are not
# comparable as they stand.  What can be compared is the ratio of the two
# samplers' times within each source: RTDR over Sun here, with Gao & Wang's
# ratio in parentheses.
print_times <- function(df, pattern_label, method_code) {
  sub <- subset(df, pattern == pattern_label & method == method_code)
  cat(sprintf("\n[%s, method=%s] time_sec_median on this machine\n",
              pattern_label, toupper(method_code)))
  for (lam in GW_LAMBDAS) {
    line <- sprintf("  lambda=%-7g | ", lam)
    for (bet in GW_BETAS) {
      r <- sub[sub$gw_lambda == lam & sub$gw_beta == bet, ]
      ours <- if (nrow(r) == 0L) NA_real_ else r$time_sec_median[1]
      line <- paste0(line, sprintf("%8s ", if (is.na(ours)) "--" else
                                     sprintf("%.3f", ours)))
    }
    cat(line, "\n")
  }
}
print_ratios <- function(df, pattern_label) {
  cat(sprintf(paste0("\n[%s] time of RTDR / time of Sun on this machine,",
                     " Gao & Wang's ratio in parentheses\n"), pattern_label))
  cell <- function(method_code, lam, bet) {
    r <- df[df$pattern == pattern_label & df$method == method_code &
              df$gw_lambda == lam & df$gw_beta == bet, ]
    if (nrow(r) == 0L) c(NA_real_, NA_real_) else
      c(r$time_sec_median[1], r$paper_value[1])
  }
  for (lam in GW_LAMBDAS) {
    line <- sprintf("  lambda=%-7g | ", lam)
    for (bet in GW_BETAS) {
      rt <- cell("rtdr", lam, bet)
      sn <- cell("sun", lam, bet)
      fmt <- function(x) if (is.na(x)) "NA" else sprintf("%.2f", x)
      line <- paste0(line, sprintf("%6s(%-6s) ", fmt(rt[1] / sn[1]),
                                   fmt(rt[2] / sn[2])))
    }
    cat(line, "\n")
  }
}
for (pat in c("T1", "T2")) {
  print_times(df, pat, "rtdr")
  print_times(df, pat, "sun")
  print_ratios(df, pat)
}

# -----------------------------------------------------------------------
# Diagnostics CSV (single-row environment / provenance record)
# -----------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
si <- sessionInfo()
diag <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  platform = R.version$platform,
  os = si$running %||% R.version$os,
  mhn_version = as.character(utils::packageVersion("mhn")),
  bench_version = as.character(utils::packageVersion("bench")),
  mode = if (QUICK) "QUICK" else "FULL",
  iterations = BENCH_ITER,
  grid_points = nrow(df),
  elapsed_min = round(elapsed_min, 2),
  stringsAsFactors = FALSE
)
diag <- with_provenance(diag, RUN)
save_csv(diag, DIAG_CSV)

if (SAVE) {
  cat(sprintf("\nResults written to: %s\n", RESULT_CSV))
  cat(sprintf("Diagnostics written to: %s\n", DIAG_CSV))
}
