#!/usr/bin/env Rscript
# rmhn_gof.R --- Distributional goodness-of-fit audit of rmhn().
#
# Role: this is a *correctness audit*, not a timing benchmark.  The
# timing scripts (replicate_gw_tables.R, auto_dispatch.R) show that
# rmhn() is fast; this script shows that the samples it produces
# actually follow the target Modified Half-Normal law, across every
# dispatch region and for each of the three methods.
#
# For each (alpha, beta, gamma) on a grid that spans the special cases
# and the four RTDR regions, and for each applicable method in
# {auto, rtdr, sun}, it draws N variates and records:
#   - the one-sample Kolmogorov-Smirnov statistic and p-value against
#     the package CDF pmhn() (an independent code path: the samplers
#     use rejection kernels in which the Fox-Wright Psi cancels, whereas
#     pmhn() evaluates the Lemma 1b series / Boost quadrature);
#   - the standardised mean error mean_z = (mean(x) - E[X]) / sqrt(Var[X]/N),
#     using the closed-form mhn_mean() / mhn_var();
#   - the relative variance error |var(x) - Var[X]| / Var[X].
# These mirror the pass/fail assertions in tests/testthat/test-rmhn.R
# (KS p > 0.001, |mean error| < 5 SE, relative variance error < 0.10),
# promoted here to a CSV so they can be summarised in a table/figure.
#
# method = "sun" is unavailable for alpha < 1 & gamma > 0 (Sun's
# Algorithm 2 is intentionally not implemented); those cells are
# recorded with error = "sun_algo2_not_implemented".  Special cases
# (gamma = 0, or alpha = 1) are intercepted before dispatch for every
# method, so all methods share the closed-form draw there.
#
# Outputs:
#   mhn/inst/benchmarks/results/rmhn_gof_<YYYYMMDD>.csv
#   mhn/inst/benchmarks/results/rmhn_gof_diagnostics_<YYYYMMDD>.csv
#
# Invocation (from repository root):
#   Rscript mhn/inst/benchmarks/rmhn_gof.R
#
# Optional environment variables:
#   MHN_GOF_QUICK=1   reduce the grid and N for a smoke test (~10 s)
#   MHN_GOF_N=N       override the sample size per cell (default 10000)
#   MHN_GOF_SEED=S    override the RNG seed (default 1)
#   MHN_GOF_OUTDIR=D  override the output directory

suppressPackageStartupMessages({
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("rmhn_gof.R requires the 'mhn' package. Run R CMD INSTALL mhn first.")
})

QUICK <- isTRUE(nchar(Sys.getenv("MHN_GOF_QUICK")) > 0L)
N <- {
  v <- suppressWarnings(as.integer(Sys.getenv("MHN_GOF_N", unset = "")))
  if (is.na(v) || v < 100L) (if (QUICK) 2000L else 10000L) else v
}
SEED <- {
  v <- suppressWarnings(as.integer(Sys.getenv("MHN_GOF_SEED", unset = "")))
  if (is.na(v)) 1L else v
}

# Grid spanning the special cases and the four RTDR regions:
#   gamma = 0            -> sqrt-Gamma special case (all alpha)
#   alpha = 1            -> truncated-normal special case (all gamma)
#   alpha >= 1, gamma!=0 -> region (a): Sun A1 (gamma>0) / A3 (gamma<0)
#   0.5 <= alpha < 1     -> region (b)
#   alpha < 0.5, gamma<=thr -> region (c); alpha < 0.5, gamma>thr -> region (d)
ALPHAS <- if (QUICK) c(0.3, 1, 5) else c(0.3, 0.7, 1, 1.5, 3, 10, 100)
GAMMAS <- if (QUICK) c(-2, 0, 5)  else c(-10, -2, 0, 2, 10)
BETA   <- 1.0
METHODS <- c("auto", "rtdr", "sun")

OUTDIR <- Sys.getenv("MHN_GOF_OUTDIR", unset = "")
if (!nzchar(OUTDIR)) OUTDIR <- file.path("mhn", "inst", "benchmarks", "results")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("rmhn_gof_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR, sprintf("rmhn_gof_diagnostics_%s.csv", TODAY))

cat(sprintf("[rmhn_gof] mode=%s N=%d alpha=%d gamma=%d methods=%d cells=%d\n",
            if (QUICK) "QUICK" else "FULL", N,
            length(ALPHAS), length(GAMMAS), length(METHODS),
            length(ALPHAS) * length(GAMMAS) * length(METHODS)))
cat(sprintf("[rmhn_gof] result_csv=%s\n", RESULT_CSV))

`%||%` <- function(a, b) if (is.null(a)) b else a

# Sun Algorithm 2 (alpha < 1 & gamma > 0) is not implemented.
sun_unavailable <- function(alpha, gamma) alpha < 1.0 && gamma > 0.0

# One-sample KS + moment summary for a single (alpha, beta, gamma, method).
gof_one <- function(alpha, beta, gamma, method, n = N, seed = SEED) {
  if (method == "sun" && sun_unavailable(alpha, gamma)) {
    return(data.frame(
      alpha = alpha, beta = beta, gamma = gamma, method = method, n = n,
      ks_stat = NA_real_, ks_pvalue = NA_real_,
      sample_mean = NA_real_, theory_mean = NA_real_, mean_z = NA_real_,
      sample_var = NA_real_, theory_var = NA_real_, var_rel_err = NA_real_,
      error = "sun_algo2_not_implemented", stringsAsFactors = FALSE))
  }
  set.seed(seed)
  x <- tryCatch(mhn::rmhn(n, alpha, beta, gamma, method = method),
                error = function(e) structure(NA_real_,
                                              gof_error = conditionMessage(e)))
  if (length(x) == 1L && is.na(x)) {
    return(data.frame(
      alpha = alpha, beta = beta, gamma = gamma, method = method, n = n,
      ks_stat = NA_real_, ks_pvalue = NA_real_,
      sample_mean = NA_real_, theory_mean = NA_real_, mean_z = NA_real_,
      sample_var = NA_real_, theory_var = NA_real_, var_rel_err = NA_real_,
      error = attr(x, "gof_error") %||% "rmhn_error",
      stringsAsFactors = FALSE))
  }

  theory_mean <- mhn::mhn_mean(alpha, beta, gamma)
  theory_var  <- mhn::mhn_var(alpha, beta, gamma)
  sample_mean <- mean(x)
  sample_var  <- stats::var(x)

  cdf <- function(q) mhn::pmhn(q, alpha = alpha, beta = beta, gamma = gamma)
  ks  <- suppressWarnings(stats::ks.test(x, cdf))

  se_mean <- sqrt(theory_var / n)
  data.frame(
    alpha = alpha, beta = beta, gamma = gamma, method = method, n = n,
    ks_stat = unname(ks$statistic), ks_pvalue = ks$p.value,
    sample_mean = sample_mean, theory_mean = theory_mean,
    mean_z = (sample_mean - theory_mean) / se_mean,
    sample_var = sample_var, theory_var = theory_var,
    var_rel_err = abs(sample_var - theory_var) / theory_var,
    error = NA_character_, stringsAsFactors = FALSE)
}

grid <- expand.grid(alpha = ALPHAS, gamma = GAMMAS, method = METHODS,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
n_cells <- nrow(grid)
rows <- vector("list", n_cells)
t0 <- Sys.time()

for (k in seq_len(n_cells)) {
  a <- grid$alpha[k]
  g <- grid$gamma[k]
  m <- grid$method[k]
  rows[[k]] <- gof_one(a, BETA, g, m)
  r <- rows[[k]]
  cat(sprintf("[%3d/%3d] alpha=%-6g gamma=%-6g %-4s ", k, n_cells, a, g, m))
  if (!is.na(r$error)) {
    cat(sprintf("%s\n", r$error))
  } else {
    cat(sprintf("KS.p=%.3f mean_z=%+.2f var_rel=%.3f\n",
                r$ks_pvalue, r$mean_z, r$var_rel_err))
  }
}

results <- do.call(rbind, rows)
write.csv(results, RESULT_CSV, row.names = FALSE)

# -----------------------------------------------------------------------
# Headline diagnostics
# -----------------------------------------------------------------------
ok <- is.na(results$error)
ks_fail   <- ok & results$ks_pvalue < 0.001
mean_fail <- ok & abs(results$mean_z) > 5
var_fail  <- ok & results$var_rel_err > 0.10

# Precompute the headline scalars once (shared by the console summary
# and the diagnostics CSV below).  `safe()` leaves its argument a
# promise, so the min()/max() over an empty set is never evaluated when
# no cell was scored.
safe <- function(v) if (any(ok)) v else NA_real_
n_eval   <- sum(ok)
min_ks_p <- safe(min(results$ks_pvalue[ok], na.rm = TRUE))
max_mean <- safe(max(abs(results$mean_z[ok]), na.rm = TRUE))
max_var  <- safe(max(results$var_rel_err[ok], na.rm = TRUE))

cat(sprintf("\n[rmhn_gof] summary:\n"))
cat(sprintf("  cells evaluated (excl. sun N/A): %d\n", n_eval))
cat(sprintf("  KS p < 0.001                   : %d  [min p = %.4f]\n",
            sum(ks_fail, na.rm = TRUE), min_ks_p))
cat(sprintf("  |mean_z| > 5                   : %d  [max = %.2f]\n",
            sum(mean_fail, na.rm = TRUE), max_mean))
cat(sprintf("  var_rel_err > 0.10             : %d  [max = %.3f]\n",
            sum(var_fail, na.rm = TRUE), max_var))

# -----------------------------------------------------------------------
# Diagnostics CSV (single-row environment / provenance record)
# -----------------------------------------------------------------------
total_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
si <- sessionInfo()
diag <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  platform = R.version$platform,
  os = si$running %||% R.version$os,
  mhn_version = as.character(utils::packageVersion("mhn")),
  mode = if (QUICK) "QUICK" else "FULL",
  n = N, seed = SEED,
  cells = n_cells,
  cells_evaluated = n_eval,
  total_secs = round(total_secs, 2),
  ks_fail_count = sum(ks_fail, na.rm = TRUE),
  min_ks_pvalue = min_ks_p,
  mean_fail_count = sum(mean_fail, na.rm = TRUE),
  max_abs_mean_z = max_mean,
  var_fail_count = sum(var_fail, na.rm = TRUE),
  max_var_rel_err = max_var,
  stringsAsFactors = FALSE
)
write.csv(diag, DIAG_CSV, row.names = FALSE)

cat(sprintf("\n[rmhn_gof] wrote %s\n", RESULT_CSV))
cat(sprintf("[rmhn_gof] wrote %s\n", DIAG_CSV))
cat(sprintf("[rmhn_gof] total time = %.1fs\n", total_secs))
