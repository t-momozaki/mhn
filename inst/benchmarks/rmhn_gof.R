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
# and the four parameter regions (a)-(d) into which the package divides the
# RTDR construction of Gao & Wang (2025), served by three setups (A, BC and
# D), and for each applicable method in
# {auto, rtdr, sun}, it draws N variates, under a seed of its own for each
# cell, and records:
#   - the one-sample Kolmogorov-Smirnov statistic and p-value against
#     the package CDF pmhn() (a different code path: the samplers use
#     rejection kernels in which the Fox-Wright Psi cancels, whereas pmhn()
#     evaluates the Lemma 1b series or a quadrature; both are the package's
#     own, and its tests check pmhn() against independent references);
#   - the standardised mean error mean_z = (mean(x) - E[X]) / sqrt(Var[X]/N),
#     using mhn_mean() / mhn_var();
#   - the relative variance error |var(x) - Var[X]| / Var[X].
# The verdicts follow the package's own sampler tests (KS p > 0.001,
# |mean error| < 5 SE) except for the variance, which is judged more loosely
# here: within four standard deviations of the sample variance, with the
# excess kurtosis from mhn_kurtosis() and a floor of 0.05, where the tests
# allow five with the kurtosis from a reference independent of the package.
# They are written to a CSV so that they can be summarised in a table.
#
# method = "sun" is unavailable for alpha < 1 & gamma > 0 (Sun's
# Algorithm 2 is intentionally not implemented); those cells are
# recorded with error = "sun_algo2_not_implemented".  Special cases
# (gamma = 0, or alpha = 1) are intercepted before dispatch for every
# method, so all methods share the closed-form draw there.
#
# Outputs:
#   rmhn_gof_<YYYYMMDD>.csv
#   rmhn_gof_diagnostics_<YYYYMMDD>.csv  (one row: headline counts and the
#       software and machine of the run, from provenance.R)
#
# Invocation:
#   Rscript <file>
# where <file> is system.file("benchmarks", "rmhn_gof.R", package = "mhn").
# The results are written to MHN_GOF_OUTDIR when it names a directory, and not
# at all otherwise.
# Next to them it records the software and machine of the run
# (provenance.R), for which it runs `R CMD config`, the compiler with
# `--version`, and `sysctl` on macOS or `lscpu` on Linux.
#
# Optional environment variables:
#   MHN_GOF_QUICK=1   reduce the grid and N for a smoke test (about a second)
#   MHN_GOF_N=N       override the sample size per cell (default 10000)
#   MHN_GOF_SEED=S    override the RNG seed (default 1)
#   MHN_GOF_OUTDIR=D  override the output directory
#   MHN_GIT_COMMIT=H  git commit of the installed package sources, recorded
#                     in the diagnostics CSV

suppressPackageStartupMessages({
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("rmhn_gof.R requires the 'mhn' package; install it first.")
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

# Grid spanning the special cases and the package's RTDR regions (a)-(d):
#   gamma = 0            -> sqrt-Gamma special case (all alpha)
#   alpha = 1            -> truncated-normal special case (all gamma)
#   alpha > 1, gamma!=0  -> region (a): Sun Algorithm 1 (gamma>0) /
#                           Algorithm 3 (gamma<0)
#   0.5 <= alpha < 1     -> region (b)
#   alpha < 0.5, gamma<=thr -> region (c); alpha < 0.5, gamma>thr -> region (d)
ALPHAS <- if (QUICK) c(0.3, 1, 5) else c(0.05, 0.3, 0.7, 1, 1.5, 3, 10, 100)
# The grid reaches well past |gamma| = 10 and beta = 1.  A grid stopping there
# would miss a failure that shows only at a large tilt, such as a region (d)
# envelope piece whose area, formed directly, overflows between gamma = 50 and
# 100, or only away from beta = 1, such as an Algorithm 3 starting point that
# drops a factor of sqrt(beta).  An audit is only as wide as its grid.
GAMMAS <- if (QUICK) c(-2, 0, 5)  else
  c(-10000, -1000, -100, -10, -2, 0, 2, 10, 100, 1000, 10000)
BETAS  <- if (QUICK) 1.0 else c(0.01, 1, 100)
METHODS <- c("auto", "rtdr", "sun")

OUTDIR <- Sys.getenv("MHN_GOF_OUTDIR", unset = "")
# The results are written only where MHN_GOF_OUTDIR names a directory; without
# it the script prints its summary and writes nothing.
SAVE <- nzchar(OUTDIR)
if (SAVE) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, path) {
  if (SAVE) utils::write.csv(x, path, row.names = FALSE)
}
TODAY <- format(Sys.Date(), "%Y%m%d")
RESULT_CSV <- file.path(OUTDIR, sprintf("rmhn_gof_%s.csv", TODAY))
DIAG_CSV   <- file.path(OUTDIR, sprintf("rmhn_gof_diagnostics_%s.csv", TODAY))

cat(sprintf("[rmhn_gof] mode=%s N=%d alpha=%d gamma=%d methods=%d cells=%d\n",
            if (QUICK) "QUICK" else "FULL", N,
            length(ALPHAS), length(GAMMAS), length(METHODS),
            length(ALPHAS) * length(BETAS) * length(GAMMAS) * length(METHODS)))
cat(if (SAVE) sprintf("[rmhn_gof] result_csv=%s\n", RESULT_CSV) else
      "[rmhn_gof] results not saved (set MHN_GOF_OUTDIR to save them)\n")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Sun Algorithm 2 (gamma > 0, any alpha > 0) is not implemented; only for
# alpha < 1 & gamma > 0 does the Sun method then have no sampler.
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

grid <- expand.grid(alpha = ALPHAS, beta = BETAS, gamma = GAMMAS,
                    method = METHODS,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
n_cells <- nrow(grid)
rows <- vector("list", n_cells)
t0 <- Sys.time()

for (k in seq_len(n_cells)) {
  a <- grid$alpha[k]
  b <- grid$beta[k]
  g <- grid$gamma[k]
  m <- grid$method[k]
  # A seed for each cell, so that the tests of different cells do not share
  # their random numbers.
  rows[[k]] <- gof_one(a, b, g, m, seed = SEED + k - 1L)
  r <- rows[[k]]
  cat(sprintf("[%4d/%4d] a=%-6g b=%-6g g=%-8g %-4s ", k, n_cells, a, b, g, m))
  if (!is.na(r$error)) {
    cat(sprintf("%s\n", r$error))
  } else {
    cat(sprintf("KS.p=%.3f mean_z=%+.2f var_rel=%.3f\n",
                r$ks_pvalue, r$mean_z, r$var_rel_err))
  }
  # The full grid runs for the better part of an hour; without this the progress
  # lines sit in R's buffer until the end whenever stdout is not a terminal.
  utils::flush.console()
}

results <- do.call(rbind, rows)

# The variance tolerance travels with the results, since it is part of the
# pass/fail record and a reader of the CSV cannot reconstruct it otherwise.
# It depends only on the parameter triple and n, not on the draws.
results$var_kurt <- mapply(function(a, b, g)
  tryCatch(mhn::mhn_kurtosis(a, b, g), error = function(e) NA_real_),
  results$alpha, results$beta, results$gamma)
results$var_tol <- pmax(4 * sqrt(pmax(results$var_kurt + 2, 2) / results$n), 0.05)

save_csv(results, RESULT_CSV)

# -----------------------------------------------------------------------
# Headline diagnostics
# -----------------------------------------------------------------------
ok <- is.na(results$error)
ks_fail   <- ok & results$ks_pvalue < 0.001
mean_fail <- ok & abs(results$mean_z) > 5

# The variance tolerance has to follow the distribution, not a fixed percentage.
# A sample variance from n draws has relative standard deviation
# sqrt((excess kurtosis + 2)/n), which at alpha = 0.05 with a negative tilt is
# 0.11 -- so a flat 10 percent flags a dozen cells that are doing nothing wrong,
# and a grid wide enough to reach them cannot use it.  Four of those standard
# deviations, with the excess kurtosis taken as at least zero, give
# 4 sqrt(2/n) = 0.057 at n = 10^4 wherever the tail is no heavier than a normal
# one, which is most of the grid and tighter than ten per cent; the bound
# reaches 0.113 at alpha = 1 under a strong negative tilt, and exceeds that only
# where the tail is heavy enough to need it.  The floor of 0.05 binds only for
# n above 12800.
var_fail  <- ok & !is.na(results$var_tol) & results$var_rel_err > results$var_tol

# Precompute the headline scalars once (shared by the console summary
# and the diagnostics CSV below).  `safe()` leaves its argument a
# promise, so the min()/max() over an empty set is never evaluated when
# no cell was scored.
safe <- function(v) if (any(ok)) v else NA_real_
n_eval   <- sum(ok)
min_ks_p <- safe(min(results$ks_pvalue[ok], na.rm = TRUE))
max_mean <- safe(max(abs(results$mean_z[ok]), na.rm = TRUE))
max_var  <- safe(max(results$var_rel_err[ok], na.rm = TRUE))
# Under a correct sampler the KS p-values of the cells, each drawn under its
# own seed, are about uniform on (0, 1); a KS test of them against that law
# shows a drift no single cell resolves.
p_uniform <- safe(suppressWarnings(stats::ks.test(results$ks_pvalue[ok],
                                                  "punif")$p.value))

cat(sprintf("\n[rmhn_gof] summary:\n"))
cat(sprintf("  cells evaluated (excl. sun N/A): %d\n", n_eval))
cat(sprintf("  KS p < 0.001                   : %d  [min p = %.4f]\n",
            sum(ks_fail, na.rm = TRUE), min_ks_p))
cat(sprintf("  |mean_z| > 5                   : %d  [max = %.2f]\n",
            sum(mean_fail, na.rm = TRUE), max_mean))
cat(sprintf("  var_rel_err over its tolerance : %d  [max err = %.3f, max err/tol = %.3f]\n",
            sum(var_fail, na.rm = TRUE), max_var,
            safe(max(results$var_rel_err[ok] / results$var_tol[ok], na.rm = TRUE))))
cat(sprintf("  KS p-values against uniform    : p = %.4f\n", p_uniform))

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
  ks_pvalue_uniformity = p_uniform,
  stringsAsFactors = FALSE
)
diag <- with_provenance(diag, RUN)
save_csv(diag, DIAG_CSV)

if (SAVE) {
  cat(sprintf("\n[rmhn_gof] wrote %s\n", RESULT_CSV))
  cat(sprintf("[rmhn_gof] wrote %s\n", DIAG_CSV))
}
cat(sprintf("[rmhn_gof] total time = %.1fs\n", total_secs))
