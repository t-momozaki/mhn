#!/usr/bin/env Rscript
# vmf_gibbs.R --- von Mises-Fisher concentration Gibbs sampler via rmhn().
#
# Role: a self-contained worked example for Section 6 of the paper.  It
# shows that a single mhn::rmhn() call performs the concentration-
# parameter update of a fully Bayesian analysis of directional data, so
# that a complete Gibbs sampler for the von Mises-Fisher (vMF) model is a
# few lines of base R once the MHN draw is available.
#
# Model.  Observed unit vectors y_1, ..., y_n on the sphere S^{p-1} are
# modelled as i.i.d. von Mises-Fisher with known mean direction mu and
# unknown concentration kappa > 0:
#
#     f_VF(y | mu, kappa) = C_p(kappa) exp(kappa * mu'y),
#     C_p(kappa) = kappa^nu / ((2*pi)^{p/2} I_nu(kappa)),   nu = p/2 - 1,
#
# with I_nu the modified Bessel function of the first kind.  The prior on
# kappa is Gamma(a0, b0) (shape a0, rate b0).  The normalising constant
# 1/I_nu(kappa)^n makes the kappa likelihood intractable, so we follow
# the data-augmentation scheme of Sun, Kong & Pal (2023), Section 2.3
# (their directional-data example).  Latent T_1, ..., T_n > 0 with
#
#     f_{T,nu}(t | kappa) = c(kappa) * sum_{k>=1} w_k * exp(-(j_{nu,k}^2 + kappa^2) t),
#     w_k = j_{nu,k}^{nu+1} / J_{nu+1}(j_{nu,k}),
#
# where {j_{nu,k}} are the positive zeros of the Bessel function J_nu, are
# introduced.  The T_i absorb the intractable 1/I_nu(kappa) factor, and
# the resulting full conditional of kappa is a Modified Half-Normal:
#
#     pi(kappa | T, Y) proportional to
#         kappa^{a0 - 1} exp(-(n*Tbar) kappa^2 + (n*mu'Ybar - b0) kappa),
#
# i.e. kappa | . ~ MHN(alpha = a0, beta = n*Tbar, gamma = n*mu'Ybar - b0),
# with Tbar = mean(T_i) and Ybar = mean(y_i)  (Sun, Kong & Pal 2023,
# their Equations (1) and (5)).  Each sweep therefore reduces to
#
#     T_i  ~ f_{T,nu}(. | kappa)            (latent augmentation draws)
#     kappa <- rmhn(1, a0, n*Tbar, n*mu'Ybar - b0)   (one MHN draw)
#
# Reference:
#   Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#   distribution: Properties and an efficient sampling scheme.
#   Communications in Statistics - Theory and Methods, 52(5), 1591-1613.
#
# What it does: simulates a small spherical data set from a vMF law with a
# known true kappa and a fixed mean direction, runs the Gibbs sampler
# above, prints the posterior mean of kappa against the truth and a mixing
# summary, and (when run standalone) writes the posterior draws to a dated
# CSV under results/.
#
# Invocation (from the repository root):
#   Rscript mhn/inst/examples/vmf_gibbs.R
#
# Optional environment variables:
#   MHN_VMF_N=N       number of retained Gibbs draws (chain length; default 1500)
#   MHN_VMF_BURN=B    burn-in sweeps discarded before recording (default 500)
#   MHN_VMF_NDATA=n   number of observed unit vectors (default 300)
#   MHN_VMF_P=p       sphere dimension p (data live on S^{p-1}; default 3)
#   MHN_VMF_KAPPA=k   true concentration used to simulate the data (default 5)
#   MHN_VMF_TRUNC=K   number of Bessel-zero terms in the T draw (default 128)
#   MHN_VMF_GRID=G    grid size for the inverse-CDF T sampler (default 512)
#   MHN_VMF_SEED=S    RNG seed (default 1)
#   MHN_VMF_QUICK=1   smaller data set and shorter chain for a smoke test
#   MHN_VMF_OUTDIR=D  override the CSV output directory

suppressPackageStartupMessages({
  if (!requireNamespace("mhn", quietly = TRUE))
    stop("vmf_gibbs.R requires the 'mhn' package. Run R CMD INSTALL mhn first.")
})

## ---- configuration -------------------------------------------------------

getenv_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
  if (is.na(v)) default else v
}
getenv_num <- function(name, default) {
  v <- suppressWarnings(as.numeric(Sys.getenv(name, unset = "")))
  if (is.na(v)) default else v
}

QUICK      <- nzchar(Sys.getenv("MHN_VMF_QUICK"))
N_DRAWS    <- getenv_int("MHN_VMF_N",     if (QUICK) 400L else 1500L)
N_BURN     <- getenv_int("MHN_VMF_BURN",  if (QUICK) 150L else 500L)
N_DATA     <- getenv_int("MHN_VMF_NDATA", if (QUICK) 120L else 300L)
P_DIM      <- getenv_int("MHN_VMF_P",     3L)
KAPPA_TRUE <- getenv_num("MHN_VMF_KAPPA", 5.0)
N_TERMS    <- getenv_int("MHN_VMF_TRUNC", 128L)
N_GRID     <- getenv_int("MHN_VMF_GRID",  512L)
SEED       <- getenv_int("MHN_VMF_SEED",  1L)

# Weakly informative Gamma(shape, rate) prior on kappa (prior mean 100).
PRIOR_SHAPE <- 1.0   # a0  (= MHN alpha in the full conditional)
PRIOR_RATE  <- 0.01  # b0

stopifnot(P_DIM >= 2L, KAPPA_TRUE > 0, N_DATA >= 2L, N_TERMS >= 8L)

## ---- Bessel-zero utilities ----------------------------------------------

# First K positive zeros of the Bessel function J_nu, found by bracketing
# around the McMahon asymptotic guess j_{nu,k} ~ (k + nu/2 - 1/4) pi and
# refining with uniroot().  Consecutive zeros are ~pi apart, so a bracket
# of half that width around each guess isolates a single sign change.
bessel_zeros <- function(nu, K) {
  z <- numeric(K)
  for (k in seq_len(K)) {
    guess <- (k + nu / 2 - 0.25) * pi
    lo <- max(1e-8, guess - pi / 2 + 1e-6)
    hi <- guess + pi / 2 - 1e-6
    # Widen the bracket if the asymptotic guess is slightly off for small k.
    tries <- 0L
    while (besselJ(lo, nu) * besselJ(hi, nu) > 0 && tries < 8L) {
      lo <- max(1e-8, lo - pi / 4)
      hi <- hi + pi / 4
      tries <- tries + 1L
    }
    z[k] <- uniroot(function(x) besselJ(x, nu), c(lo, hi),
                    tol = 1e-10)$root
  }
  z
}

## ---- von Mises-Fisher sampler (Wood 1994) -------------------------------

# Householder reflection mapping the pole e_p = (0, ..., 0, 1) onto the
# unit mean direction mu; used to rotate draws generated about the pole.
reflect_pole_to <- function(mu) {
  p <- length(mu)
  e <- c(rep(0, p - 1), 1)
  u <- e - mu
  s <- sum(u * u)
  if (s < 1e-12) diag(p) else diag(p) - 2 * (u %o% u) / s
}

# Draw n unit vectors on S^{p-1} from vMF(mu, kappa) using Ulrich's method
# as given by Wood (1994): sample the cosine W to the mean direction by
# rejection, then a uniform direction on the orthogonal S^{p-2}.
rvmf <- function(n, mu, kappa) {
  p <- length(mu)
  b  <- (-2 * kappa + sqrt(4 * kappa^2 + (p - 1)^2)) / (p - 1)
  x0 <- (1 - b) / (1 + b)
  c0 <- kappa * x0 + (p - 1) * log(1 - x0^2)
  W <- numeric(n)
  for (i in seq_len(n)) {
    repeat {
      Z  <- rbeta(1, (p - 1) / 2, (p - 1) / 2)
      Wi <- (1 - (1 + b) * Z) / (1 - (1 - b) * Z)
      if (kappa * Wi + (p - 1) * log(1 - x0 * Wi) - c0 >= log(runif(1))) {
        W[i] <- Wi
        break
      }
    }
  }
  V <- matrix(rnorm(n * (p - 1)), n, p - 1)
  V <- V / sqrt(rowSums(V^2))
  X <- cbind(sqrt(1 - W^2) * V, W)   # draws about the pole e_p
  X %*% reflect_pole_to(mu)          # rotate the pole onto mu
}

## ---- latent augmentation: sampler for T | kappa -------------------------

# Precompute the Bessel zeros and the (signed) series weights w_k.  Only
# their relative values matter: the kappa-dependent normaliser cancels in
# the inverse-CDF below, as does the global constant 2^{(2+nu)/2}.
NU      <- P_DIM / 2 - 1
J_ZERO  <- bessel_zeros(NU, N_TERMS)
J_SQ    <- J_ZERO^2
WEIGHT  <- J_ZERO^(NU + 1) / besselJ(J_ZERO, NU + 1)

# Draw n independent T_i from f_{T,nu}(. | kappa) by numerical inversion.
# The survival function of that density is
#   S(t) proportional to sum_k w_k exp(-(j_k^2 + kappa^2) t) / (j_k^2 + kappa^2),
# which converges geometrically for t > 0.  We tabulate it on a log-spaced
# grid whose ends carry negligible mass, form the CDF by self-normalising,
# and invert it by linear interpolation.
draw_T <- function(n, kappa) {
  rates <- J_SQ + kappa^2
  r_min <- rates[1]                       # slowest-decaying term
  t_max <- 40 / r_min                      # right tail exhausted here
  t_min <- min(1e-3, t_max / 1e3)          # left tail is astronomically thin
  tgrid <- exp(seq(log(t_min), log(t_max), length.out = N_GRID))
  # Survival S(t) on the grid: (K x G) matrix of exp(-rate * t).
  coef <- WEIGHT / rates
  Sgrid <- as.numeric(coef %*% exp(-outer(rates, tgrid)))
  # CDF anchored at the grid ends; force monotonicity against tiny wiggles.
  Fgrid <- (Sgrid[1] - Sgrid) / (Sgrid[1] - Sgrid[N_GRID])
  Fgrid <- cummax(pmin(pmax(Fgrid, 0), 1))
  keep  <- c(TRUE, diff(Fgrid) > 0)        # strictly increasing support
  approx(Fgrid[keep], tgrid[keep], xout = runif(n), rule = 2)$y
}

## ---- Gibbs sampler -------------------------------------------------------

# Runs the augmentation Gibbs sampler and returns the recorded draws.
run_vmf_gibbs <- function(Y, mu, prior_shape, prior_rate,
                          n_draws, n_burn) {
  n     <- nrow(Y)
  y_sum <- colSums(Y)
  p     <- ncol(Y)
  proj  <- sum(y_sum * mu)                  # n * mu'Ybar
  total <- n_draws + n_burn

  # Start at the Banerjee et al. (2005) moment estimate of kappa from the
  # mean resultant length, so the chain begins in the region of the
  # posterior rather than in the diffuse prior tail.
  r_bar <- sqrt(sum(y_sum^2)) / n
  r_bar <- min(r_bar, 1 - 1e-8)
  kappa <- max(1e-3, r_bar * (p - r_bar^2) / (1 - r_bar^2))
  draws <- data.frame(sweep = integer(n_draws), kappa = numeric(n_draws),
                      Tbar = numeric(n_draws), beta_k = numeric(n_draws),
                      gamma_k = numeric(n_draws))

  for (s in seq_len(total)) {
    T_lat   <- draw_T(n, kappa)             # latent augmentation draws
    Tbar    <- mean(T_lat)
    beta_k  <- n * Tbar                     # MHN beta
    gamma_k <- proj - prior_rate            # MHN gamma = n*mu'Ybar - b0
    kappa_new <- mhn::rmhn(1, alpha = prior_shape,
                           beta = beta_k, gamma = gamma_k)  # one MHN draw
    if (is.finite(kappa_new)) kappa <- kappa_new           # retain on failure
    if (s > n_burn) {
      i <- s - n_burn
      draws[i, ] <- list(i, kappa, Tbar, beta_k, gamma_k)
    }
  }
  draws
}

## ---- effective sample size (mixing summary) -----------------------------

effective_size <- function(x) {
  n <- length(x)
  if (stats::var(x) <= 0) return(NA_real_)
  ac  <- as.numeric(acf(x, lag.max = min(n - 1L, 200L), plot = FALSE)$acf)
  rho <- ac[-1]                              # autocorrelations at lags >= 1
  # Initial-positive-sequence truncation (Geyer 1992).
  cutoff <- which(rho < 0)[1]
  if (!is.na(cutoff) && cutoff > 1) rho <- rho[seq_len(cutoff - 1L)]
  n / (1 + 2 * sum(rho))
}

## ---- run -----------------------------------------------------------------

set.seed(SEED)

# Fixed mean direction (an arbitrary unit vector) and simulated data.
mu <- rep(1, P_DIM)
mu <- mu / sqrt(sum(mu * mu))
Y  <- rvmf(N_DATA, mu, KAPPA_TRUE)

cat(sprintf("[vmf_gibbs] %s p=%d n_data=%d kappa_true=%.3f draws=%d burn=%d terms=%d\n",
            if (QUICK) "QUICK" else "FULL", P_DIM, N_DATA, KAPPA_TRUE,
            N_DRAWS, N_BURN, N_TERMS))

post <- run_vmf_gibbs(Y, mu, PRIOR_SHAPE, PRIOR_RATE, N_DRAWS, N_BURN)

k_mean <- mean(post$kappa)
k_sd   <- sd(post$kappa)
k_ci   <- quantile(post$kappa, c(0.025, 0.975))
ess    <- effective_size(post$kappa)
lag1   <- cor(post$kappa[-1], post$kappa[-N_DRAWS])

cat("\n--- posterior summary for kappa ---\n")
cat(sprintf("  true kappa        : %.3f\n", KAPPA_TRUE))
cat(sprintf("  posterior mean    : %.3f\n", k_mean))
cat(sprintf("  posterior sd      : %.3f\n", k_sd))
cat(sprintf("  95%% credible int. : (%.3f, %.3f)\n", k_ci[1], k_ci[2]))
cat(sprintf("  covers truth      : %s\n",
            if (KAPPA_TRUE >= k_ci[1] && KAPPA_TRUE <= k_ci[2]) "yes" else "no"))
cat("\n--- mixing summary ---\n")
cat(sprintf("  retained draws    : %d\n", N_DRAWS))
cat(sprintf("  lag-1 autocorr.   : %.3f\n", lag1))
cat(sprintf("  effective size    : %.1f  (%.1f%% of draws)\n",
            ess, 100 * ess / N_DRAWS))

## ---- persist draws when run standalone ----------------------------------

# TRUE under `Rscript vmf_gibbs.R`, FALSE when source()'d interactively.
run_standalone <- !interactive() && sys.nframe() == 0L
if (run_standalone) {
  outdir <- Sys.getenv("MHN_VMF_OUTDIR", unset = "")
  if (!nzchar(outdir)) {
    here <- tryCatch(dirname(sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                error = function(e) ".")
    if (is.na(here) || !nzchar(here)) here <- "."
    outdir <- file.path(here, "results")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(outdir, sprintf("vmf_gibbs_%s.csv", format(Sys.Date(), "%Y%m%d")))
  write.csv(post, csv, row.names = FALSE)
  cat(sprintf("\n[vmf_gibbs] wrote %d draws to %s\n", nrow(post), csv))
}

invisible(post)
