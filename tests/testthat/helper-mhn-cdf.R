# References for the tests of the samplers, formed in R from the kernel of the
# law alone, k(u) = u^(alpha - 1) exp(-u^2 + z u) with u = sqrt(beta) x and
# z = gamma / sqrt(beta) (Sun et al. 2023, Eq. (1)), so that they call no
# function of the package.  The kernel is scaled by its value at a point near
# its peak, so that it neither overflows nor underflows, and integrated piece
# by piece around that point, where the mass lies.  Each piece is asked for a
# relative accuracy: the absolute 1e-10 that stats::integrate() also accepts
# by default bounds nothing for a law whose width is far below 1, and left
# the moments of gamma-like laws (z far below 0) wrong by up to 1e-6.
mhn_ref <- function(alpha, beta, gamma) {
  sb <- sqrt(beta)
  z <- gamma / sb
  log_k <- function(u) (alpha - 1) * log(u) - u^2 + z * u
  # The positive root of 2 u^2 - z u - (alpha - 1) = 0, written for z < 0
  # without the difference of two nearly equal terms.  For alpha <= 1 and
  # z <= 0 the kernel falls from u = 0 on a scale of 1 / (1 - z), and is
  # scaled at that point instead.
  peak <- if (alpha > 1) {
    d <- sqrt(z^2 + 8 * (alpha - 1))
    if (z >= 0) (z + d) / 4 else 2 * (alpha - 1) / (d - z)
  } else if (z > 0) {
    z / 2
  } else {
    1 / (1 - z)
  }
  # For alpha > 1 the peak is the mode, where (alpha - 1) / peak = 2 peak - z,
  # and log k(u) - log k(peak) = (alpha - 1) (log1p(d / peak) - d / peak) - d^2
  # with d = u - peak: no two large terms cancel, however narrow the law.
  log1pmx <- function(x) {
    ifelse(abs(x) < 1e-4, -x * x * (0.5 - x / 3 + x * x / 4), log1p(x) - x)
  }
  top <- log_k(peak)
  log_kn <- if (alpha > 1) {
    function(u) {
      d <- u - peak
      (alpha - 1) * log1pmx(d / peak) - d^2
    }
  } else {
    function(u) log_k(u) - top
  }
  # A width of the peak on the u axis, from the curvature of log k there, or
  # for alpha <= 1 and z <= 0 the scale on which the kernel falls.
  width <- if (alpha > 1) {
    1 / sqrt(2 + (alpha - 1) / peak^2)
  } else {
    1 / (sqrt(2) + max(0, -z))
  }
  # The ends of the range: from 10 widths either side of the peak, the
  # distance is doubled until log k is below -750, where what lies beyond is
  # below the smallest double, or on the left until u reaches 0.  A
  # gamma-like law (z far below 0) falls by only a few units of log k per
  # width on its right.
  right <- peak + 10 * width
  right_pts <- right
  while (log_kn(right) > -750) {
    right <- peak + 2 * (right - peak)
    right_pts <- c(right_pts, right)
  }
  left <- peak - 10 * width
  left_pts <- numeric(0)
  while (left > 0 && log_kn(left) > -750) {
    left_pts <- c(left_pts, left)
    left <- peak - 2 * (peak - left)
  }
  left <- max(0, left)
  breaks <- sort(unique(c(left, left_pts, pmax(left, peak + c(-3, 0, 3) * width),
                          right_pts)))
  # The integral of h(u) k(u), to a relative 1e-12 (1e-10 where integrate()
  # reports that rounding keeps it from that), or to 1e-16 of a width where
  # its value is smaller than that (as in the far tails of a narrow peak,
  # where the spacing of the doubles in u is a sizeable fraction of a width
  # and the kernel too ragged for a relative bound).  Where h changes
  # sign and the integral can be near 0, abs_tol is set by the caller.  For
  # alpha < 1 the piece that starts at 0 is taken in t = u^alpha, where
  # u^(alpha - 1) du = dt / alpha, so that the singularity at 0 is gone.
  integral <- function(h, lo, hi, abs_tol = 1e-16 * width) {
    pts <- c(lo, breaks[breaks > lo & breaks < hi], hi)
    sum(vapply(seq_len(length(pts) - 1L), function(i) {
      a <- pts[i]
      b <- pts[i + 1L]
      if (a == 0 && alpha < 1) {
        f <- function(t) {
          u <- t^(1 / alpha)
          h(u) * exp(-u^2 + z * u - top) / alpha
        }
        a <- 0
        b <- b^alpha
      } else {
        f <- function(u) h(u) * exp(log_kn(u))
      }
      piece <- function(rel_tol) {
        stats::integrate(f, a, b, rel.tol = rel_tol, abs.tol = abs_tol,
                         subdivisions = 1000L)$value
      }
      tryCatch(piece(1e-12), error = function(e) piece(1e-10))
    }, numeric(1)))
  }
  one <- function(u) rep(1, length(u))
  total <- integral(one, left, right)
  # Measured from the peak, so that a narrow law far from 0 keeps its digits;
  # the error allowed is 1e-13 of a width, or of a standard deviation to the
  # power p for the central moments.
  mean_u <- peak + integral(function(u) u - peak, left, right,
                            1e-13 * total * width) / total
  central <- function(p, scale) {
    integral(function(u) (u - mean_u)^p, left, right,
             1e-13 * total * scale^p) / total
  }
  var_u <- central(2, width)
  sd_u <- sqrt(var_u)
  cdf <- function(qs) {
    vapply(qs, function(q) {
      if (!(q > 0)) return(0)
      u <- q * sb
      if (u <= left) 0
      else if (u <= peak) integral(one, left, u) / total
      else if (u < right) 1 - integral(one, u, right) / total
      else 1
    }, numeric(1))
  }
  list(cdf = cdf, mean = mean_u / sb, var = var_u / beta,
       skewness = central(3, sd_u) / var_u^1.5,
       excess_kurtosis = central(4, sd_u) / var_u^2 - 3, sd_u = sd_u,
       mean_u = mean_u, sqrt_beta = sb)
}

# The reference CDF on its own, as a fixed function of q for ks.test().
mhn_cdf <- function(alpha, beta, gamma) mhn_ref(alpha, beta, gamma)$cdf

# Pearson's chi-square test of draws x against the reference, on `bins` bins
# whose edges lie near the reference quantiles j / bins.  The edges are read
# off the reference CDF on a grid and need not be exact: the expected counts
# use the CDF at the edges themselves.
chisq_ref <- function(x, ref, bins = 100L) {
  grid_u <- c(ref$mean_u * 10^seq(-6, 0, length.out = 200),
              ref$mean_u + ref$sd_u * seq(0, 12, length.out = 400)[-1])
  grid_u <- grid_u[grid_u > 0]
  grid <- grid_u / ref$sqrt_beta
  f_grid <- ref$cdf(grid)
  keep <- !duplicated(f_grid)
  edges <- stats::approx(f_grid[keep], grid[keep],
                         xout = seq_len(bins - 1L) / bins, rule = 2)$y
  edges <- unique(edges)
  p <- diff(c(0, ref$cdf(edges), 1))
  observed <- tabulate(findInterval(x, edges) + 1L, nbins = length(p))
  expected <- length(x) * p
  stat <- sum((observed - expected)^2 / expected)
  stats::pchisq(stat, df = length(p) - 1L, lower.tail = FALSE)
}

# The sample variance against the reference, to `k` of its standard errors:
# the variance of a sample variance is about var^2 (excess kurtosis + 2) / n.
expect_var_near <- function(sample_var, ref, n, k = 5, info = NULL) {
  se <- sqrt((ref$excess_kurtosis + 2) / n)
  testthat::expect_lt(abs(sample_var / ref$var - 1), k * se, label = info)
}
