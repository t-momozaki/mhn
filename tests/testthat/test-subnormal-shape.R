# A subnormal shape, alpha below 2^-1022 (src/mhn_vanishing.h).
#
# There t^alpha = 1 to double precision at every positive double t, and the
# law takes a form in which alpha only multiplies: with z = gamma / sqrt(beta)
# and u = sqrt(beta) x,
#   f(x) = alpha exp(-beta x^2 + gamma x) / (x (1 + w)),
#   P(X > x) = alpha J(u) / (1 + w),  J(u) = int_u^Inf exp(-t^2 + z t) / t dt,
#   P(X <= x) = (1 + alpha G(u)) / (1 + w),  G(u) = K0(z) - J(u),
#   E (sqrt(beta) X)^k = alpha A_k / (1 + w),
#   A_k = int_0^Inf t^(k-1) exp(-t^2 + z t) dt,
# where w = alpha K0(z) is below 1e-200 for z <= 30.  1 + alpha G(u) is
# 1 + alpha I(u), I(u) = int_0^u (exp(-t^2 + z t) - 1) / t dt, to double
# precision.  The references below are these integrals by integrate(), the
# law at gamma = 0 (Sun et al. 2023, Lemma 6a), and the law at the smallest
# normal shape, 2^-1022, which the ordinary paths serve.

u1 <- 2^-1074
# integrate() over one wide interval can miss a kernel that dies within a
# small part of it (at z = -30 over [1, 7] it is off by 5e-6), so the
# references integrate over pieces that start at the scale of the kernel,
# 1 / (1 + |z|), and double in width until they pass its peak by 40.
pieces <- function(f, lo, z) {
  sc <- 1 / (1 + abs(z))
  ends <- lo + sc * 2^(0:60)
  ends <- ends[seq_len(which(ends > max(lo, z / 2) + 40)[1])]
  starts <- c(lo, ends[-length(ends)])
  sum(mapply(function(a, b) integrate(f, a, b, rel.tol = 1e-13)$value,
             starts, ends))
}
J_ref <- function(x, z) pieces(function(t) exp(-t^2 + z * t) / t, x, z)
A_ref <- function(k, z) pieces(function(t) t^(k - 1) * exp(-t^2 + z * t), 0, z)
# J(x) for x < 1, as J(1) + int_x^1 (exp(-t^2 + z t) - 1) / t dt - log x,
# whose integrand stays bounded at the origin.
J_small <- function(x, z) {
  J_ref(1, z) - log(x) +
    integrate(function(t) ifelse(t < 1e-300, z, expm1(t * (z - t)) / t),
              x, 1, rel.tol = 1e-13)$value
}

# For a tilt above 30 the kernel passes the largest double, and the
# references are taken on the log scale with exp(z^2 / 4) out, in pieces
# about the peak at z / 2.
around <- function(f, m) {
  ends <- c(max(m - 40, 0), m - c(16, 8, 4, 2, 1, 0) , m + c(1, 2, 4, 8, 16, 40))
  ends <- sort(unique(pmax(ends, 0)))
  sum(mapply(function(a, b) integrate(f, a, b, rel.tol = 1e-13)$value,
             ends[-length(ends)], ends[-1]))
}
# log K0 = log D (the Euler term is below 1e-96 of D here).
logK0_ref <- function(z) {
  m <- z / 2
  m^2 + log(around(function(t) exp(-(t - m)^2) * (-expm1(-z * t)) / t, m))
}
logA_ref <- function(k, z) {
  m <- z / 2
  m^2 + log(around(function(t) t^(k - 1) * exp(-(t - m)^2), m))
}
# log I(u) for u below z / 2, with exp(phi(u)) out, in pieces toward u.
logI_ref <- function(u, z) {
  pu <- u * (z - u)
  ends <- u - c(256, 128, 64, 32, 16, 8, 4, 2, 1, 0.5, 0.25, 0) / (z - 2 * u)
  ends <- sort(unique(c(0, pmax(ends, 0))))
  f <- function(t) {
    ifelse(t < 1e-300, z * exp(-pu), (exp(t * (z - t) - pu) - exp(-pu)) / t)
  }
  pu + log(sum(mapply(function(a, b) integrate(f, a, b, rel.tol = 1e-13)$value,
                      ends[-length(ends)], ends[-1])))
}
log1p_exp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))

test_that("density and upper tail are proportional to a subnormal shape", {
  # Half of an odd multiple of 2^-1074 is not a double, and the paths for an
  # ordinary shape carried its rounding: in 0.1.1, at 3 * 2^-1074, the
  # density and the upper tail were 4/3 of their values at every tilt tried
  # from 0.001 to 40, and at z = 0 the density was 4/3 and the upper tail
  # 1.28 times its value.  log J, of order one, is compared with log alpha
  # taken off, which would otherwise hide its error.
  for (k in c(1, 3, 5, 101, 2^40 + 1)) {
    for (z in c(-30, -1, 0, 0.001, 1, 12, 20, 21, 40)) {
      a <- k * u1
      info <- sprintf("alpha = %g * 2^-1074, z = %g", k, z)
      expect_close(dmhn(1, a, 1, z, log = TRUE), log(a) + z - 1, 1e-13,
                   info = info)
      expect_close(pmhn(0.5, a, 1, z, lower.tail = FALSE, log.p = TRUE) -
                     log(a), log(J_ref(0.5, z)), 1e-12, abs_tol = 1e-13,
                   info = info)
    }
  }
})

test_that("the lower tail of a subnormal shape is 1 to double precision", {
  # The upper tail is alpha J(u) / (1 + w), below 1e-130 here, so F is 1 and
  # log F is -alpha J(u), to within the spacing of the doubles there.  In
  # 0.1.1, at 2^-1074, pmhn() returned 0 at every x tried under z = -30 and
  # at x = 1e-300 under z = 0, NA at every x under z from 0.001 to 12, and NA
  # below x = 6 under z = 20 and 21 and below x = 1 under z = 40; at
  # 3 * 2^-1074 under z = 21, log F was 0.
  for (k in c(1, 2, 3, 7)) {
    for (z in c(-30, 0, 0.001, 1, 12, 20, 21, 40)) {
      a <- k * u1
      x <- c(1e-300, 1e-150, 1e-10, 0.5, 1, 6)
      info <- sprintf("alpha = %g * 2^-1074, z = %g", k, z)
      expect_identical(pmhn(x, a, 1, z), rep(1, length(x)), info = info)
      J <- vapply(x, function(xi) if (xi < 1) J_small(xi, z) else J_ref(xi, z),
                  numeric(1))
      lF <- pmhn(x, a, 1, z, log.p = TRUE)
      expect_true(all(abs(lF + a * J) <= 2 * u1 + 1e-12 * (a * J)), info = info)
    }
  }
  # +0, not -0, below the support.
  expect_identical(1 / pmhn(c(0, -1), 3 * u1, 1, 1), c(Inf, Inf))
  expect_identical(1 / pmhn(c(0, -1), u1, 1, 56), c(Inf, Inf))
})

test_that("the lower tail of a subnormal shape keeps the mass by the origin", {
  # Where the peak outweighs the origin (w > 1, from z of about 53) the lower
  # tail below the peak is (1 + alpha I(u)) / (1 + w), and where alpha I(u)
  # is small it is the mass by the origin, 1 / (1 + w), itself.  Past the
  # tilt where the forms hand over (55 to 56.3, as alpha falls) the points
  # with alpha I(u) up to 2^64 are still taken from them, which the points
  # where t (z - t) is 650 to 850 straddle at z = 60, 100 and 1000.  In 0.1.1
  # pmhn(1, 2^-1074, 1, 56) was 1 against 1.04e-16, log F at z = 57, x = 1
  # was 0 against -65.03, and at 3 * 2^-1074 the mass by the origin was taken
  # at 3/2 of its size (log F -65.726 against -66.132).  The quantile of a
  # lower tail at most 1 / (1 + w) lies below the smallest positive double.
  for (a in c(u1, 3 * u1, 1e-310, 2^-1022 - u1)) {
    for (z in c(54, 55, 56, 57, 60, 100, 1000)) {
      info <- sprintf("alpha = %g, z = %g", a, z)
      l1pw <- log1p_exp(log(a) + logK0_ref(z))
      # x = u / sqrt(beta) stays a positive double at beta = 1e300.
      u <- c(1e-140, c(1e-3, 0.3, 0.6, 0.85, 0.95) * (z / 2))
      if (z >= 60) {
        phi <- c(650, 700, 750, 780, 800, 850)
        u <- c(u, 2 * phi / (z + sqrt((z - 2 * sqrt(phi)) * (z + 2 * sqrt(phi)))))
      }
      lI <- vapply(u, logI_ref, numeric(1), z = z)
      ref <- log1p_exp(log(a) + lI) - l1pw
      # At z = 54 the two smallest shapes have w < 1 and F above one half
      # everywhere, and nothing below one half is left to compare.
      keep <- ref < log(0.5)
      u <- u[keep]; lI <- lI[keep]; ref <- ref[keep]
      if (!length(u)) next
      # Both log F and the reference are differences of terms of the size of
      # log Psi, about z^2 / 4, and carry the rounding of such terms: up to
      # about 12 eps z^2 / 4 in absolute terms, the bound of ?pmhn, 6.7e-10
      # at z = 1000.  The rounding differs with the system's mathematical
      # library: at z = 1000 and beta = 1e-300 the upper tail below is 1 to
      # double precision on macOS and 5.8e-11 below 1 on Linux.
      term_err <- 3 * .Machine$double.eps * z^2
      for (beta in c(1, 1e-300, 1e300)) {
        g <- z * sqrt(beta)
        x <- u / sqrt(beta)
        expect_close(pmhn(x, a, beta, g, log.p = TRUE), ref, 1e-12,
                     abs_tol = term_err, info = info)
        # exp() passes the absolute error of log F to F as a relative one.
        lin <- ref > log(1e-300)
        expect_close(pmhn(x[lin], a, beta, g), exp(ref[lin]),
                     max(1e-11, term_err), info = info)
        # log(1 - F) and 1 - F where the forms serve the lower tail, up to
        # alpha I(u) = 2^64.  Past the hand-over 1 - F on the linear scale
        # comes from the upper tail of the paths for an ordinary shape, as at
        # 2^-1022, whose logarithm carries an absolute error of up to about
        # term_err.
        own <- log(a) + lI <= 64 * log(2)
        expect_close(pmhn(x[own], a, beta, g, lower.tail = FALSE, log.p = TRUE),
                     log1p(-exp(ref[own])), 1e-12,
                     abs_tol = max(2 * u1, term_err), info = info)
        expect_close(pmhn(x[own], a, beta, g, lower.tail = FALSE),
                     -expm1(ref[own]), 1e-12, abs_tol = max(2e-12, term_err),
                     info = info)
        # The quantile where F rises away from 1 / (1 + w), on both tails, and
        # 0 where F is 1 / (1 + w) in double precision: there the
        # probability is taken from pmhn() itself, since a reference that
        # differs from its 1 / (1 + w) in the last place would ask for a
        # quantile that is not 0.
        q_ok <- log(a) + lI > log(1e-3)
        expect_close(qmhn(ref[q_ok], a, beta, g, log.p = TRUE), x[q_ok], 1e-11,
                     info = info)
        qu <- q_ok & lin
        expect_close(qmhn(log1p(-exp(ref[qu])), a, beta, g, lower.tail = FALSE,
                          log.p = TRUE), x[qu], 1e-11, info = info)
        flat <- log(a) + lI < log(1e-17)
        expect_identical(qmhn(pmhn(x[flat], a, beta, g, log.p = TRUE), a, beta,
                              g, log.p = TRUE), rep(0, sum(flat)), info = info)
        expect_identical(qmhn(-l1pw - 1e-6, a, beta, g, log.p = TRUE), 0,
                         info = info)
        if (-l1pw > log(1e-300)) {
          expect_identical(qmhn(log1p(-exp(-l1pw - 1e-6)), a, beta, g,
                                lower.tail = FALSE, log.p = TRUE), 0,
                           info = info)
        }
      }
    }
  }
})

test_that("the lower tail of a subnormal shape agrees with the ordinary paths", {
  # A comparison that rests on neither the form nor its references: just
  # below 2^-1022 the forms serve the lower tail, and at 2^-1022 the paths
  # for an ordinary shape.  F (1 + w) - 1 is alpha I(u) at both shapes, and
  # where it is of order 1 to e^30 each gives I(u) to about 1e-12.  K0 enters
  # w at both, so that an error in its reference moves both alike.
  a0 <- 2^-1022 - u1
  a1 <- 2^-1022
  for (z in c(53.6, 54, 54.5, 54.9)) {
    lK <- logK0_ref(z)
    u <- c(seq(0.5, 0.95, by = 0.05), 0.97, 0.98) * (z / 2)
    lI <- vapply(u, logI_ref, numeric(1), z = z)
    u <- u[log(a0) + lI > 0 & log(a0) + lI < 30]
    expect_true(length(u) >= 2, info = sprintf("z = %g", z))
    I_of <- function(a) {
      lF <- pmhn(u, a, 1, z, log.p = TRUE)
      log(expm1(lF + log1p_exp(log(a) + lK))) - log(a)
    }
    expect_close(I_of(a0), I_of(a1), 1e-13, info = sprintf("z = %g", z))
  }
})

test_that("the moments of a subnormal shape follow alpha A_k", {
  # E X = alpha A_1 and Var X = alpha A_2 to a relative alpha A_1^2 / A_2;
  # the skewness is A_3 / (A_2^(3/2) sqrt(alpha)) and the excess kurtosis
  # A_4 / (A_2^2 alpha) - 3, which passes the largest double below a shape
  # that depends on the tilt (about 1.1e-308 at z = 0) and is then Inf.  In
  # 0.1.1 mhn_mean(3 * 2^-1074, 1, 12) was 1.510e-307 against 1.133e-307,
  # and mhn_skewness(2 * 2^-1074, 1, -30) -Inf against 6.34e161.
  for (a in c(3 * u1, 101 * u1, 1e-310, 1e-308)) {
    for (z in c(-30, -1, 1, 12, 21)) {
      info <- sprintf("alpha = %g, z = %g", a, z)
      A <- vapply(1:4, A_ref, numeric(1), z = z)
      if (a * A[1] > 1e-300) {
        expect_equal(mhn_mean(a, 1, z) / (a * A[1]), 1, tolerance = 1e-9,
                     info = info)
      }
      if (a * A[2] > 1e-300) {
        expect_equal(mhn_var(a, 1, z) / (a * A[2]), 1, tolerance = 1e-9,
                     info = info)
      }
      expect_equal(mhn_skewness(a, 1, z) * sqrt(a) / (A[3] / A[2]^1.5), 1,
                   tolerance = 1e-9, info = info)
      ku <- A[4] / A[2]^2 / a
      if (is.finite(ku)) {
        expect_equal(mhn_kurtosis(a, 1, z) / (ku - 3), 1, tolerance = 1e-9,
                     info = info)
      } else {
        expect_identical(mhn_kurtosis(a, 1, z), Inf, info = info)
      }
    }
  }
})

test_that("the moments of a subnormal shape hold where the peak counts", {
  # Where log w is between -37 and -14 the terms of the size of w that the
  # raw branch keeps change the variance in its eighth digit; the reference
  # forms the central moments from the raw ones, m_k = alpha A_k / (1 + w),
  # which cancel little there.  Above w = 2^-20 the law is the mixture of
  # the origin, with weight q = 1 / (1 + w), and the peak, whose moments are
  # integrated here about t = z / 2 by integrate(); from those the central
  # moments are formed without cancelling.  The raw form, which cancels by
  # up to w (557 at 2^-1074, z = 54.9), checks them more loosely.
  cases <- list(c(2^-1022 - u1, 52), c(2^-1022 - u1, 52.5),
                c(2^-1022 - u1, 52.7), c(2^-1022 - u1, 53),
                c(2^-1022 - u1, 53.2), c(2^-1022 - u1, 53.4),
                c(u1, 54.6), c(u1, 54.9))
  for (cs in cases) {
    a <- cs[1]; z <- cs[2]
    info <- sprintf("alpha = %g, z = %g", a, z)
    lw <- log(a) + logK0_ref(z)
    l1pw <- log1p_exp(lw)
    m <- vapply(1:4, function(k) exp(log(a) + logA_ref(k, z) - l1pw),
                numeric(1))
    raw <- c(m[1], m[2] - m[1]^2, m[3] - 3 * m[1] * m[2] + 2 * m[1]^3,
             m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4)
    got <- c(mhn_mean(a, 1, z), mhn_var(a, 1, z), mhn_skewness(a, 1, z),
             mhn_kurtosis(a, 1, z))
    as_moments <- function(mu) c(mu[1], mu[2], mu[3] / mu[2]^1.5,
                                 mu[4] / mu[2]^2 - 3)
    if (lw < -20 * log(2)) {
      expect_close(got, as_moments(raw), 1e-12, info = info)
    } else {
      h <- z / 2
      g <- function(t) exp(-(t - h)^2) / t
      K <- around(g, h)
      m1 <- around(function(t) t * g(t), h) / K
      v <- vapply(2:4, function(j) around(function(t) (t - m1)^j * g(t), h) / K,
                  numeric(1))
      q <- exp(-l1pw)
      p <- -expm1(-l1pw)
      c1 <- p * m1
      d <- q * m1
      mix <- c(c1, q * c1^2 + p * (v[1] + d^2),
               -q * c1^3 + p * (v[2] + 3 * v[1] * d + d^3),
               q * c1^4 + p * (v[3] + 4 * v[2] * d + 6 * v[1] * d^2 + d^4))
      expect_close(got, as_moments(mix), 1e-12, info = info)
      expect_close(got, as_moments(raw), 1e-8, info = info)
    }
  }
})

test_that("a subnormal shape joins the smallest normal one", {
  # Just below 2^-1022 the new forms serve, and at 2^-1022 the ordinary paths
  # do.  Per unit of alpha the two laws agree to double precision, at every
  # scale; the forms hand over at z = 54.97 for a0, so that z = 55 compares
  # the ordinary paths at both.
  a1 <- 2^-1022
  a0 <- a1 - u1
  for (beta in c(1e-300, 1, 1e300)) {
    for (z in c(-1e6, -30, -1, 0, 1, 12, 20, 30, 45, 53, 54.5, 55)) {
      g <- z * sqrt(beta)
      x <- c(1e-10, 0.3, 2, max(z / 2, 0) + 1) / sqrt(beta)
      info <- sprintf("beta = %g, z = %g", beta, z)
      expect_close(dmhn(x, a0, beta, g, log = TRUE) - log(a0),
                   dmhn(x, a1, beta, g, log = TRUE) - log(a1), 1e-12,
                   abs_tol = 1e-12, info = info)
      expect_close(
        pmhn(x, a0, beta, g, lower.tail = FALSE, log.p = TRUE) - log(a0),
        pmhn(x, a1, beta, g, lower.tail = FALSE, log.p = TRUE) - log(a1),
        1e-12, abs_tol = 1e-12, info = info)
      # The mean and the variance pass out of the double range at the
      # extreme scales; they are compared at beta = 1, where the moments
      # of the law at other scales follow by sqrt(beta) X ~ MHN(alpha, 1, z).
      if (beta == 1) {
        expect_close(mhn_mean(a0, beta, g) / mhn_mean(a1, beta, g), a0 / a1,
                     1e-11, info = info)
        expect_close(mhn_var(a0, beta, g) / mhn_var(a1, beta, g), a0 / a1,
                     1e-9, info = info)
      }
      expect_close(mhn_skewness(a0, beta, g) * sqrt(a0),
                   mhn_skewness(a1, beta, g) * sqrt(a1), 1e-9, info = info)
    }
  }
})

test_that("the upper quantile of a subnormal shape inverts alpha J(u)", {
  # For u far below the scale of the kernel, J(u) = K0(z) - log u, with
  # K0(0) = -0.5772.../2 and K0(1) = int_0^Inf exp(-t^2) (e^t - 1) / t dt
  # - 0.5772.../2, so the quantile of an upper tail Q is
  # exp(K0 - Q / alpha).  At 3 * 2^-1074, Q = 1e-320 and z = 1, qmhn() in
  # 0.1.1 gave 4.9e-324 against 2.6e-293.
  euler <- 0.57721566490153286
  D1 <- integrate(function(t) ifelse(t < 1e-12, 1, exp(-t^2) * expm1(t) / t),
                  0, 40, rel.tol = 1e-13)$value
  K0 <- c(`0` = -euler / 2, `1` = D1 - euler / 2)
  for (k in c(3, 5)) {
    for (z in c(0, 1)) {
      a <- k * u1
      ref <- exp(K0[[as.character(z)]] - 1e-320 / a)
      expect_equal(qmhn(1e-320, a, 1, z, lower.tail = FALSE) / ref, 1,
                   tolerance = 1e-9,
                   info = sprintf("alpha = %g * 2^-1074, z = %g", k, z))
    }
  }
  # Every probability of the lower tail short of 1 leaves the quantile below
  # the smallest positive double.
  expect_identical(qmhn(c(0, 0.5, 1 - 1e-16), 3 * u1, 1, 1), c(0, 0, 0))
  expect_identical(qmhn(1, 3 * u1, 1, 1), Inf)
})

test_that("both tails and the quantile at gamma = 0 follow Lemma 6a", {
  # At gamma = 0, beta X^2 ~ Gamma(alpha / 2, 1) (Sun et al. 2023, Lemma 6a),
  # and as alpha goes to 0, P(X > x) = Gamma(alpha / 2, beta x^2) /
  # Gamma(alpha / 2) is (alpha / 2) E1(beta x^2) to a relative O(alpha).
  euler <- 0.57721566490153286
  log_E1 <- function(y) {
    if (y < 1) {
      k <- 1:40
      return(log(-euler - log(y) - sum((-y)^k / (k * factorial(k)))))
    }
    -y + log(integrate(function(t) exp(-(t - y)) / t, y, Inf,
                       rel.tol = 1e-13)$value)
  }
  u <- c(1e-100, 1e-3, 0.5, 2, 5)
  lE <- vapply(u^2, log_E1, numeric(1))
  for (a in c(3 * u1, 1e-308, 2^-1022 - u1)) {
    lQ <- log(a) - log(2) + lE
    for (beta in c(1e-300, 1, 1e300)) {
      info <- sprintf("alpha = %g, beta = %g", a, beta)
      x <- u / sqrt(beta)
      expect_close(pmhn(x, a, beta, 0, lower.tail = FALSE, log.p = TRUE), lQ,
                   1e-13, info = info)
      # The quantile at u = 1e-100 is ill-conditioned: there log Q changes by
      # 2 / E1, under 0.005, per unit of log u.
      expect_close(qmhn(lQ, a, beta, 0, lower.tail = FALSE, log.p = TRUE), x,
                   1e-10, info = info)
      # Where Q is a normal double, the linear upper tail and the lower tail
      # on the log scale, log(1 - Q) = -Q, carry it in full.
      nq <- lQ > log(.Machine$double.xmin)
      if (any(nq)) {
        expect_close(qmhn(exp(lQ[nq]), a, beta, 0, lower.tail = FALSE), x[nq],
                     1e-8, info = info)
        expect_close(qmhn(-exp(lQ[nq]), a, beta, 0, log.p = TRUE), x[nq],
                     1e-8, info = info)
        expect_close(pmhn(x[nq], a, beta, 0, log.p = TRUE), -exp(lQ[nq]),
                     1e-12, info = info)
      }
    }
  }
})

test_that("the far upper tail of a subnormal shape follows its end point", {
  # Where the kernel falls steeply from u on, at the rate
  # lambda = u (2u - z) on the log axis, J(u) is
  # exp(-u^2 + z u) / lambda (1 + g'' / lambda^2) with g'' = u (z - 4u), the
  # next term below 1e-15 of it here (Watson's lemma).  Under z = -1e300 the
  # kernel falls by e^-4 within 4 / lambda of log u: 4e-8 at u = 1e-292,
  # which a first panel of the usual width, 1/8, would not resolve, and
  # 4e-100 at u = 1e-200, below the spacing of the doubles there.  The tails
  # are taken from log u, whose rounding, up to half a unit in its last
  # place (5.7e-14 for |log u| from 512 to 1024) or a unit for a libm that
  # rounds less closely, -u^2 + z u passes on relative to itself.
  for (cs in list(c(1e-292, -1e300), c(1e4, 0), c(1e-200, -1e300))) {
    u <- cs[1]; z <- cs[2]
    lambda <- u * (2 * u - z)
    ref <- log(3 * u1) + u * (z - u) - log(lambda) +
      log1p(u * (z - 4 * u) / lambda^2)
    expect_close(pmhn(u, 3 * u1, 1, z, lower.tail = FALSE, log.p = TRUE), ref,
                 2.5e-13, info = sprintf("u = %g, z = %g", u, z))
  }
})

test_that("a subnormal shape at a large tilt takes its law from the peak", {
  # Past the hand-over the density, the upper tail, the median and the
  # moments do not depend on alpha to double precision, and agree with those
  # at 2^-1022.  log D(z) is taken from its expansion about the peak there:
  # on the log axis the peak is 1 / (z / 2) wide, finer from z of about 1e15
  # on than the spacing of the doubles near log(z / 2), and an integral there
  # can come out -Inf, which would leave w = 0 and the log density at the mode
  # at z = 1e17 at 2.5e33 instead of -0.572.
  a1 <- 2^-1022
  for (a in c(u1, 3 * u1, 2^-1022 - u1)) {
    for (z in c(100, 1e10, 1e16, 5e16, 1e17, 2e17, 5e17, 1e20, 1e100, 1e154)) {
      info <- sprintf("alpha = %g, z = %g", a, z)
      x <- z / 2 + c(-3, -1, 0, 0.5, 2)
      expect_close(dmhn(x, a, 1, z, log = TRUE), dmhn(x, a1, 1, z, log = TRUE),
                   1e-12, abs_tol = 1e-12, info = info)
      expect_close(pmhn(x, a, 1, z, lower.tail = FALSE, log.p = TRUE),
                   pmhn(x, a1, 1, z, lower.tail = FALSE, log.p = TRUE), 1e-12,
                   abs_tol = 1e-15, info = info)
      expect_close(qmhn(c(0.1, 0.5, 0.9), a, 1, z),
                   qmhn(c(0.1, 0.5, 0.9), a1, 1, z), 1e-14, info = info)
      expect_close(mhn_mean(a, 1, z), mhn_mean(a1, 1, z), 1e-14, info = info)
      expect_close(mhn_var(a, 1, z), mhn_var(a1, 1, z), 1e-9, info = info)
    }
  }
  # At the mode the density of sqrt(beta) X is that of N(z / 2, 1 / 2) to a
  # relative 1 / z^2.
  expect_close(dmhn(5e16, u1, 1, 1e17, log = TRUE), -0.5 * log(pi), 1e-13)
})

test_that("the lower tail of a subnormal shape rises through its median", {
  # Below z / 2 the lower tail is tried first and kept at or below one half;
  # where it is above, the upper tail serves, with Q kept at or below one
  # half, so that F does not step back where the two forms, which differ by
  # their rounding (2e-14 near the median), take over from each other.
  for (cs in list(c(u1, 54.6691), c(3 * u1, 55), c(2^-1022 - u1, 53.54))) {
    a <- cs[1]; z <- cs[2]
    med <- qmhn(0.5, a, 1, z)
    x <- med * (1 + seq(-400, 400) * 2^-52)
    expect_true(all(diff(pmhn(x, a, 1, z)) >= 0),
                info = sprintf("alpha = %g, z = %g", a, z))
    expect_true(all(diff(pmhn(x, a, 1, z, log.p = TRUE)) >= 0),
                info = sprintf("alpha = %g, z = %g", a, z))
  }
})

test_that("just past the hand-over the lower tail above one half keeps its accuracy", {
  # Past the tilt where the forms hand over (w = 2^64) they still serve the
  # lower tail far below the peak, up to one half.  Above one half the
  # ordinary paths serve it, as log(1 - Q): from the forms log F would be the
  # difference of two terms of about 44, and just past the hand-over lost up
  # to 5 percent of itself near F = 1.  That happens only within about 1e-14
  # of the hand-over in z, so the tilts step through it one unit in the last
  # place at a time, from 20 below the reference's hand-over to 60 above.
  for (a in c(u1, 1e-310, 2^-1022 - u1)) {
    zh <- uniroot(function(z) log(a) + logK0_ref(z) - 64 * log(2), c(54, 57),
                  tol = 1e-13)$root
    for (z in zh * (1 + c((-20:60) * 2^-52, 1e-9, 1e-6))) {
      u <- z / 2 + c(0.5, 2, 3, 4, 4.95)
      lQ <- pmhn(u, a, 1, z, lower.tail = FALSE, log.p = TRUE)
      expect_close(pmhn(u, a, 1, z, log.p = TRUE), log1p(-exp(lQ)), 1e-10,
                   info = sprintf("alpha = %g, z = z_h + %g", a, z - zh))
    }
  }
})

test_that("vectorised parameters give the scalar results at a subnormal shape", {
  # The recycled parameters go through a cache of one state per triple
  # (mhn_param_slots.h), which a scalar call does not use.
  al <- c(3 * u1, 5 * u1, 0.5, u1, 2^-1022 - u1)
  zs <- c(1, 55, 57, 100, 1e17, 54)
  x <- c(1e-300, 1, 20, 27.4, 50, 5e16)
  n <- 60
  ia <- rep_len(seq_along(al), n); iz <- rep_len(seq_along(zs), n)
  ix <- rep_len(seq_along(x), n)
  one <- function(f, ...) vapply(seq_len(n), function(i)
    f(x[ix[i]], al[ia[i]], 1, zs[iz[i]], ...), numeric(1))
  expect_identical(dmhn(x[ix], al, 1, zs, log = TRUE), one(dmhn, log = TRUE))
  for (lt in c(TRUE, FALSE)) {
    for (lg in c(TRUE, FALSE)) {
      info <- sprintf("lower.tail = %s, log.p = %s", lt, lg)
      expect_identical(pmhn(x[ix], al, 1, zs, lower.tail = lt, log.p = lg),
                       one(pmhn, lower.tail = lt, log.p = lg), info = info)
    }
  }
  p <- c(1e-30, 0.3, 0.5, 0.9)
  ip <- rep_len(seq_along(p), n)
  for (lt in c(TRUE, FALSE)) {
    pp <- if (lt) log(p[ip]) else log1p(-p[ip])
    one_q <- vapply(seq_len(n), function(i)
      qmhn(pp[i], al[ia[i]], 1, zs[iz[i]], lower.tail = lt, log.p = TRUE),
      numeric(1))
    expect_identical(qmhn(pp, al, 1, zs, lower.tail = lt, log.p = TRUE), one_q,
                     info = sprintf("lower.tail = %s", lt))
  }
  # rmhn(n, alpha = v) draws what the element-wise loop draws, both where
  # each recycled triple keeps a setup of its own (n of at least 4 times the
  # period) and where one setup is rebuilt from triple to triple.  At
  # z = 1e17 the law is narrower than the spacing of the doubles at its
  # mode, and rmhn() refuses it, naming the shape asked for.
  zr <- zs[zs < 1e10]
  for (nn in c(n, 5)) {
    ir <- rep_len(seq_along(zr), nn)
    ja <- rep_len(seq_along(al), nn)
    set.seed(7)
    r_vec <- rmhn(nn, al, 1, zr)
    set.seed(7)
    r_one <- vapply(seq_len(nn), function(i) rmhn(1, al[ja[i]], 1, zr[ir[i]]),
                    numeric(1))
    expect_identical(r_vec, r_one, info = sprintf("n = %d", nn))
  }
  expect_error(rmhn(1, u1, 1, 1e17), "alpha = 4.94066e-324", fixed = TRUE)
  # At one beta and gamma a setup keeps log D and the sampler of the proxy
  # from one subnormal shape to the next, and gives them up for an ordinary
  # shape between them; the draws are those of the element-wise loop.  With
  # n = 7 one setup serves the shapes in turn; with n = 70 each residue of
  # the recycling keeps a setup of its own.
  for (z in c(52.5, 54.67, 57)) {
    av <- c(u1, 3 * u1, 0.3, 5 * u1, 2^-1022 - u1, 7 * u1, 1e-300)
    for (nn in c(7, 70)) {
      ia2 <- rep_len(seq_along(av), nn)
      set.seed(19)
      r_vec <- rmhn(nn, av, 1, z)
      set.seed(19)
      r_one <- vapply(seq_len(nn), function(i) rmhn(1, av[ia2[i]], 1, z),
                      numeric(1))
      expect_identical(r_vec, r_one, info = sprintf("z = %g, n = %d", z, nn))
    }
  }
  # A shape whose draws are all 0 between two that draw from the proxy, at
  # z = 53.2811, where that is so for 2^-1074 and not for 2^-1022 - 2^-1074
  # (q = 0.8), must not overwrite the proxy's sampler; and beta alone changing
  # between shapes must not reuse the setup of the other tilt.
  calls <- list(list(n = 20, a = rep(c(2^-1022 - u1, u1), 10), b = 1,
                     g = 53.2811),
                list(n = 6, a = c(u1, 3 * u1), b = c(1, 1.21), g = 54.67))
  for (cl in calls) {
    set.seed(41)
    r_vec <- rmhn(cl$n, cl$a, cl$b, cl$g)
    set.seed(41)
    r_one <- vapply(seq_len(cl$n), function(i)
      rmhn(1, rep_len(cl$a, cl$n)[i], rep_len(cl$b, cl$n)[i], cl$g),
      numeric(1))
    expect_identical(r_vec, r_one, info = sprintf("n = %d", cl$n))
  }
})

test_that("vectorised subnormal shapes at one tilt cost about one setup", {
  # log D and the sampler of the proxy depend on beta and gamma alone; built
  # for every shape they made 1e5 draws of distinct subnormal shapes seven
  # times as slow as those of tiny ordinary shapes, and with the reuse they
  # take a tenth of that time.  A timing, so not on CRAN.
  skip_on_cran()
  set.seed(43)
  a_sub <- runif(1e5) * 2^-1022
  a_ord <- 2^-1022 * (2 + runif(1e5))
  time_of <- function(a) min(vapply(1:3, function(i)
    system.time(rmhn(1e5, a, 1, 54.67))[["elapsed"]], numeric(1)))
  t_sub <- time_of(a_sub)
  t_ord <- time_of(a_ord)
  expect_true(t_sub < 0.5 * t_ord,
              info = sprintf("subnormal %.3f s, ordinary %.3f s", t_sub, t_ord))
  # A tilt that changes from draw to draw builds a setup for every draw,
  # subnormal or not; with log D from its expansion that costs about what an
  # ordinary shape's setup does, where the quadrature for log D made it seven
  # times as much.
  g <- 54.67 + runif(1e5)
  time_g <- function(a) min(vapply(1:3, function(i)
    system.time(rmhn(1e5, a, 1, g))[["elapsed"]], numeric(1)))
  t_sub_g <- time_g(a_sub)
  t_ord_g <- time_g(a_ord)
  expect_true(t_sub_g < 2 * t_ord_g,
              info = sprintf("subnormal %.3f s, ordinary %.3f s", t_sub_g,
                             t_ord_g))
})

test_that("rmhn draws a subnormal shape as its mixture of origin and peak", {
  # A draw is 0, below the smallest positive double, with probability
  # 1 / (1 + w), and otherwise comes from the peak.  In 0.1.1 rmhn() drew 0
  # for 73 percent of the draws at 2^-1074 and z = 54.669, where that
  # probability is 1/2.  Below a tilt of about 52 it is 1 in double
  # precision, and every draw is 0 without a random number.
  n <- 20000
  for (cs in list(c(u1, 54.6691), c(3 * u1, 55), c(2^-1022 - u1, 53.54),
                  c(u1, 57))) {
    a <- cs[1]; z <- cs[2]
    info <- sprintf("alpha = %g, z = %g", a, z)
    q0 <- pmhn(5e-324, a, 1, z)
    set.seed(11)
    x <- rmhn(n, a, 1, z)
    zeros <- mean(x == 0)
    expect_true(abs(zeros - q0) <= 5 * sqrt(q0 * (1 - q0) / n) + 1e-12,
                info = paste(info, "| zeros", zeros, "against", q0))
    # The draws that are not 0 follow the peak: P(X <= x | X > 0) is
    # (F(x) - q0) / (1 - q0).
    y <- sort(x[x > 0][1:2000])
    Fc <- (pmhn(y, a, 1, z) - q0) / (1 - q0)
    ks <- max(abs(Fc - seq_along(y) / length(y)),
              abs(Fc - (seq_along(y) - 1) / length(y)))
    expect_true(ks < 2.5 / sqrt(length(y)),
                info = paste(info, "| KS distance", ks))
  }
  # The draws that are not 0 have the mean of the peak: E(X | X > 0) is
  # mhn_mean() / (1 - q0), within 4 standard errors over 2e5 draws.  This
  # tells the proxy shape 2^-60 from a shape of 0.5, whose peak lies 0.007
  # further out (5 standard errors).
  for (cs in list(c(u1, 54.6691), c(u1, 57))) {
    a <- cs[1]; z <- cs[2]
    q0 <- pmhn(5e-324, a, 1, z)
    m1 <- mhn_mean(a, 1, z) / (1 - q0)
    sd1 <- sqrt((mhn_var(a, 1, z) + mhn_mean(a, 1, z)^2) / (1 - q0) - m1^2)
    set.seed(13)
    y <- rmhn(2e5, a, 1, z)
    y <- y[y > 0]
    expect_true(abs(mean(y) - m1) < 4 * sd1 / sqrt(length(y)),
                info = sprintf("alpha = %g, z = %g, t = %.2f", a, z,
                               (mean(y) - m1) / (sd1 / sqrt(length(y)))))
  }
  # At an extreme scale the test t >= 1 is made on t = sqrt(beta) x, and the
  # draws are those at beta = 1 scaled by 1 / sqrt(beta).
  set.seed(17)
  x1 <- rmhn(2000, u1, 1, 54.6691)
  for (beta in c(1e-300, 1e300)) {
    set.seed(17)
    xb <- rmhn(2000, u1, beta, 54.6691 * sqrt(beta))
    expect_false(anyNA(xb), info = sprintf("beta = %g", beta))
    expect_identical(xb == 0, x1 == 0, info = sprintf("beta = %g", beta))
    expect_close(sqrt(beta) * xb[xb > 0], x1[x1 > 0], 1e-12,
                 info = sprintf("beta = %g", beta))
  }
  set.seed(3)
  before <- .Random.seed
  expect_identical(rmhn(5, 3 * u1, 1, 30), rep(0, 5))
  expect_identical(.Random.seed, before)
})

test_that("a tiny ordinary shape keeps its lower tail finite, at any beta", {
  # Next to the origin the lower tail of a shape near 1e-300 is about
  # 1 / alpha, e^690, while the kernel at x = u / sqrt(beta) carries the
  # factor 1 / x: e^-325 or less at beta = 1e-300 (x from 1e141 here) and
  # e^-152 at beta = 1e-150.  Measured against that kernel, the part next to
  # the origin passes the largest double; the lower-tail quadrature measures
  # against that part instead wherever it is the larger by e^600.  0.1.1
  # returned -Inf at all of these points, at beta = 1 as well.  For alpha
  # below 1e-19, t^alpha = 1 to double precision, and the lower tail is
  # (1 + alpha I(u)) / (1 + w) as for a subnormal shape.  |log F| is 2.5e5
  # to 2.5e7 here, and its error is measured at up to 4.5e-16 of it.  At
  # beta = 1e-150 the part next to the origin exceeds the kernel by e^714 to
  # e^718 at alpha = 1e-246 and u = 1e-11 to 1e-9, and by e^709.9 to
  # e^711.5 at alpha = 1e-240 and u = 2e-7 to 1e-6, just past the largest
  # double (e^709.78), where the quadrature serves: from a tilt of about 900
  # at such a shape, below which the series does.
  for (a in c(2^-1022, 1e-300, 1e-250, 1e-246, 1e-240, 1e-200)) {
    for (z in c(999, 1001, 1e4)) {
      l1pw <- log1p_exp(log(a) + logK0_ref(z))
      u <- c(1e-10, 1e-6, 1e-4, 1e-2, 0.1) * (z / 2) / 50
      if (a == 1e-246) u <- c(1e-11, 1e-10, 1e-9)
      if (a == 1e-240) u <- c(2e-7, 5e-7, 1e-6)
      lI <- vapply(u, logI_ref, numeric(1), z = z)
      ref <- log1p_exp(log(a) + lI) - l1pw
      for (beta in c(1e-300, 1e-150, 1)) {
        expect_close(pmhn(u / sqrt(beta), a, beta, z * sqrt(beta), log.p = TRUE),
                     ref, 1e-14, info = sprintf("alpha = %g, z = %g, beta = %g",
                                                a, z, beta))
      }
    }
  }

  # A log probability below log F at the smallest positive double has its
  # quantile below that double, and 0 is returned; above it the quantile
  # gives the probability back.
  a <- 1e-300; b <- 1e-300; g <- 999 * sqrt(b)
  expect_true(pmhn(5e-324, a, b, g, log.p = TRUE) > -1e6)
  expect_identical(qmhn(c(-1e12, -1e8, -1e6), a, b, g, log.p = TRUE),
                   c(0, 0, 0))
  q <- qmhn(-1e5, a, b, g, log.p = TRUE)
  expect_true(q > 0)
  expect_close(pmhn(q, a, b, g, log.p = TRUE), -1e5, 1e-12)

  # At tilts of 50 to 55 pmhn() takes the series, and the quadrature is
  # reached through the internal hook.  It is checked there on its own, at
  # points where the part next to the origin mostly exceeds the kernel at x
  # by more than e^600: the terms past that part must be measured against
  # it too.  The error measured is 1.1e-13 of F, the rounding of log terms
  # of the size of log(1 / alpha).  alpha I(u) is below the rounding of 1
  # at these points, so this checks the scale, not how F depends on u.
  for (a in c(1e-300, 1e-250)) {
    for (z in c(50, 54, 55)) {
      l1pw <- log1p_exp(log(a) + logK0_ref(z))
      u <- c(0.01, 0.1, 0.5, 1)
      ref <- exp(log1p_exp(log(a) + vapply(u, logI_ref, numeric(1), z = z)) -
                   l1pw)
      for (beta in c(1, 1e-300)) {
        got <- vapply(u / sqrt(beta), function(q) mhn:::.pmhn_force(
          q, a, beta, z * sqrt(beta), "integrate"), numeric(1))
        expect_close(got, ref, 1e-12, info = sprintf(
          "quadrature, alpha = %g, z = %g, beta = %g", a, z, beta))
      }
    }
  }
})

test_that("log D of the subnormal-shape forms is rounded once", {
  # log D = m^2 + log(sqrt(pi) / m) + log1p(sum) with m = z / 2, summed with
  # m^2 split exactly by fma, is within 0.65 units in the last place.  The
  # references are D's Taylor series sum_n z^n Gamma(n/2) / (2 n!) at 256
  # bits (Rmpfr), as the nearest double hi and the rest lo, at the ten tilts
  # of a scan of 308 from 30 to 99 where the sum added left to right was
  # farthest off, 1.05 to 1.19 units as Apple clang compiles it (seven of
  # them are past 0.75 also where m^2 + 0.5 log(pi) is not fused into one
  # rounding): z = 30, 30.22, 31.35, 41.91, 42.14,
  # 46.41, 50.00, 65.74, 73.15 and 97.43; and the four besides 41.91 where
  # the sum without lo, rounded at m^2, is farthest off, 0.88 to 0.99 units:
  # z = 53.37, 40.79, 60.79 and 90.91.  They are written in hexadecimal,
  # which R reads exactly; a decimal of 17 digits it can read as the double
  # next to the one meant.
  ref <- rbind(
    c(0x1.ep+4, 0x1.bdbbac5e80c0ap+7, -0x1.8914e653a86dp-48),
    c(0x1.e398996ee7d15p+4, 0x1.c47c6ca8ce2c2p+7, -0x1.18946ac1e412ep-49),
    c(0x1.f59398996ee7dp+4, 0x1.e702948543d3p+7, -0x1.3b89b1853b5edp-49),
    c(0x1.4f4be1faff2a8p+5, 0x1.b4afa4d179cecp+8, 0x1.c366cd2341992p-49),
    c(0x1.51182eb273133p+5, 0x1.b96742d239547p+8, 0x1.2b71c794c3997p-47),
    c(0x1.7341e0500d579p+5, 0x1.0beae9aaf7a5p+9, -0x1.37f8a258244f5p-46),
    c(0x1.9006abc74be2p+5, 0x1.3737c3c04ecb9p+9, -0x1.36e2f3d96f084p-47),
    c(0x1.06f1d2f87ebfcp+6, 0x1.0d590454d6c3ep+10, 0x1.a14cc25840359p-47),
    c(0x1.249cc4cb773e8p+6, 0x1.4db44c01f78b9p+10, -0x1.4277e58acbf68p-45),
    c(0x1.85b4f37dea51bp+6, 0x1.2835aec1d0028p+11, 0x1.cde6f924dbf7fp-45),
    c(0x1.aaff2a871683cp+5, 0x1.62c015f96801bp+9, 0x1.c5082f01605c7p-50),
    c(0x1.464e6265bb9f4p+5, 0x1.9d7aa8e17f083p+8, 0x1.2c4b2f376365dp-48),
    c(0x1.e6550e2d07814p+5, 0x1.cc87d429bdbcep+9, -0x1.58b873f68b3ap-47),
    c(0x1.6ba29b19d9a44p+6, 0x1.01db75e1fd655p+11, 0x1.ecda2c14fe41ep-45))
  got <- vapply(ref[, 1], mhn:::.vanishing_log_D, numeric(1))
  units <- ((got - ref[, 2]) - ref[, 3]) / 2^(floor(log2(ref[, 2])) - 52)
  expect_true(all(abs(units) <= 0.75),
              info = paste("units off:", paste(format(units, digits = 3),
                                               collapse = " ")))
})
