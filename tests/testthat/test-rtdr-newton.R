# Unit tests for the RTDR contact-point search, exercised in isolation through
# .rtdr_contact_point_newton_test_cpp and its three fixture targets.
#
# The search solves log f(m) - log f(t) = delta for a contact point t, accepting
# early once the increment falls inside the band that Gao & Wang
# (2025, Eq. (7)) give for the transform in use -- [0.46, 2.49] under T_0 and
# [0.93, 1.99] under T_{-1/2}.  A point returned on the wrong side of the mode
# gives a tangent with the wrong sign of slope, and one far outside the band a
# hat much wider than the target, so the routine is tested here on targets
# whose contact points can be written down.

target <- function(tag) {
  switch(tag,
    log_concave_normal     = list(log_f = function(t) -0.5 * t^2,
                                  dlog_f = function(t) -t,
                                  mode = 0, left_init = -1,
                                  band = c(0.46, 2.49), delta = 1),
    log_concave_gamma_like = list(log_f = function(t) 2 * log(t) - t^2,
                                  dlog_f = function(t) 2 / t - 2 * t,
                                  # supported on t > 0, so the left search
                                  # starts at half the mode, as region A does
                                  mode = 1, left_init = 0.5,
                                  band = c(0.46, 2.49), delta = 1),
    tneghalf               = list(log_f = function(y) 0.7 * y - exp(2 * y),
                                  dlog_f = function(y) 0.7 - 2 * exp(2 * y),
                                  mode = log(0.35) / 2,
                                  left_init = log(0.35) / 2 - 1,
                                  band = c(0.93, 1.99), delta = log(4)))
}

test_that("the contact point lands inside the acceptance band", {
  for (tag in c("log_concave_normal", "log_concave_gamma_like", "tneghalf")) {
    tg <- target(tag)
    log_f_mode <- tg$log_f(tg$mode)
    for (side in c("left", "right")) {
      t_init <- if (side == "right") tg$mode + 1 else tg$left_init
      r <- mhn:::.rtdr_contact_point_newton_test_cpp(
        t_init, log_f_mode, tg$delta, tag, 100L)
      info <- sprintf("%s / %s", tag, side)
      expect_true(is.finite(r$t), info = info)
      expect_gte(r$increment, tg$band[1] - 1e-9, label = info)
      expect_lte(r$increment, tg$band[2] + 1e-9, label = info)
      # The C++ fixture must be the same function as the R one beside it.
      expect_equal(r$log_f_t, tg$log_f(r$t), tolerance = 1e-12, info = info)
      expect_equal(r$dlog_f_t, tg$dlog_f(r$t), tolerance = 1e-12, info = info)
      # Not `r$increment == log_f_mode - r$log_f_t`: the shim computes the
      # increment that way, so that comparison holds whatever the search does.
      # What can fail is convergence -- restarting from the point returned must
      # return it unchanged, since its increment is already inside the band.
      r2 <- mhn:::.rtdr_contact_point_newton_test_cpp(r$t, log_f_mode, tg$delta,
                                                      tag, 100L)
      expect_equal(r2$t, r$t, tolerance = 1e-12, info = info)
    }
  }
})

test_that("the contact point falls on the side it was started from", {
  # The envelope needs a positive slope on the left and a negative one on the
  # right; a point returned on the wrong side has the wrong sign and the region
  # setups reject it.  A left start with an absolute floor, such as
  # max(m/2, 1e-6), would cross the mode whenever the mode is below 2e-6, so
  # the starts here lie on the side they are meant for.
  for (tag in c("log_concave_normal", "log_concave_gamma_like", "tneghalf")) {
    tg <- target(tag)
    log_f_mode <- tg$log_f(tg$mode)
    right <- mhn:::.rtdr_contact_point_newton_test_cpp(
      tg$mode + 1, log_f_mode, tg$delta, tag, 100L)
    left <- mhn:::.rtdr_contact_point_newton_test_cpp(
      tg$left_init, log_f_mode, tg$delta, tag, 100L)
    expect_gt(right$t, tg$mode)
    expect_lt(right$dlog_f_t, 0)
    expect_lt(left$t, tg$mode)
    expect_gt(left$dlog_f_t, 0)
  }
})

test_that("a start outside the band takes the steps of Eq. (8) until it is inside", {
  # Gao & Wang (2025), Eq. (8): t <- t + (increment - delta) / (log f)'(t),
  # stopping at the first iterate whose increment lies in the band, which is
  # not the root of increment = delta.  For the normal target log f(t) = -t^2/2
  # with the mode at 0 the increment is t^2/2, so from 0.2 the steps visit 5.1
  # and 2.746 and stop at 1.737, inside [sqrt(0.92), sqrt(4.98)], and not at
  # sqrt(2).  Every start below lies outside the band, so at least one step is
  # taken, and none of the steps is cut short by the search's tolerances.
  eq8 <- function(tg, t, log_f_mode) {
    for (i in 1:100) {
      inc <- log_f_mode - tg$log_f(t)
      if (inc >= tg$band[1] && inc <= tg$band[2]) return(t)
      t <- t + (inc - tg$delta) / tg$dlog_f(t)
    }
    stop("Eq. (8) did not reach the band")
  }
  starts <- list(log_concave_normal = c(0.2, 10, -0.2, -5),
                 log_concave_gamma_like = c(0.05, 4),
                 tneghalf = log(0.35) / 2 + c(-4, 1.5, 3))
  for (tag in names(starts)) {
    tg <- target(tag)
    log_f_mode <- tg$log_f(tg$mode)
    for (t0 in starts[[tag]]) {
      info <- sprintf("%s from %g", tag, t0)
      inc0 <- log_f_mode - tg$log_f(t0)
      expect_true(inc0 < tg$band[1] || inc0 > tg$band[2], info = info)
      r <- mhn:::.rtdr_contact_point_newton_test_cpp(t0, log_f_mode, tg$delta,
                                                     tag, 100L)
      expect_equal(r$t, eq8(tg, t0, log_f_mode), tolerance = 1e-12, info = info)
      expect_gte(r$increment, tg$band[1], label = info)
      expect_lte(r$increment, tg$band[2], label = info)
    }
  }

  # Starting exactly at the analytic solution must return it unchanged, since
  # the increment is then delta, which is inside the band.
  r2 <- mhn:::.rtdr_contact_point_newton_test_cpp(sqrt(2), 0, 1,
                                                  "log_concave_normal", 100L)
  expect_equal(r2$t, sqrt(2), tolerance = 1e-12)
})

test_that("a start far from the mode still converges", {
  # The start Gao & Wang (2025, Section 2) recommend for their Eq. (8)
  # iteration is m +- sqrt(-2 delta / L''(m)), which for a flat mode lands
  # hundreds of units away -- 447 above the mode at alpha = 1e-5, where the
  # log-axis target's exp(2y) overflows.
  for (tag in c("log_concave_normal", "log_concave_gamma_like")) {
    tg <- target(tag)
    log_f_mode <- tg$log_f(tg$mode)
    for (t_init in c(tg$mode + 50, tg$mode + 500)) {
      r <- mhn:::.rtdr_contact_point_newton_test_cpp(
        t_init, log_f_mode, tg$delta, tag, 200L)
      expect_true(is.finite(r$t), info = sprintf("%s from %g", tag, t_init))
      expect_true(is.finite(r$dlog_f_t))
      expect_gte(r$increment, tg$band[1] - 1e-9,
                 label = sprintf("%s from %g", tag, t_init))
      expect_lte(r$increment, tg$band[2] + 1e-9,
                 label = sprintf("%s from %g", tag, t_init))
    }
  }
})

test_that("an unknown concavity tag is rejected", {
  expect_error(
    mhn:::.rtdr_contact_point_newton_test_cpp(1, 0, 1, "not_a_transform", 30L),
    "concavity")
})

test_that("a contact point short of the band still gives a valid hat", {
  # The search returns its last iterate when it does not reach the band, and
  # no caller tests the band.  Region A's left contact just above alpha = 1 is
  # such a case: log f falls towards x = 0 only like (alpha - 1) log(m / x), so
  # each Newton step overshoots below zero and is cut back to a halving of t,
  # and 30 of them gain a drop of only about 21 (alpha - 1).  At these two
  # points the drop reached is about 0.25 and 0.31, below the 0.46 that
  # Gao & Wang (2025, Eq. (7)) ask for.  At alpha = 1.0001 the band lies below
  # the smallest positive double, so that point stays short of it whatever the
  # search does; at alpha = 1.003 it lies at t of about 2e-31, and the drop is
  # checked to be short of it.  log f is concave on x > 0, so the tangent
  # there still dominates it and the draws stay exact; what is lost is the
  # Eq. (7) guarantee that the rejection constant stays below e.  The
  # constant is the envelope's area over the target's, both measured from the
  # peak; it is about 1.1 at both points.
  for (alpha in c(1.0001, 1.003)) {
    env <- mhn:::.dump_rtdr_envelope_cpp(alpha, 1, 1)
    info <- sprintf("alpha = %g", alpha)
    expect_identical(env$region, 0L)
    expect_false(env$simplified)
    log_f <- function(x) (alpha - 1) * log(x) - x^2 + x - env$log_dens_mode
    expect_lt(-log_f(env$t_l), 0.46, label = info)
    x <- seq(1e-6, 4, length.out = 4001)
    expect_true(all(vapply(x, function(xi) {
      p <- env$pieces[[if (xi <= env$p_l) 1L else if (xi <= env$p_r) 2L else 3L]]
      anchor <- if (p$type == PIECE_EXP_BOUNDED_RIGHT) p$b else p$a
      h <- if (p$type == PIECE_PLATEAU) p$base_log_dens
           else p$base_log_dens + p$slope * (xi - anchor)
      h - env$log_dens_mode >= log_f(xi) - 1e-12
    }, logical(1))), info = info)
    target <- stats::integrate(function(x) exp(log_f(x)), 0, Inf,
                               rel.tol = 1e-10)$value
    ratio <- sum(exp(env$piece_log_area - env$log_dens_mode)) / target
    expect_gt(ratio, 1, label = info)
    expect_lt(ratio, exp(1), label = info)
  }
})
