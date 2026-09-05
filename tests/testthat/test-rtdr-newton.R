# Unit tests for the RTDR contact-point search, exercised in isolation through
# .rtdr_contact_point_newton_test_cpp and its three fixture targets.
#
# The search solves log f(m) - log f(t) = delta for a contact point t, accepting
# early once the increment falls inside the band that Gao & Wang (2025, Eq. 7)
# give for the transform in use -- [0.46, 2.49] under T_0 and [0.93, 1.99] under
# T_{-1/2}.  Two of this package's sampling defects came from that routine
# returning a point it should have rejected, so it is worth testing on targets
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
  # setups reject it.  Region A used to start its left search at
  # max(m/2, 1e-6), which crossed the mode whenever the mode was below 2e-6.
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

test_that("the contact point matches the value solved for analytically", {
  # For the normal target log f(t) = -t^2/2 with the mode at 0, the increment is
  # exactly t^2/2, so the band [0.46, 2.49] is the interval
  # [sqrt(0.92), sqrt(4.98)] = [0.959, 2.232] and delta = 1 puts the exact
  # contact point at sqrt(2).
  r <- mhn:::.rtdr_contact_point_newton_test_cpp(1, 0, 1, "log_concave_normal", 100L)
  expect_gte(abs(r$t), sqrt(0.92) - 1e-9)
  expect_lte(abs(r$t), sqrt(4.98) + 1e-9)
  expect_equal(r$increment, r$t^2 / 2, tolerance = 1e-12)

  # Starting exactly at the analytic solution must return it unchanged, since
  # the increment is then delta, which is inside the band.
  r2 <- mhn:::.rtdr_contact_point_newton_test_cpp(sqrt(2), 0, 1,
                                                  "log_concave_normal", 100L)
  expect_equal(r2$t, sqrt(2), tolerance = 1e-12)
})

test_that("a start far from the mode still converges", {
  # The recommended start of Gao & Wang (2025, Eq. 8) is m +- sqrt(-2 delta /
  # L''(m)), which for a flat mode lands hundreds of units away -- 447 above the
  # mode at alpha = 1e-5, where the log-axis target's exp(2y) overflows.
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
