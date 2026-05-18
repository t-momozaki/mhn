# Tests for .check_mhn_params()

test_that(".check_mhn_params rejects invalid alpha", {
  expect_error(mhn:::.check_mhn_params(0, 1, 0), "alpha must be positive")
  expect_error(mhn:::.check_mhn_params(-1, 1, 0), "alpha must be positive")
  expect_error(mhn:::.check_mhn_params(NA, 1, 0), "alpha must be positive")
  expect_error(mhn:::.check_mhn_params("a", 1, 0), "alpha must be positive")
  expect_error(mhn:::.check_mhn_params(c(1, 2), 1, 0), "alpha must be positive")
})

test_that(".check_mhn_params rejects invalid beta", {
  expect_error(mhn:::.check_mhn_params(1, 0, 0), "beta must be positive")
  expect_error(mhn:::.check_mhn_params(1, -1, 0), "beta must be positive")
  expect_error(mhn:::.check_mhn_params(1, NA, 0), "beta must be positive")
})

test_that(".check_mhn_params rejects invalid gamma", {
  expect_error(mhn:::.check_mhn_params(1, 1, NA), "gamma must be a finite")
  expect_error(mhn:::.check_mhn_params(1, 1, "a"), "gamma must be a finite")
})

test_that(".check_mhn_params accepts valid parameters", {
  expect_silent(mhn:::.check_mhn_params(1, 1, 0))
  expect_silent(mhn:::.check_mhn_params(0.5, 2, -3))
  expect_silent(mhn:::.check_mhn_params(10, 0.001, 100))
})

test_that(".convert_to_gw returns correct values", {
  result <- mhn:::.convert_to_gw(2, 3, -1)
  expect_equal(result$lambda_gw, 2)
  expect_equal(result$alpha_gw, 3)
  expect_equal(result$beta_gw, 1)  # -(-1) = 1

  result2 <- mhn:::.convert_to_gw(1, 1, 5)
  expect_equal(result2$lambda_gw, 1)
  expect_equal(result2$alpha_gw, 1)
  expect_equal(result2$beta_gw, -5)  # -(5) = -5
})
