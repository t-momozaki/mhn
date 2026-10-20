# Tests of how the functions take their arguments: the type of each argument,
# parameters outside the parameter space, NA and NaN, zero-length arguments,
# n and the logical flags.  dmhn, pmhn, qmhn and rmhn follow base R's
# distribution functions (?dgamma: "Invalid arguments will result in return
# value NaN, with a warning"); the moment functions and mhn_mode, which take
# one triple, stop instead.

dpq <- list(
  dmhn = function(x, a = 2, b = 1, g = 1, ...) dmhn(x, a, b, g, ...),
  pmhn = function(x, a = 2, b = 1, g = 1, ...) pmhn(x, a, b, g, ...),
  qmhn = function(x, a = 2, b = 1, g = 1, ...) qmhn(x, a, b, g, ...)
)
first_arg <- c(dmhn = "x", pmhn = "q", qmhn = "p")
scalar_fns <- list(mhn_mean = mhn_mean, mhn_var = mhn_var,
                   mhn_skewness = mhn_skewness, mhn_kurtosis = mhn_kurtosis,
                   mhn_mode = mhn_mode)
non_numeric <- list(factor("3"), "3", list(3), 3 + 0i, NULL)

test_that("a non-numeric argument of dmhn, pmhn, qmhn or rmhn is an error, as in dnorm", {
  # A factor would otherwise be used through its level codes: dmhn(factor("3"))
  # was the density at 1.
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    for (v in non_numeric) {
      lab <- paste(nm, deparse(v))
      expect_error(f(v), sprintf("'%s' must be numeric", first_arg[[nm]]),
                   fixed = TRUE, label = lab)
      expect_error(f(0.5, a = v), "'alpha' must be numeric", fixed = TRUE,
                   label = lab)
      expect_error(f(0.5, b = v), "'beta' must be numeric", fixed = TRUE,
                   label = lab)
      expect_error(f(0.5, g = v), "'gamma' must be numeric", fixed = TRUE,
                   label = lab)
    }
  }
  for (v in non_numeric) {
    expect_error(rmhn(2, alpha = v), "'alpha' must be numeric", fixed = TRUE)
    expect_error(rmhn(2, beta = v), "'beta' must be numeric", fixed = TRUE)
    expect_error(rmhn(2, gamma = v), "'gamma' must be numeric", fixed = TRUE)
  }
})

test_that("integer, logical and classed numeric arguments are taken as their numbers", {
  # As in dnorm, which takes a Date as its number of days; is.numeric() is
  # FALSE for a Date, so a check made with it would refuse one.
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    expect_identical(f(1L, 2L, 1L, 1L), f(1, 2, 1, 1), label = nm)
    expect_identical(f(TRUE, 2, TRUE, FALSE), f(1, 2, 1, 0), label = nm)
  }
  day2 <- as.Date("1970-01-03")
  expect_identical(dmhn(day2, 2, 1, 1), dmhn(2, 2, 1, 1))
  expect_identical(pmhn(day2, 2, 1, 1), pmhn(2, 2, 1, 1))
  expect_identical(dmhn(1, as.difftime(2, units = "secs"), 1, 1),
                   dmhn(1, 2, 1, 1))
})

test_that("a parameter outside the parameter space gives NaN with one warning", {
  # Sun et al. (2023, Eq. (1)): alpha > 0, beta > 0, gamma real; the
  # shape and the rate must also be finite.  Each element stands alone, as
  # in dgamma(1, shape = c(1, -1)).  expect_warning() catches only the first
  # warning, so the warnings are collected and counted instead.
  with_warnings <- function(expr) {
    msgs <- character(0)
    value <- withCallingHandlers(expr, warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = msgs)
  }
  expect_nan_once <- function(r, expected, label) {
    expect_identical(r$warnings, "NaNs produced", label = label)
    expect_identical(r$value, expected, label = label)
  }
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    ok <- f(0.5)
    expect_nan_once(with_warnings(f(0.5, a = c(2, -1, 0, Inf, -Inf))),
                    c(ok, NaN, NaN, NaN, NaN), label = nm)
    expect_nan_once(with_warnings(f(0.5, b = c(1, -1, 0, Inf, -Inf))),
                    c(ok, NaN, NaN, NaN, NaN), label = nm)
    expect_nan_once(with_warnings(f(0.5, g = c(1, Inf, -Inf))),
                    c(ok, NaN, NaN), label = nm)
    # The scalar path.
    expect_nan_once(with_warnings(f(c(0.5, 0.7), a = -1)), c(NaN, NaN),
                    label = nm)
  }
  expect_nan_once(with_warnings(dmhn(0.5, -1, 1, 1, log = TRUE)), NaN,
                  label = "dmhn, log")
  expect_nan_once(with_warnings(pmhn(0.5, -1, 1, 1, lower.tail = FALSE,
                                     log.p = TRUE)), NaN, label = "pmhn, log")
  expect_nan_once(with_warnings(qmhn(log(0.3), -1, 1, 1, lower.tail = FALSE,
                                     log.p = TRUE)), NaN, label = "qmhn, log")
})

test_that("the warning names the call", {
  w <- tryCatch(dmhn(1, -1), warning = function(w) w)
  expect_identical(conditionCall(w), quote(dmhn(1, -1)))
})

test_that("an invalid element leaves the other elements as they were", {
  x <- c(0.3, 0.8, 1.5, 2.2)
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    arg <- if (nm == "qmhn") x / 3 else x
    expect_warning(out <- f(arg, a = c(2, -1, 3, NA)), "NaNs produced")
    expect_identical(out[c(1, 3)], f(arg[c(1, 3)], a = c(2, 3)), label = nm)
  }
})

test_that("an NA or NaN parameter gives NA or NaN without a warning, as in dgamma", {
  # NA wins over NaN, and over a parameter outside the space, wherever it
  # stands; no warning is given, since an input was already NA or NaN.
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    ok <- f(0.5)
    expect_silent(out <- f(0.5, a = c(2, NA, NaN)))
    expect_identical(out, c(ok, NA, NaN), label = nm)
    expect_silent(out <- f(0.5, b = c(1, NA, NaN), g = c(1, 1, NA)))
    expect_identical(out, c(ok, NA, NA), label = nm)
    expect_silent(out <- f(c(NaN, NaN, NA), a = c(-1, NA, NaN)))
    expect_identical(out, c(NaN, NA, NA), label = nm)
    expect_silent(out <- f(0.5, a = NA_real_))
    expect_identical(out, NA_real_, label = nm)
  }
  # On the probability scale too: the density is the exponential of its log,
  # and an NA must not become NaN on the way.
  expect_identical(dmhn(c(NA, NaN, 1), 2, 1, 1),
                   c(NA, NaN, dmhn(1, 2, 1, 1)))
  expect_identical(dmhn(1, c(NA, NaN), 1, 1, log = TRUE), c(NA, NaN))
})

test_that("a zero-length argument gives numeric(0), as in dnorm", {
  for (nm in names(dpq)) {
    f <- dpq[[nm]]
    expect_identical(f(0.5, a = numeric(0)), numeric(0), label = nm)
    expect_identical(f(0.5, b = numeric(0)), numeric(0), label = nm)
    expect_identical(f(0.5, g = numeric(0)), numeric(0), label = nm)
    expect_identical(f(numeric(0), a = c(2, 3)), numeric(0), label = nm)
    expect_identical(f(numeric(0), a = -1), numeric(0), label = nm)
  }
})

test_that("rmhn draws NaN with one warning for an NA or invalid parameter, as rnorm does", {
  set.seed(1)
  expect_warning(x <- rmhn(5, alpha = c(2, NA, -1, 0, Inf), gamma = 1),
                 "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_warning(x <- rmhn(5, beta = c(1, NaN, -0.5, 0, Inf), gamma = 1),
                 "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_warning(x <- rmhn(5, gamma = c(1, NA, NaN, Inf, -Inf)),
                 "NAs produced")
  expect_identical(is.nan(x), c(FALSE, TRUE, TRUE, TRUE, TRUE))
  # method = "sun" refuses alpha < 1 with gamma > 0, but a shape outside the
  # parameter space is an invalid argument, not that case.
  expect_warning(x <- rmhn(2, alpha = -1, gamma = 1, method = "sun"),
                 "NAs produced")
  expect_identical(x, c(NaN, NaN))
})

test_that("an invalid element of rmhn draws no random numbers", {
  # The other elements get the draws they get without it.
  set.seed(3)
  expect_warning(x <- rmhn(4, alpha = c(2, -1, 3, NA), gamma = 1),
                 "NAs produced")
  set.seed(3)
  expect_identical(x[c(1, 3)], rmhn(2, alpha = c(2, 3), gamma = 1))
})

test_that("a zero-length parameter of rmhn draws NA with a warning, as rnorm does", {
  expect_warning(x <- rmhn(2, alpha = numeric(0)), "NAs produced")
  expect_identical(x, c(NA_real_, NA_real_))
  expect_warning(x <- rmhn(2, gamma = numeric(0)), "NAs produced")
  expect_identical(x, c(NA_real_, NA_real_))
  expect_identical(rmhn(0, alpha = numeric(0)), numeric(0))
})

test_that("rmhn reads n as rnorm does, up to .Machine$integer.max", {
  # A vector n asks for length(n) variates, a zero-length one included.
  expect_length(rmhn(c(10, 20, 30), 2, 1, 1), 3)
  expect_identical(rmhn(numeric(0), 2, 1, 1), numeric(0))
  expect_identical(rmhn(0.5, 2, 1, 1), numeric(0))
  # NULL, NA, an infinite n and one beyond the integer range are errors, with
  # no warning on the way.
  for (n in list(NULL, NA, NaN, Inf, 2^31)) {
    warned <- FALSE
    expect_error(
      withCallingHandlers(rmhn(n, 2, 1, 1),
                          warning = function(w) warned <<- TRUE),
      "n must be a finite, non-negative number no larger than .Machine$integer.max",
      fixed = TRUE, label = deparse(n))
    expect_false(warned, label = deparse(n))
  }
  for (n in list(-1, -0.5, -Inf)) {
    expect_error(rmhn(n, 2, 1, 1), "n must be non-negative", fixed = TRUE,
                 label = deparse(n))
  }
})

test_that("a flag is read through as.logical(), and one that gives NA is an error", {
  expect_identical(dmhn(1, 2, 1, 1, log = 1), dmhn(1, 2, 1, 1, log = TRUE))
  expect_identical(dmhn(1, 2, 1, 1, log = 0), dmhn(1, 2, 1, 1))
  # NEWS 0.1.1: lower.tail = 1 gave the upper tail.
  expect_identical(pmhn(1, 2, 1, 1, lower.tail = 1), pmhn(1, 2, 1, 1))
  expect_identical(pmhn(1, 2, 1, 1, lower.tail = 0),
                   pmhn(1, 2, 1, 1, lower.tail = FALSE))
  expect_identical(pmhn(1, 2, 1, 1, log.p = 1), pmhn(1, 2, 1, 1, log.p = TRUE))
  expect_identical(pmhn(1, 2, 1, 1, lower.tail = "TRUE"), pmhn(1, 2, 1, 1))
  expect_identical(qmhn(0.3, 2, 1, 1, lower.tail = 0),
                   qmhn(0.3, 2, 1, 1, lower.tail = FALSE))
  expect_identical(qmhn(log(0.3), 2, 1, 1, log.p = 1),
                   qmhn(log(0.3), 2, 1, 1, log.p = TRUE))
  for (v in list(NA, NULL, "1", "yes")) {
    lab <- deparse(v)
    expect_error(dmhn(1, 2, 1, 1, log = v), "'log' must be TRUE or FALSE",
                 fixed = TRUE, label = lab)
    for (f in list(pmhn, qmhn)) {
      expect_error(f(0.5, 2, 1, 1, lower.tail = v),
                   "'lower.tail' must be TRUE or FALSE", fixed = TRUE,
                   label = lab)
      expect_error(f(0.5, 2, 1, 1, log.p = v), "'log.p' must be TRUE or FALSE",
                   fixed = TRUE, label = lab)
    }
  }
})

test_that("attributes of the first argument are not carried to a longer result", {
  # The names of x are kept when the result has its length; a longer result
  # would otherwise get them with NA added.
  expect_null(names(dmhn(c(a = 1), alpha = c(1, 2))))
  expect_null(names(pmhn(c(a = 1), alpha = c(1, 2))))
  expect_null(names(qmhn(c(a = 0.5), alpha = c(1, 2))))
  expect_null(dim(dmhn(matrix(1, 1, 1), alpha = c(1, 2))))
})

test_that("the moment functions and mhn_mode take three single numbers", {
  for (nm in names(scalar_fns)) {
    f <- scalar_fns[[nm]]
    for (v in c(non_numeric, list(c(1, 2), numeric(0)))) {
      lab <- paste(nm, deparse(v))
      expect_error(f(v, 1, 0), "'alpha' must be a single number",
                   fixed = TRUE, label = lab)
      expect_error(f(2, v, 0), "'beta' must be a single number",
                   fixed = TRUE, label = lab)
      expect_error(f(2, 1, v), "'gamma' must be a single number",
                   fixed = TRUE, label = lab)
    }
    expect_identical(f(2L, 1L, 0L), f(2, 1, 0), label = nm)
    expect_identical(f(TRUE, 1, FALSE), f(1, 1, 0), label = nm)
  }
})

test_that("the moment functions and mhn_mode stop on a triple outside the parameter space", {
  for (nm in names(scalar_fns)) {
    f <- scalar_fns[[nm]]
    for (a in c(NA, NaN, 0, -1, Inf, -Inf)) {
      expect_error(f(a, 1, 0), "alpha must be positive and finite",
                   fixed = TRUE, label = paste(nm, a))
    }
    for (b in c(NA, NaN, 0, -1, Inf)) {
      expect_error(f(2, b, 0), "beta must be positive and finite",
                   fixed = TRUE, label = paste(nm, b))
    }
    for (g in c(NA, NaN, Inf, -Inf)) {
      expect_error(f(2, 1, g), "gamma must be a finite numeric value",
                   fixed = TRUE, label = paste(nm, g))
    }
  }
})
