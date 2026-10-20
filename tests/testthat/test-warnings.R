# Warnings raised from the compiled code.
#
# Rcpp::warning() calls Rf_warning().  When such a warning is caught by an
# exiting handler (tryCatch(..., warning = )) or turned into an error
# (options(warn = 2)), R leaves the compiled code by a long jump, the C++
# destructors of the call do not run, and its Rcpp vectors stay protected for
# the rest of the session.  The package raises its warnings through R's
# warning() under Rcpp's unwind protection instead (src/mhn_warn.h).

# A law that dmhn() and pmhn() cannot form: NA with a warning (see "d, p and
# q are NA with a warning where the law cannot be formed" in test-pmhn.R).
unformed <- list(alpha = 1e300, beta = 1e-323, gamma = -1e-160)

test_that("a warning from the compiled code carries its message and no call", {
  w <- tryCatch(
    dmhn(c(1, 1e300), unformed$alpha, unformed$beta, unformed$gamma),
    warning = identity)
  expect_s3_class(w, "warning")
  expect_match(conditionMessage(w), "dmhn: the density could not be evaluated",
               fixed = TRUE)
  expect_null(conditionCall(w))
  as_error <- function() {
    op <- options(warn = 2)
    on.exit(options(op))
    tryCatch(pmhn(c(1, 1e300), unformed$alpha, unformed$beta, unformed$gamma),
             error = identity)
  }
  e <- as_error()
  expect_s3_class(e, "error")
  expect_match(conditionMessage(e), "pmhn: the distribution function could not",
               fixed = TRUE)
})

test_that("a warning that does not return releases the vectors of the call", {
  skip_on_cran()
  # A handler that does not return unwinds the C++ frames under Rcpp's
  # unwind protection, the default from Rcpp 1.0.10 on (DESCRIPTION
  # requires 1.0.13).  Each call below holds vectors of 1e5
  # doubles when it warns, 0.76 Mb each in the units of gc() (2^20 bytes).
  # One of them left behind at every call would add 46 Mb over 60 calls; the
  # vector heap in use after a collection (the second column of gc()) must
  # not grow by 20.
  x <- rep(1, 1e5)
  grown <- function(f) {
    for (i in 1:5) f()
    before <- gc()[2, 2]
    for (i in 1:60) f()
    gc()[2, 2] - before
  }
  expect_lt(grown(function() {
    tryCatch(dmhn(x, unformed$alpha, unformed$beta, unformed$gamma),
             warning = function(w) NULL)
  }), 20)
  expect_lt(grown(function() {
    op <- options(warn = 2)
    on.exit(options(op))
    try(pmhn(x, unformed$alpha, unformed$beta, unformed$gamma), silent = TRUE)
  }), 20)
})
