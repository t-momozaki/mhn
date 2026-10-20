# Each element against its reference, to tol of its size plus abs_tol.
# expect_equal() takes the elements that differ together: where the mean
# size of their references is below the tolerance it compares absolute
# differences, so that a small reference goes unchecked (0 passes for
# 2e-300 at tolerance 1e-8, and for 5.7e-29 beside 0.5 by default), and
# otherwise it compares their mean difference with that mean size, so that
# a small element beside a large one goes unchecked (c(100, 0.5) passes for
# c(101, 0.6) at tolerance 0.1).  An infinite reference must be met exactly,
# and an empty comparison fails.
expect_close <- function(got, ref, tol, abs_tol = 0, info = NULL) {
  ok <- (got == ref) |
    (is.finite(ref) & abs(got - ref) <= tol * abs(ref) + abs_tol)
  worst <- suppressWarnings(
    max(abs(got - ref) / pmax(abs(ref), .Machine$double.xmin)))
  expect_true(length(got) > 0 && length(got) == length(ref) &&
                isTRUE(all(ok)),
              info = paste0(info, " | worst relative difference ",
                            format(worst, digits = 3)))
}
