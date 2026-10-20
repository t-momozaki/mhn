if (requireNamespace("testthat", quietly = TRUE) &&
    utils::packageVersion("testthat") >= "3.1.5") {
  library(testthat)
  library(mhn)

  test_check("mhn")
}
