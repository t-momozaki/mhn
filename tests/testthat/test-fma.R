# The package's own fused multiply-add, used where the C library's fma does not
# round x * y + z once: mingw-w64 v11's fma on Windows, which Rtools45 links
# statically, adds partial products with a rounding at each step.  The tests
# call it through .fma_soft on every platform and compare it with x * y + z
# rounded once, formed in Rmpfr with enough bits to hold it exactly.

exact_fma <- function(x, y, z) {
  # 4200 bits hold x * y + z exactly for any finite doubles: the product of
  # two doubles spans at most 2 * 1074 + 2048 bits from its lowest possible
  # bit to its highest, and the sum adds one.
  X <- Rmpfr::mpfr(x, 4200); Y <- Rmpfr::mpfr(y, 4200); Z <- Rmpfr::mpfr(z, 4200)
  Rmpfr::asNumeric(X * Y + Z)
}

same_double <- function(a, b) {
  # Bitwise, so that -0 and 0 are told apart; any NaN matches any NaN.
  identical(writeBin(a, raw()), writeBin(b, raw())) || (is.nan(a) && is.nan(b))
}

expect_fma <- function(x, y, z, label) {
  got <- mhn:::.fma_soft(x, y, z)
  ref <- vapply(seq_along(x), function(i) exact_fma(x[i], y[i], z[i]), 0)
  ok <- vapply(seq_along(x), function(i) same_double(got[i], ref[i]), TRUE)
  first <- which(!ok)[1]
  expect_true(all(ok), label = sprintf("%s: %d of %d differ, first at x=%a y=%a z=%a",
                                       label, sum(!ok), length(ok),
                                       x[first], y[first], z[first]))
}

test_that("the residual of a product is exact", {
  skip_if_not_installed("Rmpfr")
  set.seed(1)
  x <- c(1 + 2^-30, 0.1, 1 / 3, 1.0609e60, 7.65e208, 2^1023,
         runif(200, 0.5, 2) * 2^sample(-500:500, 200, TRUE))
  y <- c(1 + 2^-30, 0.3, 3, 1.0609e60, 1.482e-23, 1.5,
         runif(200, 0.5, 2) * 2^sample(-500:500, 200, TRUE))
  expect_fma(x, y, -(x * y), "x y - fl(x y)")
})

test_that("cancellation, ties and the sign of zero follow one rounding", {
  skip_if_not_installed("Rmpfr")
  u <- 2^-52
  # (1 + u)(1 - u) = 1 - u^2 lies between two doubles; (1 + u)^2 - 1 and
  # (1 + u)(1 - u) - 1 cancel to 2 u + u^2 and -u^2; (1 + 2u)(1 + u/2) lies
  # u^2 above the midpoint 1 + 2.5 u; 3 fl(1/3) - 1 cancels to -u/4; the
  # smallest subnormal times one half is a tie, which goes to the even 0.
  # With a zero addend the product is rounded by the multiplication alone;
  # the ties met in the integer arithmetic are in the next block.
  x <- c(1 + u, 1 + u, 1 + u, 1 + 2 * u, 3, 2^-1074, 2^-1074 * 3)
  y <- c(1 - u, 1 + u, 1 - u, 1 + u / 2, 1 / 3, 0.5, 0.5)
  z <- c(0, -1, -1, 0, -1, 0, 0)
  expect_fma(x, y, z, "hand-picked")
  expect_fma(c(1, -1), c(1, 1), c(-1, 1), "exact zero sums")
  expect_fma(c(-0, 0, -0, 1), c(1, -1, -1, -2^-1074), c(-0, -0, 0, 0),
             "signed zeros")
})

test_that("each step of the integer arithmetic rounds once", {
  # Operands that reach each step of fma_soft where a wrong step would still
  # give a double near the right one, and where random operands seldom go.
  # The expected values are x * y + z rounded once (Rmpfr and a hardware fma
  # agree on them); they are written out so that the block runs without
  # Rmpfr, as on a platform where fma_soft serves every fused multiply-add.
  cases <- list(
    # The bits shifted out of the smaller term, below the sum's last place,
    # decide the rounding.
    "shifted-out bits" = c(0x1.769d5dc175517p-53, -0x1.5de2836261c2bp-1,
                           1.5, 0x1.7ffffffffffffp+0),
    # The product (1 + 2^-26)(1 + 2^-27) is a midpoint; an addend more than
    # 128 places below it, kept only as a sticky bit, rounds it up.
    "addend past 128 bits" = c(1 + 2^-26, 1 + 2^-27, 2^-300,
                               0x1.0000006000001p+0),
    # A shift of exactly one 64-bit word, whose bits shifted out are
    # gathered without shifting a word by 64.
    "shift of one word" = c(1.5, 1.5, 2^-63, 2.25),
    "shift of one word, 2" = c(1.5, 1.5, 1.25 * 2^-63, 2.25),
    # A sum that reaches 2^128 is halved, the bit shifted out kept.
    "sum reaching 2^128" = c(0x1.769d5dc175517p-52, 0x1.5de2836261c2bp+0,
                             2 - 2^-52, 0x1.0000000000001p+1),
    # A carry from the low word into the high one.
    "carry from the low word" = c(-0x1.650c4608b01ddp-761,
                                  0x1.c0f6d47fd227p+83,
                                  -0x1.c2912f8764668p-735,
                                  -0x1.3916de5d785efp-677),
    # A carry from the low word through a high word of all ones, by which
    # the sum reaches 2^128.
    "carry through all ones" = c(0x1.59a8a5548598p+0, 0x1.7b257731a897ep+0,
                                 0x1.0bceb873afc55p-12, 2),
    # 1 - 2^-54 is a midpoint; rounded to even it is 1, a mantissa of 2^53
    # that carries into the exponent.
    "rounding into the exponent" = c(1, 1, -2^-54, 1),
    # 1.5 * 2^1024 lies past the largest double.
    "overflow" = c(2^1000, 1.5 * 2^24, 1, Inf),
    # A finite product beyond the doubles, with an infinite addend.
    "infinite addend" = c(1e300, 1e300, -Inf, -Inf),
    # A zero factor with a finite one: the sum of two signed zeros.
    "zero factor, +0" = c(1, -0, 0, 0),
    "zero factor, -0" = c(-1, 0, -0, -0)
  )
  for (name in names(cases)) {
    v <- cases[[name]]
    got <- mhn:::.fma_soft(v[1], v[2], v[3])
    expect_true(same_double(got, v[4]),
                label = sprintf("%s: got %a for %a", name, got, v[4]))
  }
})

test_that("a product beyond the doubles is brought back by the addend", {
  skip_if_not_installed("Rmpfr")
  big <- .Machine$double.xmax
  x <- c(big, big, 2^600, 1e308, 1.0609e300)
  y <- c(2, 1.5, 2^500, 10, 1e10)
  z <- c(-big, -big, -big, -big, -big)
  expect_fma(x, y, z, "overflowing products")
  expect_fma(c(big, 2^1000), c(2, 2^30), c(0, 1), "overflow to infinity")
})

test_that("subnormal results and operands are rounded once", {
  skip_if_not_installed("Rmpfr")
  set.seed(2)
  tiny <- 2^-1074
  x <- c(tiny, 3 * tiny, 2^-600, 2^-537 * (1 + 2^-52), 1.5 * 2^-1022,
         runif(100, 0.5, 2) * 2^sample(-560:-500, 100, TRUE))
  y <- c(0.5, 0.5, 2^-500, 2^-537, 0.75,
         runif(100, 0.5, 2) * 2^sample(-560:-500, 100, TRUE))
  z <- c(0, tiny, 2^-1073, -tiny, 2^-1022,
         sample(c(-1, 1), 100, TRUE) * runif(100) * 2^-1060)
  expect_fma(x, y, z, "subnormal")
})

test_that("random operands across the range are rounded once", {
  skip_if_not_installed("Rmpfr")
  set.seed(3)
  n <- 2000
  e <- function() pmax(sample(-1100:1023, n, TRUE), -1074)
  s <- function() sample(c(-1, 1), n, TRUE)
  x <- s() * runif(n, 1, 2) * 2^e()
  y <- s() * runif(n, 1, 2) * 2^e()
  # Half of the addends lie next to the product and cancel it; the rest are
  # anywhere in the range.
  near <- sample(c(TRUE, FALSE), n, TRUE)
  z <- ifelse(near, -(x * y) * (1 + s() * runif(n) * 2^-40),
              s() * runif(n, 1, 2) * 2^e())
  keep <- is.finite(x * y) | is.finite(z)
  expect_fma(x[keep], y[keep], z[keep], "random")
})

test_that("infinities and NaN give what IEEE 754 fma gives", {
  got <- mhn:::.fma_soft(c(Inf, Inf, 0, 1, NaN, 1, Inf),
                         c(1, 0, Inf, 1, 1, 1, 1),
                         c(1, 1, 1, NaN, 1, Inf, -Inf))
  expect_identical(got[c(1, 6)], c(Inf, Inf))
  expect_true(all(is.nan(got[c(2, 3, 4, 5, 7)])))
})
