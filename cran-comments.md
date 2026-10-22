# CRAN submission comments — mhn 0.1.2

Maintainer: Tomotaka Momozaki <momozaki.stat@gmail.com>

This is an update of mhn 0.1.1, which has been on CRAN since 2026-09-23. It
is mainly a bug-fix release; no exported function gains or loses an argument,
and no argument's default value changes.

## Why an update so soon after 0.1.1

0.1.1 returns wrong values without a warning in several places, some of them
at ordinary parameter values, and a few of these came with 0.1.1 itself.
Run on 2026-10-21 with the published 0.1.1 and with this version, against
independent references (arbitrary-precision computations, or the exact Gamma
law the distribution approaches):

* `mhn_kurtosis(0.3, 1, sqrt(10))` returns -0.55609 where the excess kurtosis
  is -0.55730 (0.1.0 was right), and `mhn_kurtosis(1000, 1, 0)` returns
  -0.627 where it is 7.5e-7.
* `pmhn(2e-9, 3, 1, -1e8)` returns 0.0175 where the true value is 0.00115: the
  distribution is Gamma(3, rate 1e8) to within 1e-11 there, and 0.1.1 is wrong
  whenever the mode or the quantile of a strongly tilted law lies below about
  1.5e-8. `qmhn(0.2, 3, 1, -3e8)` returns 8.56e-9 instead of 5.12e-9.
* `pmhn(q, 3, 1, 10, lower.tail = FALSE, log.p = TRUE)` at the 1e-10 quantile
  returns -1.000195e-10 instead of -1e-10 (0.1.0 was right).

The sampler is affected as well: in 0.1.1, `rmhn()` returns `NaN` draws from
Sun et al.'s Algorithm 3 once `beta` exceeds about 1e24 (1998 of the 2000
draws of `rmhn(2000, 50, 1e30, -3e15, method = "sun")`), and for `alpha < 1`
from a standardised tilt of about 1e13 draws that fall on a coarse lattice
(32 distinct values in `rmhn(20000, 0.7, 1, 1e14)`, against 559 in this
version).

## Changes in this version

`NEWS.md` lists every change against 0.1.1, with the values at which each
defect appeared, the effect on seeded draws, and the speed changes. In
outline:

* **Moments and mode.** Skewness and kurtosis at a large shape and for
  `alpha < 1` with an interior peak; the variance near `alpha = 1` and `2`, at
  a strong tilt and at a large shape; the mean at a large shape and under a
  strong negative tilt; the mean and variance where they crossed the bounds
  of Sun et al.'s Lemma 4 next to `alpha = 1`, and the mean where it crossed
  them under a strong negative tilt; the mode at the boundary of Sun et al.'s
  Lemma 3c.
* **Density, distribution and quantile.** The normalising constant at a large
  shape and a strong tilt; `pmhn()` for a tiny mode or quantile, far from the
  peak, and for log tails close to zero; `qmhn()` for tiny probabilities and
  for upper-tail probabilities close to 1; the square-root-Gamma and
  truncated-normal special cases near the origin, and the shapes and tilts
  within 1.5e-8 of them, which 0.1.1 answered as the special cases.
* **Random generation.** Sun et al.'s Algorithms 1 and 3 at extreme shape,
  scale and tilt; exact zeros from underflowing Gamma variates at small
  shapes; the RTDR sampler at a steep tilt, when the mode is tiny, and from
  a shape of about 3e15 under a tilt near `-alpha / e`. The
  default method's choice between the samplers was re-derived from new
  timings, which include shapes up to 1e14 and single draws within one
  vectorised call, so for a fixed seed its draws change where it now picks
  another sampler (`NEWS.md`, Seeded draws).
* **Edges of the double range.** Where 0.1.1 stopped with an internal error
  message (`std::bad_alloc`, quadrature-library errors) or returned a wrong
  finite number at extreme parameter values, the functions now return the
  value of the exact or limiting law, or `NA`, `NaN` or an error that says
  why. The help pages of `dmhn()`, `pmhn()` and `qmhn()` state the range
  over which their results are accurate and what is returned beyond it, and
  that of `rmhn()` which laws it can sample in double precision.
* **Speed.** Most functions are as fast as in 0.1.1 or faster, several by
  orders of magnitude at a strong tilt or a large shape. Some paths are
  slower, as `NEWS.md` details with the triples measured. The log upper tail
  of `pmhn()` where the upper tail alone would lose its relative accuracy
  (below a lower-tail probability of about 1e-5 at ordinary parameters) now
  also evaluates the lower tail, and costs about what the lower tail costs.
  Where that is the series of Sun et al., whose length grows like `gamma^2 /
  beta` and with the shape, it takes per value, for lower-tail probabilities
  from 1e-6 to 1e-4, from about 2 to about 80 times as long as 0.1.1 at the
  triples measured, the ratio growing with `|gamma| / sqrt(beta)` towards 20
  and with the shape (47 against 1.1 microseconds at `(0.5, 1, 19.5)`, 90
  against 1.2 at `(85, 1, -19)`); where it is a quadrature, up to about 2
  times. 0.1.1 is faster there because it takes the logarithm of the upper
  tail alone, whose relative error reaches 0.16. The upper-tail quantile of
  a probability within about 7e-4 of 1, whose root is refined on the lower
  tail, takes 1.3 to 2.6 times as long for `|gamma| / sqrt(beta)` from 15 to
  19.5, 2 to 4 times at shapes from 50 to 85, and at most 1.15 times
  elsewhere. A shape or tilt just off `alpha = 1` or `gamma = 0` now costs
  the general case; the lower tail of `pmhn()` below the median takes up to
  8% longer at some parameters with `|gamma| / sqrt(beta)` up to 20; and a
  single `dmhn()` call with scalar arguments takes up to about 0.1
  microseconds longer. In the timings behind the default sampler choice of
  `rmhn()`, the sampler it now picks is at most 2.3% slower than the one
  0.1.1 picks, except at `alpha = 1.5` with
  `gamma / sqrt(beta) = 2 * sqrt(alpha)`, where it is up to 5.4% slower.
* **Dependencies.** `Imports` and `LinkingTo` ask for `Rcpp (>= 1.0.13)`,
  the first version whose `List::create()` takes more than 20 elements (0.1.1
  asked for 1.0.7 in `Imports` and for any version in `LinkingTo`); with
  Rcpp 1.0.12 the package does not compile. `Suggests` asks for
  `testthat (>= 3.1.5)`, which added `expect_no_error()`.
* **Documentation.** The 0.1.1 NEWS and theory vignette attributed to Sun et
  al. (2023) two misprints their article does not contain; both are
  corrected, and the vignette and NEWS list the six misprints it does
  contain in the results the package computes. `citation("mhn")` gives the
  package's DOI.
* **Scripts under `inst/`.** The result files that 0.1.1 installed under
  `benchmarks/results` and `audits/results` are no longer part of the
  package. The benchmark and audit scripts and the example write their
  results only to a directory named by an environment variable and
  otherwise only print a summary; none of them runs during `R CMD check` or
  needs the network.

## Test environments

All runs below use the tarball built from the submitted sources, or, where
marked, an earlier candidate whose compiled code is identical (the later
changes are to comments, documentation, test descriptions, `DESCRIPTION`
and `inst/CITATION`).

* macOS Tahoe 26.5.2, aarch64-apple-darwin23, R 4.6.1 (2026-06-24), Apple
  clang 17.0.0: `R CMD build`, then `R CMD check --as-cran` of
  `mhn_0.1.2.tar.gz`, with the remote incoming checks (2026-10-22). Also the
  whole test suite on a build with `-ffp-contract=fast`, which lets the
  compiler fuse a multiplication and an addition across statements, as GCC
  does by default on targets with a fused multiply-add (earlier candidate).
* macOS, R 4.5.1, with the lower bounds of the dependencies, Rcpp 1.0.13
  and BH 1.78.0-0: the whole test suite (earlier candidate).
* Linux aarch64 (Docker `rocker/r-ver:4.6.1`, g++ 13.3.0, Rcpp 1.1.2, BH
  1.90.0-1): install without compiler warnings, the whole test suite with no
  failure, and `R CMD check --no-manual --ignore-vignettes`: `Status: OK`
  (earlier candidate). Linux x86_64, the same image under emulation:
  `Status: OK` (earlier candidate).
* Linux x86_64, R 4.0.5 (the declared minimum; Docker `rocker/r-ver:4.0.5`,
  g++ 9.4.0, testthat 3.2.1.1): install without compiler warnings and the
  whole test suite, with no failure (earlier candidate).
* win-builder, R-devel (2026-10-21 r90733 ucrt) and R 4.6.1 ucrt, g++ 14.3.0:
  `Status: OK` for both (2026-10-22 and 2026-10-23), with no compiler warning
  and the tests inside the check at
  `[ FAIL 0 | WARN 0 | SKIP 30 | PASS 12420 ]`.

## R CMD check results

Local `R CMD check --as-cran` of the submitted `mhn_0.1.2.tar.gz` (702,133
bytes; installed size 1.46 MB) on 2026-10-22:

    Status: 1 NOTE

    * checking HTML version of manual ... NOTE
    Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
    Please obtain a recent version of HTML Tidy by downloading a binary
    release or compiling the source code from <https://www.html-tidy.org/>.
    Skipping checking math rendering: package 'V8' unavailable

The local machine lacks a recent HTML Tidy and the V8 package; on
win-builder the same step reports no NOTE.
Every other step is OK, including CRAN incoming feasibility, the tests, the
examples, the re-building of both vignettes and the PDF manual.

## Downstream dependencies

There are no reverse dependencies: on 2026-10-21,
`tools::package_dependencies("mhn", reverse = TRUE)` over `Depends`,
`Imports`, `LinkingTo` and `Suggests` of the CRAN package database (25,404
packages) returned none.

## Tests and vignettes

* Inside `R CMD check --as-cran` the test suite reports
  `[ FAIL 0 | WARN 0 | SKIP 30 | PASS 12420 ]` and takes 29 s. The 30
  skipped blocks, guarded with `skip_on_cran()`, are heavier goodness-of-fit
  and large-sample regression tests, two timing tests, a memory test and a
  comparison to the last bit.
* With `NOT_CRAN=true` the full suite runs: 457 blocks and 15,498
  expectations, with no failure, warning, error or skip.
* The two vignettes are static knitr/rmarkdown documents with no network
  access and seeded random numbers, and build in well under a minute.
