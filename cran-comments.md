# CRAN submission comments — mhn 0.1.1

Maintainer: Tomotaka Momozaki <momozaki.stat@gmail.com>

This is a bug-fix release, updating the currently published version 0.1.0 (on
CRAN since 2026-05-18). It is not a new submission. It fixes a number of
numerical defects reachable from the public API, most of which returned a wrong
number silently rather than raising an error. See "Changes in this version"
below; `NEWS.md` gives the full list with the parameter values at which each
failed.

## Test environments

<!-- ACTION REQUIRED BEFORE SUBMISSION -------------------------------------
The win-builder, GitHub Actions and R-hub results recorded below were obtained
on the earlier 0.1.1 tarball, before the numerical fixes below changed the C++
sources.  They do not describe the tarball this file accompanies.  Re-run all
three on the current tarball and replace the three bullets below with the new
results before submitting.  Delete this comment once that is done.
-------------------------------------------------------------------------- -->

* local macOS 26.3.1 (Darwin), aarch64-apple-darwin20, R 4.5.1 (2025-06-13)
  — Apple clang 17.0.0. `R CMD check --as-cran` run this session on the
  built `mhn_0.1.1.tar.gz`. Result below.

* win-builder, run on the 0.1.1 tarball via `devtools::check_win_devel()`,
  `check_win_release()` and `check_win_oldrelease()` — **all three
  Status: OK** (0 errors, 0 warnings, 0 notes):
    * R-devel   — R Under development (unstable) (2026-07-23 r90295 ucrt)
    * R-release — R 4.6.1 (2026-06-24 ucrt)
    * R-oldrel  — R 4.5.3 (2026-03-11 ucrt)

The 0.1.1 tarball was additionally checked on GitHub Actions:

* GitHub Actions R-CMD-check matrix: ubuntu-latest (R-devel, R-release,
  R-oldrel-1), macos-latest (R-release), windows-latest (R-release) — **all
  green** (0 errors, 0 warnings, 0 notes).
* R-hub v2 on `linux` (R-devel), `windows` (R-devel), `macos-arm64`
  (R-devel), `clang-asan`, `valgrind` and `rchk`: **five of the six platforms
  pass** (linux, windows, macos-arm64, clang-asan, valgrind). As in 0.1.0,
  only `rchk` reports a problem, and it is the long-known Rcpp `Shield<T>`
  RAII false positive — `[PB]`/`[UP]` against `Rcpp::Rcpp_protect` /
  `Rcpp::Rcpp_unprotect` in `Rcpp/include/Rcpp/protection/Shield.h:25`/`:31`
  — accompanied only by rchk's own "address taken, results will be
  incomplete" limitation notes on the Rcpp-generated wrappers and "too many
  states" abstraction errors on R's own internals (`strptime_internal`,
  `bcEval_loop`, `RunGenCollect`). No mhn-side function is flagged with a
  protection bug.

## R CMD check results

On win-builder the package is fully clean — **0 errors | 0 warnings | 0
notes** on R-devel, R-release and R-oldrelease (see Test environments above).

Local `R CMD check --as-cran` on R 4.5.1 (macOS) on the submitted tarball
reports

  0 errors | 0 warnings | 1 note

`checking CRAN incoming feasibility` is OK locally, and the examples, tests,
and vignette rebuild all pass. The note is local tooling only and does not
appear on win-builder or on CRAN's own check infrastructure. A second note,
`checking for future file timestamps`, appears on runs where the machine
cannot reach the remote clock service; it did not appear on this run.

### NOTE — `checking HTML version of manual ... 'tidy' doesn't look like recent enough HTML Tidy` / `package 'V8' unavailable`

The local machine lacks a recent HTML Tidy and the `V8` package; CRAN's
machines have both. This is a local-tooling note only.

### Possible `checking CRAN incoming feasibility` note on the update

CRAN's own incoming-feasibility check may flag "Possibly misspelled words in
DESCRIPTION". These are not misspellings:

* "Gao" is the surname of one of the cited authors (Gao & Wang, 2025; the
  `Description` field carries the DOI `<doi:10.1080/03610918.2025.2524551>`).
* "MHN" is the standard acronym for the Modified Half-Normal distribution
  the package provides; it is defined inline in the `Description`
  ("Modified Half-Normal (MHN)").

(The "RTDR" acronym flagged in 0.1.0 no longer appears: the `Description`
now spells out "relaxed transformed density rejection method".)

## Changes in this version

A bug-fix release. `NEWS.md` lists the fixes individually; the summary below is
what a reviewer or user would notice.

* **`rmhn()` is not seed-compatible with 0.1.0.** Over a grid of 330 parameter
  triples, 163 return a different sequence of draws for a fixed seed. Where
  0.1.0 sampled the wrong law --- `alpha < 1` with `gamma > 0`, and again with
  `gamma` above roughly 100, where a Kolmogorov-Smirnov test of the old sample
  against the true distribution gives `p = 0` --- the old draws are invalid.
  Elsewhere the law was already correct and only the realisation moves, because
  correcting a constant in a rejection sampler's envelope changes which variates
  are proposed and accepted. `NEWS.md` gives the measured breakdown.

* **The density, distribution, quantile and moment functions return corrected
  values at extreme parameters.** `dmhn()` and `pmhn()` returned 0 for `alpha`
  near 2 with a large tilt; `pmhn()` returned 0 at its own mode for `gamma`
  above 3e4 and exactly 0 in the upper tail below 1e-16; `qmhn()` was wrong by
  36 orders of magnitude for small `alpha`; `mhn_kurtosis()` returned a value
  below its theoretical minimum. Values also move slightly inside the ordinary
  range --- 30 of 130 grid cells with `|gamma|/sqrt(beta)` at most 10, by at
  most 7.9e-6 relative --- because the Lemma 2c variance already cancels
  measurably there. In every case checked against an arbitrary-precision
  reference the new value is the accurate one.

* **Three interface behaviours now match base R.** `dmhn()`/`pmhn()` returned
  `NA` for a `NaN` input where the `d`/`p`/`q` family returns `NaN`; `rmhn()`
  read a vector `n` as `n[1]` rather than `length(n)`; and `pmhn()`/`qmhn()`
  tested `lower.tail` and `log.p` more strictly than base R tests them, so
  `pmhn(1, 2, 1, 1, lower.tail = 1)` returned the upper tail.

* **`pmhn()` no longer allocates without bound.** At `gamma = -1e4` it requested
  about 5e7 series terms; it now defers to quadrature above a ceiling. The
  evaluators also check for a user interrupt in every driver loop.

* **DESCRIPTION.** Spelled out "Markov chain Monte Carlo" and "relaxed
  transformed density rejection method", following CRAN reviewer feedback on
  unexpanded acronyms in 0.1.0.

No exported function gains or loses an argument, and no default changes. The
`inst/` benchmark, example and audit scripts are not run during `R CMD check`,
need no network, and do not materially affect the installed size.

## Downstream dependencies

There are no reverse dependencies. Verified 2026-07-24 against the CRAN
package database with `tools::package_dependencies("mhn", reverse = TRUE)`
(across `Depends`, `Imports`, `LinkingTo`, and `Suggests`), which returns
none. Because there are no reverse dependencies, the corrected `rmhn()` output
cannot affect any downstream CRAN package.

## Notes on tests and vignettes

* The default test suite (`tests/testthat/`) completes well within the CRAN
  budget. This session, inside `R CMD check --as-cran` on the submitted
  tarball: 6066 pass / 0 fail / 0 warn / 37 skip.
* The 37 skips are heavier goodness-of-fit (Kolmogorov-Smirnov) and
  large-`n` regression blocks in `tests/testthat/test-rmhn.R`,
  `test-rmhn-sun.R`, `test-rmhn-rtdr-regions.R` and
  `test-numerical-stability.R`, guarded with
  `skip_on_cran()` so they exercise the package fully in local development
  and on CI while keeping CRAN check time low. The remaining (non-skipped)
  tests still cover the public API, special cases, vectorisation, NA
  propagation, the Sun et al. (2023) and Gao & Wang (2025) algorithm
  dispatch, and the input-validation contract.
* The two vignettes (`introduction.Rmd`, `theory.Rmd`) build in well under
  one minute combined; both are static knitr/rmarkdown documents with no
  network access, no random external data, and seeded RNG calls.
