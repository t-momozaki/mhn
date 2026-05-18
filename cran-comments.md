# CRAN submission comments — mhn 0.1.0

This is a new submission.

## Test environments

* local macOS 26.3.1 (Darwin), R 4.5.1 (2025-06-13) — Apple clang 17.0.0
* win-builder, R-devel and R-release
  (submitted via `devtools::check_win_devel()` and
  `devtools::check_win_release()`; result emails archived with the
  maintainer)
* R-hub v2 via the package's GitHub Actions workflow, on:
  `linux`, `windows`, `macos-arm64`, `clang-asan`, `valgrind`, `rchk`.
  Five of six platforms green (linux / windows / macos-arm64 /
  clang-asan / valgrind); `rchk` reports
  `[PB]` / `[UP]` against `Rcpp::Rcpp_protect` and
  `Rcpp::Rcpp_unprotect` in
  `Rcpp/include/Rcpp/protection/Shield.h:25` / `:31`, the long-known
  Rcpp `Shield<T>` RAII false positive — no mhn-side function is
  flagged. Run name: `some-birdofparadise`; results on the Actions
  tab of <https://github.com/t-momozaki/mhn>.
* GitHub Actions R-CMD-check matrix:
  ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release)

## R CMD check results

Local `R CMD check --as-cran` on R 4.5.1 (macOS) reports

  0 errors | 0 warnings | 1 note

This is the second pretest submission for v0.1.0. The first attempt
(2026-05-18) returned 2 notes; Uwe Ligges responded that the
`Versioned 'LinkingTo' value for 'BH' is only usable in R >= 3.0.2`
note should be addressed by declaring the corresponding (or stricter)
`Depends: R (>= ...)` field. We have added `Depends: R (>= 4.0.0)`,
which is the floor required by `SystemRequirements: C++17` (the
binding constraint among Rcpp >= 1.0.7, BH >= 1.78, testthat
edition 3, and the C++17 `Makevars` standard). This resolves the
LinkingTo note.

The one remaining note is:

### NOTE 1 — `checking CRAN incoming feasibility ... NOTE`

```
Maintainer: 'Tomotaka Momozaki <momozaki.stat@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Gao (...)
  MHN (...)
  RTDR (...)
```

The "New submission" line is expected for an initial submission. The
flagged words are not misspellings:

* "Gao" is the surname of one of the cited authors
  (Gao & Wang, 2025; the `Description` field carries the DOI
  `<doi:10.1080/03610918.2025.2524551>`).
* "MHN" is the standard acronym for the Modified Half-Normal
  distribution that the package provides; it is defined inline in the
  `Description` ("Modified Half-Normal (MHN)").
* "RTDR" is the standard acronym for the Relaxed Transformed Density
  Rejection sampling method (Gao & Wang, 2025), introduced inline in
  the `Description`.

(The previous `Versioned 'LinkingTo' value for 'BH' ...` NOTE is
resolved by the added `Depends: R (>= 4.0.0)`.)

A third note, `checking for future file timestamps ... NOTE` with
"unable to verify current time", is sometimes seen on the local
machine when the wall-clock verification endpoint is unreachable; it
is environment-specific and does not appear on CRAN's own check
infrastructure.

## Downstream dependencies

This is a new submission, so there are no reverse dependencies to
re-check.

## Notes on tests and vignettes

* The default test suite (`tests/testthat/`) completes well within the
  CRAN budget; on `R CMD check --as-cran` the `tests` step takes
  ~10 seconds.
* Heavier goodness-of-fit (Kolmogorov-Smirnov) and large-`n`
  regression blocks in `tests/testthat/test-rmhn*.R` are guarded with
  `skip_on_cran()` so that they exercise the package fully in local
  development and on CI while keeping CRAN check time low. The
  remaining (non-skipped) tests still cover the public API, special
  cases, vectorisation, NA propagation, the Sun et al. (2023) and
  Gao & Wang (2025) algorithm dispatch, and the input-validation
  contract.
* Two vignettes (`introduction.Rmd`, `theory.Rmd`) build in well under
  one minute combined; both are static knitr/rmarkdown documents with
  no network access, no random external data, and seeded RNG calls.
