# CRAN submission comments — mhn 0.1.0

This is a new submission.

## Test environments

* local macOS 26.3.1 (Darwin), R 4.5.1 (2025-06-13) — Apple clang 17.0.0
* win-builder, R-devel and R-release
  (submitted via `devtools::check_win_devel()` and
  `devtools::check_win_release()`; result emails archived with the
  maintainer)
* R-hub v2 via the package's GitHub Actions workflow, on:
  `linux`, `windows`, `macos-arm64`, `clang-asan`, `valgrind`, `rchk`
  (run name: `some-birdofparadise`; results on the Actions tab of
  <https://github.com/t-momozaki/mhn>)
* GitHub Actions R-CMD-check matrix:
  ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release)

## R CMD check results

Local `R CMD check --as-cran` on R 4.5.1 (macOS) reports

  0 errors | 0 warnings | 2 notes

The two notes are:

### NOTE 1 — `checking CRAN incoming feasibility ... NOTE`

```
Maintainer: 'Tomotaka Momozaki <momozaki.stat@gmail.com>'

New submission
```

Expected for an initial submission.

### NOTE 2 — `checking DESCRIPTION meta-information ... NOTE`

```
Versioned 'LinkingTo' value for 'BH' is only usable in R >= 3.0.2
```

`DESCRIPTION` declares `LinkingTo: BH (>= 1.78.0-0)` because the C++
code uses `boost::math::quadrature::gauss_kronrod` (and `tanh_sinh`
for `alpha < 1`) together with `digamma`/`trigamma`, which require a
recent-enough Boost. The version constraint is informational;
`Imports: Rcpp (>= 1.0.7)` already requires an R version far above
the 3.0.2 threshold cited by the note, so the requirement is
satisfied with ample margin.

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
