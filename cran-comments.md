# CRAN submission comments — mhn 0.1.1

Maintainer: Tomotaka Momozaki <momozaki.stat@gmail.com>

This is an update (maintenance / bug-fix) release, updating the currently
published version 0.1.0 (on CRAN since 2026-05-18). It is not a new
submission. The headline change is a sampling-correctness fix in `rmhn()`
for `alpha < 1` and `gamma > 0` (previously biased draws); see "Changes in
this version" below.

## Test environments

* local macOS 26.3.1 (Darwin), aarch64-apple-darwin20, R 4.5.1 (2025-06-13)
  — Apple clang 17.0.0. `R CMD check --as-cran` run this session on the
  built `mhn_0.1.1.tar.gz`. Result below.

* win-builder, run on the 0.1.1 tarball via `devtools::check_win_devel()`,
  `check_win_release()` and `check_win_oldrelease()` — **all three
  Status: OK** (0 errors, 0 warnings, 0 notes):
    * R-devel   — R Under development (unstable) (2026-07-23 r90295 ucrt)
    * R-release — R 4.6.1 (2026-06-24 ucrt)
    * R-oldrel  — R 4.5.3 (2026-03-11 ucrt)

The following additional checks have **not** yet been run for 0.1.1 (the
archived runs on record were for 0.1.0); they will be re-run before
submission:

* **[TO RUN]** R-hub v2 via the package's GitHub Actions workflow, on
  `linux`, `windows`, `macos-arm64`, `clang-asan`, `valgrind`, `rchk`
  (results on the Actions tab of <https://github.com/t-momozaki/mhn>). The
  0.1.0 run reported only the long-known Rcpp `Shield<T>` RAII `[PB]`/`[UP]`
  false positive under `rchk`, with no mhn-side function flagged; re-confirm
  on the 0.1.1 run.
* **[TO RUN]** GitHub Actions R-CMD-check matrix: ubuntu-latest (R-devel,
  R-release, R-oldrel-1), macos-latest (R-release), windows-latest
  (R-release).

## R CMD check results

On win-builder the package is fully clean — **0 errors | 0 warnings | 0
notes** on R-devel, R-release and R-oldrelease (see Test environments above).

Local `R CMD check --as-cran` on R 4.5.1 (macOS) reports

  0 errors | 0 warnings | 1-2 notes

`checking CRAN incoming feasibility` is OK locally, and the examples, tests,
and vignette rebuild all pass. The local notes are environment/local-only, do
not appear on win-builder or on CRAN's own check infrastructure, and vary
with the machine's network state:

### NOTE 1 — `checking for future file timestamps ... unable to verify current time`

The wall-clock time-verification endpoint (the remote clock service this
check contacts) was unreachable from the local machine, so the current time
could not be confirmed. This is environment-specific and does not appear on
CRAN's infrastructure.

### NOTE 2 — `checking HTML version of manual ... 'tidy' doesn't look like recent enough HTML Tidy` / `package 'V8' unavailable`

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

The headline change is a sampling-correctness bug fix. A concise summary
(full details in `NEWS.md`):

* **Bug fix (sampling correctness).** `rmhn(method = "rtdr")`, and hence the
  default `method = "auto"` where it routes there, drew biased samples for
  `alpha < 1` and `gamma > 0`: on the log axis the density is only
  T_{-1/2}-concave (concave under the inverse-square transform) in that
  region rather than log-concave, so the previous log-tangent envelope did
  not dominate the target and over-weighted large values. The envelope now
  follows Gao & Wang (2025, Section 3.2 and Appendix B) — an inverse-square
  (T_{-1/2}) tangent hat when `gamma > 0`, the log-tangent hat only when
  `gamma <= 0` — and drawn samples now match the target across the whole
  parameter space (verified by a Kolmogorov-Smirnov and moment
  goodness-of-fit audit). Users who called `rmhn()` with `alpha < 1` and
  `gamma > 0` under 0.1.0 could have obtained biased draws; this release
  corrects that. The density, distribution, quantile, and moment functions
  were not affected.
* **DESCRIPTION.** Spelled out "Markov chain Monte Carlo" and "relaxed
  transformed density rejection method", following CRAN reviewer feedback on
  unexpanded acronyms in 0.1.0.
* **Benchmarks** (`inst/`, not user-facing API). Refined the
  `method = "auto"` dispatch for `gamma < 0` (`alpha >= 10` now uses the Sun
  et al. (2023) Algorithm 3 in the batch regime; the `alpha < 0.1` batch
  cutoff raised from 25 to 100); fixed a unit double-conversion in
  `auto_dispatch.R` that inflated the reported times by roughly 1e6x (the
  dispatch *decisions* are ratio-based and were unaffected, as is
  `rmhn(method = "auto")` itself); added a goodness-of-fit benchmark
  `rmhn_gof.R`; the benchmarks now emit diagnostics CSVs.
* **Examples.** Added `inst/examples/vmf_gibbs.R`, a self-contained Gibbs
  sampler for the von Mises-Fisher concentration parameter whose full
  conditional is an MHN law.

The new `inst/` benchmark and example scripts are not run during
`R CMD check`, need no network, and do not materially affect the installed
size (the built tarball is 304K). None of these changes alters the exported
density / distribution / quantile / moment API behaviour; the only
user-visible change is the `rmhn` sampling-correctness fix above.

## Downstream dependencies

There are no reverse dependencies. Verified 2026-07-24 against the CRAN
package database with `tools::package_dependencies("mhn", reverse = TRUE)`
(across `Depends`, `Imports`, `LinkingTo`, and `Suggests`), which returns
none. Because there are no reverse dependencies, the intentional,
goodness-of-fit-validated change to `rmhn()` output for `alpha < 1`,
`gamma > 0` cannot affect any downstream CRAN package.

## Notes on tests and vignettes

* The default test suite (`tests/testthat/`) completes well within the CRAN
  budget. This session: 2124 pass / 0 fail / 0 warn / 27 skip.
* The 27 skips are heavier goodness-of-fit (Kolmogorov-Smirnov) and
  large-`n` regression blocks in `tests/testthat/test-rmhn*.R`, guarded with
  `skip_on_cran()` so they exercise the package fully in local development
  and on CI while keeping CRAN check time low. The remaining (non-skipped)
  tests still cover the public API, special cases, vectorisation, NA
  propagation, the Sun et al. (2023) and Gao & Wang (2025) algorithm
  dispatch, and the input-validation contract.
* The two vignettes (`introduction.Rmd`, `theory.Rmd`) build in well under
  one minute combined; both are static knitr/rmarkdown documents with no
  network access, no random external data, and seeded RNG calls.
