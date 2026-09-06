# Changelog

## mhn 0.1.1

Bug-fix and maintenance release.

This release corrects a number of numerical defects reachable from the
public API. Most returned a wrong number silently rather than raising an
error, and most require parameter values outside the ranges 0.1.0 was
tested over. Each fix below carries a regression test at the values that
failed.

### What changes relative to 0.1.0

Numbers change wherever a defect is fixed, and that includes part of the
ordinary parameter range — not only the extremes. Everything below was
measured against 0.1.0 as published, over a grid of 330 parameter
triples spanning `alpha` in \[0.05, 20\], `beta` in \[0.01, 100\] and
`gamma` in \[-1e6, 1e6\].

**Deterministic functions.** Over the 130 grid cells with
`|gamma| / sqrt(beta)` at most 10,
[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
and [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
agree with 0.1.0 to twelve significant digits in 100 cells and differ in
30. The largest change is 7.9e-6 relative, in
[`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md);
the rest are below 2e-7. In every case checked against an independent
arbitrary-precision reference the 0.1.1 value is the accurate one — the
Lemma 2c variance already cancels measurably at
`|gamma| / sqrt(beta) = 10`. Outside that range the differences are far
larger, and are the defects listed below.

**[`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) is not
seed-compatible with 0.1.0.** For a fixed seed, 163 of the 330 cells
return a different sequence of draws; 45 of those are in the ordinary
range above. Draws at `gamma = 0` are unchanged, since that case is a
closed-form transform of
[`rgamma()`](https://rdrr.io/r/stats/GammaDist.html). There are two
reasons a sequence changes, and they matter differently:

- Where 0.1.0 sampled the wrong law, the old draws are invalid and any
  analysis resting on them should be re-run. At
  `(alpha, beta, gamma) = (0.3, 1, 100)` a Kolmogorov-Smirnov test of
  the 0.1.0 sample against the true distribution gives `p = 0`; the same
  test on 0.1.1 gives `p = 0.51`.
- Elsewhere the law was already right and only the realisation moves,
  because [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)
  is a rejection sampler: correcting a constant in the envelope changes
  which variates are proposed and which are accepted, so the stream is
  consumed differently. At `(10, 1, 1)` every one of 200,000 draws from
  a fixed seed moves, and yet a Kolmogorov-Smirnov test against
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) gives
  `p = 0.78` on the 0.1.0 sample and `p = 0.70` on the 0.1.1 one, with
  both sample means within 0.05% of the true 2.4387. Analyses in this
  category remain valid; only bit-for-bit reproduction of a stored
  sample is lost.

### Bug fixes — random generation

- `rmhn(method = "rtdr")`, and hence the default `method = "auto"` where
  it routes there, drew biased samples for `alpha < 1` and `gamma > 0`:
  on the log axis the density is only `T_{-1/2}`-concave in that region,
  not log-concave, so the previous log-tangent envelope did not dominate
  it and over-weighted large values. The envelope now follows Gao & Wang
  (2025, Section 3.2 and Appendix B), using a `T_{-1/2}`
  (inverse-square) tangent hat when `gamma > 0` and the log-tangent hat
  only when `gamma <= 0`.

- RTDR drew badly biased samples for `alpha < 1` and `gamma` above
  roughly 100 under every method, since that corner routes to RTDR. An
  exponential envelope piece was evaluated by forming `exp(t)` before
  dividing it out, and region D’s left-tangent piece has `t` growing
  with the tilt — 1068 at `gamma = 100`, past the overflow at 709. The
  piece area became infinite, so piece selection always fell through to
  the last piece, and the widest piece proposed nothing inside its own
  support. At `(0.3, 1, 100)` a Kolmogorov-Smirnov test against
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) gave
  `p = 1.2e-121`. Both the area and the inverse CDF now factor the
  dominant exponential out instead of forming it.

- The RTDR envelope for `alpha <= 1` was decided by rounding once the
  standardised tilt `gamma / sqrt(beta)` passed about 1e8, and answered
  differently on different machines: at `(0.7, 1e-8, 1e8)` macOS
  returned 50 `NaN` while Linux and Windows raised an error, from the
  same input, and the outcome was not even monotone in `gamma`. Near the
  mode the log-axis density carries the constant `gamma_norm^2 / 4` —
  2.5e23 at that tilt — while every quantity the envelope reads out of
  it is an `O(1)` difference of two such values, so the contact drop was
  quantised to multiples of 3.4e7 against a target of `log 4`. Past a
  conditioning threshold the ordinate is now measured from the mode
  through an identity that only ever adds same-signed terms, and the
  accept/reject test is measured the same way — without that second half
  the sampler draws from its own hat. Draws below
  `gamma / sqrt(beta) = 1.6e7` are unchanged; above 1e8 they change and
  the new ones are the correct ones, with a Kolmogorov-Smirnov statistic
  against the exact Gaussian limit of 0.0019 where it had reached 0.384
  and a sample standard deviation within 0.07% of the truth where it had
  been 21 times too large. Regions A and D carry the same defect at the
  same tilts and are not fixed here.

- `rmhn(200, 1.5, 0.01, 3162.3)` returned 200 `NaN`, behind 200 warnings
  that the Algorithm 1 sampler had exhausted its retry budget. Sun et
  al. (2023) Theorem 1a requires the Algorithm 1 proposal scale to lie
  in `(0, beta)`, but the published expression subtracts two quantities
  of size `gamma^2` and leaves that interval for a large tilt — reaching
  0, then going negative, then snapping to `beta`. A negative scale made
  the branch test `NaN < 0` select the sqrt-Gamma proposal, and `rgamma`
  was handed a negative scale. Rationalising the expression removes
  every subtraction.

- At `alpha = 1` with `gamma = -1e8` the sampler returned 81572 negative
  values out of 2e5 — outside the `(0, Inf)` support — and collapsed the
  sample onto 13 distinct doubles. Robert (1995) proposes `y = a + e`
  with `a = -mu/sigma`, and the code returned `mu + sigma * y`; since
  `sigma * a = -mu` that is `sigma * e` with the mean removed
  analytically, but written literally it cancels. It now returns
  `sigma * e`.

- The RTDR region A left-contact search started at `max(m/2, 1e-6)`.
  Whenever the mode was below `2e-6` that floor put the start on the
  wrong side of it, and `rmhn(200, 4, 1, -1e7)` failed outright on the
  default path.

- The RTDR right-contact searches ran on the log axis, where the density
  carries `exp(2y)` and overflows past `y = 354`. The recommended
  starting offset of Gao & Wang (2025, Eq. 8) lands beyond that when the
  mode is flat — 447 above the mode at `alpha = 1e-5`. The search
  returned a non-finite slope, which the sign guard admitted because
  `NaN` fails every comparison, and the plateau piece was dropped with
  no error: the sampler’s support no longer contained the mode or the
  right tail. At `(3e-4, 1, 0.01)` it drew nothing above 0.005, where
  the true probability is 0.0015.

- Corrected the Algorithm 1 envelope constant `K_1`. The main-text
  statement of Sun et al. (2023) Theorem 1a prints the power base as
  `sqrt(beta (alpha-1))`; direct maximisation gives
  `sqrt(beta) (alpha-1)`, which is what the paper’s own proof and
  Theorem 1c carry. The sampled law was always correct — both proposals
  are valid envelopes and each acceptance test is derived independently
  — but the choice between them was not: `log K_1` came out short by
  `(alpha-1) log(alpha-1) / 2`, so for `alpha > 2` the Normal proposal
  was selected across most of the small-tilt half of the region and
  acceptance fell to 0.70, below the 0.8 that Theorem 2e guarantees for
  `alpha >= 4`. It is now 0.82 at worst on that grid.

- Sun Algorithm 1 formed `K_1 - K_2` as a difference of two large terms;
  cancelling them analytically leaves an expression with no subtraction
  of like-sized quantities.

- RTDR region A needs an interior mode in `x`, and at `alpha = 1`
  exactly there is none, so the envelope degenerated to a plateau of
  height `-Inf` spanning `[0, Inf]` and every proposal was rejected.
  `alpha = 1` now uses the region whose mode is strictly positive for
  every `alpha > 0`.

- Sun Algorithm 3’s matching-point fallback for `alpha <= 1.1` contained
  no `beta`, so it was not scale-equivariant and the acceptance rate
  drifted with the scale. When the inflection heuristic left the support
  it discarded the inflection point entirely rather than using it, which
  is the limit of Sun’s own expression.

- `rmhn(method = "sun")` refused parameter values the sun path answers
  in closed form, because the pre-scan tested the raw `alpha` and
  `gamma` while the dispatcher routes on the scale-free predicates.

### Bug fixes — density, distribution, quantile and moments

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) and
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  returned 0 everywhere for `alpha` near 2 with a large tilt — the
  projected-normal `p = 2` case of both source papers. `Psi[1, z]`
  formed `exp(z^2/4)` in real space, which overflows at
  `z = gamma/sqrt(beta) > 53.28`, so the normalising constant became
  infinite. The term is now built as a logarithm.

- The normalising constant for `gamma < 0` degraded from `|z| ~ 100` and
  collapsed to `-Inf` from `|z| ~ 600`, taking everything downstream
  with it:
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md)
  returned `Inf`,
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  returned 1 everywhere,
  [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  `NaN`. The quadrature’s upper limit extended to `|z|/2`, the peak of
  the same expression with the opposite sign, so the mass fell inside
  the first fraction of a percent of the interval. The limit now follows
  the decay the integrand actually has, and the range is split at the
  peak.

- For `gamma > 0` the normalising constant had only the Lemma 10 series,
  which raised an error once its truncation length passed a ceiling — on
  ordinary parameter sets, since that length grows like `z^2` and
  `(2.5, 1e-4, 20)` is already past it. The quadrature that serves
  `gamma < 0` now serves this too.

- The CDF quadrature integrated `[0, x]` in one piece, but for a strong
  tilt the mass lies within a unit or so of the peak, so a fixed-order
  rule returned 0 at the distribution’s own mode. In 0.1.0 that path was
  unreachable, because the series raised an error first; it became
  reachable once the series was given the quadrature fallback described
  above, and both are fixed here. The range is now anchored at the peak,
  split there, and trimmed at both ends.

- `pmhn(..., lower.tail = FALSE)` returned exactly 0, and `log.p = TRUE`
  returned `-Inf`, as soon as the survival probability fell below about
  1e-16, because it was derived as `1 - F`. The upper tail is now
  computed directly.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) at
  `alpha = 1` always computed the upper tail and derived the lower as
  `1 -` that, so `pmhn(3, 1, 1, 18)` gave 0 where the true value is
  1.076e-17. Whichever tail is smaller is now computed directly.

- `qmhn(1e-4, 0.1, 1, 1)` returned 1.4e-17, whose actual probability is
  0.019 — 23 orders of magnitude above the true 2.3e-40, silently. For
  `alpha < 1` the quantile of a small `p` sits near `p^(1/alpha)`, and
  the lower bracket bottomed out at 1.4e-17 after 30 halvings; the
  solver then returned that endpoint as though it were the root. The
  search now runs to the denormal floor.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  returned `+Inf` for every `p` at `alpha = 1` with
  `gamma <~ -11.7 sqrt(beta)`, where the truncation mass was formed as
  `1 - Phi(-mu/sigma)` in linear space and underflowed to exactly 0.

- [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
  returned -701.6 at `(3, 1, 1000)`, for a quantity that cannot fall
  below -2, and
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
  had lost its sign by `gamma = 200`. Both built central moments by
  differencing raw ones, which cancels once the standard deviation is
  small next to the mean. They are now integrated directly, in the
  variable centred on the peak so that no step forms a difference larger
  than the answer.

- [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  was 31% wrong at `(2.5, 1, 1e8)`, differencing two `log Psi` values
  that are each about `z^2/4` when their difference is about `log z`.
  Nothing about that error is systematic: it is the last bits of a total
  cancellation, so it lands at +31% for some shapes and -20% for others,
  and by `gamma = 1e9` it reaches a factor of 1.6e5. Past the point
  where the error would exceed the working tolerance the mean now comes
  from quadrature.

- `log Psi[1/2, z]` carried an absolute error of 0.73 at `z = -1e8`, in
  a quantity whose true value is -17.73: `z^2/4` and `log Phi(z/sqrt 2)`
  cancel. Substituting the scaled complementary error function cancels
  the exponential analytically.

- The mode formula cancels for `gamma < 0`; Gao & Wang (2025, Eq. 9)
  give the conjugate form for that reason and their reference code
  selects on the sign. `mhn_mode(1.1, 1, -1e8)` returned 0 instead of
  1e-9. Five sites used only the unstable form.

- Both places that asked whether the tilt is negligible compared
  `|gamma|` against `sqrt(.Machine$double.eps)` on an absolute scale,
  but the family depends on `gamma` only through
  `Delta = gamma / sqrt(beta)` (Sun et al. 2023, Theorem 1c). At
  `beta = 4e-18` a `gamma` of 1e-8 is an ordinary `Delta = 5`, and every
  function returned the untilted sqrt-Gamma answer:
  `mhn_mean(1, 4e-18, 1e-8)` gave 2.82e8 against a true 1.2503e9. Both
  tests now use `Delta`.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) with
  `gamma = -1e4` requested about 5e7 series terms, taking tens of
  seconds and hundreds of megabytes. The truncation length grows like
  `z^2` and had no ceiling; past one the series now returns the sentinel
  that already routes the caller to quadrature.

- The normalising constant for `gamma < 0` was also wrong by a factor
  that grew without bound as `alpha` fell — 1545 at `alpha = 1e-6`,
  still 0.17% at `alpha = 0.01`. Half the mass of `x^(alpha-1)` lies
  below `x = exp(-1/alpha)`, which for small `alpha` is below any
  representable double, so no quadrature can reach it. The first stretch
  of the range is now integrated in closed form. The distribution
  function and the central-moment quadrature had the same gap; it was
  invisible while the normalising constant shared it, and appeared as
  soon as that was fixed. `pmhn(1e-3, 1e-6, 1, -1000)` returned 6.3e-4
  where the true value is 1.

- `pmhn(..., lower.tail = FALSE)` returned exactly 0 for every `q` at
  small `alpha` with a strong negative tilt, and `qmhn` inverted that by
  returning its own initial bracket — the same number for every `p`,
  with no warning. The quadrature range was set from the curvature at
  the peak, which is the wrong scale when the range is anchored anywhere
  else: at `beta = 1e-4` with `gamma = -1000` the Gaussian width is 894
  while the kernel dies within 0.08. The range now follows the slope at
  the anchor. `qmhn` no longer returns an unjustified value: it warns
  and returns `NA` when it cannot bracket the root.

- [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
  was 17% high at `gamma = -1e7` and clamped to exactly 0 by
  `gamma = -1e8`. The Lemma 2c form is exact algebra, but its leading
  terms are each about `alpha/2` while the answer is about
  `alpha/gamma^2`, so it amplifies its own rounding by
  `gamma^2/(2 beta alpha)`. It is now used only where the measured
  cancellation is small enough, and the second central moment comes from
  quadrature elsewhere — the closed form remains the more accurate of
  the two wherever it is well conditioned.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  implemented its special-case dispatch twice and the two disagreed. The
  vector-parameter path derived the truncated-normal lower tail as
  `1 - upper`, which is exactly 0 once the lower tail falls below the
  rounding of 1: `pmhn(3.0184, 1, 1, 18)` was 1.35e-17 element-wise and
  0 vectorised.

### Interface

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) and
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) mapped
  a `NaN` input to `NA`, where
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) and the
  base R `d`/`p`/`q` family return `NaN`. All three now distinguish the
  two on every path.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)
  documents an `NA` draw for a non-finite parameter but guarded only
  `gamma`; an infinite `alpha` or `beta` reached the sampler.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) took
  `n[1]` where every base R `r*` function reads a vector `n` as
  `length(n)`, so `rmhn(c(10, 20, 30))` drew 10 variates instead of 3.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) tested
  `lower.tail` and `log.p` with
  [`isTRUE()`](https://rdrr.io/r/base/Logic.html), which is stricter
  than base R: `pmhn(1, 2, 1, 1, lower.tail = 1)` returned the upper
  tail, 0.5661, where `lower.tail = TRUE` gives 0.4339 — and base R’s
  `pnorm(1, lower.tail = 1)` agrees with `lower.tail = TRUE`. All four
  functions now coerce the flag as base R does, and reject one that
  resolves to `NA`.

- Errors raised when the sampler cannot build an envelope named private
  C++ routines and the sign of an intermediate slope. They now report
  what failed, at which parameter values, and that it should be
  reported.

### Performance

- Twenty of the twenty-five compiled entry points never touch the RNG
  and were still generated with an `Rcpp::RNGScope`, paying
  `GetRNGstate`/`PutRNGstate` on every call to `dmhn`, `pmhn`, `qmhn`
  and the moment helpers. Draws are unchanged.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  evaluated the normalising constant a second time for a triple whose
  CDF state had just computed it.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) skips
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) when `method`
  was not supplied; it accounted for about 55% of the cost of the
  single-variate call a Gibbs sweep makes.

- Piece selection in the RTDR sampler rebuilt a cumulative-area table on
  every proposal. The table is now built once per envelope and searched
  by bisection.

### Source

- Removed development scaffolding from the installed source: references
  to files not in the tarball, plan markers describing a half-finished
  file, and stale claims of mirroring an R implementation that no longer
  exists. Removed three unused functions and one unused struct field.

- Added `src/mhn_stable.h`, which holds the better-conditioned algebraic
  forms in one place rather than repeating them at the call sites that
  need them.

- Added interrupt checking to the four driver loops and the two
  rejection loops in the truncated-normal sampler; a long
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) call
  was uninterruptible.

- Moved the retry-budget warnings out of the sampler frames.
  `Rcpp::warning` is a bare `Rf_warning`, which under
  `options(warn = 2)` leaves the frame by a long jump without running
  destructors while the cached envelope holds five vectors.

### Tests

- The suite grows from about 2,100 expectations to about 6,300, of which
  6,066 run under `R CMD check` and 37 blocks are skipped there, with a
  regression block for each defect above at the parameter values that
  failed. New coverage: non-unit `beta` throughout the sampler grids and
  a direct check of scale equivariance; the contact-point search in
  isolation; the `NA`/`NaN` contract against base R; and the accuracy
  audit compared against a reference that shares no code with the
  package.

### DESCRIPTION

- Spelled out “Markov chain Monte Carlo” and “relaxed transformed
  density rejection method” in the `Description` field, following CRAN
  reviewer feedback on unexpanded acronyms.

### Benchmarks

- Refined the `method = "auto"` dispatch for `gamma < 0`. Previously any
  batch of 25 or more variates per setup used RTDR; benchmarking across
  three independent runs showed that for `alpha >= 10` the Sun et al.
  2023. Algorithm 3 has the lower per-proposal cost and wins in the
        batch regime too, so `auto` now keeps RTDR for large batches
        only when `alpha < 10`. This lowers the worst-case slowdown of
        `auto` relative to the per-cell optimum from about 11% to about
        4% while leaving the common Gibbs (single-variate) path
        unchanged. `inst/benchmarks/auto_dispatch.R` gained an
        `alpha < 1` grid, a setup/per-proposal cost decomposition, and a
        comparison of the shipped rule against the measured optimum.
- Raised the `gamma < 0` batch cutoff from 25 to 100 variates per setup
  when `alpha < 0.1`. The crossover between Algorithm 3 and RTDR moves
  to larger batches as the shape shrinks — it sits near 25 for `alpha`
  around 0.8 but near 100 by `alpha = 0.01` — so the old cutoff sent 25
  to 99 variates per setup to RTDR while Algorithm 3 was still up to 10%
  faster there.
- Fixed a unit double-conversion in `inst/benchmarks/auto_dispatch.R`
  that inflated the reported `median_us` / `iqr_us` times by a factor of
  about 1e6. The `method = "auto"` dispatch *decisions* are ratio-based
  and were unaffected, as is `rmhn(method = "auto")` itself.
- Added a goodness-of-fit benchmark, `inst/benchmarks/rmhn_gof.R`, that
  writes Kolmogorov-Smirnov statistics and sample-vs-theory moment
  summaries across the parameter grid to a CSV.
- The two timing benchmarks now emit a `_diagnostics_<date>.csv` with
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), hardware,
  and `mhn` version provenance, matching the audit scripts.

### Examples

- Added `inst/examples/vmf_gibbs.R`, a self-contained Gibbs sampler for
  the von Mises-Fisher concentration parameter whose full conditional is
  an MHN law. It is run with
  `source(system.file("examples", "vmf_gibbs.R", package = "mhn"))`,
  takes its sample size, chain length and true concentration from
  `MHN_VMF_*` environment variables, and reports interval coverage and
  effective sample size.

## mhn 0.1.0

CRAN release: 2026-05-27

Initial release.

### Distribution functions

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md), and
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) provide
  density, distribution, quantile, and random generation for the
  Modified Half-Normal (MHN) distribution of Sun, Kong & Pal (2023).
- All four functions are vectorised over both the evaluation argument
  and the parameters `alpha`, `beta`, `gamma`, following standard R
  recycling rules.
- A ParamCache reuses the Fox–Wright Psi normalising constant across
  consecutive elements that share an (`alpha`, `beta`, `gamma`) triple,
  so grouped inputs are evaluated significantly faster than calling the
  functions inside an R loop.

### Random generation

- `rmhn(..., method = "auto")` (default) routes each parameter triple to
  the cheapest provably-correct sampler: closed-form shortcuts for the
  special cases, Sun et al. (2023, Algorithms 1 and 3) where they win,
  and the Gao & Wang (2025) Relaxed Transformed Density Rejection (RTDR)
  sampler elsewhere.
- `method = "rtdr"` forces RTDR with its uniform 1/e acceptance bound.
- `method = "sun"` forces Sun Algorithm 1 (`gamma > 0, alpha > 1`) or
  Algorithm 3 (`gamma <= 0`); Sun Algorithm 2 is intentionally not
  implemented and an unsupported combination triggers a clear error.

### CDF and quantile

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) uses
  the Sun et al. (2023, Lemma 1b) series in log space, truncated at the
  Sun et al. (2023, Supplementary Lemma 10(d)) constructive bound K =
  max(K1, K2); the truncation residual is bounded by the user’s
  tolerance divided by `Psi`.
- For `gamma < 0` the series uses sign-separated log-sum-exp + log-
  diff-exp accumulation and a runtime cancellation guard derived from
  the double-precision precision floor: when the relative cancellation
  loss would exceed the user’s tolerance,
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) falls
  back to a peak-normalised Boost.Math quadrature (Gauss-Kronrod for
  `alpha >= 1`, tanh-sinh for `alpha < 1`) of the unnormalised density.
- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) inverts
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) via
  `boost::math::tools::toms748_solve` on the bracket
  `[sqrt(eps), E(X) + 8 sqrt(Var(X))]`, doubling the upper end as
  needed.

### Summary statistics

- [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
  [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md),
  [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md),
  and
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  evaluate the closed-form / recurrence-based expressions from Sun et
  al. (2023, Lemmas 2 and 3).

### Tests and documentation

- testthat suite with \> 1,700 expectations covering goodness-of-fit
  (Kolmogorov-Smirnov), special-case identities, vectorised recycling,
  NA / NaN propagation, and the Sun / Gao & Wang acceptance bounds.
- [`vignette("introduction", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/introduction.md)
  walks through every exported function with runnable examples.
- [`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md)
  is the theoretical companion: it covers the MHN family and its special
  cases, the Fox–Wright Psi normalising constant, Algorithms 1 and 3 of
  Sun et al. (2023), the four-region Gao & Wang RTDR construction, and
  the `rmhn(method = "auto")` decision tree.
- `citation("mhn")` returns four `bibentry` objects: the package, the
  Sun et al. (2023) paper, the Gao & Wang (2025) paper, and
  Robert (1995) for the truncated-normal sampler used at `alpha = 1`.
