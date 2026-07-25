# Changelog

## mhn 0.1.1

Bug-fix and maintenance release.

### Bug fixes

- `rmhn(method = "rtdr")`, and hence the default `method = "auto"` where
  it routes there, drew biased samples for `alpha < 1` and `gamma > 0`:
  on the log axis the density is only -concave in that region, not
  log-concave, so the previous log-tangent envelope did not dominate it
  and over-weighted large values. The envelope now follows Gao & Wang
  (2025, Section 3.2 and Appendix B), using a (inverse-square) tangent
  hat when `gamma > 0` and the log-tangent hat only when `gamma <= 0`;
  drawn samples now match the target distribution across the whole
  parameter space (verified by a Kolmogorov-Smirnov and moment
  goodness-of-fit audit). The density, distribution, quantile, and
  moment functions were not affected.

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
- `citation("mhn")` returns three `bibentry` objects: the package, the
  Sun et al. (2023) paper, and the Gao & Wang (2025) paper.
