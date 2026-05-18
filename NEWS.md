# mhn 0.1.0

Initial release.

## Distribution functions

* `dmhn()`, `pmhn()`, `qmhn()`, and `rmhn()` provide density,
  distribution, quantile, and random generation for the Modified
  Half-Normal (MHN) distribution of Sun, Kong & Pal (2023).
* All four functions are vectorised over both the evaluation argument
  and the parameters `alpha`, `beta`, `gamma`, following standard R
  recycling rules.
* A ParamCache reuses the Fox--Wright Psi normalising constant across
  consecutive elements that share an (`alpha`, `beta`, `gamma`) triple,
  so grouped inputs are evaluated significantly faster than calling the
  functions inside an R loop.

## Random generation

* `rmhn(..., method = "auto")` (default) routes each parameter triple
  to the cheapest provably-correct sampler: closed-form shortcuts for
  the special cases, Sun et al. (2023, Algorithms 1 and 3) where they
  win,
  and the Gao & Wang (2025) Relaxed Transformed Density Rejection
  (RTDR) sampler elsewhere.
* `method = "rtdr"` forces RTDR with its uniform 1/e acceptance bound.
* `method = "sun"` forces Sun Algorithm 1 (`gamma > 0, alpha > 1`) or
  Algorithm 3 (`gamma <= 0`); Sun Algorithm 2 is intentionally not
  implemented and an unsupported combination triggers a clear error.

## CDF and quantile

* `pmhn()` uses the Sun et al. (2023, Lemma 1b) series in log space,
  truncated at the Sun et al. (2023, Supplementary Lemma 10(d))
  constructive bound K = max(K1, K2); the truncation residual is
  bounded by the user's tolerance divided by `Psi`.
* For `gamma < 0` the series uses sign-separated log-sum-exp + log-
  diff-exp accumulation and a runtime cancellation guard derived from
  the double-precision precision floor: when the relative cancellation
  loss would exceed the user's tolerance, `pmhn()` falls back to a
  peak-normalised Boost.Math quadrature (Gauss-Kronrod for `alpha >= 1`,
  tanh-sinh for `alpha < 1`) of the unnormalised density.
* `qmhn()` inverts `pmhn()` via `boost::math::tools::toms748_solve` on
  the bracket `[sqrt(eps), E(X) + 8 sqrt(Var(X))]`, doubling the upper
  end as needed.

## Summary statistics

* `mhn_mean()`, `mhn_var()`, `mhn_skewness()`, `mhn_kurtosis()`, and
  `mhn_mode()` evaluate the closed-form / recurrence-based expressions
  from Sun et al. (2023, Lemmas 2 and 3).

## Tests and documentation

* testthat suite with > 1,700 expectations covering goodness-of-fit
  (Kolmogorov-Smirnov), special-case identities, vectorised recycling,
  NA / NaN propagation, and the Sun / Gao & Wang acceptance bounds.
* `vignette("introduction", package = "mhn")` walks through every
  exported function with runnable examples.
* `vignette("theory", package = "mhn")` is the theoretical companion:
  it covers the MHN family and its special cases, the Fox--Wright Psi
  normalising constant, Algorithms 1 and 3 of Sun et al. (2023),
  the four-region
  Gao & Wang RTDR construction, and the `rmhn(method = "auto")`
  decision tree.
* `citation("mhn")` returns three `bibentry` objects: the package, the
  Sun et al. (2023) paper, and the Gao & Wang (2025) paper.
