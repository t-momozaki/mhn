# mhn: The Modified Half-Normal Distribution

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![R-CMD-check](https://github.com/t-momozaki/mhn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/t-momozaki/mhn/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/t-momozaki/mhn/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/t-momozaki/mhn/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

An R package providing density, distribution, quantile, and random
generation functions for the Modified Half-Normal (MHN) distribution,
together with deterministic helpers for its moments and mode.

## Overview

The MHN($\alpha$, $\beta$, $\gamma$) distribution has support on $(0, \infty)$
and density

$$f(x \mid \alpha, \beta, \gamma) \;\propto\;
  x^{\alpha - 1} \exp(-\beta x^2 + \gamma x), \qquad x > 0,$$

where $\alpha, \beta > 0$ and $\gamma \in \mathbb{R}$. It arises as a
full conditional in Bayesian MCMC samplers for a range of models (Sun
et al. 2023, Section 2; Gao & Wang 2025, Section 1): for example, the
extra scalar of the sandwich algorithm for binary regression with a $t$
link, the latent radius of the projected normal distribution, the
concentration parameter of the von Mises-Fisher distribution, the latent
scales of robust Bayesian graphical models, and shrinkage priors such as
the Bayesian elastic net. It also generalises a number of familiar
one-sided distributions:

| Constraint                          | Reduction                                              |
|------------------------------------|--------------------------------------------------------|
| $\gamma = 0$                       | $\sqrt{\mathrm{Gamma}}$: $X^2 \sim \mathrm{Gamma}(\alpha/2, \beta)$, rate $\beta$ |
| $\alpha = 1$                       | Truncated normal on $(0, \infty)$ with location $\gamma / (2\beta)$ and scale $1/\sqrt{2\beta}$ |
| $\alpha = 1,\; \gamma = 0$         | Half-normal with scale $1/\sqrt{2\beta}$              |
| $\beta \to 0^+,\; \gamma < 0$      | $\mathrm{Gamma}(\alpha, -\gamma)$ (limit)             |

The package implements the efficient samplers of **Sun, Kong & Pal
(2023)** (Algorithms 1 / 3) and the **Gao & Wang (2025)** Relaxed
Transformed Density Rejection (RTDR) method, whose acceptance Gao & Wang
bound below by 1/*e* where its contact points lie in their band (see
`?rmhn`), and dispatches between them automatically.

## Installation

```r
install.packages("mhn")
```

The development version is on GitHub:

```r
# install.packages("remotes")
remotes::install_github("t-momozaki/mhn")
```

## Quick start

```r
library(mhn)

# Density, CDF, quantile, random generation
dmhn(c(0.5, 1, 2), alpha = 2, beta = 1, gamma = 1)
pmhn(1.5,           alpha = 2, beta = 1, gamma = 1)
qmhn(0.95,          alpha = 2, beta = 1, gamma = 1)
rmhn(10,            alpha = 2, beta = 1, gamma = 1)

# Summary statistics
mhn_mean(2, 1, 1); mhn_var(2, 1, 1); mhn_mode(2, 1, 1)
```

## Function reference

### Distribution functions

| Function | Description                                                                    |
|----------|--------------------------------------------------------------------------------|
| `dmhn()` | Density (vectorised over `x` and parameters; supports `log = TRUE`)            |
| `pmhn()` | Cumulative distribution function (`lower.tail` and `log.p` à la `pgamma`)      |
| `qmhn()` | Quantile function (TOMS 748 root-finder in the general case)                   |
| `rmhn()` | Random generation; method dispatch via `method = c("auto", "rtdr", "sun")`     |

```r
dmhn(1.5, alpha = 2, beta = 1, gamma = 1, log = TRUE)
pmhn(1.5, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
qmhn(log(0.05), alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
rmhn(5, alpha = c(1, 2, 3), beta = 1, gamma = c(0, 1, -1))  # recycled
```

### Summary statistics

| Function          | Description                                                         |
|-------------------|---------------------------------------------------------------------|
| `mhn_mean()`      | $E(X) = \Psi[(\alpha+1)/2,\, z] \,/\, (\sqrt{\beta}\, \Psi[\alpha/2,\, z])$, with $z = \gamma/\sqrt{\beta}$; by quadrature where that ratio loses accuracy |
| `mhn_var()`       | Variance, Sun et al. (2023, Lemma 2c) where that form is well conditioned, quadrature elsewhere |
| `mhn_skewness()`  | Skewness $\gamma_1$, by quadrature of the central moments           |
| `mhn_kurtosis()`  | Excess kurtosis $\gamma_2$, by quadrature of the central moments    |
| `mhn_mode()`      | Mode (`NA` for α < 1 without an interior local maximum)             |

```r
mhn_mean(2, 1, 1)       # 1.1337...
mhn_skewness(2, 1, 1)   # positive (density right-skewed)
mhn_mode(0.5, 1, -1)    # NA: monotone-decreasing density
```

## Algorithm dispatch for `rmhn(method = "auto")`

The default method routes each parameter triple to one of the exact
samplers below.  Its ordinary decision rules are benchmarked, not
theoretical: they come from two paired timing races of the two samplers,
over α from 0.01 to 10^4 and from 1.5 to 10^14 (the second also timing
single draws within one vectorised call), and three runs of
`inst/benchmarks/auto_dispatch.R` (α ≥ 0.3), all on an Apple M1 Ultra (the
races are not shipped with the package), and
`vignette("theory")` gives the cost decomposition behind them.  The rules
change only the speed: every sampler draws the exact law, and the
crossings move from machine to machine.  Because the choice depends on the
number of draws, the default method's stream under one seed can change
with `n`; `method = "rtdr"` or `"sun"` keeps one stream (see `?rmhn`).  The
overrides at the end rest on the limits of double precision instead.

```
Closed-form shortcuts:
  α = 1, γ = 0    -> half-normal via |rnorm|
  γ = 0           -> sqrt(rgamma(...)), up to α = 2^60 (beyond, the
                     Normal proposal of Sun et al. (2023, Algorithm 1))
  α = 1           -> truncated normal

Region α < 1, γ > 0:
  Gao & Wang (2025) RTDR  (Sun et al. (2023) Algorithm 2 is
                           intentionally not implemented; RTDR was
                           the faster in every such cell that
                           Gao & Wang (2025, Tables 1-2) timed)

Region γ > 0, α > 1 (z = γ/√β):
  z ≤ 2√α and samples per setup ≥ N+
                          -> RTDR (lighter per-proposal cost: at mild
                             tilts Algorithm 1 draws a Gamma variate per
                             proposal, or accepts less often)
  otherwise               -> Sun et al. (2023, Algorithm 1) (closed-form
                             setup; it wins single draws, and from about
                             z = 2.5-3√α nearly every batch)
  where N+, for z ≤ √α / for √α < z ≤ 2√α, is
       5 / 25    for α < 10           50 / 100   for 10^3 ≤ α < 10^8
      25 / 50    for 10 ≤ α < 10^3   100 / 1000  for α ≥ 10^8
  and 5 / 25 where RTDR's setup starts next to the mode (see below).
  Above α = 2^59: Algorithm 1 at any batch size.

Region γ < 0:
  samples per setup ≥ N*  -> RTDR (lighter per-proposal cost)
  samples per setup <  N* -> Sun Algorithm 3 (lighter setup cost)
  where N* = 25 for α < 0.1, 5 for 0.1 ≤ α < 1, 1 for 1 < α < 10^3,
  10 for 10^3 ≤ α < 10^5, 25 for 10^5 ≤ α < 10^8 and 50 above, and
  N* = 1 wherever RTDR's setup starts next to the mode.  Above α = 1
  Algorithm 3 finds its matching point by bisection and Newton steps,
  and RTDR won every batch size timed from α = 1.5 to 100.  From
  α = 10^3, though, RTDR's setup walks its two contact points in from
  m/2 and m + max(1, m/2) (m the mode of √β·X), which costs more than
  Algorithm 3's setup; it starts them next to the mode only at very
  large shapes, from about 2.4e12 to 3.1e12 at moderate tilts, 2.3e13
  near a tilt of -α/e and 4.5e13 at any tilt, and for z < -1e10·√α.
  Overriding these at the limits of double precision, in this order and
  whatever the batch size:
    α > 2^59                      -> RTDR (R's rgamma, which Algorithm 3
                                     draws from, is too coarse there)
    -γ/√β > 1e150                 -> Sun Algorithm 3 (it proposes the
                                     Gamma(α, -γ) limit of the law)
    α < 2.5e-324·(-γ/√β)          -> Sun Algorithm 3 (RTDR cannot place
                                     the log-axis mode, which rounds
                                     to 0)
```

Choosing a specific sampler:

```r
rmhn(1000, 2, 1, 1, method = "rtdr")   # RTDR for the general case
rmhn(1000, 2, 1, 1, method = "sun")    # Sun for the general case
                                       # (errors for α < 1, γ > 0)
```

`method` selects the sampler used for the general case. The closed-form
special cases are taken first whatever it is set to: a negligible tilt draws
from the square root of a Gamma (up to α = 2^60), and α ≈ 1 from a truncated
normal. Both are exact samplers of their special-case law, which next to
γ = 0 or α = 1 is within a total variation distance of 8.6e-9 of the law
requested, and need next to no setup.

## Documentation

Full documentation, with a searchable function reference and rendered
vignettes, is published at <https://t-momozaki.github.io/mhn/>.

* `vignette("introduction", package = "mhn")` — a 5-minute tour of
  every exported function with worked examples.
* `vignette("theory", package = "mhn")` — the mathematics behind the
  family, its normalising constant and the samplers, and how
  `method = "auto"` chooses between them.
* `?dmhn`, `?pmhn`, `?qmhn`, `?rmhn` — full argument documentation
  including the recycling rules and the `method` argument.
* `citation("mhn")` — package + underlying papers.

## Citation

If you use this package in academic work, please cite both the package
and the methodology papers (`citation("mhn")` prints all four):

* Momozaki, T. (2026). *mhn: The Modified Half-Normal Distribution.*
  R package version 0.1.2.
* Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
  distribution: Properties and an efficient sampling scheme.
  *Communications in Statistics - Theory and Methods*, 52(5),
  1591–1613. <https://doi.org/10.1080/03610926.2021.1934700>
* Gao, F., & Wang, H.-B. (2025). Generating modified-half-normal random
  variates by a relaxed transformed density rejection method.
  *Communications in Statistics - Simulation and Computation*, 1–25.
  Advance online publication.
  <https://doi.org/10.1080/03610918.2025.2524551>

## References

* Robert, C. P. (1995). Simulation of truncated normal variables.
  *Statistics and Computing*, 5(2), 121–125.
  <https://doi.org/10.1007/BF00143942> (Used by the *α* = 1
  truncated-normal special case of `rmhn`.)

## License

MIT © 2026 Tomotaka Momozaki. See [`LICENSE`](LICENSE).
