# Quantile Function of the Modified Half-Normal Distribution

Computes the quantile (inverse cumulative) function of the Modified
Half-Normal (MHN) distribution with parameters `alpha`, `beta`, and
`gamma`.

## Usage

``` r
qmhn(p, alpha = 1, beta = 1, gamma = 0, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  Numeric vector of probabilities.

- alpha:

  Shape parameter (\\\alpha \> 0\\). Scalar or numeric vector. Default:
  1.

- beta:

  Scale parameter (\\\beta \> 0\\). Scalar or numeric vector. Default:
  1.

- gamma:

  Location parameter (\\\gamma \in R\\). Scalar or numeric vector.
  Default: 0.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le q)\\,
  otherwise \\P(X \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are provided on the log scale.
  Default: `FALSE`.

## Value

A numeric vector. The output length equals
`max(length(p), length(alpha), length(beta), length(gamma))`; each input
is recycled to that length following standard R recycling rules.
`qmhn(0) = 0` and `qmhn(1) = Inf`. Probabilities outside \\\[0, 1\]\\
yield `NaN`.

## Details

For the general case, \\q = F^{-1}(p)\\ is obtained by a TOMS 748
root-finder applied to the series CDF (Sun et al., 2023, Lemma 1b). The
initial bracket is \\\[\sqrt{\epsilon},\\ E(X) + 8
\sqrt{\mathrm{Var}(X)}\]\\ and is doubled on the right (up to 30 times)
until it brackets the target probability.

Special cases are detected and dispatched to standard R primitives:

- \\\gamma = 0\\: `sqrt(qgamma(p, alpha/2, scale = 1/beta))`

- \\\alpha = 1\\: truncated-normal inverse via `qnorm`

When any of `alpha`, `beta`, `gamma` is a vector, the quantile is
evaluated element-wise. The Fox-Wright \\\Psi\\ normalizing constant and
moments \\E(X)\\, \\\mathrm{Var}(X)\\ (used to size the root-finder
bracket) are recomputed only when consecutive elements present a
different \\(\alpha, \beta, \gamma)\\ triple.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)

## Examples

``` r
# Basic evaluation
qmhn(c(0.1, 0.5, 0.9), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.4717434 1.0906276 1.8480472

# Round-trip: F(F^-1(p)) ~ p
p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
all.equal(pmhn(qmhn(p, alpha = 2, beta = 1, gamma = 1),
               alpha = 2, beta = 1, gamma = 1),
          p, tolerance = 1e-6)
#> [1] TRUE

# Tail / log forms
qmhn(0.95, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#> [1] 0.3394014
qmhn(log(0.05), alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#> [1] 0.3394014
```
