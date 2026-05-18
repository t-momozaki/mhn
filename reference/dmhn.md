# Density of the Modified Half-Normal Distribution

Computes the probability density function (or log-density) of the
Modified Half-Normal (MHN) distribution with parameters `alpha`, `beta`,
and `gamma`.

## Usage

``` r
dmhn(x, alpha = 1, beta = 1, gamma = 0, log = FALSE)
```

## Arguments

- x:

  Numeric vector of evaluation points.

- alpha:

  Shape parameter (\\\alpha \> 0\\). Scalar or numeric vector. Default:
  1.

- beta:

  Scale parameter (\\\beta \> 0\\). Scalar or numeric vector. Default:
  1.

- gamma:

  Location parameter (\\\gamma \in R\\). Scalar or numeric vector.
  Default: 0.

- log:

  Logical; if `TRUE`, log-density is returned. Default: `FALSE`.

## Value

A numeric vector. The output length equals
`max(length(x), length(alpha), length(beta), length(gamma))`; each input
is recycled to that length following standard R recycling rules. For
`x < 0`, the density is 0 (`-Inf` if `log = TRUE`).

## Details

The MHN density is \$\$f(x \mid \alpha, \beta, \gamma) = \frac{2
\beta^{\alpha/2} x^{\alpha-1} \exp(-\beta x^2 + \gamma
x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]} \quad (x \> 0)\$\$ where
\\\Psi\[a, z\]\\ is the Fox-Wright Psi function (Sun et al., 2023, Lemma
1a).

The default parameters `alpha = 1, beta = 1, gamma = 0` correspond to
the half-normal distribution \\\mathrm{HN}(1/\sqrt{2})\\.

Special cases are detected and dispatched to closed-form solutions:

- \\\gamma = 0\\: sqrt-Gamma distribution

- \\\alpha = 1\\: truncated normal distribution

Computation is performed in log-space to avoid numerical
underflow/overflow.

When any of `alpha`, `beta`, `gamma` is a vector, the density is
evaluated element-wise. The Fox-Wright \\\Psi\\ normalizing constant is
recomputed only when consecutive elements present a different \\(\alpha,
\beta, \gamma)\\ triple, so passing grouped parameters is significantly
faster than calling `dmhn` inside an R loop.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.

## See also

[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
[`mhn_mode`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)

## Examples

``` r
x <- seq(0, 5, length.out = 100)
plot(x, dmhn(x, alpha = 2, beta = 1, gamma = 1), type = "l")


# Log-density
dmhn(1, alpha = 2, beta = 1, gamma = 1, log = TRUE)
#> [1] -0.3112403
```
