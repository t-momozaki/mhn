# Excess Kurtosis of the Modified Half-Normal Distribution

Computes the excess kurtosis \\\gamma_2 = E\[(X - \mu)^4\] / \sigma^4 -
3\\ for \\X \sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_kurtosis(alpha, beta, gamma)
```

## Arguments

- alpha:

  Shape parameter (\\\alpha \> 0\\).

- beta:

  Scale parameter (\\\beta \> 0\\).

- gamma:

  Location parameter (\\\gamma \in R\\).

## Value

A numeric scalar.

## Details

Uses the moment recurrence (Sun et al., 2023, Lemma 2b) to compute raw
moments up to fourth order, then converts to central moments.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.
(Lemma 2b)

## See also

[`mhn_skewness`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)

## Examples

``` r
mhn_kurtosis(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.2450893
```
