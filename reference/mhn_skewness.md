# Skewness of the Modified Half-Normal Distribution

Computes the skewness \\\gamma_1 = E\[(X - \mu)^3\] / \sigma^3\\ for \\X
\sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_skewness(alpha, beta, gamma)
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
moments up to third order, then converts to central moments.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.
(Lemma 2b)

## See also

[`mhn_kurtosis`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)

## Examples

``` r
mhn_skewness(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.6311107
```
