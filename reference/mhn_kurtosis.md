# Excess Kurtosis of the Modified Half-Normal Distribution

Computes the excess kurtosis \\\gamma_2 = E\[(X - \mu)^4\] / \sigma^4 -
3\\ for \\X \sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_kurtosis(alpha, beta, gamma)
```

## Arguments

- alpha:

  Shape parameter (\\\alpha \> 0\\), a single number.

- beta:

  Rate-type parameter (\\\beta \> 0\\), the coefficient of the Gaussian
  tail \\\exp(-\beta x^2)\\; larger values pull mass toward zero. A
  single number.

- gamma:

  Tilt parameter (\\\gamma \in R\\), the coefficient of the exponential
  tilt \\\exp(\gamma x)\\. A single number.

## Value

A numeric scalar.

## Details

The excess kurtosis is formed from the second and fourth central
moments. In exact arithmetic these follow from the mean of Sun et al.
(2023, Lemma 2a) and the recurrence for the raw moments of their Lemma
2b; the function obtains them instead, wherever it can, by numerical
integration of the unnormalised density, which needs no normalising
constant and stays accurate under a large tilt.

Against independent quadratures, for shapes from 1 to 60 and tilts
\\\gamma / \sqrt{\beta}\\ from -20 to 12, the absolute error stayed
below 1e-9; the largest found, 6e-10, is for a shape between 2 and 3
under a positive tilt.

For a subnormal shape, below the smallest normal double 2.2e-308, the
excess kurtosis depends on the tilt \\z = \gamma / \sqrt{\beta}\\ as
follows. While the mass at the peak of the density is small against the
mass by the origin (up to \\z\\ of about 53 at the largest subnormal
shape and 54.5 at the smallest), it grows like \\1 / \alpha\\ as the
shape falls, by a factor that depends on the tilt: 2 at \\z = 0\\ and
\\8.2 \times 10^{-16}\\ at \\z = 12\\. It is `Inf` where it passes the
largest double: below a shape of about 1.1e-308 at \\z = 0\\, at every
subnormal shape for \\z = -30\\, and at none for \\z = 12\\. Where the
two masses are about equal (\\z\\ from about 53.3 to 54.7) the peak
takes over, and above a tilt of 55 to 56.3 the kurtosis no longer
depends on \\\alpha\\.

The function is not vectorised in its parameters: `alpha`, `beta` and
`gamma` must each be a single number. To evaluate several parameter
triples, use [`mapply`](https://rdrr.io/r/base/mapply.html), as in the
examples.

The numerical methods are described in Appendix A of
[`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md).

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
(Lemmas 2a and 2b)

## See also

[`mhn_skewness`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)

## Examples

``` r
mhn_kurtosis(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.2450893
mapply(mhn_kurtosis, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#> [1]  4.26550287  0.24508930 -0.04785938
```
