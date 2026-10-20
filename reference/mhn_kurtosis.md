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

The central moments are obtained by numerical integration against the
unnormalised density, in a variable centred on its peak (or, where the
density decreases from the origin, on a logarithmic scale), so no
normalising constant is needed. Building them from raw moments by the
Lemma 2b recurrence from the mean of Lemma 2a is exact algebra but
cancels once the standard deviation is small next to the mean, which is
what a large tilt produces; that expansion is used only where the
integration cannot be set up.

Against independent quadratures, for shapes from 1 to 60 and tilts
\\\gamma / \sqrt{\beta}\\ from -20 to 12, the absolute error stayed
below 1e-9; the largest found, 6e-10, is for a shape between 2 and 3
under a positive tilt.

For a subnormal shape, below the smallest normal double 2.2e-308, the
moments are taken from forms in which the shape only multiplies (see the
Details of
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)). While the
mass at the peak is small against the mass by the origin (\\w\\ well
below 1, up to a tilt \\z = \gamma / \sqrt{\beta}\\ of about 53 at the
largest subnormal shape and 54.5 at the smallest), the excess kurtosis
then grows like \\1 / \alpha\\ as the shape falls, by a factor that
depends on the tilt: 2 at \\z = 0\\ and \\8.2 \times 10^{-16}\\ at \\z =
12\\. It is `Inf` where it passes the largest double: below a shape of
about 1.1e-308 at \\z = 0\\, at every subnormal shape for \\z = -30\\,
and at none for \\z = 12\\. Where \\w\\ reaches about 1 (\\z\\ from
about 53.3 to 54.7) the peak takes over, and past the tilt where the
forms hand over (55 to 56.3) the kurtosis no longer depends on
\\\alpha\\.

The function is not vectorised in its parameters: `alpha`, `beta` and
`gamma` must each be a single number. To evaluate several parameter
triples, use [`mapply`](https://rdrr.io/r/base/mapply.html), as in the
examples.

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
