# Skewness of the Modified Half-Normal Distribution

Computes the skewness \\\gamma_1 = E\[(X - \mu)^3\] / \sigma^3\\ for \\X
\sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_skewness(alpha, beta, gamma)
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
below 1e-9; the largest found, 3e-10, is for a shape between 2 and 3
under a positive tilt.

For a subnormal shape, below the smallest normal double 2.2e-308, the
moments are taken from forms in which the shape only multiplies (see the
Details of
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)). While the
mass at the peak is small against the mass by the origin (\\w\\ well
below 1, up to a tilt \\z = \gamma / \sqrt{\beta}\\ of about 53 at the
largest subnormal shape and 54.5 at the smallest), the skewness then
grows like \\1 / \sqrt{\alpha}\\ as the shape falls, by a factor that
depends on the tilt: 1.25 at \\z = 0\\ and \\2.8 \times 10^{-8}\\ at \\z
= 12\\. It turns negative where \\w\\ is about 1 (\\z\\ from about 53.3
to 54.7), and past the tilt where the forms hand over (55 to 56.3) it no
longer depends on \\\alpha\\.

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

[`mhn_kurtosis`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)

## Examples

``` r
mhn_skewness(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.6311107
mapply(mhn_skewness, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#> [1] 1.8814298 0.6311107 0.1579617
```
