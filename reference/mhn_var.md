# Variance of the Modified Half-Normal Distribution

Computes \\\mathrm{Var}(X)\\ for \\X \sim \mathrm{MHN}(\alpha, \beta,
\gamma)\\.

## Usage

``` r
mhn_var(alpha, beta, gamma)
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

Uses the formula (Sun et al. 2023, Lemma 2c) \$\$\mathrm{Var}(X) =
\frac{\alpha}{2\beta} + E(X)\left(\frac{\gamma}{2\beta} -
E(X)\right)\$\$ where it is accurate. Its terms cancel under a strong
tilt \\z = \gamma/\sqrt{\beta}\\ of either sign, where it also magnifies
any error in \\E(X)\\ by about \\z^2/2\\, and it magnifies that error by
about \\4\alpha\\ at a large \\\alpha\\. Wherever the combined error
would exceed about \\10^{-8}\\ of the result, the second central moment
is instead obtained by numerical integration of the unnormalised density
about its peak (or, where the density decreases from the origin, on a
logarithmic scale); so it is for \\z \> 12\\, where the mean itself is
obtained by numerical integration (see
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)),
for \\z \< -30\\, and for \\\alpha \ge 1000\\ unless \\\gamma = 0\\.

For \\\alpha \ge 1\\, the variance satisfies \\\mathrm{Var}(X) \le
1/(2\beta)\\ (Sun et al., 2023, Lemma 4c).

The function is not vectorised in its parameters: `alpha`, `beta` and
`gamma` must each be a single number. To evaluate several parameter
triples, use [`mapply`](https://rdrr.io/r/base/mapply.html), as in the
examples.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
(Lemmas 2c and 4c)

## See also

[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)

## Examples

``` r
mhn_var(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.2146018
mapply(mhn_var, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#> [1] 0.07828276 0.21460184 0.35382586
```
