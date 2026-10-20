# Mean of the Modified Half-Normal Distribution

Computes \\E(X)\\ for \\X \sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_mean(alpha, beta, gamma)
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

By Sun et al. (2023, Lemma 2a) the mean is a ratio of Fox-Wright
\\\Psi\\ functions, \$\$E(X) = \frac{\Psi\[(\alpha+1)/2,\\
z\]}{\sqrt{\beta}\\ \Psi\[\alpha/2,\\ z\]}, \qquad z =
\gamma/\sqrt{\beta}.\$\$ The ratio is used for \\z \le 12\\ as long as
the two logarithms it is formed from are small enough to leave it
accurate. For a larger positive tilt, or a very large \\\alpha\\, both
logarithms are large while their difference is not, and the mean is
instead obtained by numerical integration of the unnormalised density
about its peak, which needs no \\\Psi\\ value. (For \\z \> 12\\ and a
shape so small that most of the mass lies by the origin, that
integration holds the mean as a small difference; where it finds the
mean below a thousandth of the peak, the ratio is used instead.) Under a
strong negative tilt, \\z \< -30\\, the ratio serves only for \\\alpha
\< 6\\, and there it is capped at the upper bound \\(\gamma +
\sqrt{\gamma^2 + 8\alpha\beta})/(4\beta)\\ of Sun et al. (2023, Lemma
4a), which holds for every \\\alpha \> 0\\ by the remark after the
lemma; from \\\alpha = 6\\ on the mean is integrated there. At \\\gamma
= 0\\ the ratio is \\\Gamma((\alpha+1)/2) / (\sqrt{\beta}\\
\Gamma(\alpha/2))\\ (Sun et al. 2023, Lemma 6a), which for \\\alpha \ge
1\\ is computed as \\\sqrt{\pi} / (\sqrt{\beta}\\ B(\alpha/2, 1/2))\\,
with \\B\\ the beta function, whose logarithm R's `lbeta` forms without
the cancellation between the two log-gamma values that loses digits at a
large \\\alpha\\.

The function is not vectorised in its parameters: `alpha`, `beta` and
`gamma` must each be a single number. To evaluate several parameter
triples, use [`mapply`](https://rdrr.io/r/base/mapply.html), as in the
examples.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
(Lemmas 2a, 4a and 6a)

## See also

[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)

## Examples

``` r
mhn_mean(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.8862269
mapply(mhn_mean, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#> [1] 0.2339600 0.8862269 2.3958050
```
