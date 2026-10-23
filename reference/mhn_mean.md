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
\gamma/\sqrt{\beta}.\$\$ At \\\gamma = 0\\ the ratio is
\\\Gamma((\alpha+1)/2) / (\sqrt{\beta}\\ \Gamma(\alpha/2))\\ (Sun et al.
2023, Lemma 6a). For every \\\alpha \> 0\\ the mean is at most
\\(\gamma + \sqrt{\gamma^2 + 8\alpha\beta})/(4\beta)\\ (Sun et al. 2023,
Lemma 4a and the remark after it).

The ratio is used where it can be evaluated accurately. Where it cannot,
mainly under a strong tilt or at a very large \\\alpha\\, the mean is
instead obtained by numerical integration of the unnormalised density,
which needs no \\\Psi\\ value.

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
