# Mean of the Modified Half-Normal Distribution

Computes \\E(X)\\ for \\X \sim \mathrm{MHN}(\alpha, \beta, \gamma)\\.

## Usage

``` r
mhn_mean(alpha, beta, gamma)
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

The mean is computed as a ratio of Fox-Wright Psi functions: \$\$E(X) =
\frac{\Psi\[(\alpha+1)/2,\\ \gamma/\sqrt{\beta}\]}{ \sqrt{\beta}\\
\Psi\[\alpha/2,\\ \gamma/\sqrt{\beta}\]}\$\$

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.
(Lemma 2a)

## See also

[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)

## Examples

``` r
mhn_mean(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.8862269
```
