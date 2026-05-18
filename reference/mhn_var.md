# Variance of the Modified Half-Normal Distribution

Computes \\\mathrm{Var}(X)\\ for \\X \sim \mathrm{MHN}(\alpha, \beta,
\gamma)\\.

## Usage

``` r
mhn_var(alpha, beta, gamma)
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

Uses the formula (Sun et al., 2023, Lemma 2c): \$\$\mathrm{Var}(X) =
\frac{\alpha}{2\beta} + E(X)\left(\frac{\gamma}{2\beta} -
E(X)\right)\$\$

For \\\alpha \geq 1\\, the variance satisfies \\\mathrm{Var}(X) \leq
1/(2\beta)\\ (Sun et al., 2023, Lemma 4c).

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.
(Lemma 2c)

## See also

[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)

## Examples

``` r
mhn_var(alpha = 2, beta = 1, gamma = 0)
#> [1] 0.2146018
```
