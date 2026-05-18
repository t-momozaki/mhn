# Mode of the Modified Half-Normal Distribution

Computes the mode (most probable value) of the MHN distribution.

## Usage

``` r
mhn_mode(alpha, beta, gamma)
```

## Arguments

- alpha:

  Shape parameter (\\\alpha \> 0\\).

- beta:

  Scale parameter (\\\beta \> 0\\).

- gamma:

  Location parameter (\\\gamma \in R\\).

## Value

A numeric scalar. Returns `NA` when no interior mode exists (density is
monotonically decreasing on \\(0, \infty)\\).

## Details

The mode depends on \\\alpha\\:

- \\\alpha \> 1\\:

  \\(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) / (4\beta)\\ (Sun et
  al., 2023, Lemma 3b).

- \\\alpha = 1\\:

  \\\max(0, \gamma / (2\beta))\\, obtained as the mode of the truncated
  normal \\\mathrm{TN}(\gamma/(2\beta), 1/\sqrt{2\beta}, 0, \infty)\\
  that the MHN reduces to in this case (Sun et al., 2023, Lemma 6b).

- \\0 \< \alpha \< 1\\:

  An interior mode exists only when \\\gamma \> 0\\ and \\\alpha \geq
  1 - \gamma^2 / (8\beta)\\ (Sun et al., 2023, Lemma 3c); otherwise the
  density is monotonically decreasing (Sun et al., 2023, Lemma 3d) and
  `NA` is returned.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.
(Lemma 3b–d, Lemma 6b)

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)

## Examples

``` r
mhn_mode(alpha = 2, beta = 1, gamma = 1)
#> [1] 1
mhn_mode(alpha = 1, beta = 1, gamma = 2)
#> [1] 1
mhn_mode(alpha = 0.5, beta = 1, gamma = -1)  # NA
#> [1] NA
```
