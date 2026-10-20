# Mode of the Modified Half-Normal Distribution

Computes the mode of the MHN distribution. For \\\alpha \ge 1\\ this is
the maximiser of the density. For \\0 \< \alpha \< 1\\ the density is
unbounded as \\x \to 0^+\\, so a non-`NA` result is an interior local
maximum, not the most probable value.

## Usage

``` r
mhn_mode(alpha, beta, gamma)
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

A numeric scalar. For \\0 \< \alpha \< 1\\ with no interior local
maximum (density strictly decreasing on \\(0, \infty)\\) the result is
`NA`; at \\\alpha = 1\\ with \\\gamma \le 0\\ the mode sits on the
boundary and 0 is returned. See Details.

## Details

The mode depends on \\\alpha\\:

- \\\alpha \> 1\\:

  \\(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) / (4\beta)\\ (Sun et
  al. 2023, Lemma 3b).

- \\\alpha = 1\\:

  \\\max(0, \gamma / (2\beta))\\, obtained as the mode of the truncated
  normal \\\mathrm{TN}(\gamma/(2\beta), 1/\sqrt{2\beta}, 0, \infty)\\
  that the MHN reduces to in this case (Sun et al. 2023, Lemma 6b). For
  \\\gamma \le 0\\ this is the boundary point 0. Sun et al. (2023, Lemma
  3d) state that no mode exists for \\\gamma \< 0\\ and \\\alpha \le
  1\\; at \\\alpha = 1\\ the density is bounded and decreases from its
  supremum at \\x = 0\\, and that point is returned. Only \\\alpha = 1\\
  itself takes this branch. The Lemma 3b value tends to it as \\\alpha
  \to 1\\ from above, and for \\\gamma \> 0\\ so does the Lemma 3c value
  from below; for \\\gamma \le 0\\ and \\\alpha \< 1\\ there is no mode
  (`NA`).

- \\0 \< \alpha \< 1\\:

  An interior local maximum exists only when \\\gamma \> 0\\ and
  \\\alpha \> 1 - \gamma^2 / (8\beta)\\, and it is then located at
  \\(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) / (4\beta)\\ (Sun et
  al. 2023, Lemma 3c). Lemma 3c states the condition with \\\ge\\; the
  strict inequality is used here because at \\\alpha = 1 -
  \gamma^2/(8\beta)\\ the two stationary points coincide at
  \\\gamma/(4\beta)\\, where \\\frac{d}{dx} \log f(x) =
  -\frac{2\beta}{x} \left(x - \frac{\gamma}{4\beta}\right)^2 \le 0\\, so
  the density is still strictly decreasing and has an inflection point,
  not a maximum, there. In every other case the density is strictly
  decreasing on \\(0, \infty)\\ (Sun et al. 2023, Lemma 3d, for \\\gamma
  \ne 0\\; for \\\gamma = 0\\ the derivative \\-(2\beta x^2 + 1 -
  \alpha)/x\\ of \\\log f\\ is negative) and `NA` is returned.

The function is not vectorised in its parameters: `alpha`, `beta` and
`gamma` must each be a single number. To evaluate several parameter
triples, use [`mapply`](https://rdrr.io/r/base/mapply.html), as in the
examples.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
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
mhn_mode(alpha = 0.5, beta = 1, gamma = 5)   # a local maximum only
#> [1] 2.395644
mapply(mhn_mode, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#> [1]        NA 0.7071068 2.3507811
```
