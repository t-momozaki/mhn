# Density of the Modified Half-Normal Distribution

Computes the probability density function (or log-density) of the
Modified Half-Normal (MHN) distribution with parameters `alpha`, `beta`,
and `gamma`.

## Usage

``` r
dmhn(x, alpha = 1, beta = 1, gamma = 0, log = FALSE)
```

## Arguments

- x:

  Numeric vector of evaluation points.

- alpha:

  Shape parameter (\\\alpha \> 0\\). Scalar or numeric vector. Default:
  1.

- beta:

  Rate-type parameter (\\\beta \> 0\\), the coefficient of the Gaussian
  tail \\\exp(-\beta x^2)\\; larger values pull mass toward zero. Scalar
  or numeric vector. Default: 1.

- gamma:

  Tilt parameter (\\\gamma \in R\\), the coefficient of the exponential
  tilt \\\exp(\gamma x)\\. Scalar or numeric vector. Default: 0.

- log:

  Logical; if `TRUE`, log-density is returned. Default: `FALSE`. Read
  with [`as.logical()`](https://rdrr.io/r/base/logical.html), so that 1
  is `TRUE`; a value that gives `NA` is an error.

## Value

A numeric vector. The output length equals
`max(length(x), length(alpha), length(beta), length(gamma))`, or 0 if
any of them has length zero; each input is recycled to that length
following standard R recycling rules. For `x < 0`, the density is 0
(`-Inf` if `log = TRUE`). At `x = 0` it is `Inf` for `alpha < 1`, the
truncated normal density at 0 for `alpha = 1` (\\2\sqrt{\beta/\pi}\\ at
`gamma = 0`), and 0 for `alpha > 1`. As in
[`dgamma`](https://rdrr.io/r/stats/GammaDist.html), an element whose
parameters lie outside the parameter space (`alpha` or `beta` not
positive and finite, or `gamma` not finite) is `NaN`, with a warning,
and an element with an `NA` or `NaN` among its inputs is `NA` or `NaN`
(`NA` if any of them is `NA`), without one; the other elements are
unaffected. The arguments must be numeric (double, integer or logical):
a factor, a character string, a complex number, a list or `NULL` is an
error, as in [`dnorm`](https://rdrr.io/r/stats/Normal.html). As with
[`dnorm`](https://rdrr.io/r/stats/Normal.html), the attributes of `x`
(names, `dim`, `dimnames`) are kept when the result has the length of
`x`; a classed `x`, such as a Date, is taken as its numbers and gives a
plain vector. Unlike `dnorm`, attributes of a longer parameter vector
are not carried over.

## Details

The MHN density (Sun et al. 2023, Lemma 1a) is \$\$f(x \mid \alpha,
\beta, \gamma) = \frac{2 \beta^{\alpha/2} x^{\alpha-1} \exp(-\beta x^2 +
\gamma x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]} \quad (x \> 0)\$\$
where \\\Psi\[a, z\] = \sum\_{n \ge 0} \Gamma(a + n/2) z^n / n!\\ is the
Fox-Wright function \\{}\_1\Psi_1\[(a, 1/2); (1, 0); z\]\\ (Sun et al.
2023, Supplementary Eq. (1)).

The default parameters `alpha = 1, beta = 1, gamma = 0` correspond to
the half-normal distribution \\\mathrm{HN}(1/\sqrt{2})\\.

At \\\gamma = 0\\ and at \\\alpha = 1\\ exactly the density is computed
from the closed form of the law it reduces to; a shape or a tilt merely
close to those values is a different law, and is computed as in the
general case:

- \\\gamma = 0\\: sqrt-Gamma distribution, for \\\alpha \le 10^6\\; a
  larger shape is computed as in the general case.

- \\\alpha = 1\\: truncated normal distribution (Sun et al., 2023, Lemma
  6b), the normal law of mean \\\mu = \gamma/(2\beta)\\ and standard
  deviation \\\sigma = 1/\sqrt{2\beta}\\ restricted to \\(0, \infty)\\.
  However far 0 lies from the mean, the density holds to the Exponential
  and normal limits it tends to; where \\\gamma/(2\beta)\\ passes the
  largest double under a negative tilt, it is the Exponential density
  \\-\gamma e^{\gamma x}\\.

Otherwise the log density is computed from the formula above; for a
large shape or tilt, where the terms of that formula grow too large, the
law is measured instead from the peak of the density of \\\log X\\.

**Accuracy and range.** With \\z = \gamma/\sqrt{\beta}\\ the scale-free
tilt, for \\\alpha \le 10^{15}\\, \\\|z\| \le 10^{15}\\ and \\\beta\\ in
\\\[10^{-300}, 10^{300}\]\\ the log density is accurate to the working
tolerance `sqrt(.Machine$double.eps)`. Outside that range the law can be
narrower than the spacing of the doubles near its peak (for \\z\\ above
about \\10^{16}\\, or \\\alpha\\ above about \\10^{31}\\), and the
density is still that of the exact law at the double `x` given, almost
always 0 there. For a positive tilt whose mode \\\gamma/(2\beta)\\ lies
beyond the largest double, which needs a subnormal \\\beta\\,
`log = TRUE` returns `-Inf` where the log density is finite. Two cases
give `NA` with a warning instead: within 45 standard deviations of the
peak of a law whose peak cannot be placed to within \\10^{-11}\\ of its
standard deviation, which needs a law narrower than the spacing by a
factor of order \\10^{19}\\ or more and then happens for some parameters
only (as at \\\alpha = \beta = \gamma = 10^{200}\\); and where the whole
law lies beyond the largest double without a positive tilt to place it
(a shape near the largest double against a subnormal \\\beta\\). Where
\\z\\ itself overflows, the law is its limit: Gamma(\\\alpha\\, rate
\\\|\gamma\|\\) under a negative tilt and normal about its mode with
variance \\1/(2\beta)\\ under a positive one.

A subnormal shape, below the smallest normal double \\2.2 \times
10^{-308}\\, is computed apart, at \\\gamma = 0\\ too, from a form of
the law in which \\\alpha\\ only multiplies and \\X\\ lies below the
smallest positive double with probability \\1 / (1 + w)\\, where \\w\\,
defined in Appendix A.4 of the vignette cited below, grows with the
tilt. Below a tilt \\z\\ of about 52 (53 at the smallest shapes) that
probability is 1 in double precision:
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md) is 1 at
every positive \\x\\,
[`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md) is 0 for
every lower-tail probability below 1 given on the probability scale (the
upper tail, or `log.p = TRUE`, still gives the quantiles of the far
upper tail), and
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md) returns 0.
Above a tilt between about 55 and 56.3, as \\\alpha\\ falls, the mass by
the origin is at most \\2^{-64}\\, about \\5.4 \times 10^{-20}\\, of the
whole, and the density, the upper tail, the moments and the law of the
draws no longer depend on \\\alpha\\.

When any of `alpha`, `beta`, `gamma` is a vector, the density is
evaluated element-wise. The Fox-Wright \\\Psi\\ normalizing constant is
computed once for each position of the recycling cycle while that cycle
is at most 256 long, and otherwise once for each run of consecutive
elements that share the \\(\alpha, \beta, \gamma)\\ triple. So recycled
parameters, and parameter vectors in which equal triples are adjacent,
are much faster than triples that alternate along a long vector, and
than calling `dmhn` inside an R loop.

The numerical methods, and the reasons for the ranges above, are
described in Appendix A of
[`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md).

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.

## See also

[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
[`mhn_mode`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)

## Examples

``` r
x <- seq(0, 5, length.out = 100)
plot(x, dmhn(x, alpha = 2, beta = 1, gamma = 1), type = "l")


# Log-density
dmhn(1, alpha = 2, beta = 1, gamma = 1, log = TRUE)
#> [1] -0.3112403
```
