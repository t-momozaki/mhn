# Distribution Function of the Modified Half-Normal Distribution

Computes the cumulative distribution function (CDF) of the Modified
Half-Normal (MHN) distribution with parameters `alpha`, `beta`, and
`gamma`.

## Usage

``` r
pmhn(q, alpha = 1, beta = 1, gamma = 0, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- q:

  Numeric vector of quantiles.

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

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le q)\\,
  otherwise \\P(X \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.
  Default: `FALSE`. Both flags are read with
  [`as.logical()`](https://rdrr.io/r/base/logical.html), so that 1 is
  `TRUE`; a value that gives `NA` is an error.

## Value

A numeric vector. The output length equals
`max(length(q), length(alpha), length(beta), length(gamma))`, or 0 if
any of them has length zero; each input is recycled to that length
following standard R recycling rules. For `q <= 0` the CDF is 0; for
`q = Inf` it is 1. As in
[`pgamma`](https://rdrr.io/r/stats/GammaDist.html), an element whose
parameters lie outside the parameter space (`alpha` or `beta` not
positive and finite, or `gamma` not finite) is `NaN`, with a warning,
and an element with an `NA` or `NaN` among its inputs is `NA` or `NaN`
(`NA` if any of them is `NA`), without one; the other elements are
unaffected. The arguments must be numeric (double, integer or logical):
a factor, a character string, a complex number, a list or `NULL` is an
error, as in [`pnorm`](https://rdrr.io/r/stats/Normal.html). As with
[`pnorm`](https://rdrr.io/r/stats/Normal.html), the attributes of `q`
(names, `dim`, `dimnames`) are kept when the result has the length of
`q`; a classed `q`, such as a Date, is taken as its numbers and gives a
plain vector. Unlike `pnorm`, attributes of a longer parameter vector
are not carried over.

## Details

In the general case the lower tail is computed, where that is accurate
(see below), from the series representation (Sun et al. 2023, Lemma 1b)
\$\$F(x \mid \alpha, \beta, \gamma) = \frac{1}{\Psi\[\alpha/2,
\gamma/\sqrt{\beta}\]} \sum\_{i=0}^{\infty} \frac{z^i}{i!}\\
\Gamma(s_i)\\ P(s_i, \beta x^2)\$\$ where \\z = \gamma/\sqrt{\beta}\\,
\\s_i = (\alpha + i)/2\\, and \\P(s, y)\\ is the regularized lower
incomplete gamma function.

The default parameters `alpha = 1, beta = 1, gamma = 0` correspond to
the half-normal distribution \\\mathrm{HN}(1/\sqrt{2})\\.

At \\\gamma = 0\\ and at \\\alpha = 1\\ exactly the distribution
function is computed from the law it reduces to; a shape or a tilt
merely close to those values is a different law, and is computed as in
the general case:

- \\\gamma = 0\\: \\F(q) = P(\alpha/2, \beta q^2)\\, the regularized
  incomplete gamma function of the square-root-Gamma law (Sun et al.
  2023, Lemma 6a), for \\\alpha \le 10^6\\; a larger shape is computed
  as in the general case.

- \\\alpha = 1\\: the truncated-normal CDF (Sun et al. 2023, Lemma 6b).
  Near the truncation point the lower tail keeps its relative accuracy
  however small it is. Where \\-\gamma/\sqrt{2\beta}\\ overflows the law
  is Exponential with rate \\\|\gamma\|\\.

The series is used, where it is accurate, for \\\alpha \le 100\\ with
\\\|z\| \le 20\\, and for \\\alpha \< 1\\ up to \\\|z\|\\ of about 830
(for \\\|z\| \> 20\\ it costs milliseconds per value); elsewhere the
lower tail is computed by numerical integration of the density. The
upper tail \\P(X \> q)\\ is computed directly, by numerical integration
of the density over \\(q, \infty)\\, and not as \\1 - F\\, which would
round to 0 once it falls below about \\10^{-16}\\. Each tail is computed
directly on its own side of the peak of the density of \\\log X\\, which
lies at or above the median, and as one minus the other tail beyond it.
For a large shape or tilt, where the terms of these forms grow too
large, both tails are measured from that peak instead. With
`log.p = TRUE` both tails are formed on the log scale throughout, so
that \\\log F\\ is still resolved where \\F\\ rounds to 1 or is below
the smallest positive double, and a log probability close to 0, about
minus the other tail, keeps the relative accuracy of that other tail to
within about \\10^{-9}\\ (Appendix A.4 of the vignette cited below).

**Accuracy and range.** With \\z = \gamma/\sqrt{\beta}\\ as above, for
\\\alpha \le 10^{15}\\, \\\|z\| \le 10^{15}\\ and \\\beta\\ in
\\\[10^{-300}, 10^{300}\]\\ both tails are accurate to the working
tolerance `sqrt(.Machine$double.eps)` within 10 standard deviations of
the peak, and to about twice that beyond. Outside that range the law can
be narrower than the spacing of the doubles near its peak (for \\z\\
above about \\10^{16}\\, or \\\alpha\\ above about \\10^{31}\\); the
probabilities are then still those of the exact law at the doubles `q`
given, a step from 0 to 1 between two adjacent doubles, whose value at
the double nearest the peak is set by where the peak lies within its
spacing. Two cases give `NA` with a warning instead: within 45 standard
deviations of the peak of a law whose peak cannot be placed to within
\\10^{-11}\\ of its standard deviation, which needs a law narrower than
the spacing by a factor of order \\10^{19}\\ or more and then happens
for some parameters only (as at \\\alpha = \beta = \gamma = 10^{200}\\);
and where the whole law lies beyond the largest double without a
positive tilt to place it (a shape near the largest double against a
subnormal \\\beta\\). A third case lies at the other end of the range:
where \\\beta\\ and the tilt are both close to the largest double
(\\\beta\\ from about \\2 \times 10^{307}\\ with \\\gamma\\ near the
largest double, or from about \\5 \times 10^{307}\\ with \\\gamma\\ from
about \\10^{306}\\), at some points \\q\\ above the peak where \\\beta
q\\ is above about half the largest double, the log upper tail is `-Inf`
where it is finite (the upper tail itself is then 0, which is its value
in double precision), and next to the point where \\\beta q\\ overflows
both are `NA` with a warning. Where \\z\\ itself overflows, the law is
its limit: Gamma(\\\alpha\\, rate \\\|\gamma\|\\) under a negative tilt,
and normal about its mode with variance \\1/(2\beta)\\ under a positive
one. A subnormal shape, below the smallest normal double, is computed
apart (see
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)), and there
\\\log F\\ keeps its relative accuracy too; below a tilt \\z\\ of about
52 (53 at the smallest shapes) `pmhn` is then 1 at every positive `q`.

When any of `alpha`, `beta`, `gamma` is a vector, the CDF is evaluated
element-wise. The Fox-Wright \\\Psi\\ normalizing constant is computed
once for each position of the recycling cycle while that cycle is at
most 256 long, and otherwise once for each run of consecutive elements
that share the \\(\alpha, \beta, \gamma)\\ triple. So recycled
parameters, and parameter vectors in which equal triples are adjacent,
are much faster than triples that alternate along a long vector, and
than calling `pmhn` inside an R loop.

The numerical methods, and the reasons for the ranges above, are
described in Appendix A of
[`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md).

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)

## Examples

``` r
# Basic evaluation
pmhn(c(0.5, 1, 1.5), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.1129101 0.4338796 0.7614259

# Tail / log forms
pmhn(2, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#> [1] 0.06369619
pmhn(2, alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#> [1] -0.06581527

# Special case: gamma = 0 reduces to sqrt-Gamma
all.equal(pmhn(1.5, alpha = 2, beta = 1, gamma = 0),
          pgamma(1.5^2, shape = 1, rate = 1))
#> [1] TRUE
```
