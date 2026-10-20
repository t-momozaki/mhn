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

In the general case the lower tail is computed, over the part of the
parameter space given below, from the series representation \$\$F(x \mid
\alpha, \beta, \gamma) = \frac{1}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]}
\sum\_{i=0}^{\infty} \frac{z^i}{i!}\\ \Gamma(s_i)\\ P(s_i, \beta
x^2)\$\$ where \\z = \gamma/\sqrt{\beta}\\, \\s_i = (\alpha + i)/2\\,
and \\P(s, y)\\ is the regularized lower incomplete gamma function (Sun
et al. 2023, Lemma 1b: there each term
\\2\beta^{\alpha/2}\gamma^i\beta^{-(\alpha+i)/2}/(2\\ i!)\\ reduces to
\\z^i/i!\\, and \\\Gamma(s)\\ P(s, y) = \gamma(s, y)\\, where
\\\gamma(s, y)\\ is the lower incomplete gamma function used in the
paper). The infinite sum is truncated at an upper bound on the length
\\K = \max\\K_1, K_2\\\\ of Sun et al. (2023), Supplementary Lemma
10(d), from the geometric decay of the terms that Lemma 10(a, b)
establish, which makes the truncation residual bounded by the working
tolerance `sqrt(.Machine$double.eps)` divided by \\\Psi\\; for \\\gamma
\> 0\\, where every term is positive, it is truncated earlier where that
bound is below the rounding of the sum, so that the residual is below
the larger of that and \\2^{-53} F\\. The series is tried first for
\\\alpha \< 1\\ at every tilt up to \\z = 1000\\, and for \\\alpha \ge
1\\ only where \\\|z\| \le 20\\ and \\\alpha \le 100\\. There a
numerical integration of the density on \\\[0, q\]\\ replaces it where,
for \\\gamma \< 0\\, the rounding of the alternating terms, set against
the cancellation between them, could make its relative error exceed
`sqrt(.Machine$double.eps)`; where it would need more than \\10^6\\
terms, which happens for \\\|z\|\\ above about 830, and so only for
\\\alpha \< 1\\; and where \\\beta q^2\\ is below the smallest normal
double. Elsewhere, outside the ranges measured from the peak (see
below), the integration is used directly, being the faster of the two
there. It is Gauss-Kronrod for \\\alpha \ge 1\\, and for \\\alpha \<
1\\, where the density is unbounded at the origin, a closed form next to
the origin, then tanh-sinh up to \\q\\ or up to the interior peak that a
positive tilt can give the density, whichever comes first, and
Gauss-Kronrod past that peak. For \\z \> 1000\\ the integral and
\\\Psi\\ are both taken with \\\gamma^2/(4\beta)\\ removed from their
logarithms, which are otherwise too large for their difference to be
accurate.

The upper tail \\P(X \> q)\\ is computed directly, by Gauss-Kronrod
quadrature of the density over \\(q, \infty)\\ on a logarithmic axis,
and not as \\1 - F\\, which would round to 0 once it falls below about
\\10^{-16}\\. Each tail is computed directly on its own side of the peak
of the density of \\\log X\\, the positive root of \\2\beta x^2 - \gamma
x - \alpha = 0\\, which lies at or above the median; above it the lower
tail is returned as one minus the upper tail. With `log.p = TRUE` both
tails are formed on the log scale throughout, so that \\\log F\\ is
still resolved where \\F\\ rounds to 1 or is below the smallest positive
double. At or below the median \\\log P(X \> q)\\ is about \\-F\\, and
the logarithm of the quadrature carries an absolute error of about
\\\epsilon\\ times the size of the terms it is assembled from,
\\\|\log\Psi\| + (\alpha/2)\|\log\beta\| + m(1 + \|\log v^\*\|)\\ (\\m\\
and \\v^\*\\ as in **Range**), with \\\gamma^2/(4\beta)\\ removed from
\\\log\Psi\\ for \\z \> 1000\\, as in the quadratures. Wherever that
could pass \\10^{-9}\\ of \\F\\, it is formed instead as \\\log(1 - F)\\
from \\\log F\\, which keeps the relative accuracy of \\F\\: for
parameters of order one only below \\F\\ of about \\10^{-5}\\. For
\\\alpha \< 1\\ with \\\|z\| \> 20\\, where the series costs
milliseconds, \\F\\ is taken for this from the numerical integration
wherever that agrees with the upper tail to within the error of the
latter. The lower tail has the matching switch. Below that peak \\F\\
can still exceed one half, and for a small shape it is close to 1 over
most of the range below the peak, one less a tail proportional to
\\\alpha\\. There \\\log F\\, about \\-P(X \> q)\\, is assembled from
terms of the size of \\\log\Psi\\, \\(\alpha/2)\log\beta\\ and
\\(\alpha - 1)\log q\\ that cancel, and carries an absolute error of
about \\\epsilon\\ times their size. Wherever that could pass
\\10^{-9}\\ of \\\log F\\ and the error of the upper tail, carried into
\\\log F\\, is estimated at less than half as much, it is formed instead
as \\\log(1 - P(X \> q))\\ from the upper tail, which keeps the relative
accuracy of \\P(X \> q)\\. Without that switch the error would pass
\\10^{-9}\\ of \\\log F\\ for a shape below about \\10^{-6}\\ under a
small tilt, and below a smaller shape under a large one. A subnormal
shape, below the smallest normal double, is served apart (see
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)), and there
\\\log F\\ keeps its relative accuracy too.

At \\\gamma = 0\\ and at \\\alpha = 1\\ exactly the distribution
function is computed from the law it reduces to; a shape or a tilt
merely close to those values is a different law, and is computed as in
the general case:

- \\\gamma = 0\\: \\F(q) = P(\alpha/2, \beta q^2)\\, the regularized
  incomplete gamma function of the square-root-Gamma law (Sun et al.
  2023, Lemma 6a), by `pgamma`, or as its leading term \\(\beta
  q^2)^{\alpha/2} / \Gamma(\alpha/2 + 1)\\, exact to double precision,
  where \\\beta q^2\\ is below the smallest normal double; for \\\alpha
  \le 10^6\\. Above it the Gamma law is so narrow that rounding \\\beta
  q^2\\ moves it by a visible part of its width (up to about \\10^{-6}\\
  of \\F\\ at \\\alpha = 10^{20}\\), and \\\gamma = 0\\ is measured from
  the peak as below.

- \\\alpha = 1\\: the truncated-normal CDF (Sun et al. 2023, Lemma 6b),
  by `pnorm`, or through the scaled complementary error function where
  the truncation point 0 lies three or more standard deviations above
  the mean \\\gamma/(2\beta)\\ of the untruncated normal,
  \\-\gamma/\sqrt{2\beta} \ge 3\\. Near the truncation point the lower
  tail is summed from a short series, so that it keeps its relative
  accuracy however small it is.

For a shape \\\alpha\\ or a scale-free tilt \\z = \gamma/\sqrt{\beta}\\
above \\10^6\\, and wherever the terms of the ordinary forms grow too
large (see **Range**), both tails are measured from the peak \\v^\*\\ of
the density of \\\log X\\, the positive root of \\2\beta x^2 - \gamma
x - \alpha = 0\\, in the variable \\w = \log(x/v^\*)\\: the tail on the
side of \\\log(q/v^\*)\\ away from the peak is the integral of
\\\exp(-\alpha(e^w - 1 - w) - \beta v^{\*2}(e^w - 1)^2)\\ over that
side, taken by Gauss-Kronrod quadrature in offsets from `q`, over the
same integral on the whole line, and the other tail is its complement.
Nothing of the size of \\\alpha \log x\\ or \\\gamma^2/(4\beta)\\ is
formed, whose rounding would otherwise leave the tails \\10^{-7}\\ off
at \\\alpha = 10^8\\ and farther off as \\z\\ grows, and no abscissa is
rounded to the doubles next to the peak, which are a fraction
\\z\\2^{-53}\\ of the law's standard deviation apart. \\v^\*\\ is held
as a double next to it plus its offset from that double.

**Range.** The ordinary forms are used only while the larger \\m\\ of
\\\alpha\\ and \\z\\ is at most \\10^6\\ and \\m (1 + \|\log v^\*\|)\\,
the size of their largest terms (\\\alpha \log x\\, \\(\alpha/2)
\log\beta\\, and the rounding of quadrature abscissae placed in \\\log
x\\), is at most \\10^7\\: at \\\beta = 1\\ up to a shape of \\10^6\\ or
a tilt of about \\7.3 \times 10^5\\, and at \\\beta = 10^{\pm 300}\\ up
to about \\3 \times 10^4\\ for either. That keeps their error within 10
standard deviations of the peak below \\10^{-8}\\; measured from the
peak it is about \\10^{-14}\\. Farther out in the upper tail the error
of the ordinary forms grows with the distance from the peak, and close
to that bound it reaches about \\3 \times 10^{-8}\\ at 36 standard
deviations, where the tail is of order \\10^{-284}\\. So for \\\alpha
\le 10^{15}\\, \\\|z\| \le 10^{15}\\ and \\\beta\\ in \\\[10^{-300},
10^{300}\]\\ both tails are accurate to the working tolerance
`sqrt(.Machine$double.eps)` within 10 standard deviations of the peak,
and to about twice that beyond. Outside that range the law can be
narrower than the spacing of the doubles near its peak (for \\z\\ above
about \\10^{16}\\, or \\\alpha\\ above about \\10^{31}\\); the
probabilities are then still those of the exact law at the doubles `q`
given, a step from 0 to 1 between two adjacent doubles, whose value at
the double nearest the peak is set by where the peak lies within its
spacing. Two cases give `NA` with a warning instead: within 45 standard
deviations of the peak of a law narrower than the spacing by a factor of
order \\10^{19}\\ (\\z\\ beyond about \\10^{35}\\, \\\alpha\\ beyond
about \\10^{70}\\), where the position of the peak cannot be resolved
finely enough; and where the whole law lies beyond the largest double
without a positive tilt to place it (a shape near the largest double
against a subnormal \\\beta\\). A third case lies at the other end of
the range: where \\\beta\\ and the tilt are both close to the largest
double (\\\beta\\ from about \\2 \times 10^{307}\\ with \\\gamma\\ near
the largest double, or from about \\5 \times 10^{307}\\ with \\\gamma\\
from about \\10^{306}\\), at points \\q\\ above the peak where \\\beta
q\\ passes about half the largest double, the upper tail is `NA` with a
warning or, with `log.p = TRUE`, `-Inf` where the log probability is
finite. Where \\z\\ itself overflows, the law is its limit:
Gamma(\\\alpha\\, rate \\\|\gamma\|\\) under a negative tilt, and normal
about its mode with variance \\1/(2\beta)\\ under a positive one. For
\\\alpha = 1\\ with the truncation point \\z_0 = -\gamma/\sqrt{2\beta}\\
far beyond 3, the same forms hold up to \\z_0 = \infty\\, where the law
is Exponential with rate \\\|\gamma\|\\; below \\-10^3\\ the
standardised offset is formed from \\2\beta q - \gamma\\ in one rounding
(as \\\beta q - \gamma/2\\, so also where \\2\beta\\ overflows), not
from a rounded mean \\\gamma/(2\beta)\\; the lower tail within
\\\|z_0\|/2\\ standard deviations of 0, where that offset is resolved
only to about \\\epsilon \|z_0\|\\ of one, is taken from the ratio of
the two normal lower tails written with the exact distance
\\q\sqrt{2\beta}\\ from the truncation point.

When any of `alpha`, `beta`, `gamma` is a vector, the CDF is evaluated
element-wise. The Fox-Wright \\\Psi\\ normalizing constant is computed
once for each position of the recycling cycle while that cycle is at
most 256 long, and otherwise once for each run of consecutive elements
that share the \\(\alpha, \beta, \gamma)\\ triple. So recycled
parameters, and parameter vectors in which equal triples are adjacent,
are much faster than triples that alternate along a long vector, and
than calling `pmhn` inside an R loop.

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
