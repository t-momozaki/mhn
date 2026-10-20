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

- \\\gamma = 0\\: sqrt-Gamma distribution, for \\\alpha \le 10^6\\.
  Above it the Gamma law of \\\beta X^2\\ is so narrow that rounding
  \\\beta x^2\\ moves it by a visible part of its width (up to about
  \\10^{-5}\\ of the log density at \\\alpha = 10^{20}\\), and \\\gamma
  = 0\\ is measured from the peak as below.

- \\\alpha = 1\\: truncated normal distribution (Sun et al., 2023, Lemma
  6b), the normal law of mean \\\mu = \gamma/(2\beta)\\ and standard
  deviation \\\sigma = 1/\sqrt{2\beta}\\ restricted to \\(0, \infty)\\:
  `dnorm` divided by its upper tail at 0 where 0 lies within 1000
  standard deviations of the mean (\\\|\gamma\|/\sqrt{2\beta} \le
  1000\\). Farther out the same ratio is written with the normal hazard
  where 0 lies above the mean, so that two logarithms of about
  \\-\gamma^2/(4\beta)\\ are not subtracted, and from the standardised
  point \\(x - \mu)/\sigma\\, formed directly, where it lies below.
  Where \\\gamma/(2\beta)\\ passes the largest double under a negative
  tilt, the density is the exponential one, \\-\gamma e^{\gamma x}\\.

A subnormal shape, below the smallest normal double \\2.2 \times
10^{-308}\\, is computed apart, at \\\gamma = 0\\ too. There half of
\\\alpha\\ is not always a double and \\1/\alpha\\ can overflow, while
\\t^\alpha = 1\\ to double precision at every positive double \\t\\; the
density is then \\\alpha e^{-\beta x^2 + \gamma x} / (x (1 + w))\\, with
\\w = \alpha K_0(z)\\ and \\K_0(z) = \int_0^1 (e^{-t^2 + z t} - 1) / t
\\ dt + \int_1^\infty e^{-t^2 + z t} / t \\ dt\\. \\1 / (1 + w)\\ is the
probability that \\X\\ lies below the smallest positive double. Below a
tilt of about 52 (53 at the smallest shapes) it is 1 in double
precision: [`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
is 1 at every positive \\x\\,
[`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md) is 0 for
every lower-tail probability below 1 given on the probability scale (the
upper tail, or `log.p = TRUE`, still gives the quantiles of the far
upper tail), and
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md) returns 0.
For a large tilt \\w\\ is the mass at the peak of the kernel against the
mass by the origin; it is below \\10^{-200}\\ in size for \\z \le 30\\,
and negative below \\z\\ of about 0.3.
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md) and the
moment functions take the same form. It serves the whole law up to a
tilt between about 55 and 56.3, as \\\alpha\\ falls. Above it the mass
by the origin is at most \\2^{-64}\\, about \\5.4 \times 10^{-20}\\, of
the whole, and the density, the upper tail, the moments and the law of
the draws no longer depend on \\\alpha\\; the lower tail far below the
peak, of which that mass is more than \\2^{-64}\\, is still taken from
the same form up to where it reaches one half.

Computation is performed in log-space to avoid numerical
underflow/overflow. For a shape \\\alpha\\ or a scale-free tilt \\z =
\gamma/\sqrt{\beta}\\ above \\10^6\\, and wherever the terms of the
ordinary form grow too large (see **Range**), the log density is
measured instead from the peak \\v^\*\\ of the density of \\\log X\\,
the positive root of \\2\beta x^2 - \gamma x - \alpha = 0\\, as
\\(\alpha - 1)(\log(1+\rho) - \rho) - \rho - \beta d^2\\ less the
logarithm of \\v^\*\\ times the peak-normalised integral of the law,
with \\d = x - v^\*\\ and \\\rho = d/v^\*\\. The ordinary form adds
terms of size \\\alpha \log x\\ and \\(\alpha/2) \log\beta\\, and for
\\z \> 10^3\\ measures `x` from the centre \\\gamma/(2\beta)\\ rounded
to a double; their rounding passes \\10^{-7}\\ at \\\alpha = 10^8\\, and
\\10^{-6}\\ at \\z = 10^{10}\\ unless \\\beta\\ is a power of 2. Here
neither is formed, and \\d\\ is exact: \\v^\*\\ is held as a double next
to it plus its offset from that double.

**Range.** The ordinary form is used only while the larger \\m\\ of
\\\alpha\\ and \\z\\ is at most \\10^6\\ and \\m (1 + \|\log v^\*\|)\\,
the size of its largest terms, is at most \\10^7\\: at \\\beta = 1\\ up
to a shape of \\10^6\\ or a tilt of about \\7.3 \times 10^5\\, and at
\\\beta = 10^{\pm 300}\\ up to about \\3 \times 10^4\\ for either. That
keeps its rounding to about \\2 \times 10^{-9}\\; measured from the peak
the error is about \\10^{-12}\\, relative to the size of the log density
where that is large. So for \\\alpha \le 10^{15}\\, \\\|z\| \le
10^{15}\\ and \\\beta\\ in \\\[10^{-300}, 10^{300}\]\\ the log density
is accurate to the working tolerance `sqrt(.Machine$double.eps)`.
Outside that range the law can be narrower than the spacing of the
doubles near its peak (for \\z\\ above about \\10^{16}\\, or \\\alpha\\
above about \\10^{31}\\), and the density is still that of the exact law
at the double `x` given, almost always 0 there. For a positive tilt
whose mode \\\gamma/(2\beta)\\ lies beyond the largest double, which
needs a subnormal \\\beta\\, `log = TRUE` returns `-Inf` where the log
density is finite. Two cases give `NA` with a warning instead: within 45
standard deviations of the peak of a law narrower than the spacing by a
factor of order \\10^{19}\\ (\\z\\ beyond about \\10^{35}\\, \\\alpha\\
beyond about \\10^{70}\\), where the position of the peak cannot be
resolved finely enough; and where the whole law lies beyond the largest
double without a positive tilt to place it (a shape near the largest
double against a subnormal \\\beta\\). Where \\z\\ itself overflows, the
law is its limit: Gamma(\\\alpha\\, rate \\\|\gamma\|\\) under a
negative tilt and normal about its mode with variance \\1/(2\beta)\\
under a positive one. For \\\alpha = 1\\ and a truncation point \\z_0 =
-\gamma/\sqrt{2\beta}\\ beyond \\10^3\\ the truncated-normal density is
written with the normal hazard \\\phi(z_0)/\bar\Phi(z_0)\\, and below
\\-10^3\\ with the offset \\2\beta x - \gamma\\ formed in one rounding
(as \\\beta x - \gamma/2\\, so also where \\2\beta\\ overflows), so that
it holds to the Exponential and normal limits it tends to.

When any of `alpha`, `beta`, `gamma` is a vector, the density is
evaluated element-wise. The Fox-Wright \\\Psi\\ normalizing constant is
computed once for each position of the recycling cycle while that cycle
is at most 256 long, and otherwise once for each run of consecutive
elements that share the \\(\alpha, \beta, \gamma)\\ triple. So recycled
parameters, and parameter vectors in which equal triples are adjacent,
are much faster than triples that alternate along a long vector, and
than calling `dmhn` inside an R loop.

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
