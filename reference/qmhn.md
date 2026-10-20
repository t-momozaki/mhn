# Quantile Function of the Modified Half-Normal Distribution

Computes the quantile (inverse cumulative) function of the Modified
Half-Normal (MHN) distribution with parameters `alpha`, `beta`, and
`gamma`.

## Usage

``` r
qmhn(p, alpha = 1, beta = 1, gamma = 0, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  Numeric vector of probabilities.

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

  Logical; if `TRUE` (default), probabilities are \\P(X \le x)\\,
  otherwise \\P(X \> x)\\, \\x\\ being the quantile returned.

- log.p:

  Logical; if `TRUE`, probabilities are provided on the log scale.
  Default: `FALSE`. Both flags are read with
  [`as.logical()`](https://rdrr.io/r/base/logical.html), so that 1 is
  `TRUE`; a value that gives `NA` is an error.

## Value

A numeric vector. The output length equals
`max(length(p), length(alpha), length(beta), length(gamma))`, or 0 if
any of them has length zero; each input is recycled to that length
following standard R recycling rules. As in
[`qgamma`](https://rdrr.io/r/stats/GammaDist.html), an element whose
parameters lie outside the parameter space (`alpha` or `beta` not
positive and finite, or `gamma` not finite) is `NaN`, with a warning,
and an element with an `NA` or `NaN` among its inputs is `NA` or `NaN`
(`NA` if any of them is `NA`), without one; the other elements are
unaffected. The arguments must be numeric (double, integer or logical):
a factor, a character string, a complex number, a list or `NULL` is an
error, as in [`qnorm`](https://rdrr.io/r/stats/Normal.html).
`qmhn(0) = 0` and `qmhn(1) = Inf`. Probabilities outside \\\[0, 1\]\\
(log probabilities above 0 when `log.p = TRUE`) yield `NaN`, with a
warning as in [`qnorm`](https://rdrr.io/r/stats/Normal.html). If the
distribution function cannot be evaluated where the root is sought, a
warning is issued and `NA` is returned for that element. The attributes
of `p` (names, `dim`, `dimnames`) are kept when the result has the
length of `p`; a classed `p`, such as a Date, is taken as its numbers
and gives a plain vector. Unlike `qnorm`, attributes of a longer
parameter vector are not carried over.

## Details

For the general case, \\q = F^{-1}(p)\\ is obtained by a TOMS 748
root-finder applied to the distribution function computed as in
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md), or for
`lower.tail = FALSE` to the logarithm of its directly computed upper
tail. Near 1 that logarithm is about \\-(1 - p)\\, and its error, set by
the size of the terms it is assembled from (see `pmhn`), is a relative
error in \\1 - p\\. Where that could pass \\10^{-9}\\, within about \\3
\times 10^{-6}\\ of 1 at the least and further out for a large shape,
tilt or \\\|\log\beta\|\\, the distribution function is solved for one
minus the probability instead, which is exact. Where `pmhn` takes that
distribution function from the series of Sun et al. (2023), which can
cost far more than the upper tail (\\\alpha \< 1\\ with \\\|z\|\\ up to
about 830, and \\\alpha \le 100\\ with \\\|z\| \le 20\\, \\z =
\gamma/\sqrt{\beta}\\), the root of the upper tail is refined instead by
Newton steps on the logarithm of the distribution function (taken as in
`pmhn` for the log upper tail). A lower-tail `p` there is solved on the
series itself, and for \\\alpha \< 1\\ with \\\|z\|\\ from about 100 one
quantile can take from a tenth of a second to seconds (about 1 s at
\\\alpha = 0.8\\, \\z = 829\\ and \\p = 0.5\\); where \\1 - p\\ is known
exactly, `lower.tail = FALSE` gives the same quantile at once. Within
\\10^{-7}\\ of 1, where the upper tail may not resolve one minus the
probability at all, that root is taken for a probability it does resolve
where needed, and the steps start from it. Likewise a lower-tail `p`
close to 1 (within \\10^{-7}\\ of 1 at least; see **Range**), where the
distribution function is resolved only to the rounding of 1, is solved
on the logarithm of the upper tail for one minus it, which is exact. A
lower-tail `log.p` probability below `log(.Machine$double.xmin)`, about
-708.4, whose exponential is no longer a normal double, is solved on the
logarithm of the distribution function. The initial bracket is
\\\[s\sqrt{\epsilon},\\ s\]\\ with \\s = E(X) + 8
\sqrt{\mathrm{Var}(X)}\\, which scales with the distribution. Its left
end is halved (down to the smallest positive double) and its right end
doubled (up to the largest double) until it brackets the target
probability, each step keeping the bracket a factor of two wide. Where
the root-finder does not reach its tolerance on that bracket, as for a
tiny lower-tail `p` with a large \\\alpha\\, the equation is solved
again for \\\log q\\ on the logarithm of the distribution function, and
where that fails too, or the bracket lies near the largest double, the
bracket is bisected. The root is found to a relative \\2^{-46}\\, and
where `pmhn` measures the law from its peak (see there) also to within
\\2^{-46}\\ of the law's standard deviation, which can be narrower. A
quantile below `.Machine$double.xmin` is resolved to the spacing of the
subnormal doubles there, about \\4.9 \times 10^{-324}\\, and one below
that spacing is returned as 0, as
[`qgamma`](https://rdrr.io/r/stats/GammaDist.html) returns a quantile
that underflows; a quantile beyond the largest double is returned as
`Inf`, as `qgamma` returns one that overflows.

**Range.** The quantile inverts
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md), whose
range of accuracy it inherits: for \\\alpha \le 10^{15}\\, \\\|z\| \le
10^{15}\\ and \\\beta\\ in \\\[10^{-300}, 10^{300}\]\\ it is the root,
to the tolerance above or to the spacing of the doubles near it,
whichever is larger, of an equation that resolves the smaller of the two
tail probabilities to about \\2 \times 10^{-9}\\ of itself or better
(see there), which places it within about that many of the law's
standard deviations. Next to the origin, where \\F(x)\\ is about a
multiple of \\x^\alpha\\, its relative error is that of the smaller tail
probability divided by about \\\alpha\\. A lower-tail `p` close to 1 is
solved on the upper tail, as the root of \\P(X \> x) = 1 - p\\ with
\\1 - p\\ formed exactly, wherever the error of \\F\\ near 1 could pass
\\10^{-9}\\ of \\1 - p\\ and the upper tail resolves it more accurately:
\\F\\ is resolved there only to the rounding of 1, and for a small
shape, where it is close to 1 over most of the range, to the error of
the lower tail (see
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)). At
parameters of order one that is so within about \\10^{-5}\\ of 1, and
farther for a small or large shape, a large tilt or an extreme
\\\beta\\. Where the law is narrower than the spacing of the doubles
(for \\z\\ above about \\10^{16}\\, or \\\alpha\\ above about
\\10^{31}\\) its distribution function is a step between two adjacent
doubles, and the quantile of every `p` strictly inside \\(0, 1)\\ is the
one of the two nearer to the true quantile, placed from the exact offset
of the law's peak. Where `pmhn` is `NA` (see there), so is the quantile,
with a warning. Where the log upper tail of `pmhn` is `-Inf` in place of
a finite value (\\\beta\\ and the tilt both close to the largest double,
see there), the quantile of a `log.p` probability beyond that point is
the point at which it turns `-Inf`, without a warning.

At \\\gamma = 0\\ and at \\\alpha = 1\\ exactly the quantile is computed
from the law it reduces to; a shape or a tilt merely close to those
values is a different law, and is computed as in the general case:

- \\\gamma = 0\\: \\q = \sqrt{y / \beta}\\, where \\P(\alpha/2, y) = p\\
  for the regularized incomplete gamma function \\P\\ (Sun et al. 2023,
  Lemma 6a); \\y\\ is taken from `qgamma`, or from the leading term of
  \\P\\ where it is below the smallest normal double; for \\\alpha \le
  10^6\\. Above it \\\gamma = 0\\ is solved on the distribution function
  measured from the peak, as for any other tilt (see
  [`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)).

- \\\alpha = 1\\: the truncated-normal inverse (Sun et al. 2023, Lemma
  6b), by `qnorm`, refined by one Newton step on the scaled
  complementary error function where the truncation point 0 lies three
  or more standard deviations above the mean \\\gamma/(2\beta)\\ of the
  untruncated normal, \\-\gamma/\sqrt{2\beta} \ge 3\\, and by Newton's
  method on that function alone where it lies more than 30 standard
  deviations above it. Near the truncation point it is solved by
  Halley's method on a short series, so that the quantile of a small
  lower-tail probability keeps its relative accuracy. Where
  \\-\gamma/\sqrt{2\beta}\\ overflows the law is Exponential with rate
  \\\|\gamma\|\\, and its quantile is that one's.

When any of `alpha`, `beta`, `gamma` is a vector, the quantile is
evaluated element-wise. The Fox-Wright \\\Psi\\ normalizing constant and
moments \\E(X)\\, \\\mathrm{Var}(X)\\ (used to size the root-finder
bracket) are computed once for each position of the recycling cycle
while that cycle is at most 256 long, and otherwise once for each run of
consecutive elements that share the \\(\alpha, \beta, \gamma)\\ triple.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)

## Examples

``` r
# Basic evaluation
qmhn(c(0.1, 0.5, 0.9), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.4717434 1.0906276 1.8480472

# Round-trip: F(F^-1(p)) ~ p
p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
all.equal(pmhn(qmhn(p, alpha = 2, beta = 1, gamma = 1),
               alpha = 2, beta = 1, gamma = 1),
          p, tolerance = 1e-6)
#> [1] TRUE

# Tail / log forms
qmhn(0.95, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#> [1] 0.3394014
qmhn(log(0.05), alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#> [1] 0.3394014
```
