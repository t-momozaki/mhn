# Random Generation from the Modified Half-Normal Distribution

Draws random variates from the Modified Half-Normal (MHN) distribution
with parameters `alpha`, `beta`, and `gamma`.

## Usage

``` r
rmhn(n, alpha = 1, beta = 1, gamma = 0, method = c("auto", "rtdr", "sun"))
```

## Arguments

- n:

  Number of variates to draw: a non-negative number no larger than
  `.Machine$integer.max`, truncated to an integer, so that `n = 0`
  returns `numeric(0)`. Following the base R convention of
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) and
  [`rgamma()`](https://rdrr.io/r/stats/GammaDist.html), if `length(n)`
  is not 1 the number required is taken to be `length(n)`, so that a
  zero-length `n` also returns `numeric(0)`; `NULL` is an error.

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

- method:

  Sampling algorithm. One of `"auto"` (default), `"rtdr"`, or `"sun"`.
  See Details.

## Value

A numeric vector of length `n`. As in
[`rnorm`](https://rdrr.io/r/stats/Normal.html), an element whose
parameters (after recycling to length `n`) are `NA` or `NaN`, or lie
outside the parameter space (`alpha` or `beta` not positive and finite,
or `gamma` not finite), is `NaN`, with a single warning "NAs produced",
and draws no random numbers. A zero-length `alpha`, `beta` or `gamma`
makes every element `NA`, with the same warning. The parameters must be
numeric (double, integer or logical): a factor, a character string, a
complex number, a list or `NULL` is an error. If a rejection sampler
accepts none of 1000 proposals for an element (for a subnormal shape,
none of 1000 draws of the proxy below, each of up to 1000 proposals),
that element is `NaN`, and a single warning reports how many draws were
affected; the acceptance rates of the samplers make this vanishingly
unlikely. Where the distribution lies beyond the largest double
precision number, and at some other extreme values, `rmhn` instead stops
with an error for the whole call (see Details).

## Details

The MHN density is \$\$f(x \mid \alpha, \beta, \gamma) = \frac{2
\beta^{\alpha/2} x^{\alpha-1} \exp(-\beta x^2 + \gamma
x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]} \quad (x \> 0)\$\$ where
\\\Psi\[a, z\] = \sum\_{n \ge 0} \Gamma(a + n/2) z^n / n!\\ is the
Fox-Wright function \\{}\_1\Psi_1\[(a, 1/2); (1, 0); z\]\\. `rmhn` does
not evaluate \\\Psi\\; the rejection-sampling kernels cancel it out. The
one exception is the weight of the origin, \\1 / (1 + w)\\, that it
draws for a subnormal shape (see Subnormal shapes below).

The default parameters `alpha = 1, beta = 1, gamma = 0` correspond to
the half-normal distribution \\\mathrm{HN}(1/\sqrt{2})\\.

### Choice of sampler

The `method` argument selects the rejection sampler:

- `"auto"` (default): Uses closed-form special cases when applicable
  (\\\gamma \approx 0\\ -\> sqrt-Gamma, \\\alpha \approx 1\\ -\>
  truncated normal, where \\\approx\\ means within
  `sqrt(.Machine$double.eps)`, for the tilt on
  \\\|\gamma\|/\sqrt{\beta}\\; away from \\\gamma = 0\\ and \\\alpha =
  1\\ themselves the law drawn is then within a total variation distance
  of \\8.6 \times 10^{-9}\\ of the one requested). Otherwise it selects
  the faster of the provably correct samplers, as measured by timing
  them against each other on an Apple M1 Ultra, for the parameter region
  and the number of variates drawn per setup. The choice changes only
  the speed, not the law drawn, and the crossings move from machine to
  machine. For \\\gamma \> 0\\ it uses RTDR (Gao & Wang 2025) when
  \\\alpha \< 1\\. When \\\alpha \> 1\\ it uses RTDR if \\\alpha \le
  2^{59}\\, the tilt \\z = \gamma/\sqrt{\beta}\\ is at most
  \\2\sqrt{\alpha}\\ and enough variates share a setup, and Sun et
  al. (2023) Algorithm 1 otherwise. Enough means 5, 25, 50 or 100
  variates for \\z \le \sqrt{\alpha}\\, and 25, 50, 100 or 1000 above
  it, for \\\alpha \< 10\\, \\10 \le \alpha \< 1000\\, \\1000 \le \alpha
  \< 10^8\\ and \\\alpha \ge 10^8\\ respectively, and 5 or 25 where
  RTDR's setup starts next to the mode (see below). For \\\gamma \< 0\\
  it uses RTDR above \\\alpha = 2^{59}\\ (see below). Up to that shape
  it uses Sun et al. Algorithm 3 when fewer than 25 (\\\alpha \< 0.1\\),
  5 (\\0.1 \le \alpha \< 1\\), 10 (\\1000 \le \alpha \< 10^5\\), 25
  (\\10^5 \le \alpha \< 10^8\\) or 50 (\\\alpha \ge 10^8\\) variates
  share a setup, for \\-\gamma/\sqrt{\beta} \> 10^{150}\\ and where RTDR
  cannot place the distribution (see below), and RTDR otherwise, which
  for \\1 \< \alpha \< 1000\\, and wherever RTDR's setup starts next to
  the mode, means at any number of variates per setup. For \\\alpha \>
  1\\ that setup searches for the two contact points of RTDR's envelope.
  It starts the search next to the mode \\m\\ of \\\sqrt{\beta}X\\ at
  very large shapes, where it builds the envelope relative to the mode
  (from about \\2.4 \times 10^{12}\\ to \\3.1 \times 10^{12}\\ at
  moderate tilts, from about \\2.3 \times 10^{13}\\ near a tilt of
  \\-\alpha/e\\, and from about \\4.5 \times 10^{13}\\ at any tilt), and
  under a very strong negative tilt, \\z \< -10^{10}\sqrt{\alpha}\\; the
  default method decides this with the setup's own test. Elsewhere the
  search walks the two points in from \\m/2\\ and \\m + \max(1, m/2)\\,
  and from a shape of about 1000 on that makes RTDR's setup dearer than
  Algorithm 3's. The number of variates per setup is counted as `n` for
  scalar parameters, and as `n` divided by the least common multiple of
  the lengths of `alpha`, `beta` and `gamma`, the period after which the
  recycled triples repeat, rounded down, where a setup is kept for each
  position within that period (see below). Otherwise it is 1: parameters
  given as vectors of length `n` count one per setup, as in a Gibbs
  sampler, even when equal triples are grouped into runs, and such a
  call draws what the loop of single draws would. A long run of one
  triple is drawn fastest by its own scalar call.

- `"rtdr"`: Use the Relaxed Transformed Density Rejection method of Gao
  & Wang (2025) for the general case. Gao & Wang bound the acceptance
  probability below by \\1/e \approx 0.368\\ uniformly over the
  parameter space where the contact points of the envelope lie in their
  band (their Eq. (7)); the search for them does not always reach that
  band, and the acceptance measured stays above the bound there too. The
  closed-form special cases are still taken first: a negligible tilt
  draws from the square root of a Gamma, and \\\alpha \approx 1\\ from a
  truncated normal. Both are exact samplers of their special-case law,
  which next to \\\gamma = 0\\ or \\\alpha = 1\\ is within a total
  variation distance of \\8.6 \times 10^{-9}\\ of the law requested, and
  faster than any rejection scheme, so no sampler choice overrides them,
  except that a negligible tilt above \\\alpha = 2^{60}\\ is drawn by
  RTDR (see below). Where RTDR cannot place a law of tiny shape under a
  strong negative tilt, Sun et al. Algorithm 3 is used instead (see
  below). Note: Gao & Wang (2025) use the parameterization \\(\lambda,
  \alpha, \beta)\\ with density proportional to \\x^{\lambda - 1}
  \exp(-\alpha x^2 - \beta x)\\; the mapping to the Sun et al.
  parameterization used here is \\\lambda \leftrightarrow \alpha\\,
  \\\alpha \leftrightarrow \beta\\, \\\beta \leftrightarrow -\gamma\\
  (sign flip on the linear term).

- `"sun"`: Use the Sun et al. (2023) algorithms for the general case,
  with the same closed-form special cases taken first. Algorithm 1 is
  used when \\\gamma \> 0\\ and \\\alpha \> 1\\; Algorithm 3 is used
  when \\\gamma \le 0\\, up to \\\alpha = 2^{59}\\, beyond which the
  call stops with an error (see below). The combination \\\alpha \< 1\\
  with \\\gamma \> 0\\ is unsupported and triggers an error, unless a
  special case answers it first.

### Vector parameters and setups

Vector parameters are recycled to length `n` following standard R rules:
only the first `n` elements of each parameter are used, and any further
elements are silently ignored, matching the convention of `rnorm`.

Internally the setup of the chosen sampler is built once for each run of
equal consecutive \\(\alpha, \beta, \gamma)\\ triples and, when short
parameter vectors whose triples repeat after a period of at most 256
draws are recycled over at least four such periods, at most once for
each position within that period, so passing parameters grouped by
triple, or as short recycled vectors, is faster than calling `rmhn`
inside an R loop. Grouping saves setups but does not change the number
of variates per setup counted above, and so not the choice of sampler.

### Reproducibility

Under one seed, `method = "rtdr"` and `method = "sun"` draw one stream
whatever `n`: the first `k` draws of `rmhn(n, ...)` are those of
`rmhn(k, ...)`, and draws split over several calls are those of one
call. The default method does not keep this, as it chooses its sampler
from the number of variates per setup (see above): `rmhn(4, 5, 1, 2)`
and the first four draws of `rmhn(100, 5, 1, 2)` come from different
samplers. Its stream therefore changes with `n`, with how the parameters
are arranged, and with a version of the package that moves a threshold,
although the law drawn does not. Where the draws must not depend on `n`,
as when a pilot run is compared with the main one, set `method = "rtdr"`
or `method = "sun"`.

### Limits of double precision

Double precision limits which laws can be sampled, and how finely. Every
sampler draws from the exact law and rounds the draw to double
precision, and a draw below the smallest positive double is returned as
0. As \\\alpha\\ or \\\gamma/\sqrt{\beta}\\ grows, the distribution
narrows around its mode relative to the mode itself, and a rounding
error that moved every draw alike, such as that of a computed mode,
would become a visible part of its standard deviation. Wherever such a
shift could exceed a few times \\10^{-7}\\ of a standard deviation, the
draws are placed from the mode computed to about twice double precision,
so that each is the exact draw rounded once. R's
[`rgamma`](https://rdrr.io/r/stats/GammaDist.html), which the
square-root-Gamma closed form and Sun et al. Algorithm 3 draw from,
rounds \\\sqrt{a - 1/2}\\ and the scale, and so shifts every variate of
shape \\a\\ alike by up to about \\3 \times 2^{-53}\sqrt{a}\\ standard
deviations, \\2.5 \times 10^{-7}\\ at \\a = 2^{59}\\. Algorithm 3's own
setup, which rounds the shape, power and rate of its Gamma proposal,
adds a shift of the same kind, and its draws were measured to shift by
up to about \\4.5 \times 2^{-53}\sqrt{\alpha}\\ standard deviations in
all (\\3.8 \times 10^{-7}\\ at \\\alpha = 2^{59}\\). Both are used only
up to a Gamma shape of \\2^{59}\\ (about \\5.8 \times 10^{17}\\), at
\\\alpha\\ up to \\2^{60}\\ and \\2^{59}\\ respectively. Above \\\alpha
= 2^{60}\\ a negligible tilt of either sign is drawn by the Normal
proposal of Algorithm 1 centred at the mode (by RTDR under
`method = "rtdr"`): Sun et al. (2023) state Algorithm 1 for \\\gamma \>
0\\ and Theorem 1 for \\\gamma \ge 0\\, but centred at the mode the
Normal envelope holds for either sign. Above \\\alpha = 2^{59}\\ a
negative tilt that is not negligible is drawn by RTDR under the default
method, while `method = "sun"` stops with an error.

The RTDR sampler stops with an error for the whole call once the
standard deviation is no larger than the spacing of the double precision
numbers at the mode: from \\\gamma/\sqrt{\beta}\\ of about \\9 \times
10^{15}\\ (\\2^{53}\\ at \\\beta = 1\\), and, for \\\gamma \le 0\\, from
\\\alpha\\ between \\2^{103}\\ and \\2^{106}\\ (\\1.0 \times 10^{31}\\
to \\8.1 \times 10^{31}\\). The standard deviation of the law over its
mode \\m\\ is about \\1/\sqrt{2\alpha}\\ at \\\gamma = 0\\ and
\\1/\sqrt{\alpha}\\ under a strong negative tilt, where the law is close
to the Gamma distribution with shape \\\alpha\\ and rate \\-\gamma\\,
and the doubles at \\m\\ are between \\2^{-53} m\\ and \\2^{-52} m\\
apart; so the stop comes between \\2^{103}\\ and \\2^{105}\\ (\\4.1
\times 10^{31}\\) at \\\gamma = 0\\ and between \\2^{104}\\ and
\\2^{106}\\ under a strong negative tilt, depending on where \\m\\ falls
between two powers of two. There the law is normal to far below double
precision, with mean \\m\\ = `mhn_mode(alpha, beta, gamma)` and variance
\\1/((\alpha - 1)/m^2 + 2\beta)\\, and can be drawn from directly. Sun
et al. Algorithm 1, which the default method uses for such laws with
\\\gamma \> 0\\ and \\\alpha \> 1\\ (up to \\\alpha = 2^{59}\\ their
tilt exceeds \\2\sqrt{\alpha}\\, and above that shape Algorithm 1 serves
at any tilt), and the truncated normal of \\\alpha \approx 1\\ do not
stop there: their draws are placed at the mode to within a few units in
the last place, so that they take only a few distinct values near it.
For a strongly negative tilt the law approaches the Gamma distribution
with shape \\\alpha\\ and rate \\-\gamma\\, which Algorithm 3 draws at
any tilt up to its shape limit. RTDR works with \\\sqrt{\beta}\\x\\, or
for \\\alpha \le 1\\ with its logarithm \\y = \log(\sqrt{\beta}\\x)\\,
and cannot place the law in three cases. The first is where
\\-\gamma/\sqrt{\beta}\\ passes the largest double. The second is for
\\\alpha \le 1\\ and a finite \\-\gamma/\sqrt{\beta}\\: the mode of
\\\log(\sqrt{\beta}X)\\ is the logarithm of about
\\\alpha\sqrt{\beta}/\|\gamma\|\\, and where \\\alpha\\ is below about
\\2.5 \times 10^{-324}\\ (half the smallest positive double) times
\\-\gamma/\sqrt{\beta}\\, which needs \\\alpha\\ below \\4.4 \times
10^{-16}\\, that number rounds to 0; almost all the mass of \\X\\ then
lies below the smallest positive double. The third is for a shape of a
few, where \\-\gamma/\sqrt{\beta}\\ comes within a small factor of the
largest double. In the second case, and in the first for \\\alpha \<
10^{-15}\\, but for a subnormal shape (see below), `method = "rtdr"`
draws from Algorithm 3, which samples the law exactly there; otherwise
it stops and names that Gamma law. The default method uses Algorithm 3
in all three cases, except where \\-\gamma/\sqrt{\beta}\\ passes the
largest double at a shape above \\2^{59}\\; there it stops as RTDR does.
Where the distribution lies beyond the largest double precision number,
no draw can represent it, and every sampler, the closed forms included,
stops with an error.

### Subnormal shapes

A subnormal shape, below the smallest normal double \\2.2 \times
10^{-308}\\, is drawn as the mixture its law is (see the Details of
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)): 0 with
probability \\1 / (1 + w)\\, the probability of a value below the
smallest positive double, and otherwise from the peak, as the law of the
shape \\2^{-60}\\ above \\\sqrt{\beta} x = 1\\, which has the peak's
density there to within a relative \\4 \times 10^{-18}\\, by the sampler
`method` selects for that shape. The choice between 0 and the peak
compares one uniform variate with \\1 / (1 + w)\\, and so has the
resolution of the generator, \\2^{-32}\\ for R's default. Where \\1 /
(1 + w)\\ is 1 in double precision, below a tilt of about 52 (53 at the
smallest shapes), every draw is 0 and uses no random number; a shape
from \\2^{-1022}\\ up is drawn by the samplers above, which use random
numbers there too, so that the stream after the draws moves as the shape
crosses \\2^{-1022}\\ although every draw is 0 on either side. A refusal
names the shape asked for.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
[doi:10.1080/03610926.2021.1934700](https://doi.org/10.1080/03610926.2021.1934700)

Gao, F., & Wang, H.-B. (2025). Generating modified-half-normal random
variates by a relaxed transformed density rejection method.
*Communications in Statistics - Simulation and Computation*, 1–25.
Advance online publication.
[doi:10.1080/03610918.2025.2524551](https://doi.org/10.1080/03610918.2025.2524551)

Robert, C. P. (1995). Simulation of truncated normal variables.
*Statistics and Computing*, 5(2), 121–125.
[doi:10.1007/BF00143942](https://doi.org/10.1007/BF00143942)

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)

## Examples

``` r
set.seed(1)
rmhn(10, alpha = 2, beta = 1, gamma = 0.5)
#>  [1] 0.4882149 1.7209460 0.2702322 0.5039321 1.0199839 0.8549756 0.5065696
#>  [8] 1.1410072 0.7866116 1.0854880

# Vector parameters are recycled to length n.
set.seed(1)
rmhn(5, alpha = c(1, 2, 3, 4, 5))
#> [1] 0.4429697 0.7989284 0.5821857 1.4320656 2.2998471
```
