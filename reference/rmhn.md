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
none of 1000 draws of the proxy described under Subnormal shapes in
Details, each of up to 1000 proposals), that element is `NaN`, and a
single warning reports how many draws were affected; the acceptance
rates of the samplers make this vanishingly unlikely. Where the
distribution lies beyond the largest double precision number, and at
some other extreme values, `rmhn` instead stops with an error for the
whole call (see Details).

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

Two closed-form special cases are taken first, whatever `method`: a
negligible tilt, \\\|\gamma\|/\sqrt{\beta}\\ within
`sqrt(.Machine$double.eps)` of 0, is drawn as the square root of a Gamma
variate, and a shape \\\alpha\\ within `sqrt(.Machine$double.eps)` of 1
from a truncated normal. Both are exact samplers of their special-case
law, which next to \\\gamma = 0\\ or \\\alpha = 1\\ is within a total
variation distance of \\8.6 \times 10^{-9}\\ of the law requested. The
one exception is a negligible tilt above \\\alpha = 2^{60}\\, which is
drawn by the Normal proposal of Sun et al. (2023) Algorithm 1 centred at
the mode, or by RTDR under `method = "rtdr"` (Sun et al. state Algorithm
1 for \\\gamma \> 0\\; centred at the mode, its Normal envelope holds
for either sign). Otherwise `method` selects the rejection sampler:

- `"auto"` (default): the faster of the provably correct samplers, as
  measured by timing them against each other. For \\\alpha \< 1\\ and
  \\\gamma \> 0\\ it uses RTDR. Otherwise it chooses between Sun et al.
  Algorithm 1 (\\\gamma \> 0\\) or Algorithm 3 (\\\gamma \< 0\\) and
  RTDR by the shape, the tilt \\\gamma/\sqrt{\beta}\\ and the number of
  variates per setup (see Vector parameters and setups); Section 7 of
  [`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md)
  gives the decision tree. The choice changes the speed and the stream
  of draws (see Reproducibility), not the law drawn.

- `"rtdr"`: the Relaxed Transformed Density Rejection method of Gao &
  Wang (2025). Gao & Wang bound its acceptance probability below by
  \\1/e \approx 0.368\\ for every parameter value, provided the contact
  points of the envelope lie in the band of their Eq. (7); the search
  for them does not always reach that band, and the acceptance measured
  stays above the bound there too. Where RTDR cannot place a law of tiny
  shape under a strong negative tilt, Sun et al. Algorithm 3 is used
  instead (see Limits of double precision). Note: Gao & Wang (2025) use
  the parameterization \\(\lambda, \alpha, \beta)\\ with density
  proportional to \\x^{\lambda - 1} \exp(-\alpha x^2 - \beta x)\\; the
  mapping to the Sun et al. parameterization used here is \\\lambda
  \leftrightarrow \alpha\\, \\\alpha \leftrightarrow \beta\\, \\\beta
  \leftrightarrow -\gamma\\ (sign flip on the linear term).

- `"sun"`: the algorithms of Sun et al. (2023): Algorithm 1 for \\\gamma
  \> 0\\ and \\\alpha \> 1\\, and Algorithm 3 for \\\gamma \< 0\\, up to
  \\\alpha = 2^{59}\\, beyond which the call stops with an error (see
  Limits of double precision). The combination \\\alpha \< 1\\ with
  \\\gamma \> 0\\ is unsupported: unless a special case answers it
  first, the call stops with the error
  `Sun method not available for alpha<1 and gamma>0; use method="rtdr"`.

### Vector parameters and setups

Vector parameters are recycled to length `n` following standard R rules:
only the first `n` elements of each parameter are used, and any further
elements are silently ignored, matching the convention of `rnorm`.

The setup of the chosen sampler is built once for each run of equal
consecutive \\(\alpha, \beta, \gamma)\\ triples. Short parameter vectors
repeat their triples after a period, the least common multiple of the
lengths of `alpha`, `beta` and `gamma`; when that period is at most 256
draws and they are recycled over at least four such periods, a setup is
kept for each position within the period. So passing parameters grouped
by triple, or as short recycled vectors, is faster than calling `rmhn`
inside an R loop, and a long run of one triple is drawn fastest by its
own scalar call.

The number of variates per setup, from which `method = "auto"` chooses
its sampler, and which therefore changes the draws for a fixed seed, is
`n` for scalar parameters, and `n` divided by the period, rounded down,
where a setup is kept for each position within the period. Otherwise it
is 1: parameters given as vectors of length `n` count one per setup, as
in a Gibbs sampler, even when equal triples are grouped into runs, and
such a call draws what the loop of single draws would. Grouping saves
setups but does not change this number, and so not the choice of
sampler.

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
0. A rounding error that would shift every draw alike is kept to a few
times \\10^{-7}\\ of a standard deviation. For this reason the
square-root-Gamma closed form is used only up to \\\alpha = 2^{60}\\,
and Sun et al. Algorithm 3 only up to \\\alpha = 2^{59}\\ (about \\5.8
\times 10^{17}\\). Above \\\alpha = 2^{59}\\ a negative tilt that is not
negligible is drawn by RTDR under the default method, while
`method = "sun"` stops with an error; a positive tilt is drawn by
Algorithm 1 under both.

The RTDR sampler stops with an error for the whole call once the
standard deviation is no larger than the spacing of the double precision
numbers at the mode, so that its draws could take only a few distinct
values: from \\\gamma/\sqrt{\beta}\\ of about \\9 \times 10^{15}\\
(\\2^{53}\\ at \\\beta = 1\\), and, for \\\gamma \le 0\\, from
\\\alpha\\ between \\2^{103}\\ and \\2^{106}\\ (\\1.0 \times 10^{31}\\
to \\8.1 \times 10^{31}\\), depending on the tilt and on where the mode
falls between two powers of two. There the law is normal to far below
double precision, with mean \\m\\ = `mhn_mode(alpha, beta, gamma)` and
variance \\1/((\alpha - 1)/m^2 + 2\beta)\\, and can be drawn from
directly, as the error message suggests. Sun et al. Algorithm 1, which
the default method uses for such laws with \\\gamma \> 0\\ and \\\alpha
\> 1\\, and the truncated normal of \\\alpha \approx 1\\ do not stop
there: their draws are placed at the mode to within a few units in the
last place, so that they take only a few distinct values near it.

RTDR cannot place the law in three cases, all under a negative tilt:
where \\-\gamma/\sqrt{\beta}\\ passes the largest double; for \\\alpha
\le 1\\, where \\\alpha\\ is below about \\2.5 \times 10^{-324}\\ times
\\-\gamma/\sqrt{\beta}\\, which needs \\\alpha\\ below \\4.4 \times
10^{-16}\\, and almost all the mass of \\X\\ lies below the smallest
positive double; and for a shape of a few, where
\\-\gamma/\sqrt{\beta}\\ comes within a small factor of the largest
double. In the second case, and in the first for \\\alpha \< 10^{-15}\\,
`method = "rtdr"` draws from Algorithm 3, which samples the law exactly
there, except for a subnormal shape, which is drawn as described under
Subnormal shapes; otherwise it stops with an error that names the Gamma
distribution with shape \\\alpha\\ and rate \\-\gamma\\, which the law
is to double precision. The default method uses Algorithm 3 in all three
cases, except where \\-\gamma/\sqrt{\beta}\\ passes the largest double
at a shape above \\2^{59}\\; there it stops as RTDR does. Where the
distribution lies beyond the largest double precision number, no draw
can represent it, and every sampler, the closed forms included, stops
with an error.

### Subnormal shapes

A subnormal shape, below the smallest normal double \\2.2 \times
10^{-308}\\, is drawn as the mixture its law is (see the Details of
[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)): 0 with
probability \\1 / (1 + w)\\, the probability of a value below the
smallest positive double, and otherwise from the peak, through a proxy:
the law of the shape \\2^{-60}\\ above \\\sqrt{\beta} x = 1\\, drawn by
the sampler `method` selects for that shape. The choice between 0 and
the peak compares one uniform variate with \\1 / (1 + w)\\, and so has
the resolution of the generator, \\2^{-32}\\ for R's default. Where \\1
/ (1 + w)\\ is 1 in double precision, below a tilt of about 52 (53 at
the smallest shapes), every draw is 0 and uses no random number; a shape
from \\2^{-1022}\\ up is drawn by the samplers above, which use random
numbers there too, so that the stream after the draws moves as the shape
crosses \\2^{-1022}\\ although every draw is 0 on either side. A refusal
names the shape asked for.

The algorithms, the timings behind the choice of `method = "auto"`, and
the handling of extreme parameters are described in
[`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md)
(Section 7 and Appendix A).

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
