# Changelog

## mhn 0.1.2

Mainly a bug-fix release, written for users of 0.1.1; the default method
of [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) also
chooses its sampler by re-derived thresholds (see Seeded draws below and
Performance).

0.1.1 returns wrong values without a warning in several places, some of
them at ordinary parameter values: the skewness and kurtosis at a large
shape (`mhn_kurtosis(1000, 1, 0)` is -0.627 where the true value is
7.5e-7), the variance for a shape close to 1, the distribution function
under a strong negative tilt when the mode or the quantile lies below
about 1.5e-8, the logarithm of an upper tail close to one and the
quantile of such a tail, and the density, distribution function and
quantile for a shape or a tilt within 1.5e-8 of `alpha = 1` or
`gamma = 0` but not equal to it. A few of these came with 0.1.1 itself,
and 0.1.0 is right there. The true values quoted below come from
independent references: arbitrary-precision computations, or the Gamma
and normal laws the distribution approaches in the limit. Each fix
carries a regression test at the values that failed.

### What changes relative to 0.1.1

**Deterministic functions.** On the grid of 330 parameter triples
described under 0.1.1 below, in the 130 cells with
`|gamma| / sqrt(beta)` at most 10,
[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
and [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
agree with 0.1.1 to twelve significant digits in 128 cells and differ by
at most 5.3e-11 relative in the other two;
[`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) at `p` =
0.001, 0.1, 0.5, 0.9 and 0.999 differs in 13 cells, by at most 5.9e-9;
[`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
does not change.
[`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
and
[`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
change by up to 2.2e-5 and 2.8e-4 relative, for `alpha <= 1/2` at
`gamma / sqrt(beta) = 10`, where 0.1.1 is the inaccurate one (see
below). Outside these cells values change by more, and those changes are
the fixes listed below.

**Seeded draws.** For a fixed seed
[`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) returns
the draws of 0.1.1 over most of the parameter space wherever it uses the
sampler 0.1.1 used. The default method now uses RTDR in parts of the
space where 0.1.1 used Sun et al.’s Algorithm 1 or 3 (see Performance),
and there its draws follow a different sequence. On the 330-cell grid
(200 draws per cell under `set.seed(20260903)`, default method), 300
cells give the same draws and 30 a different sequence: 3 (`alpha < 1/2`
at a tilt of 1e7) where RTDR itself changed, and 27 where the default
method changed sampler (12 with `gamma > 0` and `alpha > 1` at a tilt up
to `2 sqrt(alpha)`, 15 with `gamma < 0` and `alpha >= 10`). In more
detail:

- Unchanged: RTDR below a standardised tilt `gamma / sqrt(beta)` of
  about 5e6 for `alpha < 1/2`, 1.6e7 for `1/2 <= alpha < 1` and 7.7e6
  for `alpha > 1` (6.1e6 at a shape of 1e12), and at every negative tilt
  checked (down to -1e8), both on the default path and under
  `method = "rtdr"`, at shapes up to 1e12 (apart from a rare change in
  the trailing digits at an extreme scale); the square-root-Gamma case
  `gamma = 0` for shapes from 0.05 to 1e16; and the truncated normal of
  `alpha = 1` up to a tilt of about 1e8, or at any tilt when `beta = 1`.
- The same sequence, with each draw changed only slightly: Sun et al.
  Algorithm 1 (`gamma > 0`, `alpha > 1`), in the last bit or two, for
  shapes below 1e6 at every `beta` checked (1e-20 to 1e20), as does the
  truncated normal of `alpha = 1` beyond a tilt of 1e8 with `beta != 1`;
  and Algorithm 3, by at most about 1e-12 relative for `beta <= 1`, a
  few 1e-12 for `beta` up to 10, 1e-9 at `beta = 1e6` and 1e-7 at
  `beta = 1e10`, and at `beta = 1e20`, where the search of 0.1.1 for the
  matching point had lost precision, by up to 4e-3 relative (0.005
  standard deviations). Algorithm 1 supplies every draw with `gamma > 0`
  and `alpha > 1` under `method = "sun"`, and those of the default
  method there at few draws per setup or a tilt above `2 sqrt(alpha)`;
  Algorithm 3 supplies every draw with `gamma < 0` under
  `method = "sun"`, and those of the default method for `gamma < 0` with
  `alpha < 1` or `alpha >= 1000` at few draws per setup, past a tilt of
  1e150 and where RTDR cannot place the law
  ([`?rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md) gives
  the rule). An accept-reject decision that fell within such a change
  would shift the sequence from there on; none did in the draws
  compared.
- A different sequence: RTDR, and so the default method under a positive
  tilt, at a shape within a few times 2^-1022 and a tilt
  `gamma / sqrt(beta)` of about 52.5 to 55, from the first draw whose
  log-axis value would pass the largest double and that the corrected
  test then rejects; such a draw has a chance of about
  `exp(-1.8e308 alpha)` per proposal from the left tail (e^-4 at
  2^-1022, 1.6e-8 at 1e-307). Where every draw is 0 the values stay and
  the stream after them can move (see the RTDR entry under the
  random-generation fixes).
- A different sequence: every subnormal shape, below 2.2e-308. A draw is
  0 with probability `1 / (1 + w)`, decided by one uniform variate, and
  otherwise drawn as from the shape 2^-60 (see the subnormal-shape entry
  below); below a tilt of about 52 (53 at the smallest shapes) every
  draw is 0 and uses no random number, so that the draws after it in the
  stream change too.
- A different sequence: RTDR wherever it now works in the offset from
  the mode, above the tilts named above and, for `alpha > 1`, at the
  large shapes where it builds its envelope relative to the mode (from
  about 2.4e12 to 3.1e12 at moderate tilts, 2.3e13 near a tilt of
  `-alpha / e` and 4.5e13 at any tilt); `gamma = 0`, and Algorithm 3, at
  a shape below about 0.05, from the first Gamma variate that underflows
  to zero onwards (at `gamma = 0` this happens with probability about
  1.5e-5 per draw at `alpha = 0.03` and 1e-8 at `alpha = 0.05`);
  Algorithm 1 at a shape of about 1e6 or more with `beta != 1`; shapes
  above about 1e14 under a negative tilt, above about 1e16 under a
  positive one and above `2^60` (about 1.2e18) at `gamma = 0`; and every
  parameter set at which 0.1.1 returned `NaN`, stopped with an error, or
  sampled the wrong law (see below).
- A different sequence also where the default method now takes RTDR and
  0.1.1 took Algorithm 1 or 3 (see Performance): `gamma < 0` with
  `1 < alpha < 1000` at fewer than 25 draws per setup or with
  `alpha >= 10`, and with `alpha >= 1000` from 10, 25 or 50 draws per
  setup, and at any number where RTDR’s setup starts next to the mode;
  `gamma < 0` with `0.1 <= alpha < 1` at 5 to 24 draws per setup and
  with `alpha < 0.1` at 25 to 99; and `gamma > 0` with `alpha > 1` at a
  tilt up to `2 sqrt(alpha)` from 5 to 1000 draws per setup, by shape
  and tilt. This includes the single draws of a Gibbs sampler with
  `gamma < 0` and `1 < alpha < 1000`. With recycled parameter vectors
  the draws per setup are also counted differently (see Performance),
  which can change the sampler, and with it the sequence.

**Invalid arguments.**
[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) and
[`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) now treat
a parameter outside the parameter space, `NA` and a zero-length argument
as the base R distribution functions do, element by element, where 0.1.1
stopped the whole call; and a factor or a string is an error, where
0.1.1 used it as a number (see Interface). Valid arguments give the
values and the draws they gave before.

### Bug fixes — moments and mode

- [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
  and
  [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
  were wrong at a large shape under a non-positive tilt, because their
  quadrature range ended where the kernel had fallen by only a few units
  of log density: `mhn_kurtosis(1000, 1, 0)` returned -0.627 and
  `mhn_skewness(1000, 1, 0)` 0.00106, where the true values are 7.5e-7
  and 0.0224; at `alpha = 1e4` the kurtosis was -1.14, and
  `mhn_kurtosis(200, 1, -2)` had the wrong sign (-0.0040 against
  0.00044).

- For `alpha < 1` with a positive tilt strong enough to give the density
  an interior peak, the quadrature panel next to the origin contained
  the singularity of `x^(alpha - 1)`: `mhn_skewness(0.01, 1, 3)` and
  `mhn_kurtosis(0.01, 1, 3)` returned 1.20 and 0.46 against 3.47 and
  12.15, and `mhn_kurtosis(0.3, 1, sqrt(10))` was -0.55609 against
  -0.55730. 0.1.0 was right at these values.

- The normalising constant was taken from its `alpha = 1` or `alpha = 2`
  closed form for any shape within `sqrt(.Machine$double.eps)` of those
  values, which put it out of step with the rest of the computation:
  `mhn_var(1 + 1e-9, 1, -20)` was 7.4e-7 relative off,
  `mhn_var(1 + 1.4e-8, 1, -1e6)` returned 1.0e-7 where the variance is
  1.0e-12, and the density at `alpha = 2 + 1e-8`, `gamma = 5` integrated
  to 1 + 9.6e-9. The closed forms are now used only at `alpha = 1` and
  `2` exactly.
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  likewise treated a shape near 1 as 1: `mhn_mode(1 + 1e-9, 1, -1)` was
  0, where the mode is 1.0e-9.

- [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
  used the closed form of Sun et al. (2023, Lemma 2c) where its
  cancellation was larger than its guard allowed for, and returned
  variances above the bound `1/(2 beta)` that Sun et al. (2023, Lemma
  4c) prove for `alpha >= 1`: 0.50127 at `(2, 1, 900)` against a true
  0.4999988, 6.17 at `(1e8, 1, 1)` and 1074 at `(1e9, 1, -1)`, where the
  variance is about 0.25. Each moment helper now bounds the error of its
  closed form and integrates directly where that bound exceeds the
  working tolerance.

- `mhn_mean(alpha, beta, 0)` lost accuracy at a large shape: at
  `alpha = 1e12` it was 1.9e-4 low, and at `alpha = 1e15` 21.7 times too
  large. It is now accurate to the rounding of its logarithm, about
  `log(alpha) / 2`: within 4e-15 relative up to `alpha = 1e15` and 4e-14
  at any shape up to the largest double.

- [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  returned the inflection point `gamma / (4 beta)` at the boundary
  `alpha = 1 - gamma^2 / (8 beta)` of Sun et al. (2023, Lemma 3c), where
  the density has no interior maximum: `mhn_mode(0.5, 1, 2)` was 0.5. It
  now returns `NA` there, as for any `alpha < 1` without an interior
  mode. (Lemma 3c states its condition with `<=` and places a local
  maximum there; at equality the stationary points merge into an
  inflection point, so the package departs from the lemma at this
  boundary.)

- `mhn_kurtosis(2, 1e300, 1)` was `NaN`; it is now 0.2451.

- For `alpha` just below 1 and a small positive tilt, just past the tilt
  `gamma / sqrt(beta) = sqrt(8 (1 - alpha))` at which the density
  acquires its interior peak,
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
  and
  [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
  were up to 1e-7 off: `mhn_kurtosis(0.99998, 1, 0.013)` returned
  0.85932105 where the true value is 0.85932098. The quadrature now
  takes the bend of `x^(alpha - 1)` right of the peak on its own, and on
  a grid of shapes from 1e-6 to 1 - 1e-8 with an interior peak they are
  within 1.2e-11 of arbitrary-precision references.

- At `alpha = 1` under a positive tilt `gamma / sqrt(beta)` from about
  11.1 on,
  [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  returned values just below the mode `gamma / (2 beta)` and
  [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
  values just above `1/(2 beta)`, against `E(X) >= mode`, which the
  truncated normal satisfies and Sun et al. (2023, Lemma 4a) prove for
  `alpha > 1`, and the bound `Var(X) <= 1/(2 beta)` of their Lemma 4c:
  up to a tilt of 12 by up to 6e-15 of the mean and 2.2e-13, and beyond
  by more (8e-12 of the mean at tilts up to 1000; for the variance see
  the item on the Lemma 2c form above). 0.1.1 did the same for a shape
  within 1.5e-8 of 1, which it treated as 1. Where the closed forms of
  Lemmas 2a and 2c would still cross these bounds next to `alpha = 1`,
  the quadrature is now taken instead, and for `alpha >= 1` a variance
  from the quadrature is capped at `1/(2 beta)`, which is then the
  nearer value.

- Under a strong negative tilt
  [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  could exceed the upper bound of Lemma 4a. At a shape of 6 or more it
  was off by up to 2.8e-9 (`mhn_mean(3e6, 1, -35)` returned
  1216.02602953 where the mean is 1216.02602617), and
  `mhn_mean(1e5, 1, -1e9)` returned 1.0000000000387e-4, where the mean
  and the bound both round to 0.9999999999998e-4. There the mean now
  comes from the quadrature, which keeps to the bound to within an ulp
  and is within 3e-15 of arbitrary-precision references for shapes from
  6 to 3e6. Below a shape of 6 it was off by up to 3.6e-12 for tilts
  down to -1e14, and above the bound at about a third of the points
  checked. It is still the Lemma 2a ratio there, now within 5e-14 of the
  references down to that tilt, and where that error would carry it
  above the bound, from a tilt of about -1e7 on, the bound is returned
  instead, which is closer to the mean.

### Bug fixes — density and normalising constant

- For a small shape under a positive tilt the normalising constant came
  from Sun et al.’s Supplementary Lemma 10 series, whose second term was
  formed from `(alpha / 2 + 1) - 1`: that keeps only the leading bits of
  `alpha / 2`, and none below 1.1e-16, and every later term of the even
  half of the series inherited the error, about `2.2e-16 / alpha`
  relative. Wherever the peak carries the mass, from a tilt of about 12
  at such a shape, the density and the moments were off with it: at
  `alpha = 1e-16` and a tilt from 15 up the density was twice its value
  and `mhn_mean(1e-16, 1, 20)` 19.90 against 9.949; at `1e-15` and
  `1e-14` the density was 6 and 1 percent off at a tilt of 20, and at
  `1e-12` the mean 4.4e-5 off. The term is now formed as
  `alpha / 2 + (k - 1)`, and past a tilt of 20 the constant comes from
  quadrature (see below).

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) treated
  a shape within `sqrt(.Machine$double.eps)` (1.5e-8) of 1 as 1 and a
  tilt `|gamma| / sqrt(beta)` below it as 0, and answered for the
  truncated normal or the square root of a Gamma law; the normalising
  constant, and with it the moments, did the same for the tilt. The law
  requested differs from those by about
  `(alpha - 1) (log(sqrt(beta) x) - E log U)` or
  `gamma / sqrt(beta) (sqrt(beta) x - E U)` in the log density,
  `U = sqrt(beta) X`: `dmhn(1e-300, 1 + 1.4e-8, 1, -1)` was 9.7e-6
  relative off, `pmhn(6, 2.5, 1, 1.4e-8, lower.tail = FALSE)` 7.1e-8 and
  `mhn_mean(2.5, 1, 1.4e-8)` 3.1e-9. The closed forms are now used at
  `alpha = 1` and `gamma = 0` exactly, and the general forms, which keep
  their accuracy up to those values, everywhere else. On 126 triples
  next to them the density, both tails and the quantile agree with an
  arbitrary-precision reference to within 3.3e-11 up to a shape of 50
  (0.1.1: up to 9.7e-6 in the density), and above it to within the
  ordinary accuracy of the general forms, 9e-11 at `alpha = 1e5`, which
  there exceeds the error of taking a tilt below about 1e-12 as 0.
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) still
  draws from those laws within the same window, which its help page
  describes as `alpha ~ 1` and `gamma ~ 0`: the law drawn is within a
  total variation distance of 8.6e-9 of the one requested, and the draws
  for a fixed seed do not change. (The general samplers draw the exact
  law there, but of Sun et al.’s samplers only Algorithm 2, which the
  package does not implement, covers `alpha < 1` with `gamma > 0`, so
  `method = "sun"` would stop next to `alpha = 1` from below under a
  positive tilt, and next to `gamma = 0` from above for any
  `alpha < 1`.)

- For `alpha` above about 1.0e7 with a nonzero tilt, the normalising
  constant lost a factor of 2 and
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md)
  returned twice the density: `dmhn(mhn_mode(1.5e7, 1, 1), 1.5e7, 1, 1)`
  was 1.5958, where the density is 0.79785.
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) were
  wrong in the same range: `pmhn(mhn_mode(1.5e7, 1, 1), 1.5e7, 1, 1)`
  was 1 and `pmhn(mhn_mode(1e8, 1, 1), 1e8, 1, 1)` 0, where both are
  0.5, and `qmhn(0.5, 1e8, 1, 1)` was 7071.66, where the median is
  7071.32.

- Under a negative tilt at a large shape the quadrature for the
  normalising constant could miss part of the mass: at
  `(2795.013, 1, -9.112)` the density was 9.6e-5 relative too low. For a
  shape below 1 its logarithm was off by up to 4.4e-10, at
  `(0.7148, 1, -0.6437)`, and by more than 1e-12 at one in a hundred of
  9.4e4 such points checked; the density and the distribution function
  carry that error. Under a negative tilt it is now within 4.4
  `eps (1 + |log Psi|)` of an arbitrary-precision reference at 1.27e5
  points with shapes from 1e-4 to 2000 and `gamma / sqrt(beta)` from
  -1e5 to -1e-4.

- For `gamma / sqrt(beta)` above 20, and for `alpha` above 100, the
  normalising constant now comes from quadrature, which there is both
  faster and more accurate than the Lemma 10 series:
  `dmhn(500, 2.5, 1, 1000)` was 4.8e-9 relative off and took more than
  10 ms; it now takes a few microseconds.

- The quadratures held the error of each panel to a tolerance that grew
  with the panel’s width, so their accuracy and cost depended on `beta`,
  on which the shape of the law does not depend. At a large `beta` the
  tolerance could not be met: `mhn_skewness(1.5, 1e12, -3e6)` took 38
  ms, against about 45 microseconds for the same law at `(1.5, 1, -3)`.
  Both now take about 9 microseconds.

- A subnormal shape, below the smallest normal double 2.2e-308, was not
  served correctly. There half of `alpha` is not always a double (half
  of `3 * 2^-1074` rounds to `2 * 2^-1074`, and half of `2^-1074` to 0),
  and `1 / alpha` passes the largest double below 5.6e-309. The series
  and the square-root-Gamma case that start from `alpha / 2`, and the
  quadratures that hold the mass by the origin, about `1 / alpha`, went
  wrong with them: `dmhn(1, 2^-1074, 1, 1, log = TRUE)` returned `Inf`
  and `pmhn(0.5, 2^-1074, 1, 1)` `NA`; at `3 * 2^-1074` the log density
  was -743.054 against -743.341 (a third too large as a density),
  `pmhn(1e-300, 3 * 2^-1074, 1, 1)` was 0 against 1,
  `mhn_mean(3 * 2^-1074, 1, 12)` 1.510e-307 against 1.133e-307 and
  `qmhn(1e-320, 3 * 2^-1074, 1, 1, lower.tail = FALSE)` 4.9e-324 against
  2.56e-293; and `mhn_skewness(2 * 2^-1074, 1, -30)` was `-Inf` against
  6.34e161. At such a shape `t^alpha` is 1 to double precision at every
  positive double `t`, and the law takes a form in which `alpha` only
  multiplies: the density is
  `alpha exp(-beta x^2 + gamma x) / (x (1 + w))`, the upper tail
  `alpha J(sqrt(beta) x) / (1 + w)` with
  `J(u) = int_u^Inf exp(-t^2 + z t) / t dt`, the lower tail
  `(1 + alpha (K0(z) - J(sqrt(beta) x))) / (1 + w)`, and the moments
  `E X^k = alpha A_k / ((1 + w) beta^(k/2))` with
  `A_k = int_0^Inf t^(k-1) exp(-t^2 + z t) dt`,
  `z = gamma / sqrt(beta)`, `w = alpha K0(z)` and
  `K0(z) = int_0^1 (exp(-t^2 + z t) - 1) / t dt + int_1^Inf exp(-t^2 + z t) / t dt`.
  `1 / (1 + w)` is the probability that `X` lies below the smallest
  positive double; it is 1 in double precision below a tilt of about 52
  (53 at the smallest shapes), and for a large tilt `w` is the mass at
  the peak of the kernel against the mass by the origin.
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) and the
  four moment functions now use these forms for a subnormal shape up to
  a tilt between about 55 and 56.3, as `alpha` falls. Above it the mass
  by the origin is at most 2^-64 (5.4e-20) of the whole, the density,
  the upper tail and the moments are those of the ordinary paths, and
  the draws follow the same law; the lower tail far below the peak, of
  which that mass is more than 2^-64, and its quantiles still come from
  these forms up to one half (`pmhn(1, 2^-1074, 1, 57, log.p = TRUE)`
  was 0 against -65.03).
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) draws 0
  with probability `1 / (1 + w)`: at `2^-1074` and a tilt of 54.669,
  where that is 1/2, it drew 0 for 73 percent of the draws. While the
  mass at the peak is small (`w` well below 1, up to a tilt of about 53
  to 54.5 as `alpha` falls), the skewness and the excess kurtosis grow
  like `1 / sqrt(alpha)` and `1 / alpha`, by factors that depend on the
  tilt, and the kurtosis is `Inf` where it passes the largest double:
  below a shape of about 1.1e-308 at `z = 0`, at every subnormal shape
  for `z = -30`, and at none for `z = 12`.
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  was already right there.

### Bug fixes — distribution and quantile

- At `gamma = 0`,
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) with
  `log.p = TRUE` lost the far tail opposite to the one given. That case
  hands the probability to R’s
  [`qgamma()`](https://rdrr.io/r/stats/GammaDist.html) (`X^2` follows a
  Gamma law), which returns a wrong quantile for a log probability close
  to 0 at some shapes, of either tail:
  `qmhn(log1p(-1e-300), 1e-10, 1, 0, log.p = TRUE)` returned 10.92,
  where the quantile is 25.70, and at `alpha = 10` the quantile was 0.18
  to 0.33 percent off;
  `qmhn(-1e-300, 500, 1, 0, lower.tail = FALSE, log.p = TRUE)` returned
  2.4267, where the quantile is 2.4561, and at `alpha` from about 200 to
  3e4 it was up to 8.4 percent off. A log probability above `log(1/2)`
  is now passed to [`qgamma()`](https://rdrr.io/r/stats/GammaDist.html)
  as the logarithm of the other tail.

- The lower-tail quadrature of
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md), which
  serves a strong negative tilt, normalised its kernel at a point
  floored at `sqrt(.Machine$double.eps)`, so it was wrong whenever the
  mode or `q` lay below about 1.5e-8: `pmhn(2e-9, 3, 1, -1e8)` returned
  0.0175 and `pmhn(1e-9, 3, 1, -5e6)` 4.3e-6, where the law is Gamma(3,
  rate `-gamma`) to within 1e-11 and the true values are 0.00115 and
  2.1e-8. [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  inherited the error: `qmhn(0.2, 3, 1, -3e8)` returned 8.56e-9 against
  5.12e-9, and `qmhn(0.5, 3, 1, -3e8)` failed with a warning.

- Far from the peak the lower-tail quadrature could fall back towards
  `P(X <= mode)` and stop being monotone:
  `pmhn(c(0.01, 1, 6), 2.5, 1, -3e4)` was `1 1 0.300`, and
  `pmhn(2.47e8, 0.7, 1, -20)` was 0.9595 where it is 1.

- `pmhn(q, lower.tail = FALSE, log.p = TRUE)` took the logarithm of the
  upper-tail quadrature, so where `P(X > q)` is close to one its
  logarithm, about `-F(q)`, kept only the absolute accuracy of the
  quadrature: at the 1e-10 quantile of `(2, 1, 1)` it returned
  -9.9999953e-11 for -1e-10, and at that of `(3, 1, 10)` -1.000195e-10.
  That absolute error is about `eps` times the size of the terms the
  logarithm is assembled from,
  `|log Psi| + (alpha/2) |log beta| + m (1 + |log v*|)`, with
  `gamma^2 / (4 beta)` removed from `log Psi` for `gamma / sqrt(beta)`
  above 1000 (see
  [`?pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)).
  Wherever it could pass 1e-9 of `F`, below the median, the logarithm is
  now formed as `log1p(-F)` from the lower tail: at parameters of order
  one below `F` of about 1e-5, and at every such `q` where the terms are
  largest (at `(2e4, 1e300, 5e151)` the logarithm of the quadrature is
  2.6e-8 off at `F = 0.05`). Elsewhere it stays the logarithm of the
  quadrature, whose relative error there is held below 1e-9 by the error
  bound it is gated on: four times the size-of-terms estimate (about
  twice the largest error measured), plus `8 eps (1 + |log Psi|)` for
  the error of the normalising constant itself (at most 4.4
  `eps (1 + |log Psi|)` on scans of 1.55e5 points). Against an
  arbitrary-precision reference it is at most 1.4e-10 at 1005 such
  points on 390 triples, and 1.1e-10 at 945 points placed from just
  below the gate to ten times it on 244 random triples with `alpha < 1`
  and `gamma / sqrt(beta)` from -3e4 to -20, where 0.1.1 is up to 4.0e-7
  off, and 6.1e-11 at 600 points from just above the gate to ten times
  it on 150 random triples with `alpha` from 5 to 200 and
  `gamma / sqrt(beta)` from -40 to -1. 0.1.0 was right here.

- In the same way `pmhn(q, log.p = TRUE)` is now formed in logarithms
  throughout. Outside the special cases 0.1.1 took the logarithm of `F`:
  near `F = 1`, at the 1e-10 upper quantile of `(2, 1, 1)`, it returned
  -1.00000008e-10 for -1e-10, and wherever `F` fell below the smallest
  positive double, about 4.9e-324, it returned `-Inf`:
  `pmhn(75, 2, 1, 300, log.p = TRUE)` was `-Inf` and is -5631.28. Where
  `F` is a subnormal double that logarithm kept only part of its digits.

- For a shape near 1e-300 at a small `beta`, that lower tail needed one
  more change. Next to the origin the integral of the density is about
  `1 / alpha`, and the lower-tail quadrature measured it against the
  kernel at `q`, which carries the factor `1 / q`: at
  `pmhn(1e145, 1e-300, 1e-300, 1001e-150, log.p = TRUE)` their ratio is
  `exp(1024)`, past the largest double, so that the integral would
  overflow. Where that part exceeds the kernel by more than `exp(600)`
  it is now the measure instead, which leaves the integral as it is; the
  call returns -249803.83 (0.1.1: `-Inf`).

- For a small shape `pmhn(q, log.p = TRUE)` lost the relative accuracy
  of the lower tail below the peak, where `F` is close to 1: its
  logarithm, about `-P(X > q)`, came from terms of the size of
  `log Psi`, `(alpha / 2) log(beta)` and `(alpha - 1) log(q)` that
  cancel. Under a small tilt the error passed 1e-9 below a shape of
  about 1e-6 (7.5e-11 at 1e-5, 1.1e-9 at 1e-6), and under a large tilt
  below a smaller shape. `pmhn(2.5, 1e-100, 1, 5, log.p = TRUE)` was 0
  where it is -1.52e-98, `pmhn(2.5, 1e-12, 1, 5, log.p = TRUE)` 1.1e-4
  off, and `pmhn(2.67e-308, 2^-1022, 1, -0.5, log.p = TRUE)` -1.1e-13
  where it is -1.57e-305. Wherever that error could pass 1e-9 of the
  logarithm and the error of the upper tail, carried into it, is
  estimated at less than half as much, it is now `log(1 - P(X > q))`
  from the upper tail, as the logarithm of the upper tail is formed from
  the lower one on the other side (see above). On scans with shapes from
  2.2e-308 to 1e-4, tilts from -1e5 to 45 and `beta` from 1e-300 to
  1e300 it is within 1e-11 of independent quadratures where it is taken
  (see the next entry for a peak below the smallest normal double).
  Where it is taken, the upper tail is evaluated as well as the lower
  one: 5.0 microseconds a point against 3.0 at `(1e-10, 1, 0.5)`, and
  4.0 against 2.6 at `(1e-100, 1, 5)`.

- Where the peak of the density of `log X` lies below the smallest
  normal double, 2.2e-308, as it can for a tiny shape under a strong
  negative tilt at a large `beta`, the upper-tail quadrature took the
  slope of its integrand at the peak from the peak’s stationarity, which
  holds only at the exact peak, while a subnormal peak is resolved only
  to the spacing of the doubles there:
  `pmhn(1.482e-323, 1.41e-217, 7.65e208, -37.47 * sqrt(7.65e208), lower.tail = FALSE, log.p = TRUE)`
  was 3.5e-7 off. The slope is now formed directly there. On 200 such
  points the logarithm of the upper tail is within 2.2e-16 of an
  independent quadrature, and that of the lower tail within 1e-13.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) solved
  a lower-tail `p` on `F(x) = p` unless `1 - p` was below 1e-7, while
  `F` near 1 is resolved only to the rounding of 1 and, for a small
  shape, to the error of the lower tail above: a relative error of that
  divided by `1 - p`. `qmhn(log1p(-1e-6), 1e-8, 1, 0.5, log.p = TRUE)`
  was 1.3e-7 of itself off, and at `(2, 1, 0.5)` the quantile of
  `p = 1 - 1e-7` left `P(X > q)` 2.1e-9 of `1 - p` away from it, `1 - p`
  taken as the double it is. Wherever that error could pass 1e-9 of
  `1 - p` and the upper tail resolves it more accurately, the quantile
  is now the root of `P(X > x) = 1 - p`, with `1 - p` formed exactly:
  within about 1e-5 of 1 at parameters of order one, and farther for a
  small or large shape, a large tilt or an extreme `beta`. At
  `(2, 1, 0.5)` that difference is now 4e-15.

- `qmhn(1e-100, 2, 1, 3)` returned 1.62e-40, whose probability is
  5.2e-82; it now returns 7.11e-50.

- `qmhn(p, lower.tail = FALSE)` for `p` close to 1 solved on the
  logarithm of the upper tail, which is about `-(1 - p)` there and
  carries the absolute error described above. That is a relative error
  in `1 - p`, and at `1 - p = 1e-6` it moved the quantile by 8.2e-5 of a
  standard deviation at `(9e5, 1, -500)`, 7.3e-5 at
  `(2e4, 1e-300, 5e-151)`, 4.6e-5 at `(0.3, 1, 300)` and 7.4e-4 at
  `(0.5, 1, 900)`. At the last two most of that error is the error of
  the normalising constant under a strong positive tilt (2.3e-10 and
  3.6e-9 in its logarithm, fixed above), which this equation turns into
  a relative error in `1 - p`; with the constant fixed, the equation
  alone leaves 8.9e-7 and 3.4e-6. Where that relative error could pass
  1e-9, the quantile is now the root of `F(x) = 1 - p`, found directly
  where `F` comes from a quadrature, and by Newton steps from the
  upper-tail root where it comes from the series of Sun et al. (within
  1e-7 of 1, where needed, from the root for a `1 - p` the upper tail
  still resolves). On 141 triples at upper-tail `p` from 0.7 to
  `1 - 1e-12` it is within 1.3e-9 of a standard deviation of an
  arbitrary-precision reference, where 0.1.1 is up to 8.8e-3 off down to
  `1 - 1e-7` (at `(0.5, 1e-300, 1e-147)`, again mostly the error of the
  normalising constant; the equation alone left up to 2.6e-3, at
  `(2e4, 1e300, 5e151)`). Its relative error is that of `1 - p` divided
  by `x f(x) / F(x)`, about `alpha` next to the origin: on 178 other
  triples it is up to 6.8e-10, at `alpha` of about 0.015, against up to
  2.9e-5 in 0.1.1 down to `1 - 1e-7`; closer to 1, 0.1.1 missed by
  2.7e-3 or more at a tenth of the points. Where the quantile lies below
  the smallest positive double, as it can for `alpha < 1` under a
  negative tilt, 0.1.1 returned a positive value
  (`qmhn(log1p(-1e-35), 0.0677, 1, -267.08, lower.tail = FALSE, log.p = TRUE)`
  was 3.9e-229); it now returns 0, as the lower-tail quantile of `1 - p`
  does.

- At `gamma = 0` the square-root-Gamma forms squared `x`, which
  underflows below 1.5e-154: `dmhn(1e-170, 1.5, 1, 0)` was `Inf`,
  `pmhn(1e-200, 0.01, 1, 0)` was 0 where it is 0.0100, and
  `qmhn(0.01, 0.01, 1, 0)` was 0 where it is 7.5e-201.

- At `alpha = 1` the lower tail and the quantile cancelled:
  `pmhn(1e-16, 1, 1, -3)` was 0 where it is 3.5e-16, and
  `qmhn(1e-17, 1, 1, -3)` was 0 where it is 2.85e-18.

- When [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  cannot bracket a root it now returns `NA`, as its warning says; 0.1.1
  returned `NaN`. A
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) call
  with scalar parameters can now be interrupted.

- For a shape below about 1e-11 under a positive tilt, the upper tail at
  a very small `q` left out the stretch next to the origin, which holds
  most of the tail there:
  `pmhn(1e-200, 1e-20, 1, 1, lower.tail = FALSE)` returned 3.29e-19
  where the value is 4.61e-18, and
  `pmhn(1e-30, 1e-16, 1, 1, lower.tail = FALSE)` 6.954e-15 where it is
  7.003e-15.
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) with
  `lower.tail = FALSE` inherited the error:
  `qmhn(1e-12, 1e-14, 1, 1, lower.tail = FALSE)` returned 9.660e-44 and
  now returns 9.616e-44.

### Bug fixes — random generation

- RTDR drew 0 too often at a shape just above the smallest normal double
  where the law puts a fair share of its mass below the smallest
  positive double: at `alpha = 2^-1022` and a tilt of 53.333, where that
  share is 1/2,
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) drew 0
  for 50.8 percent of the draws (15 standard errors over 1e6 draws), and
  at `1.25 * 2^-1022` and 53.325, where it is 1/2 too, for 50.3 percent.
  On the log axis the leftmost piece of the envelope falls like
  `exp(alpha y)`, and a draw from it, `a + log(u) / alpha`, passed the
  largest double and became `-Inf` once `|log u|` exceeded `alpha` times
  the largest double, which at `2^-1022` happens for `u` below `e^-4`;
  the acceptance test then compared `-Inf` with `-Inf` and accepted
  every such draw. Such a draw is now tested at the largest negative
  double, where the ratio of the law to the envelope has already reached
  the value it keeps below `y = -745`, and is accepted with that
  probability, so that the acceptance rate is again the ratio of the two
  areas. Under R’s default generator, whose uniform variates are not
  below 2^-33, this can happen only below a shape of about 1.3e-307;
  under another generator it can reach somewhat larger shapes, and is
  handled the same way.

- Sun et al. Algorithm 3 searched for its matching point with an
  absolute stopping width, so once `beta` passed about 1e24 its
  acceptance rate collapsed and draws came back as `NaN`:
  `rmhn(2000, 50, 1e30, -3e15, method = "sun")` returned `NaN` for all
  but a handful of its draws, with a warning. The default method used
  Algorithm 3 for many of these laws, so it was affected too:
  `rmhn(500, 10, 1e30, -5e15)` returned `NaN` for about 70% of its
  draws. The search now runs at `beta = 1` and is rescaled.

- Algorithm 1 drew from too wide a law at a very large shape: its
  standard deviation was about 22% too large at `(1e20, 1, 1e10)` and
  46% at `(1e28, 1, 1e10)`, and a Kolmogorov-Smirnov test of 20,000
  draws against the normal limit gave `p` below 1e-40. The standard
  deviation of 20,000 draws is now within 1% of the limit’s. For `beta`
  above about 1e154 it returned only `NaN`, as in
  `rmhn(2000, 2, 1e160, 1e80)`.

- For small shapes a Gamma variate that underflowed became an exact
  zero: at `(0.01, 1, 0)` 2.4% of the draws were 0, where the law puts
  0.058% of its mass below the point at which a double rounds to 0, and
  at `alpha = 0.005` 15.7% against 2.41%. Such a variate is now replaced
  by a draw from the exact tail of the law.

- For `alpha < 1` with a strong positive tilt, RTDR built its envelope
  on an absolute log axis. From `gamma / sqrt(beta)` of about 1e13 its
  draws fell on a lattice: at `(0.7, 1, 1e14)` 20,000 draws took only
  about 30 distinct values, and a Kolmogorov-Smirnov test against the
  normal limit gave `p` below 1e-100. They now take some 570 values, the
  doubles `2^-7` apart near the mode, with the frequencies of the normal
  limit rounded to them. For `1/2 <= alpha < 1` it stopped from about
  3e14 with an error that called the parameter values a defect, and for
  `alpha < 1/2` from about 3.3e14. All regions now build and sample
  their pieces in the offset from the mode, and stop only where the law
  is narrower than the spacing of the doubles at its mode (see below).

- On the default path for `gamma < 0`, RTDR drew from the wrong law,
  without a warning, once the mode of `sqrt(beta) X` fell below about
  1e-27: `rmhn(20000, 2, 1e-56, -1)` had mean 1.17 where the law,
  Gamma(2, 1) to double precision, has mean 2 (Kolmogorov-Smirnov
  `p = 0`; the draws now have mean 2.0), and `rmhn(20000, 4, 1, -1e30)`
  stopped with an error that called itself a defect.

- For `alpha > 1` under a tilt `gamma / sqrt(beta)` near `-alpha / e`,
  where the mode of `sqrt(beta) X` lies near e, RTDR drew from the wrong
  law without a warning, or stopped with an error that called itself a
  defect, from a shape of about 3e15:
  `rmhn(20000, 1e17, 1, -1e17 / exp(1), method = "rtdr")` had a mean
  about 0.75 standard deviations below the law’s (Kolmogorov-Smirnov
  `p = 0` against Algorithm 3), and
  `rmhn(50, 3e17, 1, -3e17 / exp(1), method = "rtdr")` stopped. Its
  setup chose whether to build the envelope relative to the mode from
  the size of the log-density at the mode, which there is the small
  difference of terms of order `alpha` and carries their rounding error;
  it now chooses from the size of those terms. Over 81 parameter sets
  with shapes from 1e14 to 5e17 and tilts within 1% of `-alpha / e`,
  0.1.1 drew the wrong law in at least 21 and stopped in 8; RTDR’s draws
  now agree with those of Algorithm 3 in all 81 (two-sample
  Kolmogorov-Smirnov tests of 20,000 draws each, smallest `p` 0.0013).
  The default method of 0.1.1 used Algorithm 3 there; the default method
  now uses RTDR there at any number of draws (see Performance).

- `rmhn(20, 20, 1e-300, -1e300)` returned only `NaN`.

- Under `method = "rtdr"`, a law of tiny shape under a strong negative
  tilt, whose mass lies almost wholly below the smallest positive
  double, stopped with an error that called itself a defect, as in
  `rmhn(10, 1e-17, 1, -1e300, method = "rtdr")`. It is now drawn by
  Algorithm 3, which returns zeros there.

### Behaviour at the edges of the double range

- The help pages of
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) now
  state the range over which their results are accurate and what is
  returned beyond it, and that of
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) which
  laws double precision lets it sample and what it returns or stops with
  where it cannot. For
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) that
  range is `alpha` and `|gamma| / sqrt(beta)` up to 1e15 with `beta` in
  `[1e-300, 1e300]`. Beyond it the law can be narrower than the spacing
  of the doubles near its peak, and the probabilities are those of the
  exact law at the doubles given:
  `pmhn(mhn_mode(2, 1, 1e16), 2, 1, 1e16)` is now 0.5 (0.1.1: 0.585),
  and `pmhn(mhn_mode(2, 1e26, 1e17), 2, 1e26, 1e17)` 0.5 (0.1.1: 0).
  Under a nonzero tilt, from a shape of about 1e22, 0.1.1’s
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) also
  returned `NA` or stopped, with an error from the quadrature library, a
  length error or `std::bad_alloc`; at `gamma = 0` it did not stop.

- Past `|gamma| / sqrt(beta)` of about 1e154, where `gamma^2` overflows,
  the moment helpers stopped with an error from the quadrature library,
  and `mhn_mode(2, 1, 1e200)` was `Inf`. They now return the values of
  the limit laws: `mhn_mean(2, 1, 1e200)` and `mhn_var(2, 1, 1e200)` are
  5e199 and 0.5 (the normal limit), and `mhn_mean(2, 1, -1e200)` and
  `mhn_skewness(2, 1, -1e200)` are 2e-200 and 1.414 (the Gamma limit).

- At `alpha = 1` the density was `NaN` where `2 beta` overflows or the
  tilt is close to the largest double: `dmhn(1e-309, 1, 1, -1e308)` and
  `dmhn(0.5, 1, 1e308, 1e308)` are now 9.05e307 and 5.64e153.

- RTDR stops with an error, saying why, only where the standard
  deviation of the law is no larger than the spacing of the doubles at
  its mode: from `gamma / sqrt(beta)` of about 9e15 (`2^53` at
  `beta = 1`), and under `method = "rtdr"` with `gamma <= 0` from a
  shape between `2^103` and `2^106`. So `rmhn(3, 0.3, 1, 8e15)`, which
  0.1.1 refused, now samples, while
  `rmhn(3, 1.5e31, 1, 0, method = "rtdr")`, whose draws 0.1.1 centred on
  the mode rounded to a double (0.035 standard deviations from the true
  centre), now stops; the default method still samples it.

- `method = "sun"` now stops with an error for a negative tilt
  (`gamma < 0`, other than a negligible one, which is drawn as at
  `gamma = 0`) above a shape of `2^59` (about 5.8e17), where R’s
  [`rgamma()`](https://rdrr.io/r/stats/GammaDist.html) can shift the
  variates it supplies by more than 1e-7 of a standard deviation; the
  default method uses RTDR there. Above `alpha = 2^60` a negligible tilt
  is drawn from the Normal proposal of Algorithm 1 centred at the mode
  (by RTDR under `method = "rtdr"`) instead of as the square root of a
  Gamma variate.

- With `beta` far outside `[1e-300, 1e300]`, two intermediate quantities
  of the closed forms for the mode and for the points the distribution
  function starts from left the double range although the results do
  not. Where `alpha` and `beta` both exceed about 2.2e307,
  `sqrt(8 beta (alpha - 1))` overflows: `mhn_mode(1e308, 1e308, -1e308)`
  is now 0.5 (0.1.1: `NaN`), and
  `pmhn(0.9 * (1 + sqrt(5)) / 2, 1e308, 5e307, 1e308)`, a point below a
  law that sits on its mode 1.618, is 0 (0.1.1 stopped with an error
  from the quadrature library). Where `gamma^2` and `8 beta |alpha - 1|`
  both fall below the smallest normal double (for `alpha = 2`, a `beta`
  below about 2.8e-309 with `|gamma|` below 1.5e-154), they lose digits:
  `mhn_mode(2, 1e-322, 3e-161)` is now 1.799249e161 (0.1.1:
  1.799002e161), and
  `pmhn(0.9 / sqrt(1e-320), 150, 1e-320, 1e-170, log.p = TRUE)` is
  -268.49 (0.1.1: `NaN`). On 1287 parameter triples with `beta` in
  `[1e-300, 1e300]` every function and every seeded draw is unchanged by
  this, except
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  at one triple (`alpha = 1 + 1e-9`, `beta = 1e-300`), in its last
  digit.

### Interface

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) and
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) treat
  invalid arguments as the base R distribution functions do
  ([`?dgamma`](https://rdrr.io/r/stats/GammaDist.html): “Invalid
  arguments will result in return value NaN, with a warning”), element
  by element:

  - An element whose parameters lie outside the parameter space (`alpha`
    or `beta` not positive and finite, or `gamma` not finite) is `NaN`,
    with the warning “NaNs produced”
    ([`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md):
    “NAs produced”). 0.1.1 stopped the whole call, in
    [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) for a
    zero or negative `alpha` or `beta`.
  - An `NA` or `NaN` parameter gives `NA` or `NaN` without a warning
    from [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
    [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
    [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
    where 0.1.1 stopped the whole call, and `NaN` with the warning “NAs
    produced” from
    [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md),
    where 0.1.1 gave `NA` without one. An `NA` input gives `NA` and a
    `NaN` input `NaN`;
    [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) had
    relied on the exponential to keep an `NA` as `NA`, which R does not
    promise.
  - A zero-length `alpha`, `beta` or `gamma` gives `numeric(0)`, and in
    [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) `NA`
    draws with a warning, where 0.1.1 stopped.

  The other elements are unaffected, and an invalid element of
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) draws
  no random numbers. The moment functions and
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md),
  which take a single triple, still stop on an invalid one.

- A factor, a character string, a complex number, a list or `NULL` given
  to [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) or
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) is an
  error, as in [`dnorm()`](https://rdrr.io/r/stats/Normal.html); 0.1.1
  used a factor through its level codes, so that `dmhn(factor("3"))` was
  the density at 1, and a string through
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html). A Date or
  another vector stored as numbers is taken as its numbers, as before.
  The moment functions and
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  say which argument is not a single number, where 0.1.1 gave Rcpp’s
  “Expecting a single value” or “Not compatible with requested type”.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) now
  warns “NaNs produced” for a probability outside `[0, 1]`, as the base
  R quantile functions do.

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) keep
  the attributes of their first argument, such as `names` and `dim`, as
  the base R `d`/`p`/`q` functions do.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) uses
  only the elements of `alpha`, `beta` and `gamma` that its draws use,
  the first `n`, and ignores the rest, as
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) does:
  `rmhn(1, alpha = c(2, -1))` stopped with “alpha must be positive” and
  now draws from `alpha = 2`.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) reads
  `n` as [`rnorm()`](https://rdrr.io/r/stats/Normal.html) does: a
  zero-length `n` asks for no variates, where 0.1.1 stopped, and `NULL`
  and a negative fraction such as `-0.5` are errors, where 0.1.1
  returned `numeric(0)` for the latter. An `n` that is `NA`, `NaN`,
  infinite, beyond the integer range or a string that is not a number is
  reported as such rather than as negative, and one beyond the integer
  range without a coercion warning on the way.

- Messages from internal routines (`std::bad_alloc`, quadrature-library
  errors) no longer reach the user, and the errors of
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) say
  what went wrong at which parameter values instead of calling a
  representability limit a defect. The help page of
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)
  documents the `NaN` draws and the errors it can return.

- The warnings of the compiled code
  ([`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)) are
  raised through R’s [`warning()`](https://rdrr.io/r/base/warning.html)
  under Rcpp’s unwind protection (its default from Rcpp 1.0.10), and
  carry no call. A warning that is caught by
  `tryCatch(..., warning = )`, or turned into an error by
  `options(warn = 2)`, does not return to the compiled code; raised by
  `Rf_warning()`, as `Rcpp::warning()` raises it, it would leave by a
  long jump that skips the C++ destructors, and the vectors of the call
  would stay allocated for the rest of the session.

- mhn now requires Rcpp 1.0.13: `Imports` and `LinkingTo` ask for
  `Rcpp (>= 1.0.13)`, where 0.1.1 asked for 1.0.7 in `Imports` and for
  any version in `LinkingTo`. Rcpp turns on the unwind protection above
  by default from 1.0.10, and from 1.0.13 its `List::create()` takes
  more than 20 elements, as an internal routine of mhn needs; mhn does
  not compile against Rcpp 1.0.12.

### Performance

Measured on an Apple M4 Max against 0.1.1, as medians of repeated calls:

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) at the
  default parameters is about 4 times as fast (1e5 values in 0.98 ms
  against 4.3 ms), and 1.7 times as fast with recycled parameter
  vectors.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) for `q`
  above the median is 3 to 12 times as fast (5 against 17 microseconds
  per value at `(2.5, 1, -8)`), and about 160 times at `(2.5, 1, 30)`;
  for `|gamma| / sqrt(beta)` above 20 it is faster on both sides of the
  median (0.36 against 520 microseconds per value below it at
  `(2, 1, -50)`); and at a large shape `pmhn(71, 1e4, 1, 1)` takes a few
  microseconds against about 8 milliseconds.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) is 1.8
  times as fast at `(2.5, 1, 1.5)`, and at `alpha = 1` 5.6 times with
  recycled parameters and about 150 times for `qmhn(p, 1, 1, -5)`;
  `qmhn(0.5, 1.5, 1, 300)` takes about 30 microseconds against 0.2 s,
  and `qmhn(p, 2, 1, 2000)` about 20 microseconds per probability
  against 14 ms.

- [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  and
  [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
  at `(2, 1, 900)` take 5 microseconds together, against 17 ms.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) with
  short recycled parameter vectors builds one sampler setup per distinct
  triple instead of one per draw: `rmhn(1e4, c(2.5, 3), 1, -2)` takes
  0.54 ms against 3.4 ms. A single draw with `method = "sun"` or
  `method = "rtdr"` given exactly takes half as long, and a batch of
  draws from Algorithm 3 takes about 13% less time.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) with
  the default method chooses between RTDR and Sun et al.’s Algorithms 1
  and 3 by thresholds re-derived from timings of the current samplers
  (Apple M1 Ultra). For `gamma < 0` it uses RTDR from 25 draws per setup
  for `alpha < 0.1` (0.1.1: from 100), from 5 for `0.1 <= alpha < 1`
  (0.1.1: from 25), at any number of draws for `1 < alpha < 1000`, and
  from 10, 25 or 50 for `alpha` from 1000 to 1e5, from 1e5 to 1e8 and
  above (0.1.1: from 25, and never for `alpha >= 10`). RTDR’s setup
  walks its two contact points in from half and about 1.5 times the
  mode, which from a shape of about 1000 costs more than Algorithm 3’s
  setup; where it starts them next to the mode instead, under a very
  strong negative tilt (`gamma / sqrt(beta) < -1e10 * sqrt(alpha)`) or,
  at moderate tilts, at a shape beyond about 2.4e12 to 3.1e12, near a
  tilt of `-alpha / e` beyond about 2.3e13, and at any tilt beyond about
  4.5e13, the default method uses RTDR at any number of draws. For
  `gamma > 0` and `alpha > 1` it uses RTDR when
  `gamma / sqrt(beta) <= 2 * sqrt(alpha)` and at least 5, 25, 50 or 100
  draws (for a tilt up to `sqrt(alpha)`) or 25, 50, 100 or 1000 (above
  it) share a setup, for `alpha` below 10, from 10 to 1000, from 1000 to
  1e8 and above, or 5 or 25 where RTDR’s setup starts next to the mode,
  and Algorithm 1 otherwise (0.1.1: Algorithm 1 always). Over a paired
  timing race of the two samplers in 560 cells the chosen sampler is
  0.025% slower than the faster one on average and 7.1% at worst,
  against 6.7% and 72% under the rule of 0.1.1; over a second race of
  1458 cells, with shapes up to 1e14, tilts up to `3 * sqrt(alpha)` and
  single draws within one vectorised call (one setup per draw), it is
  0.11% slower on average and 13% at worst, and never slower in the
  vectorised calls, against 7.9% and 49%; in three runs of
  `inst/benchmarks/auto_dispatch.R` it is at most 3.8% slower, against
  66%. The draws per setup are counted as `n` for scalar parameters and,
  where short recycled parameter vectors keep a setup for each position
  within the period after which their triples repeat (a period of at
  most 256 draws, recycled at least four times), as `n` divided by that
  period, rounded down; otherwise as 1. 0.1.1 divided `n` by the length
  of the longest parameter vector, which with recycled parameters could
  count several draws per setup where every draw rebuilt one.

- Slower. The figures below are per value against 0.1.1, from three
  rounds in which the two versions ran side by side in persistent R
  sessions, alternating measure by measure, on a machine carrying other
  load; the ratios agree between rounds to within 10% at most triples
  and 20% at all.

  - `pmhn(q, lower.tail = FALSE, log.p = TRUE)` also evaluates the lower
    tail where the logarithm of the upper tail would lose its relative
    accuracy (see above): below `F` of about 1e-5 at parameters of order
    one, and of about 7e-4 at `(50, 1, 20)`. There it costs about what
    the lower tail costs, in either version. For `|gamma| / sqrt(beta)`
    up to 20 and `alpha` up to 100 that lower tail is the series of Sun
    et al., whose length grows like `gamma^2 / beta` with the tilt and
    also with the shape, and so does the factor. With `q` spread evenly
    in `log F` from 1e-6 to 1e-4 it takes 2.2 and 2.6 times as long as
    0.1.1 at `(0.5, 1, -5)` and `(2.5, 1, 1.5)`, 6.5 and 8.7 times at
    `(10, 1, 5)` and `(3, 1, 10)`, 11 and 15 times at `(0.9, 1, 12)` and
    `(5, 1, 13)`, 26 to 44 times for `gamma / sqrt(beta)` from 15 to
    19.5 at shapes from 0.01 to 10 (47 against 1.1 microseconds at
    `(0.5, 1, 19.5)`), and 45 to 77 times at shapes from 50 to 85 (90
    against 1.2 at `(85, 1, -19)`); with `q` spread evenly in `log F`
    from 1e-6 to 0.5, 1.5 to 23 times and 18 to 33 times. At the same
    `|gamma| / sqrt(beta)` a negative tilt costs less: 7 and 12 times at
    `(0.5, 1, -19.5)` and `(3, 1, -18)` from 1e-6 to 1e-4. Where the
    lower tail comes from a quadrature (`alpha > 100`, or
    `|gamma| / sqrt(beta) > 20` with `alpha >= 1`), both tails are
    quadratures, and from 1e-6 to 1e-4 it takes 1.1 to 1.9 times as long
    (1.8 against 1.0 microseconds at `(400, 1, -3)`, 1.5 times at
    `(1.5, 1, 25)` and 1.1 at `(5, 1, -40)`). For `alpha < 1` with a
    larger tilt, where the series would cost milliseconds and the lower
    tail is taken from its quadrature, checked against the upper tail,
    it takes 1.0 to 1.1 times as long at `(0.3, 1, -30)`,
    `(0.3, 1, -300)` and three triples with `alpha` from 0.16 to 0.39
    and `gamma / sqrt(beta)` from -120 to -1993, and less time under a
    strong positive tilt, where the normalising constant is now cheaper
    (0.9 to 3.4 against 13 to 25 microseconds at `(0.3, 1, 300)`). With
    `q` spread evenly in `F` below the median, where the upper tail
    alone serves, it takes as long as 0.1.1 to within 13% at all these
    triples.

  - `qmhn(p, lower.tail = FALSE)` for `p` close to 1 now works from the
    lower tail (see above). Where that is the series, which the Newton
    steps that refine the root sum a few times per value, it costs more:
    for 20 values of `1 - p` spread evenly in their logarithm, from
    1e-12 to 9e-8, from 1e-7 to 1e-5 or from 1e-5 to 7e-4, it takes 1.3
    to 2.6 times as long as 0.1.1 for `gamma / sqrt(beta)` from 15 to
    19.5 (118 against 46 microseconds at `(0.5, 1, 19.5)` from 1e-12 to
    9e-8), 2.0 to 4.0 times at shapes from 50 to 85 (236 against 59 at
    `(85, 1, -19)`), and 1.15 times at `(3, 1, -18)` from 1e-7 to 1e-5.
    At the other triples measured it takes at most 7% longer, and mostly
    less time: 74 times less under a strong positive tilt (92 against
    6800 microseconds at `(0.5, 1, 900)`) and 58 times at
    `(0.16, 1, -1993)` (16 against 964).

  - A single call of
    [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) with
    scalar arguments takes up to about 0.1 microseconds longer, 2 to 11%
    in the measurements, for keeping the attributes of `x` (see
    Interface).

  - A shape within `sqrt(.Machine$double.eps)` of 1 or a tilt below it,
    other than exactly 1 or 0, now costs what the general case costs
    next to it:
    [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) takes
    4.2 microseconds per value at `(1 + 1e-12, 1, -1)` and
    [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) 40,
    against 0.03 for both at `alpha = 1`.

  - For `|gamma| / sqrt(beta)` up to 20, the lower tail of
    [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) below
    the median, which is still summed from the series, takes up to 8%
    longer at some of these triples (8% at `(15.9, 1, -1.69)`, 7% at
    `(50, 1, -15)` and `(85, 1, -19)`), and up to a third less at
    others.

  - [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) for
    `gamma > 0` and `alpha > 1` with a parameter triple that changes at
    every draw (one Algorithm 1 setup per draw, as in a vectorised Gibbs
    update) takes 9 to 25% longer per draw than 0.1.1 at shapes from 2.5
    to 1e6, 7 to 14 nanoseconds more, 3 to 18% longer at shapes of 1.1
    and 1.5, and up to 13% longer at shapes from 1e8 to 1e12, where near
    a tilt `gamma / sqrt(beta)` of `sqrt(alpha)` it takes up to 20% less
    time instead (vectors of 2000 triples, 8 to 12 rounds). The setup of
    Algorithm 1 now forms the constants of the acceptance test of its
    proposal and checks that they stay within the double range.

  - Where 0.1.1 was wrong, the right answer can also cost more:
    `mhn_kurtosis(1000, 1, 0)` takes 11 microseconds against 6.

  - In the two timing races of the default method’s choice (see above),
    which time both samplers of this version, the sampler
    [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) now
    picks is at most 2.3% slower than the one the rule of 0.1.1 picks,
    except at `alpha = 1.5` with `gamma / sqrt(beta) = 2 * sqrt(alpha)`,
    where from 25 draws per setup it picks RTDR, 4.6 to 5.4% slower than
    Algorithm 1.

  Nothing else was measured to be slower beyond run-to-run noise: 20
  further workloads covering the density, both tails, the quantile, the
  moments and the samplers (single draws, batches of 1e5 from each
  sampler and a call with recycled parameters) took at most 7% longer.

### Documentation

- The NEWS entry of 0.1.1 attributed to Sun et al. (2023) a misprint in
  the constant `K_1` of the first part of Theorem 1 (printed “a.”), and
  the theory vignette of 0.1.1 one in the Gamma argument of the mixture
  weights of Lemma 7. The printed article contains neither: that part of
  Theorem 1 prints the numerator of the base of `K_1` as
  `sqrt(beta) (alpha - 1)` and Lemma 7 prints `Gamma(alpha/2 + i/2)`.
  The NEWS entry is corrected below, and the vignette now lists the six
  misprints the article does contain in the results it draws on that
  affect what is computed from them (it has other slips there, which do
  not); the package uses the correct form of each result it implements
  (the weights of Lemma 7 serve only Algorithm 2, which it does not
  implement):

  - Theorem 1 prints the density of the square root of a Gamma variable
    without its factor 2, so that it integrates to 1/2;
  - Theorem 1 prints the indicator that picks the proposal as
    `I(K_1(delta, ...) <= K_2(mu, ...))`, with the first arguments of
    `K_1` and `K_2` exchanged;
  - the weights of Lemma 7 carry a spurious factor
    `1 / (2 beta^(alpha/2))`, so they do not sum to one;
  - Supplementary Lemma 10 prints `alpha x` and `(alpha + 1) x` inside
    the discriminants of `C_1` and `C_2`, where the derivation gives
    `alpha x^2` and `(alpha + 1) x^2`;
  - step 2 of Algorithm 3 prints the Gamma shape as
    `(beta m + |gamma| alpha) / (2 beta m + |gamma|)`, where the display
    that introduces the proposal has
    `alpha (beta m + |gamma|) / (2 beta m + |gamma|)`;
  - the curvature `l''` of Supplementary Section 2.17 prints its first
    term without the minus sign that differentiating `l'` gives.

- The entry for 0.1.1 below is corrected where it did not describe 0.1.1
  or 0.1.0 as released: the counts of values and draws that changed, the
  scope of the RTDR overflow, the acceptance rate and the behaviour of
  RTDR at a steep tilt, four failures it attributed to 0.1.0 that 0.1.0
  did not show, the cost of the series at `gamma = -1e4`, the value
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) returns
  when it cannot bracket a root, the loops that check for an interrupt,
  and the size of the test suite.

- The help page of
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) said
  that a `-Inf` in `alpha` or `beta` gives an `NA` draw, whereas 0.1.1
  stopped the call for a zero or negative `alpha` or `beta`, `-Inf`
  included; 0.1.2 gives `NaN` with a warning (see Interface). The help
  page now also says how the number of variates per setup that steers
  the default method is counted, and that the switch to RTDR above
  `alpha = 2^59` concerns a negative tilt only.

- The help pages describe `beta` as the rate-type coefficient of
  `exp(-beta x^2)` and `gamma` as the coefficient of the tilt
  `exp(gamma x)` (0.1.1 called them a scale and a location parameter),
  and the theory vignette describes the normalising constant, the
  distribution function and the quantile as they are now computed.

- The help page of
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) gives
  the accuracy of the ordinary forms as measured: below 1e-8 within 10
  standard deviations of the peak, and about 3e-8 at 36 standard
  deviations in the upper tail close to the bound at which the forms
  measured from the peak take over. It names a third case outside the
  stated range, `beta` from about 2e307, in which a point far above the
  peak gives `NA` with a warning or `-Inf`.

- The help page of
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) says
  that a lower-tail `p` is solved on the series where
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) uses
  it, which for `alpha < 1` with `|z|` from about 100 can take from a
  tenth of a second to seconds for one quantile, and that `1 - p` with
  `lower.tail = FALSE` gives the same quantile at once. It says what is
  returned where the log upper tail of
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) is
  `-Inf` for `beta` from about 2e307.

- The help pages of
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) say
  when the normalising constant is reused across a vector of parameters:
  for each position of a recycling cycle of at most 256, and otherwise
  for each run of consecutive elements with the same parameters. The
  help page of
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) also
  says that `log = TRUE` returns `-Inf` in place of a finite log density
  for a positive tilt whose mode lies beyond the largest double.

- The help pages of
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
  and
  [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
  state the accuracy measured against independent quadratures (an
  absolute error below 1e-9 for shapes from 1 to 60 and tilts from -20
  to 12, the largest, 3e-10 and 6e-10, for a shape between 2 and 3 under
  a positive tilt), and how a subnormal shape is served.

- `citation("mhn")` gives the package’s DOI, 10.32614/CRAN.package.mhn,
  and the reference to Gao & Wang (2025) in the help pages, the README
  and the vignettes says it is an advance online publication, pages
  1-25.

### Tests

- The suite grows from about 9,300 expectations to about 15,500, of
  which about 12,400 run under `R CMD check`, where 26 heavier blocks,
  two timing tests, a memory test and a comparison to the last bit are
  skipped.

### Benchmark and audit scripts

- The result files that 0.1.1 installed under `benchmarks/results` and
  `audits/results` are no longer part of the package:
  `system.file("benchmarks", "results", package = "mhn")` and
  `system.file("audits", "results", package = "mhn")` are now `""`. The
  scripts under `benchmarks/` and `audits/` and the example under
  `examples/` remain.

- The benchmark and audit scripts and the example write their results
  only to the directory named by their output-directory variable
  (`MHN_BENCH_OUTDIR`, `MHN_GW_OUTDIR`, `MHN_GOF_OUTDIR`,
  `MHN_VMF_OUTDIR`), and otherwise only print a summary, where 0.1.1
  created a results directory inside the working directory. Next to the
  results they record the mhn version and source commit, R, the
  platform, CPU, compiler flags, BLAS and the run times
  (`benchmarks/provenance.R`). `auto_dispatch.R` takes its seed from
  `MHN_BENCH_SEED`; `replicate_gw_tables.R` prints the ratio of the two
  samplers’ times beside Gao & Wang’s ratio, their times not being
  comparable with the package’s; and `vmf_gibbs.R` no longer sums every
  autocorrelation into the effective sample size when the one at lag 1
  is negative.

- `series_breakdown_audit.R` forms the incomplete gamma values of its
  500-bit reference from their logarithms. 0.1.1 took them in double
  precision, where under a strong positive tilt they underflow to zero,
  so that the audit reported breakdowns of the double series that were
  errors of the reference. It also measures the error of a logarithm as
  the absolute difference, which is the relative error of the value.

## mhn 0.1.1

CRAN release: 2026-09-23

Bug-fix and maintenance release.

This release corrects a number of numerical defects reachable from the
public API. Most returned a wrong number silently rather than raising an
error, and most require parameter values outside the ranges 0.1.0 was
tested over. Each fix below carries a regression test at the values that
failed.

### What changes relative to 0.1.0

Numbers change wherever a defect is fixed, and that includes part of the
ordinary parameter range — not only the extremes. Everything below was
measured against 0.1.0 as published, over a grid of 330 parameter
triples spanning `alpha` in \[0.05, 20\], `beta` in \[0.01, 100\] and
`gamma` in \[-1e6, 1e6\].

**Deterministic functions.** Over the 130 grid cells with
`|gamma| / sqrt(beta)` at most 10,
[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
[`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
and [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
agree with 0.1.0 to twelve significant digits in 100 cells and differ in
30. The largest change is 8.2e-6 relative, in
[`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md);
[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) and
[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) move by
at most 3.1e-9 and
[`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
by 1.6e-7. In every case checked against an independent
arbitrary-precision reference the 0.1.1 value is the accurate one — the
Lemma 2c variance already cancels measurably at
`|gamma| / sqrt(beta) = 10`. Outside that range the differences are far
larger, and are the defects listed below.

**[`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) is not
seed-compatible with 0.1.0.** For a fixed seed, 163 of the 330 cells
return a different sequence of draws; 51 of those are in the ordinary
range above. Draws at `gamma = 0` are unchanged, since that case is a
closed-form transform of
[`rgamma()`](https://rdrr.io/r/stats/GammaDist.html). There are two
reasons a sequence changes, and they matter differently:

- Where 0.1.0 sampled the wrong law, the old draws are invalid and any
  analysis resting on them should be re-run. At
  `(alpha, beta, gamma) = (0.3, 1, 100)` a Kolmogorov-Smirnov test of
  the 0.1.0 sample against the true distribution gives `p = 0`; the same
  test on 0.1.1 gives `p = 0.51`.
- Elsewhere the law was already right and only the realisation moves,
  because [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)
  is a rejection sampler: correcting a constant in the envelope changes
  which variates are proposed and which are accepted, so the stream is
  consumed differently. At `(10, 1, 1)` every one of 200,000 draws from
  a fixed seed moves, and yet a Kolmogorov-Smirnov test against
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) gives
  `p = 0.78` on the 0.1.0 sample and `p = 0.70` on the 0.1.1 one, with
  both sample means within 0.05% of the true 2.4387. Analyses in this
  category remain valid; only bit-for-bit reproduction of a stored
  sample is lost.

### Bug fixes — random generation

- `rmhn(method = "rtdr")`, and hence the default `method = "auto"` where
  it routes there, drew biased samples for `alpha < 1` and `gamma > 0`
  outside the corner `alpha < 1/2`,
  `gamma / sqrt(beta) > 2 * (1 - sqrt(1 - 2 * alpha))`: on the log axis
  the density is only `T_{-1/2}`-concave in that region, not
  log-concave, so the previous log-tangent envelope did not dominate it
  and over-weighted large values. The envelope now follows Gao & Wang
  (2025, Section 3.2 and Appendix B), using a `T_{-1/2}`
  (inverse-square) tangent hat when `gamma > 0` and the log-tangent hat
  only when `gamma <= 0`; the inflection-point envelope used in that
  corner was not affected.

- RTDR drew badly biased samples wherever an exponential envelope piece
  rose by more than `t = 709` in log scale, because the piece was
  evaluated by forming `exp(t)` before dividing it out. For
  `alpha < 1/2` this is region D’s left-tangent piece once
  `gamma / sqrt(beta)` exceeds about 76 (`t = 1068` at 100), under both
  the default method and `method = "rtdr"`. For `alpha > 1` it is region
  A’s left-tail piece under `method = "rtdr"` once `gamma / sqrt(beta)`
  exceeds about 530. Once `alpha` exceeds about 5e4 it also happens at
  small and negative tilts, which the default method reached for
  `gamma < 0` whenever 25 or more draws shared one parameter set. The
  piece area became infinite, so piece selection always fell through to
  the last piece, and the widest piece proposed nothing inside its own
  support. At `(0.3, 1, 100)` a Kolmogorov-Smirnov test against
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) gave
  `p = 1.2e-121`. Both the area and the inverse CDF now factor the
  dominant exponential out instead of forming it.

- The RTDR envelope was decided by rounding once the standardised tilt
  `gamma / sqrt(beta)` passed about 1e8, and answered differently on
  different machines: at `(0.7, 1e-8, 1e8)` macOS returned 50 `NaN`
  while Linux and Windows raised an error, from the same input, and the
  outcome was not even monotone in `gamma`. Near the mode the density
  carries the constant `gamma_norm^2 / 4` — 2.5e23 at that tilt — while
  every quantity the envelope reads out of it is an `O(1)` difference of
  two such values, so the contact drop was quantised to multiples of
  3.4e7 against a target of `log 4`. All three envelope regions carried
  it. The ordinate is now measured from the mode, past a conditioning
  threshold, through an identity that only ever adds same-signed terms,
  and the accept/reject test is measured the same way — without that
  second half the sampler draws from its own hat.

  Fixing the arithmetic exposed three further defects that it does not
  itself cure. The contact searches started a fixed distance from the
  mode and needed more Newton steps than their budget allowed — 31
  against 30 at `gamma / sqrt(beta) = 1e10` — so they returned a contact
  point whose drop was hundreds instead of one; they now start where Gao
  & Wang (2025, Sec. 2) recommend starting their Eq. (8) iteration. The
  secant count `ceil(rho)` was narrowed into an `int`, which is
  undefined behaviour past a tilt of about 1.8e9. And the plateau piece
  was dropped without complaint whenever its edges crossed, which is how
  the sampler’s support came to exclude the mode.

  This fix leaves draws below `gamma / sqrt(beta) = 5e6` unchanged.
  Above it they can change, and the new ones are the correct ones:
  against the distribution’s exact Gaussian limit the Kolmogorov-Smirnov
  statistic is 0.004 where it had reached 1.000, and the sample standard
  deviation is within 0.5% of the truth where it had been 1593 times too
  large. The acceptance rate — which Gao & Wang bound below by `1/e` and
  which no goodness-of-fit test can see, because an envelope that is
  merely too wide still samples the right law — now holds steady from a
  tilt of 1e6 to 1e13: between 0.82 and 0.89 for `alpha > 1` and for
  `alpha < 1/2`, where it had fallen to 0.001, and about 0.76 for
  `1/2 <= alpha < 1` with a positive tilt, the rate of the
  inverse-square hat described above.

  Once the law is narrower than the spacing of the doubles at its mode,
  no sampler working in these coordinates can proceed, and RTDR stops
  with an error instead of returning draws that are mostly `NaN`. Short
  of that point this release does not yet sample every law correctly:
  for `alpha < 1` the draws fall on a coarse lattice from a tilt of
  about 1e13 (27 distinct values in 2000 draws at `(0.7, 1, 1e14)`)
  until the sampler stops, near 3.3e14 for `alpha < 1/2`, and for
  `1/2 <= alpha < 1` the error raised from about 3e14 calls the
  parameter values a defect.

- `rmhn(200, 1.5, 0.01, 3162.3)` returned 200 `NaN`, behind 200 warnings
  that the Algorithm 1 sampler had exhausted its retry budget. Sun et
  al. (2023) require the Algorithm 1 proposal rate to lie in `(0, beta)`
  (Theorem 1, first part printed “a.”), and the `delta_opt` printed in
  its second part “a.” does so for every `gamma > 0`. Evaluated as
  printed in double precision, however, it subtracts two quantities of
  size `gamma^2`, and for a large tilt the rounded result leaves the
  interval — reaching 0, then going negative, then snapping to `beta`. A
  negative rate made the branch test `NaN < 0` select the sqrt-Gamma
  proposal, and `rgamma` was handed a negative scale. Rationalising the
  expression removes every subtraction.

- At `alpha = 1` with `gamma = -1e8` the sampler returned 81572 negative
  values out of 2e5 — outside the `(0, Inf)` support — and collapsed the
  sample onto 13 distinct doubles. Robert (1995) proposes `y = a + e`
  with `a = -mu/sigma`, and the code returned `mu + sigma * y`; since
  `sigma * a = -mu` that is `sigma * e` with the mean removed
  analytically, but written literally it cancels. It now returns
  `sigma * e`.

- The RTDR region A left-contact search started at `max(m/2, 1e-6)`.
  Whenever the mode was below `2e-6` that floor put the start on the
  wrong side of it, and `rmhn(200, 4, 1, -1e7)` failed outright on the
  default path.

- The RTDR right-contact searches ran on the log axis, where the density
  carries `exp(2y)` and overflows past `y = 354`. The starting offset
  Gao & Wang (2025, Sec. 2) recommend for their Eq. (8) iteration lands
  beyond that when the mode is flat — at `alpha = 1e-5`, 447 above the
  mode under a negative tilt and 372 under a small positive one — and in
  region D, which started one unit above the mode, a single Newton step
  could land there. The search returned a non-finite slope, which the
  sign guard admitted because `NaN` fails every comparison, and the
  plateau piece was dropped with no error: the sampler’s support no
  longer contained the mode or the right tail. At `(3e-4, 1, 0.01)` it
  drew nothing above 0.005, where the true probability is 0.0015.

- Corrected the Algorithm 1 envelope constant `K_1`. Version 0.1.0 had
  implemented the numerator of its power base as `sqrt(beta (alpha-1))`,
  where the first part of Sun et al. (2023) Theorem 1 (printed “a.”)
  prints `sqrt(beta) (alpha-1)`. The sampled law was always correct —
  both proposals are valid envelopes and each acceptance test is derived
  independently — but the choice between them was not: `log K_1` came
  out short by `(alpha-1) log(alpha-1) / 2`, so for `alpha > 2` the
  Normal proposal was selected across most of the small-tilt half of the
  region. Over `alpha` in {4, 10} and `gamma` in {0.05, 0.5, 2, 20},
  with `beta = 1`, the worst acceptance probability was 0.69, below the
  0.8 that Theorem 2e guarantees for `alpha >= 4`; it is now 0.82.

- Sun Algorithm 1 formed `K_1 - K_2` as a difference of two large terms;
  cancelling them analytically leaves an expression with no subtraction
  of like-sized quantities.

- RTDR region A needs an interior mode in `x`, and at `alpha = 1`
  exactly there is none, so the envelope degenerated to a plateau of
  height `-Inf` spanning `[0, Inf]` and every proposal was rejected.
  `alpha = 1` now uses the region whose mode is strictly positive for
  every `alpha > 0`.

- The matching point Sun et al. (2023, Section 4.2) print for
  `alpha <= 1.1`, `alpha^2 / (1 + alpha)`, contains no `beta`; used as
  printed it was not scale-equivariant, and the acceptance rate of Sun
  Algorithm 3 drifted with the scale. It is now read as the value at
  `beta = 1` and divided by `sqrt(beta)`. When the inflection heuristic
  left the support, the setup discarded the inflection point entirely
  rather than using it, which is the limit of Sun’s own expression.

- `rmhn(method = "sun")` refused parameter values the sun path answers
  in closed form, because the pre-scan tested the raw `alpha` and
  `gamma` while the dispatcher routes on the scale-free predicates.

### Bug fixes — density, distribution, quantile and moments

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) and
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  returned 0 everywhere for `alpha` near 2 with a large tilt — the
  two-dimensional projected-normal case (`alpha = p = 2`) of Sun et
  al. (2023, Section 2.2). `Psi[1, z]` formed `exp(z^2/4)` in real
  space, which overflows at `z = gamma/sqrt(beta) > 53.28`, so the
  normalising constant became infinite. The term is now built as a
  logarithm.

- The normalising constant for `gamma < 0` degraded from `|z| ~ 100` and
  collapsed to `-Inf` from `|z| ~ 600`, taking everything downstream
  with it:
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md)
  returned `Inf`,
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  returned 1 everywhere,
  [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md)
  `NaN`. The quadrature’s upper limit extended to `|z|/2`, the peak of
  the same expression with the opposite sign, so the mass fell inside
  the first fraction of a percent of the interval. The limit now follows
  the decay the integrand actually has, and the range is split at the
  peak.

- For `gamma > 0` the normalising constant had only the Lemma 10 series,
  which raised an error once its truncation length passed a ceiling — on
  ordinary parameter sets, since that length grows like `z^2` and
  `(2.5, 1e-4, 20)` is already past it. The quadrature that serves
  `gamma < 0` now serves this too. Where the closed-form mean would
  difference two `log Psi` values of size `z^2/4`, it now comes from
  quadrature instead.

- The CDF quadrature integrated `[0, x]` in one piece, but for a strong
  tilt the mass lies within a unit or so of the peak, so a fixed-order
  rule returned 0 at the distribution’s own mode. In 0.1.0 that path was
  unreachable, because the series raised an error first; it became
  reachable once the series was given the quadrature fallback described
  above, and both are fixed here. The range is now anchored at the peak,
  split there, and trimmed at both ends.

- `pmhn(..., lower.tail = FALSE)` returned exactly 0, and `log.p = TRUE`
  returned `-Inf`, as soon as the survival probability fell below about
  1e-16, because it was derived as `1 - F`. The upper tail is now
  computed directly.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) at
  `alpha = 1` always computed the upper tail and derived the lower as
  `1 -` that, so `pmhn(3, 1, 1, 18)` gave 0 where the true value is
  1.076e-17. Whichever tail is smaller is now computed directly.

- `qmhn(1e-4, 0.1, 1, 1)` returned 1.4e-17, whose actual probability is
  0.019 — 23 orders of magnitude above the true 2.3e-40, silently. For
  `alpha < 1` the quantile of a small `p` sits near `p^(1/alpha)`, and
  the lower bracket bottomed out at 1.4e-17 after 30 halvings; the
  solver then returned that endpoint as though it were the root. The
  search now runs to the denormal floor.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  returned `+Inf` for every `p` at `alpha = 1` with
  `gamma <~ -11.7 sqrt(beta)`, where the truncation mass was formed as
  `1 - Phi(-mu/sigma)` in linear space and underflowed to exactly 0.

- [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md)
  returned -701.6 at `(3, 1, 1000)`, for a quantity that cannot fall
  below -2, and
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md)
  had lost its sign by `gamma = 200`. Both built central moments by
  differencing raw ones, which cancels once the standard deviation is
  small next to the mean. They are now integrated directly, in the
  variable centred on the peak so that no step forms a difference larger
  than the answer.

- `log Psi[1/2, z]` carried an absolute error of 0.73 at `z = -1e8`, in
  a quantity whose true value is -17.73: `z^2/4` and `log Phi(z/sqrt 2)`
  cancel. Substituting the scaled complementary error function cancels
  the exponential analytically.

- The mode formula cancels for `gamma < 0`; Gao & Wang (2025, Eq. 9)
  print the conjugate form alongside it, and their reference code
  selects between the two on the sign. `mhn_mode(1.1, 1, -1e8)` returned
  0 instead of 1e-9. Five sites used only the unstable form.

- Both places that asked whether the tilt is negligible compared
  `|gamma|` against `sqrt(.Machine$double.eps)` on an absolute scale,
  but the shape depends on `gamma` only through
  `Delta = gamma / sqrt(beta)`, since `sqrt(beta) X` is MHN(alpha, 1,
  Delta). At `beta = 4e-18` a `gamma` of 1e-8 is an ordinary
  `Delta = 5`, and every function returned the untilted sqrt-Gamma
  answer: `mhn_mean(1, 4e-18, 1e-8)` gave 2.82e8 against a true
  1.2503e9. Both tests now use `Delta`.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) with
  `gamma = -1e4` requested about 1.4e8 series terms, taking over ten
  seconds and more than a gigabyte of memory. The truncation length
  grows like `z^2` and had no ceiling; past one the series now returns
  the sentinel that already routes the caller to quadrature.

- The normalising constant for `gamma < 0` was also wrong by a factor
  that grew without bound as `alpha` fell — 1545 at `alpha = 1e-6`,
  still 0.17% at `alpha = 0.01`. A fraction 1/e of the mass of
  `x^(alpha-1)` on (0, 1) lies below `x = exp(-1/alpha)`, which for
  small `alpha` is below any representable double, so no quadrature can
  reach it. The first stretch of the range is now integrated in closed
  form. The distribution function and the central-moment quadrature had
  the same gap; it was invisible while the normalising constant shared
  it, and appeared as soon as that was fixed.

- `pmhn(..., lower.tail = FALSE)` returned exactly 0 for every `q` at
  small `alpha` with a strong negative tilt, and `qmhn` inverted that by
  returning its own initial bracket — the same number for every `p`,
  with no warning. The quadrature range was set from the curvature at
  the peak, which is the wrong scale when the range is anchored anywhere
  else: at `beta = 1e-4` with `gamma = -1000` the Gaussian width is 894
  while the kernel dies within 0.08. The range now follows the slope at
  the anchor. `qmhn` no longer returns an unjustified value: it warns,
  and returns `NaN`, when it cannot bracket the root.

- [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)
  multiplies any error in the mean by about `gamma^2/(2 beta)`: the
  Lemma 2c form is exact algebra, but its leading terms are each about
  `alpha/(2 beta)` while the answer is about `alpha/gamma^2`. With the
  degraded normalising constant it returned 3.6e-5 at `(1.5, 1, -300)`
  against a true 1.7e-5, and beyond `gamma / sqrt(beta)` of about -600
  it returned 0 for `alpha > 1` and `alpha/(2 beta)` for `alpha <= 1`.
  It is now used only where the measured cancellation is small enough,
  and the second central moment comes from quadrature elsewhere — the
  closed form remains the more accurate of the two wherever it is well
  conditioned.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md)
  implemented its special-case dispatch twice and the two disagreed. The
  vector-parameter path formed the truncated-normal upper tail as
  `1 - (1 - upper)`, which is exactly 0 once the upper tail falls below
  the rounding of 1: `pmhn(3, 1, 4, -5, lower.tail = FALSE)` was
  1.49e-23 element-wise and 0 vectorised.

### Interface

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) and
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) mapped
  a `NaN` input to `NA`, where
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) and the
  base R `d`/`p`/`q` family return `NaN`. All three now distinguish the
  two on every path.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md)
  documents an `NA` draw for a non-finite parameter but guarded only
  `gamma`; an infinite `alpha` or `beta` reached the sampler.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) took
  `n[1]` where every base R `r*` function reads a vector `n` as
  `length(n)`, so `rmhn(c(10, 20, 30))` drew 10 variates instead of 3.

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) tested
  `lower.tail` and `log.p` with
  [`isTRUE()`](https://rdrr.io/r/base/Logic.html), which is stricter
  than base R: `pmhn(1, 2, 1, 1, lower.tail = 1)` returned the upper
  tail, 0.5661, where `lower.tail = TRUE` gives 0.4339 — and base R’s
  `pnorm(1, lower.tail = 1)` agrees with `lower.tail = TRUE`. All three
  functions
  ([`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)) now
  read the flag with
  [`as.logical()`](https://rdrr.io/r/base/logical.html), which agrees
  with base R for `TRUE`, `FALSE`, `1` and `0`, and reject one that
  resolves to `NA`.

- Errors raised when the sampler cannot build an envelope named private
  C++ routines and the sign of an intermediate slope. They now report
  what failed, at which parameter values, and that it should be
  reported.

### Performance

- Twenty of the twenty-five compiled entry points never touch the RNG
  and were still generated with an `Rcpp::RNGScope`, paying
  `GetRNGstate`/`PutRNGstate` on every call to `dmhn`, `pmhn`, `qmhn`
  and the moment helpers. Draws are unchanged.

- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  evaluated the normalising constant a second time for a triple whose
  CDF state had just computed it.

- [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) skips
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) when `method`
  was not supplied; it accounted for about 55% of the cost of the
  single-variate call a Gibbs sweep makes.

- Piece selection in the RTDR sampler rebuilt a cumulative-area table on
  every proposal. The table is now built once per envelope and searched
  by bisection.

### Source

- Removed development scaffolding from the installed source: references
  to files not in the tarball, plan markers describing a half-finished
  file, and stale claims of mirroring an R implementation that no longer
  exists. Removed three unused functions and one unused struct field.

- Added `src/mhn_stable.h`, which holds the better-conditioned algebraic
  forms in one place rather than repeating them at the call sites that
  need them.

- Added interrupt checking to the driver loops of
  [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) and
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md), to the
  loop of [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md)
  over recycled parameters, and to the two rejection loops in the
  truncated-normal sampler; a long
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) call
  with recycled parameters was uninterruptible.

- Moved the retry-budget warnings out of the sampler frames.
  `Rcpp::warning` is a bare `Rf_warning`, which under
  `options(warn = 2)` leaves the frame by a long jump without running
  destructors while the cached envelope holds five vectors.

### Tests

- The suite grows from about 2,200 expectations to about 9,300, of which
  about 6,100 run under `R CMD check`, where 40 heavier blocks are
  skipped, with a regression block for each defect above at the
  parameter values that failed. New coverage: non-unit `beta` throughout
  the sampler grids and a direct check of scale equivariance; the
  contact-point search in isolation; the `NA`/`NaN` contract against
  base R; and the accuracy audit compared against a reference that
  shares no code with the package.

### DESCRIPTION

- Spelled out “Markov chain Monte Carlo” and “relaxed transformed
  density rejection method” in the `Description` field, following CRAN
  reviewer feedback on unexpanded acronyms.

### Benchmarks

- Refined the `method = "auto"` dispatch for `gamma < 0`. Previously any
  batch of 25 or more variates per setup used RTDR; benchmarking across
  three independent runs showed that for `alpha >= 10` the Sun et al.
  2023. Algorithm 3 has the lower per-proposal cost and wins in the
        batch regime too, so `auto` now keeps RTDR for large batches
        only when `alpha < 10`. This lowers the worst-case slowdown of
        `auto` relative to the per-cell optimum from about 11% to about
        4% while leaving the common Gibbs (single-variate) path
        unchanged. `inst/benchmarks/auto_dispatch.R` gained an
        `alpha < 1` grid, a setup/per-proposal cost decomposition, and a
        comparison of the shipped rule against the measured optimum.
- Raised the `gamma < 0` batch cutoff from 25 to 100 variates per setup
  when `alpha < 0.1`. The crossover between Algorithm 3 and RTDR moves
  to larger batches as the shape shrinks — it sits near 25 for `alpha`
  around 0.8 but near 100 by `alpha = 0.01` — so the old cutoff sent 25
  to 99 variates per setup to RTDR while Algorithm 3 was still up to 10%
  faster there.
- Fixed a unit double-conversion in `inst/benchmarks/auto_dispatch.R`
  that inflated the reported `median_us` / `iqr_us` times by a factor of
  about 1e6. The `method = "auto"` dispatch *decisions* are ratio-based
  and were unaffected, as is `rmhn(method = "auto")` itself.
- Added a goodness-of-fit benchmark, `inst/benchmarks/rmhn_gof.R`, that
  writes Kolmogorov-Smirnov statistics and sample-vs-theory moment
  summaries across the parameter grid to a CSV.
- The two timing benchmarks now emit a `_diagnostics_<date>.csv` with
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), hardware,
  and `mhn` version provenance, matching the audit scripts.

### Examples

- Added `inst/examples/vmf_gibbs.R`, a self-contained Gibbs sampler for
  the von Mises-Fisher concentration parameter whose full conditional is
  an MHN law. It is run with
  `source(system.file("examples", "vmf_gibbs.R", package = "mhn"))`,
  takes its sample size, chain length and true concentration from
  `MHN_VMF_*` environment variables, and reports whether the 95%
  credible interval covers the true concentration, and the effective
  sample size.

## mhn 0.1.0

CRAN release: 2026-05-27

Initial release.

### Distribution functions

- [`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md), and
  [`rmhn()`](https://t-momozaki.github.io/mhn/reference/rmhn.md) provide
  density, distribution, quantile, and random generation for the
  Modified Half-Normal (MHN) distribution of Sun, Kong & Pal (2023).
- All four functions are vectorised over both the evaluation argument
  and the parameters `alpha`, `beta`, `gamma`, following standard R
  recycling rules.
- A ParamCache reuses the Fox–Wright Psi normalising constant across
  consecutive elements that share an (`alpha`, `beta`, `gamma`) triple,
  so grouped inputs are evaluated significantly faster than calling the
  functions inside an R loop.

### Random generation

- `rmhn(..., method = "auto")` (default) routes each parameter triple to
  the cheapest provably-correct sampler: closed-form shortcuts for the
  special cases, Sun et al. (2023, Algorithms 1 and 3) where they win,
  and the Gao & Wang (2025) Relaxed Transformed Density Rejection (RTDR)
  sampler elsewhere.
- `method = "rtdr"` forces RTDR with its uniform 1/e acceptance bound.
- `method = "sun"` forces Sun Algorithm 1 (`gamma > 0, alpha > 1`) or
  Algorithm 3 (`gamma <= 0`); Sun Algorithm 2 is intentionally not
  implemented and an unsupported combination triggers a clear error.

### CDF and quantile

- [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) uses
  the Sun et al. (2023, Lemma 1b) series in log space, truncated at the
  Sun et al. (2023, Supplementary Lemma 10(d)) constructive bound K =
  max(K1, K2); the truncation residual is bounded by the user’s
  tolerance divided by `Psi`.
- For `gamma < 0` the series uses sign-separated log-sum-exp + log-
  diff-exp accumulation and a runtime cancellation guard derived from
  the floor set by double precision: when the relative cancellation loss
  would exceed the user’s tolerance,
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) falls
  back to a peak-normalised Boost.Math quadrature (Gauss-Kronrod for
  `alpha >= 1`, tanh-sinh for `alpha < 1`) of the unnormalised density.
- [`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) inverts
  [`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) via
  `boost::math::tools::toms748_solve` on the bracket
  `[sqrt(eps), E(X) + 8 sqrt(Var(X))]`, doubling the upper end as
  needed.

### Summary statistics

- [`mhn_mean()`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
  [`mhn_var()`](https://t-momozaki.github.io/mhn/reference/mhn_var.md),
  [`mhn_skewness()`](https://t-momozaki.github.io/mhn/reference/mhn_skewness.md),
  [`mhn_kurtosis()`](https://t-momozaki.github.io/mhn/reference/mhn_kurtosis.md),
  and
  [`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
  evaluate the closed-form / recurrence-based expressions from Sun et
  al. (2023, Lemmas 2 and 3).

### Tests and documentation

- testthat suite with \> 1,700 expectations covering goodness-of-fit
  (Kolmogorov-Smirnov), special-case identities, vectorised recycling,
  NA / NaN propagation, and the Sun / Gao & Wang acceptance bounds.
- [`vignette("introduction", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/introduction.md)
  walks through every exported function with runnable examples.
- [`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md)
  is the theoretical companion: it covers the MHN family and its special
  cases, the Fox–Wright Psi normalising constant, Algorithms 1 and 3 of
  Sun et al. (2023), the four-region Gao & Wang RTDR construction, and
  the `rmhn(method = "auto")` decision tree.
- `citation("mhn")` returns four `bibentry` objects: the package, the
  Sun et al. (2023) paper, the Gao & Wang (2025) paper, and
  Robert (1995) for the truncated-normal sampler used at `alpha = 1`.
