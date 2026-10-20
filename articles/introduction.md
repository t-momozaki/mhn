# Introduction to the mhn Package

The **mhn** package provides density, distribution, quantile, and random
generation functions for the Modified Half-Normal (MHN) distribution,
plus helpers for its moments and mode. This vignette walks through the
basic usage of every exported function.

``` r

library(mhn)
```

## The MHN distribution

The MHN(\\\alpha\\, \\\beta\\, \\\gamma\\) distribution has support on
\\(0, \infty)\\ and density

\\f(x \mid \alpha, \beta, \gamma) = \frac{2 \beta^{\alpha/2} \\
x^{\alpha-1} \\ \exp(-\beta x^2 + \gamma x)} {\Psi\[\alpha/2,\\
\gamma/\sqrt{\beta}\]}, \qquad x \> 0,\\

where \\\alpha \> 0\\, \\\beta \> 0\\, \\\gamma \in \mathbb{R}\\, and
\\\Psi\[a, z\]\\ is the Fox–Wright \\{}\_1\Psi_1\\ function used as the
normalising constant (Sun, Kong & Pal 2023). The three parameters
control the polynomial factor \\x^{\alpha-1}\\, the Gaussian tail
\\\exp(-\beta x^2)\\, and the exponential tilt \\\exp(\gamma x)\\
respectively.

## Density: `dmhn()`

[`dmhn()`](https://t-momozaki.github.io/mhn/reference/dmhn.md) evaluates
the density (or log-density) at one or more points:

``` r

dmhn(c(0.5, 1, 2), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.4702986 0.7325378 0.1982764
dmhn(1, alpha = 2, beta = 1, gamma = 1, log = TRUE)
#> [1] -0.3112403
```

It is vectorised over both `x` and the parameters, with standard R
recycling:

``` r

dmhn(c(0.5, 1, 1.5), alpha = c(1, 2), beta = 1, gamma = c(0, 1, -1))
#> [1] 0.87878258 0.73253783 0.04310111
```

``` r

x <- c(seq(0.001, 0.3, length.out = 60), seq(0.3, 5, length.out = 140))
plot(x, dmhn(x, alpha = 2, beta = 1, gamma = 1),
     type = "l", lwd = 2, ylab = "density", ylim = c(0, 1),
     main = expression("MHN(" * alpha * ", " * beta * ", " * gamma * ")"))
lines(x, dmhn(x, alpha = 0.5, beta = 1, gamma = 1), lwd = 2, col = "tomato")
lines(x, dmhn(x, alpha = 0.3, beta = 1, gamma = 4), lwd = 2, col = "steelblue")
legend("topright", bty = "n",
       legend = c("(2, 1, 1)    log-concave, mode near 1",
                  "(0.5, 1, 1)  monotone decreasing",
                  "(0.3, 1, 4)  boundary spike + interior bump"),
       col = c("black", "tomato", "steelblue"), lwd = 2)
```

![MHN densities for three parameter triples; the steelblue curve sits in
the alpha \> 0 regime where the density combines a boundary divergence
at x to 0+ with an interior local maximum near x = 1.8 (Sun et al. 2023,
Lemma 3c). The y-axis is clipped at 1; the divergent left tails of the
tomato and steelblue curves continue upward beyond the
plot.](introduction_files/figure-html/dmhn-plot-1.png)

MHN densities for three parameter triples; the steelblue curve sits in
the alpha \< 1, gamma \>\> 0 regime where the density combines a
boundary divergence at x to 0+ with an interior local maximum near x =
1.8 (Sun et al. 2023, Lemma 3c). The y-axis is clipped at 1; the
divergent left tails of the tomato and steelblue curves continue upward
beyond the plot.

## Distribution function: `pmhn()`

[`pmhn()`](https://t-momozaki.github.io/mhn/reference/pmhn.md) returns
\\P(X \le q)\\ (or its complement / log-probability via `lower.tail` /
`log.p`). In the general case each tail is computed directly on its own
side of the peak of the density of \\\log X\\, which lies at or above
the median.

Below that peak the lower tail comes from one of two routes, chosen by
measuring their cost and accuracy. One is the Sun et al. (2023) Lemma 1b
series, summed in log space and truncated at an upper bound on the
length \\K = \max\\K_1, K_2\\\\ of their Supplementary Lemma 10(d), from
the geometric decay of the terms that Lemma 10(a, b) establish. The
other is a Boost.Math quadrature of the unnormalised density:
Gauss–Kronrod, and for \\\alpha \< 1\\, where the density is unbounded
at the origin, a closed form next to the origin, then tanh–sinh up to
\\q\\ or up to the interior peak that a positive tilt can give the
density, whichever comes first, and Gauss–Kronrod past that peak. The
series is tried first for \\\alpha \< 1\\ up to \\\gamma/\sqrt{\beta} =
10^3\\, where the quadrature misses the working tolerance
`sqrt(.Machine$double.eps)` at scattered points, and for \\\alpha \ge
1\\ while \\\|\gamma\|/\sqrt{\beta} \le 20\\ and \\\alpha \le 100\\.
Where it is tried, the quadrature still replaces it when, for \\\gamma
\< 0\\, the cancellation between its alternating terms could push its
error past that tolerance, when it would need more than \\10^6\\ terms,
and when \\\beta q^2\\ is below the smallest normal double. Elsewhere
the quadrature is used directly: there it is the faster of the two, and
within the tolerance.

The upper tail \\P(X \> q)\\ is computed directly as well, by quadrature
on a logarithmic axis, rather than as \\1 - F\\, which would round to 0
once it falls below about \\10^{-16}\\; above the peak, where it is the
smaller tail, the lower tail is returned as one minus it. With
`log.p = TRUE` both tails stay on the log scale, with two exceptions. At
or below the median, where \\F \le 1/2\\, \\\log P(X \> q)\\ is about
\\-F\\, and the logarithm of the quadrature carries an absolute error of
about \\\epsilon\\ times the size of the terms it is assembled from.
Wherever that could pass \\10^{-9}\\ of \\F\\, below \\F\\ of about
\\10^{-5}\\ for parameters of order one, \\\log P(X \> q)\\ is formed
instead as \\\log(1 - F)\\ from \\\log F\\, which keeps the relative
accuracy of a small \\F\\. Likewise, below the peak of a small shape,
where \\F\\ is close to 1 and \\\log F\\ is about \\-P(X \> q)\\, \\\log
F\\ is formed as \\\log(1 - P(X \> q))\\ from the upper tail wherever
that is the more accurate; without that switch its error would pass
\\10^{-9}\\ for a shape below about \\10^{-6}\\ under a small tilt. At a
very large shape or tilt both tails are measured from the peak instead
(see [`?pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)).

``` r

pmhn(c(0.5, 1, 1.5, 2), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.1129101 0.4338796 0.7614259 0.9363038
pmhn(2, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#> [1] 0.06369619
pmhn(2, alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#> [1] -0.06581527
```

A direct cross-check against `integrate(dmhn, 0, q)` confirms accuracy:

``` r

q <- 1.5
ref <- integrate(function(x) dmhn(x, 2, 1, 1), 0, q,
                 rel.tol = 1e-10)$value
all.equal(pmhn(q, 2, 1, 1), ref)
#> [1] TRUE
```

``` r

plot(x, dmhn(x, 2, 1, 1), type = "l", lwd = 2, ylab = "f(x)", main = "Density")
plot(x, pmhn(x, 2, 1, 1), type = "l", lwd = 2, ylab = "F(x)", main = "CDF",
     ylim = c(0, 1))
```

![Density and matching CDF for MHN(2, 1,
1).](introduction_files/figure-html/pmhn-plot-1.png)

Density and matching CDF for MHN(2, 1, 1).

## Quantile function: `qmhn()`

[`qmhn()`](https://t-momozaki.github.io/mhn/reference/qmhn.md) inverts
the CDF. Write \\p\_\ell\\ and \\p_u = 1 - p\_\ell\\ for the lower- and
upper-tail probabilities: the argument `p` (its logarithm with
`log.p = TRUE`) is \\p\_\ell\\ with `lower.tail = TRUE` and \\p_u\\ with
`lower.tail = FALSE`, and the other one is formed from it without
cancellation where it is needed. In the general case a TOMS 748
root-finder solves \\F(x) = p\_\ell\\, or \\\log P(X \> x) = \log p_u\\
on the directly computed upper tail where that tail resolves the
probability better: for an upper-tail `p` whose \\p\_\ell\\ that
logarithm resolves to \\10^{-9}\\ of itself (at least \\3 \times
10^{-6}\\ or so from 1, more for a large shape, tilt or
\\\|\log\beta\|\\), and for a lower-tail `p` close enough to 1 that
\\F\\ no longer resolves \\p_u\\ to \\10^{-9}\\ of itself (within about
\\10^{-5}\\ of 1 for parameters of order one, farther for a small or
large shape, a large tilt or an extreme \\\beta\\; see
[`?qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md)). Closer
to 1 an upper-tail `p` is solved as \\F(x) = p\_\ell\\, or, where \\F\\
comes from the series and costs far more, its upper-tail root is refined
by Newton steps on \\\log F\\; within \\10^{-7}\\ of 1 that root is
taken, where needed, for a \\p\_\ell\\ the upper tail still resolves.
With `log.p = TRUE`, a lower-tail probability whose logarithm is below
`log(.Machine$double.xmin)`, about \\-708.4\\, is solved as \\\log F(x)
= \log p\_\ell\\ instead. The initial bracket is \\\[s\sqrt{\epsilon},\\
s\]\\ with \\s = E(X) + 8 \sqrt{\mathrm{Var}(X)}\\ and \\\epsilon\\ the
machine epsilon; its left end is halved and its right end doubled until
it brackets the target. Where the root-finder does not reach its
tolerance on that bracket, as for a tiny lower-tail `p` at a large
shape, it solves for \\\log x\\ instead, on \\\log F\\ for a lower-tail
equation, and failing that the bracket is bisected. The special cases
invert the laws they reduce to (see
[`?qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md)).

``` r

qmhn(c(0.1, 0.5, 0.9), alpha = 2, beta = 1, gamma = 1)
#> [1] 0.4717434 1.0906276 1.8480472
```

The round-trip identity holds within the inverter’s tolerance:

``` r

p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
all.equal(pmhn(qmhn(p, 2, 1, 1), 2, 1, 1), p, tolerance = 1e-6)
#> [1] TRUE
```

`lower.tail` and `log.p` follow the same conventions as
`qnorm`/`qgamma`:

``` r

qmhn(0.95, 2, 1, 1, lower.tail = FALSE)        # = qmhn(0.05)
#> [1] 0.3394014
qmhn(log(0.05), 2, 1, 1, log.p = TRUE)         # same value, log-input
#> [1] 0.3394014
```

## Random generation: `rmhn()`

`rmhn(n, alpha, beta, gamma)` draws `n` variates. The default
`method = "auto"` chooses between the special-case shortcuts (sqrt-Gamma
for \\\gamma = 0\\; truncated normal for \\\alpha = 1\\), Sun et al.
(2023) Algorithms 1 / 3, and the Gao & Wang (2025) Relaxed Transformed
Density Rejection (RTDR) sampler. `method = "rtdr"` or `method = "sun"`
selects which of the two general-case samplers is used (with a few
exceptions at the limits of double precision, listed in
[`?rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)); the
special-case shortcuts are taken first in either case, being exact and
needing next to no setup. The default method chooses its sampler from
the number of variates per setup, so under one seed its stream can
change with `n` where that choice does (at `(2, 1, 1)` the first four
draws of `rmhn(100, ...)` are not those of `rmhn(4, ...)`);
`method = "rtdr"` or `"sun"` keeps one stream whatever `n` (see
[`?rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)).

``` r

rmhn(5, alpha = 2, beta = 1, gamma = 1)
#> [1] 0.5497233 1.9097224 0.3042778 1.0601083 1.1484879

# Vector parameters are recycled to length n.
rmhn(5, alpha = c(1, 2, 3, 4, 5))
#> [1] 0.5649849 0.1332783 1.1693196 1.3508566 0.9682530
```

``` r

set.seed(42)
draws <- rmhn(10000, alpha = 2, beta = 1, gamma = 1)
hist(draws, breaks = 60, probability = TRUE, col = "grey90", border = "white",
     xlab = "x", main = "rmhn(10000, 2, 1, 1)")
lines(x, dmhn(x, 2, 1, 1), lwd = 2, col = "tomato")
```

![10,000 draws from MHN(2, 1, 1) with the true density
overlaid.](introduction_files/figure-html/rmhn-overlay-1.png)

10,000 draws from MHN(2, 1, 1) with the true density overlaid.

Switching `method` is useful for benchmarking. Both samplers draw the
exact law, so a sample from either is consistent with it; here two
samples, drawn one after the other from a single seed, are each compared
with the distribution function by a one-sample Kolmogorov-Smirnov test:

``` r

set.seed(1)
s_rtdr <- rmhn(5000, 2, 1, 1, method = "rtdr")
s_sun  <- rmhn(5000, 2, 1, 1, method = "sun")
F <- function(q) pmhn(q, 2, 1, 1)
c(rtdr = ks.test(s_rtdr, F)$p.value, sun = ks.test(s_sun, F)$p.value)
#>      rtdr       sun 
#> 0.3785647 0.6279249
```

## Moments and mode

The package provides helpers for the common summary statistics. The mean
is the Fox–Wright ratio of Sun et al. (2023, Lemma 2a) and the variance
follows from it by their Lemma 2c wherever these closed forms are
accurate; elsewhere, and for the skewness and kurtosis, the moments are
integrated numerically against the density. The mode is the closed form
of Sun et al. (2023, Lemma 3).

``` r

data.frame(
  Quantity = c("mean", "variance", "skewness", "excess kurtosis", "mode"),
  Function = c("mhn_mean", "mhn_var", "mhn_skewness",
               "mhn_kurtosis", "mhn_mode"),
  Value = c(
    mhn_mean(2, 1, 1),
    mhn_var(2, 1, 1),
    mhn_skewness(2, 1, 1),
    mhn_kurtosis(2, 1, 1),
    mhn_mode(2, 1, 1)
  )
)
#>          Quantity     Function       Value
#> 1            mean     mhn_mean 1.133731086
#> 2        variance      mhn_var 0.281519367
#> 3        skewness mhn_skewness 0.463887428
#> 4 excess kurtosis mhn_kurtosis 0.006868473
#> 5            mode     mhn_mode 1.000000000
```

For \\\alpha \< 1\\ without an interior local maximum (for instance with
\\\gamma \le 0\\),
[`mhn_mode()`](https://t-momozaki.github.io/mhn/reference/mhn_mode.md)
returns `NA`; at \\\alpha = 1\\ with \\\gamma \le 0\\ it returns the
boundary mode 0:

``` r

mhn_mode(0.5, 1, -1)
#> [1] NA
mhn_mode(1, 1, -1)
#> [1] 0
```

## Special cases

The MHN family contains several familiar distributions:

| Constraint | Reduction |
|----|----|
| \\\gamma = 0\\ | \\X^2 \sim \mathrm{Gamma}(\alpha/2, \beta)\\, rate \\\beta\\ (sqrt-Gamma) |
| \\\alpha = 1\\ | Truncated normal on \\(0, \infty)\\ |
| \\\alpha = 1, \gamma = 0\\ | Half-normal \\\|Z\|, Z \sim N(0, 1/(2\beta))\\ |

At \\\gamma = 0\\ and at \\\alpha = 1\\ exactly, `dmhn` / `pmhn` /
`qmhn` work through the corresponding closed form, the Gamma law of
\\\beta X^2\\ or the truncated normal, rather than the general series; a
shape or a tilt merely close to those values is a different law, and
goes through the general route, which keeps its accuracy there. `rmhn`
draws from the closed forms also for \\\|\alpha - 1\|\\ or
\\\|\gamma\|/\sqrt{\beta}\\ below `sqrt(.Machine$double.eps)`, a law
within a total variation distance of \\8.6 \times 10^{-9}\\ of the one
requested. A zero tilt goes through the general route above a shape of
\\10^6\\ (\\2^{60}\\ for `rmhn`), where the Gamma law is too narrow for
the rounding of \\\beta x^2\\ or of `rgamma`.

``` r

# gamma = 0: dmhn matches the change-of-variable sqrt-Gamma density.
xx <- c(0.5, 1, 1.5, 2)
mhn_d  <- dmhn(xx, alpha = 2, beta = 1, gamma = 0)
ref_d  <- dgamma(xx^2, shape = 1, rate = 1) * 2 * xx
all.equal(mhn_d, ref_d)
#> [1] TRUE
```

``` r

mu_tn    <- 1 / (2 * 1)
sigma_tn <- 1 / sqrt(2 * 1)
tn_dens  <- function(x) {
  dnorm(x, mu_tn, sigma_tn) / pnorm(0, mu_tn, sigma_tn, lower.tail = FALSE)
}
plot(x, dmhn(x, 1, 1, 1), type = "l", lwd = 2, ylab = "density",
     main = "alpha = 1: MHN reduces to truncated normal")
points(x, tn_dens(x), pch = ".", col = "tomato", cex = 2)
legend("topright", bty = "n",
       legend = c("dmhn(x, 1, 1, 1)", "TN reference"),
       col = c("black", "tomato"), lwd = c(2, NA), pch = c(NA, 19))
```

![MHN(1, 1, 1) and the equivalent truncated normal
density.](introduction_files/figure-html/special-overlay-1.png)

MHN(1, 1, 1) and the equivalent truncated normal density.

## Further reading

- See [`?dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
  [`?pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md),
  [`?qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md),
  [`?rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md) for full
  argument documentation, including the recycling rules and the `method`
  argument of `rmhn`.
- `citation("mhn")` returns four entries: the package itself and the
  three underlying papers.
- Sun, Kong & Pal (2023) develop the parametric family and Algorithms 1
  and 3 used here. Gao & Wang (2025) introduce the RTDR sampler, whose
  acceptance they bound below by \\1/e\\ where its contact points lie in
  a band they give.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1591–1613.
<https://doi.org/10.1080/03610926.2021.1934700>

Gao, F., & Wang, H.-B. (2025). Generating modified-half-normal random
variates by a relaxed transformed density rejection method.
*Communications in Statistics - Simulation and Computation*, 1–25.
Advance online publication.
<https://doi.org/10.1080/03610918.2025.2524551>

Robert, C. P. (1995). Simulation of truncated normal variables.
*Statistics and Computing*, 5(2), 121–125.
<https://doi.org/10.1007/BF00143942>
