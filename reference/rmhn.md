# Random Generation from the Modified Half-Normal Distribution

Draws random variates from the Modified Half-Normal (MHN) distribution
with parameters `alpha`, `beta`, and `gamma`.

## Usage

``` r
rmhn(n, alpha = 1, beta = 1, gamma = 0, method = c("auto", "rtdr", "sun"))
```

## Arguments

- n:

  Non-negative integer giving the number of variates to draw. `n = 0`
  returns `numeric(0)`.

- alpha:

  Shape parameter (\\\alpha \> 0\\). Scalar or numeric vector. Default:
  1.

- beta:

  Scale parameter (\\\beta \> 0\\). Scalar or numeric vector. Default:
  1.

- gamma:

  Location parameter (\\\gamma \in R\\). Scalar or numeric vector.
  Default: 0.

- method:

  Sampling algorithm. One of `"auto"` (default), `"rtdr"`, or `"sun"`.
  See Details.

## Value

A numeric vector of length `n`. If any of `alpha`, `beta`, `gamma`
(after recycling to length `n`) is `NA` or non-finite (`Inf`, `-Inf`,
`NaN`), the corresponding output element is `NA`.

## Details

The MHN density is \$\$f(x \mid \alpha, \beta, \gamma) = \frac{2
\beta^{\alpha/2} x^{\alpha-1} \exp(-\beta x^2 + \gamma
x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]} \quad (x \> 0)\$\$ where
\\\Psi\[a, z\]\\ is the Fox-Wright Psi function. `rmhn` does not
evaluate \\\Psi\\; the rejection-sampling kernels cancel it out.

The default parameters `alpha = 1, beta = 1, gamma = 0` correspond to
the half-normal distribution \\\mathrm{HN}(1/\sqrt{2})\\.

The `method` argument selects the rejection sampler:

- `"auto"`: Special-case shortcuts when applicable (\\\gamma \approx 0\\
  -\> sqrt-Gamma, \\\alpha \approx 1\\ -\> truncated normal). Otherwise
  dispatches to RTDR (Gao & Wang, 2025).

- `"rtdr"`: Force the Relaxed Transformed Density Rejection method of
  Gao & Wang (2025). The acceptance probability is bounded below by
  \\1/e \approx 0.368\\ uniformly over the parameter space. Note: Gao &
  Wang (2025) use the parameterization \\(\lambda, \alpha, \beta)\\ with
  density proportional to \\x^{\lambda - 1} \exp(-\alpha x^2 - \beta
  x)\\; the mapping to the Sun et al. parameterization used here is
  \\\lambda \leftrightarrow \alpha\\, \\\alpha \leftrightarrow \beta\\,
  \\\beta \leftrightarrow -\gamma\\ (sign flip on the linear term).

- `"sun"`: Force the Sun et al. (2023) algorithms. Algorithm 1 is used
  when \\\gamma \> 0\\ and \\\alpha \> 1\\; Algorithm 3 is used when
  \\\gamma \le 0\\. The combination \\\alpha \< 1\\ with \\\gamma \> 0\\
  is unsupported and triggers an error.

Vector parameters are recycled to length `n` following standard R rules.
Trailing parameter elements beyond index `n - 1` are silently ignored,
matching the convention of `rnorm`.

Internally the setup state of the chosen sampler is reused as long as
consecutive \\(\alpha, \beta, \gamma)\\ triples are equal, so passing
parameters grouped by triple is faster than calling `rmhn` inside an R
loop.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.

Gao, F. & Wang, H.-B. (2025). Generating modified-half-normal random
variates by a relaxed transformed density rejection method.
*Communications in Statistics - Simulation and Computation*.

Robert, C. P. (1995). Simulation of truncated normal variables.
*Statistics and Computing*, 5(2), 121–125.

## See also

[`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md),
[`mhn_mean`](https://t-momozaki.github.io/mhn/reference/mhn_mean.md),
[`mhn_var`](https://t-momozaki.github.io/mhn/reference/mhn_var.md)

## Examples

``` r
set.seed(1)
rmhn(10, alpha = 2, beta = 1, gamma = 0.5)
#>  [1] 0.3277054 0.9876283 0.2750659 0.7977315 1.2334905 1.0218017 0.6370495
#>  [8] 0.5679219 0.2958218 1.9631443

# Vector parameters are recycled to length n.
set.seed(1)
rmhn(5, alpha = c(1, 2, 3, 4, 5))
#> [1] 0.4429697 0.7989284 0.5821857 1.4320656 2.2998471
```
