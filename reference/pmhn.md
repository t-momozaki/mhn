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

  Scale parameter (\\\beta \> 0\\). Scalar or numeric vector. Default:
  1.

- gamma:

  Location parameter (\\\gamma \in R\\). Scalar or numeric vector.
  Default: 0.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le q)\\,
  otherwise \\P(X \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.
  Default: `FALSE`.

## Value

A numeric vector. The output length equals
`max(length(q), length(alpha), length(beta), length(gamma))`; each input
is recycled to that length following standard R recycling rules. For
`q <= 0` the CDF is 0; for `q = Inf` it is 1.

## Details

The CDF is computed via the series representation \$\$F(x \mid \alpha,
\beta, \gamma) = \frac{1}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]}
\sum\_{i=0}^{\infty} \frac{z^i}{i!}\\ \Gamma(s_i)\\ P(s_i, \beta
x^2)\$\$ where \\z = \gamma/\sqrt{\beta}\\, \\s_i = (\alpha + i)/2\\,
and \\P(s, y)\\ is the regularized lower incomplete gamma function (Sun
et al., 2023, Lemma 1b; equivalent to the paper's form via the identity
\\\Gamma(s)\\ P(s, y) = \gamma(s, y)\\, where \\\gamma(s, y)\\ is the
lower incomplete gamma function used in the paper). The infinite sum is
truncated at the constructive bound \\K = \max\\K_1, K_2\\\\ from Sun et
al. (2023), Supplementary Lemma 10(d), which makes the truncation
residual bounded by the user's tolerance divided by \\\Psi\\. When
double-precision cancellation in the alternating-sign accumulator for
\\\gamma \< 0\\ would exceed that tolerance, the series is replaced by a
Gauss-Kronrod (or tanh-sinh for \\\alpha \< 1\\) numerical integration
of the density on \\\[0, q\]\\.

Special cases are detected and dispatched to standard R primitives:

- \\\gamma = 0\\: `pgamma(q^2, alpha/2, scale = 1/beta)`

- \\\alpha = 1\\: truncated-normal CDF via `pnorm`

When any of `alpha`, `beta`, `gamma` is a vector, the CDF is evaluated
element-wise. The Fox-Wright \\\Psi\\ normalizing constant is recomputed
only when consecutive elements present a different \\(\alpha, \beta,
\gamma)\\ triple, so passing grouped parameters is significantly faster
than calling `pmhn` inside an R loop.

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507–1536.

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
