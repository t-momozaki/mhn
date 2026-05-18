# mhn: The Modified Half-Normal Distribution

Provides density
([`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)),
distribution function
([`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)), quantile
function ([`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md)),
random generation
([`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)), and
moments / mode of the Modified Half-Normal (MHN) distribution.

## MHN Distribution

The MHN distribution has probability density function \$\$f(x \mid
\alpha, \beta, \gamma) = \frac{2\beta^{\alpha/2} x^{\alpha-1}
\exp(-\beta x^2 + \gamma x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]}\$\$
for \\x \> 0\\, where \\\Psi\[a, z\]\\ is the Fox-Wright Psi function.

## Parameters

- `alpha`:

  Shape parameter (\\\alpha \> 0\\). Controls the \\x^{\alpha-1}\\ term.

- `beta`:

  Scale (rate) parameter (\\\beta \> 0\\). Controls the \\\exp(-\beta
  x^2)\\ term.

- `gamma`:

  Location (skewness) parameter (\\\gamma \in R\\). Controls the
  \\\exp(\gamma x)\\ term.

## Special Cases

- \\\gamma = 0\\:

  Square-root Gamma distribution (\\X^2 \sim \textrm{Gamma}(\alpha/2,
  \beta)\\) (Sun et al., 2023, Lemma 6a).

- \\\alpha = 1\\:

  Truncated normal distribution \\\textrm{TN}(\gamma/(2\beta),
  1/\sqrt{2\beta}, 0, \infty)\\ (Sun et al., 2023, Lemma 6b).

- \\\alpha = 1, \gamma = 0\\:

  Half-normal distribution \\\textrm{HN}(1/\sqrt{2\beta})\\ (Sun et al.,
  2023, Lemma 6c).

## References

Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
distribution: Properties and an efficient sampling scheme.
*Communications in Statistics - Theory and Methods*, 52(5), 1507-1536.
[doi:10.1080/03610926.2021.1934700](https://doi.org/10.1080/03610926.2021.1934700)

Gao, F. & Wang, H.-B. (2025). Generating modified-half-normal random
variates by a relaxed transformed density rejection method.
*Communications in Statistics - Simulation and Computation*.
[doi:10.1080/03610918.2025.2524551](https://doi.org/10.1080/03610918.2025.2524551)

## See also

Useful links:

- <https://github.com/t-momozaki/mhn>

- <https://t-momozaki.github.io/mhn/>

- Report bugs at <https://github.com/t-momozaki/mhn/issues>

## Author

**Maintainer**: Tomotaka Momozaki <momozaki.stat@gmail.com>
