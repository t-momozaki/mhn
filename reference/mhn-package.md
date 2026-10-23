# mhn: The Modified Half-Normal Distribution

Provides density
([`dmhn`](https://t-momozaki.github.io/mhn/reference/dmhn.md)),
distribution function
([`pmhn`](https://t-momozaki.github.io/mhn/reference/pmhn.md)), quantile
function ([`qmhn`](https://t-momozaki.github.io/mhn/reference/qmhn.md)),
random generation
([`rmhn`](https://t-momozaki.github.io/mhn/reference/rmhn.md)), and
moments / mode of the Modified Half-Normal (MHN) distribution.

## Details

[`vignette("introduction", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/introduction.md)
shows the functions in use, and
[`vignette("theory", package = "mhn")`](https://t-momozaki.github.io/mhn/articles/theory.md)
sets out the theory they implement; its Appendix A describes the
numerical methods.

## MHN Distribution

The MHN distribution has probability density function \$\$f(x \mid
\alpha, \beta, \gamma) = \frac{2\beta^{\alpha/2} x^{\alpha-1}
\exp(-\beta x^2 + \gamma x)}{\Psi\[\alpha/2, \gamma/\sqrt{\beta}\]}\$\$
for \\x \> 0\\, where \\\Psi\[a, z\] = \sum\_{n \ge 0} \Gamma(a + n/2)
z^n / n!\\ is the Fox-Wright function \\{}\_1\Psi_1\[(a, 1/2); (1, 0);
z\]\\.

## Parameters

- `alpha`:

  Shape parameter (\\\alpha \> 0\\). Controls the \\x^{\alpha-1}\\ term.

- `beta`:

  Rate-type parameter (\\\beta \> 0\\). Controls the Gaussian tail
  \\\exp(-\beta x^2)\\; larger values pull mass toward zero.

- `gamma`:

  Tilt parameter (\\\gamma \in R\\). Controls the exponential tilt
  \\\exp(\gamma x)\\; larger values move mass away from zero, but
  \\\gamma\\ does not translate the law, whose support stays \\(0,
  \infty)\\.

## Special Cases

- \\\gamma = 0\\:

  Square-root Gamma distribution (\\X^2 \sim \textrm{Gamma}(\alpha/2,
  \beta)\\, with rate \\\beta\\) (Sun et al. 2023, Lemma 6a).

- \\\alpha = 1\\:

  Truncated normal distribution \\\textrm{TN}(\gamma/(2\beta),
  1/\sqrt{2\beta}, 0, \infty)\\ (Sun et al. 2023, Lemma 6b).

- \\\alpha = 1, \gamma = 0\\:

  Half-normal distribution \\\textrm{HN}(1/\sqrt{2\beta})\\ (Sun et al.
  2023, Lemma 6c).

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

Useful links:

- <https://github.com/t-momozaki/mhn>

- <https://t-momozaki.github.io/mhn/>

- Report bugs at <https://github.com/t-momozaki/mhn/issues>

## Author

**Maintainer**: Tomotaka Momozaki <momozaki.stat@gmail.com>
([ORCID](https://orcid.org/0000-0001-6735-0804))
