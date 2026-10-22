# Moment functions for the Modified Half-Normal distribution.
# Thin wrappers around .mhn_mean_cpp / .mhn_var_cpp / .mhn_skewness_cpp /
# .mhn_kurtosis_cpp (src/mhn_moments.cpp).

#' Mean of the Modified Half-Normal Distribution
#'
#' Computes \eqn{E(X)} for
#' \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}{X ~ MHN(alpha, beta,
#' gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}), a single number.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. A single number.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   A single number.
#'
#' @return A numeric scalar.
#'
#' @details
#' By Sun et al. (2023, Lemma 2a) the mean is a ratio of Fox-Wright
#' \eqn{\Psi}{Psi} functions,
#' \deqn{E(X) = \frac{\Psi[(\alpha+1)/2,\, z]}{\sqrt{\beta}\,
#'   \Psi[\alpha/2,\, z]}, \qquad z = \gamma/\sqrt{\beta}.}{%
#'   E(X) = Psi[(alpha+1)/2, z] / (sqrt(beta) Psi[alpha/2, z]),
#'   z = gamma/sqrt(beta).}
#' At \eqn{\gamma = 0} the ratio is
#' \eqn{\Gamma((\alpha+1)/2) / (\sqrt{\beta}\,
#' \Gamma(\alpha/2))}{Gamma((alpha+1)/2) / (sqrt(beta) Gamma(alpha/2))}
#' (Sun et al. 2023, Lemma 6a). For every \eqn{\alpha > 0}{alpha > 0} the
#' mean is at most
#' \eqn{(\gamma + \sqrt{\gamma^2 + 8\alpha\beta})/(4\beta)}{(gamma +
#' sqrt(gamma^2 + 8 alpha beta)) / (4 beta)} (Sun et al. 2023, Lemma 4a and
#' the remark after it).
#'
#' The ratio is used where it can be evaluated accurately. Where it cannot,
#' mainly under a strong tilt or at a very large \eqn{\alpha}{alpha}, the
#' mean is instead obtained by numerical integration of the unnormalised
#' density, which needs no \eqn{\Psi}{Psi} value.
#'
#' The function is not vectorised in its parameters: \code{alpha},
#' \code{beta} and \code{gamma} must each be a single number. To evaluate
#' several parameter triples, use \code{\link{mapply}}, as in the examples.
#'
#' The numerical methods are described in Appendix A of
#' \code{vignette("theory", package = "mhn")}.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. (Lemmas 2a, 4a and 6a)
#'
#' @seealso \code{\link{mhn_var}}, \code{\link{dmhn}}
#'
#' @examples
#' mhn_mean(alpha = 2, beta = 1, gamma = 0)
#' mapply(mhn_mean, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#'
#' @export
mhn_mean <- function(alpha, beta, gamma) {
  .check_single(alpha, beta, gamma)
  .mhn_mean_cpp(alpha, beta, gamma)
}


#' Variance of the Modified Half-Normal Distribution
#'
#' Computes \eqn{\mathrm{Var}(X)}{Var(X)} for
#' \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}{X ~ MHN(alpha, beta,
#' gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}), a single number.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. A single number.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   A single number.
#'
#' @return A numeric scalar.
#'
#' @details
#' Uses the formula (Sun et al. 2023, Lemma 2c)
#' \deqn{\mathrm{Var}(X) = \frac{\alpha}{2\beta} +
#'   E(X)\left(\frac{\gamma}{2\beta} - E(X)\right)}{Var(X) =
#'   alpha/(2 beta) + E(X) (gamma/(2 beta) - E(X))}
#' where it is accurate. Where it would lose accuracy, under a strong tilt
#' \eqn{z = \gamma/\sqrt{\beta}}{z = gamma/sqrt(beta)} of either sign and at
#' a large \eqn{\alpha}{alpha} unless \eqn{\gamma = 0}, the variance is
#' instead obtained by numerical integration of the unnormalised density.
#'
#' For \eqn{\alpha \ge 1}{alpha >= 1}, the variance satisfies
#' \eqn{\mathrm{Var}(X) \le 1/(2\beta)}{Var(X) <= 1/(2 beta)} (Sun et al.,
#' 2023, Lemma 4c).
#'
#' The function is not vectorised in its parameters: \code{alpha},
#' \code{beta} and \code{gamma} must each be a single number. To evaluate
#' several parameter triples, use \code{\link{mapply}}, as in the examples.
#'
#' The numerical methods are described in Appendix A of
#' \code{vignette("theory", package = "mhn")}.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. (Lemmas 2c and 4c)
#'
#' @seealso \code{\link{mhn_mean}}, \code{\link{dmhn}}
#'
#' @examples
#' mhn_var(alpha = 2, beta = 1, gamma = 0)
#' mapply(mhn_var, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#'
#' @export
mhn_var <- function(alpha, beta, gamma) {
  .check_single(alpha, beta, gamma)
  .mhn_var_cpp(alpha, beta, gamma)
}


#' Skewness of the Modified Half-Normal Distribution
#'
#' Computes the skewness
#' \eqn{\gamma_1 = E[(X - \mu)^3] / \sigma^3}{gamma_1 = E[(X - mu)^3] / sigma^3}
#' for
#' \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}{X ~ MHN(alpha, beta,
#' gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}), a single number.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. A single number.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   A single number.
#'
#' @return A numeric scalar.
#'
#' @details
#' The skewness is formed from the second and third central moments. In
#' exact arithmetic these follow from the mean of Sun et al. (2023,
#' Lemma 2a) and the recurrence for the raw moments of their Lemma 2b; the
#' function obtains them instead, wherever it can, by numerical integration
#' of the unnormalised density, which needs no normalising constant and
#' stays accurate under a large tilt.
#'
#' Against independent quadratures, for shapes from 1 to 60 and tilts
#' \eqn{\gamma / \sqrt{\beta}}{gamma / sqrt(beta)} from -20 to 12, the
#' absolute error stayed below 1e-9; the largest found, 3e-10, is for a shape
#' between 2 and 3 under a positive tilt.
#'
#' For a subnormal shape, below the smallest normal double 2.2e-308, the
#' skewness depends on the tilt
#' \eqn{z = \gamma / \sqrt{\beta}}{z = gamma / sqrt(beta)} as follows.
#' While the mass at the peak of the density is small against the mass by
#' the origin (up to \eqn{z}{z} of about 53 at the largest subnormal shape
#' and 54.5 at the smallest), it grows like
#' \eqn{1 / \sqrt{\alpha}}{1 / sqrt(alpha)} as the shape falls, by a factor
#' that depends on the tilt: 1.25 at \eqn{z = 0} and
#' \eqn{2.8 \times 10^{-8}}{2.8e-8} at \eqn{z = 12}.  It turns negative
#' where the two masses are about equal (\eqn{z}{z} from about 53.3 to
#' 54.7), and above a tilt of 55 to 56.3 it no longer depends on
#' \eqn{\alpha}{alpha}.
#'
#' The function is not vectorised in its parameters: \code{alpha},
#' \code{beta} and \code{gamma} must each be a single number. To evaluate
#' several parameter triples, use \code{\link{mapply}}, as in the examples.
#'
#' The numerical methods are described in Appendix A of
#' \code{vignette("theory", package = "mhn")}.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. (Lemmas 2a and 2b)
#'
#' @seealso \code{\link{mhn_kurtosis}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_skewness(alpha = 2, beta = 1, gamma = 0)
#' mapply(mhn_skewness, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#'
#' @export
mhn_skewness <- function(alpha, beta, gamma) {
  .check_single(alpha, beta, gamma)
  .mhn_skewness_cpp(alpha, beta, gamma)
}


#' Excess Kurtosis of the Modified Half-Normal Distribution
#'
#' Computes the excess kurtosis
#' \eqn{\gamma_2 = E[(X - \mu)^4] / \sigma^4 - 3}{gamma_2 = E[(X - mu)^4] /
#' sigma^4 - 3}
#' for
#' \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}{X ~ MHN(alpha, beta,
#' gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}), a single number.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. A single number.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   A single number.
#'
#' @return A numeric scalar.
#'
#' @details
#' The excess kurtosis is formed from the second and fourth central
#' moments. In exact arithmetic these follow from the mean of Sun et al.
#' (2023, Lemma 2a) and the recurrence for the raw moments of their
#' Lemma 2b; the function obtains them instead, wherever it can, by
#' numerical integration of the unnormalised density, which needs no
#' normalising constant and stays accurate under a large tilt.
#'
#' Against independent quadratures, for shapes from 1 to 60 and tilts
#' \eqn{\gamma / \sqrt{\beta}}{gamma / sqrt(beta)} from -20 to 12, the
#' absolute error stayed below 1e-9; the largest found, 6e-10, is for a shape
#' between 2 and 3 under a positive tilt.
#'
#' For a subnormal shape, below the smallest normal double 2.2e-308, the
#' excess kurtosis depends on the tilt
#' \eqn{z = \gamma / \sqrt{\beta}}{z = gamma / sqrt(beta)} as follows.
#' While the mass at the peak of the density is small against the mass by
#' the origin (up to \eqn{z}{z} of about 53 at the largest subnormal shape
#' and 54.5 at the smallest), it grows like \eqn{1 / \alpha}{1 / alpha} as
#' the shape falls, by a factor that depends on the tilt: 2 at \eqn{z = 0} and
#' \eqn{8.2 \times 10^{-16}}{8.2e-16} at \eqn{z = 12}.  It is \code{Inf}
#' where it passes the largest double: below a shape of about 1.1e-308 at
#' \eqn{z = 0}, at every subnormal shape for \eqn{z = -30}, and at none for
#' \eqn{z = 12}.  Where the two masses are about equal (\eqn{z}{z} from
#' about 53.3 to 54.7) the peak takes over, and above a tilt of 55 to 56.3
#' the kurtosis no longer depends on \eqn{\alpha}{alpha}.
#'
#' The function is not vectorised in its parameters: \code{alpha},
#' \code{beta} and \code{gamma} must each be a single number. To evaluate
#' several parameter triples, use \code{\link{mapply}}, as in the examples.
#'
#' The numerical methods are described in Appendix A of
#' \code{vignette("theory", package = "mhn")}.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. (Lemmas 2a and 2b)
#'
#' @seealso \code{\link{mhn_skewness}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_kurtosis(alpha = 2, beta = 1, gamma = 0)
#' mapply(mhn_kurtosis, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#'
#' @export
mhn_kurtosis <- function(alpha, beta, gamma) {
  .check_single(alpha, beta, gamma)
  .mhn_kurtosis_cpp(alpha, beta, gamma)
}
