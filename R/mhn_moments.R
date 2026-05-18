# Moment functions for the Modified Half-Normal distribution.
# Thin wrappers around .mhn_mean_cpp / .mhn_var_cpp / .mhn_skewness_cpp /
# .mhn_kurtosis_cpp (src/mhn_moments.cpp).

#' Mean of the Modified Half-Normal Distribution
#'
#' Computes \eqn{E(X)} for \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}).
#' @param beta Scale parameter (\eqn{\beta > 0}).
#' @param gamma Location parameter (\eqn{\gamma \in R}).
#'
#' @return A numeric scalar.
#'
#' @details
#' The mean is computed as a ratio of Fox-Wright Psi functions:
#' \deqn{E(X) = \frac{\Psi[(\alpha+1)/2,\, \gamma/\sqrt{\beta}]}{
#'   \sqrt{\beta}\, \Psi[\alpha/2,\, \gamma/\sqrt{\beta}]}}
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536. (Lemma 2a)
#'
#' @seealso \code{\link{mhn_var}}, \code{\link{dmhn}}
#'
#' @examples
#' mhn_mean(alpha = 2, beta = 1, gamma = 0)
#'
#' @export
mhn_mean <- function(alpha, beta, gamma) {
  .mhn_mean_cpp(alpha, beta, gamma)
}


#' Variance of the Modified Half-Normal Distribution
#'
#' Computes \eqn{\mathrm{Var}(X)} for
#' \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}).
#' @param beta Scale parameter (\eqn{\beta > 0}).
#' @param gamma Location parameter (\eqn{\gamma \in R}).
#'
#' @return A numeric scalar.
#'
#' @details
#' Uses the formula (Sun et al., 2023, Lemma 2c):
#' \deqn{\mathrm{Var}(X) = \frac{\alpha}{2\beta} +
#'   E(X)\left(\frac{\gamma}{2\beta} - E(X)\right)}
#'
#' For \eqn{\alpha \geq 1}, the variance satisfies
#' \eqn{\mathrm{Var}(X) \leq 1/(2\beta)} (Sun et al., 2023, Lemma 4c).
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536. (Lemma 2c)
#'
#' @seealso \code{\link{mhn_mean}}, \code{\link{dmhn}}
#'
#' @examples
#' mhn_var(alpha = 2, beta = 1, gamma = 0)
#'
#' @export
mhn_var <- function(alpha, beta, gamma) {
  .mhn_var_cpp(alpha, beta, gamma)
}


#' Skewness of the Modified Half-Normal Distribution
#'
#' Computes the skewness \eqn{\gamma_1 = E[(X - \mu)^3] / \sigma^3}
#' for \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}).
#' @param beta Scale parameter (\eqn{\beta > 0}).
#' @param gamma Location parameter (\eqn{\gamma \in R}).
#'
#' @return A numeric scalar.
#'
#' @details
#' Uses the moment recurrence (Sun et al., 2023, Lemma 2b) to compute
#' raw moments up to third order, then converts to central moments.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536. (Lemma 2b)
#'
#' @seealso \code{\link{mhn_kurtosis}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_skewness(alpha = 2, beta = 1, gamma = 0)
#'
#' @export
mhn_skewness <- function(alpha, beta, gamma) {
  .mhn_skewness_cpp(alpha, beta, gamma)
}


#' Excess Kurtosis of the Modified Half-Normal Distribution
#'
#' Computes the excess kurtosis
#' \eqn{\gamma_2 = E[(X - \mu)^4] / \sigma^4 - 3}
#' for \eqn{X \sim \mathrm{MHN}(\alpha, \beta, \gamma)}.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}).
#' @param beta Scale parameter (\eqn{\beta > 0}).
#' @param gamma Location parameter (\eqn{\gamma \in R}).
#'
#' @return A numeric scalar.
#'
#' @details
#' Uses the moment recurrence (Sun et al., 2023, Lemma 2b) to compute
#' raw moments up to fourth order, then converts to central moments.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536. (Lemma 2b)
#'
#' @seealso \code{\link{mhn_skewness}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_kurtosis(alpha = 2, beta = 1, gamma = 0)
#'
#' @export
mhn_kurtosis <- function(alpha, beta, gamma) {
  .mhn_kurtosis_cpp(alpha, beta, gamma)
}
