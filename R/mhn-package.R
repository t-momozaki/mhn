#' mhn: The Modified Half-Normal Distribution
#'
#' Provides density (\code{\link{dmhn}}), distribution function
#' (\code{\link{pmhn}}), quantile function (\code{\link{qmhn}}),
#' random generation (\code{\link{rmhn}}), and moments / mode of the
#' Modified Half-Normal (MHN) distribution.
#'
#' @section MHN Distribution:
#' The MHN distribution has probability density function
#' \deqn{f(x \mid \alpha, \beta, \gamma) =
#'   \frac{2\beta^{\alpha/2} x^{\alpha-1}
#'   \exp(-\beta x^2 + \gamma x)}{\Psi[\alpha/2, \gamma/\sqrt{\beta}]}}
#' for \eqn{x > 0}, where \eqn{\Psi[a, z]} is the Fox-Wright Psi function.
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{alpha}}{Shape parameter (\eqn{\alpha > 0}).
#'     Controls the \eqn{x^{\alpha-1}} term.}
#'   \item{\code{beta}}{Scale (rate) parameter (\eqn{\beta > 0}).
#'     Controls the \eqn{\exp(-\beta x^2)} term.}
#'   \item{\code{gamma}}{Location (skewness) parameter (\eqn{\gamma \in R}).
#'     Controls the \eqn{\exp(\gamma x)} term.}
#' }
#'
#' @section Special Cases:
#' \describe{
#'   \item{\eqn{\gamma = 0}}{Square-root Gamma distribution
#'     (\eqn{X^2 \sim \textrm{Gamma}(\alpha/2, \beta)})
#'     (Sun et al., 2023, Lemma 6a).}
#'   \item{\eqn{\alpha = 1}}{Truncated normal distribution
#'     \eqn{\textrm{TN}(\gamma/(2\beta), 1/\sqrt{2\beta}, 0, \infty)}
#'     (Sun et al., 2023, Lemma 6b).}
#'   \item{\eqn{\alpha = 1, \gamma = 0}}{Half-normal distribution
#'     \eqn{\textrm{HN}(1/\sqrt{2\beta})}
#'     (Sun et al., 2023, Lemma 6c).}
#' }
#'
#' @section References:
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal distribution:
#' Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5), 1507-1536.
#' \doi{10.1080/03610926.2021.1934700}
#'
#' Gao, F. & Wang, H.-B. (2025). Generating modified-half-normal random
#' variates by a relaxed transformed density rejection method.
#' \emph{Communications in Statistics - Simulation and Computation}.
#' \doi{10.1080/03610918.2025.2524551}
#'
#' @importFrom stats dgamma pgamma qgamma rgamma
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @importFrom stats integrate uniroot
#' @importFrom Rcpp sourceCpp
#' @useDynLib mhn, .registration = TRUE
#' @keywords internal
"_PACKAGE"
