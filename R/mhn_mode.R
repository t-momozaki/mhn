# Mode of the Modified Half-Normal distribution.
# Thin wrapper around .mhn_mode_cpp (src/mhn_mode.cpp).

#' Mode of the Modified Half-Normal Distribution
#'
#' Computes the mode of the MHN distribution. For
#' \eqn{\alpha \ge 1}{alpha >= 1} this is the maximiser of the density. For
#' \eqn{0 < \alpha < 1} the density is unbounded as
#' \eqn{x \to 0^+}{x -> 0+}, so a non-\code{NA} result is an interior local
#' maximum, not the most probable value.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}), a single number.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. A single number.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   A single number.
#'
#' @return A numeric scalar. For \eqn{0 < \alpha < 1} with no interior
#'   local maximum (density strictly decreasing on \eqn{(0, \infty)}) the
#'   result is \code{NA}; at \eqn{\alpha = 1} with
#'   \eqn{\gamma \le 0}{gamma <= 0} the mode sits on the boundary and 0 is
#'   returned. See Details.
#'
#' @details
#' The mode depends on \eqn{\alpha}:
#' \describe{
#'   \item{\eqn{\alpha > 1}}{
#'     \eqn{(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) /
#'     (4\beta)}{(gamma + sqrt(gamma^2 + 8 beta (alpha - 1))) / (4 beta)}
#'     (Sun et al. 2023, Lemma 3b).}
#'   \item{\eqn{\alpha = 1}}{
#'     \eqn{\max(0, \gamma / (2\beta))}{max(0, gamma / (2 beta))}, obtained
#'     as the mode of the truncated normal
#'     \eqn{\mathrm{TN}(\gamma/(2\beta), 1/\sqrt{2\beta}, 0,
#'     \infty)}{TN(gamma / (2 beta), 1 / sqrt(2 beta), 0, Inf)} that the MHN
#'     reduces to in this case (Sun et al. 2023, Lemma 6b). For
#'     \eqn{\gamma \le 0}{gamma <= 0} this is the boundary point 0. Sun et
#'     al. (2023, Lemma 3d) state that no mode exists for \eqn{\gamma < 0}
#'     and \eqn{\alpha \le 1}{alpha <= 1}; at \eqn{\alpha = 1} the density
#'     is bounded and decreases from its supremum at \eqn{x = 0}, and that
#'     point is returned. Only
#'     \eqn{\alpha = 1} itself takes this branch. The Lemma 3b value tends
#'     to it as \eqn{\alpha \to 1}{alpha -> 1} from above, and for
#'     \eqn{\gamma > 0} so does the Lemma 3c value from below; for
#'     \eqn{\gamma \le 0}{gamma <= 0} and \eqn{\alpha < 1} there is no mode
#'     (\code{NA}).}
#'   \item{\eqn{0 < \alpha < 1}}{
#'     An interior local maximum exists only when \eqn{\gamma > 0} and
#'     \eqn{\alpha > 1 - \gamma^2 / (8\beta)}{alpha > 1 - gamma^2 / (8 beta)},
#'     and it is then located at
#'     \eqn{(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) /
#'     (4\beta)}{(gamma + sqrt(gamma^2 + 8 beta (alpha - 1))) / (4 beta)}
#'     (Sun et al. 2023, Lemma 3c). Lemma 3c states the condition with
#'     \eqn{\ge}{>=}; the strict inequality is used here because at
#'     \eqn{\alpha = 1 - \gamma^2/(8\beta)}{alpha = 1 - gamma^2 / (8 beta)}
#'     the two stationary points coincide at
#'     \eqn{\gamma/(4\beta)}{gamma / (4 beta)}, where
#'     \eqn{\frac{d}{dx} \log f(x) = -\frac{2\beta}{x}
#'     \left(x - \frac{\gamma}{4\beta}\right)^2 \le 0}{d/dx log f(x) =
#'     -(2 beta / x) (x - gamma / (4 beta))^2 <= 0}, so the density is still
#'     strictly decreasing and has an inflection point, not a maximum,
#'     there. In every other case the density is strictly decreasing on
#'     \eqn{(0, \infty)} (Sun et al. 2023, Lemma 3d, for
#'     \eqn{\gamma \ne 0}{gamma != 0}; for \eqn{\gamma = 0} the derivative
#'     \eqn{-(2\beta x^2 + 1 - \alpha)/x}{-(2 beta x^2 + 1 - alpha) / x} of
#'     \eqn{\log f}{log f} is negative) and \code{NA} is returned.}
#' }
#'
#' The function is not vectorised in its parameters: \code{alpha},
#' \code{beta} and \code{gamma} must each be a single number. To evaluate
#' several parameter triples, use \code{\link{mapply}}, as in the examples.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. (Lemma 3b--d, Lemma 6b)
#'
#' @seealso \code{\link{dmhn}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_mode(alpha = 2, beta = 1, gamma = 1)
#' mhn_mode(alpha = 1, beta = 1, gamma = 2)
#' mhn_mode(alpha = 0.5, beta = 1, gamma = -1)  # NA
#' mhn_mode(alpha = 0.5, beta = 1, gamma = 5)   # a local maximum only
#' mapply(mhn_mode, alpha = c(0.5, 2, 5), beta = 1, gamma = c(-1, 0, 3))
#'
#' @export
mhn_mode <- function(alpha, beta, gamma) {
  .check_single(alpha, beta, gamma)
  .mhn_mode_cpp(alpha, beta, gamma)
}
