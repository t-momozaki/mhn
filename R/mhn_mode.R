# Mode of the Modified Half-Normal distribution.
# Thin wrapper around .mhn_mode_cpp (src/mhn_mode.cpp).

#' Mode of the Modified Half-Normal Distribution
#'
#' Computes the mode (most probable value) of the MHN distribution.
#'
#' @param alpha Shape parameter (\eqn{\alpha > 0}).
#' @param beta Scale parameter (\eqn{\beta > 0}).
#' @param gamma Location parameter (\eqn{\gamma \in R}).
#'
#' @return A numeric scalar. Returns \code{NA} when no interior mode exists
#'   (density is monotonically decreasing on \eqn{(0, \infty)}).
#'
#' @details
#' The mode depends on \eqn{\alpha}:
#' \describe{
#'   \item{\eqn{\alpha > 1}}{
#'     \eqn{(\gamma + \sqrt{\gamma^2 + 8\beta(\alpha - 1)}) / (4\beta)}
#'     (Sun et al., 2023, Lemma 3b).}
#'   \item{\eqn{\alpha = 1}}{
#'     \eqn{\max(0, \gamma / (2\beta))}, obtained as the mode of the
#'     truncated normal \eqn{\mathrm{TN}(\gamma/(2\beta), 1/\sqrt{2\beta}, 0,
#'     \infty)} that the MHN reduces to in this case (Sun et al., 2023,
#'     Lemma 6b).}
#'   \item{\eqn{0 < \alpha < 1}}{
#'     An interior mode exists only when \eqn{\gamma > 0} and
#'     \eqn{\alpha \geq 1 - \gamma^2 / (8\beta)} (Sun et al., 2023,
#'     Lemma 3c); otherwise the density is monotonically decreasing
#'     (Sun et al., 2023, Lemma 3d) and \code{NA} is returned.}
#' }
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536. (Lemma 3b--d, Lemma 6b)
#'
#' @seealso \code{\link{dmhn}}, \code{\link{mhn_mean}}
#'
#' @examples
#' mhn_mode(alpha = 2, beta = 1, gamma = 1)
#' mhn_mode(alpha = 1, beta = 1, gamma = 2)
#' mhn_mode(alpha = 0.5, beta = 1, gamma = -1)  # NA
#'
#' @export
mhn_mode <- function(alpha, beta, gamma) {
  .mhn_mode_cpp(alpha, beta, gamma)
}
