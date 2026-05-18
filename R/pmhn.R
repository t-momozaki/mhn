# Cumulative distribution function for the Modified Half-Normal distribution.
# Thin wrapper around .pmhn_cpp (src/mhn_pmhn.cpp).

#' Distribution Function of the Modified Half-Normal Distribution
#'
#' Computes the cumulative distribution function (CDF) of the
#' Modified Half-Normal (MHN) distribution with parameters \code{alpha},
#' \code{beta}, and \code{gamma}.
#'
#' The CDF is computed via the series representation
#' \deqn{F(x \mid \alpha, \beta, \gamma) =
#'   \frac{1}{\Psi[\alpha/2, \gamma/\sqrt{\beta}]}
#'   \sum_{i=0}^{\infty} \frac{z^i}{i!}\, \Gamma(s_i)\, P(s_i, \beta x^2)}
#' where \eqn{z = \gamma/\sqrt{\beta}}, \eqn{s_i = (\alpha + i)/2}, and
#' \eqn{P(s, y)} is the regularized lower incomplete gamma function
#' (Sun et al., 2023, Lemma 1b; equivalent to the paper's form via the
#' identity \eqn{\Gamma(s)\, P(s, y) = \gamma(s, y)}, where
#' \eqn{\gamma(s, y)} is the lower incomplete gamma function used in
#' the paper).  The infinite sum is truncated at the constructive bound
#' \eqn{K = \max\{K_1, K_2\}} from Sun et al. (2023), Supplementary
#' Lemma 10(d), which makes the truncation residual bounded by the
#' user's tolerance divided by \eqn{\Psi}.  When
#' double-precision cancellation in the alternating-sign accumulator
#' for \eqn{\gamma < 0} would exceed that tolerance, the series is
#' replaced by a Gauss-Kronrod (or tanh-sinh for \eqn{\alpha < 1})
#' numerical integration of the density on \eqn{[0, q]}.
#'
#' @param q Numeric vector of quantiles.
#' @param alpha Shape parameter (\eqn{\alpha > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param beta Scale parameter (\eqn{\beta > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param gamma Location parameter (\eqn{\gamma \in R}). Scalar or numeric
#'   vector. Default: 0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the
#'   log scale. Default: \code{FALSE}.
#'
#' @return A numeric vector. The output length equals
#'   \code{max(length(q), length(alpha), length(beta), length(gamma))}; each
#'   input is recycled to that length following standard R recycling rules.
#'   For \code{q <= 0} the CDF is 0; for \code{q = Inf} it is 1.
#'
#' @details
#' Special cases are detected and dispatched to standard R primitives:
#' \itemize{
#'   \item \eqn{\gamma = 0}: \code{pgamma(q^2, alpha/2, scale = 1/beta)}
#'   \item \eqn{\alpha = 1}: truncated-normal CDF via \code{pnorm}
#' }
#'
#' When any of \code{alpha}, \code{beta}, \code{gamma} is a vector, the CDF
#' is evaluated element-wise.  The Fox-Wright \eqn{\Psi} normalizing
#' constant is recomputed only when consecutive elements present a
#' different \eqn{(\alpha, \beta, \gamma)} triple, so passing grouped
#' parameters is significantly faster than calling \code{pmhn} inside an
#' R loop.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1507--1536.
#'
#' @seealso \code{\link{dmhn}}, \code{\link{qmhn}}, \code{\link{rmhn}}
#'
#' @examples
#' # Basic evaluation
#' pmhn(c(0.5, 1, 1.5), alpha = 2, beta = 1, gamma = 1)
#'
#' # Tail / log forms
#' pmhn(2, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#' pmhn(2, alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#'
#' # Special case: gamma = 0 reduces to sqrt-Gamma
#' all.equal(pmhn(1.5, alpha = 2, beta = 1, gamma = 0),
#'           pgamma(1.5^2, shape = 1, rate = 1))
#'
#' @export
pmhn <- function(q, alpha = 1, beta = 1, gamma = 0,
                 lower.tail = TRUE, log.p = FALSE) {
  .pmhn_cpp(as.numeric(q),
            as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
            isTRUE(lower.tail), isTRUE(log.p))
}

# Diagnostic hook for inst/audits/cdf_series_accuracy.R: evaluates
# the general-case CDF via a single, caller-specified code path,
# bypassing both the special-case shortcuts and the automatic
# series-to-integration fallback used by pmhn().  Scalar arguments
# only.  Internal helper: leading dot follows the package convention
# (cf. .check_mhn_params, .dmhn_cpp), and @noRd suppresses Rd output.
#' @noRd
.pmhn_force <- function(q, alpha, beta, gamma,
                        method = c("series", "integrate")) {
  method <- match.arg(method)
  .pmhn_force_cpp(as.numeric(q),
                  as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
                  method)
}
