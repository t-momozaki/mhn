# Quantile function for the Modified Half-Normal distribution.
# Thin wrapper around .qmhn_cpp (src/mhn_qmhn.cpp).

#' Quantile Function of the Modified Half-Normal Distribution
#'
#' Computes the quantile (inverse cumulative) function of the
#' Modified Half-Normal (MHN) distribution with parameters \code{alpha},
#' \code{beta}, and \code{gamma}.
#'
#' For the general case, \eqn{q = F^{-1}(p)}{q = F^-1(p)} is obtained by a
#' root-finder applied to the distribution function computed as in
#' \code{\link{pmhn}}, or for \code{lower.tail = FALSE} to the logarithm of
#' its directly computed upper tail.
#'
#' @param p Numeric vector of probabilities.
#' @param alpha Shape parameter (\eqn{\alpha > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of
#'   the Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values
#'   pull mass toward zero. Scalar or numeric vector. Default: 1.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   Scalar or numeric vector. Default: 0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}, \eqn{x} being the
#'   quantile returned.
#' @param log.p Logical; if \code{TRUE}, probabilities are provided on the
#'   log scale. Default: \code{FALSE}. Both flags are read with
#'   \code{as.logical()}, so that 1 is \code{TRUE}; a value that gives
#'   \code{NA} is an error.
#'
#' @return A numeric vector. The output length equals
#'   \code{max(length(p), length(alpha), length(beta), length(gamma))}, or 0
#'   if any of them has length zero; each input is recycled to that length
#'   following standard R recycling rules.
#'   As in \code{\link[stats]{qgamma}}, an element whose parameters lie
#'   outside the parameter space (\code{alpha} or \code{beta} not positive
#'   and finite, or \code{gamma} not finite) is \code{NaN}, with a warning,
#'   and an element with an \code{NA} or \code{NaN} among its inputs is
#'   \code{NA} or \code{NaN} (\code{NA} if any of them is \code{NA}),
#'   without one; the other elements are unaffected.  The arguments must be
#'   numeric (double, integer or logical): a factor, a character string, a
#'   complex number, a list or \code{NULL} is an error, as in
#'   \code{\link[stats]{qnorm}}.
#'   \code{qmhn(0) = 0} and \code{qmhn(1) = Inf}. Probabilities outside
#'   \eqn{[0, 1]} (log probabilities above 0 when \code{log.p = TRUE})
#'   yield \code{NaN}, with a warning as in \code{\link[stats]{qnorm}}.  If
#'   the distribution function cannot be evaluated where the root is sought,
#'   a warning is issued and \code{NA} is returned for that element.  The
#'   attributes of \code{p} (names, \code{dim}, \code{dimnames}) are kept
#'   when the result has the length of \code{p}; a classed \code{p}, such
#'   as a Date, is taken as its numbers and gives a plain vector.
#'   Unlike \code{qnorm}, attributes of a longer parameter vector are not
#'   carried over.
#'
#' @details
#' The default parameters \code{alpha = 1, beta = 1, gamma = 0} correspond
#' to the half-normal distribution
#' \eqn{\mathrm{HN}(1/\sqrt{2})}{HN(1/sqrt(2))}.
#'
#' At \eqn{\gamma = 0}{gamma = 0} and at \eqn{\alpha = 1}{alpha = 1}
#' exactly the quantile is computed from the law it reduces to; a shape or a
#' tilt merely close to those values is a different law, and is computed as
#' in the general case:
#' \itemize{
#'   \item \eqn{\gamma = 0}: \eqn{q = \sqrt{y / \beta}}{q = sqrt(y / beta)},
#'     where \eqn{P(\alpha/2, y) = p} for the regularized incomplete gamma
#'     function \eqn{P} (Sun et al. 2023, Lemma 6a), for
#'     \eqn{\alpha \le 10^6}{alpha <= 1e6}; a larger shape is computed as in
#'     the general case.
#'   \item \eqn{\alpha = 1}: the truncated-normal inverse (Sun et al. 2023,
#'     Lemma 6b).  The quantile of a small lower-tail probability keeps its
#'     relative accuracy.  Where \eqn{-\gamma/\sqrt{2\beta}}{-gamma /
#'     sqrt(2 beta)} overflows the law is Exponential with rate
#'     \eqn{|\gamma|}{|gamma|}, and its quantile is that one's.
#' }
#'
#' A probability close to 1 is solved, where that is more accurate, on the
#' other tail for one minus it, which is formed exactly: a lower-tail
#' \code{p} on the upper tail, as the root of
#' \eqn{P(X > x) = 1 - p}{P(X > x) = 1 - p}, since the distribution
#' function is resolved there only to the rounding of 1 (and, for a small
#' shape, to the error of the lower tail; see \code{\link{pmhn}}), and an
#' upper-tail \code{p} on the distribution function.  Where \code{pmhn}
#' takes the distribution function from its series
#' (\eqn{\alpha < 1}{alpha < 1} with \eqn{|z|}{|z|} up to about 830, and
#' \eqn{\alpha \le 100}{alpha <= 100} with \eqn{|z| \le 20}{|z| <= 20},
#' \eqn{z = \gamma/\sqrt{\beta}}{z = gamma/sqrt(beta)}), a lower-tail
#' \code{p} is solved on the series itself, and for
#' \eqn{\alpha < 1}{alpha < 1} with \eqn{|z|}{|z|} from about 100 one
#' quantile can take from a tenth of a second to seconds (about 1 s at
#' \eqn{\alpha = 0.8}{alpha = 0.8}, \eqn{z = 829}{z = 829} and
#' \eqn{p = 0.5}{p = 0.5}); where \eqn{1 - p}{1 - p} is known exactly,
#' \code{lower.tail = FALSE} gives the same quantile at once.
#'
#' \strong{Accuracy and range.}  The root is found to a relative
#' \eqn{2^{-46}}{2^-46}, and where \code{pmhn} measures the law from its
#' peak (see there) also to within \eqn{2^{-46}}{2^-46} of the law's
#' standard deviation, which can be narrower.  The quantile inverts
#' \code{\link{pmhn}}, whose range of accuracy it inherits: for
#' \eqn{\alpha \le 10^{15}}{alpha <= 1e15},
#' \eqn{|z| \le 10^{15}}{|z| <= 1e15} and \eqn{\beta}{beta} in
#' \eqn{[10^{-300}, 10^{300}]}{[1e-300, 1e300]} it is the root, to that
#' tolerance or to the spacing of the doubles near it, whichever is larger,
#' of an equation whose tail probabilities are as accurate as those of
#' \code{pmhn} (to the working tolerance within 10 standard deviations of
#' the peak, and to about twice that beyond; see there), which places it
#' within about that many of the law's standard deviations.
#' Next to the origin, where \eqn{F(x)}{F(x)} is about a multiple of
#' \eqn{x^\alpha}{x^alpha}, its relative error is that of the smaller tail
#' probability divided by about \eqn{\alpha}{alpha}.  A quantile below
#' \code{.Machine$double.xmin} is resolved to the spacing of the subnormal
#' doubles there, about \eqn{4.9 \times 10^{-324}}{4.9e-324}, and one below
#' that spacing is returned as 0, as \code{\link[stats]{qgamma}} returns a
#' quantile that underflows; a quantile beyond the largest double is
#' returned as \code{Inf}, as \code{qgamma} returns one that overflows.
#' Where the law is narrower than the spacing of the doubles (for
#' \eqn{z}{z} above about \eqn{10^{16}}{1e16}, or \eqn{\alpha}{alpha} above
#' about \eqn{10^{31}}{1e31}) its distribution function is a step between
#' two adjacent doubles, and the quantile of every \code{p} strictly inside
#' \eqn{(0, 1)} is the one of the two nearer to the true quantile, placed
#' from the exact offset of the law's peak.  Where \code{pmhn} is \code{NA}
#' (see there), so is the quantile, with a warning.  Where the log upper
#' tail of \code{pmhn} is \code{-Inf} in place of a finite value
#' (\eqn{\beta}{beta} and the tilt both close to the largest double, see
#' there), the quantile of a \code{log.p} probability beyond that point is
#' the point at which it turns \code{-Inf}, without a warning.  For a
#' subnormal shape (see \code{\link{dmhn}}) and a tilt \eqn{z} below about
#' 52 (53 at the smallest shapes), \code{qmhn} is 0 for every lower-tail
#' probability below 1 given on the probability scale; the upper tail, or
#' \code{log.p = TRUE}, still gives the quantiles of the far upper tail.
#'
#' When any of \code{alpha}, \code{beta}, \code{gamma} is a vector, the
#' quantile is evaluated element-wise.  The Fox-Wright \eqn{\Psi}
#' normalizing constant and moments \eqn{E(X)},
#' \eqn{\mathrm{Var}(X)}{Var(X)} (used to size the root-finder bracket) are
#' computed once for each position of the recycling cycle while that cycle
#' is at most 256 long, and otherwise once for each run of consecutive
#' elements that share the \eqn{(\alpha, \beta, \gamma)} triple.
#'
#' The numerical methods, and the reasons for the ranges above, are
#' described in Appendix A of \code{vignette("theory", package = "mhn")}.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613.
#'
#' @seealso \code{\link{dmhn}}, \code{\link{pmhn}}, \code{\link{rmhn}}
#'
#' @examples
#' # Basic evaluation
#' qmhn(c(0.1, 0.5, 0.9), alpha = 2, beta = 1, gamma = 1)
#'
#' # Round-trip: F(F^-1(p)) ~ p
#' p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' all.equal(pmhn(qmhn(p, alpha = 2, beta = 1, gamma = 1),
#'                alpha = 2, beta = 1, gamma = 1),
#'           p, tolerance = 1e-6)
#'
#' # Tail / log forms
#' qmhn(0.95, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
#' qmhn(log(0.05), alpha = 2, beta = 1, gamma = 1, log.p = TRUE)
#'
#' @export
qmhn <- function(p, alpha = 1, beta = 1, gamma = 0,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!(is.double(p) && is.double(alpha) && is.double(beta) && is.double(gamma))) {
    .check_numeric(alpha, beta, gamma, p, "p")
  }
  out <- .qmhn_cpp(as.numeric(p),
                   as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
                   .as_flag(lower.tail, "lower.tail"), .as_flag(log.p, "log.p"))
  if (anyNA(out)) .warn_nan(out, p, alpha, beta, gamma)
  # Keep names and dimensions the way qnorm does (see dmhn).
  if (!is.null(attributes(p)) && !is.object(p) && length(out) == length(p)) {
    attributes(out) <- attributes(p)
  }
  out
}
