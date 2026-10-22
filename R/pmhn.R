# Cumulative distribution function for the Modified Half-Normal distribution.
# Thin wrapper around .pmhn_cpp (src/mhn_pmhn.cpp).

#' Distribution Function of the Modified Half-Normal Distribution
#'
#' Computes the cumulative distribution function (CDF) of the
#' Modified Half-Normal (MHN) distribution with parameters \code{alpha},
#' \code{beta}, and \code{gamma}.
#'
#' In the general case the lower tail is computed, where that is accurate
#' (see below), from the series representation (Sun et al. 2023,
#' Lemma 1b)
#' \deqn{F(x \mid \alpha, \beta, \gamma) =
#'   \frac{1}{\Psi[\alpha/2, \gamma/\sqrt{\beta}]}
#'   \sum_{i=0}^{\infty} \frac{z^i}{i!}\, \Gamma(s_i)\, P(s_i, \beta x^2)}{
#'   F(x | alpha, beta, gamma) = (1 / Psi[alpha/2, gamma/sqrt(beta)])
#'   sum_(i >= 0) z^i / i! Gamma(s_i) P(s_i, beta x^2)}
#' where \eqn{z = \gamma/\sqrt{\beta}}{z = gamma/sqrt(beta)},
#' \eqn{s_i = (\alpha + i)/2}, and \eqn{P(s, y)} is the regularized lower
#' incomplete gamma function.
#'
#' @param q Numeric vector of quantiles.
#' @param alpha Shape parameter (\eqn{\alpha > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of
#'   the Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values
#'   pull mass toward zero. Scalar or numeric vector. Default: 1.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   Scalar or numeric vector. Default: 0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the
#'   log scale. Default: \code{FALSE}. Both flags are read with
#'   \code{as.logical()}, so that 1 is \code{TRUE}; a value that gives
#'   \code{NA} is an error.
#'
#' @return A numeric vector. The output length equals
#'   \code{max(length(q), length(alpha), length(beta), length(gamma))}, or 0
#'   if any of them has length zero; each input is recycled to that length
#'   following standard R recycling rules.
#'   For \code{q <= 0} the CDF is 0; for \code{q = Inf} it is 1.
#'   As in \code{\link[stats]{pgamma}}, an element whose parameters lie
#'   outside the parameter space (\code{alpha} or \code{beta} not positive
#'   and finite, or \code{gamma} not finite) is \code{NaN}, with a warning,
#'   and an element with an \code{NA} or \code{NaN} among its inputs is
#'   \code{NA} or \code{NaN} (\code{NA} if any of them is \code{NA}),
#'   without one; the other elements are unaffected.  The arguments must be
#'   numeric (double, integer or logical): a factor, a character string, a
#'   complex number, a list or \code{NULL} is an error, as in
#'   \code{\link[stats]{pnorm}}.
#'   As with \code{\link[stats]{pnorm}}, the attributes of \code{q} (names,
#'   \code{dim}, \code{dimnames}) are kept when the result has the length of
#'   \code{q}; a classed \code{q}, such as a Date, is taken as its numbers
#'   and gives a plain vector.
#'   Unlike \code{pnorm}, attributes of a longer parameter vector are not
#'   carried over.
#'
#' @details
#' The default parameters \code{alpha = 1, beta = 1, gamma = 0} correspond
#' to the half-normal distribution
#' \eqn{\mathrm{HN}(1/\sqrt{2})}{HN(1/sqrt(2))}.
#'
#' At \eqn{\gamma = 0}{gamma = 0} and at \eqn{\alpha = 1}{alpha = 1}
#' exactly the distribution function is computed from the law it reduces to;
#' a shape or a tilt merely close to those values is a different law, and is
#' computed as in the general case:
#' \itemize{
#'   \item \eqn{\gamma = 0}: \eqn{F(q) = P(\alpha/2, \beta q^2)}{F(q) =
#'     P(alpha/2, beta q^2)}, the regularized incomplete gamma function of
#'     the square-root-Gamma law (Sun et al. 2023, Lemma 6a), for
#'     \eqn{\alpha \le 10^6}{alpha <= 1e6}; a larger shape is computed as in
#'     the general case.
#'   \item \eqn{\alpha = 1}: the truncated-normal CDF (Sun et al. 2023,
#'     Lemma 6b).  Near the truncation point the lower tail keeps its
#'     relative accuracy however small it is.  Where
#'     \eqn{-\gamma/\sqrt{2\beta}}{-gamma / sqrt(2 beta)} overflows the law
#'     is Exponential with rate \eqn{|\gamma|}{|gamma|}.
#' }
#'
#' The series is used, where it is accurate, for
#' \eqn{\alpha \le 100}{alpha <= 100} with \eqn{|z| \le 20}{|z| <= 20}, and
#' for \eqn{\alpha < 1}{alpha < 1} up to \eqn{|z|}{|z|} of about 830 (for
#' \eqn{|z| > 20}{|z| > 20} it costs milliseconds per value); elsewhere
#' the lower tail is computed by numerical integration of the density.  The
#' upper tail \eqn{P(X > q)} is computed directly, by numerical integration
#' of the density over \eqn{(q, \infty)}{(q, Inf)}, and not as
#' \eqn{1 - F}, which would round to 0 once it falls below about
#' \eqn{10^{-16}}{1e-16}.  Each tail is computed directly on its own side of
#' the peak of the density of \eqn{\log X}{log X}, which lies at or above
#' the median, and as one minus the other tail beyond it.  For a large shape
#' or tilt, where the terms of these forms grow too large, both tails are
#' measured from that peak instead.  With \code{log.p = TRUE} both tails are
#' formed on the log scale throughout, so that \eqn{\log F}{log F} is still
#' resolved where \eqn{F} rounds to 1 or is below the smallest positive
#' double, and a log probability close to 0, about minus the other tail,
#' keeps the relative accuracy of that other tail to within about
#' \eqn{10^{-9}}{1e-9} (Appendix A.4 of the vignette cited below).
#'
#' \strong{Accuracy and range.}  With
#' \eqn{z = \gamma/\sqrt{\beta}}{z = gamma/sqrt(beta)} as above, for
#' \eqn{\alpha \le 10^{15}}{alpha <= 1e15},
#' \eqn{|z| \le 10^{15}}{|z| <= 1e15} and \eqn{\beta}{beta} in
#' \eqn{[10^{-300}, 10^{300}]}{[1e-300, 1e300]} both tails are accurate to
#' the working tolerance \code{sqrt(.Machine$double.eps)} within 10 standard
#' deviations of the peak, and to about twice that beyond.  Outside that
#' range the law can be narrower than the spacing of the doubles near its
#' peak (for \eqn{z}{z} above about \eqn{10^{16}}{1e16}, or
#' \eqn{\alpha}{alpha} above about \eqn{10^{31}}{1e31}); the probabilities
#' are then still those of the exact law at the doubles \code{q} given, a
#' step from 0 to 1 between two adjacent doubles, whose value at the double
#' nearest the peak is set by where the peak lies within its spacing.  Two
#' cases give \code{NA} with a warning instead: within 45 standard
#' deviations of the peak of a law whose peak cannot be placed to within
#' \eqn{10^{-11}}{1e-11} of its standard deviation, which needs a law
#' narrower than the spacing by a factor of order \eqn{10^{19}}{1e19} or
#' more and then happens for some parameters only (as at
#' \eqn{\alpha = \beta = \gamma = 10^{200}}{alpha = beta = gamma = 1e200});
#' and where the whole law lies beyond the largest
#' double without a positive tilt to place it (a shape near the largest
#' double against a subnormal \eqn{\beta}{beta}).  A third case lies at the
#' other end of the range: where \eqn{\beta}{beta} and the tilt are both
#' close to the largest double (\eqn{\beta}{beta} from about
#' \eqn{2 \times 10^{307}}{2e307} with \eqn{\gamma}{gamma} near the
#' largest double, or from about \eqn{5 \times 10^{307}}{5e307} with
#' \eqn{\gamma}{gamma} from about \eqn{10^{306}}{1e306}), at some points
#' \eqn{q}{q} above the peak where \eqn{\beta q}{beta q} is above about half
#' the largest double, the log upper tail is \code{-Inf} where it is finite
#' (the upper tail itself is then 0, which is its value in double
#' precision), and next to the point where \eqn{\beta q}{beta q} overflows
#' both are \code{NA} with a warning.
#' Where \eqn{z}{z} itself overflows, the law is its limit:
#' Gamma(\eqn{\alpha}{alpha}, rate \eqn{|\gamma|}{|gamma|}) under a negative
#' tilt, and normal about its mode with variance
#' \eqn{1/(2\beta)}{1/(2 beta)} under a positive one.  A subnormal shape,
#' below the smallest normal double, is computed apart (see
#' \code{\link{dmhn}}), and there \eqn{\log F}{log F} keeps its relative
#' accuracy too; below a tilt \eqn{z} of about 52 (53 at the smallest
#' shapes) \code{pmhn} is then 1 at every positive \code{q}.
#'
#' When any of \code{alpha}, \code{beta}, \code{gamma} is a vector, the CDF
#' is evaluated element-wise.  The Fox-Wright \eqn{\Psi} normalizing
#' constant is computed once for each position of the recycling cycle while
#' that cycle is at most 256 long, and otherwise once for each run of
#' consecutive elements that share the \eqn{(\alpha, \beta, \gamma)}
#' triple.  So recycled parameters, and parameter vectors in which equal
#' triples are adjacent, are much faster than triples that alternate along a
#' long vector, and than calling \code{pmhn} inside an R loop.
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
  if (!(is.double(q) && is.double(alpha) && is.double(beta) && is.double(gamma))) {
    .check_numeric(alpha, beta, gamma, q, "q")
  }
  out <- .pmhn_cpp(as.numeric(q),
                   as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
                   .as_flag(lower.tail, "lower.tail"), .as_flag(log.p, "log.p"))
  if (anyNA(out)) .warn_nan(out, q, alpha, beta, gamma)
  # Keep names and dimensions the way pnorm does (see dmhn).
  if (!is.null(attributes(q)) && !is.object(q) && length(out) == length(q)) {
    attributes(out) <- attributes(q)
  }
  out
}

# Diagnostic hook for the tests and for accuracy audits: evaluates
# the general-case CDF via a single, caller-specified code path,
# bypassing both the special-case shortcuts and the automatic
# series-to-integration fallback used by pmhn().  Scalar arguments
# only.  Internal helper: leading dot follows the package convention
# (cf. .dmhn_cpp, .qmhn_cpp), and @noRd suppresses Rd output.
#' @noRd
.pmhn_force <- function(q, alpha, beta, gamma,
                        method = c("series", "integrate")) {
  method <- match.arg(method)
  .pmhn_force_cpp(as.numeric(q),
                  as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
                  method)
}
