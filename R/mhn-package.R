#' mhn: The Modified Half-Normal Distribution
#'
#' Provides density (\code{\link{dmhn}}), distribution function
#' (\code{\link{pmhn}}), quantile function (\code{\link{qmhn}}),
#' random generation (\code{\link{rmhn}}), and moments / mode of the
#' Modified Half-Normal (MHN) distribution.
#'
#' \code{vignette("introduction", package = "mhn")} shows the functions in
#' use, and \code{vignette("theory", package = "mhn")} sets out the theory
#' they implement; its Appendix A describes the numerical methods.
#'
#' @section MHN Distribution:
#' The MHN distribution has probability density function
#' \deqn{f(x \mid \alpha, \beta, \gamma) =
#'   \frac{2\beta^{\alpha/2} x^{\alpha-1}
#'   \exp(-\beta x^2 + \gamma x)}{\Psi[\alpha/2, \gamma/\sqrt{\beta}]}}{
#'   f(x | alpha, beta, gamma) = 2 beta^(alpha/2) x^(alpha-1)
#'   exp(-beta x^2 + gamma x) / Psi[alpha/2, gamma/sqrt(beta)]}
#' for \eqn{x > 0}, where
#' \eqn{\Psi[a, z] = \sum_{n \ge 0} \Gamma(a + n/2) z^n / n!}{Psi[a, z] =
#' sum_{n >= 0} Gamma(a + n/2) z^n / n!} is the Fox-Wright function
#' \eqn{{}_1\Psi_1[(a, 1/2); (1, 0); z]}{1Psi1[(a, 1/2); (1, 0); z]}.
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{alpha}}{Shape parameter (\eqn{\alpha > 0}).
#'     Controls the \eqn{x^{\alpha-1}}{x^(alpha-1)} term.}
#'   \item{\code{beta}}{Rate-type parameter (\eqn{\beta > 0}).
#'     Controls the Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)};
#'     larger values pull mass toward zero.}
#'   \item{\code{gamma}}{Tilt parameter (\eqn{\gamma \in R}{gamma in R}).
#'     Controls the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)};
#'     larger values move mass away from zero, but \eqn{\gamma} does not
#'     translate the law, whose support stays \eqn{(0, \infty)}{(0, Inf)}.}
#' }
#'
#' @section Special Cases:
#' \describe{
#'   \item{\eqn{\gamma = 0}}{Square-root Gamma distribution
#'     (\eqn{X^2 \sim
#'     \textrm{Gamma}(\alpha/2, \beta)}{X^2 ~ Gamma(alpha/2, beta)}, with rate
#'     \eqn{\beta}{beta})
#'     (Sun et al. 2023, Lemma 6a).}
#'   \item{\eqn{\alpha = 1}}{Truncated normal distribution
#'     \eqn{\textrm{TN}(\gamma/(2\beta),
#'     1/\sqrt{2\beta}, 0, \infty)}{TN(gamma/(2 beta), 1/sqrt(2 beta), 0, Inf)}
#'     (Sun et al. 2023, Lemma 6b).}
#'   \item{\eqn{\alpha = 1, \gamma = 0}}{Half-normal distribution
#'     \eqn{\textrm{HN}(1/\sqrt{2\beta})}{HN(1/sqrt(2 beta))}
#'     (Sun et al. 2023, Lemma 6c).}
#' }
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal distribution:
#' Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613. \doi{10.1080/03610926.2021.1934700}
#'
#' Gao, F., & Wang, H.-B. (2025). Generating modified-half-normal random
#' variates by a relaxed transformed density rejection method.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--25.
#' Advance online publication. \doi{10.1080/03610918.2025.2524551}
#'
#' Robert, C. P. (1995). Simulation of truncated normal variables.
#' \emph{Statistics and Computing}, 5(2), 121--125.
#' \doi{10.1007/BF00143942}
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib mhn, .registration = TRUE
#' @keywords internal
"_PACKAGE"

# Read a d/p/q logical flag through as.logical(), so that log = 1 behaves as
# log = TRUE rather than silently selecting the other branch.  Base R's
# distribution functions read the flag as an integer instead, so they take
# lower.tail = 0.5 as FALSE where this takes it as TRUE, and NA as TRUE.  A
# flag that gives NA, from NA, NULL or a string such as "1", is an error here:
# a flag that cannot be resolved should not quietly pick a branch.
.as_flag <- function(x, name) {
  v <- as.logical(x)[1L]
  if (is.na(v)) stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
  v
}

# Stop unless each argument is stored as numbers, as base R's distribution
# functions do: a double, integer or logical vector is taken, and so is a
# classed one such as a Date, as its numbers; a factor, a string, a complex
# number, a list or NULL is an error.  A factor would otherwise be used through
# its level codes.  is.numeric() would not do, being FALSE for a Date.  x is
# the first argument of dmhn, pmhn and qmhn, named x_name; rmhn has none.
# The callers test the usual case, every argument a double, inline, and call
# this only when that test fails: the inline test costs a small part of a
# call to this function.
.check_numeric <- function(alpha, beta, gamma, x = 0, x_name = "x") {
  stored_as_numbers <- function(v) is.double(v) || is.integer(v) || is.logical(v)
  if (!stored_as_numbers(x)) {
    stop("'", x_name, "' must be numeric", call. = FALSE)
  }
  if (!stored_as_numbers(alpha)) stop("'alpha' must be numeric", call. = FALSE)
  if (!stored_as_numbers(beta)) stop("'beta' must be numeric", call. = FALSE)
  if (!stored_as_numbers(gamma)) stop("'gamma' must be numeric", call. = FALSE)
  invisible()
}

# The same for the moment functions and mhn_mode, which take one triple.
.check_single <- function(alpha, beta, gamma) {
  if (is.double(alpha) && is.double(beta) && is.double(gamma) &&
      length(alpha) == 1L && length(beta) == 1L && length(gamma) == 1L) {
    return(invisible())
  }
  single_number <- function(v) {
    (is.double(v) || is.integer(v) || is.logical(v)) && length(v) == 1L
  }
  if (!single_number(alpha)) {
    stop("'alpha' must be a single number", call. = FALSE)
  }
  if (!single_number(beta)) stop("'beta' must be a single number", call. = FALSE)
  if (!single_number(gamma)) {
    stop("'gamma' must be a single number", call. = FALSE)
  }
  invisible()
}

# Warn, as base R's d/p/q functions do, when an element of the result is NaN
# although none of the inputs it was computed from is NA or NaN: a parameter
# outside the parameter space, or for qmhn a probability outside [0, 1].  The
# warning names the user's call.  Called only when anyNA(out), a scan without
# allocation, so that a result without NA costs nothing more.
.warn_nan <- function(out, x, alpha, beta, gamma) {
  n <- length(out)
  produced <- is.nan(out)
  for (v in list(x, alpha, beta, gamma)) {
    produced <- produced & !is.na(rep_len(v, n))
  }
  if (any(produced)) warning(simpleWarning("NaNs produced", sys.call(-1L)))
}
