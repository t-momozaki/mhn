# Random generation for the Modified Half-Normal distribution.
# Thin wrapper around .rmhn_cpp (src/mhn_rmhn.cpp).

#' Random Generation from the Modified Half-Normal Distribution
#'
#' Draws random variates from the Modified Half-Normal (MHN) distribution
#' with parameters \code{alpha}, \code{beta}, and \code{gamma}.
#'
#' The MHN density is
#' \deqn{f(x \mid \alpha, \beta, \gamma) =
#'   \frac{2 \beta^{\alpha/2} x^{\alpha-1}
#'   \exp(-\beta x^2 + \gamma x)}{\Psi[\alpha/2, \gamma/\sqrt{\beta}]}
#'   \quad (x > 0)}
#' where \eqn{\Psi[a, z]} is the Fox-Wright Psi function. \code{rmhn} does
#' not evaluate \eqn{\Psi}; the rejection-sampling kernels cancel it out.
#'
#' @param n Non-negative integer giving the number of variates to draw.
#'   \code{n = 0} returns \code{numeric(0)}. Following the base R convention
#'   of \code{rnorm()} and \code{rgamma()}, if \code{length(n) > 1} the
#'   number required is taken to be \code{length(n)}.
#' @param alpha Shape parameter (\eqn{\alpha > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param beta Scale parameter (\eqn{\beta > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param gamma Location parameter (\eqn{\gamma \in R}). Scalar or numeric
#'   vector. Default: 0.
#' @param method Sampling algorithm. One of \code{"auto"} (default),
#'   \code{"rtdr"}, or \code{"sun"}. See Details.
#'
#' @return A numeric vector of length \code{n}. If any of \code{alpha},
#'   \code{beta}, \code{gamma} (after recycling to length \code{n}) is
#'   \code{NA} or non-finite (\code{Inf}, \code{-Inf}, \code{NaN}), the
#'   corresponding output element is \code{NA}.
#'
#' @details
#' The default parameters \code{alpha = 1, beta = 1, gamma = 0} correspond
#' to the half-normal distribution \eqn{\mathrm{HN}(1/\sqrt{2})}.
#'
#' The \code{method} argument selects the rejection sampler:
#' \itemize{
#'   \item \code{"auto"} (default): Uses closed-form special cases when
#'     applicable (\eqn{\gamma \approx 0} -> sqrt-Gamma,
#'     \eqn{\alpha \approx 1} -> truncated normal). Otherwise it selects the
#'     fastest provably correct sampler for the parameter region and the
#'     number of variates drawn per setup: for \eqn{\gamma > 0} it uses Sun
#'     et al. (2023) Algorithm 1 when \eqn{\alpha > 1} and RTDR (Gao & Wang,
#'     2025) when \eqn{\alpha < 1}; for \eqn{\gamma < 0} it uses Sun et al.
#'     Algorithm 3 for small batches and for \eqn{\alpha \ge 10}, and RTDR
#'     for larger batches with \eqn{\alpha < 10}. A batch counts as large at
#'     25 variates per setup, raised to 100 for \eqn{\alpha < 0.1} where the
#'     crossover between the two samplers occurs later. These thresholds were
#'     fixed by benchmarking (see \code{inst/benchmarks/auto_dispatch.R}).
#'   \item \code{"rtdr"}: Use the Relaxed Transformed Density Rejection
#'     method of Gao & Wang (2025) for the general case. The acceptance
#'     probability is bounded below by \eqn{1/e \approx 0.368} uniformly over
#'     the parameter space. The closed-form special cases are still taken
#'     first: a negligible tilt draws from the square root of a Gamma, and
#'     \eqn{\alpha \approx 1} from a truncated normal. Both are exact and
#'     faster than any rejection scheme, so no sampler choice overrides them.
#'     Note: Gao & Wang (2025) use the parameterization
#'     \eqn{(\lambda, \alpha, \beta)} with density proportional to
#'     \eqn{x^{\lambda - 1} \exp(-\alpha x^2 - \beta x)}; the mapping to
#'     the Sun et al. parameterization used here is
#'     \eqn{\lambda \leftrightarrow \alpha},
#'     \eqn{\alpha \leftrightarrow \beta},
#'     \eqn{\beta \leftrightarrow -\gamma} (sign flip on the linear term).
#'   \item \code{"sun"}: Use the Sun et al. (2023) algorithms for the general
#'     case, with the same closed-form special cases taken first.
#'     Algorithm 1 is used when \eqn{\gamma > 0} and \eqn{\alpha > 1};
#'     Algorithm 3 is used when \eqn{\gamma \le 0}. The combination
#'     \eqn{\alpha < 1} with \eqn{\gamma > 0} is unsupported and triggers
#'     an error, unless a special case answers it first.
#' }
#'
#' Vector parameters are recycled to length \code{n} following standard R
#' rules: only the first \code{n} elements of each parameter are used, and any
#' further elements are silently ignored, matching the convention of
#' \code{rnorm}.
#'
#' Internally the setup state of the chosen sampler is reused as long as
#' consecutive \eqn{(\alpha, \beta, \gamma)} triples are equal, so passing
#' parameters grouped by triple is faster than calling \code{rmhn} inside
#' an R loop.
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
#' \emph{Communications in Statistics - Theory and Methods}, 52(5),
#' 1591--1613.
#'
#' Gao, F. & Wang, H.-B. (2025). Generating modified-half-normal random
#' variates by a relaxed transformed density rejection method.
#' \emph{Communications in Statistics - Simulation and Computation}.
#'
#' Robert, C. P. (1995). Simulation of truncated normal variables.
#' \emph{Statistics and Computing}, 5(2), 121--125.
#'
#' @seealso \code{\link{dmhn}}, \code{\link{mhn_mean}}, \code{\link{mhn_var}}
#'
#' @examples
#' set.seed(1)
#' rmhn(10, alpha = 2, beta = 1, gamma = 0.5)
#'
#' # Vector parameters are recycled to length n.
#' set.seed(1)
#' rmhn(5, alpha = c(1, 2, 3, 4, 5))
#'
#' @export
rmhn <- function(n, alpha = 1, beta = 1, gamma = 0,
                 method = c("auto", "rtdr", "sun")) {
  # match.arg is roughly 55% of the cost of rmhn(1, alpha, beta, gamma), the
  # single-variate call a Gibbs sampler makes once per sweep.  Skipping it when
  # the argument was not supplied costs nothing and validates exactly as before
  # for any value that was.
  method <- if (missing(method)) "auto" else match.arg(method)
  # Follow the base R convention shared by rnorm(), runif() and rgamma():
  # a vector n asks for length(n) variates.
  n <- if (length(n) > 1L) length(n) else as.integer(n)[1L]
  .rmhn_cpp(n,
            as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
            method)
}
