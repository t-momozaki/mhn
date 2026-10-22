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
#'   \quad (x > 0)}{f(x | alpha, beta, gamma) = 2 beta^(alpha/2) x^(alpha-1)
#'   exp(-beta x^2 + gamma x) / Psi[alpha/2, gamma/sqrt(beta)],  x > 0}
#' where \eqn{\Psi[a, z] = \sum_{n \ge 0} \Gamma(a + n/2) z^n / n!}{Psi[a, z] =
#' sum_{n >= 0} Gamma(a + n/2) z^n / n!} is the Fox-Wright function
#' \eqn{{}_1\Psi_1[(a, 1/2); (1, 0); z]}{1Psi1[(a, 1/2); (1, 0); z]}.
#' \code{rmhn} does not evaluate \eqn{\Psi}; the rejection-sampling kernels
#' cancel it out.  The one exception is the weight of the origin,
#' \eqn{1 / (1 + w)}, that it draws for a subnormal shape (see Subnormal
#' shapes below).
#'
#' @param n Number of variates to draw: a non-negative number no larger than
#'   \code{.Machine$integer.max}, truncated to an integer, so that
#'   \code{n = 0} returns \code{numeric(0)}. Following the base R convention
#'   of \code{rnorm()} and \code{rgamma()}, if \code{length(n)} is not 1 the
#'   number required is taken to be \code{length(n)}, so that a zero-length
#'   \code{n} also returns \code{numeric(0)}; \code{NULL} is an error.
#' @param alpha Shape parameter (\eqn{\alpha > 0}). Scalar or numeric vector.
#'   Default: 1.
#' @param beta Rate-type parameter (\eqn{\beta > 0}), the coefficient of the
#'   Gaussian tail \eqn{\exp(-\beta x^2)}{exp(-beta x^2)}; larger values pull
#'   mass toward zero. Scalar or numeric vector. Default: 1.
#' @param gamma Tilt parameter (\eqn{\gamma \in R}{gamma in R}), the
#'   coefficient of the exponential tilt \eqn{\exp(\gamma x)}{exp(gamma x)}.
#'   Scalar or numeric vector. Default: 0.
#' @param method Sampling algorithm. One of \code{"auto"} (default),
#'   \code{"rtdr"}, or \code{"sun"}. See Details.
#'
#' @return A numeric vector of length \code{n}. As in
#'   \code{\link[stats]{rnorm}}, an element whose parameters (after
#'   recycling to length \code{n}) are \code{NA} or \code{NaN}, or lie
#'   outside the parameter space (\code{alpha} or \code{beta} not positive
#'   and finite, or \code{gamma} not finite), is \code{NaN}, with a single
#'   warning "NAs produced", and draws no random numbers.  A zero-length
#'   \code{alpha}, \code{beta} or \code{gamma} makes every element
#'   \code{NA}, with the same warning.  The parameters must be numeric
#'   (double, integer or logical): a factor, a character string, a complex
#'   number, a list or \code{NULL} is an error.  If a rejection sampler
#'   accepts none of 1000 proposals for an element (for a subnormal shape,
#'   none of 1000 draws of the proxy described under Subnormal shapes in
#'   Details, each of up to 1000 proposals),
#'   that element is \code{NaN}, and a single warning reports how many draws
#'   were affected;
#'   the acceptance rates of the samplers make this vanishingly unlikely.
#'   Where the distribution lies beyond the largest double precision number,
#'   and at some other extreme values, \code{rmhn} instead stops with an
#'   error for the whole call (see Details).
#'
#' @details
#' The default parameters \code{alpha = 1, beta = 1, gamma = 0} correspond
#' to the half-normal distribution \eqn{\mathrm{HN}(1/\sqrt{2})}{HN(1/sqrt(2))}.
#'
#' \subsection{Choice of sampler}{
#' Two closed-form special cases are taken first, whatever \code{method}:
#' a negligible tilt, \eqn{|\gamma|/\sqrt{\beta}}{|gamma|/sqrt(beta)} within
#' \code{sqrt(.Machine$double.eps)} of 0, is drawn as the square root of a
#' Gamma variate, and a shape \eqn{\alpha}{alpha} within
#' \code{sqrt(.Machine$double.eps)} of 1 from a truncated normal. Both are
#' exact samplers of their special-case law, which next to
#' \eqn{\gamma = 0}{gamma = 0} or \eqn{\alpha = 1}{alpha = 1} is within a
#' total variation distance of \eqn{8.6 \times 10^{-9}}{8.6e-9} of the law
#' requested. The one exception is a negligible tilt above
#' \eqn{\alpha = 2^{60}}{alpha = 2^60}, which is drawn by the Normal
#' proposal of Sun et al. (2023) Algorithm 1 centred at the mode, or by RTDR
#' under \code{method = "rtdr"} (Sun et al. state Algorithm 1 for
#' \eqn{\gamma > 0}{gamma > 0}; centred at the mode, its Normal envelope
#' holds for either sign). Otherwise \code{method} selects the
#' rejection sampler:
#' \itemize{
#'   \item \code{"auto"} (default): the faster of the provably correct
#'     samplers, as measured by timing them against each other. For
#'     \eqn{\alpha < 1}{alpha < 1} and \eqn{\gamma > 0}{gamma > 0} it uses
#'     RTDR. Otherwise it chooses between Sun et al. Algorithm 1
#'     (\eqn{\gamma > 0}{gamma > 0}) or Algorithm 3
#'     (\eqn{\gamma < 0}{gamma < 0}) and RTDR by the shape, the tilt
#'     \eqn{\gamma/\sqrt{\beta}}{gamma/sqrt(beta)} and the number of variates
#'     per setup (see Vector parameters and setups); Section 7 of
#'     \code{vignette("theory", package = "mhn")} gives the decision tree.
#'     The choice changes the speed and the stream of draws
#'     (see Reproducibility), not the law drawn.
#'   \item \code{"rtdr"}: the Relaxed Transformed Density Rejection method of
#'     Gao & Wang (2025). Gao & Wang bound its acceptance probability below
#'     by \eqn{1/e \approx 0.368}{1/e ~ 0.368} for every parameter value,
#'     provided the contact points of the envelope lie in the band of their
#'     Eq. (7); the search for them does not always reach that band,
#'     and the acceptance measured stays above the bound there too. Where
#'     RTDR cannot place a law of tiny shape under a strong negative tilt,
#'     Sun et al. Algorithm 3 is used instead (see Limits of double
#'     precision). Note: Gao & Wang (2025) use the parameterization
#'     \eqn{(\lambda, \alpha, \beta)} with density proportional to
#'     \eqn{x^{\lambda - 1} \exp(-\alpha x^2 - \beta x)}{x^(lambda - 1)
#'     exp(-alpha x^2 - beta x)}; the mapping to the Sun et al.
#'     parameterization used here is
#'     \eqn{\lambda \leftrightarrow \alpha}{lambda <-> alpha},
#'     \eqn{\alpha \leftrightarrow \beta}{alpha <-> beta},
#'     \eqn{\beta \leftrightarrow -\gamma}{beta <-> -gamma} (sign flip on
#'     the linear term).
#'   \item \code{"sun"}: the algorithms of Sun et al. (2023): Algorithm 1
#'     for \eqn{\gamma > 0} and \eqn{\alpha > 1}, and Algorithm 3 for
#'     \eqn{\gamma < 0}, up to \eqn{\alpha = 2^{59}}{alpha = 2^59}, beyond
#'     which the call stops with an error (see Limits of double precision).
#'     The combination \eqn{\alpha < 1} with \eqn{\gamma > 0} is
#'     unsupported: unless a special case answers it first, the call stops
#'     with the error \code{Sun method not available for alpha<1 and
#'     gamma>0; use method="rtdr"}.
#' }
#' }
#'
#' \subsection{Vector parameters and setups}{
#' Vector parameters are recycled to length \code{n} following standard R
#' rules: only the first \code{n} elements of each parameter are used, and any
#' further elements are silently ignored, matching the convention of
#' \code{rnorm}.
#'
#' The setup of the chosen sampler is built once for each run of equal
#' consecutive \eqn{(\alpha, \beta, \gamma)} triples. Short parameter
#' vectors repeat their triples after a period, the least common multiple
#' of the lengths of \code{alpha}, \code{beta} and \code{gamma}; when that
#' period is at most 256 draws and they are recycled over at least four
#' such periods, a setup is kept for each position within the period. So
#' passing parameters grouped by triple, or as short recycled vectors, is
#' faster than calling \code{rmhn} inside an R loop, and a long run of one
#' triple is drawn fastest by its own scalar call.
#'
#' The number of variates per setup, from which \code{method = "auto"}
#' chooses its sampler, and which therefore changes the draws for a fixed
#' seed, is \code{n} for scalar parameters, and \code{n} divided by the
#' period, rounded down, where a setup is kept for each position within the
#' period. Otherwise it is 1: parameters given as vectors of length
#' \code{n} count one per setup, as in a Gibbs sampler, even when equal
#' triples are grouped into runs, and such a call draws what the loop of
#' single draws would. Grouping saves setups but does not change this
#' number, and so not the choice of sampler.
#' }
#'
#' \subsection{Reproducibility}{
#' Under one seed, \code{method = "rtdr"} and \code{method = "sun"} draw one
#' stream whatever \code{n}: the first \code{k} draws of \code{rmhn(n, ...)}
#' are those of \code{rmhn(k, ...)}, and draws split over several calls are
#' those of one call.  The default method does not keep this, as it chooses
#' its sampler from the number of variates per setup (see above):
#' \code{rmhn(4, 5, 1, 2)} and the first four draws of
#' \code{rmhn(100, 5, 1, 2)} come from different samplers.  Its stream
#' therefore changes with \code{n}, with how the parameters are arranged,
#' and with a version of the package that moves a threshold, although the law
#' drawn does not.  Where the draws must not depend on \code{n}, as when a
#' pilot run is compared with the main one, set \code{method = "rtdr"} or
#' \code{method = "sun"}.
#' }
#'
#' \subsection{Limits of double precision}{
#' Double precision limits which laws can be sampled, and how finely. Every
#' sampler draws from the exact law and rounds the draw to double precision,
#' and a draw below the smallest positive double is returned as 0. A
#' rounding error that would shift every draw alike is kept to a few times
#' \eqn{10^{-7}}{1e-7} of a standard deviation. For this reason the
#' square-root-Gamma closed form is used only up to
#' \eqn{\alpha = 2^{60}}{alpha = 2^60}, and Sun et al. Algorithm 3 only up
#' to \eqn{\alpha = 2^{59}}{alpha = 2^59} (about
#' \eqn{5.8 \times 10^{17}}{5.8e17}). Above
#' \eqn{\alpha = 2^{59}}{alpha = 2^59} a negative tilt that is not
#' negligible is drawn by RTDR under the default method, while
#' \code{method = "sun"} stops with an error; a positive tilt is drawn by
#' Algorithm 1 under both.
#'
#' The RTDR sampler stops with an error for the whole call once the standard
#' deviation is no larger than the spacing of the double precision numbers
#' at the mode, so that its draws could take only a few distinct values:
#' from \eqn{\gamma/\sqrt{\beta}}{gamma/sqrt(beta)} of about
#' \eqn{9 \times 10^{15}}{9e15}
#' (\eqn{2^{53}}{2^53} at \eqn{\beta = 1}{beta = 1}),
#' and, for \eqn{\gamma \le 0}{gamma <= 0}, from \eqn{\alpha} between
#' \eqn{2^{103}}{2^103} and \eqn{2^{106}}{2^106}
#' (\eqn{1.0 \times 10^{31}}{1.0e31} to \eqn{8.1 \times 10^{31}}{8.1e31}),
#' depending on the tilt and on where the mode falls between two powers of
#' two. There the law is normal to far below double precision, with mean
#' \eqn{m} = \code{mhn_mode(alpha, beta, gamma)} and variance
#' \eqn{1/((\alpha - 1)/m^2 + 2\beta)}{1/((alpha - 1)/m^2 + 2 beta)}, and
#' can be drawn from directly, as the error message suggests. Sun et al.
#' Algorithm 1, which the default method uses for such laws with
#' \eqn{\gamma > 0} and \eqn{\alpha > 1}, and the truncated normal of
#' \eqn{\alpha \approx 1}{alpha ~ 1} do not stop there: their draws
#' are placed at the mode to within a few units in the last place, so that
#' they take only a few distinct values near it.
#'
#' RTDR cannot place the law in three cases, all under a negative tilt:
#' where \eqn{-\gamma/\sqrt{\beta}}{-gamma/sqrt(beta)} passes the largest
#' double; for \eqn{\alpha \le 1}{alpha <= 1}, where \eqn{\alpha} is below
#' about \eqn{2.5 \times 10^{-324}}{2.5e-324} times
#' \eqn{-\gamma/\sqrt{\beta}}{-gamma/sqrt(beta)}, which needs
#' \eqn{\alpha}{alpha} below \eqn{4.4 \times 10^{-16}}{4.4e-16}, and almost
#' all the mass of \eqn{X} lies below the smallest positive double; and for
#' a shape of a few, where \eqn{-\gamma/\sqrt{\beta}}{-gamma/sqrt(beta)}
#' comes within a small factor of the largest double. In the second case,
#' and in the first for \eqn{\alpha < 10^{-15}}{alpha < 1e-15},
#' \code{method = "rtdr"} draws from Algorithm 3, which samples the law
#' exactly there, except for a subnormal shape, which is drawn as described
#' under Subnormal shapes; otherwise it stops with
#' an error that names the Gamma distribution with shape \eqn{\alpha} and
#' rate \eqn{-\gamma}{-gamma}, which the law is to double precision. The
#' default method uses Algorithm 3 in all three cases, except where
#' \eqn{-\gamma/\sqrt{\beta}}{-gamma/sqrt(beta)} passes the largest double
#' at a shape above \eqn{2^{59}}{2^59}; there it stops as RTDR does. Where
#' the distribution lies beyond the largest double precision number, no draw
#' can represent it, and every sampler, the closed forms included, stops
#' with an error.
#' }
#'
#' \subsection{Subnormal shapes}{
#' A subnormal shape, below the smallest normal double
#' \eqn{2.2 \times 10^{-308}}{2.2e-308}, is drawn as the mixture its law is
#' (see the Details of \code{\link{dmhn}}): 0 with probability
#' \eqn{1 / (1 + w)}, the probability of a value below the smallest positive
#' double, and otherwise from the peak, through a proxy: the law of the
#' shape \eqn{2^{-60}}{2^-60} above
#' \eqn{\sqrt{\beta} x = 1}{sqrt(beta) x = 1}, drawn by the sampler
#' \code{method} selects for that shape.  The choice between 0 and the peak
#' compares one uniform variate with \eqn{1 / (1 + w)}, and so has the
#' resolution of the generator, \eqn{2^{-32}}{2^-32} for R's default.  Where
#' \eqn{1 / (1 + w)} is 1 in double precision, below a tilt of about 52 (53
#' at the smallest shapes), every draw is 0 and uses no random number; a
#' shape from \eqn{2^{-1022}}{2^-1022} up is drawn by the samplers above,
#' which use random numbers there too, so that the stream after the draws
#' moves as the shape crosses \eqn{2^{-1022}}{2^-1022} although every draw
#' is 0 on either side.  A refusal names the shape asked for.
#' }
#'
#' The algorithms, the timings behind the choice of \code{method = "auto"},
#' and the handling of extreme parameters are described in
#' \code{vignette("theory", package = "mhn")} (Section 7 and Appendix A).
#'
#' @references
#' Sun, J., Kong, M., & Pal, S. (2023). The Modified-Half-Normal
#' distribution: Properties and an efficient sampling scheme.
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
  # match.arg costs more than the sampling itself in a single-variate call,
  # the call a Gibbs sampler makes once per sweep, so a missing argument and
  # an exact choice skip it; anything else, including a partial match, still
  # goes through match.arg.
  method <- if (missing(method)) {
    "auto"
  } else if (is.character(method) && length(method) == 1L && !is.na(method) &&
             (method == "auto" || method == "rtdr" || method == "sun")) {
    method
  } else {
    match.arg(method)
  }
  if (!(is.double(alpha) && is.double(beta) && is.double(gamma))) {
    .check_numeric(alpha, beta, gamma)
  }
  # Follow the base R convention shared by rnorm(), runif() and rgamma():
  # a vector n asks for length(n) variates, one of length zero included, and
  # NULL is an error.  n goes on as a double, so that one beyond the integer
  # range is refused without a coercion warning on the way.
  n <- if (is.null(n) || length(n) == 1L) as.double(n)[1L] else length(n)
  .rmhn_cpp(n,
            as.numeric(alpha), as.numeric(beta), as.numeric(gamma),
            method)
}
