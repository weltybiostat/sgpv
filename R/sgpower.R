################################################################
##	Purpose: 	Compute power/type I error
##            for Second-Generation p-values approach
##
##	Function:	sgpower
##	Version:	1.2
##
##	Author:		Valerie F. Welty and Jeffrey D. Blume
##	Date:		December 6, 2018
################################################################
#
#' Power functions for Second-Generation p-values.
#'
#' @description
#' @param true The true value for the parameter of interest at which to calculate power. Note that this is on the absolute scale of the parameter, and not the standard deviation or standard error scale.
#' @param null.lo The lower bound of the indifference zone (null interval) upon which the second-generation \emph{p}-value is based
#' @param null.hi The upper bound for the indifference zone (null interval) upon which the second-generation \emph{p}-value is based
#' @param std.err Standard error for the distribution of the estimator for the parameter of interest. Note that this is the standard deviation for the estimator, not the standard deviation parameter for the data itself. This will be a function of the sample size(s).
#' @param interval.type Class of interval estimate used for calculating the SGPV. Options are \code{confidence} for a \eqn{(1-\alpha)100}\% confidence interval and \code{likelihood} for a \eqn{1/k} likelihood support interval (\code{credible} not yet supported)
#' @param interval.level Level of interval estimate. If \code{interval.type} is \code{confidence}, the level is \eqn{\alpha}. If \code{interval.type} is \code{likelihood}, the level is \eqn{1/k} (not \eqn{k}).

#' @details
#' @return A list containing the following components:
#' \describe{
#' \item{\code{power.alt}}{Probability of SGPV = 0 calculated assuming the parameter is equal to \code{true}. That is, \code{power.alt}\eqn{ = P(SGPV = 0 | \theta = } \code{true}).}
#' \item{\code{power.inc}}{Probability of 0 < SGPV < 1 calculated assuming the parameter is equal to \code{true}. That is, \code{power.inc}\eqn{ = P(0 < SGPV < 1 | \theta = } \code{true}).}
#' \item{\code{power.null}}{Probability of SGPV = 1 calculated assuming the parameter is equal to \code{true}. That is, \code{power.null}\eqn{ = P(SGPV = 1 | \theta = } \code{true}).}
#' \item{\code{`type I error summaries`}}{Named vector that includes different ways the type I error may be summarized for an interval null hypothesis. \code{min} is the minimum type I error over the range (\code{null.lo}, \code{null.hi}), which occurs at the midpoint of (\code{null.lo}, \code{null.hi}). \code{max} is the maximum type I error over the range (\code{null.lo}, \code{null.hi}), which occurs at the boundaries of the null hypothesis, \code{null.lo} and \code{null.hi}. \code{mean} is the average type I error (unweighted) over the range (\code{null.lo}, \code{null.hi}). If \eqn{0} is included in the null hypothesis region, then \code{`type I error summaries`} also contains \code{at 0}, the type I error calculated assuming the true parameter value \eqn{\theta} is equal to \eqn{0}.}
#' }
#' @seealso \code{\link{fdrisk}, \link{sgpvalue}}
#' @keywords
#' @export
#' @examples
#' sgpower(true=2, null.lo=-1, null.hi=1, std.err=1, interval.type='confidence', 'interval.level'=0.05)
#'
#' sgpower(true=0, null.lo=-1, null.hi=1, std.err=1, interval.type='confidence', 'interval.level'=0.05)
#'
#' # plot the power curve
#' sigma = 5
#' n = 20
#' theta = seq(-10, 10, by=0.1)
#' power = sgpower(true=theta, null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), interval.type='confidence', 'interval.level'=0.05)$power.alt
#' plot(theta, power, type='l', ylab='power')
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. In press. https://doi.org/10.1080/00031305.2018.1537893
#'

sgpower <- function(true, null.lo, null.hi, std.err=1, interval.type, interval.level) {

  if (!(interval.type %in% c('confidence', 'likelihood'))) stop("Parameter `interval.type` must be one of the following: \n  * confidence \n  * likelihood \n  (credible not currently supported for power calculations)")

  if(interval.type == 'confidence') Z = qnorm(1-interval.level/2)
  if(interval.type == 'likelihood') Z = qnorm(1-2*pnorm(-sqrt(2*log(1/interval.level)))/2)

  ## P(SGPV=0 | true )    (see Blume et al. (2018) eq.(S4) for CI/LSI)
  power.0 = pnorm(null.lo/std.err - true/std.err - Z) + pnorm(-null.hi/std.err + true/std.err - Z)

  ## P(SGPV=1 | true )    (see Blume et al. (2018) eq.(S7) for CI/LSI)
  # -> only for symmetric null hypothesis
  if((null.hi-null.lo) >= 2*Z*std.err) {
    power.1 = pnorm(null.hi/std.err - true/std.err - Z) - pnorm(null.lo/std.err - true/std.err + Z)
  }
  if((null.hi-null.lo) < 2*Z*std.err) {
    power.1 = 0
  }

  ## P(0<SGPV<1 | true)   (see Blume et al. (2018) eq.(S8, S9) for CI/LSI)
  # -> only for symmetric null hypothesis
  if((null.hi-null.lo) <= 2*Z*std.err) {
    power.inc = 1 - pnorm(null.lo/std.err - true/std.err - Z) - pnorm(-null.hi/std.err + true/std.err - Z)
  }
  if((null.hi-null.lo) > 2*Z*std.err) {
    power.inc = 1 - (pnorm(null.lo/std.err - true/std.err - Z) + pnorm(-null.hi/std.err + true/std.err - Z)) - (pnorm(null.hi/std.err - true/std.err - Z) - pnorm(null.lo/std.err - true/std.err + Z))
  }


  ## check
  if(any(round(power.0+power.inc+power.1,7) != 1)) warning(paste0('error: power.0+power.inc+power.1 != 1 for indices ', paste(which(round(power.0+power.inc+power.1,7) != 1), collapse=', ')))

  ## bonus: type I error summaries
  pow0 = function(x) pnorm(null.lo/std.err - x/std.err - Z) + pnorm(-null.hi/std.err + x/std.err - Z)

  minI = pow0((null.lo+null.hi)/2)
  maxI = pow0(null.lo)
  avgI = 1/(null.hi-null.lo)*integrate(f=pow0, lower=null.lo, upper=null.hi)$value

  typeI = c('min'=minI, 'max'=maxI, 'mean'=avgI)
  if(null.lo<=0 & 0<=null.hi) {
    typeI = c('at 0'=pow0(0), 'min'=minI, 'max'=maxI, 'mean'=avgI)
  }

  ## more bonus: P(inc | null) (basically analogous to type II error but for confirmation)
  ## (TBD)

  return(list('power.alt' = power.0, 'power.inc' = power.inc, 'power.null' = power.1, 'type I error summaries' = typeI))

}


