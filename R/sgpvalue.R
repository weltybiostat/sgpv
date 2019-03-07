################################################################
##	Purpose: 	Compute Second-Generation p-values and delta-gaps
##
##	Function:	sgpvalue
##	Version:	2.1
##
##	Author:		Jeffrey D. Blume and Valerie F. Welty
##	Date:		November 1, 2018
################################################################
#
#' Second-Generation p-values
#'
#' @description This function computes the second-generation \emph{p}-value (SGPV) and its associated delta gaps, as introduced in Blume et al. (2018).
#'
#' @param est.lo A numeric vector of lower bounds of interval estimates. Values may be finite or \code{-Inf} or \code{Inf}. Must be of same length as \code{est.hi}.
#' @param est.hi A numeric vector of upper bounds of interval estimates. Values may be finite or \code{-Inf} or \code{Inf}. Must be of same length as \code{est.lo}.
#' @param null.lo A numeric vector of lower bounds of null intervals. Values may be finite or \code{-Inf} or \code{Inf}. Must be of same length as \code{null.hi}.
#' @param null.hi A numeric vector of upper bounds of null intervals. Values may be finite or \code{-Inf} or \code{Inf}. Must be of same length as \code{null.lo}.
#' @param inf.correction A small scalar to denote a positive but infinitesimally small SGPV. Default is 1e-5. SGPVs that are infinitesimally close to 1 are assigned \code{1-inf.correction}. This option can only be invoked when one of the intervals has infinite length.
#' @param warnings Warnings toggle. Warnings are on by default.
#'
#' @details Values of \code{NA} or \code{NaN} for \code{est.lo}, \code{est.hi}, \code{null.lo}, or \code{null.lo} will yield a warning and result in a SGPV of \code{NA} or \code{NaN}.
#'
#'When \code{null.hi} and \code{null.lo} are of length 1, the same null interval is used for every interval estimate of [\code{est.lo}, \code{est.hi}]. If \code{null.hi} is not of length 1, its length must match that of \code{est.hi}.
#'
#'When possible, one should compute the second-generation \emph{p}-value on a scale that is symmetric about the null hypothesis. For example, if the parameter of interest is an odds ratio, computations are typically done on the log scale. This keeps the magnitude of positive and negative delta-gaps comparable. Also, recall that the delta-gaps magnitude is not comparable across different null intervals.
#'
#' @return A list containing the following components:
#' \describe{
#' \item{\code{p.delta}}{Vector of second-generation p-values}
#' \item{\code{delta.gap}}{Vector of delta-gaps. Reported as \code{NA} when the corresponding second-generation p-value is not zero.}
#' }
#' @seealso \code{\link{fdrisk}, \link{sgpower}, \link{plotsgpv}}
#' @keywords
#' @export
#' @examples
#'
#' ## Simple example for three estimated log odds ratios but the same null interval
#' lb <- c(log(1.05), log(1.3), log(0.97))
#' ub <- c(log(1.8), log(1.8), log(1.02))
#' sgpv <- sgpvalue(est.lo = lb, est.hi = ub, null.lo = log(1/1.1), null.hi = log(1.1))
#' sgpv$p.delta
#' [1] 0.1220227 0.0000000 1.0000000
#' sgpv$delta.gap
#' [1]       NA 1.752741       NA
#'
#' ## Works with infinte interval bounds
#' sgpvalue(est.lo = log(1.3), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1))
#'
#' $p.delta
#' [1] 0
#' $delta.gap
#' [1] 0.1670541
#'
#' Warning message:
#' In sgpvalue(est.lo = log(1.3), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1))
#' At least one interval has infinite length
#'
#' sgpvalue(est.lo = log(1.05), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1))
#'
#' $p.delta
#' [1] 1e-05
#' $delta.gap
#' [1] NA
#'
#' Warning message:
#' In sgpvalue(est.lo = log(1.05), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1)) :
#' At least one interval has infinite length
#'
#' ## Example t-test with simulated data
#' ## (not run)
#' set.seed(1776)
#' x1 <- rnorm(15,mean=0,sd=2) ; x2 <- rnorm(15,mean=3,sd=2)
#' ci <- t.test(x1,x2)$conf.int[1:2]
#' sgpvalue(est.lo = ci[1], est.hi = ci[2], null.lo = -1, null.hi = 1)
#'
#' set.seed(2019)
#' x1 <- rnorm(15,mean=0,sd=2) ; x2 <- rnorm(15,mean=3,sd=2)
#' ci <- t.test(x1,x2)$conf.int[1:2]
#' sgpvalue(est.lo = ci[1], est.hi = ci[2], null.lo = -1, null.hi = 1)
#'
#' ## Simulated two-group dichotomous data for different parameters
#' set.seed(1492)
#' n1 <- n2 <- 30
#' x1 <- rbinom(1,size=n1,p=0.15) ; x2 <- rbinom(1,size=n2,p=0.50)
#'
#' # On the difference in proportions
#' ci.p  <- prop.test(c(x1,x2),n=c(n1,n2))$conf.int[1:2]
#' sgpvalue(est.lo = ci.p[1], est.hi = ci.p[2], null.lo = -0.2, null.hi = 0.2)
#'
#' # On the log odds ratio scale
#' a <- x1 ; b <- x2 ; c <- n1-x1 ; d <- n2-x2
#' ci.or <- log(a*d/(b*c)) + c(-1,1)*1.96*sqrt(1/a+1/b+1/c+1/d)	# Delta-method SE for log odds ratio
#' sgpvalue(est.lo = ci.or[1], est.hi = ci.or[2], null.lo = log(1/1.5), null.hi = log(1.5))
#'
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. In press. https://doi.org/10.1080/00031305.2018.1537893
#'



sgpvalue <- function(est.lo, est.hi, null.lo, null.hi,
			inf.correction=1e-5, warnings=TRUE){

  #### Errors
  if(length(null.lo)!=length(null.hi)) {
    stop('null.lo and null.hi of different lengths')
  }

  if(length(est.lo)!=length(est.hi)) {
    stop('est.lo and est.hi of different lengths')
  }

  ####
  if(length(null.lo)==1) {
    null.lo = rep(null.lo, length(est.lo))
    null.hi = rep(null.hi, length(est.lo))
  }

	#### Compute Interval Lengths ####
	est.len   <- est.hi-est.lo
	null.len  <- null.hi-null.lo

	#### Warnings ####
	na.any = (any(is.na(est.lo))|any(is.na(est.hi))|any(is.na(null.lo))|any(is.na(null.hi)))

	if ((na.any==TRUE)&warnings){
		warning('At least one input is NA')}

	if ((na.any==FALSE)&any(est.len<0)&any(null.len<0)&warnings) {
		warning('At least one interval length is negative')}

	if ((na.any==FALSE)&any(is.infinite(abs(est.len)+abs(null.len)))&warnings) {
		warning('At least one interval has infinite length') }

	if ((na.any==FALSE)&any(est.len==0,null.len==0)&warnings) {
		warning('At least one interval has zero length') }

	#### SGPV computation ####

	overlap <- pmin(est.hi,null.hi)-pmax(est.lo,null.lo)
	overlap <- pmax(overlap,0)		## Negative implies disjoint

	bottom 	<- pmin(2*null.len,est.len)

	p.delta <- overlap/bottom

	#### Zero-length & Infinite-length intervals ####

	## Overwrite NA and NaN due to bottom = Inf
	p.delta[overlap==0]   <- 0

	## Overlap finite & non-zero but bottom = Inf
	p.delta[overlap!=0&is.finite(overlap)&is.infinite(bottom)] <- inf.correction

	## Interval estimate is a point (overlap=zero) but can be in null or equal null pt
	p.delta[est.len==0&null.len>=0&est.lo>=null.lo&est.hi<=null.hi] <- 1

	## Null interval is a point (overlap=zero) but is in interval estimate
	p.delta[est.len>0&null.len==0&est.lo<=null.lo&est.hi>=null.hi] <- 1/2

	## One-sided intervals with overlap; overlap == Inf & bottom==Inf
	p.delta[is.infinite(overlap)&is.infinite(bottom)&((est.hi<=null.hi)|(est.lo>=null.lo))] <- 1
	p.delta[is.infinite(overlap)&is.infinite(bottom)&((est.hi>null.hi)|(est.lo<null.lo))] <- 1-inf.correction

	## Interval estimate is entire real line and null interval is NOT entire real line
	p.delta[est.lo==-Inf&est.hi==Inf] <- 1/2

	## Null interval is entire real line
	p.delta[null.lo==-Inf&null.hi==Inf] <- NA

	if (any(null.lo==-Inf)&any(null.hi==Inf)&warnings) {
		warning('at least one null interval is entire real line') }

	## Return NA for nonsense intervals
	p.delta[(est.lo>est.hi)|(null.lo>null.hi)] <- NA

	if ((any(est.lo>est.hi)|any(null.lo>null.hi))&warnings) {
		warning('Some interval limits likely reversed') }

	## Calculate delta gap
	delta.gap <- rep(NA, length(p.delta))
	delta.gap[!is.na(p.delta)&(p.delta==0)] = 0

	gap = (pmax(est.lo, null.lo) - pmin(null.hi, est.hi))

	delta = null.len/2

	# Report unscaled delta gap if null has infinite length
  	delta[null.len==Inf] = 1

  	# Report unscaled delta gap if null has length zero
  	delta[null.len==0] = 1

	dg = gap/delta

  	delta.gap[!is.na(p.delta)&(p.delta==0)] = dg[!is.na(p.delta)&(p.delta==0)]

	return(list("p.delta"=p.delta, "delta.gap"=delta.gap))

}

###
##
#
