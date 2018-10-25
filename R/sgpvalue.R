################################################################
##	Purpose: 	Compute Second-Generation p-values and delta-gaps
##
##	Function:	sgpvalue
##	Version:	2.1
##
##	Author:		Jeffrey D. Blume and Valerie F. Welty
##	Date:		Oct 16, 2018
################################################################
#
#' Second-Generation p-values
#'
#' @description This function calculates the Second-Generation p-value (SGPV) and associated delta gaps, as introduced in Blume et al. (2018).
#'
#' @param est.lo A (non-empty) numeric vector of lower bounds of interval estimates. These may be finite or they may be \code{-Inf} or \code{Inf}.
#' @param est.hi A (non-empty) numeric vector of upper bounds of interval estimates. These may be finite or they may be \code{-Inf} or \code{Inf}. Must be of same length as \code{est.lo}.
#' @param null.lo A (non-empty) numeric vector of lower bounds of null intervals. These may be finite or they may be \code{-Inf} or \code{Inf}.
#' @param null.hi A (non-empty) numeric vector of upper bounds of null intervals. These may be finite or they may be \code{-Inf} or \code{Inf}. Must be of same length as \code{null.lo}. If not of length 1, must be of same length as \code{est.lo} and \code{est.hi}.
#' @param inf.correction A (usually) small scalar to denote positive but infinitesimally small SGPV. Default = 1e-5. Values that are infinitesimally close to 1 are assigned \code{1-inf.correction}.
#' @param length.warn Warnings toggle. Default is TRUE.
#'
#' @details
#' [Details]
#'
#' Note that the second-generation \emph{p}-value should be calculated on the "symmetric" scale. For example, if the parameter of interest is an odds ratio, the inputs \code{est.lo}, \code{est.hi}, \code{null.lo}, and \code{null.hi} should be on the scale of log odds ratios.
#' @return A list containing the following components:
#' \describe{
#' \item{\code{p.delta}}{the vector of calculated second-generation p-values}
#' \item{\code{delta.gap}}{the vector of calculated delta-gaps. Will be \code{NA} for those with corresponding second-generation p-values not equal to zero.}
#' }
#' @seealso \code{\link{fdrisk}}
#' @keywords
#' @export
#' @examples
#' ## Single inputs
#' > sgpvalue(est.lo = log(1.05), est.hi = log(1.8), null.lo = log(1/1.1), null.hi = log(1.1))
#'
#' $p.delta
#' [1] 0.1220227
#'
#' $delta.gap
#' [1] NA
#'
#' > sgpvalue(est.lo = log(1.3), est.hi = log(1.8), null.lo = log(1/1.1), null.hi = log(1.1))
#'
#' $p.delta
#' [1] 0
#'
#' $delta.gap
#' [1] 1.752741
#'
#' > sgpvalue(est.lo = log(0.97), est.hi = log(1.02), null.lo = log(1/1.1), null.hi = log(1.1))
#'
#' $p.delta
#' [1] 1
#'
#' $delta.gap
#' [1] NA
#'
#' ## Vectors of inputs
#' > lb = c(log(1.05), log(1.3), log(0.97))
#' > ub = c(log(1.8), log(1.8), log(1.02))
#' > sgpvalue(est.lo = lb, est.hi = ub, null.lo = rep(log(1/1.1), 3), null.hi = rep(log(1.1), 3))
#'
#' $p.delta
#' [1] 0.1220227 0.0000000 1.0000000
#'
#' $delta.gap
#' [1] NA        1.752741  NA
#'
#' ## To get only the vector of second-generation p-values or only the vector of delta-gaps:
#' > sgpvalue(est.lo = lb, est.hi = ub, null.lo = rep(log(1/1.1), 3), null.hi = rep(log(1.1), 3))$p.delta
#'
#' [1] 0.1220227 0.0000000 1.0000000
#'
#' > sgpvalue(est.lo = lb, est.hi = ub, null.lo = rep(log(1/1.1), 3), null.hi = rep(log(1.1), 3))$delta.gap
#'
#' [1] NA        1.752741  NA
#'
#' ## Interval bounds may be infinite:
#' > sgpvalue(est.lo = log(1.3), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1))
#'
#' $p.delta
#' [1] 0
#'
#' $delta.gap
#' [1] 0.1670541
#'
#' Warning message:
#'   In sgpvalue(est.lo = log(1.3), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1)) :
#'     At least one interval has infinite length
#'
#' > sgpvalue(est.lo = log(1.05), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1))
#'
#' $p.delta
#' [1] 1e-05
#'
#' $delta.gap
#' [1] NA
#'
#' Warning message:
#'     In sgpvalue(est.lo = log(1.05), est.hi = Inf, null.lo = -Inf, null.hi = log(1.1)) :
#'       At least one interval has infinite length
#'
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr (2018) Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Welty VF, Dupont WD, Greevy RA Jr. An Introduction to Second-Generation \emph{p}-values. The American Statistician. doi.
#'



sgpvalue <- function(est.lo, est.hi, null.lo, null.hi,
			inf.correction=1e-5, length.warn=TRUE){

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

	if ((na.any==TRUE)&length.warn){
		warning('At least one input is NA')}

	if ((na.any==FALSE)&any(est.len<0)&any(null.len<0)&length.warn) {
		warning('At least one interval length is negative')}

	if ((na.any==FALSE)&any(is.infinite(abs(est.len)+abs(null.len)))&length.warn) {
		warning('At least one interval has infinite length') }

	if ((na.any==FALSE)&any(est.len==0,null.len==0)&length.warn) {
		warning('At least one interval has zero length') }

	#### SGPV computation ####
	## assume ordering correc. Fix at end.

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

	## Null Interval is a point (overlap=zero) but is in interval
	p.delta[est.len>0&null.len==0&est.lo<=null.lo&est.hi>=null.hi] <- 1/2

	## One-sided intervals with overlap; overlap == Inf & bottom==Inf
	p.delta[is.infinite(overlap)&is.infinite(bottom)&((est.hi<=null.hi)|(est.lo>=null.lo))] <- 1
	p.delta[is.infinite(overlap)&is.infinite(bottom)&((est.hi>null.hi)|(est.lo<null.lo))] <- 1-inf.correction

	## Interval Estimate entire real line and null is NOT entire real line
	p.delta[est.lo==-Inf&est.hi==Inf] <- 1/2

	## Null Interval entire real line
	p.delta[null.lo==-Inf&null.hi==Inf] <- NA

	    ## ! to do: add warning ""

	## Return NA for nonsense intervals
	p.delta[(est.lo>est.hi)|(null.lo>null.hi)] <- NA

	    ## ! to do: add warning ""

	## Calculate delta gap
	delta.gap <- rep(NA, length(p.delta))
	delta.gap[!is.na(p.delta)&(p.delta==0)] = 0

	gap = (pmax(est.lo, null.lo) - pmin(null.hi, est.hi))
	# gap[!is.na(p.delta)&(p.delta==0)]
	delta = null.len/2 # ifelse(null.len==0, 1/inf.correction, null.len/2)
	  # unscaled delta gap if null has infinite length:
	  delta[null.len==Inf] = 1

	  # unscaled delta gap if null has length zero:
	  delta[null.len==0] = 1

	dg = gap/delta

  delta.gap[!is.na(p.delta)&(p.delta==0)] = dg[!is.na(p.delta)&(p.delta==0)]

	return(list("p.delta"=p.delta, "delta.gap"=delta.gap))

}

###
##
#
