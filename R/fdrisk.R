################################################################
##	Purpose: 	Compute False Discovery or Confirmation risk
##            for Second-Generation p-values approach
##
##	Function:	fdrisk
##	Version:	1.1
##
##	Author:		Valerie F. Welty and Jeffrey D. Blume
##	Date:		November 8, 2018
################################################################
#
#' False Discovery Risk for Second-Generation p-values
#'
#' @description This function computes the false discovery risk for a second-generation \emph{p}-value of 0, or the false confirmation risk for a second-generation \emph{p}-value of 1.
#'
#' @param sgpvalue The observed second-generation \emph{p}-value
#' @param null.lo The lower bound of the indifference zone (null interval) upon which the second-generation \emph{p}-value was based
#' @param null.hi The upper bound for the indifference zone (null interval) upon which the second-generation \emph{p}-value was based
#' @param pt.est Point estimate of the unknown parameter
#' @param std.err Standard error of the point estimate
#' @param interval.type Class of interval estimate used. This determines the functional form of the power function. Options are \code{confidence} for a (1-$\alpha$)100\% confidence interval and \code{likelihood} for a 1/k likelihood support interval (\code{credible} not yet supported)
#' @param interval.level Level of interval estimate. If \code{interval.type} is \code{confidence}, the level is $\alpha$. If \code{interval.type} is \code{likelihood}, the level is 1/k.
#' @param null.weights Probability distribution for the null parameter space. Options are currently \code{point} and \code{uniform}.
#' @param null.space Support of the null probability distribution. If \code{null.weights} is \code{point}, then \code{null.space} is a scalar. If \code{null.weights} is \code{uniform}, then \code{null.space} is a vector of length two.
#' @param alt.weights Probability distribution for the alternative parameter space. Options are currently \code{point} and \code{uniform}.
#' @param alt.space Support for the alternative probability distribution. If \code{alt.weights} is \code{point}, then \code{alt.space} is a scalar. If \code{alt.weights} is \code{uniform}, then \code{alt.space} is a vector of length two.
#' @param pi0 Prior probability of the null hypothesis.
#'
#' @details
#' [Details]
#'
#' When possible, one should compute the second-generation \emph{p}-value and FDR/FCR on a scale that is symmetric about the null hypothesis. For example, if the parameter of interest is an odds ratio, inputs \code{pt.est}, \code{std.err}, \code{null.lo},  \code{null.hi}, \code{null.space}, and \code{alt.space} are typically on the log scale.
#'
#' @return Numeric scalar representing the False discovery risk (FDR) or false confirmation risk (FCR) for the observed second-generation \emph{p}-value. If \code{sgpvalue} = 0, the function returns false discovery risk (FDR). If \code{sgpvalue} = 1, the function returns false confirmation risk (FCR).
#' @seealso \code{\link{sgpvalue}}
#' @keywords
#' @export
#' @examples
#' # false discovery risk with 95% confidence level
#' fdrisk(sgpvalue = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  pt.est = 2.5,  std.err = 0.8,  null.weights = 'uniform',  null.space = c(log(1/1.1), log(1.1)),  alt.weights = 'uniform',  alt.space = 2.5 + c(-1,1)*qnorm(1-0.05/2)*0.8,  interval.type = 'confidence',  interval.level = 0.05)
#' [1] 0.04879258
#'
#' # false discovery risk with 1/8 likelihood support level
#' fdrisk(sgpvalue = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  pt.est = 2,  std.err = 0.8,  null.weights = 'point',  null.space = 0,  alt.weights = 'uniform',  alt.space = 2.5 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.04119162
#'
#' # false discovery risk with LSI and wider null hypothesis
#' fdrisk(sgpvalue = 0,  null.lo = log(1/1.5), null.hi = log(1.5),  pt.est = 2.5,  std.err = 0.8,  null.weights = 'point',  null.space = 0,  alt.weights = 'uniform',  alt.space = 2.5 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.01688343
#'
#' # false confirmation risk example
#' fdrisk(sgpvalue = 1,  null.lo = log(1/1.5), null.hi = log(1.5),  pt.est = 0.01,  std.err = 0.15,  null.weights = 'uniform',  null.space = 0.01 + c(-1,1)*qnorm(1-0.041/2)*0.15,  alt.weights = 'uniform',  alt.space = c(log(1.5), 1.25*log(1.5)),  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.03059522
#'
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Welty VF, Dupont WD, Greevy RA Jr. (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. In press. https://doi.org/10.1080/00031305.2018.1537893
#'



fdrisk <- function(sgpvalue, null.lo, null.hi, pt.est, std.err,
                   interval.type = 'likelihood', interval.level = 1/8, pi0 = 0.5,
                   null.weights, null.space, alt.weights, alt.space) {


  ### NOTE: need the indifference zone bounds for the power function
  ### NOTE: need to specify that it is a confidence interval for the power function

  if(0 < sgpvalue & sgpvalue < 1) stop('sgpvalue must take a value of 0 or 1 to use fdrisk')

  Fdr = NULL
  Fcr = NULL
  pi1 = 1 - pi0

  if(interval.type == 'confidence') {
    Z = qnorm(1-interval.level/2)
    est.lo = pt.est - Z*std.err
    est.hi = pt.est + Z*std.err
  }

  if(interval.type == 'likelihood') {
    Z = qnorm(1-2*pnorm(-sqrt(2*log(1/interval.level)))/2)
    est.lo = pt.est - Z*std.err
    est.hi = pt.est + Z*std.err
  }

  if(interval.type != 'confidence' & interval.type != 'likelihood') {
    stop('false discovery risk calculations for interval estimates other than confidence or likelihood support intervals are currently not supported')
  }

  # if sgpv != 0 or sgpv != 1
  if( (sgpvalue == 0 & ((est.lo > null.lo & est.lo < null.hi) | (est.hi > null.lo & est.hi < null.hi))) ) {
    stop('indifference zone must not intersect with the observed interval when SGPV = 0')
  }

  if(interval.type == 'confidence' | interval.type == 'likelihood') {

      ## `power0` = P(SGPV = 0 | theta = x) ##
      power0 = function(x) {
        pow.i = pnorm(null.lo/std.err - x/std.err - Z)
        pow.ii = pnorm(-null.hi/std.err + x/std.err - Z)
        return(pow.i + pow.ii)
      }
      ## `power1` = P(SGPV = 1 | theta = x) ##
      power1 = function(x) {
        pow.i = pnorm(null.hi/std.err - x/std.err - Z)
        pow.ii = pnorm(null.lo/std.err - x/std.err + Z)
        return(pow.i - pow.ii)
      }

  }

  ## `power` = P(SGPV = sgpvalue | theta = x)
  if(sgpvalue == 0) {power = function(x) {power0(x)}}
  if(sgpvalue == 1) {power = function(x) {power1(x)}}

  ## calculate P.sgpv.H0
  ## `P.sgpv.H0` = P(SGPV = sgpvalue | H0 )
  P.sgpv.H0 = NULL

  # If the indifference zone was a point, don't allow alternate specification of null probability distribution, just calculate type I at the point
  if(null.lo == null.hi)  {

    P.sgpv.H0 = power(x = null.lo)
    if(any(null.lo != null.space)) {warning(paste0('for a point indifference zone, specification of a different `null.space` not permitted; `null.space` set to be ',round(null.lo, 2),'.'))}

  }

  # Indifference zone not a point
  if(null.lo != null.hi) {

    # point
    if(null.weights == "point")  {

      if(length(null.space)!=1) stop('null space must be a vector of length 1 when using a point null probability distribution')

      P.sgpv.H0 = power(x = null.space)

      }

    # uniform
    if(null.weights == "uniform") {

      if(length(null.space)<2)  stop('null space must not be a point to use uniform averaging')

      if(length(null.space)==2) {

        if(max(null.space)>null.hi|min(null.space)<null.lo) {
          warning('null space must be inside originally specified indifference zone; bounds have been truncated')
          if(max(null.space)>null.hi) null.space[which.max(null.space)]=null.hi
          if(min(null.space)<null.lo) null.space[which.min(null.space)]=null.lo
        }

        P.sgpv.H0 = 1/(max(null.space) - min(null.space)) * integrate(f=power, lower = min(null.space), upper = max(null.space))$value
      }

    }

    # generalized beta
    if(null.weights == "GBeta") {
      warning('placeholder for future implementation of Generalized Beta null probability distribution')
    }

    # truncated normal
    if(null.weights == "TruncNormal") {
      warning('placeholder for future implementation of truncated Normal null probability distribution')
    }

  }

  ## calculate P.sgpv.H1
  ## `P.sgpv.H1` = P(SGPV = sgpvalue | H1 )
  P.sgpv.H1 = NULL

  # point
  if(alt.weights == "point")  {

    if(length(alt.space)!=1) stop('alternative space must be a vector of length 1 when using a point alternative probability distribution')

    if(length(alt.space)==1) {

      if( ((alt.space>=null.lo)&(alt.space<=null.hi)) ) stop('alternative space must be outside of the originally specified indifference zone')

      P.sgpv.H1 = power(x = alt.space)

    }

  }

  # uniform
  if(alt.weights == "uniform") {

    if(length(alt.space)<2)   stop('alternative space must not be a point to use uniform averaging')

    if(length(alt.space)==2)  {

      if(all(alt.space > null.lo) & all(alt.space < null.hi)) stop("alternative space may not be contained inside indifference zone; `null.space` and `alt.space` may be flipped")
      if(any(alt.space > null.lo & alt.space < null.hi)) stop("alternative space may not intersect indifference zone")

      P.sgpv.H1 = 1/(max(alt.space) - min(alt.space)) * integrate(f=power, lower = min(alt.space), upper = max(alt.space))$value

    }

  }

  # generalized beta
  if(alt.weights == "GBeta") {
    warning('placeholder for future implementation of Generalized Beta alternative probability distribution')
  }

  # truncated normal
  if(alt.weights == "TruncNormal") {
    warning('placeholder for future implementation of truncated Normal alternative probability distribution')
  }

  ### Calculate FDR or FCR
  if(sgpvalue == 0) Fdr = (1 + P.sgpv.H1 / P.sgpv.H0 * pi1 / pi0 ) ^ (-1)
  if(sgpvalue == 1) Fcr = (1 + P.sgpv.H0 / P.sgpv.H1 * pi0 / pi1 ) ^ (-1)

  return(c(Fdr, Fcr))   ## if FDR, Fcr is null, so c(Fdr, Fcr)=Fdr, and if FCR, Fdr is null, so c(Fdr, Fcr)=Fcr

}



