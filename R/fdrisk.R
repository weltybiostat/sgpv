################################################################
##	Purpose: 	Compute False Discovery or Confirmation risk
##            for Second-Generation p-values approach
##
##	Function:	fdrisk
##	Version:	1.2
##
##	Author:		Valerie F. Welty and Jeffrey D. Blume
##	Date:		December 6, 2018
################################################################
#
#' False Discovery Risk for Second-Generation p-values
#'
#' @description This function computes the false discovery risk (sometimes called the "empirical bayes FDR") for a second-generation \emph{p}-value of 0, or the false confirmation risk for a second-generation \emph{p}-value of 1.
#'
#' @param sgpval The observed second-generation \emph{p}-value. Default is \eqn{0}, which gives the false discovery risk.
#' @param null.lo The lower bound of the indifference zone (null interval) upon which the second-generation \emph{p}-value was based
#' @param null.hi The upper bound for the indifference zone (null interval) upon which the second-generation \emph{p}-value was based
#' @param std.err Standard error of the point estimate
#' @param interval.type Class of interval estimate used. This determines the functional form of the power function. Options are \code{confidence} for a \eqn{(1-\alpha)100}\% confidence interval and \code{likelihood} for a \eqn{1/k} likelihood support interval (\code{credible} not yet supported).
#' @param interval.level Level of interval estimate. If \code{interval.type} is \code{confidence}, the level is \eqn{\alpha}. If \code{interval.type} is \code{likelihood}, the level is \eqn{1/k} (not \eqn{k}).
#' @param null.weights Probability distribution for the null parameter space. Options are currently \code{Point}, \code{Uniform}, and \code{TruncNormal}.
#' @param null.space Support of the null probability distribution. If \code{null.weights} is \code{Point}, then \code{null.space} is a scalar. If \code{null.weights} is \code{Uniform}, then \code{null.space} is a vector of length two.
#' @param alt.weights Probability distribution for the alternative parameter space. Options are currently \code{Point}, \code{Uniform}, and \code{TruncNormal}.
#' @param alt.space Support for the alternative probability distribution. If \code{alt.weights} is \code{Point}, then \code{alt.space} is a scalar. If \code{alt.weights} is \code{Uniform}, then \code{alt.space} is a vector of length two.
#' @param pi0 Prior probability of the null hypothesis. Default is \eqn{0.5}.
#'
#' @details
#'
#' When possible, one should compute the second-generation \emph{p}-value and FDR/FCR on a scale that is symmetric about the null hypothesis. For example, if the parameter of interest is an odds ratio, inputs \code{pt.est}, \code{std.err}, \code{null.lo},  \code{null.hi}, \code{null.space}, and \code{alt.space} are typically on the log scale.
#'
#' If \code{TruncNormal} is used for \code{null.weights}, then the distribution used is a truncated Normal distribution with mean equal to the midpoint of \code{null.space}, and standard deviation equal to \code{std.err}, truncated to the support of \code{null.space}. If \code{TruncNormal} is used for \code{alt.weights}, then the distribution used is a truncated Normal distribution with mean equal to the midpoint of \code{alt.space}, and standard deviation equal to \code{std.err}, truncated to the support of \code{alt.space}. Further customization of these parameters for the truncated Normal are currently not possible, although they may be implemented in future versions.
#'
#' @return Numeric scalar representing the False discovery risk (FDR) or false confirmation risk (FCR) for the observed second-generation \emph{p}-value. If \code{sgpval} = \eqn{0}, the function returns false discovery risk (FDR). If \code{sgpval} = \eqn{1}, the function returns false confirmation risk (FCR).
#' @seealso \code{\link{sgpvalue}, \link{sgpower}}
#' @keywords
#' @export
#' @examples
#'
#' # false discovery risk with 95% confidence level
#' fdrisk(sgpval = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  std.err = 0.8,  null.weights = 'Uniform',  null.space = c(log(1/1.1), log(1.1)),  alt.weights = 'Uniform',  alt.space = 2 + c(-1,1)*qnorm(1-0.05/2)*0.8,  interval.type = 'confidence',  interval.level = 0.05)
#' [1] 0.059499
#'
#' # false discovery risk with 1/8 likelihood support level
#' fdrisk(sgpval = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  std.err = 0.8,  null.weights = 'Point',  null.space = 0,  alt.weights = 'Uniform',  alt.space = 2 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.050555
#'
#' ## with truncated normal weighting distribution
#' fdrisk(sgpval = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  std.err = 0.8,  null.weights = 'Point',  null.space = 0,  alt.weights = 'TruncNormal',  alt.space = 2 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.049026
#'
#' # false discovery risk with LSI and wider null hypothesis
#' fdrisk(sgpval = 0,  null.lo = log(1/1.5), null.hi = log(1.5),  std.err = 0.8,  null.weights = 'Point',  null.space = 0,  alt.weights = 'Uniform',  alt.space = 2.5 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.016883
#'
#' # false confirmation risk example
#' fdrisk(sgpval = 1,  null.lo = log(1/1.5), null.hi = log(1.5),  std.err = 0.15,  null.weights = 'Uniform',  null.space = 0.01 + c(-1,1)*qnorm(1-0.041/2)*0.15,  alt.weights = 'Uniform',  alt.space = c(log(1.5), 1.25*log(1.5)),  interval.type = 'likelihood',  interval.level = 1/8)
#' [1] 0.030595
#'
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. In press. https://doi.org/10.1080/00031305.2018.1537893
#'


fdrisk <- function(sgpval=0, null.lo, null.hi, std.err,
                   interval.type, interval.level, pi0 = 0.5,
                   null.weights, null.space, alt.weights, alt.space) {

  ## Warnings
  if(!(sgpval %in% c(0, 1))) stop('sgpval must take a value of 0 or 1 to use fdrisk')

  if (!(interval.type %in% c('confidence', 'likelihood'))) stop("Parameter `interval.type` must be one of the following: \n  * confidence \n  * likelihood \n  (credible/bayesian not currently supported for fdrisk)")

  ## Relevant quantities
  Fdr = NULL  # FDR = (1 + P(SGPV=0 | H1 ) / P(SGPV=0 | H0 ) *  P(H1) / P(H0) ) ^ (-1)
  Fcr = NULL  # FCR = (1 + P(SGPV=1 | H0 ) / P(SGPV=1 | H1 ) *  P(H0) / P(H1) ) ^ (-1)

  P.sgpv.H0 = NULL  # `P.sgpv.H0` = P(SGPV=0 | H0 )
  P.sgpv.H1 = NULL  # `P.sgpv.H1` = P(SGPV=1 | H1 )

  ## Power functions
  if(sgpval==0) {power = function(x) {sgpower(true = x, null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.alt}}
  if(sgpval==1) {power = function(x) {sgpower(true = x, null.lo, null.hi, std.err, interval.type, interval.level)$power.null}}
    # power(x) = P(SGPV=`sgpval` | theta = x)
      # if `sgpval`=0, power (probability) to get SGPV=0 (Blume et al. (2018) eq.(S4) for CI/LSI w/ symmetric null)
      # if `sgpval`=1, power (probability) to get SGPV=1 (Blume et al. (2018) eq.(S7) for CI/LSI w/ symmetric null)

  ### calculate P.sgpv.H0

    # point null
    if(null.lo == null.hi)  {
      if(any(null.lo != null.space)) {warning(paste0('for a point indifference zone, specification of a different `null.space` not permitted; `null.space` set to be ',round(null.lo, 2),'.'))}
      P.sgpv.H0 = power(x = null.lo)
    }

    # interval null
    if(null.lo != null.hi) {

      # P.sgpv.H0 @ point (=type I error at null.space)
      if(null.weights == "Point")  {
        if(length(null.space)!=1) stop('null space must be a vector of length 1 when using a point null probability distribution')
        P.sgpv.H0 = power(x = null.space)
      }

      # P.sgpv.H0 averaged: check `null.space` input
      if(null.weights %in% c('Uniform', 'GBeta', 'TruncNormal')) {
        if(length(null.space)<2)  stop('null space must not be a point to use averaging methods')
        if(length(null.space)==2) {
          # truncate bounds to edge of null if null.spcae falls outside indifference zone
          if(max(null.space)>null.hi|min(null.space)<null.lo) {
            warning('null space must be inside originally specified null hypothesis; at least one null space bound has been truncated')
            if(max(null.space)>null.hi) null.space[which.max(null.space)]=null.hi
            if(min(null.space)<null.lo) null.space[which.min(null.space)]=null.lo
          }
        }
      }

      # P.sgpv.H0 averaged uniformly
      if(null.weights == 'Uniform') {
        P.sgpv.H0 = 1/(max(null.space) - min(null.space)) * integrate(f=power, lower = min(null.space), upper = max(null.space))$value
      }

      # P.sgpv.H0 averaged using generalized beta as weighting distribution function
      if(null.weights == "GBeta") {
        warning('placeholder for future implementation of Generalized Beta null probability distribution')

        P.sgpv.H0 = NULL
      }

      # P.sgpv.H0 averaged using truncated normal as weighting distribution function
      if(null.weights == "TruncNormal") {

        # default: mean of Normal distr at midpoint of null.space
        truncNorm.mu = mean(null.space)
        # default: std. dev of Normal distr same as assumed for estimator
        truncNorm.sd = std.err

        integrand = function(x) {
          power(x) * ( dnorm(x, truncNorm.mu, truncNorm.sd) * (pnorm(max(null.space), truncNorm.mu, truncNorm.sd) - pnorm(min(null.space), truncNorm.mu, truncNorm.sd))^(-1) )
        }

        P.sgpv.H0 = integrate(f=integrand, lower = min(null.space), upper = max(null.space))$value

      }

    }

  ### calculate P.sgpv.H1

    # P.sgpv.H1 @ point
    if(alt.weights == "Point")  {
      if(length(alt.space)!=1) stop('alt space must be a vector of length 1 when using a point alternative probability distribution')
      if( ((alt.space>=null.lo)&(alt.space<=null.hi)) ) stop('alternative space must be outside of the originally specified indifference zone')

      P.sgpv.H1 = power(x = alt.space)
    }

    # P.sgpv.H1 averaged: check `alt.space` input
    if(alt.weights %in% c('Uniform', 'GBeta', 'TruncNormal')) {
      if(length(alt.space)<2)   stop('alt space must not be a point to use averaging methods')
      if(length(alt.space)==2)  {
        if(all(alt.space > null.lo) & all(alt.space < null.hi)) stop("alternative space can not be contained inside indifference zone; `null.space` and `alt.space` might be flipped")
        if(any(alt.space > null.lo & alt.space < null.hi)) stop("alternative space can not intersect indifference zone")
      }
    }

    # P.sgpv.H1 averaged uniformly
    if(alt.weights == 'Uniform') {
      P.sgpv.H1 = 1/(max(alt.space) - min(alt.space)) * integrate(f=power, lower = min(alt.space), upper = max(alt.space))$value
    }

    # P.sgpv.H1 averaged using generalized beta as weighting distribution function
    if(alt.weights == "GBeta") {
      warning('placeholder for future implementation of Generalized Beta null probability distribution')
      P.sgpv.H1 = NULL
    }

    # P.sgpv.H1 averaged using truncated normal as weighting distribution function
    if(alt.weights == "TruncNormal") {

      # default: mean of Normal distr at midpoint of alt.space
      truncNorm.mu = mean(alt.space)
      # default: std. dev of Normal distr same as assumed for estimator
      truncNorm.sd = std.err

      if(any(is.na(c(truncNorm.mu, truncNorm.sd)))) stop('error: `trunNorm.mu` and `truncNorm.sd` must be NULL or numeric; may not be NA')

      integrand = function(x) {
        power(x) * ( dnorm(x, truncNorm.mu, truncNorm.sd) * (pnorm(max(alt.space), truncNorm.mu, truncNorm.sd) - pnorm(min(alt.space), truncNorm.mu, truncNorm.sd))^(-1) )
      }

      P.sgpv.H1 = integrate(f=integrand, lower = min(alt.space), upper = max(alt.space))$value

    }

  ## Calculate FDR or FCR
  if(sgpval == 0) Fdr = (1 + P.sgpv.H1 / P.sgpv.H0 *  (1-pi0) / pi0     ) ^ (-1)
  if(sgpval == 1) Fcr = (1 + P.sgpv.H0 / P.sgpv.H1 *  pi0     / (1-pi0) ) ^ (-1)

  return(c(Fdr, Fcr))   ## if FDR, Fcr is null, so c(Fdr, Fcr)=Fdr, and if FCR, Fdr is null, so c(Fdr, Fcr)=Fcr

}



