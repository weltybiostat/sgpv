################################################################
##	Purpose: 	Plot power curve from
##            Second-Generation p-values approach
##
##	Function:	plotsgpower
##	Version:	1.1
##
##	Author:		Rebecca T. Irlmeier, Valerie F. Welty and Jeffrey D. Blume
##	Date:		November 23, 2020
################################################################
#
#' Plot power curves for Second-Generation p-Values
#'
#' @description This function calculates power and type I error values from significance testing based on second-generation p-values as the inferential metric and plots the power curve to visualize the operating charateristics of the inferential procedure.
#'
#'
#' @param null.lo A scalar representing the lower bound of the null interval hypothesis (indifference zone) upon which the second-generation \emph{p}-value is based.
#' @param null.hi A scalar representing the upper bound of the null interval hypothesis (indifference zone) upon which the second-generation \emph{p}-value is based.
#' @param std.err Standard error for the distribution of the estimator for the parameter of interest. Note that this is the standard deviation for the estimator, not the standard deviation parameter for the data itself. This will be a function of the sample size(s).
#' @param alt Optional scalar or vector of alternative value(s) for the parameter of interest. Default is \code{NA}. If provided, a blue dotted line (or one at each point) will be plotted and the power will be printed.
#' @param x.lim Optional numeric vector of length two giving the lower and upper bounds of the x-axis for the power curve. Default is \code{NA}, where the x-axis range will be optimized to fit the entirety of the power curve (which is dependent upon the width of the null zone and the standard error of the estimator).
#' @param interval.type Class of interval estimate used for calculating the SGPV. Options are \code{"confidence"} for a \eqn{(1-\alpha)100}\% confidence interval and \code{"likelihood"} for a \eqn{1/k} likelihood support interval (\code{credible} not yet supported).
#' @param interval.level Level of interval estimate. If \code{interval.type = "confidence"} is used, the level is \eqn{\alpha}. If \code{interval.type = "likelihood"} is used, the level is \eqn{1/k} (not \eqn{k}).
#' @param plot.option Used to specify the type of plot desired. If \code{plot.option = 1}, the classical power curve and its corresponding SGPV power curve are shown. If \code{plot.option = 2}, the three power curves provided by \code{sgpower} are shown. Default is \code{plot.option = 1}.
#' @param title.lab Title text.
#' @param x.lab x-axis label.
#' @param y.lab y-axis label.
#' @param legend.on Toggle for plotting the legend. Default is \code{TRUE}.
#' @param null.col Coloring of shading for the null interval hypothesis (indifference zone) region. Default is Hawkes Blue: \code{null.col = rgb(208, 216, 232, maxColorValue = 255)}.
#' @param pow.col  Vector of length three specifying the colors for the the three power curves given when \code{plot.option = 2}. The first color option corresponds to the \eqn{Pr(SGPV = 0 | \theta)} line, the second color option corresponds to the \eqn{Pr(0 < SGPV < 1 | \theta)} line, and the third color option corresponds to the \eqn{Pr(SGPV = 1 | \theta)} line. Default is \code{pow.col = c("cornflowerblue", "firebrick3", "green4")}.
#' @param pow.lty Vector of length three specifying the line types (\code{lty}) for the three power curves given when \code{plot.option = 2}. The first line type option corresponds to the \eqn{Pr(SGPV = 0 | \theta)} line, the second line type option corresponds to the \eqn{Pr(0 < SGPV < 1 | \theta)} line, and the third line type option corresponds to the \eqn{Pr(SGPV = 1 | \theta)} line. Default is \code{pow.lty = c(1,1,1)} for solid lines.
#' @param null.pt Optional numeric scalar representing a point null hypothesis. Default is \code{NA}. If a value is given, it will be plotted as a black dashed line and the type I error at that point will be printed.
#' @param acc Optional parameter specifying the resolution of the x-axis. Default is \code{acc = 100} for plotting the power curve as a sequence of 100 (x, y) points.
#'
#' @seealso \code{\link{fdrisk}, \link{sgpvalue}, \link{plotsgpv}}
#' @export
#' @examples
#'
#' sigma = 5
#' n = 20
#'
#' plotsgpower(alt = NA, null.lo = -1, null.hi = 1,
#'             std.err = sigma/sqrt(n), x.lim = c(-8,8),
#'            interval.type = 'confidence', interval.level = 0.05,
#'            plot.option = 2, null.pt = 0)
#'
#' plotsgpower(alt = c(-4,2),
#'             null.lo = -1, null.hi = 1, std.err = sigma/sqrt(n),
#'             x.lim = NA, interval.type = 'confidence',
#'             interval.level = 0.05, plot.option = 2)
#'
#' plotsgpower(alt = NA, null.lo = -1, null.hi = 1,
#'             std.err = sigma/sqrt(n), x.lim = NA,
#'             interval.type = 'confidence', interval.level = 0.05,
#'             plot.option = 1, null.pt = NA)
#'
#' plotsgpower(alt = c(-4,2), null.lo = -1, null.hi = 1,
#'             std.err = 1, x.lim = NA, interval.type = 'likelihood',
#'             interval.level = 0.05, plot.option = 1, null.pt = 0)
#'
#'
#' @references
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. 73:sup1, 157-167, DOI: https://doi.org/10.1080/00031305.2018.1537893
#'
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' @importFrom graphics lines clip
#'

plotsgpower <-  function (null.lo, null.hi, std.err,
                          alt = NA, x.lim = NA,
                          interval.type, interval.level = 0.05, plot.option = 1,
                          null.col = rgb(208, 216, 232, maxColorValue = 255),
                          pow.col=c("cornflowerblue","firebrick3","green4"),
                          pow.lty=c(1,1,1),
                          title.lab = "", x.lab = "Parameter", y.lab = "Probability", legend.on = TRUE,
                          null.pt = NA, acc = 100) {

  #### Errors
  if (length(x.lim)!=2 & !is.na(x.lim[1])) {
    stop('x.lim needs a lower and upper bound')
  }

  #### Get range of truths
  if (is.na(x.lim[1])) {
    x.lim <- c(null.lo-5*std.err, null.hi+5*std.err)
    x.range <- seq(x.lim[1], x.lim[2], by=0.1/acc)
  }
  if (!is.na(x.lim[1])) {
    x.lim <- sort(x.lim)
    x.range <- seq(x.lim[1], x.lim[2], by=0.1/acc)
  }

  #### Null point
  if (is.na(null.pt)) {null.pt.plot <- mean(c(null.lo, null.hi))}
  if (!is.na(null.pt)) {null.pt.plot <- null.pt}

  ### Compute classical power
  if (interval.type == "confidence") {Z = qnorm(1 - interval.level/2)}
  if (interval.type == "likelihood") {Z = qnorm(1 - 2 * pnorm(-sqrt(2 * log(1/interval.level)))/2)}
  power.classic <- pnorm(null.pt.plot/std.err - x.range/std.err - Z) + pnorm(-null.pt.plot/std.err + x.range/std.err - Z)

  ### Compute SGPV power
  power.alt = sgpower(true=x.range, null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.alt
  power.inc = sgpower(true=x.range, null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.inc
  power.null = sgpower(true=x.range, null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.null

  #### Set plot limits
  x.limits <- c(x.lim[1], x.lim[2])
  y.limits <- c(0,1)

  ### Plot
  plot(1, null.pt.plot, ylim=y.limits, xlim=x.limits, type="n",
       ylab=y.lab, xlab=x.lab, main=title.lab)
  if (!is.na(alt[1])) {
    axis(side = 1, at = alt)
  }
  #axis(side = 1, at = c(ceiling((x.limits[1])), ceiling((x.limits[1]+null.pt.plot)/2), null.pt.plot, floor((null.pt.plot+x.limits[2])/2), floor(x.limits[2]), alt))
  #axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1))
  rect(null.lo, 0, null.hi, 1, col = null.col, border = NA)

  # if (plot.axis[1]=="TRUE") {axis(side=1)}
  # if (plot.axis[2]=="TRUE") {axis(side=2)}

  if (plot.option == 1) {
    # power.alt
    lines(x.range, power.alt, type="l", col=pow.col[1], lty=pow.lty[1], lwd=1.5)
    # power.classic
    lines(x.range, power.classic, type="l", col=pow.col[2], lty=pow.lty[1], lwd=1.5)
  }

  if (plot.option == 2) {
    # power.alt
    lines(x.range, power.alt, type="l", col=pow.col[1], lty=pow.lty[1], lwd=1.5)
    # power.inc
    lines(x.range, power.inc, col=pow.col[2], lty=pow.lty[2], lwd=1.5)
    # power.null
    lines(x.range, power.null, col=pow.col[3], lty=pow.lty[3], lwd=1.5)
  }

  ### Add legend
  if (legend.on == TRUE) {
    if (plot.option == 1) {
      if (!is.na(alt[1])) {
        legend("right", c("Interval Null", "SGPV Power", "Classic Power", "Alternative"),
               col = c(null.col, pow.col[1], pow.col[2], "black"), lwd = c(8, 1.5, 1.5, 1), lty = c(1, pow.lty[1], pow.lty[2], 3),
               bty = "n", cex=0.7, pch=c())
      }
      if (is.na(alt[1])) {
        legend("right", c("Interval Null", "SGPV Power", "Classic Power"),
               col = c(null.col, pow.col[1], pow.col[2]), lwd = c(8, 1.5, 1.5), lty = c(1, pow.lty[1], pow.lty[2]),
               bty = "n", cex=0.7)
      }
    }
    if (plot.option == 2) {
      if (!is.na(alt[1])) {
        legend("right", c("Interval Null", expression("Pr(p"[delta] *" = 1)"), expression("Pr(0 < p"[delta]*" < 1)"), expression("Pr(p"[delta] *" = 1)"), "Alternative"),
               col = c(null.col, pow.col[1], pow.col[2], pow.col[3], "black"), lwd = c(8, 1.5, 1.5, 1.5, 1), lty = c(1, pow.lty[1], pow.lty[2], pow.lty[3], 3),
               bty = "n", cex=0.7)
      }
      if (is.na(alt[1])) {
        legend("right", c("Interval Null", expression("Pr(p"[delta] *" = 1)"), expression("Pr(0 < p"[delta]*" < 1)"), expression("Pr(p"[delta] *" = 1)")),
               col = c(null.col, pow.col[1], pow.col[2], pow.col[3]), lwd = c(8, 1.5, 1.5, 1.5), lty = c(1, pow.lty[1], pow.lty[2], pow.lty[3]),
               bty = "n", cex=0.7)
      }
    }
  }

  ### Add null.pt
  if (!is.na(null.pt)) {
    clip(x1=x.lim[1], x2=x.lim[2], y1=0, y2=1)
    abline(v=null.pt.plot, lty=2, lwd=0.5)
  }

  if (!is.na(null.pt) & !is.na(alt[1])) {
    pow.pt <- c(null.pt, alt)
  }
  if (!is.na(null.pt) & is.na(alt[1])) {
    pow.pt <- null.pt
  }
  if (is.na(null.pt) & !is.na(alt[1])) {
    pow.pt <- alt
  }
  if (is.na(null.pt) & is.na(alt[1])) {
    pow.pt <- NA
  }
  pow.pt <- sort(pow.pt)

  ### Add vertical line for alternative of interest
  if (!is.na(pow.pt[1])) {
    pow.pt.df <- matrix(NA, nrow=length(pow.pt), ncol=4)
    for (i in 1:length(pow.pt)) {
      abline(v=pow.pt[i], lty=3, col="black")
      pow.pt.df[i,1] <- pnorm(null.pt.plot/std.err - pow.pt[i]/std.err - Z) + pnorm(-null.pt.plot/std.err + pow.pt[i]/std.err - Z)
      pow.pt.df[i,2] <- sgpower(true=pow.pt[i], null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.alt
      pow.pt.df[i,3] <- sgpower(true=pow.pt[i], null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.inc
      pow.pt.df[i,4] <- sgpower(true=pow.pt[i], null.lo=null.lo, null.hi=null.hi, std.err=std.err, interval.type=interval.type, interval.level=interval.level)$power.null
    }
    #pow.pt.df <- round(as.data.frame(pow.pt.df), 4)
    pow.pt.df <- round(pow.pt.df, 4)
    rownames(pow.pt.df) <- ifelse(pow.pt < 0, paste0("At ", pow.pt, ":"), paste0("At  ", pow.pt, ":"))
    colnames(pow.pt.df) <- c("   Classic Power", "     SGPV Power*", "     Pr(0 < SGPV < 1)", "     Pr(SGPV = 1)")
    print(pow.pt.df)
    cat("--\n")
    cat("*SGPV Power: Pr(SGPV = 0)")
  }
}



## Examples

# sigma = 5
# n = 20
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=c(-8,8), interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=0)
# plotsgpower(alt=c(-4,2), null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=NA, interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=NA)
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=NA, interval.type='confidence', 'interval.level'=0.05, plot.option=1, null.pt=NA, y.lab="Power")
# plotsgpower(alt=c(-4,2), null.lo=-1, null.hi=1, std.err=1, x.lim=NA, interval.type='likelihood', 'interval.level'=0.05, plot.option=1, null.pt=0)
#
# a=plotsgpower(alt=c(-8,2,4), null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=c(-10,10), interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=NA)
# plotsgpower(alt=c(-1,2,4), null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=NA, interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=NA)
# plotsgpower(alt=2, null.lo=-1, null.hi=1, std.err=sigma/sqrt(n), x.lim=c(-10,10), interval.type='confidence', 'interval.level'=0.05, plot.option=1, null.pt=0)
#
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=0.5, x.lim=c(-10,10), interval.type='likelihood', 'interval.level'=0.05, plot.option=2, null.pt=0)
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=1, x.lim=NA, interval.type='likelihood', 'interval.level'=0.05, plot.option=1, null.pt=0)
#
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=0.3, x.lim=NA, interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=0)
# plotsgpower(alt=NA, null.lo=-1, null.hi=1, std.err=0.3, x.lim=c(-10,10), interval.type='confidence', 'interval.level'=0.05, plot.option=1, null.pt=NA)
#
# plotsgpower(alt=NA, null.lo=5, null.hi=7, std.err=sigma/sqrt(n), x.lim=c(1,11), interval.type='confidence', 'interval.level'=0.05, plot.option=2, null.pt=NA)

