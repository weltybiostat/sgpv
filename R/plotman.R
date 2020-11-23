################################################################
##	Purpose: 	Plot a "modified Manhattan-style plot" colored
##          according to Second-Generation p-value status
##
##	Function:	plotman
##	Version:	1.1
##
##	Author:		Rebecca T. Irlmeier, Valerie F. Welty and Jeffrey D. Blume
##	Date:		November 23, 2020
################################################################
#
#' Second-Generation p-Value Plotting
#'
#' @description This function displays a modified Manhattan-style plot colored according to second-generation p-value status. There are several variations of this plot that can be made depending upon user input for \code{type} as well as the \code{set.order} and \code{x.show} options. These plots allow the user to visualize the overall result of a large-scale analysis succintly and to visually assess the differences in the results using second-generation p-value techniques as opposed to classical p-value techniques.
#'
#'
#' @param est.lo A numeric vector of lower bounds of interval estimates. Must be of same length as \code{est.hi}.
#' @param est.hi A numeric vector of upper bounds of interval estimates.  Must be of same length as \code{est.lo}.
#' @param null.lo A scalar representing the lower bound of the null interval hypothesis (indifference zone). Value must be finite.
#' @param null.hi A scalar representing the upper bound of the null interval hypothesis (indifference zone). Value must be finite.
#' @param set.order A numeric vector giving the desired order along the x-axis. Alternatively, if \code{set.order} is set to \code{"sgpv"}, the second-generation \emph{p}-value ranking is used. The default option is \code{NA}, which uses the original input ordering.
#' @param x.show A numeric scalar representing the maximum ranking on the x-axis that is displayed. Default is to display all rankings.
#' @param type A string specifying the desired Manhattan-style plot to be graphed. This argument specifies the variable on the y-axis. If \code{type = "delta-gap"}, the delta-gaps are ranked. If \code{type = "p-value"}, the classic p-values are ranked. If \code{type = "comparison"}, the classic p-values are ranked by SGPV. Default is \code{type = "delta-gap"}.
#' @param p.values A numeric vector giving the classic \emph{p}-values. This is required when \code{type = "p-value"} or \code{type = "comparison"}, and is not required when \code{type = "delta-gap"}. The \code{p.values} input may be any desired transformation of the p-values. For example, if the desired transformation is \eqn{-log10(p-value)} as in a traditional Manhattan plot, the \eqn{-log10(p-values)} should be provided for \code{p.values}. The corresponding x or y axis label(s) should be updated to reflect any transformations.
#' @param ref.lines A numeric scalar or vector giving the points on the y-axis at which to add a horizontal reference line. For example, if \code{p.values} is set to \eqn{-log10(p-values)} and the type of plot selected shows the (transformed) p-values on the y-axis, possible locations for the reference lines could be at the \eqn{-log10(0.05)}, \eqn{-log10(Bonferroni)} and \eqn{-log10(FDR)} significance levels.
#' @param int.col Vector of length three specifing the colors of the points according to SGPV result. The first color option corresponds to the \eqn{SGPV = 0} results, the second color option corresponds to the \eqn{0 < SGPV < 1} results, and the third color option corresponds to the \eqn{SGPV = 1} results. Default is \code{int.col = c("cornflowerblue","firebrick3","darkslateblue")}.
#'
#' @param int.pch Plotting symbol for points. Default is \code{16} for small points.
#' @param int.cex Size of plotting symbol for points. Default is \code{0.4}.
#' @param null.pt An optional numeric scalar representing a point null hypothesis. Default is \code{NA}.
#' @param title.lab Title text.
#' @param x.lab A title for the x-axis. Default is the generic \code{"Position (by set.order)"}.
#' @param y.lab A title for the y-axis. Default is the generic \code{"Outcome label"}.
#' @param legend.on Toggle for plotting the legend. Default is \code{TRUE}.
#'
#'
#'
#' @details
#'
#' Use \code{set.order} to provide the classical p-value ranking. For example, if \code{pvalue.vector} is a vector of classical p-values, then set \code{set.order=order(pvalue.vector)} to sort the x-axis according to p-value rank.
#'
#' Use \code{type} and \code{p.values} to provide the \eqn{-log10(p-values)} for the y-axis. For example, if \code{pvalue.vector} is a vector of classical p-values, then set \code{type="p-value"} (or \code{type="comparison"}) and \code{p.values=-log10(pvalue.vector)} to set the y-axis. Then, set the y-axis title to something like \code{y.lab="-log10(p)"}.
#'
#'
#' @seealso \code{\link{sgpvalue}, \link{plotsgpv}, \link{sgpower}, \link{plotsgpower}}
#' @export
#' @examples
#'
#'
#' #  Use leukstats data
#' data(leukstats)
#'
#' # ID number on the x-axis, delta-gap on the y-axis, using an interval null hypothesis of (-0.3, 0.3) for the log mean difference in expression levels (fold change).
#' plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#'        null.lo=-0.3, null.hi=0.3,
#'        set.order=NA,
#'        type="delta-gap",
#'        ref.lines=NA,
#'        int.pch=16, int.cex=0.4,
#'        title.lab="Leukemia Example",
#'        y.lab="Delta-gap",
#'        x.lab="Position (ID)",
#'        legend.on=TRUE)
#'
#' # ID number on the x-axis, -log10(classical p-value) on the y-axis, using an interval null hypothesis of (-0.3, 0.3) for the log mean difference in expression levels (fold change).
#' plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#'        null.lo=-0.3, null.hi=0.3,
#'        set.order=NA,
#'        type="p-value",
#'        p.values=-log10(leukstats$p.value),
#'        ref.lines=-log10(0.05),
#'        int.pch=16, int.cex=0.4,
#'        title.lab="Leukemia Example",
#'        y.lab=expression("-log"[10]*"(p-value)"),
#'        x.lab="Position (ID)",
#'        legend.on=TRUE)
#'
#' # Second-generation p-value (SGPV) on the x-axis, -log10(classical p-value) on the y-axis, using an interval null hypothesis of (-0.3, 0.3) for the log mean difference in expression levels (fold change).
#' plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#'        null.lo=-0.3, null.hi=0.3,
#'        set.order="sgpv",
#'        type="comparison",
#'        p.values=-log10(leukstats$p.value),
#'        ref.lines=c(-log10(0.05), -log10(0.001)),
#'        int.pch=16, int.cex=0.4,
#'        title.lab="Leukemia Example",
#'        y.lab=expression("-log"[10]*"(p-value)"),
#'        x.lab="Second-generation p-value ranking",
#'        legend.on=TRUE)
#'
#' @references
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. 73:sup1, 157-167, DOI: https://doi.org/10.1080/00031305.2018.1537893
#'
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#'
#'
#' @importFrom grDevices rgb
#' @importFrom graphics abline axis legend plot points rect segments par mtext text


plotman <- function( est.lo, est.hi, null.lo, null.hi,
                     set.order=NA, x.show=NA, type="delta-gap",
                     p.values=NA, ref.lines=NA, null.pt=NA,
                     int.col = c("cornflowerblue", "firebrick3", "darkslateblue"),
                     int.pch=16, int.cex=0.4,
                     title.lab=NA, x.lab="Position (by set.order)", y.lab="Outcome label",
                     legend.on=TRUE ) {

  #### Errors
  if(length(null.lo)!=length(null.hi)) {
    stop('null.lo and null.hi of different lengths')
  }

  if(length(est.lo)!=length(est.hi)) {
    stop('est.lo and est.hi of different lengths')
  }

  if(length(null.lo)!=1|length(null.lo)!=1) {
    stop('null.lo and null.hi must be scalars')
  }

  if(type %in% c("delta-gap", "p-value", "comparison") == FALSE) {
    stop('type must be "delta-gap", "p-value" or "comparison"')
  }

  if(type %in% c("p-value", "comparison") & length(p.values)==1) {
    stop('input vector of p-values')
  }

  #### Compute SGPVs
  sgpv <- sgpvalue(est.lo=est.lo, est.hi=est.hi, null.lo=null.lo, null.hi=null.hi)
  sgpv.combo <- ifelse(sgpv$p.delta==0,-sgpv$delta.gap,sgpv$p.delta)

  #### y-axis
  if (type=="delta-gap") {
    y.var.plot <- sgpv$delta.gap
    y.var.plot <- ifelse(is.na(y.var.plot), 0, y.var.plot)
  }
  if (type %in% c("p-value", "comparison")) {
    y.var.plot <- p.values
  }

  #### Set plot limits
  x.max 	<- length(est.lo)
  x 		<- seq(1,x.max,1)

  if (is.na(x.show)) { x.show <- x.max}
  x.limits <- c(1,min(x.show,x.max,na.rm=TRUE))
  y.limits <- c(0,ceiling(max(y.var.plot, na.rm=TRUE)))

  #### Set order of x-axis (by input ordering)
  if (is.na(set.order[1]))  {
    set.order <- x
    set.order <- set.order[1:x.show]

    #### Plot for type = "delta-gap"
    if (type=="delta-gap") {
      par(mar = c(6, 5, 5, 5))
      plot(1, null.pt, ylim=y.limits, xlim=x.limits, type="n",
         ylab=y.lab, xlab=x.lab, main=title.lab, las=1)
      mtext("*", at=0, side=2, las=2, adj=1, line=0.65)
      # Footnote
      mtext("* Includes instances where the delta-gap is not defined", side=1, at=0, adj=-0.23, cex=0.8, line=4.5)
    }
    #### Plot for type = "p-value"
    if (type=="p-value") {
      par(mar = c(5, 5, 5, 5))
      plot(1, null.pt, ylim=y.limits, xlim=x.limits, type="n",
           ylab=y.lab, xlab=x.lab, main=title.lab, las=1)
    }
  }

  #### Set order of x-axis (by SGPV)
  if (set.order[1]=="sgpv") {
    set.order <- order(sgpv.combo)
    set.order <- set.order[1:x.show]

    # If including SGPV = 1
    if (x.show > length(sgpv$p.delta)-sum(sgpv$p.delta==1)) {
      x.ord <- sum(sgpv$p.delta!=1) # number ordered by p-value/delta-gap
      x.NA <- x.show-x.ord # number ordered by input ordering since delta-gap is NA
      par(mar = c(6, 5, 5, 5))
      plot(1, null.pt, ylim=y.limits, xlim=x.limits, type="n", xaxt="n",
           ylab=y.lab, xlab=x.lab, main=title.lab, las=1)
      # SGPV = 0 and 0 < SGPV < 1 x-axis labeling
      round.digits <- ifelse(x.ord > 2500, -3, ifelse(x.ord > 1000, -2, -1))
      axis(side=1, at=floor(seq(x.limits[1], x.ord, by=round(x.ord/4, round.digits)))-c(0,1,1), labels=FALSE)
      text(x=floor(seq(x.limits[1], x.ord, by=round(x.ord/4, round.digits))), y=-1.5, labels=floor(seq(x.limits[1], x.ord, by=round(x.ord/4, round.digits)))-c(0,1,1), srt=35, pos=1, xpd=TRUE)
      # SGPV = 1 x-axis labeling
      segments(x0=x.ord+1, y0=y.limits[1], x1=x.ord+1, y1=mean(y.limits), lty=2, lwd=1)
      rep.NA <- ifelse(x.NA/x.show < 0.1, 1, ifelse(x.NA/x.show < 0.2, 2, 3))
      axis(side=1, at=seq(x.ord+1, x.show, length.out=rep.NA), labels=FALSE)
      text(x=seq(x.ord+1, x.show, length.out=rep.NA), y=-1.5, labels=rep(x.ord+1,rep.NA), srt=35, pos=1, xpd=TRUE)
      # Footnote
      mtext("Delta-gap rankings used when SGPV = 0; when SGPV = 1, rankings are tied", side=1, at=1, adj=-0.05, cex=0.8, line=4.5)
    }

    # SGPV = 0 and 0 < SGPV < 1 (with no SGPV = 1) x-axis labeling
    if (x.show <= length(sgpv$p.delta)-sum(sgpv$p.delta==1)) {
      par(mar = c(6, 5, 5, 5))
      plot(1, null.pt, ylim=y.limits, xlim=x.limits,
           ylab=y.lab, xlab=x.lab, main=title.lab, las=1)
      # Footnote
      mtext("Delta-gap rankings used when SGPV = 0; when SGPV = 1, rankings are tied", side=1, at=1, adj=-0.05, cex=0.8, line=4.5)
    }
  }

  #### Subset intervals by SGPV value for coloring
  gap.marker <- 1*is.na(sgpv$delta.gap[set.order])
  x <- x[1:x.show]

  #### Plot points for type = "delta-gap"
  if (type=="delta-gap") {
    set.out  <- x[gap.marker[x]==0]
    set.in.both   <- x[(gap.marker[x]==1)]

    ## Intervals where 0<SGPV<=1
    points(x[set.in.both],y.var.plot[set.order][set.in.both],cex=int.cex,pch=int.pch,col=int.col[1])

    ## Intervals where SGPV==0
    points(x[set.out],y.var.plot[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])

    #### Legend for type = "delta-gap"
    if (legend.on==TRUE) {
      legend("topright", c(expression("p"[delta]*" = 0"),
                           expression("p"[delta]*" > 0")),
             lty = 1, col = c(int.col[2], int.col[1]),
             lwd = c(1.5, 1.5), bty = "n")
    }
  }

  #### Plot points for type = "p-value" and type = "comparison"
  if (type!="delta-gap") {
    set.out  <- x[gap.marker[x]==0]
    set.in <- x[(gap.marker[x] == 1) & (sgpv$p.delta[set.order] == 1)]
    set.both <- x[(gap.marker[x] == 1) & (sgpv$p.delta[set.order] < 1)]

    ## Intervals where SGPV==0
    points(x[set.out],y.var.plot[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])

    ## Intervals where 0<SGPV<1
    points(x[set.both],y.var.plot[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])

    ## Intervals where SGPV==1
    points(x[set.in],y.var.plot[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])

    ## Line for y-axis cutoffs
    if (!is.na(ref.lines[1])) {
      abline(h=ref.lines, lty=2)
    }

    #### Legend for "p-value"
    if (legend.on==TRUE & type=="p-value") {
      legend("topright", c(expression("p"[delta]*" = 0"),
                           expression("0 < p"[delta]*" < 1"),
                           expression("p"[delta] *" = 1")),
             lty = 1, col = c(int.col[2], int.col[1], int.col[3]),
             lwd = c(1.5, 1.5, 1.5), bty = "n")
    }

    #### Legend for type = "comparison" and showing SGPV = 1
    if (legend.on==TRUE & type=="comparison" & x.show > length(sgpv$p.delta)-sum(sgpv$p.delta==1)) {
      legend("topright", c(expression("p"[delta]*" = 0"),
                           expression("0 < p"[delta]*" < 1"),
                           expression("p"[delta] *" = 1")),
             lty = 1, col = c(int.col[2], int.col[1], int.col[3]),
             lwd = c(1.5, 1.5, 1.5), bty = "n")
    }

    #### Legend for type = "comparison" and not showing SGPV = 1
    if (legend.on==TRUE & type=="comparison" & x.show <= length(sgpv$p.delta)-sum(sgpv$p.delta==1)) {
      legend("topright", c(expression("p"[delta]*" = 0"),
                           expression("0 < p"[delta]*" < 1")),
             lty = 1, col = c(int.col[2], int.col[1]),
             lwd = c(1.5, 1.5), bty = "n")
    }
  }
}




# # Examples
# library(sgpv)
# data(leukstats)
#
# # x-axis ID, y-axis delta-gap
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#          null.lo=-0.3, null.hi=0.3,
#          set.order=NA,
#          type="delta-gap",
#          ref.lines=NA,
#          int.pch=16, int.cex=0.4,
#          title.lab="Leukemia Example",
#          y.lab="Delta-gap",
#          x.lab="Position (ID)",
#          legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order=NA,
#         x.show=3000,
#         type="delta-gap",
#         ref.lines=NA,
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab="Delta-gap",
#         x.lab="Position (ID)",
#         legend.on=TRUE)
#
# # x-axis ID, y-axis -log10(classical p-value)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order=NA,
#         type="p-value",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=-log10(0.05),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Position (ID)",
#         legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order=NA,
#         x.show=3500,
#         type="p-value",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=-log10(0.05),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Position (ID)",
#         legend.on=TRUE)
#
# # x-axis sgpv, y-axis -log10(classical p-value)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order="sgpv",
#         type="comparison",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=c(-log10(0.05), -log10(0.001)),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Second-generation p-value ranking",
#         legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order="sgpv",
#         x.show=4000,
#         type="comparison",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=c(-log10(0.05)),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Second-generation p-value ranking",
#         legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order="sgpv",
#         x.show=3300,
#         type="comparison",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=c(-log10(0.05)),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Second-generation p-value ranking",
#         legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order="sgpv",
#         x.show=3100,
#         type="comparison",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=c(-log10(0.05)),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Second-generation p-value ranking",
#         legend.on=TRUE)
# plotman(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
#         null.lo=-0.3, null.hi=0.3,
#         set.order="sgpv",
#         x.show=2000,
#         type="comparison",
#         p.values=-log10(leukstats$p.value),
#         ref.lines=c(-log10(0.05)),
#         int.pch=16, int.cex=0.4,
#         title.lab="Leukemia Example",
#         y.lab=expression("-log"[10]*"(p-value)"),
#         x.lab="Second-generation p-value ranking",
#         legend.on=TRUE)
