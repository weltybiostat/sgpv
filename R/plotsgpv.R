################################################################
##	Purpose: 	Plot interval estimates according to
##            	Second-Generation p-value rankings
##
##	Function:	plotsgpv
##	Version:	1.1
##
##	Author:		Valerie F. Welty and Jeffrey D. Blume
##	Date:		March 7, 2019
################################################################
#
#' Second-Generation \emph{p}-value plotting
#'
#' @description This function displays user supplied interval estimates (support intervals, confidence intervals, credible intervals, etc.) according to its associated second-generation \emph{p}-value ranking.
#'
#' @param est.lo A numeric vector of lower bounds of interval estimates. Values must be finite for interval to be drawn. Must be of same length as \code{est.hi}.
#' @param est.hi A numeric vector of upper bounds of interval estimates. Values must be finite for interval to be drawn. Must be of same length as \code{est.lo}.
#' @param null.lo A scalar representing the lower bound of null interval (indifference zone). Value must be finite.
#' @param null.hi A scalar representing the upper bound of null interval (indifference zone). Value must be finite.
#' @param set.order A numeric vector giving the desired order along the x-axis. If \code{set.order} is set to \code{sgpv}, the second-generation \emph{p}-value ranking is used. If \code{set.order} is set to \code{NA}, the original input ordering is used.
#' @param x.show A scalar representing the maximum ranking on the x-axis that is displayed. Default is to display all intervals.
#' @param null.col Coloring of the null interval (indifference zone). Default is Hawkes Blue: \code{rgb(208,216,232,max=255)}.
#' @param int.col Coloring of the intervals according to SGPV ranking. Default is \code{c("cornflowerblue","firebrick3","darkslateblue")} for SGPVs of \eqn{0}, in \eqn{(0,1)}, and \eqn{1} respectively.
#' @param int.pch Plotting symbol for interval endpoints. Default is \code{NA}, no symbol. Use \code{16} for small endpoints.
#' @param int.cex Size of plotting symbol for interval endpoints. Default is \eqn{0.4}.
#' @param plot.axis Toggle for default axis plotting. Default is \code{c(TRUE,TRUE)} for \eqn{(x-axis,y-axis)} respectively.
#' @param null.pt A scalar representing a point null hypothesis. Default is \code{NA}. If set, the function will draw a horizontal dashed black line at this location.
#' @param outline.zone Toggle for drawing a slim white outline around the null zone. Helpful visual aid when plotting many intervals. Default is \code{TRUE}.
#' @param title.lab Title text.
#' @param x.lab x-axis label.
#' @param y.lab y-axis label.
#' @param legend.on Toggle for plotting the legend. Default is \code{TRUE}.
#'
#' @details
#'
#' Use \code{set.order} to provide the classical p-value ranking. For example, if \code{pvalue.vector} is a vector of classical p-values, then set \code{set.order=order(pvalue.vector)} to sort the x-axis according to p-value rank.
#'
#' Interval estimates with infinite or undefined limits should be manually truncated or avoided altogether. While the sgpvalue funciton will handle these cases, this function assumes they have been truncated or removed because there is no standard way to plot them.
#'
#' @seealso \code{\link{sgpvalue}, \link{sgpower}, \link{fdrisk}}
#' @keywords
#' @export
#' @examples
#'
#' # Need leukstats data loaded locally
#' plotsgpv(est.lo=leuk.stats$ci.lo, est.hi=leuk.stats$ci.hi,
#'			null.lo=-0.3, null.hi=0.3,
#'			set.order=order(leuk.stats$p.value),
#'			x.show=7000,
#'			plot.axis=c("TRUE","FALSE"),
#'			null.pt=0, outline.zone=TRUE,
#'			title.lab="Leukemia Example", y.lab="Fold Change (base 10)",
#'			x.lab="Classical p-value ranking",
#'			legend.on=TRUE)
#' axis(side=2,at=round(log(c(1/1000,1/100,1/10,1/2,1,2,10,100,1000),
#'		base=10),2),labels=c("1/1000","1/100","1/10","1/2",1,2,10,100,1000),
#'		las=2)
#'
#'
#' @references
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation \emph{p}-values. \emph{The American Statistician}. In press. https://doi.org/10.1080/00031305.2018.1537893
#'

plotsgpv <- function( est.lo, est.hi, null.lo, null.hi,
		set.order="sgpv", x.show=NA, null.col=rgb(208,216,232,max=255),
		int.col=c("cornflowerblue","firebrick3","darkslateblue"),
		int.pch=NA, int.cex=0.4, plot.axis=c(TRUE,TRUE),
		null.pt=NA, outline.zone=TRUE,
		title.lab="Title", x.lab="Position (by set.order)", y.lab="Outcome label",
		legend.on=TRUE ){

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

	#### Set plot limits
	x.max 	<- length(est.lo)
	x 		<- seq(1,x.max,1)

	x.limits <- c(1,min(x.show,x.max,na.rm=TRUE))
	y.limits <- c(floor(min(est.lo,est.hi)),ceiling(max(est.lo,est.hi)))

	#### Compute SGPVs
	sgpv <- sgpvalue(est.lo=est.lo, est.hi=est.hi, null.lo=null.lo, null.hi=null.hi)
	sgpv.combo <- ifelse(sgpv$p.delta==0,-sgpv$delta.gap,sgpv$p.delta)

	#### Set order of x-axis
	if (is.na(set.order[1]))  {set.order <- x}
	if (set.order[1]=="sgpv") {set.order <- order(sgpv.combo)}

	#### Subset intervals by SGPV value for coloring
	gap.marker <- 1*is.na(sgpv$delta.gap[set.order])

	set.out  <- x[gap.marker[x]==0]
	set.in   <- x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]==1)]
	set.both <- x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]<1)]

	#### Plotting
	plot(1, null.pt, ylim=y.limits, xlim=x.limits, type="n", yaxt="n", xaxt="n",
		ylab=y.lab, xlab=x.lab, main=title.lab)

	if (plot.axis[1]=="TRUE") {axis(side=1)}
	if (plot.axis[2]=="TRUE") {axis(side=2,las=2)}

	rect(1, null.lo, x.max, null.hi, col=null.col, border=NA)

	## Intervals where 0<SGPV<1
	points(x[set.both],est.lo[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
	points(x[set.both],est.hi[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
	segments(x[set.both], est.lo[set.order][set.both], x[set.both], est.hi[set.order][set.both], lty=1, col=int.col[1])

	## Intervals where SGPV==1
	points(x[set.in],est.lo[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
	points(x[set.in],est.hi[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
	segments(x[set.in], est.lo[set.order][set.in], x[set.in], est.hi[set.order][set.in], lty=1, col=int.col[3])

	## Intervals where SGPV==0
	points(x[set.out],est.lo[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
	points(x[set.out],est.hi[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
	segments(x[set.out], est.lo[set.order][set.out], x[set.out], est.hi[set.order][set.out], lty=1, col=int.col[2])

	## Detail indifference zone
	abline(h=null.pt,lty=2)

	if (outline.zone==TRUE) {
			abline(h=null.lo, col="white", lty=1, lwd=0.8)
			abline(h=null.hi, col="white", lty=1, lwd=0.8)
							}
	#### Legend
	if (legend.on==TRUE) {
	legend("topright",c("Interval Null", expression("p"[delta]*" = 0"),
						expression("0 < p"[delta]*" < 1"), expression("p"[delta]*" = 1")),
			lty=1,col=c(null.col,int.col[2],int.col[1],int.col[3]),lwd=c(6,1.5,1.5,1.5),bty="n")
						}
}

###
##
#
