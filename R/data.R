#' Test Statistics from Gloub (1999) Leukemia data set
#'
#' Data are from 7218 gene specific t-tests for a difference in mean expression (on the log scale; AML versus ALL) in the Gloub data set (1999). Data are from 72 patients using a pooled t-test (df=70). Included in the dataframe are the following: t-statistic (\code{t.stat}), p-value (\code{p.value}), CI lower limit (\code{ci.lo}), CI upper limit (\code{ci.hi}), estimate (\code{estimate}), standard error (\code{se}).
#'
#' @docType data
#'
#' @usage data(leukstats)
#'
#' @format An object of class \code{data.frame}. Includes the following: t-statistic (\code{t.stat}), p-value (\code{p.value}), CI lower limit (\code{ci.lo}), CI upper limit (\code{ci.hi}), estimate (\code{estimate}), standard error (\code{se}).
#'
#' @keywords datasets
#'
#' @references Gloub (1999) and used in Blume et. al. (2018) PlosONE.
#'
#' Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation \emph{p}-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE} 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
#'
#' @source https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)
#'
#' @examples
#' data(leukstats)
#' order(leukstats$p.value)
"leukstats"
