#' Plot upper confidence limit of total bycatch given none observed
#' 
#' \code{plot_uclnegobs} plots upper confidence limit of total bycatch vs 
#'   observer coverage when no bycatch is observed, given total fishery effort, 
#'   dispersion index, and confidence level. The function returns (1) minimum 
#'   observer coverage needed to fall within user-specified upper confidence 
#'   limit for bycatch when none was observed, and/or (2) the upper confidence 
#'   limit for bycatch given specified observer coverage and no observed bycatch. 
#'   
#' @param te an integer greater than 1. Total effort in fishery (e.g., trips 
#'   or sets).
#' @param d a number greater than or equal to 1. Dispersion index. The dispersion 
#'   index corresponds to the variance-to-mean ratio of effort-unit-level bycatch, 
#'   so \code{d = 1} corresponds to Poisson-distributed bycatch, and \code{d > 1} 
#'   to overdispersed bycatch.
#' @param cl a non-negative number less than or equal to 100. Confidence level
#'   for upper confidence limit of total bycatch (as percentage), given no bycatch 
#'   observed. 
#' @param targetucl a non-negative number. Target maximum upper confidence 
#'   limit for total bycatch given zero bycatch observed. If set to 0, no 
#'   corresponding minimum observer coverage will be highlighted or returned.
#' @param fixedoc a non-negative number between 0 and 100. Percent observer coverage 
#'   for which to return ucl value.
#' @param ymax a positive number. Upper limit for y-axis of plot.
#' @param silent logical. If \code{TRUE}, print output to terminal is suppressed.
#' @param showplot logical. If \code{FALSE}, plotting is suppressed.
#' @param ... additional arguments for compatibility with Shiny. 
#' 
#' @details
#' Upper confidence limits are based on the probability density function for 
#' the corresponding Poisson or negative binomial distribution. Upper confidence 
#' limits based on \code{d}+/-1 (as allowed by specification of d) are also plotted. 
#' If \code{fixedoc} specified, corresponding upper confidence limit is provided 
#' in printed output and returned object, but not in plot.
#' 
#' \strong{Caveat:} \code{plot_uclnegobs} assumes that (1) observer coverage is 
#' representative, (2) bycatch is in terms of individuals (not weight) per unit 
#' effort,and (3) the specified dispersion index reflects the highest level of 
#' any hierarchical variance (e.g., using dispersion index at trip level if 
#' greater than that at set level). Violating these assumptions will likely 
#' result in negatively biased projections of the observer coverage needed to 
#' meet a specified objective. More conservative (higher) projections can be 
#' obtained by using a higher dispersion index \code{d}. Users may want to 
#' explore uncertainty in dispersion index and in bycatch per unit effort by 
#' varying those inputs.
#' 
#' @return A list with components:
#'   \item{ucldat}{a data frame with the following fields for each coverage level included: 
#'   number of observed effort units (\code{nobs}), 
#'   proportion observer coverage (\code{pobs}), 
#'   upper confidence limit of total bycatch given none observed (\code{ucl}),
#'   and finite population correction (\code{fpc}) used in calculating \code{ucl}.}
#'   \item{targetucl}{specified target maximum upper confidence limit of bycatch.}
#'   \item{targetoc}{minimum observer coverage (as percentage) for which upper 
#'   confidence limit of bycatch is \code{targetucl} when none observed.}
#'   \item{fixedoc}{specified percentage observer coverage for which upper 
#'   confidence limit of bycatch is returned.}
#'   \item{fixednoc}{observer coverage (as effort) corresponding to \code{fixedoc}.}
#'   \item{fixedoc.ucl}{upper confidence limit of total bycatch corresponding 
#'   to zero bycatch observed in \code{fixedoc} coverage.}
#'   \item{te}{specified total effort.}
#'   \item{d}{specified dispersion index.}
#'   \item{cl}{specified confidence level.} 
#'   
#' @return Returned invisibly. 
#' 
#' @export 
plot_uclnegobs <- function(te, d = 2, cl = 95, targetucl = 0, fixedoc = 0, 
                           ymax = 100, showplot = TRUE, silent = FALSE, ...) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<=1) stop("te must be a positive integer > 1")
  if (d<1) stop("d must be >= 1")
  if (targetucl<0) stop("targetucl must be >= 0.")
  if (fixedoc)
    if (fixedoc<0 || fixedoc>100) stop("fixedoc must be >= 0 and <= 100.")
  if (ymax<=0) stop("ymax must be > 0")
  
  # get shiny flag if specified or set to FALSE
  myArgs <- match.call()
  if (!("as.shiny" %in% names(myArgs))) as.shiny <- FALSE
  else as.shiny <- myArgs$as.shiny
  
  # upper confidence limit of bycatch given none observed
  a <- 1 - cl/100
  dv <- c(d-1, d, d+1) # vary d
  if (te<1000) { oc <- 1:te
  } else { oc <- round(seq(0.001,1,0.001)*te) }
  df <- data.frame(nobs = oc, pobs = oc/te)
  df$ucl <- df$ucl.dl <- df$ucl.dh <- NA
  df$fpc <- sqrt((te - df$nobs)/(te-1))
  
  ucl.dl <- df$fpc * te * solveucl(a=a, d=dv[1], n=df$nobs)
  df$ucl <- df$fpc * te * solveucl(a=a, d=dv[2], n=df$nobs)
  ucl.dh <- df$fpc * te * solveucl(a=a, d=dv[3], n=df$nobs)
  
  # get minimum required observer coverage if target ucl specified
  if (targetucl) {
    itarget <- min(which(df$ucl <= targetucl))
    targetoc <- 100*ceiling_dec(df$pobs[itarget], 3)
  }
  
  if (fixedoc) {
    fixednoc <- round(fixedoc/100 * te)
    fixedoc.new <- 100*round(fixednoc/te, 3)
    # should return provided fixedoc or is rounded recalculated one closer?
    fixedoc <- ifelse(abs(fixedoc/100*te-fixednoc) < abs(fixedoc.new/100*te-fixednoc),
                      fixedoc, fixedoc.new)
    # how many decimal points in fixedoc?
    ndecpt.foc <- nchar(strsplit(as.character(fixedoc),".",fixed=TRUE)[[1]][2])
    ndecpt.foc <- ifelse(is.na(ndecpt.foc), 1, ndecpt.foc)
    # finite population correction
    fixedoc.fpc <- sqrt((te - fixednoc)/(te-1))
    # upper confidence limit
    fixedoc.ucl <- fixedoc.fpc * te * solveucl(a=a, d=d, n=fixednoc)
  }
  
  # plot
  if (showplot) {
    graphics::plot(100*df$pobs, log10(df$ucl), type="l", lty=1, lwd=2,
                   xlim=c(0,100), ylim=log10(c(utils::tail(ucl.dl,2)[1],min(max(ucl.dh),ymax))), 
                   xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxt="n",
                   xlab="Observer Coverage (%)", ylab="Upper Confidence Limit of Bycatch",
                   main=paste0("One-Tailed ", cl, "% UCL of Bycatch Given None Observed"))
    graphics::lines(100*df$pobs, log10(ucl.dl), lty=2, lwd=2)
    graphics::lines(100*df$pobs, log10(df$ucl), lty=2, lwd=2)
    graphics::lines(100*df$pobs, log10(ucl.dh), lty=3, lwd=2)
    graphics::axis(side=2, at=log10(c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000)), 
                   labels=c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000))
    if (targetucl) {
      graphics::lines(c(0,100),log10(rep(targetucl,2)), col=2, lwd=2, lty=4)
      graphics::points(targetoc, log10(df$ucl[itarget]), pch=8, col=2, cex=1.5, lwd=2)
      graphics::legend("topright", lty=c(2,1,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
                       legend=c(paste0("d=",dv[1]), paste0("d=",dv[2]), paste0("d=",dv[3]),
                                "target UCL", "min coverage"))
    } else {
      graphics::legend("topright", lty=c(2,1,3), lwd=2, col=1, 
                       legend=c(paste0("d=",dv[1]), paste0("d=",dv[2]), paste0("d=",dv[3])))
    }
  }
  
  # print recommended minimum observer coverage
  if (targetucl) {
    rec1 <- paste0("Minimum observer coverage to ensure that the upper confidence",
                 " limit of ", targetucl, " is not \nexceeded when no bycatch is ",
                 "observed is ", format(targetoc, nsmall=1), "%.\n")
  } else { rec1 <- "" }
  if (fixedoc) {
    rec2 <- paste0("Upper confidence limit for bycatch given none observed in ",
                 format(fixedoc, nsmall=ndecpt.foc), "% \n(", fixednoc, " effort units)",
                 " coverage is ", my_ceiling(fixedoc.ucl,3), ".\n")
  } else { rec2 <- "" }
  rec <- paste0(rec1, rec2)
  if (!as.shiny & !silent) 
    cat(paste0(rec, "Please review the caveats in the associated documentation.\n"))
  
  # return recommended minimum observer coverage
  if (as.shiny) {
    return(invisible(list(rec = rec)))
  } else { 
    return(invisible(list(ucldat=df, targetucl=ifelse(targetucl, targetucl, NA), 
                        targetoc=ifelse(targetucl, targetoc, NA), 
                        fixedoc=ifelse(fixedoc, fixedoc, NA), 
                        fixednoc=ifelse(fixedoc, fixednoc, NA),
                        fixedoc.ucl=ifelse(fixedoc, my_ceiling(fixedoc.ucl,3), NA),
                        te=te, d=d, cl=cl)))
  }
}
