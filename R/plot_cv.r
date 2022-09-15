#' Plot bycatch estimation CV vs. observer coverage
#' 
#' \code{plot_cv} plots projected bycatch estimation CVs vs observer 
#' coverage, and returns minimum observer coverage needed to achieve 
#' user-specified target CV and percentile. 
#'   
#' @param te an integer greater than 1. Total effort in fishery (e.g., trips 
#'   or sets).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion index. The dispersion 
#'   index corresponds to the variance-to-mean ratio of effort-unit-level bycatch, 
#'   so \code{d = 1} corresponds to Poisson-distributed bycatch, and \code{d > 1} 
#'   to overdispersed bycatch.
#' @param targetcv a non-negative number less than 1. Target CV (as a proportion). 
#'   If set to 0, no corresponding minimum observer coverage will be 
#'   highlighted or returned.
#' @param silent logical. If \code{TRUE}, print output to terminal is suppressed.
#' @param showplot logical. If \code{FALSE}, plotting is suppressed.
#' @param ... additional arguments for compatibility with Shiny. 
#' 
#' @details
#' \strong{Caveat:} \code{plot_cv} assumes that (1) observer coverage is 
#' representative, (2) bycatch (\code{bpue}) is in terms of individuals (not 
#' weight) per unit effort, and (3) the specified dispersion index reflects 
#' the highest level of any hierarchical variance (e.g., using dispersion index 
#' at trip level if greater than that at set level). Violating these assumptions 
#' will likely result in negatively biased projections of bycatch estimation CV 
#' for a given level of observer coverage. More conservative projections can be 
#' obtained by using a higher dispersion index \code{d}. Users may want to 
#' explore uncertainty in dispersion index and in bycatch per unit effort by 
#' varying those inputs.
#'   
#' @return If \code{targetcv} is non-zero, a list with one component:
#'   \item{targetoc}{minimum observer coverage in terms of percentage.} 
#' @return Returned invisibly. 
#'   
#' @export 
plot_cv <- function(te, bpue, d = 2, targetcv = 0.3, 
                           showplot = TRUE, silent = FALSE, ...) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<2) stop("te must be a positive integer > 1")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  if(targetcv<0 || targetcv>=1) stop("targetcv must be >= 0 and < 1")
  
  # get shiny flag if specified or set to FALSE
  myArgs <- match.call()
  if (!("as.shiny" %in% names(myArgs))) as.shiny <- FALSE
  else as.shiny <- myArgs$as.shiny
  
  # get CV of total bycatch estimates for range of observer coverage
  if (te<1000) { oc <- 1:te 
  } else { oc <- round(seq(0.001,1,0.001)*te) }
  df <- data.frame(nobs = oc, pobs = oc/te)
  fpc <- (te-df$nobs)/(te-1)   # finite population correction, Cochran 1973
  s2 <- d*bpue   # variance of NB or Poisson distribution
  df$cv <- sqrt((fpc * s2)/df$nobs) / bpue
  
  # get minimum required observer coverage if targetcv provided
  if (targetcv) {
    itarget <- min(which(df$cv <= targetcv))
    targetoc <- 100*ceiling_dec(df$pobs[itarget], 3)
  }
  
  # plot 
  if (showplot) {
    with(df, graphics::plot(100*pobs, cv, type="l", lty=1, lwd=2,
                            xlim=c(0,100), ylim=c(0,1), 
                            xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,1,10),
                            xlab="Observer Coverage (%)", ylab="CV of Bycatch Estimate",
                            main="CV of Bycatch Estimate vs Observer Coverage"))
    with(df, graphics::lines(100*pobs, cv))
    graphics::abline(h=1, v=100)
    # add minimum required observer coverage
    if (targetcv) {
      graphics::abline(h=targetcv, col=2, lwd=2, lty=2)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc, targetcv, pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      legpos <- ifelse(any(df$pobs > 0.7 & df$cv > 0.5), 
                       "bottomleft", "topright")
      graphics::legend(legpos, lty=c(2,0), pch=c(NA,8), col=c(2,2), lwd=c(2,2), 
                       pt.cex=1.5, y.intersp=1.1, legend=c("target CV", "min coverage"))
    }
  }
  
  # print recommended minimum observer coverage
  if (targetcv)  {
    if (!is.na(targetoc)) {
      rec <- paste0("Minimum observer coverage to achieve CV \u2264 ", targetcv, 
                     " is ", format(targetoc, nsmall=1), "%.\n")
    } else {
      rec <- paste0("Simulated observer coverage levels do not include range corresponding to ",
                 "minimum observer coverage to achieve CV \u2264 ", targetcv, ".\n")
    }
  } else { rec <- "" }
  if (!as.shiny & !silent) 
    cat(paste0(rec, "Please review the caveats in the documentation.\n"))
  
  # return recommended minimum observer coverage
  if (as.shiny) {
    return(invisible(list(rec = rec)))
  } else { 
    if (targetcv)
      return(invisible(list(targetoc = targetoc)))
  }
}
