#' Plot bycatch estimation CV vs. observer coverage
#' 
#' \code{plot_cv_obscov} plots projected bycatch estimation CVs vs observer 
#' coverage, and returns minimum observer coverage needed to achieve 
#' user-specified target CV and percentile. 
#'   
#' @param simlist list output from \code{sim_cv_obscov}.
#' @param targetcv a non-negative number less than 1. Target CV 
#'   (as a proportion). If \eqn{targetcv = 0}, no corresponding minimum observer 
#'   coverage will be highlighted.
#' @param silent logical. If silent = TRUE, print output to terminal is suppressed.
#' @param showplot logical. If plot = FALSE, plotting is suppressed.
#' 
#' @details
#' \strong{Caveat:} \code{plot_cv_obscov} assumes that (1) observer coverage is 
#' representative, (2) bycatch specified for \code{sim_cv_obscov}is in terms of 
#' individuals (not weight) per unit effort, and (3) the dispersion index 
#' specified for \code{sim_cv_obscov} reflects the highest level of any 
#' hierarchical variance (e.g., using dispersion index at trip level if greater 
#' than that at set level). Violating these assumptions will likely result in 
#' negatively biased projections of bycatch estimation CV for a given level of 
#' observer coverage. Users may want to explore uncertainty in dispersion index 
#' and in bycatch per unit effort by varying those inputs.See documentation for 
#' \code{sim_cv_obscov} for additional details.
#'   
#' @return A list with one component:
#'   \item{pobs}{minimum observer coverage in terms of percentage.} 
#' @return Returned invisibly. 
#'   
#' @export 
plot_cv_obscov <- function(simlist = simlist, targetcv = 0.3, 
                           showplot = TRUE, silent = FALSE) {
  
  # check input values
  if(targetcv<0 || targetcv>=1) stop("targetcv must be >= 0 and < 1")
  
  # get minimum required observer coverage
  if (targetcv)
    targetoc <- ifelse(simlist$te <= 20,
                       with(simlist$simsum, pobs[min(which(cvsim<=targetcv))]),
                       stats::approx(simlist$simsum$cvsim, simlist$simsum$pobs, targetcv)$y)
  
  # plot 
  if (showplot) {
    with(simlist$simsum, graphics::plot(100*pobs, cvsim, xlim=c(0,100), ylim=c(0,1), 
                                        xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,1,10),
                                        xlab="Observer Coverage (%)", ylab="CV of Bycatch Estimate",
                                        main="CV of Bycatch Estimate vs Observer Coverage"))
    #with(simsum, graphics::points(100*pobs, qcv))
    with(simlist$simsum, graphics::lines(100*pobs, cvsim))
    graphics::abline(h=1, v=100)
    # add minimum required observer coverage
    if (targetcv) {
      graphics::abline(h=targetcv, col=2, lwd=2, lty=2)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc*100, targetcv, pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      legpos <- ifelse(any(simlist$simsum$pobs > 0.7 & simlist$simsum$cvsim > 0.5), 
                       "bottomleft", "topright")
      graphics::legend(legpos, lty=c(2,0), pch=c(NA,8), col=c(2,2), lwd=c(2,2), 
                       pt.cex=1.5, y.intersp=1.1, legend=c("target CV", "min coverage"))
    }
  }
  
  # print recommended minimum observer coverage
  if (!silent) {
    if (targetcv)  {
      if (!is.na(targetoc)) {
        cat(paste0("Minimum observer coverage to achieve CV \u2264 ", targetcv, " is ", 
                   my_ceiling(targetoc*100,2), "%.\n"))
      } else {
        cat(paste0("Simulated observer coverage levels do not include range corresponding to ",
                   "minimum observer coverage to achieve CV \u2264 ", targetcv, ".\n"))
      }
    }
    cat(paste0("Please review the caveats in the associated documentation.\n"))
    cat(paste0("Note that results are interpolated from simulation-based projections and may vary slightly \n",
               "with repetition.\n"))
  }
  
  # return recommended minimum observer coverage
  if (targetcv) 
    return(invisible(list(pobs = my_ceiling(targetoc*100,2))))
}
