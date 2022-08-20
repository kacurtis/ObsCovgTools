#' Plot probability of positive bycatch vs observer coverage
#' 
#' \code{plot_probposobs} plots (1) probability of observing at least one bycatch
#'   event vs observer coverage and (2) probability of any bycatch occurring in 
#'   total fishery effort, given total fishery effort, bycatch per unit effort, 
#'   and dispersion index. 
#'   
#' @param te an integer greater than 1. Total effort in fishery (e.g., trips 
#'   or sets).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of effort-unit-level bycatch, so \eqn{d = 1} corresponds to Poisson-
#'   distributed bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param targetppos a non-negative number less than or equal to 100. Target 
#'   probability of positive observed bycatch (as percentage), given positive 
#'   bycatch in total effort. If 0, no corresponding minimum observer coverage 
#'   will be highlighted.
#' @param silent logical. If silent = TRUE, print output to terminal is suppressed.
#' @param showplot logical. If plot = FALSE, plotting is suppressed.
#' 
#' @details  
#' Probabilities are based on the probability density function for the 
#' corresponding Poisson or negative binomial distribution.
#' 
#' The probability that any bycatch occurs in the given total effort is shown
#' by the horizontal black dotted line. The conditional probability of observing 
#' any bycatch if it occurs is shown by the solid black line.  The product of 
#' these first two probabilities gives the absolute probability of observing any
#' bycatch (dashed black line).The minimum observer coverage to achieve the target 
#' probability of observing bycatch if it occurs (x-axis value of red star) is 
#' where the conditional bycatch detection probability (solid black line) 
#' intersects with the target probability (red dash-dot line).
#' 
#' \strong{Caveat:} \code{plot_probposobs} assumes that (1) observer coverage is 
#' representative, (2) bycatch (\code{bpue}) is in terms of individuals (not 
#' weight) per unit effort, and (3) the specified dispersion index reflects 
#' the highest level of any hierarchical variance (e.g., using dispersion index 
#' at trip level if greater than that at set level). Violating these assumptions 
#' will likely result in positively biased projections of the probability of 
#' observing bycatch at a given level of observer coverage. More conservative 
#' projections can be obtained by using a higher dispersion index \code{d}. Users 
#' may want to explore uncertainty in dispersion index and in bycatch per unit 
#' effort by varying those inputs.
#' 
#' @return A list with components:
#'   \item{pobs}{minimum observer coverage in terms of percentage.} 
#'   \item{nobs}{corresponding observed effort.}
#'   \item{ppos.te}{probability of any bycatch occurring in total effort}
#' @return Returned invisibly. 
#' 
#' @export 
plot_probposobs <- function(te, bpue, d = 2, targetppos = 95, showplot = TRUE, 
                            silent = FALSE) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<=1) stop("te must be a positive integer > 1")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  if (targetppos<0 || targetppos>100) stop("targetppos must be >= 0 and <= 100")
  
  # percent probablity of positive observed bycatch
  if (te<1000) { oc <- 1:te 
  } else { oc <- round(seq(0.001,1,0.001)*te) }
  df <- data.frame(nobs = oc, pobs = oc/te)
  df$pp <- 1-probnzeros(df$nobs, bpue, d)   # probability of positive observed bycatch
  ppt <- utils::tail(df$pp,1)   # probability of positive bycatch in total effort
  df$ppc <- df$pp/ppt
  if (targetppos) {
    itarget <- min(which(df$ppc >= targetppos/100))
    targetoc <- df$pobs[itarget]
    targetnoc <- df$nobs[itarget]
  }
  
  # plot
  if (showplot) {
    opar <- graphics::par(no.readonly = TRUE)
    graphics::par(xpd=TRUE)
    graphics::plot(100*df$pobs, 100*(df$ppc), type="l", lty=1, lwd=2,
                   xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
                   xlab="Observer Coverage (%)", ylab="Probability of Positive Bycatch (%)",
                   main="Probability of Positive Bycatch")
    graphics::lines(x=c(0,100),y=rep(100*ppt,2),lwd=3, lty=3)
    graphics::lines(100*df$pobs, 100*df$pp, lwd=3, lty=2)
    graphics::lines(100*df$pobs, 100*df$ppc, lwd=2)
    legpos <- ifelse(any(df$pobs > 0.6 & df$pp < 0.3 ), "topleft", "bottomright")
    if (targetppos) {
      graphics::lines(c(0,100),rep(targetppos,2), col=2, lwd=2, lty=4)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc*100, 100*df$ppc[itarget], pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      graphics::legend(legpos, lty=c(1,2,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
                       legend=c("in observed effort if total bycatch > 0", "in observed effort",
                                "in total effort", "in target coverage if total bycatch > 0", "min coverage"))
    } else {
      graphics::legend(legpos, lty=c(1,2,3), lwd=2, col=1, 
                       legend=c("in observed effort if total bycatch > 0","in observed effort","in total effort"))
    }
    graphics::par(opar)
  }
  
  # print recommended minimum observer coverage
  if (!silent) {
    cat(paste0("The probability that any bycatch occurs in the given total effort is ", 
               signif(100*ppt,3), "%.\n"))
    if (targetppos) 
      cat(paste0("Minimum observer coverage to achieve at least ", targetppos, 
                 "% probability of observing \nbycatch when total bycatch is positive is ", 
                 my_ceiling(targetoc*100,3), "% (", targetnoc, " trips or sets).\n"))
    cat(paste0("Please review the caveats in the associated documentation.\n"))
  }
  
  # return recommended minimum observer coverage
  return(invisible(list(pobs = ifelse(targetppos, my_ceiling(targetoc*100,3), NA), 
                        nobs = ifelse(targetppos, my_ceiling(targetoc*te,3), NA),
                        ppos.te = signif(100*ppt,3))))
}
