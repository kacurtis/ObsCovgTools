# Hidden function to solve for upper confidence limit of bpue at given confidence 
# level when zero in observed
solveucl <- function(bpue, n, a, d) {
  pz <- a^(1/n)
  if(d==1) { return(abs(stats::ppois(0, bpue) - pz))
  } else { return(abs(stats::pnbinom(0, size=(bpue/(d-1)), prob=1/d) - pz)) }
}
  

#' Plot upper confidence limit of total bycatch for none observed
#' 
#' \code{plot_uclnegobs} plots upper confidence limit of total bycatch vs 
#'   observer coverage when no bycatch observed, given total effort in 
#'   sets/hauls, negative binomial dispersion index, and confidence level.
#'   
#' @param te an integer greater than 1. Total effort in fishery (sets/hauls).
#' @param d a number greater than or equal to 1. Negative binomial dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param cl a non-negative number less than or equal to 100. Confidence level
#'   for upper confidence limit of bycatch (as percentage), given no bycatch 
#'   observed. 
#' @param target.ucl a non-negative number. Maximum allowable upper confidence 
#'   limit for bycatch given zero bycatch observed. If 0, no corresponding 
#'   minimum observer coverage will be highlighted.
#' @param silent logical. If silent = TRUE, print output to terminal is suppressed.
#' @param showplot logical. If plot = FALSE, plotting is suppressed.
#' 
#' @details
#' Upper confidence limits are based on the probability density function for the 
#' corresponding Poisson or negative binomial distribution.
#' 
#' Note that unlike \code{plot_cv_obscov}, \code{plot_uclnegobs} is designed 
#' as a one-step tool. 
#'   
#' \strong{Caveat:} \code{plot_uclnegobs} assumes representative observer coverage 
#' and no hierarchical sources of variance (e.g., vessel- or trip-level variation). 
#' Violating these assumptions will likely result in negatively biased projections 
#' of the upper confidence limit of total bycatch given zero observed. More 
#' conservative projections can be obtained by using higher-level units of effort 
#' (e.g., \code{te} as number of trips instead of number of sets/hauls).
#' 
#' @return A list with components:
#'   \item{ucldat}{a tibble with the following fields: 
#'   proportion observer coverage (\code{pobs}), number of observed trips/sets
#'   (\code{nobs}), and upper confidence limit of total bycatch given none 
#'   observed (\code{ucl}).}
#'   \item{d}{the negative binomial dispersion index used.}
#'   \item{cl}{specified confidence level.} 
#'   \item{target.ucl}{maximum upper confidence limit of bycatch specified.}
#'   \item{target.oc}{minimum observer coverage for which upper confidence 
#'   limit of bycatch is (\code{target.ucl}) when none observed.}
#' @return Returned invisibly. 
#' 
#' @export 
plot_uclnegobs <- function(te, d = 2, cl = 95, target.ucl = 0, showplot = TRUE, 
                            silent = FALSE) {
  # check input values
  if ((ceiling(te) != floor(te)) || te<=1) stop("te must be a positive integer > 1")
  if (d<1) stop("d must be >= 1")
  if (target.ucl<0) stop("target.ucl must be >= 0")
  # upper confidence limit of bycatch given none observed
  a <- 1 - cl/100
  if (te<1000) { oc <- data.frame(nobs = 1:te, pobs = (1:te)/te)
  } else { oc <- data.frame(nobs = round(seq(0.001,1,0.001)*te), 
                            pobs = round(seq(0.001,1,0.001)*te)/te)
  }
  oc$ucl <- NA
  for (i in 1:nrow(oc)) oc$ucl[i] <- te * optim(0.1, fn=solveucl, n=oc$nobs[i], a=a, d=d, method="BFGS")$par
  if (target.ucl) {
    itarget <- min(which(oc$ucl <= target.ucl))
    targetoc <- oc$pobs[itarget]
    targetnos <- oc$nobs[itarget]
  }
  
  # plot ###****
  if (showplot) {
    opar <- graphics::par(no.readonly = TRUE)
    graphics::par(xpd=TRUE)
    graphics::plot(100*oc$pobs, 100*(oc$ucl), type="l", lty=1, lwd=2,
                   xlim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), 
                   xlab="Observer Coverage (%)", ylab="Upper Confidence Limit of Bycatch",
                   main="Upper One-Tailed ", cl, "% Confidence Limit of Bycatch Given None Observed")
    if (target.ppos) {
      graphics::lines(c(0,100),rep(target.ppos,2), col=2, lwd=2, lty=4)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc*100, 100*oc$ppc[itarget], pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      graphics::legend("bottomleft", lty=c(1,2,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
                       legend=c("in observed effort if total bycatch > 0", "in observed effort",
                                "in total effort", "in target coverage if total bycatch > 0", "min coverage"))
    } else {
      graphics::legend("bottomleft", lty=c(1,2,3), lwd=2, col=1, 
                       legend=c("in observed effort if total bycatch > 0","in observed effort","in total effort"))
    }
    graphics::par(opar)
  }
  # return recommended minimum observer coverage
  if (target.ppos) {
    if (!silent) {
      cat(paste0("The probability that any bycatch occurs in the given total effort is ", 
                 signif(100*ppt,3), "%.\n",
                 "Minimum observer coverage to achieve at least ", target.ppos, 
                 "% probability of observing \nbycatch when total bycatch is positive is ", 
                 my_ceiling(targetoc*100,3), "% (", targetnos, " sets).\n",
                 "Please review the caveats in the associated documentation.\n"))
    }
  } else {
    if (!silent) {
      cat(paste0("The probability that any bycatch occurs in the given total effort",
                 " is ", signif(100*ppt,2), "%.\n",
                 "Please review the caveats in the associated documentation.\n"))
    }
  }
  return(invisible(list(oc=oc, target.ucl=ifelse(target.ucl, target.ucl, NA), 
              targetoc=ifelse(target.ucl, targetoc, NA), 
              targetnos=ifelse(target.ucl, targetnos, NA),
              te=te, d=d, cl=cl)))
}
