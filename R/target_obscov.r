#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL 


## Quiets concerns of R CMD check re: the .'s that appear in pipelines
## and the "n" that is produced by dplyr::count() in a pipeline
if (getRversion() >= "2.15.1") utils::globalVariables(c("n"))


# Hidden functions to execute progress bar in Shiny
progress_init <- function(shiny.progress = FALSE) {
  if (shiny.progress) { return(NULL)
  } else { return(utils::txtProgressBar(style=3)) }
} 
  
progbar <-  function(i, total, pb, shiny.progress = FALSE) {
  if (shiny.progress) {
    shiny::incProgress(500 / total)
    return(NULL)
  } else {
    utils::setTxtProgressBar(pb, i/total)
    return(pb)
  }
}


# Hidden function to round up to specified significant digits
# Extended from code by JasonWang on stackoverflow at 
# https://stackoverflow.com/questions/37583715/round-up-values-to-a-specific-significant-figure-in-r
my_ceiling <- function(x, s){
     num_string <- format(x, scientific=TRUE)
     n <- strsplit(num_string, "e")
     n1 <- sapply(n, function(x) as.numeric(x[1]))
     n2 <- sapply(n, function(x) as.numeric(x[2]))
     ceiling(n1*10^(s-1))/(10^(s-1)) * 10^(n2)
}


#' Simulate CV response to observer coverage
#'
#' \code{sim_cv_obscov} simulates bycatch estimation CVs for a range 
#' of observer coverage levels, given bycatch rate, dispersion 
#' index, and total fishery effort. 
#' 
#' \code{sim_cv_obscov} runs \code{nsim} simulations per level of observer 
#' coverage, from the larger of 0.1\% or one set/haul to 100\%. Simulated 
#' bycatch estimates are calculated as mean observed bycatch per unit effort. 
#' CV at each observer coverage level is calculated as the square root of mean 
#' square estimation error divided by "true" bycatch (\code{bpue}).
#' 
#' Warning: Large total effort (>100K sets/hauls) may require several minutes 
#' of execution time. Increasing \code{nsim} from the default of 1000 will 
#' also increase execution time. 
#' 
#' \strong{Caveat:} \code{sim_cv_obscov} assumes representative observer coverage 
#' and no hierarchical sources of variance (e.g., vessel- or trip-level variation). 
#' Violating these assumptions will likely result in negatively biased projections of 
#' bycatch estimation CV for a given level of observer coverage. More conservative 
#' projections can be obtained by using higher-level units of effort (e.g., 
#' \code{bpue} as mean bycatch per trip instead of bycatch per set/haul, and 
#' \code{te} as number of trips instead of number of sets/hauls).
#' 
#' @param te an integer greater than one. Total effort in fishery (sets/hauls).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param nsim a positive integer. Number of simulations to run.
#' @param ...  additional arguments for compatibility with Shiny.
#'   
#' @return A list with components:
#'   \item{simsum}{a tibble with one row per observer coverage level and the 
#'   following fields: simulated proportion observer coverage (\code{pobs}), 
#'   number of observed trips/sets (\code{nobs}), and bycatch estimation CV 
#'   (\code{cvsim}).} 
#'   \item{simdat}{a tibble with one row per simulation and the following fields: 
#'   simulated proportion observer coverage (\code{pobs}), number of observed 
#'   sets (\code{nobs}), true (realized) bycatch per unit effort (\code{tbpue}), 
#'   observed bycatch per unit effort (\code{obpue}), and error of observed bycatch 
#'   per unit effort (\code{oberr} = \code{obpue} - \code{tbpue}).}
#'   \item{bpue}{the bycatch per unit effort used.}
#'   \item{d}{the dispersion index used.}
#'   
#' @export 
sim_cv_obscov <- function(te, bpue, d = 2, nsim = 1000, ...) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<2) stop("te must be a positive integer > 1")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  if ((ceiling(nsim) != floor(nsim)) || nsim<=0) stop("nsim must be a positive integer")
  
  # simulate observer coverage and bycatch estimation
  if (te<20) { oc <- 1:te/te 
  } else { oc <- c(seq(0.001,0.005,0.001), seq(0.01,0.05,0.01), seq(0.10,1,0.05)) }
  simdat <- tibble::tibble(pobs = rep(oc, nsim), 
                           nobs = round(.data$pobs * te)) %>% 
    dplyr::filter(.data$nobs > 0) %>% 
    dplyr::mutate(tbpue=NA, obpue=NA, oberr=NA)
  set.seed(Sys.time())
  pb <- progress_init(...)
  
  for (i in 1:nrow(simdat)) {
    sets <- if(d==1) { Runuran::urpois(te, bpue) 
    } else { Runuran::urnbinom(te, size=(bpue/(d-1)), prob=1/d) }
    obsets <- sample(sets, simdat$nobs[i])
    simdat$tbpue[i] <- mean(sets)
    simdat$obpue[i] <- mean(obsets)
    simdat$oberr[i] <- simdat$obpue[i] - simdat$tbpue[i]

    if (i %% 500 == 0) {
      pb <- progbar(i, nrow(simdat), pb, ...)
    }
  }
  
  simsum <- simdat %>% dplyr::group_by(.data$pobs) %>% 
    dplyr::summarize(nobs=unique(.data$nobs),
                     cvsim=sqrt(mean(.data$oberr^2))/bpue)
  return(list(simsum=simsum, simdat=simdat, te=te, bpue=bpue, d=d))
}


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
#' \strong{Caveat:} \code{sim_cv_obscov} assumes representative observer 
#' coverage and no hierarchical sources of variance (e.g., vessel- or trip-level 
#' variation). Violating these assumptions will likely result in negatively biased 
#' projections of bycatch estimation CV for a given level of observer coverage. 
#' See documentation for \code{sim_obs_cov} for additional details.
#'   
#' @return A list with components:
#'   \item{pobs}{minimum observer coverage in terms of percentage.} 
#'   \item{nobs}{corresponding observed effort.}
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
              my_ceiling(targetoc*100,2), "% (", my_ceiling(targetoc*simlist$te,2), " hauls).\n"))
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
    return(invisible(list(pobs = my_ceiling(targetoc*100,2), nobs=my_ceiling(targetoc*simlist$te,2))))
}


#' Get probability of zero bycatch given effort
#' 
#' \code{probnzeros} returns probability of zero bycatch in a specified number 
#' of sets/hauls, given bycatch per unit effort and dispersion index. 
#' 
#' @param n a vector of positive integers. Observed effort levels (in terms of 
#'   sets/hauls) for which to calculate probability of zero bycatch.
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#'   
#' @details
#' Calculated from the probability density at zero of the corresponding Poisson
#' (\eqn{d = 1}) or negative binomial (\eqn{d < 1}) distribution.
#' 
#' \strong{Caveat:} \code{probnzeros} assumes representative observer coverage 
#' and no hierarchical sources of variance (e.g., vessel- or trip-level variation). 
#' Violating these assumptions will likely result in negatively biased projections 
#' of the probability of observing zero bycatch at a given level of observer coverage. 
#' More conservative projections can be obtained by using higher-level units of effort 
#' (e.g., \code{bpue} as mean bycatch per trip instead of bycatch per set/haul, and 
#' \code{n} as number of trips instead of number of sets/hauls).
#'   
#' @return Vector of same length as \code{n} with probabilities of zero bycatch. 
#' @return Returned invisibly
#' 
#' @export
probnzeros <- function(n, bpue, d) {
  
  # check input values
  if (any((ceiling(n) != floor(n)) | n<1)) stop("n must be a vector of positive integers")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  
  # calculate probability of observing zero bycatch in n units of effort
  ## one unit
  pz <- if(d==1) { stats::ppois(0, bpue)
  } else { stats::pnbinom(0, size=(bpue/(d-1)), prob=1/d) }
  ## n units
  pnz <- pz^n 
  
  return(invisible(pnz))
}


#' Plot probability of positive bycatch vs observer coverage
#' 
#' \code{plot_probposobs} plots (1) probability of observing at least one bycatch
#'   event vs observer coverage and (2) probability of any bycatch occurring in 
#'   total effort, given total effort in sets/hauls, bycatch per unit effort, and 
#'   dispersion index. 
#'   
#' @param te an integer greater than 1. Total effort in fishery (sets/hauls).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param target.ppos a non-negative number less than or equal to 100. Target 
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
#' obability of observing bycatch if it occurs (x-axis value of red star) is 
#' where the conditional bycatch detection probability (solid black line) 
#' intersects with the target probability (red dash-dot line).
#' 
#' Note that unlike \code{plot_cv_obscov}, \code{plot_probposobs} is designed 
#' as a one-step tool, and does not take output from user calls to 
#' \code{probnzeros}. 
#'   
#' \strong{Caveat:} \code{plot_probposobs} assumes representative observer coverage 
#' and no hierarchical sources of variance (e.g., vessel- or trip-level variation). 
#' Violating these assumptions will likely result in positively biased projections 
#' of the probability of observing bycatch at a given level of observer coverage. 
#' More conservative projections can be obtained by using higher-level units of effort 
#' (e.g., \code{bpue} as mean bycatch per trip instead of bycatch per set/haul, and 
#' \code{te} as number of trips instead of number of sets/hauls).
#' 
#' @return A list with components:
#'   \item{pobs}{minimum observer coverage in terms of percentage.} 
#'   \item{nobs}{corresponding observed effort.}
#'   \item{ppos.te}{probability of any bycatch occurring in total effort}
#' @return Returned invisibly. 
#' 
#' @export 
plot_probposobs <- function(te, bpue, d = 2, target.ppos = 80, showplot = TRUE, 
                            silent = FALSE) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<=1) stop("te must be a positive integer > 1")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  if (target.ppos<0 || target.ppos>100) stop("target.ppos must be >= 0 and <= 100")
  
  # percent probablity of positive observed bycatch
  if (te<1000) { oc <- 1:te 
  } else { oc <- round(seq(0.001,1,0.001)*te) }
  df <- data.frame(nobs = oc, pobs = oc/te)
  df$pp <- 1-probnzeros(df$nobs, bpue, d)   # probability of positive observed bycatch
  ppt <- utils::tail(df$pp,1)   # probability of positive bycatch in total effort
  df$ppc <- df$pp/ppt
  if (target.ppos) {
    itarget <- min(which(df$ppc >= target.ppos/100))
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
    if (target.ppos) {
      graphics::lines(c(0,100),rep(target.ppos,2), col=2, lwd=2, lty=4)
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
    if (target.ppos) 
      cat(paste0("Minimum observer coverage to achieve at least ", target.ppos, 
                "% probability of observing \nbycatch when total bycatch is positive is ", 
                my_ceiling(targetoc*100,3), "% (", targetnoc, " sets).\n"))
    cat(paste0("Please review the caveats in the associated documentation.\n"))
  }
  
  # return recommended minimum observer coverage
  return(invisible(list(pobs = ifelse(target.ppos, my_ceiling(targetoc*100,3), NA), 
                        nobs = ifelse(target.ppos, my_ceiling(targetoc*te,3), NA),
                        ppos.te = signif(100*ppt,3))))
}


# Hidden function to solve for upper one-tailed confidence limit (1-a) of bpue 
# when zero in observed effort n, given overdispersion d
solveucl <- function(a, d, n) {
  if (d==1) {
    bpue <- -1 * log(a)/n
  } else {
    bpue <- (log(a)*(d-1))/(n*log(1/(d-1)/(1/(d-1)+1)))
  }
  return(bpue)
}


#' Plot upper confidence limit of total bycatch for none observed
#' 
#' \code{plot_uclnegobs} plots upper confidence limit of total bycatch vs 
#'   observer coverage when no bycatch is observed, given total effort in 
#'   sets/hauls, dispersion index, and confidence level.
#'   
#' @param te an integer greater than 1. Total effort in fishery (sets/hauls).
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param cl a non-negative number less than or equal to 100. Confidence level
#'   for upper confidence limit of bycatch (as percentage), given no bycatch 
#'   observed. 
#' @param target.ucl a non-negative number. Maximum allowable upper confidence 
#'   limit for bycatch given zero bycatch observed. If 0, no corresponding 
#'   minimum observer coverage will be highlighted.
#' @param fixed.oc a positive number between 1 and 100. Fixed observer coverage 
#'   (as percentage) for which ucl value is desired.
#' @param ymax a positive number. Maximum bycatch level to display in plot.
#' @param silent logical. If silent = TRUE, print output to terminal is suppressed.
#' @param showplot logical. If plot = FALSE, plotting is suppressed.
#' 
#' @details
#' Upper confidence limits are based on the probability density function for 
#' the corresponding Poisson or negative binomial distribution. Upper confidence 
#' limits for dispersion +/-1 are also plotted. If \code{fixed.oc} specified, 
#' corresponding upper confidence limit is provided in printed output only.
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
#'   \item{ucldat}{a tibble with the following fields for each coverage level included: 
#'   number of observed trips/sets (\code{nobs}), 
#'   proportion observer coverage (\code{pobs}), 
#'   upper confidence limit of total bycatch given none observed (\code{ucl}),
#'   and finite population correction (\code{fpc}) used in calculating \code{ucl}.}
#'   \item{target.ucl}{maximum upper confidence limit of bycatch specified.}
#'   \item{targetoc}{minimum observer coverage (as proportion) for which upper 
#'   confidence limit of bycatch is \code{target.ucl} when none observed.}
#'   \item{targetnoc}{minimum observer coverage (as effort) for which upper 
#'   confidence limit of bycatch is \code{target.ucl} when none observed.}
#'   \item{te}{specified total effort.}
#'   \item{d}{specified dispersion index.}
#'   \item{cl}{specified confidence level.} 
#'   
#' @return Returned invisibly. 
#' 
#' @export 
plot_uclnegobs <- function(te, d = 2, cl = 95, target.ucl = 0, fixed.oc = NA, 
                           ymax = 100, showplot = TRUE, silent = FALSE) {
  
  # check input values
  if ((ceiling(te) != floor(te)) || te<=1) stop("te must be a positive integer > 1")
  if (d<1) stop("d must be >= 1")
  if (target.ucl<0) stop("target.ucl must be >= 0.")
  if (!is.na(fixed.oc))
    if (fixed.oc<=0 || fixed.oc>100) stop("fixed.oc must be > 0 and <= 100.")
  if (ymax<=0) stop("ymax must be > 0")
  
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
  
  # identify target observer coverage if target ucl specified
  if (target.ucl) {
    itarget <- min(which(df$ucl <= target.ucl))
    targetoc <- df$pobs[itarget]
    targetnoc <- df$nobs[itarget]
  }
  
  if (!is.na(fixed.oc)) {
    fixed.noc <- round(fixed.oc/100 * te)
    fixed.oc <- fixed.noc/te
    fixed.fpc <- sqrt((te - fixed.noc)/(te-1))
    fixed.oc.ucl <- fixed.fpc * te * solveucl(a=a, d=d, n=fixed.noc)
  }
      
  # plot
  if (showplot) {
    graphics::plot(100*df$pobs, log10(df$ucl), type="l", lty=1, lwd=2,
                   xlim=c(0,100), ylim=log10(c(utils::tail(ucl.dl,2)[1],min(max(ucl.dh),ymax))), 
                   xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxt="n",
                   xlab="Observer Coverage (%)", ylab="Upper Confidence Limit of Bycatch",
                   main=paste0("Upper One-Tailed ", cl, "% Confidence Limit of Bycatch Given None Observed"))
    graphics::lines(100*df$pobs, log10(ucl.dl), lty=2, lwd=2)
    graphics::lines(100*df$pobs, log10(df$ucl), lty=2, lwd=2)
    graphics::lines(100*df$pobs, log10(ucl.dh), lty=3, lwd=2)
    graphics::axis(side=2, at=log10(c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000)), 
         labels=c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000))
    if (target.ucl) {
      graphics::lines(c(0,100),log10(rep(target.ucl,2)), col=2, lwd=2, lty=4)
      graphics::points(targetoc*100, log10(df$ucl[itarget]), pch=8, col=2, cex=1.5, lwd=2)
      graphics::legend("topright", lty=c(2,1,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
                       legend=c(paste0("d=",dv[1]), paste0("d=",dv[2]), paste0("d=",dv[3]),
                                "target UCL given none obs.", "min coverage"))
    } else {
      graphics::legend("topright", lty=c(2,1,3), lwd=2, col=1, 
                       legend=c(paste0("d=",dv[1]), paste0("d=",dv[2]), paste0("d=",dv[3])))
    }
  }
  
  # print recommended minimum observer coverage
  if (!silent) {
    if (target.ucl) 
      cat(paste0("Minimum observer coverage to ensure that the upper confidence",
                 " limit of ", target.ucl, " is not exceeded when no bycatch is ",
                 "observed is ", my_ceiling(targetoc*100,3), "% (", targetnoc, 
                 " sets).\n"))
    if (!is.na(fixed.oc)) 
      cat(paste0("Upper confidence limit for bycatch given none observed in ",
                 my_ceiling(fixed.oc*100,3), "% (", fixed.noc, " sets) coverage is ",
                 my_ceiling(fixed.oc.ucl,3),".\n"))
    cat(paste0("Please review the caveats in the associated documentation.\n"))
  }
  
  # return recommended minimum observer coverage
  return(invisible(list(ucldat=df, target.ucl=ifelse(target.ucl, target.ucl, NA), 
                        targetoc=ifelse(target.ucl, targetoc, NA), 
                        targetnoc=ifelse(target.ucl, targetnoc, NA),
                        te=te, d=d, cl=cl)))
}
