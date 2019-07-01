#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL 


## Quiets concerns of R CMD check re: the .'s that appear in pipelines
## and the "n" that is produced by dplyr::count() in a pipeline
if (getRversion() >= "2.15.1") utils::globalVariables(c("n"))


# Hidden functions to execute progress bar in Shiny
progress_init <- function(shiny.progress = FALSE) {
  if (shiny.progress) return(NULL)
  else return(utils::txtProgressBar(style=3))
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
#' \code{sim_cv_obscov} simulates bycatch estimation CVs resulting from a range 
#' of observer coverage levels, given bycatch rate, negative binomial dispersion 
#' index, and total fishery effort. 
#' 
#' \code{sim_cv_obscov} runs \code{nsim} simulations per level of observer 
#' coverage, from the larger of 0.1\% or two sets/hauls to 100\%. Simulated 
#' bycatch estimates use a simple mean-per-unit approach with finite population 
#' correction. Projections of bycatch estimation CV exclude simulations with 
#' zero observed bycatch. Since estimating variance requires at least two numbers, 
#' total effort must be at least three sets/hauls to evaluate observer coverage 
#' less than 100\%. 
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
#' @param te an integer greater than 2. Total effort in fishery (sets/hauls).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Negative binomial dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#' @param nsim a positive integer. Number of simulations to run.
#' @param ...  additional arguments for compatibility with Shiny.
#'   
#' @return A list with components:
#'   \item{simdat}{a tibble with one row per simulation and the following fields: 
#'   simulated proportion observer coverage (\code{simpoc}), number of observed sets 
#'   (\code{nobsets}), total observed bycatch (\code{ob}), variance of observed 
#'   bycatch (\code{obvar}), mean observed bycatch per unit effort (\code{xsim}), 
#'   finite population correction (\code{fpc}), standard error of observed bycatch 
#'   per unit effort (\code{sesim}), and CV of observed bycatch per unit effort 
#'   (\code{cvsim}). 
#'   For simulations with zero observed bycatch, \code{cvsim} will be NaN.}
#'   \item{bpue}{the bycatch per unit effort used.}
#'   \item{d}{the negative binomial dispersion index used.}
#'   
#' @export 
sim_cv_obscov <- function(te, bpue, d = 2, nsim = 1000, ...) {
  # check input values
  if ((ceiling(te) != floor(te)) || te<3) stop("te must be a positive integer >=3")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  if ((ceiling(nsim) != floor(nsim)) || nsim<=0) stop("nsim must be a positive integer")
  # simulate observer coverage and bycatch estimation
  if (te<20) obscov <- 1:te/te 
  else obscov <- c(seq(0.001,0.005,0.001), seq(0.01,0.05,0.01), seq(0.10,1,0.05))
  simdat <- tibble::tibble(simpoc = rep(obscov, nsim), 
                           nobsets = round(.data$simpoc * te)) %>% 
    dplyr::filter(.data$nobsets > 1) %>% 
    dplyr::mutate(ob=NA, obvar=NA)
  set.seed(Sys.time())
  pb <- progress_init(...)
  
  for (i in 1:nrow(simdat)) {
    obsets <- if(d==1) Runuran::urpois(simdat$nobsets[i], bpue) 
    else Runuran::urnbinom(simdat$nobsets[i], size=(bpue/(d-1)), prob=1/d)
    simdat$ob[i] <- sum(obsets)
    simdat$obvar[i] <- stats::var(obsets)

    if (i %% 500 == 0) {
      pb <- progbar(i, nrow(simdat), pb, ...)
    }
  }
  
  simdat <- simdat %>% 
    dplyr::mutate(xsim=.data$ob/.data$nobsets, fpc=1-.data$nobsets/te, 
                  sesim=sqrt(.data$fpc*.data$obvar/.data$nobsets), 
                  cvsim=.data$sesim/.data$xsim)
  return(list(simdat=simdat, te=te, bpue=bpue, d=d))
}


#' Plot bycatch estimation CV vs. observer coverage
#' 
#' \code{plot_cv_obscov} plots projected bycatch estimation CVs vs observer 
#' coverage for user-specified percentile (i.e., probability of achieving CV) 
#' and several default percentiles, and returns minimum observer coverage needed 
#' to achieve user-specified target CV and percentile. 
#'   
#' @param simlist list output from \code{sim_cv_obscov}.
#' @param targetcv a non-negative number less than 1. Target CV 
#'   (as a proportion). If \eqn{targetcv = 0}, no corresponding minimum observer 
#'   coverage will be highlighted.
#' @param q a positive number less than 100. Custom percentile to 
#'   be plotted, desired probability of achieving the target CV or lower. 
#'   Ignored if \eqn{targetcv = 0}.
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
#'   \item{pobscov}{minimum observer coverage in terms of percentage.} 
#'   \item{nobsets}{corresponding observed effort.}
#' @return Returned invisibly. 
#'   
#' @export 
plot_cv_obscov <- function(simlist = simlist, targetcv = 0.3, q = 80, 
                           showplot = TRUE, silent = FALSE) {
  # check input values
  if(targetcv<0 || targetcv>=1) stop("targetcv must be >= 0 and < 1")
  if(q<=0 || q>=100) stop("q must be > 0 and < 100")
  # get quantiles of bycatch estimation CVs
  simsum <- simlist$simdat %>% 
    dplyr::filter(.data$ob>0) %>% 
    dplyr::group_by(.data$simpoc) %>% 
    dplyr::summarize(nsim=n(), meanob=mean(.data$ob), nobsets=mean(.data$nobsets), 
                     qcv=stats::quantile(.data$cvsim,q/100,na.rm=T), 
                     q50=stats::quantile(.data$cvsim,0.5,na.rm=T), 
                     q80=stats::quantile(.data$cvsim,0.8,na.rm=T), 
                     q95=stats::quantile(.data$cvsim,0.95,na.rm=T), 
                     min=min(.data$ob), max=max(.data$ob))
  # get minimum required observer coverage
  if (targetcv)
    targetoc <- stats::approx(simsum$qcv, simsum$simpoc, targetcv)$y
  # plot 
  if (showplot) {
    with(simsum, graphics::plot(100*simpoc, qcv, 
                      xlim=c(0,100), ylim=c(0,1), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,1,10),
                      xlab="Observer Coverage (%)", ylab="CV of Bycatch Estimate",
                      main="CV of Bycatch Estimate vs Observer Coverage"))
    with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(1,q50,1,1),col="gray90", lty=0))
    with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(1,q80,1,1),col="gray75", lty=0))
    with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(1,q95,1,1),col="gray60", lty=0))
    with(simsum, graphics::points(100*simpoc, qcv))
    with(simsum, graphics::lines(100*simpoc, qcv))
    graphics::abline(h=1, v=100)
    legpos <- ifelse(any(simsum$simpoc > 0.7 & simsum$q95 > 0.5), "bottomleft", "topright")
    # add minimum required observer coverage
    if (targetcv) {
      graphics::abline(h=targetcv, col=2, lwd=2, lty=2)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc*100, targetcv, pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      graphics::legend(legpos, lty=c(2,0,1,0,0,0), pch=c(NA,8,1,rep(15,3)), col=c(2,2,1,"gray90","gray75","gray60"), 
             lwd=c(2,2,rep(1,4)), pt.cex=1.5, y.intersp=1.1,
             legend=c("target CV", "min coverage", paste0(q,"th percentile"),
                      ">50th percentile",">80th percentile",">95th percentile"))
    } else {
      graphics::legend(legpos, lty=c(1,0,0,0), pch=c(1,rep(15,3)), col=c(1,"gray90","gray75","gray60"), 
             lwd=c(rep(1,4)), pt.cex=1.5, y.intersp=1.1,
             legend=c(paste0(q,"th percentile"), 
                      ">50th percentile",">80th percentile",">95th percentile"))
    }
  }
  # return recommended minimum observer coverage
  if (!silent) {
    if (targetcv)  {
      cat(paste0("Minimum observer coverage to achieve CV \u2264 ", targetcv, " with ", q, "% probability is ", 
              my_ceiling(targetoc*100,2), "% (", my_ceiling(targetoc*simlist$te,2), " hauls).\n", 
              "Please review the caveats in the associated documentation.\n"))
    }
    cat(paste0("Note that results are interpolated from simulation-based projections and may vary slightly \n",
              "with repetition (see plot of sample sizes underlying projections of bycatch estimation CV).\n"))
  }
  if (targetcv) 
    return(invisible(list(pobscov = my_ceiling(targetoc*100,2), nobsets=my_ceiling(targetoc*simlist$te,2))))
}


#' Get probability of zero bycatch given effort, bycatch rate, and dispersion
#' 
#' \code{get_probzero} returns probability of zero bycatch in a specified number 
#' of sets/hauls, given bycatch per unit effort and negative binomial dispersion 
#' index. 
#' 
#' @param n a vector of positive integers. Observed effort levels (in terms of 
#'   sets/hauls) for which to calculate probability of zero bycatch.
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Negative binomial dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of set-level bycatch, so \eqn{d = 1} corresponds to Poisson-distributed 
#'   bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#'   
#' @details
#' Calculated from the probability density at zero of the corresponding Poisson
#' (\eqn{d = 1}) or negative binomial (\eqn{d < 1}) distribution.
#' 
#' \strong{Caveat:} \code{get_probzero} assumes representative observer coverage 
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
get_probzero <- function(n, bpue, d) {
  # check input values
  if (any((ceiling(n) != floor(n)) | n<1)) stop("n must be a vector of positive integers")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  # calculate probability of observing zero bycatch in n sets
  pz <- if(d==1) stats::ppois(0, bpue)^n 
  else stats::pnbinom(0, size=(bpue/(d-1)), prob=1/d)^n
  return(invisible(pz))
}


#' Plot sample sizes underlying projections of bycatch estimation CV
#' 
#' \code{plot_cvsim_samplesize} plots (1) the number of simulations with positive
#' observed bycatch underlying each projected bycatch estimation CV vs corresponding
#' observer coverage level, which is inversely proprtional to (2) the probability 
#' of no bycatch being observed (based on effort and the probability density at zero,
#' given bycatch rate and negative binomial dispersion). The latter includes the 
#' probability of no bycatch occurring in the given total effort. The smaller the 
#' sample size, the less precise the projection (the more it will vary among repeated 
#' runs of the simulator). Small sample sizes can be countered by increasing the 
#' number of simulations in \code{sim_cv_obscov}. 
#' 
#' Note that the projected probability of observing zero bycatch assumes 
#' representative observer coverage and no hierarchical sources of variance 
#' (e.g., vessel- or trip-level variation). Violating these assumptions will
#' likely result negatively biased projections of the probability of observing zero 
#' bycatch at a given level of observer coverage. 
#' 
#' @param simlist List output from sim_cv_obscov.
#' 
#' @return None
#' 
#' @export 
plot_cvsim_samplesize <- function(simlist = simlist) {
  s <- simlist$simdat %>% 
    dplyr::filter(.data$ob>0) %>% 
    dplyr::group_by(.data$simpoc, .data$nobsets) %>% 
    dplyr::summarize(npos=n())
  pz <- get_probzero(s$nobsets, simlist$bpue, simlist$d)
  omar <- graphics::par()$mar
  graphics::par(mar = c(5.1,4.1,4.1,4.1))
  with(s, graphics::plot(100*simpoc, npos, pch=22,
               xlim=c(0,100), ylim=c(0,round(max(npos),-1)+10), xaxs="i", yaxs="i",
               xaxp=c(0,100,10), yaxp=c(0,1000,10),
               xlab="Observer Coverage (%)", ylab="Simulations with Positive Bycatch",
               main="Sample Sizes Underlying Projections of Bycatch Estimation CV"))
  graphics::par(new=T)
  graphics::plot(100*s$simpoc, 100*pz, type="l", lwd=2, xaxs="i", yaxs="i", xlim=c(0,100), ylim=c(0,100),
       axes=F, xlab=NA, ylab=NA, col=2)
  graphics::axis(side = 4, col=2, col.axis=2)
  graphics::mtext(side = 4, line = 3, "Probability of Zero Bycatch (%)", col=2)
  graphics::abline(h=100*utils::tail(pz,1), lty=3, lwd=2, col=2)
  legpos <- ifelse(any(s$simpoc > 0.7 & pz < 0.2), "right", "bottomright")
  graphics::legend(legpos, lty=c(1,3), col=2, lwd=2, text.col=2, bty="n", legend=c("in observed effort","in total effort"))
  graphics::par(mar=omar)
}  


#' Plot probability of positive observed bycatch vs observer coverage
#' 
#' \code{plot_probposobs} plots probability of observing at least one bycatch 
#'   event vs observer coverage, given total effort in sets/hauls, bycatch per 
#'   unit effort, and negative binomial dispersion index. 
#'   
#' @param te an integer greater than 1. Total effort in fishery (sets/hauls).
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Negative binomial dispersion 
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
#' \code{get_probzero}. 
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
#'   \item{pobscov}{minimum observer coverage in terms of percentage.} 
#'   \item{nobsets}{corresponding observed effort.}
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
  if (te<1000) obscov <- 1:te/te 
  else obscov <- seq(0.001,1,0.001)
  oc <- tibble::tibble(obscov = obscov,
               nobsets = round(.data$obscov * te)) %>% 
    dplyr::filter(.data$nobsets>0) %>% as.data.frame()
  oc$pp <- 1-get_probzero(oc$nobsets, bpue, d)   # probability of positive observed bycatch
  ppt <- utils::tail(oc$pp,1)   # probability of positive bycatch in total effort
  targetoc <- log(1-(target.ppos/100)*ppt)/log(get_probzero(1,bpue,d))/te   # target observer coverage
  # plot
  if (showplot) {
    graphics::plot(100*oc$obscov, 100*(oc$pp/ppt), type="l", lty=1, lwd=2,
         xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
         xlab="Observer Coverage (%)", ylab="Probability of Positive Bycatch (%)",
         main="Probability of Positive Bycatch")
    graphics::abline(h=100*ppt,lwd=2, lty=3)
    graphics::lines(100*oc$obscov, 100*oc$pp, lwd=2, lty=2)
    graphics::lines(100*oc$obscov, 100*(oc$pp/ppt), lwd=2)
    legpos <- ifelse(any(oc$obscov > 0.6 & oc$pp < 0.3 ), "topleft", "bottomright")
    if (target.ppos) {
      graphics::abline(h=target.ppos, col=2, lwd=2, lty=4)
      graphics::par(xpd=TRUE)
      graphics::points(targetoc*100, target.ppos, pch=8, col=2, cex=1.5, lwd=2)
      graphics::par(xpd=FALSE)
      graphics::legend(legpos, lty=c(1,2,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
             legend=c("in observed effort if total bycatch > 0", "in observed effort",
                      "in total effort", "in target coverage if total bycatch > 0", "min coverage"))
    } else {
      graphics::legend(legpos, lty=c(1,2,3), lwd=2, col=1, 
             legend=c("in observed effort if total bycatch > 0","in observed effort","in total effort"))
    }
  }
  # return recommended minimum observer coverage
  if (target.ppos) {
    if (!silent) {
      cat(paste0("Minimum observer coverage to achieve at least ", target.ppos, 
                "% probability of observing \nbycatch when total bycatch is positive is ", 
                my_ceiling(targetoc*100,2), "% (", my_ceiling(targetoc*te,2), " sets). The probability\n",
                "that any bycatch occurs in the given total effort is ", signif(100*ppt,2), "%.\n",
                "Please review the caveats in the associated documentation.\n"))
    }
    return(invisible(list(pobscov=my_ceiling(targetoc*100,2), 
                          nobsets=my_ceiling(targetoc*te,2),
                          ppos.te=signif(100*ppt,2))))
  }
}
