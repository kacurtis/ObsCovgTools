#' @importFrom magrittr %>%
#' @importFrom graphics abline legend points
#' @importFrom stats quantile var
#' @importFrom utils tail
NULL 


# Hidden function to execute progress bar
progbar = function(it, total, shiny.progress=FALSE) {
  if (shiny.progress) {
    shiny::incProgress(500 / total)
  } else {
    svMisc::progress(it/total*100)
  }
}


#' Simulate CV response to observer coverage
#'
#' \code{sim_obscov_cv} simulates bycatch estimation CVs resulting from a range 
#' of observer coverage levels, given bycatch rate, negative binomial dispersion 
#' parameter, and total fishery effort. The function runs 1000 simulations per
#' level of observer coverage, for observer coverage levels ranging from 0.1%
#' or two sets/hauls, whichever is greater, to 100%. Simulated bycatch estimates
#' use a simple mean-per-unit approach with finite population correction.
#' WARNING: Calls specifying large (>10K sets/hauls) total effort may take 
#' several minutes to simulate. 
#' 
#' @param te Integer scalar greater than 10. Total effort in fishery (sets).
#' @param bpue Numeric greater than zero. Bycatch per unit effort.
#' @param d Numeric >= 1. Negative binomial dispersion parameter. The dispersion
#'   parameter corresponds to the variance-to-mean ratio of set-level bycatch, 
#'   so d=1 corresponds to Poisson-distributed bycatch, and d>1 corresponds to
#'   overdispersed bycatch.
#' @param nsim Integer scalar >= 1. Number of simulations.
#' @param ...  Additional arguments for compatibility with Shiny.
#'   
#' @return A tibble with one row per simulation and the following fields: 
#'   simulated percent observer coverage (simpoc), number of observed sets 
#'   (nobsets), total observed bycatch (ob), variance of observed bycatch 
#'   (obvar), mean observed bycatch per unit effort (xsim), finite population 
#'   correction (fpc), standard error of observed bycatch per unit effort 
#'   (sesim), and CV of observed bycatch per unit effort (cvsim).
#'    
#'   For simulations with zero observed bycatch, cvsim will be NaN.
#'   
#' @export 
sim_obscov_cv <- function(te, bpue, d=2, nsim=1000, ...) {  
  obscov <- c(seq(0.001,0.005,0.001), seq(0.01,0.05,0.01), seq(0.10,1,0.05))
  simdat <- tibble::tibble(simpoc = rep(obscov, nsim), 
                           nobsets = round(simpoc * te)) %>% 
    dplyr::filter(nobsets > 1) %>% dplyr::mutate(ob=NA, obvar=NA)
  
  for (i in 1:nrow(simdat)) {
    obsets <- if(d==1) Runuran::urpois(simdat$nobsets[i], bpue) 
    else Runuran::urnbinom(simdat$nobsets[i], size=(bpue/(d-1)), prob=1/d)
    simdat$ob[i] <- sum(obsets)
    simdat$obvar[i] <- var(obsets)
    
    if (i %% 500 == 0) progbar(i, nrow(simdat), ...)
  }
  
  simdat <- simdat %>% 
    dplyr::mutate(xsim=ob/nobsets, fpc=1-nobsets/te, sesim=sqrt(fpc*obvar/nobsets), cvsim=sesim/xsim)
  return(simdat)
}


#' Plot CV vs. observer coverage
#' 
#' \code{plot_obscov_cv} plots CV of bycatch estimates vs observer coverage for
#'   user-specified quantile (i.e., probability of achieving CV) and several 
#'   default quantiles, and prints minimum observer coverage needed to achieve 
#'   user targets. 
#'
#' @param simdat Tibble output from sim_obscov_cv.
#' @param targetcv Numeric, 0 < targetcv <=100. Target CV (as percentage). 
#'    If 0, no corresponding minimum observer coverage will be highlighted.
#' @param q Numeric, 0 < q <=0.95. Desired probability (as a proportion) of 
#'   achieving at least target CV or lower.
#'   
#' @return A list with minimum observer coverage in terms of percentage ($pobscov) 
#'   and effort ($nobsets) corresponding to user specifications.
#' @return Returned invisibly. 
#'   
#' @export 
plot_obscov_cv <- function(simdat=simdat, targetcv=30, q=0.8) {
  # get quantiles of bycatch estimation CVs
  simsum <- simdat %>% 
    dplyr::filter(ob>0) %>% 
    dplyr::group_by(simpoc) %>% 
    dplyr::summarize(nsim=n(), meanob=mean(ob), nobsets=mean(nobsets), qcv=quantile(cvsim,q,na.rm=T), 
                     q50=quantile(cvsim,0.5,na.rm=T), q80=quantile(cvsim,0.8,na.rm=T), 
                     q95=quantile(cvsim,0.95,na.rm=T), min=min(ob), max=max(ob))
  # plot 
  with(simsum, plot(100*simpoc, 100*qcv, 
                    xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
                    xlab="Observer Coverage (%)", ylab="CV of Bycatch Estimate (%)",
                    main="CV of Bycatch Estimate vs Observer Coverage"))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q50,100,100),col="gray90", lty=0))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q80,100,100),col="gray80", lty=0))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q95,100,100),col="gray70", lty=0))
  with(simsum, points(100*simpoc, 100*qcv))
  with(simsum, lines(100*simpoc, 100*qcv))
  abline(h=100, v=100)
  legpos <- ifelse(any(simsum$simpoc > 0.7 & simsum$q95 > 0.5), "bottomleft", "topright")
  # get (and add to plot) minimum required observer coverage
  if (targetcv) {
    abline(h=targetcv, col=2, lwd=2, lty=2)
    targetoc <- simsum %>% dplyr::filter(qcv <= targetcv/100) %>% dplyr::filter(simpoc==min(simpoc))
    points(targetoc$simpoc*100, targetoc$qcv*100, pch=8, col=2, cex=1.5, lwd=2)
    legend(legpos, lty=c(2,0,1,0,0,0), pch=c(NA,8,1,rep(15,3)), col=c(2,2,1,"gray90","gray80","gray70"), 
           lwd=c(2,2,rep(1,4)), pt.cex=1.5,
           legend=c("target CV", "min coverage", paste(q*100,"th Percentile", sep=""),
                    ">50th Percentile",">80th Percentile",">95th Percentile"), y.intersp=1.1)
  } else {
    legend(legpos, lty=c(1,0,0,0), pch=c(1,rep(15,3)), col=c(1,"gray90","gray80","gray70"), 
           lwd=c(rep(1,4)), pt.cex=1.5,
           legend=c(paste(q*100,"th Percentile", sep=""),
                    ">50th Percentile",">80th Percentile",">95th Percentile"), y.intersp=1.1)
  }
  # return recommended minimum observer coverage
  if (targetcv)
    cat(paste("Minimum observer coverage to achieve ", targetcv, "% CV with ", q*100, "% probability is ", 
            targetoc$simpoc*100, "% (", targetoc$nobsets, " hauls).\n", sep=""))
  cat("Note that results are simulation-based and may vary slightly with repetition.\n")
  if (targetcv) 
    return(invisible(list(pobscov = targetoc$simpoc*100, nobsets=targetoc$nobsets)))
}


#' Get probability of zero bycatch given effort, bycatch rate, and dispersion
#' 
#' \code{get_probzero} Returns probability of zero bycatch in a given number of
#' sets, calculated from the probability density of the Poisson or negative 
#' binomial distribution given bycatch per unit effort and negative binomial
#' dispersion parameter.
#' 
#' @param n Integer vector. Observed effort levels (sets) for which to calculate
#'  probability of zero bycatch.
#' @param bpue Numeric greater than zero. Bycatch per unit effort.
#' @param d Numeric >= 1. Negative binomial dispersion parameter. The dispersion
#'   parameter corresponds to the variance-to-mean ratio of set-level bycatch, 
#'   so d=1 corresponds to Poisson-distributed bycatch, and d>1 corresponds to
#'   overdispersed bycatch.
#'   
#' @return Vector of same length as n with probabilities of zero bycatch. 
#' @return Returned invisibly
#' 
#' @export
get_probzero <- function(n, bpue, d) {
  pz <- if(d==1) ppois(0, bpue)^n 
  else pnbinom(0, size=(bpue/(d-1)), prob=1/d)^n
  return(invisible(pz))
}


#' Plot sample size for CV estimates vs observer coverage
#' 
#' \code{plot_samplesize} Plots sample size (simulations with positive
#' observed bycatch) vs observer coverage level, along with probability of 
#' observing zero bycatch based on effort and the probability density at zero
#' given bycatch rate and negative binomial dispersion. 
#' 
#' @param simdat Tibble output from sim_obscov_cv.
#' @param bpue Numeric greater than zero. Bycatch per unit effort.
#' @param d Numeric >= 1. Negative binomial dispersion parameter. The dispersion
#'   parameter corresponds to the variance-to-mean ratio of set-level bycatch, 
#'   so d=1 corresponds to Poisson-distributed bycatch, and d>1 corresponds to
#'   overdispersed bycatch.
#' 
#' @return None
#' 
#' @export 
plot_samplesize <- function(simdat=simdat, bpue, d) {
  s <- simdat %>% 
    dplyr::filter(ob>0) %>% 
    dplyr::group_by(simpoc, nobsets) %>% 
    dplyr::summarize(npos=n())
  pz <- get_probzero(s$nobsets, bpue, d)
  omar <- par()$mar
  par(mar = c(4.1,4.1,3,4.1))
  with(s, plot(100*simpoc, npos, pch=22,
               xlim=c(0,100), ylim=c(0,round(max(npos),-1)+10), xaxs="i", yaxs="i",
               xaxp=c(0,100,10), yaxp=c(0,1000,10),
               xlab="Observer Coverage (%)", ylab="Simulations with Positive Bycatch",
               main="Sample Size for CV Estimates"))
  par(new=T)
  plot(100*s$simpoc, 100*pz, type="l", lwd=2, xaxs="i", yaxs="i", xlim=c(0,100), ylim=c(0,100),
       axes=F, xlab=NA, ylab=NA, col=2)
  axis(side = 4, col=2, col.axis=2)
  mtext(side = 4, line = 3, "Probability of Zero Bycatch (%)", col=2)
  abline(h=100*tail(pz,1), lty=3, lwd=2, col=2)
  legpos <- ifelse(any(s$simpoc > 0.7 & pz < 0.2), "right", "bottomright")
  legend(legpos, lty=c(1,3), col=2, lwd=2, text.col=2, bty="n", legend=c("in observed effort","in total effort"))
  par(mar=omar)
}  


#' Plot probability of positive observed bycatch vs observer coverage ### add options to this and plotcv to solve or not
#' 
#' @param te Integer scalar greater than 10. Total effort in fishery (sets).
#' @param bpue Numeric greater than zero. Bycatch per unit effort.
#' @param d Numeric >= 1. Negative binomial dispersion parameter. The dispersion
#'   parameter corresponds to the variance-to-mean ratio of set-level bycatch, 
#'   so d=1 corresponds to Poisson-distributed bycatch, and d>1 corresponds to
#'   overdispersed bycatch.
#' @param target.ppos Numeric, 0 < target.ppos <=100. Target probability of
#'   positive observed bycatch (as percentage), given positive bycatch in total 
#'   effort. If 0, no corresponding minimum observer coverage will be highlighted.
#' 
#' @return A list with minimum observer coverage in terms of percentage ($pobscov) 
#'   and effort ($nobsets) corresponding to user specifications.
#' @return Returned invisibly. 
#' 
#' @export 
plot_probposobs <- function(te, bpue, d, target.ppos=80) {
  # percent probablity of positive observed bycatch
  oc <- tibble::tibble(obscov = c(seq(0.001,0.005,0.001), seq(0.01,1,0.01)),
               nobsets = round(obscov * te)) %>% 
    dplyr::filter(nobsets>0) %>% as.data.frame()
  oc$pp <- 1-get_probzero(oc$nobsets, bpue, d)   # probability of positive observed bycatch
  ppt <- tail(oc$pp,1)   # probability of positive bycatch in total effort
  plot(100*oc$obscov, 100*(oc$pp/ppt), type="l", lty=1, lwd=2,
       xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
       xlab="Observer Coverage (%)", ylab="Probability of Positive Bycatch (%)",
       main="Probability of Positive Bycatch")
  abline(h=100*ppt,lwd=2, lty=3)
  lines(100*oc$obscov, 100*oc$pp, lwd=2, lty=2)
  lines(100*oc$obscov, 100*(oc$pp/ppt), lwd=2)
  legpos <- ifelse(any(oc$obscov > 0.6 & oc$pp < 0.3 ), "topleft", "bottomright")
  if (target.ppos) {
    abline(h=target.ppos, col=2, lwd=2, lty=4)
    itargetoc <- which.max(100*oc$pp/ppt >= target.ppos)
    points(oc$obscov[itargetoc]*100, (oc$pp/ppt)[itargetoc]*100, pch=8, col=2, cex=1.5, lwd=2)
    legend(legpos, lty=c(1,2,3,4,NA), pch=c(NA,NA,NA,NA,8), lwd=2, col=c(1,1,1,2,2), pt.cex=1.5, 
           legend=c("in observed effort if total bycatch > 0", "in observed effort",
                    "in total effort", "in target observer coverage", "min coverage"))
  } else {
    legend(legpos, lty=c(1,2,3), lwd=2, col=1, 
           legend=c("in observed effort if total bycatch > 0","in observed effort","in total effort"))
  }
  # return recommended minimum observer coverage
  if (target.ppos) {
    cat(paste("Minimum observer coverage to achieve ", target.ppos, "% probability ",
              "of observing bycatch when\ntotal bycatch is positive is ", 
              oc$obscov[itargetoc]*100, "% (", oc$nobsets[itargetoc], " hauls).\n", 
              sep=""))
    return(invisible(list(pobscov=oc$obscov[itargetoc]*100, nobsets=oc$nobsets[itargetoc])))
  }
}