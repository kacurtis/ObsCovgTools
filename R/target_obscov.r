# Shiny notes:
# Approximately this layout: https://shiny.rstudio.com/gallery/submitbutton-demo.html
# Keyboard entry for te; default value 500; 
#   label: "Total effort (e.g., hauls) in fishery (Larger effort takes longer: ~40 s for 50K, ~3 min for 500K)"
# Keyboard entry for bpue; default value 0.05; label: "Bycatch per unit effort"
# Keyboard entry for d; default value 2; label: "Dispersion (d ~ Var/Mean; usually d < 3 for bycatch of conservation concern)"
# then an actionButton ("submit"), elicits reactive sim_obscov_cv
# sliders for targetcv; from 10 to 100; default value 30; label: "Target CV for bycatch estimates"
#   and q; from 0.5 to 0.95; default value 0.8; label: "Probability of achieving target CV"
# Interactve updating of plot of CV vs observer coverage and observer coverage recommendation based on sliders (plot_obscov_cv)
# Interactive updating of text under plot: 
#   "Minimum coverage to achieve targetcv% CV with q*100% probability is x$pobscov% (x$nobsets hauls).",
#   where targetcv and q are user-specified and x is returned from plot_obscov_cv.
# Static plot of Probability of zero observed bycatch vs observer coverage (plot_probzeroobs)
# dimmed plot while recalculating

# Add following limits on inputs:
# d limits: numeric, 1 <= d 
# bpue limits: numeric, 0 < bpue; no upper limit
# te limits: integer, te >= 10 

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
sim_obscov_cv <- function(te, bpue, d=2, ...) {
  nsim <- 1000
  obscov <- c(seq(0.001,0.005,0.001), seq(0.01,0.05,0.01), seq(0.10,1,0.05))
  simdat <- tibble::tibble(simpoc = rep(obscov, each=nsim), 
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
  abline(h=targetcv, col=2, lwd=2, lty=2)
  abline(h=100, v=100)
  # get (and add to plot) minimum required observer coverage
  targetoc <- simsum %>% dplyr::filter(qcv < targetcv/100) %>% dplyr::filter(simpoc==min(simpoc))
  points(targetoc$simpoc*100, targetoc$qcv*100, pch=8, col=2, cex=1.5, lwd=2)
  legpos <- ifelse(any(simsum$simpoc > 0.7 & simsum$q95 > 0.5), "bottomleft", "topright")
  legend(legpos, lty=c(2,0,1,0,0,0), pch=c(NA,8,1,rep(15,3)), col=c(2,2,1,"gray90","gray80","gray70"), 
         lwd=c(2,2,rep(1,4)), pt.cex=1.5,
         legend=c("target CV", "min coverage", paste(q*100,"th Percentile", sep=""),
                  ">50th Percentile",">80th Percentile",">95th Percentile"), y.intersp=1.1)
  # return recommended minimum observer coverage
  cat(paste("Minimum observer coverage to achieve ", targetcv, "% CV with ", q*100, "% probability is ", 
            targetoc$simpoc*100, "% (", targetoc$nobsets, " hauls).", sep=""))
  return(invisible(list(pobscov = targetoc$simpoc*100, nobsets=targetoc$nobsets)))
}
  

#' Plot probability of observing zero bycatch vs observer coverage
#' 
#' @param simdat Tibble output from sim_obscov_cv.
#' 
#' @return None
#' 
#' @export 
plot_probzeroobs <- function(simdat=simdat) {
  nd <- simdat %>% 
    dplyr::group_by(simpoc) %>% 
    dplyr::summarize(n=n(), ndp = sum(ob==0)/n)
  with(nd, plot(100*simpoc, 100*ndp, 
                xlim=c(0,100), ylim=c(0,round(max(100*ndp),-1)+10), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
                xlab="Observer Coverage (%)", ylab="Probability of Zero Bycatch (%)",
                main="Probability of Zero Bycatch"))
  abline(h=tail(nd$ndp,1)*100,col=2)
  legpos <- ifelse(any(nd$simpoc > 0.75 & nd$ndp > 0.8), "bottomleft", "topright")
  legend(legpos, lty=c(0,1), pch=c(1,NA), lwd=1, col=c(1,2), legend=c("in observed effort","in total effort"))
}
