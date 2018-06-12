# Shiny notes:
# Approximately this layout: https://shiny.rstudio.com/gallery/submitbutton-demo.html
# Keyboard entry for te ("Total effort (e.g.., sets or hauls) in fishery", default value 500)
# Keyboard entry for bpue ("Bycatch per unit effort", default value 0.05)
# Keyboard entry for d ("dispersion", default value 2)
# then an actionButton ("submit"), elicits reactive sim_obscov_cv
# sliders for targetcv ("Target CV", 10 to 100, default at 30) and q ("Probability of achieving target CV", 0.5 to 0.95, default at 0.5)
# interactve updating of plot of CV vs observer coverage and observer coverage recommendation based on sliders (plot_obscov_cv)
# Static plot of Probability of zero observed bycatch vs observer coverage (plot_probzeroobs)
# text under plot with result from plot_obscov_cv: 
# dimmed plot while recalculating

# Add following limits on inputs:
# d limits: numeric, 1 <= d (Note: var(x)/mean(x); most bycatch of conservation concern will have d between 1 and 3) 
# bpue limits: numeric, 0 < bpue; no upper limit (Note: e.g., 0.01 if one interaction in every hundred sets/hauls)
# te limits: integer, te >= 10 (Note: large te may take several minutes, e.g., 40 s for 50K, 3 min for 500K)


#' @importFrom magrittr %>%
#' @importFrom graphics abline legend points
#' @importFrom stats quantile var
#' @importFrom utils tail
NULL 


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
sim_obscov_cv <- function(te, bpue, d=2) {
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
  }
  
  simdat <- simdat %>% 
    dplyr::mutate(xsim=ob/nobsets, fpc=1-nobsets/te, sesim=sqrt(fpc*obvar/nobsets), cvsim=sesim/xsim)
  return(simdat)
}


#' Plot CV vs. observer coverage
#'
#' @param simdat Tibble output from sim_obscov_cv.
#' @param targetcv Numeric, 0 < targetcv <=100. Target CV (as percentage).
#' @param q Numeric, 0 < q <=0.95. Desired probability (as a proportion) of 
#'   achieving at least target CV or lower.
#'   
#' @return Numeric scalar. Estimate of minimum observer coverage (percentage) 
#'   required to meet target CV.
#'   
#' @export 
plot_obscov_cv <- function(simdat=simdat, targetcv=30, q=0.8) {
  # get quantiles of bycatch estimation CVs
  simsum <- simdat %>% 
    dplyr::filter(ob>0) %>% 
    dplyr::group_by(simpoc) %>% 
    dplyr::summarize(nsim=n(), meanob=mean(ob), nobsets=mean(nobsets), qcv=quantile(cvsim,q,na.rm=T), 
                     q50=quantile(cvsim,0.5,na.rm=T), q80=quantile(cvsim,0.8,na.rm=T), 
                     q90=quantile(cvsim,0.9,na.rm=T), q95=quantile(cvsim,0.95,na.rm=T),
                     min=min(ob), max=max(ob))
  # plot 
  with(simsum, plot(100*simpoc, 100*qcv, 
                    xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i", xaxp=c(0,100,10), yaxp=c(0,100,10),
                    xlab="Observer Coverage (%)", ylab="CV of Bycatch Estimate (%)",
                    main="CV of Bycatch Estimate vs Observer Coverage"))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q50,100,100),col="gray90", lty=0))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q80,100,100),col="gray80", lty=0))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q90,100,100),col="gray70", lty=0))
  with(simsum, polygon(c(100*simsum$simpoc[1],100*simsum$simpoc,100,0), c(100,100*q95,100,100),col="gray60", lty=0))
  with(simsum, points(100*simpoc, 100*qcv))
  with(simsum, lines(100*simpoc, 100*qcv))
  abline(h=targetcv, col=2, lwd=2, lty=2)
  # get (and add to plot) minimum required observer coverage
  targetoc <- simsum %>% dplyr::filter(qcv < targetcv/100) %>% dplyr::filter(simpoc==min(simpoc))
  points(targetoc$simpoc*100, targetoc$qcv*100, pch=8, col=2, cex=1.5, lwd=2)
  legend("topright", lty=c(2,0,0,0,1,0,0,0,0), pch=c(NA,8,NA,NA,1,rep(15,4)), col=c(2,2,2,2,1,"gray90","gray80","gray70","gray60"), 
         lwd=c(2,2,rep(1,7)), pt.cex=1.5,
         legend=c("target CV", "min coverage", paste("to achieve ",targetcv,"% CV",sep=""), paste("with ",q*100,"% probability",sep=""),
                  paste(q*100,"th Percentile", sep=""),
                  ">50th Percentile",">80th Percentile",">90th Percentile",">95th Percentile"), y.intersp=1.1)
  # return recommended minimum observer coverage
  return(targetoc$simpoc*100)
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
                main="Annual Probability of Zero Bycatch"))
  abline(h=tail(nd$ndp,1)*100,col=2)
  legend("topright", lty=c(0,1), pch=c(1,NA), lwd=1, col=c(1,2), legend=c("in observed effort","in total effort"))
}
