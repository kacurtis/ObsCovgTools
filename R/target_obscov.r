# Shiny notes:
# Approximately this layout: https://shiny.rstudio.com/gallery/submitbutton-demo.html
# Keyboard entry for te ("Total effort (e.g.., sets or hauls) in fishery", default value 500)
# Keyboard entry for bpue ("Bycatch per unit effort", default value 0.05)
# Keyboard entry for d ("dispersion", default value 2)
# then an actionButton ("submit"), elicits reactive sim_obscov_cv
# sliders for targetcv ("Target CV", 10 to 100, default at 30) and q ("Probability of achieving target CV", 0.5 to 0.95, default at 0.5)
# interactve updating of plot and observer coverage recommendation based on sliders (plot_obscov_cv)
# text under plot with result from plot_obscov_cv: 
# dimmed plot while recalculating


# Simulate CV response to observer coverage
#
# Arguments: 
# te:    total effort in fishery
# bpue:  bycatch per unit effort
# d:     dispersion (>=1)
#
# Value:
# a data frame containing summarized results for each simulation
#
sim_obscov_cv <- function(te, bpue, d=2) {
  nsim <- 1000
  obscov <- c(seq(0.001,0.005,0.001), seq(0.01,0.05,0.01), seq(0.10,1,0.05))
  simdat <- tibble::tibble(simpoc = rep(obscov, each=nsim), 
                           nobsets = round(simpoc * te)) %>% 
    dplyr::filter(nobsets > 1) %>% dplyr::mutate(ob=NA, obvar=NA)   # need nobsets > 1 to calculate a variance
  
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


# Plot CV response to observer coverage
#
# Arguments:
# cv:    target CV (as proportion, i.e., 0 < cv < 1)
# q:     desired probability (as a proportion) of achieving at least target CV or lower
#
# Value:
# a scalar estimate of minimum observer coverage (percentage) required to meet target CV

plot_obscov_cv <- function(simdat=simdat, targetcv=30, q=0.8) {
  simsum <- simdat %>% 
    dplyr::filter(ob>0) %>% 
    dplyr::group_by(simpoc) %>% 
    dplyr::summarize(nsim=n(), meanob=mean(ob), nobsets=mean(nobsets), qcv=quantile(cvsim,q,na.rm=T), 
                     q50=quantile(cvsim,0.5,na.rm=T), q80=quantile(cvsim,0.8,na.rm=T),
                     q90=quantile(cvsim,0.9,na.rm=T), q95=quantile(cvsim,0.95,na.rm=T),
                     min=min(ob), max=max(ob))
  
  # plot CV vs. observer coverage
  opar <- par()
  par(mfrow=c(2,1))
  layout(matrix(1:2,2,1))
  with(simsum, plot(100*simpoc, 100*qcv, ylim=c(0,100), 
                    xlab="Observer Coverage (%)", ylab=paste(q*100,"th Percentile CV (%)", sep=""),
                    main="Bycatch Estimate CV with Observer Coverage"))
  with(simsum, lines(100*simpoc, 100*qcv))
  abline(h=targetcv, col=2)
  # get (and plot) minimum required observer coverage
  targetoc <- simsum %>% dplyr::filter(qcv < targetcv/100) %>% dplyr::filter(simpoc==min(simpoc))
  points(targetoc$simpoc*100, targetoc$qcv*100, pch=8, lwd=2, col=2)
  legend("topright", lty=c(1,0,0), pch=c(NA,8,NA), lwd=2, col=2, legend=c("target CV", "minimum observer","coverage needed"))
  
  # plot probability of zero (observed and total) bycatch vs observer coverage
  nd <- simdat %>% 
    dplyr::group_by(simpoc) %>% 
    dplyr::summarize(n=n(), ndp = sum(ob==0)/n)
  with(nd, plot(100*simpoc, 100*ndp, ylim=c(0,100), yaxs="i",
                xlab="Observer Coverage (%)", ylab="Probability of Zero Bycatch (%)",
                main="Annual Probability of Zero Bycatch"))
  abline(h=tail(nd$ndp,1)*100,col=2)
  legend("topright", lty=c(0,1), pch=c(1,NA), lwd=1, col=c(1,2), legend=c("in observed effort","in total effort"))
  
  par <- opar
  return(targetoc$simpoc*100)
}