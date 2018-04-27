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
# te:    total effort in fishery (should add option to use landings instead)
# bpue:  bycatch per unit effort
# d:     dispersion (>=1)
#
# Value:
# a data frame containing summarized results for each simulation
#
sim_obscov_cv <- function(te, bpue, d=2) {  
  nsim = 1000
  simpoc <- rep(c(seq(0.01,0.05,0.01), seq(0.10,0.95,0.05)), each=nsim)   # proportion observer coverage 
  simdat <- tibble::tibble(simpoc, nobsets = round(simpoc * te)) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(obsets = list(if(d==1) rpois(nobsets, bpue) else rnbinom(nobsets, size=(bpue/(d-1)), mu=bpue)), 
           ob=sum(obsets), obvar=var(obsets)) %>%  
    dplyr::filter(ob>0 & ob<nobsets) %>% 
    dplyr::mutate(xsim = ob/nobsets, sesim = sqrt(obvar/nobsets), cvsim = sesim/xsim) %>% 
    dplyr::select(-obsets) %>% 
    dplyr::ungroup()
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

plot_obscov_cv <- function(simdat=simdat, targetcv=30, q=0.5) {
  simsum <- simdat %>% dplyr::group_by(simpoc) %>% 
    dplyr::summarize( nsim=n(), meanob=mean(ob), nobsets=mean(nobsets), qcv=quantile(cvsim,q,na.rm=T), min=min(ob), max=max(ob))
  # plot 
  with(simsum, plot(100*simpoc, 100*qcv, xlab="Observer Coverage (%)", ylab=paste(q*100,"th Percentile CV (%)", sep="")))
  abline(h=targetcv, col=2)
  # get minimum required observer coverage
  targetoc <- simsum %>% dplyr::filter(qcv < targetcv/100) %>% dplyr::select(simpoc) %>% dplyr::unlist() %>% min()
  return(targetoc*100)
}