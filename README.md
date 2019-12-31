ObsCovgTools
==================================


### Description

ObsCovgTools provides tools for evaluating fishery observer coverage, 
particularly with respect to documenting and estimating rare bycatch. Current 
functionality includes evaluating observer coverage in terms of (1) probabilities 
of observing a bycatch event and of any bycatch occurring in total effort, given 
mean bycatch rate, dispersion index (variance to mean ratio in the bycatch rate), 
and total fishery effort; (2) upper confidence limit of bycatch when none was 
observed, given total fishery effort and dispersion index; and (3) bycatch 
estimation CV (coefficient of variation), given bycatch rate, dispersion 
index, and total fishery effort. Estimates in all cases are based directly on or
simulated from the corresponding Poisson or negative binomial probability 
distribution.


### Caveat

The current implementation of ObsCovgTools assumes observer coverage is 
representative of the fishery, and does not account for hierarchical sources of 
variance (e.g., vessel- or trip-level variation). Violating these assumptions 
may result in negatively biased projections of observer coverage required to meet 
specific objectives. Unless hierarchical sources of variance can be ruled out as 
potentially important, using higher-level units of effort is advised (e.g., mean 
bycatch per trip and number of trips, instead of mean bycatch per set and 
number of sets). 


### Shiny app

This package has been implemented as a Shiny web application, coauthored by 
Howard Coleman, which can be accessed at https://kacurtis.shinyapps.io/obscov/.

