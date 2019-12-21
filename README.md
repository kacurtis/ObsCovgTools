ObsCovgTools
==================================


### Description

ObsCovgTools provides tools for evaluating fishery observer coverage, 
particularly with respect to documenting and estimating rare bycatch. Current 
functionality includes (1) estimating probabilities of observing a bycatch event 
as a function of observer coverage and of any bycatch occurrring in total
effort, given mean bycatch rate, dispersion index (variance to mean ratio in the 
bycatch rate), and total fishery effort in user-defined units (e.g., trips, 
sets, or hooks); (2) estimating upper confidence limit of bycatch when none was 
observed, as a function of observer coverage, given total fishery effort and 
dispersion index; and (3) estimating bycatch estimation CV (coefficient of 
variation) as a function of observer coverage, given bycatch rate, dispersion 
index, and total fishery effort. Probability of at least one 
bycatch event occurring or being observed is based on the corresponding Poisson 
or negative binomial probability mass function. Bycatch estimation CVs are based 
on the standardized root mean square error of bycatch rate based on simulated 
bycatch in total and observed effort.


### Caveat

The current implementation of ObsCovgTools assumes observer coverage is 
representative of the fishery, and does not account for hierarchical sources of 
variance (e.g., vessel- or trip-level variation). Violating these assumptions 
will likely to result in underestimates of observer coverage required to meet 
specific objectives. More conservative estimates can be obtained by using 
higher-level units of effort (e.g., mean bycatch per trip instead of bycatch 
per set/haul, and number of trips instead of number of sets/hauls). 


### Shiny app

This package has been implemented as a Shiny app, coauthored by Howard Coleman,
which can be accessed at https://kacurtis.shinyapps.io/obscov/

