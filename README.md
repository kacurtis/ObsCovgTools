ObsCovgTools
==================================


### Description

ObsCovgTools estimates observer coverage required to achieve a specific 
objective, given mean bycatch rate, negative binomial dispersion parameter 
(characterizing variance in the bycatch rate), and total fishery effort in 
user-defined units (e.g., trips, sets, or hooks). Currently, potential 
objectives implemented are (1) a target probability of observing a bycatch 
event (given that it occurred in total effort), and (2) a target bycatch 
estimation CV (coefficient of variation). The user can also specify the desired 
probability of achieving the target CV. Probability of observing a bycatch 
event, and of any bycatch occurring in a given amount of total effort, is based 
on the corresponding Poisson or negative binomial probability distribution 
function. Bycatch estimation CVs are based on simulated mean-per-unit bycatch 
estimates with finite population correction, and exclude simulations with zero 
observed bycatch.


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

