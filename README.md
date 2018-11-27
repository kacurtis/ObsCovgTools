ObsCovgTools
==============================================================

Description

ObsCovgTools estimates observer coverage required to achieve a specific
objective, given bycatch rate, negative binomial dispersion parameter, and
total fishery effort. Currently, potential objectives implemented are 
(1) target probability of observing a bycatch event (given that it occurred in 
total effort), and (2) target bycatch estimation CV (coefficient of variation). 
The user can also specify the desired probability of achieving the target CV. 
Probability of observing a bycatch event is based on the corresponding Poisson
or negative binomial probability distribution function. Bycatch estimation CV 
is based on simulated mean-per-unit bycatch estimates with finite population 
correction, which assume representative observer coverage. CV estimates exclude 
simulations with zero observed bycatch. 