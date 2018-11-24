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


Disclaimer

*This software package is developed and maintained by scientists at the NOAA Fisheries Southwest Fisheries Science Center and should be considered a fundamental research communication. The recommendations and conclusions presented here are those of the authors and this software should not be construed as official communication by NMFS, NOAA, or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the highest quality, tools such as this are under constant development and are subject to change.*