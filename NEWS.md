# ObsCovgTools 3.3.0

* Implemented an analytical solution for CV corresponding to observer coverage, replacing sim_cv_obscov.r and plot_cv_obscov.r with a single function analogous to those for the other two management objectives. Elimination of Monte Carlo simulations provides major efficiency gains and increases precision of results.
* Calculation of the precision of bycatch estimates:
Given bycatch per unit effort $r$, dispersion index $d$, and total effort $N$, variance is given by $s^2=dr$, 
and the expected coefficient of variation of bycatch estimates for observed effort $n$, corrected for sampling from a finite population per Cochran (1973), is 
$$CV={ \sqrt{{ s^2 {{N-n} \over {N-1}}} \over n} \over r}$$
* This also resulted in removing the previous dependency on the Runuran package.


# ObsCovgTools 3.2.1-2

* Added a `NEWS.md` file to track changes to the package.
* Added a tab to shiny app to show new url (only activated for shinyapps.io deployment)


# ObsCovgTools 3.2.1

* Shifted to R 4.0.0 build and associated packages to deploy on noaa server
* Additional corrections to significant digit reporting 


# ObsCovgTools 3.2.0

* Corrected significant digit reporting 
* Removed reporting of target coverage in terms of effort units (only provided in percentage now) to avoid significant digit issues when greater number significant digits specified for total effort than available for observer coverage target
* Added cautionary note regarding need to consider uncertainty in BPUE input


# ObsCovgTools 3.0.0

* Added management objective of not exceeding a specified upper confidence limit of total bycatch  when none observed



# ObsCovgTools 2.0.0

* Added management objective of specified probability of observing bycatch if it occurs 
* Changed package name from obscov4CV to ObsCovgTools
