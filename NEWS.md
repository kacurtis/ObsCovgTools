# ObsCovgTools 3.3.0

* Implemented an analytical solution for CV corresponding to observer coverage, combining sim_cv_obscov.r and plot_cv_obscov.r into one function analogous to those for the other two management objectives. Elimination of Monte Carlo simulations provides major efficiency gains and increases precision of result. 


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
