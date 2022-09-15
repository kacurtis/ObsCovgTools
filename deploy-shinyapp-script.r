# deploy shiny app with rsconnect
# first restart R session for clean slate
devtools::install_github("kacurtis/ObsCovgTools")
rsconnect::deployApp(appDir="inst/shinyapp/obscov", appName="obscov")

# may need to add server info if more than one to choose from (account, server variables)