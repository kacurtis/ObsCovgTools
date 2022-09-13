# deploy shiny app with rsconnect
# first restart R session for clean slate
devtools::install_github("kacurtis/ObsCovgTools")
rsconnect::deployApp(appDir="inst/shinyapp/obscov", appName="obscov")

# add two separate servers (account, server variables)