#.libPaths(c("/usr/lib64/R/shiny_library",.libPaths()))
library(ObsCovgTools)

no.msg <- paste("")

total.effort.msg <- paste("Total effort should be a positive integer",
                          "greater than 1.")
total.effort.title <- "Total effort value"

bpue.msg <- paste("BPUE should be a positive number.")
bpue.title <- "Bycatch per Unit Effort (BPUE) value"

dispersion.msg <- "Dispersion index should be a number greater than or equal to one."
dispersion.title <- "Dispersion index value"

target.ucl.msg <- "Target maximum upper confidence limit should be a number greater than or equal to zero."

fixedoc.ucl.msg <- "Percent observer coverage for which to return UCL should be a number greater than or equal to zero and less than or equal to 100."

ymax.ucl.msg <- "Upper limit for y-axis should be a number greater than zero."

check.te.inst <- function(input) {
  if (is.na(as.numeric(input)) || 
      ((as.numeric(input)) <= 1) ||
      (input != as.integer(input)) )
    total.effort.msg
  else NULL
}
check.te.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) || 
      ((as.numeric(input)) <= 1) ||
      (input != as.integer(input)) )
    no.msg
  else NULL
}

check.bpue.inst <- function(input) {
  if (is.na(as.numeric(input)) || 
      (as.numeric(input) <= 0) )
    bpue.msg
  else NULL
}
check.bpue.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) || 
      (as.numeric(input) <= 0) )
    no.msg
  else NULL
}

check.d.inst <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 1) )
    dispersion.msg
  else NULL
}
check.d.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 1) )
    no.msg
  else NULL
}

check.target.ucl.inst <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 0) )
    target.ucl.msg
  else NULL
}
check.target.ucl.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 0) )
    no.msg
  else NULL
}

check.fixedoc.ucl.inst <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 0) ||
      (as.numeric(input) > 100) )
    fixedoc.ucl.msg
  else NULL
}
check.fixedoc.ucl.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 0) ||
      (as.numeric(input) > 100) )
    no.msg
  else NULL
}

check.ymax.ucl.inst <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) <= 0) )
    ymax.ucl.msg
  else NULL
}
check.ymax.ucl.inst.lab <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) <= 0) )
    no.msg
  else NULL
}


server <- function(input, output, session) {
  
  observe({
    if (session$clientData$url_hostname=="kacurtis.shinyapps.io") {
      prependTab(inputId = "alltabs", 
                 tabPanel("New URL", 
                          tags$h3("The observer coverage simulator has a new home!"),
                          HTML("<p>Please bookmark and use the new URL:</p><p><a href='https://connect.fisheries.noaa.gov/ObsCovg/'>https://connect.fisheries.noaa.gov/ObsCovg/</a></p>")
                          
                 ),
                 select = TRUE)
    }
  })

  
  plotlabels.ppos <- reactiveValues(ppos='')
  
  ## validation functions for user inputs to ppos tab
  
  output$ppos_obscov_plot <- renderPlot({
    
    ## validation of user inputs
    
    validate(check.te.inst(input$te.ppos), 
             check.bpue.inst(input$bpue.ppos), 
             check.d.inst(input$d.ppos))
    
    plotlabels.ppos$ppos <- capture.output(
        plot_probposobs(te = as.numeric(input$te.ppos), 
                        bpue = as.numeric(input$bpue.ppos), 
                        d = as.numeric(input$d.ppos),
                        targetppos = as.numeric(input$target.ppos)))
    
  })
  
  output$ppos_obscov_plot_label <- renderText({
    
    validate(check.te.inst.lab(input$te.ppos), 
             check.bpue.inst.lab(input$bpue.ppos), 
             check.d.inst.lab(input$d.ppos))
    
    oc.ppos.out <- plot_probposobs(te = as.numeric(input$te.ppos), 
                    bpue = as.numeric(input$bpue.ppos), 
                    d = as.numeric(input$d.ppos),
                    targetppos = as.numeric(input$target.ppos),
                    silent = TRUE, as.shiny = TRUE)
    
    rec2 <- paste0(" The conditional probability of observing any bycatch if it occurs ", 
                   "(solid black line) is obtained by dividing the absolute probability",
                   " of observing any bycatch (black dashed line) by the probability ",
                   "that any bycatch occurs in the given total effort.")
    
    HTML(paste0("<p></p><ul><li>", oc.ppos.out$rec, "</li><li>", rec2, "</li><li>Please review the caveats in the About tab.</li></ul>"))
  })
  
    
  plotlabels.ucl <- reactiveValues(ucl='')
  
  ## validation functions for user inputs to ucl tab
  
  output$ucl_obscov_plot <- renderPlot({
    
    ## validation of user inputs
    
    validate(check.te.inst(input$te.ucl), 
             check.d.inst(input$d.ucl),
             check.target.ucl.inst(input$target.ucl),
             check.fixedoc.ucl.inst(input$fixedoc.ucl),
             check.ymax.ucl.inst(input$ymax.ucl))
    
    plotlabels.ucl$ucl <- capture.output(
      plot_uclnegobs(te = as.numeric(input$te.ucl), 
                     d = as.numeric(input$d.ucl),
                     cl = as.numeric(input$cl.ucl),
                     targetucl = as.numeric(input$target.ucl),
                     fixedoc = as.numeric(input$fixedoc.ucl),
                     ymax = as.numeric(input$ymax.ucl)))
    
  })
  
  output$ucl_obscov_plot_label <- renderText({
    validate(check.te.inst.lab(input$te.ucl), 
             check.d.inst.lab(input$d.ucl),
             check.target.ucl.inst.lab(input$target.ucl),
             check.fixedoc.ucl.inst.lab(input$fixedoc.ucl),
             check.ymax.ucl.inst.lab(input$ymax.ucl))
    
    oc.ucl.out <- plot_uclnegobs(te = as.numeric(input$te.ucl), 
                                 d = as.numeric(input$d.ucl),
                                 cl = as.numeric(input$cl.ucl),
                                 targetucl = as.numeric(input$target.ucl),
                                 fixedoc = as.numeric(input$fixedoc.ucl),
                                 ymax = as.numeric(input$ymax.ucl),
                                 silent = TRUE, as.shiny = TRUE)
    if (oc.ucl.out$rec=="")
      HTML(paste0("<p></p><ul><li>Please review the caveats in the About tab.</li></ul>"))
    else 
      HTML(paste0("<p></p><ul><li>", oc.ucl.out$rec, 
                  "</li><li>Please review the caveats in the About tab.</li></ul>"))
  })
  
  
  plotlabels.cv <- reactiveValues(cv = '')
  
  output$cv_obscov_plot <- renderPlot({
    
    ## validation of user inputs
    
    validate(check.te.inst(input$te.cv), 
             check.bpue.inst(input$bpue.cv), 
             check.d.inst(input$d.cv))
    
    plotlabels.cv$cv <- capture.output(
      plot_cv(te = as.numeric(input$te.cv), 
              bpue = as.numeric(input$bpue.cv), 
              d = as.numeric(input$d.cv),
              targetcv = as.numeric(input$target.cv)))
    
  })
  
  output$cv_obscov_plot_label <- renderText({
    
    validate(check.te.inst.lab(input$te.cv), 
             check.bpue.inst.lab(input$bpue.cv), 
             check.d.inst.lab(input$d.cv))
    
    oc.cv.out <- plot_cv(te = as.numeric(input$te.cv), 
                                   bpue = as.numeric(input$bpue.cv), 
                                   d = as.numeric(input$d.cv),
                                   targetcv = as.numeric(input$target.cv),
                                   silent = TRUE, as.shiny = TRUE)
    
    if (oc.cv.out$rec=="")
      HTML(paste0("<p></p><ul><li>Please review the caveats in the About tab.</li></ul>"))
    else
      HTML(paste0("<p></p><ul><li>", oc.cv.out$rec, 
                  "</li><li>Please review the caveats in the About tab.</li></ul>"))
  })
  
  
  outputOptions(output, suspendWhenHidden=FALSE)
}

