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

processing.title <- HTML(paste0("<center>Computing - ",
                                "please wait...</center>"))
processing.msg <- HTML(paste0("<center>Reload your browser to cancel processing.", 
                              "<br>(To the left of the address bar, press Reload and then Stop Loading).",
                              "<center>"))

progress.caption <- "Simulation progress"

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
    
    if (as.logical(input$target.ppos)) {
      rec <- paste0("The probability that any bycatch occurs in the given total ",
                    "effort (horizontal black dotted line) is ",oc.ppos.out$ppos.te, 
                    "%. Minimum observer coverage to achieve at least ", 
                    input$target.ppos, "% probability of observing bycatch when",
                    " total bycatch is positive is ", oc.ppos.out$pobs, "%.")
    } else rec <- paste0("The probability that any bycatch occurs in the given total ",
                         "effort (horizontal black dotted line) is ",oc.ppos.out$ppos.te, 
                         "%.")
    HTML(paste0(rec,
                " The conditional probability of observing any bycatch if it occurs ", 
                "(solid black line) is obtained by dividing the absolute probability",
                " of observing any bycatch (black dashed line) by the probability ",
                "that any bycatch occurs in the given total effort. Please review",
                " the caveats in the About tab."))
  })
    
  plotlabels.ucl <- reactiveValues(ucl='')
  
  ## validation functions for user inputs to ppos tab
  
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
    
    if (as.logical(as.numeric(input$target.ucl))) {
      rec1 <- paste0("Minimum observer coverage to ensure that the upper confidence",
                    " limit of ", input$target.ucl, " is not exceeded when no bycatch is ",
                    " observed is ", oc.ucl.out$targetoc, "%.\n")
    } else { rec1 <- "" }
    if (as.logical(as.numeric(input$fixedoc.ucl))) {
      rec2 <- paste0("Upper confidence limit for bycatch given none observed in ",
                     oc.ucl.out$fixedoc, "% (", oc.ucl.out$fixednoc, " trips or sets)",
                     " coverage is ", oc.ucl.out$fixedoc.ucl, ".\n")
    } else { rec2 <- "" }
    HTML(paste0(rec1, rec2, "Please review the caveats in the About tab."))
  })
  
  plotlabels.cv <- reactiveValues(cv = '')
  
  simlist <- reactive({
    if (input$submit.cv > 0) {
      isolate({

        ## validation of user inputs
        
        if (is.na(as.numeric(input$te)) ||
              ((te <- as.numeric(input$te)) < 2) ||
              (te != as.integer(te))) {
          showModal(modalDialog(title=total.effort.title,
                                 total.effort.msg,
                                 easyClose=TRUE))
          return()
        }

        if (is.na(as.numeric(input$bpue)) ||
              ((bpue <- as.numeric(input$bpue)) <= 0)) {
          showModal(modalDialog(title=bpue.title,
                                 bpue.msg,
                                 easyClose=TRUE))
          return()
        }

        if (is.na(as.numeric(input$d)) ||
              ((d <- as.numeric(input$d)) < 1)) {
          showModal(modalDialog(title=dispersion.title,
                                 dispersion.msg,
                                 easyClose=TRUE))
          return()
        }

        showModal(modalDialog(
          title=processing.title,
          processing.msg,
          footer=NULL))
        withProgress(message=progress.caption,
                     value=0,
                     {
                       x <- sim_cv_obscov(te=te, bpue=bpue, d=d,
                                          shiny.progress=TRUE)
                     })
        removeModal()
        x
      })
    }
  })
  
  output$cv_obscov_plot <- renderPlot({
    if (!is.null(simlist())) {
      plot_cv_obscov(simlist = simlist(),
                     targetcv = as.numeric(input$targetcv),
                     silent = TRUE)
    }
  })

  output$cv_obscov_plot_label <- renderText({
    if (is.null(simlist())) {
      ''
    } else {
      oc.cv.out <- plot_cv_obscov(simlist = simlist(),
                                  targetcv = as.numeric(input$targetcv),
                                  silent = TRUE, as.shiny = TRUE)
      if (as.logical(input$targetcv)) {
        if (!is.na(oc.cv.out$pobs)) {
          rec <- paste0("Minimum observer coverage to achieve CV <= ",
                      input$targetcv, " is ",oc.cv.out$pobs,"%.")
        } else {
          rec <- paste0("Simulated observer coverage levels do not include ", 
                        "range corresponding to minimum observer coverage to ",
                        "achieve CV <= ", input$targetcv, ".\n")
        }
      } else rec <- ""
      HTML(paste0("<p></p>",rec," Results are interpolated from ",
                  "simulation-based projections and may vary slightly",
                  " with repetition.", 
                  " Please review the caveats in the About tab."))
    }
  })

  
  output$plotsavailable <- reactive(
    !is.null(simlist())
  )
  outputOptions(output, 'plotsavailable', suspendWhenHidden=FALSE)
}

