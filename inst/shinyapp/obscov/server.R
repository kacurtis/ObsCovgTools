#.libPaths(c("/usr/lib64/R/shiny_library",.libPaths()))
library(ObsCovgTools)

total.effort.msg <- paste("Total effort should be a positive integer",
                          "greater than 1.")
total.effort.title <- "Total effort value"

bpue.msg <- paste("BPUE should be a positive number.")
bpue.title <- "Bycatch per Unit Effort (BPUE) value"

dispersion.msg <- "Dispersion index should be a number greater than or equal to one."
dispersion.title <- "Dispersion index value"

processing.title <- HTML(paste0("<center>Computing - ",
                                "please wait...</center>"))
processing.msg <- HTML(paste0("<center>Reload your browser to cancel processing.", 
                              "<br>(To the left of the address bar, press Reload and then Stop Loading).",
                              "<center>"))

progress.caption <- "Simulation progress"

check.te.ppos <- function(input) {
  if (is.na(as.numeric(input)) || 
      ((as.numeric(input)) <= 1) ||
      (input != as.integer(input)) )
    total.effort.msg
  else NULL
}

check.bpue.ppos <- function(input) {
  if (is.na(as.numeric(input)) || 
      (as.numeric(input) <= 0) )
    bpue.msg
  else NULL
}

check.d.ppos <- function(input) {
  if (is.na(as.numeric(input)) ||
      (as.numeric(input) < 1) )
    dispersion.msg
  else NULL
}


server <- function(input, output, session) {
  
  plotlabels.ppos <- reactiveValues(ppos='')
  
  ## validation functions for user inputs to ppos tab
  
  output$ppos_obscov_plot <- renderPlot({
    
    ## validation of user inputs
    
    validate(check.te.ppos(input$te.ppos), 
             check.bpue.ppos(input$bpue.ppos), 
             check.d.ppos(input$d.ppos))
    
    plotlabels.ppos$ppos <- capture.output(
        plot_probposobs(te=as.numeric(input$te.ppos), 
                        bpue=as.numeric(input$bpue.ppos), 
                        d=as.numeric(input$d.ppos),
                        as.numeric(input$target.ppos)))
    
  })
  
  output$ppos_obscov_plot_label <- renderText({
    validate(check.te.ppos(input$te.ppos), 
             check.bpue.ppos(input$bpue.ppos), 
             check.d.ppos(input$d.ppos))
    
    oc.ppos.out <- plot_probposobs(te = as.numeric(input$te.ppos), 
                    bpue = as.numeric(input$bpue.ppos), 
                    d = as.numeric(input$d.ppos),
                    target.ppos = as.numeric(input$target.ppos),
                    silent = TRUE)
    
    if (as.logical(input$target.ppos)) {
      rec <- paste0("The probability that any bycatch occurs in the given total ",
                    "effort (horizontal black dotted line) is ",oc.ppos.out$ppos.te, 
                    "%. Minimum observer coverage to achieve at least ", 
                    input$target.ppos, "% probability of observing bycatch when",
                    " total bycatch is positive is ", oc.ppos.out$pobscov, "% (", 
                    oc.ppos.out$nobsets, " sets). ")
    } else rec <- paste0("The probability that any bycatch occurs in the given total ",
                         "effort (horizontal black dotted line) is ",oc.ppos.out$ppos.te, 
                         "%.")
    HTML(paste0(rec,
                "The conditional probability of observing any bycatch if it occurs ", 
                "(solid black line) is obtained by dividing the absolute probability",
                " of observing any bycatch (black dashed line) by the probability ",
                "that any bycatch occurs in the given total effort. Please review",
                " the caveat in the About tab."))
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
                                  silent = TRUE)
      if (as.logical(input$targetcv)) {
        if (!is.na(oc.cv.out$pobscov)) {
          rec <- paste0("Minimum observer coverage to achieve CV <= ",
                      input$targetcv, " is ",oc.cv.out$pobscov,"% (",
                      oc.cv.out$nobsets, " hauls).")
        } else {
          rec <- paste0("Simulated observer coverage levels do not include ", 
                        "range corresponding to minimum observer coverage to ",
                        "achieve CV <= ", input$targetcv, ".\n")
        }
      } else rec <- ""
      HTML(paste0("<p></p>",rec," Results are interpolated from ",
                  "simulation-based projections and may vary slightly",
                  " with repetition. Please review the caveats in the",
                  " About tab."))
    }
  })

  
  output$plotsavailable <- reactive(
    !is.null(simlist())
  )
  outputOptions(output, 'plotsavailable', suspendWhenHidden=FALSE)
}

