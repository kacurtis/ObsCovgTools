.libPaths(c("/usr/lib64/R/shiny_library",.libPaths()))
## require("appFrame")

ui <- fluidPage(
  tags$head(
    includeCSS("www/style.css")
  ),
  ## appFrameHeaderFixed(),
  titlePanel(title="Observer coverage simulator - 2.3",
             windowTitle="Observer coverage simulator"),
  tabsetPanel(id="alltabs",
              tabPanel("Objective: Probability of Observing Bycatch", value="ppos", 
                       fluidRow(
                         column(3, wellPanel(
                           textInput(inputId="te.ppos",
                                     label=HTML(paste("Total effort in fishery (e.g., hauls)")),
                                     value="500",
                                     placeholder="Integer > 1"),
                           textInput(inputId="bpue.ppos", label="Bycatch per Unit Effort (BPUE)",
                                     value="0.01",
                                     placeholder="Number > 0"),
                           textInput(inputId="d.ppos",
                                     label=HTML(paste("Dispersion index for BPUE",
                                                      "<div class='extext'>",
                                                      "(Dispersion index d is approximated by variance:mean",
                                                      " ratio of BPUE, with more skewed data having ",
                                                      "higher d. d=2 is a relatively conservative ",
                                                      "default.)",
                                                      "</div>")),
                                     value="2",
                                     placeholder="Number \u2265 1"),
                           sliderInput(inputId="target.ppos",
                                       label=HTML(paste("Target probability of observing bycatch",
                                       "<div class='extext'>",
                                       "(Given positive bycatch in total effort. Set to ",
                                       "zero to omit.)",
                                       "</div>")),
                                       min=0, max=100, value=80, step=1)
                           )),
                         column(8,
                                plotOutput("ppos_obscov_plot", width=700, height = 400),
                                htmlOutput("ppos_obscov_plot_label")
                         )
                       )
              ),
              tabPanel("Objective: Bycatch Estimation CV", value="cv",
  fluidRow(
    column(3, wellPanel(
      textInput(inputId="te",
                label=HTML(paste("Total effort in fishery (e.g., hauls)",
                                 "<div class='extext'>",
                                 "(Larger effort takes longer: ~30s for 5K,",
                                 "~10 min for 500K)",
                                 "</div>")),
                value="500",
                placeholder="Integer > 1"),
      textInput(inputId="bpue", label="Bycatch per Unit Effort (BPUE)",
                value="0.01",
                placeholder="Number > 0"),
      textInput(inputId="d",
                label=HTML(paste("Dispersion index for BPUE",
                                 "<div class='extext'>",
                                 "(Dispersion index d is approximated by variance:mean",
                                 " ratio of BPUE, with more skewed data having ",
                                 "higher d. d=2 is a relatively conservative ",
                                 "default.)",
                                 "</div>")),
                value="2",
                placeholder="Number \u2265 1"),
      actionButton(inputId="submit.cv", label="Submit"),
      conditionalPanel(condition="output.plotsavailable",
                       hr(),
                       sliderInput(inputId="targetcv",
                                   label=HTML(paste("Target CV for bycatch estimates",
                                                    "<div class='extext'>",
                                                    "(Set to zero to omit.)",
                                                    "</div>")),
                                   min=0, max=.99, value=0.30, step=0.01)))),
    column(8,
      plotOutput("cv_obscov_plot", width=700, height = 400),
      htmlOutput("cv_obscov_plot_label")
    )
  )
  ),
  tabPanel("About",
           h4("About the observer coverage simulator"),
           includeHTML("html/about.html")
           )
  )
  ## br(),
  ## br(),
  ## appFrameFooterFixed(displayAppsURL="../..")
)

