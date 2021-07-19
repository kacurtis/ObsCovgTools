.libPaths(c("/usr/lib64/R/shiny_library",.libPaths()))
## require("appFrame")

ui <- fluidPage(
  tags$head(
    includeCSS("www/style.css")
  ),
  ## appFrameHeaderFixed(),
  titlePanel(title=div(img(src="logo_120.png"), "Observer coverage simulator - 3.1.1"),
             windowTitle="Observer coverage simulator"),
  tabsetPanel(id="alltabs",
              tabPanel("Objective: Probability of Positive Bycatch", value="ppos", 
                       fluidRow(
                         column(3, wellPanel(
                           textInput(inputId="te.ppos",
                                     label=HTML(paste("Total effort in fishery (e.g., trips)")),
                                     value="500",
                                     placeholder="Integer > 1"),
                           textInput(inputId="bpue.ppos", label="Bycatch per Unit Effort (BPUE)",
                                     value="0.01",
                                     placeholder="Number > 0"),
                           textInput(inputId="d.ppos",
                                     label=HTML(paste("Dispersion index for BPUE",
                                                      "<div class='extext'>",
                                                      "(Dispersion index d is variance:mean",
                                                      " ratio of BPUE, with more skewed data having ",
                                                      "higher d. d=2 is a relatively conservative ",
                                                      "default for rare event bycatch of marine mammals",
                                                      " and turtles. Seabirds and fishes tend to be more",
                                                      " skewed.)",
                                                      "</div>")),
                                     value="2",
                                     placeholder="Number \u2265 1"),
                           sliderInput(inputId="target.ppos",
                                       label=HTML(paste("Target probability of observing bycatch",
                                       "<div class='extext'>",
                                       "(Given positive bycatch in total effort. Set to ",
                                       "zero to omit.)",
                                       "</div>")),
                                       min=0, max=100, value=95, step=1)
                           )),
                         column(8,
                                plotOutput("ppos_obscov_plot", width=700, height = 400),
                                htmlOutput("ppos_obscov_plot_label")
                         )
                       )
              ),
              tabPanel("Objective: Upper Confidence Limit Given No Bycatch Observed", value="ucl", 
                       fluidRow(
                         column(3, wellPanel(
                           textInput(inputId="te.ucl",
                                     label=HTML(paste("Total effort in fishery (e.g., trips)")),
                                     value="500",
                                     placeholder="Integer > 1"),
                           textInput(inputId="d.ucl",
                                     label=HTML(paste("Dispersion index for BPUE",
                                                      "<div class='extext'>",
                                                      "(Dispersion index d is variance:mean",
                                                      " ratio of BPUE, with more skewed data having ",
                                                      "higher d. d=2 is a relatively conservative ",
                                                      "default for rare event bycatch of marine mammals",
                                                      " and turtles. Seabirds and fishes tend to be more",
                                                      " skewed.)",
                                                      "</div>")),
                                     value="2",
                                     placeholder="Number \u2265 1"),
                           sliderInput(inputId="cl.ucl",
                                       label=HTML(paste("Confidence level for upper confidence limit (UCL)")),
                                       min=50, max=100, value=95, step=1),
                           textInput(inputId="target.ucl",
                                     label=HTML(paste("Target maximum UCL given zero bycatch observed",
                                                      "<div class='extext'>",
                                                      "(Set to zero to omit.)",
                                                      "</div>")),
                                     value="0",
                                     placeholder="Number \u2265 0"),
                           textInput(inputId="fixedoc.ucl",
                                     label=HTML(paste("Percent observer coverage for which to return UCL")),
                                     value="0",
                                     placeholder="Number between 0 and 100"),
                           textInput(inputId="ymax.ucl",
                                     label=HTML(paste("Upper limit for y-axis")),
                                     value="100",
                                     placeholder="Number > 0")
                         )),
                         column(8,
                                plotOutput("ucl_obscov_plot", width=700, height = 400),
                                htmlOutput("ucl_obscov_plot_label")
                         )
                       )
              ),
              tabPanel("Objective: Bycatch Estimation CV", value="cv",
                       fluidRow(
                         column(3, wellPanel(
                           textInput(inputId="te",
                                     label=HTML(paste("Total effort in fishery (e.g., trips)",
                                                      "<div class='extext'>",
                                                      "(Larger effort takes longer: ~20s for 10K,",
                                                      "~75s for 100K)",
                                                      "</div>")),
                                     value="500",
                                     placeholder="Integer > 1"),
                           textInput(inputId="bpue", label="Bycatch per Unit Effort (BPUE)",
                                     value="0.01",
                                     placeholder="Number > 0"),
                           textInput(inputId="d",
                                     label=HTML(paste("Dispersion index for BPUE",
                                                      "<div class='extext'>",
                                                      "(Dispersion index d is variance:mean",
                                                      " ratio of BPUE, with more skewed data having ",
                                                      "higher d. d=2 is a relatively conservative ",
                                                      "default for rare event bycatch of marine mammals",
                                                      " and turtles. Seabirds and fishes tend to be more",
                                                      " skewed.)",
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

