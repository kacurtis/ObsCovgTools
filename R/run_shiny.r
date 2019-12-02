#' Run Shiny interface with ObsCovgTools
#' 
#' \code{run_shiny} runs a shiny application for the main functions in ObsCovgTools.
#'   
#' @return None 
#' 
#' @export
run_shiny <- function() {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package \"shiny\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
    
  appDir <- system.file("shinyapp", "obscov", package = "ObsCovgTools")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ObsCovgTools`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
  
}