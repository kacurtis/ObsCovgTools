# Hidden functions to execute progress bar in Shiny

progress_init <- function(shiny.progress = FALSE) {
  if (shiny.progress) { return(NULL)
  } else { return(utils::txtProgressBar(style=3)) }
} 

progbar <-  function(i, total, pb, shiny.progress = FALSE) {
  if (shiny.progress) {
    shiny::incProgress(500 / total)
    return(NULL)
  } else {
    utils::setTxtProgressBar(pb, i/total)
    return(pb)
  }
}
