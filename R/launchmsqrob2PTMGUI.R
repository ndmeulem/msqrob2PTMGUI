#' Launches the msqrob2PTMGUI Shiny App
#'
#' @param maxSize maximum memory size that input files are allowed to have in Mb
#'
#' @export launchmsqrob2PTMGUI
#'
#' @return shiny application object
#'
#' @example
#' \dontrun{launchmsqrob2PTMGUI()}
#'
#' @import shiny shinymeta shinyjs BiocManager
#'


# wrapper for shiny::shinyApp()
launchmsqrob2PTMGUI <- function(maxSize=500) {
  shinyjs::useShinyjs()
  onStart = shinyhelper::observe_helpers(help_dir = system.file("helpfiles", package="msqrob2PTMGUI"))
  options(shiny.maxRequestSize=maxSize*1024^2)
  shinyApp(ui = ui, server = server)
}
