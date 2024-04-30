#' Run domainArch app
#' @export
#' @return A shiny application
#' @import BiocStyle
#' @import shinyBS
#' @import DT
#' @import PhyloProfile
#' @import RColorBrewer
#' @importFrom data.table tstrsplit fread  
#' @import dplyr
#' @import ggplot2
#' @import grid
#' @importFrom gridExtra arrangeGrob
#' @import jsonlite
#' @import pbapply
#' @import shinyFiles
#' @import shinythemes
#' @import stringr
#' @rawNamespace import(shinyjs, except = colourInput)

runDomainArch <- function(){
    appDir <- system.file("domainArch", package = "domainArch")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `domainArch`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
