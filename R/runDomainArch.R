#' Run domainArch app
#' @export
#' @return A shiny application
#' @import BiocStyle
#' @import shinyBS
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
