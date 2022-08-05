#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c(
    "data.table", "DT", "ggplot2", "gridExtra", "pbapply", "RColorBrewer",
    "shiny", "shinyBS", "shinyFiles", "shinyjs", "shinyalert", "shinythemes", "plyr", 
    "PhyloProfile", "stringr"
)

# Load packages
lapply(packages, library, character.only = TRUE)