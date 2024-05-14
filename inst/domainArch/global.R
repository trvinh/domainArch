#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c(
    "BiocStyle", "data.table", "DT", "dplyr", "extrafont","ggplot2","gridExtra", 
    "grid", "jsonlite", "pbapply", "PhyloProfile", "RColorBrewer", "shiny", 
    "shinyBS", "shinyFiles", "shinyjs", "shinythemes",  "stringr"
)

# Load packages
lapply(packages, library, character.only = TRUE)