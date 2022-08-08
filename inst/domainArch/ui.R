shinyUI(
    fluidPage(

        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        useShinyjs(),

        # Application title
        titlePanel("", windowTitle = "domainArch"),
        navbarPage(
            em(strong("domainArch v0.0.1")),
            id = "tabs",
            collapsible = TRUE,
            inverse = TRUE,
            fluid = TRUE,
            position = "fixed-top",
            tabPanel(
                "Playground",
                sidebarPanel(
                    width = 3,
                    # radioButtons(
                    selectInput(
                        "inputType", "Input type",
                        c(
                            "Domain file" = "File",
                            "Domain folder" = "Folder",
                            "Annotation folder" = "Anno"
                        ),
                        width = 150
                    ),
                    conditionalPanel(
                        condition = "input.inputType == 'File'",
                        shinyFilesButton(
                            "domainFile", "Input file" ,
                            title = "Please provide domain file:",
                            multiple = FALSE,
                            buttonType = "default", class = NULL
                        ),
                        uiOutput("domainFile.ui")
                    ),
                    conditionalPanel(
                        condition = "input.inputType == 'Folder'",
                        shinyDirButton(
                            "domainDir", "Input directory" ,
                            title = "Please select a folder",
                            buttonType = "default", class = NULL
                        ),
                        uiOutput("domainDir.ui"),
                        br()
                    ),
                    conditionalPanel(
                        condition = "input.inputType == 'Anno'",
                        shinyDirButton(
                            "annoDir", "Input directory" ,
                            title = "Please select a folder",
                            buttonType = "default", class = NULL
                        ),
                        uiOutput("annoDir.ui"),
                        br()
                    ),
                    br(),
                    uiOutput("seedID.ui"),
                    # conditionalPanel(
                    #     condition = "input.inputType == 'Anno'",
                    #     bsButton("loadJson", "Load annotation file(s)")
                    # ),
                    uiOutput("seqID.ui"),
                    bsButton("doPlot", "Plot", style = "info")
                ),
                mainPanel(
                    fluidRow(
                        column(
                            2,
                            createPlotSize("archiHeight", "Plot height(px)",400)
                        ),
                        column(
                            2,
                            createPlotSize("archiWidth", "Plot width(px)", 800)
                        ),
                        column(
                            2,
                            createTextSize(
                                "titleArchiSize", "Title size(px)", 14, 100
                            )
                        ),
                        column(
                            2,
                            createTextSize(
                                "labelArchiSize", "DomainID size(px)", 12, 150
                            )
                        ),
                        column(
                            4,
                            selectInput(
                                "feature",
                                "Exclude features",
                                choices = c("flps","seg","coils","signalp","tmhmm","smart","pfam"),
                                multiple = TRUE
                            )
                        )
                    ),
                    uiOutput("domainPlot.ui"),
                    downloadButton("archiDownload", "Download plot", class = "butDL"),
                    hr(),
                    tableOutput("domainTable")
                )
            )
        )
    )
)
