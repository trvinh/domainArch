shinyUI(
    fluidPage(

        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        useShinyjs(),

        # Application title
        titlePanel("", windowTitle = "domainArch"),
        navbarPage(
            em(strong("domainArch v0.0.7")),
            id = "tabs",
            collapsible = TRUE,
            inverse = TRUE,
            fluid = TRUE,
            position = "fixed-top",
            tabPanel(
                "Playground",
                sidebarPanel(
                    width = 3,
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
                    selectizeInput(
                        "seq1",
                        "Protein 1",
                        choices = NULL
                    ),
                    uiOutput("seedProtID.ui"),
                    selectizeInput(
                        "seq2",
                        "Protein 2",
                        choices = NULL
                    ),
                    br(),
                    selectInput(
                        "font","Font", choices = extrafont::fonts(),
                        selected = "Arial"
                    ),
                    bsButton("doPlot", "Plot", style = "info")
                ),
                mainPanel(
                    # createArchitecturePlotUI("archiPlot")
                    fluidRow(
                        column(
                            3,
                            radioButtons(
                                "resolveOverlap",
                                "Merge non-overlapped features",
                                choices = c("Yes","No"), selected = "Yes",
                                inline = TRUE
                            ),
                            checkboxGroupInput(
                                "namePosition",
                                "Display feature names",
                                choices = c(
                                    "On the plot" = "plot",
                                    "As a legend" = "legend",
                                    "On the y-axis" = "axis"
                                ),
                                selected = c("plot","axis")
                            )
                        ),
                        column(
                            3,
                            selectInput(
                                "feature",
                                "Exclude features",
                                choices = c(
                                    "flps","seg","coils","signalp","tmhmm",
                                    "smart","pfam",
                                    "without E-value" = "noEvalue",
                                    "without Bit-score" = "noBitscore"
                                ),
                                multiple = TRUE
                            ),
                            checkboxInput(
                                "featureOpt", "Other feature options", value = FALSE
                            ),
                            checkboxInput("plotConfig", "Plot configuration", value = FALSE)
                        ),
                        column(
                           3,
                           strong("Show information"),
                           checkboxGroupInput(
                               "showWeight",
                               "",
                               choices = "Weight"
                           ),
                           checkboxGroupInput(
                               "showScore",
                               "",
                               choices = c(
                                   "E-value", "Bit-score"
                               )
                           ),
                           uiOutput("filterEvalue.ui"),
                           uiOutput("filterBitscore.ui")
                        ),
                        column(
                            3,
                            checkboxGroupInput(
                                "showInstance",
                                "Show only instances with",
                                choices = c(
                                    "Best E-value" = "evalue",
                                    "Best Bit-score" = "bitscore",
                                    "Paths" = "path"
                                )
                            )
                        )
                    ),
                    br(),
                    fluidRow(
                        conditionalPanel(
                            condition = "input.featureOpt == 1",
                            column(
                                3,
                                radioButtons(
                                    "nameType","Type of feature names", inline = TRUE,
                                    choices = c("Labels","Texts"), selected = "Labels"
                                )
                                
                            ),
                            column(
                                4,
                                conditionalPanel(
                                    condition = "input.nameType == 'Labels'",
                                    radioButtons(
                                        "labelPos","Label position", inline = TRUE,
                                        choices = c("Above","Inside","Below"),
                                        selected = "Above"
                                    )
                                ),
                                conditionalPanel(
                                    condition = "input.nameType == 'Texts'",
                                    colourpicker::colourInput(
                                        "nameColor",
                                        "Feature name color",
                                        value = "#000000"
                                    )
                                )
                            ),
                            column(
                                5,
                                selectInput(
                                    "excludeNames",
                                    "Exclude feature names of",
                                    choices = c(
                                        "flps","seg","coils","signalp","tmhmm",
                                        "smart","pfam"
                                    ),
                                    selected = c("tmhmm","signalp","seg","coils"),
                                    multiple = TRUE
                                )
                            ),
                            column(
                                3,
                                radioButtons(
                                    "featureClassSort",
                                    "Sort feature classes by shared features",
                                    choices = c("Yes","No"), selected = "Yes",
                                    inline = TRUE
                                )
                            ),
                            column(
                                5,
                                conditionalPanel(
                                    condition = "input.featureClassSort == 'No'",
                                    selectInput(
                                        "featureClassOrder",
                                        "Feature class order",
                                        choices = c(
                                            "pfam", "smart", "tmhmm", "coils", "signalp",
                                            "seg", "flps"
                                        ),
                                        selected = c(
                                            "pfam", "smart", "tmhmm", "coils", "signalp",
                                            "seg", "flps"
                                        ),
                                        multiple = TRUE
                                    )
                                )
                            )
                        )
                    ),
                    br(),
                    fluidRow(
                        conditionalPanel(
                            condition = "input.plotConfig == 1",
                            column(
                                3,
                                createPlotSize("archiHeight", "Plot height(px)",400, 200),
                                createPlotSize("archiWidth", "Plot width (px)", 800, 200)
                            ),
                            column(
                                3,
                                createTextSize(
                                    "titleArchiSize", "Title/Seq ID size (px)", 14, 200
                                ),
                                createTextSize(
                                    "labelArchiSize", "Axis label size(px)", 12, 200
                                )
                            ),
                            column(
                                6,
                                column(
                                    6,
                                    createTextSize(
                                        "segmentSize", "Feature segment size (mm)", 5, 200
                                    )
                                ),
                                column(
                                    6,
                                    createTextSize(
                                        "nameSize", "Feature ID size (mm)", 3, 200
                                    )
                                ),
                                column(
                                    12,
                                    sliderInput(
                                        "firstDist", "Distance between plot title and the 1st feature",
                                        min = 0, max = 5, value = 0.5, step = 0.1, width = 400
                                    )
                                )
                            ),
                            column(
                                6,
                                radioButtons(
                                    "colorType","Color feature instances", inline = TRUE,
                                    choices = c("Shared","Unique","All","Feature class"), selected = "All"
                                ),
                                checkboxInput(
                                    "ignoreInstanceNo", "Ignore number of instances", value = FALSE
                                ),
                            ),
                            column(
                                6,
                                selectInput(
                                    "colorPallete",
                                    "Color pallete",
                                    choices = c("Paired", "Set1", "Set2", "Set3", "Accent", "Dark2"),
                                    selected = "Paired"
                                ),
                                selectInput(
                                    "font",
                                    "Font",
                                    choices = fonts(),
                                    selected = "Arial"
                                )
                            )
                        )
                    ),
                    hr(),
                    uiOutput("domainPlot.ui"),
                    verbatimTextOutput("hover_info"),
                    br(),
                    downloadButton("archiDownload", "Download plot", class = "butDL"),
                    hr(),
                    tableOutput("linkTable"),
                    checkboxInput(
                        "showDomainTable", "Show detailed feature table", value = FALSE
                    ),
                    conditionalPanel(
                        condition = "input.showDomainTable == 1",
                        DT::dataTableOutput("domainTable")
                    )
                )
            ),
            navbarMenu(
                "Functions",
                tabPanel(
                    "Split domain file",
                    shinyFilesButton(
                        "domainFileIn", "Input domain file" ,
                        title = "Please provide domain file:",
                        multiple = FALSE,
                        buttonType = "default", class = NULL
                    ),
                    uiOutput("domainFileIn.ui"),
                    br(),
                    shinyDirButton(
                        "domainDirOut",
                        "Select output directory" ,
                        title = paste(
                            "Please select output directory"
                        ),
                        buttonType = "default", class = NULL
                    ),
                    br(),
                    uiOutput("domainDirOut.ui"),
                    br(),
                    bsButton(
                        "doSplitDomain",
                        "Split domain file",
                        style = "warning",
                        icon("file-export")
                    ),
                    hr(),
                    verbatimTextOutput("splitDomainFileStatus")
                )
            )
        )
    )
)
