#' set size limit for input (9999mb)
options(
    shiny.maxRequestSize = 9999 * 1024 ^ 2, # size limit for input 9999mb
    scipen = 999 # disabling scientific notation
)

#' MAIN SERVER =================================================================
shinyServer(function(input, output, session) {
    # Automatically stop a Shiny app when closing the browser tab
    session$allowReconnect(TRUE)
    homePath = c(wd='~/') # for shinyFileChoose
    
    # input file ===============================================================
    getDomainFile <- reactive({
        shinyFileChoose(
            input, "domainFile", roots = homePath, session = session,
            filetypes = c('', 'domains')
        )
        fileSelected <- parseFilePaths(homePath, input$domainFile)
        return(replaceHomeCharacter(as.character(fileSelected$datapath)))
    })
    output$domainFile.ui <- renderUI({
        req(getDomainFile())
        if (length(getDomainFile()) > 0) {
            outString <- getDomainFile()
            if (nchar(outString) > 30)
                outString <- paste0(
                    substrLeft(outString, 15), "...", substrRight(outString, 15)
                )
            em(outString)
        }
    })
    
    # input folder =============================================================
    getDomainDir <- reactive({
        shinyDirChoose(
            input, "domainDir", roots = homePath, session = session
        )
        domainPath <- parseDirPath(homePath, input$domainDir)
        return(replaceHomeCharacter(as.character(domainPath)))
    })
    output$domainDir.ui <- renderUI({
        req(getDomainDir())
        if (length(getDomainDir()) > 0) {
            outString <- getDomainDir()
            if (nchar(outString) > 30)
                outString <- paste0(
                    substrLeft(outString, 15), "...", substrRight(outString, 15)
                )
            em(outString)
        }
    })
    
    # render seq IDs ==================================================================
    output$seedID.ui <- renderUI({
        if (input$inputType == "File") req(getDomainFile())
        if (input$inputType == "Folder") req(getDomainDir())
        selectInput(
            "seed",
            "Seed/Group ID",
            choices = c(getGroupIds(input$inputType, getDomainFile(), getDomainDir())[[1]])
        )
    })
    
    output$seqID.ui <- renderUI({
        req(input$seed)
        domainFile <- NULL
        if (input$inputType == "File") domainFile <- getDomainFile()
        if (input$inputType == "Folder") 
            domainFile <- paste0(getDomainDir(),"/",input$seed,".domains")
        if(!file.exists(domainFile)) stop(paste(domainFile, "not found!"))
        
        list(
            selectInput(
                "seq1",
                "Protein 1",
                getOrthoIDs(input$seed, domainFile),
            ),
            selectInput(
                "seq2",
                "Protein 2",
                getOrthoIDs(input$seed, domainFile),
            )
        )
    })
    
    # domain plot ==============================================================
    # * get domain info ========================================================
    getDomainInformation <- reactive({
        req(input$doPlot)
        withProgress(message = 'Reading domain input...', value = 0.5, {
            if (input$inputType == "File") {
                domainDf <- parseDomainInput(
                    input$seed,
                    getDomainFile(),
                    "file"
                )
            } else if (input$inputType == "Folder") {
                domainDf <- parseDomainInput(
                    input$seed,
                    getDomainDir(),
                    "folder"
                )
            }
            orthoIDtmp <- gsub("\\|",":",c(input$seq1, input$seq2))
            return(domainDf[domainDf$orthoID %in% orthoIDtmp,])
        })
    })
    
    output$domainPlot <- renderPlot({
        req(getDomainInformation())
        if (input$doPlot > 0) {
            if (is.null(getDomainInformation())) {
                msgPlot()
            } else {
                g <- createArchiPlot2(
                    c(input$seed, input$seq2), 
                    getDomainInformation(), 
                    input$labelArchiSize, input$titleArchiSize
                )
                if (any(g == "No domain info available!")) {
                    msgPlot()
                } else {
                    grid::grid.draw(g)
                }
            }
        }
    })
    output$domainPlot.ui <- renderUI({
        plotOutput(
            "domainPlot",
            height = input$archiHeight,
            width = input$archiWidth
        )
    })
    
    output$domainTable <- renderTable({
        req(getDomainInformation())
        req(input$seq2)
        features <- getDomainLink(c(input$seed, input$seq2), getDomainInformation())
        features
    }, sanitize.text.function = function(x) x)
    
    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            g <- createArchiPlot2(
                c(input$seed, input$seq2), 
                getDomainInformation(), 
                input$labelArchiSize, input$titleArchiSize
            )
            grid.draw(g)
            ggsave(
                file, plot = g,
                width = input$archiWidth * 0.056458333,
                height = input$archiHeight * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )
})