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
    
    # input anno folder ========================================================
    getAnnoDir <- reactive({
        shinyDirChoose(
            input, "annoDir", roots = homePath, session = session
        )
        annoPath <- parseDirPath(homePath, input$annoDir)
        return(replaceHomeCharacter(as.character(annoPath)))
    })
    output$annoDir.ui <- renderUI({
        req(getAnnoDir())
        if (length(getAnnoDir()) > 0) {
            outString <- getAnnoDir()
            if (nchar(outString) > 30)
                outString <- paste0(
                    substrLeft(outString, 15), "...", substrRight(outString, 15)
                )
            em(outString)
        }
    })
    
    # render seed/species IDs ==================================================
    output$seedID.ui <- renderUI({
        if (input$inputType == "File") req(getDomainFile())
        if (input$inputType == "Folder") req(getDomainDir())
        type <- 1
        if (input$inputType == "Anno") {
            req(getAnnoDir())
            type <- 2
        }
        if (type == 1) {
            selectInput(
                "seed1",
                "Seed/Group ID",
                choices = c(getGroupIds(input$inputType, getDomainFile(), getDomainDir()))
            )
        } else if (type == 2) {
            list(
                selectInput(
                    "seed1",
                    "Species 1",
                    choices = c(getSpecList(getAnnoDir()))
                ),
                selectInput(
                    "seed2",
                    "Species 2",
                    choices = c("none", getSpecList(getAnnoDir()))
                )
            )
        }
    })
    
    # render protein IDs =======================================================
    getJonsList <- reactive({
        req(input$seed1)
        req(input$seed2)
        jsonList <- list()
        if (input$inputType == "Anno") {
            domainFile <- paste0(getAnnoDir(),"/",input$seed1,".json")
            jsonList [[1]] <- json2list(domainFile)
            if (input$seed2 != "none") {
                domainFile2 <- paste0(getAnnoDir(),"/",input$seed2,".json")
                jsonList[[2]] <- json2list(domainFile2)
            } else {
                jsonList[[2]] <- jsonList[[1]]
            }
        }
        return(jsonList)
    })
    
    output$seqID.ui <- renderUI({
        req(input$seed1)
        domainFile <- domainFile2 <- NULL
        
        if (input$inputType == "File") domainFile <- getDomainFile()
        if (input$inputType == "Folder") 
            domainFile <- paste0(getDomainDir(),"/",input$seed1,".domains")
        if (input$inputType == "Anno") {
            domainFile <- paste0(getAnnoDir(),"/",input$seed1,".json")
            if (input$seed2 != "none")
                domainFile2 <- paste0(getAnnoDir(),"/",input$seed2,".json")
        }
            
        if(!file.exists(domainFile)) stop(paste(domainFile, "not found!"))
        
        if (grepl(".domains", domainFile)) {
            list(
                selectInput(
                    "seq1",
                    "Protein 1",
                    getOrthoIDs(input$seed1, domainFile),
                ),
                selectInput(
                    "seq2",
                    "Protein 2",
                    c("none", getOrthoIDs(input$seed1, domainFile)),
                    selected = "none"
                )
            )
        } else {
            withProgress(message = 'Reading JSON input...', value = 0.5, {
                jsonList <- getJonsList()
                list(
                    selectInput(
                        "seq1",
                        "Protein 1",
                        names(jsonList[[1]]$feature)
                    ),
                    selectInput(
                        "seq2",
                        "Protein 2",
                        c("none", names(jsonList[[2]]$feature))
                    )
                )
            })
        }
    })
    
    # domain plot ==============================================================
    # * get domain info ========================================================
    getDomainInformation <- reactive({
        req(input$doPlot)
        withProgress(message = 'Reading domain input...', value = 0.5, {
            outDf <- NULL
            if (input$inputType != "Anno") {
                if (input$inputType == "File") {
                    domainDf <- parseDomainInput(
                        input$seed1,
                        getDomainFile(),
                        "file"
                    )
                } else if (input$inputType == "Folder") {
                    domainDf <- parseDomainInput(
                        input$seed1,
                        getDomainDir(),
                        "folder"
                    )
                } 
                # filter domain df by orthoIDs
                seq2 <- input$seq2
                if (input$seq2 == "none") seq2 <- input$seq1
                orthoIDtmp <- gsub("\\|",":",c(input$seq1, seq2))
                outDf <- domainDf[domainDf$orthoID %in% orthoIDtmp,]
            } else if (input$inputType == "Anno") {
                jsonList <- getJonsList()
                domainDf1 <- parseDomainFromJson(
                    input$seed1,
                    input$seq1,
                    jsonList[[1]]
                )
                if (input$seed2 == "none" || input$seed2 == input$seed1) {
                    seed2 <- input$seed1
                    if (input$seq2 == "none" || input$seq2 == input$seq1) {
                        outDf <- domainDf1
                    } else {
                        domainDf2 <- parseDomainFromJson(
                            seed2,
                            input$seq2,
                            jsonList[[2]]
                        )
                        outDf <- rbind(domainDf1, domainDf2)
                    }
                } else {
                    if (input$seq2 == "none") {
                        outDf <- domainDf1
                    } else {
                        domainDf2 <- parseDomainFromJson(
                            input$seed2,
                            input$seq2,
                            jsonList[[2]]
                        )
                        outDf <- rbind(domainDf1, domainDf2)
                    }
                }
            }
            # filter domain df by features
            outDf[c("feature_type","feature_id")] <- str_split_fixed(outDf$feature, '_', 2)
            outDf <- outDf[!(outDf$feature_type %in% input$feature),]
            return(outDf)
        })
    })
    
    output$domainPlot <- renderPlot({
        req(getDomainInformation())
        if (input$doPlot > 0) {
            if (is.null(getDomainInformation())) {
                msgPlot()
            } else {
                seq2 <- input$seq2
                if (input$seq2 == "none") seq2 <- input$seq1
                g <- createArchiPlot2(
                    c(input$seed1, seq2), 
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
        if (is.null(getDomainInformation())) {
            msg <- paste0(
                "<p><em>No domain found for this protein.</em></p>"
            )
            HTML(msg)
        } else {
            plotOutput(
                "domainPlot",
                height = input$archiHeight,
                width = input$archiWidth
            )
        }
    })
    
    output$domainTable <- renderTable({
        req(getDomainInformation())
        req(input$seq1)
        req(input$seq2)
        seq2 <- input$seq2
        if (input$seq2 == "none") seq2 <- input$seq1
        features <- getDomainLink(c(input$seed1, input$seq1), getDomainInformation())
        if(!(is.null(input$seed2))) {
            seed2 <- input$seed2
            if(input$seed2 == "none") seed2 <- input$seed1
            features2 <- getDomainLink(c(seed2, seq2), getDomainInformation())
            features <- rbind(features, features2)
        }
        features <- features[!duplicated(features),]
    }, sanitize.text.function = function(x) x)
    
    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            seq2 <- input$seq2
            if (input$seq2 == "none") seq2 <- input$seq1
            g <- createArchiPlot2(
                c(input$seed1, seq2), 
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