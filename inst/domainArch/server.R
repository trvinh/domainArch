#' set size limit for input (9999mb)
options(
    shiny.maxRequestSize = 9999 * 1024 ^ 2, # size limit for input 9999mb
    scipen = 999, # disabling scientific notation
    htmlwidgets.TOJSON_ARGS = list(na = 'string') # to show NA values in table
)
if (is.null(fonts())) extrafont::font_import()

#' MAIN SERVER =================================================================
shinyServer(function(input, output, session) {
    # Automatically stop a Shiny app when closing the browser tab
    session$allowReconnect(TRUE)
    homePath = c(wd='~/') # for shinyFileChoose
    
    # load ncbi taxonomy db from PhyloProfile
    defaultTaxDB <- system.file(
        "PhyloProfile", "data", package = "PhyloProfile", mustWork = TRUE
    )
    currentNCBIinfo <- NULL
    if (file.exists(paste0(defaultTaxDB, "/preProcessedTaxonomy.txt"))) {
        currentNCBIinfo <- as.data.frame(data.table::fread(paste0(defaultTaxDB, "/preProcessedTaxonomy.txt")))
    }

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
    
    # disable input after plotting =============================================
    observeEvent(input$doPlot,{
        shinyjs::disable("inputType")
        if (input$inputType == "File") shinyjs::disable("domainFile")
        if (input$inputType == "Folder") shinyjs::disable("domainDir")
        if (input$inputType == "Anno") shinyjs::disable("annoDir")
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
            choices = c(getGroupIds(input$inputType, getDomainFile(), getDomainDir()))
            selectInput(
                "seed1",
                "Seed/Group ID",
                choices = choices[nzchar(choices)]
            )
        } else if (type == 2) {
            choices = c(getSpecList(getAnnoDir()))
            list(
                selectInput(
                    "seed1",
                    "Species 1",
                    choices = choices[nzchar(choices)]
                ),
                selectInput(
                    "seed2",
                    "Species 2",
                    choices = c("none", choices[nzchar(choices)])
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
    
    observe({
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
            updateSelectizeInput(
                session, "seq1",
                "Protein 1",
                c("seed", getOrthoIDs(input$seed1, domainFile, currentNCBIinfo)),
                selected = "seed",
                server = TRUE
            )
            updateSelectizeInput(
                session, "seq2",
                "Protein 2",
                c("none", getOrthoIDs(input$seed1, domainFile, currentNCBIinfo)),
                selected = "none",
                server = TRUE
            )
                
        } else {
            withProgress(message = 'Reading JSON input...', value = 0.5, {
                jsonList <- getJonsList()
                updateSelectizeInput(
                    session, "seq1",
                    "Protein 1",
                    names(jsonList[[1]]$feature),
                    selected = names(jsonList[[1]]$feature)[1],
                    server = TRUE
                )
                updateSelectizeInput(
                    session, "seq2",
                    "Protein 2",
                    c("none", names(jsonList[[2]]$feature)),
                    selected = "none",
                    server = TRUE
                )
            })
        }
    })
    
    # render seed ID ===========================================================
    output$seedProtID.ui <- renderUI({
        req(getDomainInformation())
        if (input$seq1 == "seed") {
            df <- getDomainInformation()
            if (nrow(df) > 0) {
                tmp <- unlist(lapply(
                    df$orthoID, function(x) ifelse(grepl(x, df$seedID) == FALSE, x, NA)
                ))
                seed <- unique(tmp[!is.na(tmp)])
                seedId <- sapply(str_split(seed, "\\:"), "[", 3)
                seedTaxId <- sapply(str_split(seed, "@"), "[", 2)
                seedTaxName <- PhyloProfile::id2name(seedTaxId, currentNCBIinfo)
                list(
                    em(paste(seedTaxName$fullName, seedId, collapse = " - ")),
                    br(), br()
                )
            }
        }
    })
    
    # remove option to show best path if no seed protein =======================
    observe({
        if (input$seq1 != "seed") {
            updateCheckboxGroupInput(
                session, "showInstance",
                "Show only instances with",
                choices = c(
                    "Best E-value" = "evalue", 
                    "Best Bit-score" = "bitscore"
                )
            )
        } else {
            updateCheckboxGroupInput(
                session, "showInstance",
                "Show only instances with",
                choices = c(
                    "Best E-value" = "evalue", 
                    "Best Bit-score" = "bitscore",
                    "Paths" = "path"
                )
            )
        }
    })
    
    # update excludeNames if no feature type on the y-axis =====================
    observe({
        req(input$showName)
        if (!("axis" %in% input$showName) & !("legend" %in% input$showName)) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm",
                    "smart","pfam"
                )
            )
        } else if ("axis" %in% input$showName | "legend" %in% input$showName) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm",
                    "smart","pfam"
                ),
                selected = c("tmhmm","signalp","flps","seg","coils")
            )
        }
    })
    
    # domain plot ==============================================================
    # * get domain info ========================================================
    getDomainInformation <- reactive({
        req(input$doPlot)
        # withProgress(message = 'Reading domain input...', value = 0.5, {
            outDf <- NULL
            if (input$inputType != "Anno") {
                if (input$seq1 == "seed" & input$seq2 == "none")
                    stop("Please specify protein 1 and/or protein 2!")
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
                if (input$seq1 == "seed") {
                    seq2Tmp <- gsub("\\|",":",input$seq2)
                    orthoDf <- domainDf[domainDf$orthoID == seq2Tmp,]
                    seedDf <- domainDf[domainDf$seedID == unique(orthoDf$seedID),]
                    outDf <- rbind(seedDf, orthoDf)
                } else {
                    orthoIDtmp <- gsub("\\|",":",c(input$seq1, input$seq2))
                    outDf <- domainDf[domainDf$orthoID %in% orthoIDtmp,]
                }
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
            return(outDf)
        # })
    })
    
    # * render e-value / bitscore filter =======================================
    output$filterEvalue.ui <- renderUI({
        req(getDomainInformation())
        df <- getDomainInformation()
        maxEvalue = 1
        if ("evalue" %in% colnames(df))
            maxEvalue <- format(max(df$evalue[!is.na(df$evalue)]), scientific = TRUE, digits = 2)
        if ("E-value" %in% input$showScore) {
            numericInput(
                "minEvalue", "Filter E-value:",
                min = 0,
                max = maxEvalue,
                value = format(0.00001, scientific = TRUE, digits = 2)
            )
        }
    })
    
    output$filterBitscore.ui <- renderUI({
        req(getDomainInformation())
        df <- getDomainInformation()
        if ("Bit-score" %in% input$showScore) {
            numericInput(
                "minBitscore", "Filter Bit-score:",
                min = min(df$bitscore[!is.na(df$bitscore)]),
                max = 9999,
                value = min(df$bitscore[!is.na(df$bitscore)])
            )
        }
    })
    
    # * filter data ============================================================
    filterDomainData <- reactive({
        req(getDomainInformation())
        outDf <- getDomainInformation()
        
        # filter domain df by features
        if (is.null(outDf)) return(NULL)
        if (nrow(outDf) == 0) return(NULL)
        outDf[c("feature_type","feature_id")] <- str_split_fixed(outDf$feature, '_', 2)
        outDf <- outDf[!(outDf$feature_type %in% input$feature),]
        # filter filters without e-value and/or bitscore
        if ("evalue" %in% colnames(outDf)) {
            if ("noEvalue" %in% input$feature)
                outDf <- outDf[!is.na(outDf$evalue),]
            if ("noBitscore" %in% input$feature)
                outDf <- outDf[!is.na(outDf$bitscore),]
        }
        # modify feature IDs
        outDf$feature_id_mod <- outDf$feature_id
        outDf$feature_id_mod <- gsub("SINGLE", "LCR", outDf$feature_id_mod)
        outDf$feature_id_mod[outDf$feature_type == "coils"] <- "Coils"
        outDf$feature_id_mod[outDf$feature_type == "seg"] <- "LCR"
        outDf$feature_id_mod[outDf$feature_type == "tmhmm"] <- "TM"
        # exclude features IDs
        if (!is.null(input$excludeNames)) {
            outDf$feature_id_mod[outDf$feature_type %in% input$excludeNames] <- NA
        }
        
        # enable/disable option for showing evalue/bitscore
        if ("evalue" %in% colnames(outDf)) {
            shinyjs::enable("showScore")
        } else {
            shinyjs::disable("showScore")
        }
        
        # filter data by e-value, bit-score and feature path
        if ("evalue" %in% colnames(outDf)) {
            # filter by e-value and/or bit-score
            if ("E-value" %in% input$showScore) {
                req(input$minEvalue)
                minEvalue <- format(input$minEvalue, scientific = FALSE)
                naOutDf <- outDf[is.na(outDf$evalue),]
                outDf <- outDf[!is.na(outDf$evalue) & outDf$evalue <= input$minEvalue,]
                outDf <- rbind(outDf,naOutDf)
            }   
            if ("Bit-score" %in% input$showScore) {
                req(input$minBitscore)
                naOutDf <- outDf[is.na(outDf$bitscore),]
                outDf <- outDf[!is.na(outDf$bitscore) & outDf$bitscore >= input$minBitscore,]
                outDf <- rbind(outDf,naOutDf)
            }
            # get only best instances
            if ("evalue" %in% input$showInstance) {
                naOutDf <- outDf[is.na(outDf$evalue),]
                outDf <- outDf %>% group_by(feature, orthoID) %>% filter(evalue == min(evalue))
                outDf <- rbind(outDf,naOutDf)
            }
            if ("bitscore" %in% input$showInstance) {
                naOutDf <- outDf[is.na(outDf$bitscore),]
                outDf <- outDf %>% group_by(feature, orthoID) %>% filter(bitscore == max(bitscore))
                outDf <- rbind(outDf,naOutDf)
            }
            if ("path" %in% input$showInstance) {
                outDf <- outDf %>% group_by(feature) %>% filter(path == "Y")
            }
            # Format e-values
            outDf$evalue[!is.na(outDf$evalue)] <- 
                format(outDf$evalue[!is.na(outDf$evalue)], scientific = TRUE, digits = 2)
        }
        # return
        return(outDf[!is.na(outDf$seedID),])
    })
    
    # * create domain plot =====================================================
    # observeEvent(input$doPlot, {
    #     seq2 <- input$seq2
    #     if (input$seq2 == "none") seq2 <- input$seq1
    #     print(c(input$seed1, seq2))
    #     callModule(
    #         createArchitecturePlot, "archiPlot",
    #         pointInfo = reactive(c(input$seed1, seq2)),
    #         domainInfo = filterDomainData,
    #         currentNCBIinfo = reactive(currentNCBIinfo),
    #         font = reactive(input$font)
    #     )
    # })
    
    output$domainPlot <- renderPlot({
        req(filterDomainData())
        if (input$doPlot > 0) {
            if (is.null(filterDomainData())) {
                msgPlot()
            } else {
                seq2 <- input$seq2
                if (input$seq2 == "none") seq2 <- input$seq1
                g <- PhyloProfile::createArchiPlot(
                    c(input$seed1, seq2), filterDomainData(),
                    input$labelArchiSize, input$titleArchiSize, input$showScore,
                    input$showWeight, input$namePosition, input$firstDist,input$nameType,
                    input$nameSize, input$segmentSize, input$nameColor, input$labelPos,
                    input$colorType, input$ignoreInstanceNo, currentNCBIinfo,
                    input$featureClassSort, input$featureClassOrder, input$colorPallete,
                    input$resolveOverlap, input$font
                )
                if (any(g == "No domain info available!")) {
                    msgPlot()
                } else {
                    suppressWarnings(grid::grid.draw(g))
                }
            }
        }
    })

    output$domainPlot.ui <- renderUI({
        req(filterDomainData())
        if (is.null(filterDomainData())) {
            msg <- paste0(
                "<p><em>No domain found for this protein.</em></p>"
            )
            HTML(msg)
        } else {
            plotOutput(
                "domainPlot",
                height = input$archiHeight,
                width = input$archiWidth,
                click = "plot_click",
            )
        }
    })
    
    # output$hover_info <- renderPrint({
    #     cat("input$plot_click:\n")
    #     str(input$plot_click)
    # })
    
    output$linkTable <- renderTable({
        req(filterDomainData())
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

    output$domainTable <- DT::renderDataTable({
        req(filterDomainData())
        req(input$seq1)
        req(input$seq2)
        outDf <- getDomainInformation()
        outDf <- outDf[
            ,c(
                "orthoID", "length", "feature_id", "feature_type", "start", "end", "evalue",
                "bitscore", "pStart", "pEnd", "pLen"
            )
        ]
        outDf$orthoID <- gsub(":", "\\|", outDf$orthoID)
        colnames(outDf) <- c(
            "Sequence ID", "Length", "Feature", "Type", "Start", "End", "E-value",
            "Bit-score", "pHMM start", "pHMM end", "pHMM length"
        )
        outDf
    })

    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.svg")
        },
        content = function(file) {
            seq2 <- input$seq2
            if (input$seq2 == "none") seq2 <- input$seq1
            g <- PhyloProfile::createArchiPlot(
                c(input$seed1, seq2), filterDomainData(),
                input$labelArchiSize, input$titleArchiSize, input$showScore,
                input$showWeight, input$namePosition, input$firstDist,input$nameType,
                input$nameSize, input$segmentSize, input$nameColor, input$labelPos,
                input$colorType, input$ignoreInstanceNo, currentNCBIinfo,
                input$featureClassSort, input$featureClassOrder, input$colorPallete,
                input$resolveOverlap, input$font
            )
            suppressWarnings(ggsave(
                file, plot = g,
                width = input$archiWidth * 0.056458333,
                height = input$archiHeight * 0.056458333,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            ))
        }
    )
    
    # * Split multi OG domain file =============================================
    getDomainFileIn <- reactive({
        shinyFileChoose(
            input, "domainFileIn", roots = homePath, session = session,
            filetypes = c('', 'domains')
        )
        fileSelected <- parseFilePaths(homePath, input$domainFileIn)
        return(replaceHomeCharacter(as.character(fileSelected$datapath)))
    })
    output$domainFileIn.ui <- renderUI({
        req(getDomainFileIn())
        if (length(getDomainFileIn()) > 0) {
            outString <- getDomainFileIn()
            if (nchar(outString) > 30)
                outString <- paste0(
                    substrLeft(outString, 15), "...", substrRight(outString, 15)
                )
            em(outString)
        }
    })
    
    getDomainDirOut <- reactive({
        shinyDirChoose(
            input, "domainDirOut", roots = homePath, session = session
        )
        domainPathOut <- parseDirPath(homePath, input$domainDirOut)
        return(replaceHomeCharacter(as.character(domainPathOut)))
    })
    
    output$domainDirOut.ui <- renderUI({
        req(getDomainDirOut())
        if (length(getDomainDirOut()) > 0) {
            em(paste("Save output files to", getDomainDirOut()))
        }
    })
    
    observeEvent(input$doSplitDomain, {
        withCallingHandlers({
            shinyjs::html("splitDomainFileStatus", "")
            splitDomainFile(getDomainFileIn(), getDomainDirOut())
        },
        message = function(m) {
            shinyjs::html(
                id = "splitDomainFileStatus", html = m$message, add = TRUE
            )
        })
        updateButton(session, "doSplitDomain", disabled = TRUE)
    })
})