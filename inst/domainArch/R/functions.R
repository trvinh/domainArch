replaceHomeCharacter <- function (fullPath = NULL) {
    homeName <- system("echo $HOME", intern = TRUE)
    stringr::str_replace(fullPath, "~", homeName)
}

getFileName <- function (filePath = NULL) {
    if (is.null(filePath)) stop ("No complete path given!")
    tmp <- strsplit(filePath, "/")[[1]]
    return(tmp[length(tmp)])
}

substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n) {
    substr(x, 1, n)
}

createPlotSize <- function(id, title, value) {
    numericInput(id,
                 title,
                 min = 100,
                 max = 3200,
                 step = 50,
                 value = value,
                 width = 100)
}

createTextSize <- function(id, title, value, width) {
    numericInput(id,
                 title,
                 min = 3,
                 max = 99,
                 step = 1,
                 value = value,
                 width = width)
}


msgPlot <- function() {
    msg <- paste(
        "No information about domain architecture!",
        "Please check:","if the domain file exists; or ",
        "if the selected genes (seed & ortholog) do exist in the domain file",
        sep = "\n"
    )
    x <- c(1,2,3,4,5)
    y <- c(1,2,3,4,5)
    g <- ggplot(data.frame(x, y), aes(x,y)) +
        geom_point(color = "white") +
        annotate(
            "text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red"
        ) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank()) +
        ylim(0,1)
    return(g)
}

getGroupIds <- function (
        inputType = "File", domainFile = NULL, domainDir = NULL
) {
    if(inputType == "File") {
        df <- read.csv(
            domainFile, header = FALSE, sep = "\t", 
            stringsAsFactors = FALSE
        )
        df[c("groupID", "tmp")] <- str_split_fixed(df$V1, '#', 2)
        return(levels(as.factor(df$groupID)))
        #     list(
        #         levels(as.factor(df$groupID)), 
        #         #levels(as.factor(df$V2[df$groupID == input$seed]))
        #         levels(as.factor(df$groupID))
        #     )
        # )
    } else if (inputType == "Folder") {
        files <- list.files(domainDir, pattern = ".domains")
        return(str_replace(levels(as.factor(files)), ".domains", ""))
        #     list(
        #         str_replace(levels(as.factor(files)), ".domains", ""),
        #         NULL
        #     )
        # )
    }
}

getOrthoIDs <- function (groupID = NULL, file = NULL) {
    if (is.null(groupID) || is.null(file)) return(NULL)
    df <- read.csv(
        file, header = FALSE, sep = "\t", 
        stringsAsFactors = FALSE
    )
    df[c("groupID", "tmp")] <- str_split_fixed(df$V1, '#', 2)
    return(levels(as.factor(df$V2[df$groupID == groupID])))
}

getSpecList <- function (annoDir = NULL) {
    if (is.null(annoDir)) stop("No annotation folder given!")
    files <- list.files(annoDir, pattern = ".json")
    return(str_replace(levels(as.factor(files)), ".json", ""))
}

json2list <- function(jsonFile = NULL) {
    if (is.null(jsonFile)) return(NULL)
    jsonList <- fromJSON(jsonFile)
    return(jsonList)
}

parseDomainFromJson <- function (spec = NULL, protIds = NULL, jsonList = NULL) {
    if (is.null(protIds)) stop("No protein ID given")
    if (is.null(spec)) stop("No species ID given")
    outList <- lapply(
        protIds,
        function (id) {
            featList <- lapply(
                c("flps", "tmhmm", "signalp", "coils2", "seg", "smart", "pfam"), 
                function (x) {
                    if (length(jsonList$feature[[id]][[x]]) > 0) {
                        instanceList <- lapply(
                            names(jsonList$feature[[id]][[x]]), 
                            function (y) {
                                return(
                                    c(
                                        y, 
                                        jsonList$feature[[id]][[x]][[y]]$instance[1:2]
                                    )
                                )
                            }
                        )
                        return(data.frame(do.call(rbind, instanceList)))
                    }
                }
            )
            tmpDf <- data.frame(do.call(rbind, featList))
            if (nrow(tmpDf) > 0) {
                tmpDf$length <- jsonList$feature[[id]]$length
                tmpDf$protId <- id
                return(tmpDf)
            }
        }
    )
    outDf <- data.frame(do.call(rbind, outList))
    if (nrow(outDf) > 0) {
        colnames(outDf) <- c("feature", "start", "end", "length", "orthoID")
        outDf$seedID <- paste0(spec,"#",outDf$orthoID)
        outDf$weight <- 0
        outDf$path <- "N"
        outDf$start <- as.integer(outDf$start)
        outDf$end <- as.integer(outDf$end)
        return(outDf)
    } else {
        message("No domain found for ", protIds, " in ", spec)
        return(NULL)
    }
    
}


createArchiPlot2 <- function(
        info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12,
        showScore = NULL, showName = "plot"
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    # subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    subdomainDf <- domainDf
    subdomainDf$feature <- as.character(subdomainDf$feature)
    orthoID <- NULL
    
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # get minStart and maxEnd
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))
        
        if (nrow(orthoDf) > 0) {
            # change order of one df's features based on order of other df's
            if (length(orthoDf$feature) < length(seedDf$feature)) {
                orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
                orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
                orderedOrthoDf <- sortDomains(orderedSeedDf, orderedOrthoDf)
            } else {
                orderedSeedDf <- seedDf[order(seedDf$feature), ]
                orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
                orderedSeedDf <- sortDomains(orderedOrthoDf, orderedSeedDf)
            }
            # join weight values and feature names
            if ("weight" %in% colnames(orderedOrthoDf)) {
                NULL
                # orderedOrthoDf$yLabel <- paste0(
                #     orderedOrthoDf$feature," (",round(orderedOrthoDf$weight, 2),")")
            } else orderedOrthoDf$yLabel <- orderedOrthoDf$feature
            if ("weight" %in% colnames(orderedSeedDf)) {
                NULL
                # orderedSeedDf$yLabel <- paste0(
                #     orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
            } else orderedSeedDf$yLabel <- orderedSeedDf$feature
            # plotting
            g <- pairDomainPlotting(
                seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, showScore, showName)
        } else {
            orderedSeedDf <- seedDf[order(seedDf$feature), ]
            if ("weight" %in% colnames(orderedSeedDf)) {
                orderedSeedDf$yLabel <- paste0(
                    orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
            } else orderedSeedDf$yLabel <- orderedSeedDf$feature
            # plotting
            g <- pairDomainPlotting(
                seed, seed, orderedSeedDf, orderedSeedDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, showScore, showName)
        }
        return(g)
    }
}

singleDomainPlotting <- function(
        df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
        minStart = NULL, maxEnd = NULL, colorScheme = NULL, 
        showScore = NULL, showName = "plot"
){
    feature <- feature_id_mod <- end <- start <- NULL
    
    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if (is.null(colorScheme)) {
        colorScheme <- structure(
            getQualColForVector(levels(as.factor(df$feature))),
            .Names = levels(as.factor(df$feature)))}
    # initiate ggplot object
    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        scale_color_manual(values = colorScheme)
    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, size = 1, color = "#b2b2b2", alpha = 0.0,
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw features
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature),
        size = 3)
    # add feature names
    if ("plot" %in% showName) {
        gg <- gg + geom_label(
            aes(label = str_wrap(feature_id_mod),
                x = (start+end)/2),
            color = "black", vjust = -0.25
        )
    }
    # add scores if selected
    if ("Bit-score" %in% showScore & "E-value" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(evalue == "NA" & bitscore == "NA", "", str_wrap(
                    paste0(
                        "E-value: ", evalue, "; Bitscore: ", bitscore
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    } else if ("Bit-score" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(bitscore == "NA", "", str_wrap(
                    paste0(
                        "Bitscore: ", bitscore
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    } else if ("E-value" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(evalue == "NA", "", str_wrap(
                    paste0(
                        "E-value: ", evalue
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    }
    # theme format
    gg <- gg + labs(
        title = paste0(gsub(":", sep, geneID)), color = "Feature"
    )
    gg <- gg + theme_minimal() + theme(panel.border = element_blank())
    # gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(plot.title = element_text(face = "bold", size = titleSize))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
    # add feature names on the axis (if required)
    if ("axis" %in% showName) {
        if ("plot" %in% showName) {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature, labels = df$feature_type)
        } else {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature, labels = df$feature)
        }
        gg <- gg + theme(axis.text.y = element_text(size = labelSize))
    } else {
        gg <- gg + theme(axis.text.y = element_blank())
    }
    # add legend (if required)
    if ("legend" %in% showName) {
        gg <- gg + theme(legend.position = "bottom")
    } else {
        gg <- gg + theme(legend.position = "none")
    }
    # add space on the top of the plot (for feature name)
    gg <- gg + coord_cartesian(
        clip = 'off', ylim = c(1, nlevels(as.factor(df$feature)) + 0.5)
    )
    return(gg)
}

pairDomainPlotting <- function(
        seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
        minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12,
        showScore = NULL, showName = "plot"
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")
    # create color scheme, so that the same features in seed & ortholog will
    # have the same colors
    featureSeed <- levels(as.factor(seedDf$feature))
    featureOrtho <- levels(as.factor(orthoDf$feature))
    allFeatures <- c(featureSeed, featureOrtho)
    allColors <- getQualColForVector(allFeatures)
    colorScheme <- structure(allColors, .Names = allFeatures)
    # plot
    sep <- "|"
    plotOrtho <- singleDomainPlotting(
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd,colorScheme, showScore, showName)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorScheme, showScore, showName)
    if (ortho == seed) {
        g <- plotSeed
    } else {
        seedHeight <- length(levels(as.factor(seedDf$feature)))
        orthoHeight <- length(levels(as.factor(orthoDf$feature)))
        if ("legend" %in% showName) {
            g <- grid_arrange_shared_legend(
                plotSeed, plotOrtho,
                ncol = 1, nrow = 2, position = "bottom"
            )
        } else {
            g <- gridExtra::arrangeGrob(plotSeed, plotOrtho, ncol = 1, nrow = 2)
        }
    }
    return(g)
}

sortDomains <- function(seedDf, orthoDf){
    if (is.null(seedDf) | is.null(orthoDf))
        stop("Domain data for seed & ortholog cannot be NULL!")
    orderNo <- NULL
    # get list of features in seedDf
    featureList <- as.data.frame(levels(as.factor(seedDf$feature)))
    colnames(featureList) <- c("feature")
    # and add order number to each feature
    featureList$orderNo <- seq(length(featureList$feature))

    # merge those info to orthoDf
    orderedOrthoDf <- merge(orthoDf, featureList, all.x = TRUE)

    # sort orthoDf
    index <- with(orderedOrthoDf, order(orderNo))
    orderedOrthoDf <- orderedOrthoDf[index, ]
    
    #turn feature column into a character vector
    orderedOrthoDf$feature <- as.character(orderedOrthoDf$feature)
    #then turn it back into an ordered factor (to keep this order when plotting)
    orderedOrthoDf$feature <- factor(
        orderedOrthoDf$feature, levels = unique(orderedOrthoDf$feature)
    )
    #return sorted df
    return(orderedOrthoDf)
}


#' Join multiple plots and merge legends
grid_arrange_shared_legend <- function (
    ..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")
) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
        position,
        "bottom" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 1,
            heights = grid::unit.c(unit(1, "npc") - lheight, lheight)
        ),
        "right" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 2,
            widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)
        )
    )
    
    return(combined)
}

#' get pfam and smart domain links
#' @return dataframe with domain IDs and their database links
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainLink <- function(info, domainDf) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    subdomainDf$feature <- as.character(subdomainDf$feature)
    orthoID <- NULL
    feature <- NULL
    if (nrow(subdomainDf) < 1) {
        return(data.frame("ID" = character(), "PFAM" = character(), "SMART"= character()))
    }
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        feature <- c(
            levels(as.factor(orthoDf$feature)),
            levels(as.factor(seedDf$feature))
        )
    }
    # get URLs
    featurePfam <- unique(feature[grep("pfam", feature)])
    pfamDf <- data.frame(ID = character(), PFAM = character())
    if (length(featurePfam) > 0)
        pfamDf <- createLinkTable(featurePfam, "pfam")
    
    featureSmart <- unique(feature[grep("smart", feature)])
    smartDf <- data.frame(ID = character(), SMART = character())
    if (length(featureSmart) > 0)
        smartDf <- createLinkTable(featureSmart, "smart")
    
    featDf <- merge(pfamDf, smartDf, by = "ID", all = TRUE)
    colnames(featDf) <- c("ID", "PFAM", "SMART")
    return(featDf)
}

#' plot error message
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
createLinkTable <- function(featureList, featureType) {
    feature <- sub("_","@", featureList)
    featDf <- NULL
    if (length(feature) > 0) {
        tmpDf <- data.frame(
            do.call(
                'cbind',
                data.table::tstrsplit(as.character(feature), '@', fixed = TRUE)
            )
        )
        
        featDf <- data.frame("ID" = levels(as.factor(tmpDf$X2)))
        if (featureType == "pfam") {
            # featDf$type <- "PFAM"
            featDf$link <- paste0(
                "<a href='https://pfam.xfam.org/family/", featDf$ID,
                "' target='_blank'>", featDf$ID, "</a>"
            )
        } else {
            # featDf$type <- "SMART"
            featDf$link <- paste0(
                "<a href='http://smart.embl-heidelberg.de/smart/",
                "do_annotation.pl?BLAST=DUMMY&DOMAIN=",
                featDf$ID, "' target='_blank'>",
                featDf$ID, "</a>"
            )
        }
    }
    return(featDf)
}