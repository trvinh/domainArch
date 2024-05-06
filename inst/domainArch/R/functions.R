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

createPlotSize <- function(id, title, value, width) {
    numericInput(id,
                 title,
                 min = 100,
                 max = 3200,
                 step = 50,
                 value = value,
                 width = width)
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

getOrthoIDs <- function (groupID = NULL, file = NULL, currentNCBIinfo = NULL) {
    if (is.null(groupID) || is.null(file)) return(NULL)
    df <- read.csv(
        file, header = FALSE, sep = "\t",
        stringsAsFactors = FALSE
    )
    df[c("groupID", "tmp")] <- str_split_fixed(df$V1, '#', 2)
    idDf <- data.frame(oriID = levels(as.factor(df$V2[df$groupID == groupID])))
    idDf$seqID <- sapply(str_split(idDf$oriID, "\\|"), "[", 3)
    idDf$ncbiID <- sapply(str_split(idDf$oriID, "@"), "[", 2)
    if (is.null(currentNCBIinfo)) {
        return(setNames(idDf$oriID, idDf$oriID))
    } else {
        id2nameDf <- PhyloProfile::id2name(idDf$ncbiID, currentNCBIinfo)
        idDf <- merge(idDf, id2nameDf, by = "ncbiID", all.x = TRUE)
        return(setNames(idDf$oriID, paste(idDf$fullName, idDf$seqID, sep = " - ")))
    }
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
        showScore = NULL, showName = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000", labelPos = "Above",
        colorType = "Unique", ignoreInstanceNo = FALSE, currentNCBIinfo = NULL,
        featureTypeSort = "Yes", featureTypeOrder = NULL, colorPallete = "Paired",
        resolveOverlap = "Yes"
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

        # simplify seed/ortho seq IDs if they are in bionf format
        if (!is.null(currentNCBIinfo)) {
            if (str_count(seed, ":") >= 2 & str_count(seed, "@") >= 2) {
                seedTmp <- strsplit(as.character(seed),':', fixed = TRUE)[[1]]
                seedSpec <-
                    strsplit(as.character(seedTmp[2]),'@', fixed = TRUE)[[1]][2]
                seed <- paste0(
                    id2name(seedSpec, currentNCBIinfo)[,2], " - ", seedTmp[3]
                )
                if (ortho != seed) {
                    orthoTmp <- strsplit(as.character(ortho),':', fixed = TRUE)[[1]]
                    orthoSpec <-
                        strsplit(as.character(orthoTmp[2]),'@',fixed = TRUE)[[1]][2]
                    ortho <- paste0(
                        id2name(orthoSpec, currentNCBIinfo)[,2], " - ", orthoTmp[3]
                    )
                }
            }
        }

        # add feature colors
        featureColorDf <- addFeatureColors(seedDf, orthoDf, colorType, colorPallete, ignoreInstanceNo)
        seedDf <- featureColorDf[[1]]
        orthoDf <- featureColorDf[[2]]
        # resolve (non)overlapped features
        if (resolveOverlap == "Yes") {
            seedDf <- resolveOverlapFeatures(seedDf)
            orthoDf <- resolveOverlapFeatures(orthoDf)
        } else {
            seedDf$featureOri <- seedDf$feature
            seedDf$featureOri <- as.character(seedDf$featureOri)
            orthoDf$featureOri <- orthoDf$feature
            orthoDf$featureOri <- as.character(orthoDf$featureOri)
        }

        if (nrow(orthoDf) > 0) {
            if (all.equal(seedDf, orthoDf)[1] == TRUE) featureTypeSort <- "No"
            # sort features
            if (featureTypeSort == "Yes") {
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
            } else {
                # change order based on list of feature types
                orderedSeedDf <- sortDomainsByList(seedDf, featureTypeOrder)
                orderedOrthoDf <- sortDomainsByList(orthoDf, featureTypeOrder)
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
                labelArchiSize, titleArchiSize, showScore, showName, firstDist,
                nameType, nameSize, nameColor, labelPos, colorPallete)
        } else {
            # orderedSeedDf <- seedDf[order(seedDf$feature), ]
            orderedSeedDf <- sortDomainsByList(seedDf, featureTypeOrder)
            if ("weight" %in% colnames(orderedSeedDf)) {
                NULL
                # orderedSeedDf$yLabel <- paste0(
                #     orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
            } else orderedSeedDf$yLabel <- orderedSeedDf$feature
            # plotting
            g <- pairDomainPlotting(
                seed, seed, orderedSeedDf, orderedSeedDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, showScore, showName, firstDist,
                nameType, nameSize, nameColor, labelPos, colorPallete)
        }
        return(g)
    }
}

singleDomainPlotting <- function(
        df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
        minStart = NULL, maxEnd = NULL, colorPallete = "Set2",
        showScore = NULL, showName = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000", labelPos = "Above"
){
    feature <- feature_id_mod <- end <- start <- NULL

    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if ("color" %in% colnames(df)) {
        colorScheme <- structure(
            df$color, .Names = df$featureOri
        )
    } else {
        colorScheme <- structure(
            head(
                suppressWarnings(
                    RColorBrewer::brewer.pal(
                        nlevels(as.factor(df$feature)), colorPallete
                    )
                ),
                levels(as.factor(df$feature))
            ),
            .Names = levels(as.factor(df$feature)))
    }

    # initiate ggplot object
    gg <- ggplot(df, aes(y = feature, x = end))

    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, linewidth = 1, color = "#b2b2b2", alpha = 0.0,
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw features
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature, color = as.factor(featureOri)),
        linewidth = nameSize, lineend = "round", alpha = 0.7) +
        scale_color_manual(values = colorScheme)

    # add feature names
    if ("plot" %in% showName) {
        if (nameType == "Labels") {
            if (labelPos == "Above") {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", vjust = -0.5
                )
            } else if (labelPos == "Below") {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", vjust = 1.5
                )
            } else {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", size = nameSize - 2
                )
            }
        } else {
            gg <- gg + geom_text(
                aes(label = str_wrap(feature_id_mod),
                    x = (start+end)/2),
                color = nameColor, check_overlap = TRUE, size = nameSize - 2
            )
        }
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
        if ("plot" %in% showName | "legend" %in% showName) {
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
        clip = 'off', ylim = c(1, nlevels(as.factor(df$feature)) + firstDist)
    )
    return(gg)
}

pairDomainPlotting <- function(
        seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
        minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12,
        showScore = NULL, showName = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000", labelPos = "Above",
        colorPallete = "Paired"
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")

    sep <- "|"
    plotOrtho <- singleDomainPlotting(
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd, colorPallete,
        showScore, showName, firstDist, nameType, nameSize, nameColor, labelPos)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorPallete,
        showScore, showName, firstDist, nameType, nameSize, nameColor, labelPos)
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

sortDomainsByList <- function(domainDf, featureTypeOrder) {
    if (is.null(domainDf) | is.null(featureTypeOrder))
        stop("Domain data or feature type order is NULL!")
    featureTypeOrder <- rev(featureTypeOrder[
        featureTypeOrder %in% levels(as.factor(domainDf$feature_type))
    ])
    orderedDomainDf <- left_join(
        data.frame(feature_type = featureTypeOrder),
        domainDf, by = "feature_type"
    )
    orderedDomainDf$feature <- factor(
        orderedDomainDf$feature, levels = unique(orderedDomainDf$feature)
    )
    return(orderedDomainDf)
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

#' Check if a color pallete has enough colors for a list of items
checkColorPallete <- function(items, pallete) {
    colorDf <- data.frame(RColorBrewer::brewer.pal.info)
    colorDf$name <- row.names(colorDf)
    if (length(items) < colorDf$maxcolors[colorDf$name == pallete])
        return(TRUE)
    return(FALSE)
}

#' Identify feature type(s) containing overlapped features
checkOverlapDomains <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    df <- domainDf[, c("feature", "feature_type", "start", "end")]
    df <- data.frame(df %>% group_by(feature_type) %>% add_count(feature_type))
    df$tmp <- paste(df$feature_type, df$end, sep = "_")
    overlappedType <- lapply(
        df$tmp[df$n > 1],
        function (x) {
            if (length(df$feature[df$tmp == x]) > 1)
                return(df$feature_type[df$tmp == x])
            ed <- df$end[df$tmp == x][1]
            st <- df$start[df$tmp == x][1]
            type <- df$feature_type[df$tmp == x][1]

            subDf <- df[
                !(df$tmp == x) & df$feature_type == type & df$start < ed & df$end > st,
            ]
            if (nrow(subDf) > 0) {
                if(!(subDf$feature[1] == df$feature[df$tmp == x]))
                    return(subDf$feature_type)
            }
        }
    )
    return(unique(unlist(overlappedType)))
}

#' Modify domain df to resolve overlapped features
resolveOverlapFeatures <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    overlappedType <- checkOverlapDomains(domainDf)
    domainDf$featureOri <- domainDf$feature
    domainDf$featureOri <- as.character(domainDf$featureOri)
    if (length(overlappedType) > 0) {
        domainDf$feature <- ifelse(
            !(domainDf$feature_type %in% overlappedType),
            domainDf$feature_type,
            domainDf$feature
        )
    } else {
        domainDf$feature <- domainDf$feature_type
    }
    return(domainDf)
}


#' Create color scheme, so that the same features in seed & ortholog will
#' have the same colors

addFeatureColors <- function(
        seedDf = NULL, orthoDf = NULL, colorType = "all",
        colorPallete = "Paired", ignoreInstanceNo = FALSE
) {
    if (is.null(seedDf) | is.null(orthoDf)) stop("Domain Df cannot be null!")

    featureSeedCount <- seedDf %>% dplyr::count(feature)
    featureOrthoCount <- orthoDf %>% dplyr::count(feature)
    featureCount <- merge(featureSeedCount, featureOrthoCount, by = "feature", all = TRUE)
    featureCount$type <- "unique"
    featureCount$type[featureCount$n.x == featureCount$n.y] <- "shared"
    if (ignoreInstanceNo == TRUE)
        featureCount$type[!is.na(featureCount$n.x) & !is.na(featureCount$n.y)] <- "shared"

    featureCount$feature <- as.character(featureCount$feature)
    sharedFeatures <- unique(featureCount$feature[featureCount$type == "shared"])
    uniqueFeatures <- unique(featureCount$feature[featureCount$type == "unique"])
    allFeatures <- c(sharedFeatures, uniqueFeatures)

    if (colorType == "Unique" & length(uniqueFeatures) == 0)
        colorType = "All"
    if (colorType == "Shared" & length(sharedFeatures) == 0)
        colorType = "All"

    if (colorType == "Unique") {
        sharedFeaturesColors <- rep("#C9C9C9", length(sharedFeatures))
        uniqueFeaturesColors <- getQualColForVector(uniqueFeatures)
        if (checkColorPallete(uniqueFeatures, colorPallete) == TRUE) {
            uniqueFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(uniqueFeatures), colorPallete),
                    length(uniqueFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "Shared") {
        sharedFeaturesColors <- getQualColForVector(sharedFeatures)
        uniqueFeaturesColors <- rep("#C9C9C9", length(uniqueFeatures))
        if (checkColorPallete(sharedFeatures, colorPallete) == TRUE) {
            sharedFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(sharedFeatures), colorPallete),
                    length(sharedFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "All") {
        allColors <- getQualColForVector(allFeatures)
        if (checkColorPallete(allFeatures, colorPallete) == TRUE) {
            allColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(allFeatures), colorPallete),
                    length(allFeatures)
                ))
        }
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else {
        tmpDf <- data.frame(str_split_fixed(allFeatures, "_", 2))
        tmpDf$name <- paste(tmpDf$X1, tmpDf$X2, sep = "_")
        tmpDf$name[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        tmpDf$X2[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        typeColorDf <- data.frame(
            colors = head(
                suppressWarnings(RColorBrewer::brewer.pal(nlevels(as.factor(tmpDf$X1)), colorPallete)),
                nlevels(as.factor(tmpDf$X1))
            ),
            X1 = levels(as.factor(tmpDf$X1))
        )
        tmpDf <- merge(tmpDf, typeColorDf, all.x = TRUE)
        colorScheme <- data.frame(color = tmpDf$colors, feature = tmpDf$name)
    }

    # add color to seedDf and orthoDf
    seedDf <- merge(seedDf, colorScheme, by = "feature", all.x = TRUE)
    orthoDf <- merge(orthoDf, colorScheme, by = "feature", all.x = TRUE)

    return(list(seedDf, orthoDf))
}

#' Split a multi ortholo group file into single files
splitDomainFile <- function(domainFile = NULL, outPath = NULL) {
    if (is.null(domainFile)) stop("Domain file cannot be NULL")
    if (is.null(outPath)) stop("Output path cannot be NULL")

    df <- fread(
        domainFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"
    )
    setDT(df)[, paste0("# pairID", 1:2) := tstrsplit(`# pairID`, "#")]

    dir.create(file.path(outPath), showWarnings = FALSE)
    outList <- lapply(
        levels(as.factor(df$`# pairID1`)),
        function (x) {
            outDf <- df %>% filter(`# pairID1` == x) %>% select(
                `# pairID`, orthoID, seqLen, feature, fStart, fEnd, fWeight, fPath,
                interProID, `e-value`, bitScore, pStart, pEnd, pLen
            )
            write.table(
                outDf, file = paste0(outPath,"/",x,".domains"), quote = FALSE,
                row.names = FALSE, sep = "\t"
            )
            return(x)
        }
    )
    outList <- unlist(outList)
    outFiles <- gsub(".domains", "", list.files(outPath))
    if (length(setdiff(outList,outFiles)) == 0) {
        message(paste0(length(outList), " following files have been saved in ", outPath,"/"))
        message(paste(paste(head(outList, 10), ".domains", sep = ""), collapse = "\n"))
        if (length(outList) > 10) message("... (and some more)")
        message("Done!")
    } else {
        message("Something went wrong...")
        message("Please check the R terminal and the output folder!")
    }

}
