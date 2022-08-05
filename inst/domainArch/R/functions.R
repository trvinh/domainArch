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

getGroupIds <- function (inputType, domainFile, domainDir) {
    if(inputType == "File") {
        df <- read.csv(
            domainFile, header = FALSE, sep = "\t", 
            stringsAsFactors = FALSE
        )
        df[c("groupID", "tmp")] <- str_split_fixed(df$V1, '#', 2)
        return(
            list(
                levels(as.factor(df$groupID)), 
                #levels(as.factor(df$V2[df$groupID == input$seed]))
                levels(as.factor(df$groupID))
            )
        )
    } else if (inputType == "Folder") {
        files <- list.files(domainDir, pattern = ".domains")
        return(
            list(
                str_replace(levels(as.factor(files)), ".domains", ""),
                NULL
            )
        )
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


createArchiPlot2 <- function(
        info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12
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
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))
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
            orderedOrthoDf$yLabel <- paste0(
                orderedOrthoDf$feature," (",round(orderedOrthoDf$weight, 2),")")
        } else orderedOrthoDf$yLabel <- orderedOrthoDf$feature
        if ("weight" %in% colnames(orderedSeedDf)) {
            orderedSeedDf$yLabel <- paste0(
                orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
        } else orderedSeedDf$yLabel <- orderedSeedDf$feature
        # plotting
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        g <- pairDomainPlotting(
            seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
            labelArchiSize, titleArchiSize)
        return(g)
    }
}

singleDomainPlotting <- function(
        df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
        minStart = NULL, maxEnd = NULL, colorScheme = NULL
){
    feature <- end <- start <- NULL
    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if (is.null(colorScheme)) {
        colorScheme <- structure(
            getQualColForVector(levels(as.factor(df$feature))),
            .Names = levels(as.factor(df$feature)))}
    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        geom_segment(
            data = df, color = "white", size = 0,
            aes(y = feature, yend = feature, x = minStart, xend = maxEnd)) +
        scale_color_manual(values = colorScheme)
    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, size = 1, color = "#b2b2b2",
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw line and points
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature),
        size = 1.5) +
        geom_point(data = df, aes(y = feature, x = start),
                   color = "#b2b2b2", size = 3, shape = 3) +
        geom_point(data = df, aes(y = feature, x = end),
                   color = "#edae52", size = 3, shape = 5)
    # draw dashed line for domain path
    gg <- gg + geom_segment(
        data = df[df$path == "Y", ], size = 3, linetype = "dashed",
        aes(x = start, xend = end, y = feature, yend = feature))
    # theme format
    gg <- gg + scale_y_discrete(
        expand = c(0.075, 0), breaks = df$feature, labels = df$yLabel)
    gg <- gg + labs(title = paste0(gsub(":", sep, geneID)), y = "Feature")
    gg <- gg + theme_minimal() + theme(panel.border = element_blank())
    gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(plot.title = element_text(face = "bold", size = titleSize))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_text(size = labelSize),
        axis.title.y = element_text(size = labelSize),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
    return(gg)
}

pairDomainPlotting <- function(
        seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
        minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12
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
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd,colorScheme)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorScheme)
    if (ortho == seed) {
        g <- gridExtra::arrangeGrob(plotSeed, ncol = 1)
    } else {
        seedHeight <- length(levels(as.factor(seedDf$feature)))
        orthoHeight <- length(levels(as.factor(orthoDf$feature)))
        g <- gridExtra::arrangeGrob(
            plotSeed, plotOrtho, ncol = 1, heights = c(seedHeight, orthoHeight)
        )
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
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
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
