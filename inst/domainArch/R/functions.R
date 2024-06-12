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
    } else if (inputType == "Folder") {
        files <- list.files(domainDir, pattern = ".domains")
        return(str_replace(levels(as.factor(files)), ".domains", ""))
    }
}

getOrthoIDs <- function (groupID = NULL, file = NULL, currentNCBIinfo = NULL) {
    if (is.null(groupID) || is.null(file)) return(NULL)
    df <- read.csv(
        file, header = FALSE, sep = "\t",
        stringsAsFactors = FALSE
    )
    df[c("groupID", "tmp")] <- str_split_fixed(df$V1, '#', 2)
    
    idDf <- data.frame(oriID = levels(as.factor(df$tmp[df$groupID == groupID])))
    idDf$seqID <- sapply(str_split(idDf$oriID, "\\|"), "[", 3)
    idDf$ncbiID <- sapply(str_split(idDf$oriID, "@"), "[", 2)
    if (nrow(idDf) > 0) {
        if (is.null(currentNCBIinfo)) {
            return(setNames(idDf$oriID, idDf$oriID))
        } else {
            id2nameDf <- PhyloProfile::id2name(idDf$ncbiID, currentNCBIinfo)
            idDf <- merge(idDf, id2nameDf, by = "ncbiID", all.x = TRUE)
            return(setNames(idDf$oriID, paste(idDf$fullName, idDf$seqID, sep = " - ")))
        }
    } else return(NULL)
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

#' Split a multi ortholo group file into single files
splitDomainFile <- function(domainFile = NULL, outPath = NULL) {
    if (is.null(domainFile)) stop("Domain file cannot be NULL")
    if (is.null(outPath)) stop("Output path cannot be NULL")

    df <- fread(
        domainFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"
    )
    setDT(df)[, paste0("# pairID", 1:2) := tstrsplit(`# pairID`, "#")]

    dir.create(file.path(outPath), showWarnings = FALSE)
    message("Please wait...")
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
