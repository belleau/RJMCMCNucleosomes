#' @title Nucleosome positioning mapping on a segment
#'
#' @description Use of a fully Bayesian hierarchical model for chromosome-wide
#' profiling of nucleosome positions based on high-throughput short-read
#' data (MNase-Seq data). Beware that for a genome-wide profiling, each
#' chromosome must be treated separatly. This function is optimized to run
#' on segments that are smaller sections of the chromosome.
#'
#' @param reads a \code{GRanges} containing forward and
#' reverse reads. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
#' Default: \code{NULL}.
#'
#' @param nbrIterations a positive \code{integer} or \code{numeric}, the
#' number of iterations. Non-integer values of
#' \code{nbrIterations} will be casted to \code{integer} and truncated towards
#' zero.
#'
#' @param kMax a positive \code{integer} or \code{numeric}, the maximum number
#' of degrees of freedom per region. Non-integer values
#' of \code{kMax} will be casted to \code{integer} and truncated towards zero.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution. Default: 3.
#'
#' @param minInterval a \code{numeric}, the minimum distance between two
#' nucleosomes.
#'
#' @param maxInterval a \code{numeric}, the maximum distance between two
#' nucleosomes.
#'
#' @param minReads a positive \code{integer} or \code{numeric}, the minimum
#' number of reads in a potential canditate region. Non-integer values
#' of \code{minReads} will be casted to \code{integer} and truncated towards
#' zero. Default: 5.
#'
#' @param adaptIterationsToReads a \code{logical} indicating if the number
#' of iterations must be modified in function of the number of reads.
#' Default: \code{TRUE}.
#'
#' @param vSeed a \code{integer}. A seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: -1.
#'
#' @param saveAsRDS a \code{logical}. When \code{TRUE}, a RDS file containing
#' the complete output of the c++ rjmcmc() function is created.
#' Default : \code{FALSE}.
#'
#' @return a \code{list} of \code{class} "rjmcmcNucleosomes" containing:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{k} a \code{integer}, the final estimation of the number
#' of nucleosomes. \code{0} when no nucleosome is detected.
#' \item \code{mu} a \code{GRanges} containing the positions of the
#' nucleosomes and '*' as strand. The \code{seqnames} of the \code{GRanges}
#' correspond to the \code{seqName} input value. \code{NA} when no nucleosome
#' is detected.
#' \item \code{k_max} a \code{integer}, the maximum number of nucleosomes
#' obtained during the iteration process. \code{NA} when no nucleosome is
#' detected.
#' }
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_01)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(reads = reads_demo_01, seqName = "chr_SYNTHETIC",
#'             nbrIterations = 1000, lambda = 2, kMax = 30,
#'             minInterval = 146, maxInterval = 292, minReads = 5,
#'             vSeed = 10113, saveAsRDS = FALSE)
#'
#' ## Print the final estimation of the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' ## Print the maximum number of nucleosomes obtained during the iteration
#' ## process
#' result$k_max
#'
#' @author Rawane Samb, Pascal Belleau, Astrid Deschenes
#' @importFrom stats aggregate
#' @export
rjmcmc <- function(reads, seqName = NULL,
                    nbrIterations, kMax, lambda = 3,
                    minInterval, maxInterval, minReads = 5,
                    adaptIterationsToReads = TRUE, vSeed = -1,
                    saveAsRDS = FALSE) {

    # Get call information
    cl <- match.call()

    # Parameters validation
    validateRJMCMCParameters(reads = reads,
                            seqName = seqName,
                            nbrIterations = nbrIterations,
                            kMax = kMax,
                            lambda = lambda,
                            minInterval = minInterval,
                            maxInterval = maxInterval,
                            minReads = minReads,
                            adaptIterationsToReads = adaptIterationsToReads,
                            vSeed = vSeed)

    resultRJMCMC <- NULL

    if (length(reads) > 0) {

        ## Only keep reads associated to the specified chromosome
        if (!is.null(seqName)) {
            reads <- reads[seqnames(reads) == seqName]
        }

        if (length(reads) > 0) {

            startPosForwardReads <- start(reads[strand(reads) == "+"])

            startPosReverseReads <- end(reads[strand(reads) == "-"])

            # Find nucleosome positions
            if(length(startPosForwardReads) > 0 &
                            length(startPosReverseReads) > 0){
                resultRJMCMC <- rjmcmcNucleo(startPosForwardReads,
                                        startPosReverseReads,
                                        nbrIterations, kMax, lambda,
                                        minInterval, maxInterval, minReads,
                                        adaptIterationsToReads, vSeed)
            }
        }
    }

    # Save output in a RDS file
    if (saveAsRDS) {
        options(digits.secs = 2)
        file_name <- gsub(Sys.time(), pattern = "[:. ]", replacement = "_",
                            perl = TRUE)
        saveRDS(object = resultRJMCMC,
                file = paste0("RJMCMCNucleosomes_output_", file_name, ".RDS"))
    }

    if (is.null(resultRJMCMC)) {
        ## Set values when no nucleosome can be found
        k = 0
        mu = NA
        k_max = NA
    } else {
        # Find k value with the maximum of iterations
        iterPerK <- data.frame(k = resultRJMCMC$k, it = resultRJMCMC$it)
        sumIterPerK <- aggregate(it ~ k, data = iterPerK, sum)
        maxRow <- which.max( sumIterPerK[,"it"])
        k <- sumIterPerK$k[maxRow]

        # Find mu values associated to the k value
        if(is.null(seqName)){
            seqName <- seqnames(reads)[1]
        }
        # Set mu as a GRanges
        mu <- GRanges(seqnames = rep(seqName, k),
                    ranges = IRanges(start=round(resultRJMCMC$muHat[k,][1:k]),
                    end=round(resultRJMCMC$muHat[k,][1:k])),
                    strand = rep('*',k))
        # Get the k_max value
        k_max <- resultRJMCMC$k_max
    }

    # Format output
    result <- list(
        call = cl,
        k = k,
        mu = mu,
        k_max = k_max
    )

    class(result)<-"rjmcmcNucleosomes"

    return(result)
}


#' @title Merge nucleosome information from all RDS files present
#' in a same directory. Beware that only nucleosome information from same
#' chromosome should be merged together.
#'
#' @description Merge nucleosome information, from all RDS files present
#' in a same directory, into one object
#' of \code{class} "rjmcmcNucleosomesMerge".
#'
#' @param directory a \code{character}, the
#' name of the directory (relative or absolute path) containing RDS files. The
#' RDS files must
#' contain R object of \code{class} "rjmcmcNucleosomes" or
#' "rjmcmcNucleosomesMerge".
#'
#' @return a \code{list} of \code{class} "rjmcmcNucleosomesMerge" containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item \code{mu} a \code{GRanges} containing the positions of the
#' nucleosomes.
#' }
#'
#' @examples
#'
#' ## Use a directory present in the RJMCMC package
#' directoryWithRDSFiles <- system.file("extdata",
#' package = "RJMCMCNucleosomes")
#'
#' ## Merge nucleosomes info from RDS files present in directory
#' ## It is assumed that all files present in the directory are nucleosomes
#' ## result for the same chromosome
#' result <- mergeAllRDSFilesFromDirectory(directoryWithRDSFiles)
#'
#' ## Print the number and the position of the nucleosomes
#' result$k
#' result$mu
#'
#' ## Class of the output object
#' class(result)
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
mergeAllRDSFilesFromDirectory <- function(directory) {

    ## Validate that the directory exist
    validateDirectoryParameters(directory)

    ## Get the list of all RDS files present in the directory
    fileList <- dir(directory, pattern = ".rds", full.names = TRUE,
                        ignore.case = TRUE)

    ## Extract information from each file
    return(mergeAllRDSFiles(fileList))
}


#' @title Merge nucleosome information from selected RDS files.
#'
#' @description Merge nucleosome information present in RDS files into one
#' object of \code{class} "rjmcmcNucleosomesMerge".
#'
#' @param RDSFiles a \code{array}, the
#' names of all RDS used to merge nucleosome information. The files must
#' contain R object of \code{class} "rjmcmcNucleosomes" or
#' "rjmcmcNucleosomesMerge".
#'
#' @return a \code{list} of \code{class} "rjmcmcNucleosomesMerge" containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item \code{mu} a \code{GRanges} containing the positions of the
#' nucleosomes.
#' }
#'
#' @examples
#'
#' ## Use RDS files present in the RJMCMC package
#' RDSFiles <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#' full.names = TRUE, pattern = "*RDS")
#'
#' ## Merge nucleosomes info from RDS files present in directory
#' result <- mergeRDSFiles(RDSFiles)
#'
#' ## Print the number and the position of the nucleosomes
#' result$k
#' result$mu
#'
#' ## Class of the output object
#' class(result)
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
mergeRDSFiles <- function(RDSFiles) {

    ## Validate parameters
    validateRDSFilesParameters(RDSFiles)

    ## Return merge information provided by each file
    return(mergeAllRDSFiles(RDSFiles))
}


#' @title A post-treatment function to merge closely positioned nucleosomes,
#' from the same chromosome, identified by the \code{\link{rjmcmc}} function.
#'
#' @description A helper function which merges closely positioned nucleosomes
#' to rectify the over splitting and provide a more conservative approach.
#' Beware that each chromosome must be treated separatly.
#'
#' @param reads a \code{GRanges} containing forward and
#' reverse reads. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
#' Default: \code{NULL}.
#'
#' @param resultRJMCMC an object of \code{class}
#' "rjmcmcNucleosomes" or "rjmcmcNucleosomesMerge", the information
#' about nucleosome positioning for an entire chromosome or a region that must
#' be treated as one unit.
#'
#' @param extendingSize a positive \code{numeric} or a positive \code{integer}
#' indicating the size of the consensus region used to group closeley
#' positioned nucleosomes.The minimum size of the consensus region is equal to
#' twice the value of the \code{extendingSize} parameter. The numeric will
#' be treated as an integer. Default: 74.
#'
#' @param chrLength a positive \code{numeric} or a positive \code{integer}
#' indicating the length of the current chromosome. The length of the
#' chromosome is used to ensure that the consensus positions are all
#' located inside the chromosome.
#'
#' @return a \code{GRanges}, the updated nucleosome positions.
#' When no nucleosome is present, \code{NULL} is returned.
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_02)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(reads = reads_demo_02,
#'             seqName = "chr_SYNTHETIC", nbrIterations = 1000,
#'             lambda = 2, kMax = 30, minInterval = 146,
#'             maxInterval = 490, minReads = 3, vSeed = 11)
#'
#' ## Before post-treatment
#' result
#'
#' ##Post-treatment function which merged closely positioned nucleosomes
#' postResult <- postTreatment(reads = reads_demo_02,
#'                     seqName = "chr_SYNTHETIC", result, 100, 73500)
#'
#' ## After post-treatment
#' postResult
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
postTreatment <- function(reads, seqName = NULL,
                            resultRJMCMC, extendingSize = 74L, chrLength) {

    ## Validate parameters
    validatePrepMergeParameters(reads, seqName,
                                    resultRJMCMC, extendingSize, chrLength)

    ## Only keep reads associated to the specified chromosome
    if (!is.null(seqName)) {
        reads <- reads[seqnames(reads) == seqName]
    }

    ## Run post merging function and return results
    return(postMerge(reads, resultRJMCMC, extendingSize, chrLength))
}


#' @title Generate a graph of nucleosome positions with read coverage
#'
#' @description Generate a graph for
#' a \code{GRanges} or a \code{GRangesList} of nucleosome positions. In
#' presence of only one prediction (with multiples nucleosome positions), a
#' \code{GRanges} is used. In presence of more thant one predictions (as
#' example, before and after post-treatment or results from
#' different software), a \code{GRangesList} with
#' one entry per prediction is used. All predictions must have been obtained
#' using the same reads.
#'
#' @param nucleosomePositions a \code{GRanges} or a \code{GRangesList}
#' containing the nucleosome positions for one or
#' multiples predictions obtained using the same reads. In presence of
#' only one prediction (with multiples nucleosome positions), a \code{GRanges}
#' is used. In presence of more thant one predictions (as example, before and
#' after post-treatment or results from different software), a
#' \code{GRangesList} with one entry per prediction is used.
#'
#' @param reads a \code{GRanges} containing forward and
#' reverse reads. The \code{GRanges} should contain at least one read.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
#' Default: \code{NULL}.
#'
#' @param xlab a \code{character} string containing the label of the x-axis.
#'
#' @param ylab a \code{character} string containing the label of the y-axis.
#'
#' @param names a \code{vector} of a \code{character} string containing the
#' label of each prediction set. The \code{vector} must be the same length of
#' the \code{nucleosomePositions} \code{list} or 1 in presence of a
#' \code{vector}. When \code{NULL}, the name of the elements of the \code{list}
#' are used or the string "Nucleosome" for a \code{vector} are used.
#' Default: \code{NULL}.
#'
#' @return a graph containing the nucleosome positions and the read coverage
#'
#' @examples
#'
#' ## Load reads dataset
#' data(reads_demo_01)
#'
#' ## Run RJMCMC method
#' result <- rjmcmc(reads = reads_demo_01,
#'             seqName = "chr_SYNTHETIC",
#'             nbrIterations = 4000, lambda = 2, kMax = 30,
#'             minInterval = 146, maxInterval = 292, minReads = 5,
#'             vSeed = 10213)
#'
#' ## Create graph using the synthetic map
#' plotNucleosomes(nucleosomePositions = result$mu, seqName = "chr_SYNTHETIC",
#'             reads = reads_demo_01)
#'
#' @author Astrid Deschenes
#' @importFrom IRanges coverage
#' @importFrom graphics plot lines abline points legend polygon
#' @importFrom grDevices rainbow
#' @importFrom BiocGenerics start end
#' @export
plotNucleosomes <- function(nucleosomePositions, reads,
                                seqName = NULL, xlab = "position",
                                ylab = "coverage", names=NULL) {

    validatePlotNucleosomesParameters(nucleosomePositions,
                        reads, seqName, xlab, ylab, names)

    ## Only keep reads associated to the specified chromosome
    if (!is.null(seqName)) {
        reads <- reads[seqnames(reads) == seqName]
        nucleosomePositions <- nucleosomePositions[
                    seqnames(nucleosomePositions) == seqName]
    } else {
        seqName <- as.character(seqnames(reads))[1]
    }

    ## Set variables differently if vector or list
    if (class(nucleosomePositions)=="GRangesList") {
        nbrItems <-length(nucleosomePositions)
        posColors <- c(rainbow(nbrItems), "gray")
        if (is.null(names) && !is.null(names(nucleosomePositions))) {
            extraNames <- names(nucleosomePositions)
        } else if(!is.null(names)){
            extraNames <- names
        } else {
            extraNames <- 1:nbrItems
        }
    } else {
        nbrItems <-1
        posColors <- c("green", "gray")
        if (is.null(names)) {
            extraNames <- "Nucleosome"
        } else {
            extraNames <- names
        }
    }

    posNames <- c(extraNames, "Coverage")

    coverageSeqName <- coverage(reads)[[seqName]]

    ## Set Y axis maximum range
    y_max <- max(coverageSeqName, na.rm = TRUE) + 10

    ## Step in between each result, when more than one result
    step = ceiling(y_max / 80)

    ## Always set Y axis minimum to zero
    y_min <- -1 - (step * nbrItems)

    ## Set X axis minimum ans maximum
    x_min <- min(c(unlist(start(nucleosomePositions)),
                            start(reads), end(reads)))
    x_min <- floor(x_min)
    x_max <- max(c(unlist(start(nucleosomePositions)),
                            start(reads), end(reads)))
    x_max <- ceiling(x_max)

    # Plot coverage
    coverage <- c(0, as.integer(coverageSeqName), 0)
    position <- c(0, 1:(length(coverage) - 1))
    plot(coverageSeqName, type = "l", col = "gray",
            ylim = c(y_min, y_max), xlim = c(x_min, x_max), xlab = xlab,
            ylab = ylab)
    polygon(c(x_min, position, 0), c(0, coverage, 0), col="gray",
            border = "gray", ylim = c(y_min, y_max), xlim = c(x_min, x_max))

    # Plot nucleosome positions
    if (nbrItems > 1) {
        for (i in 1:nbrItems) {
            y_pos = (-(step)) * i
            nucl <- start(nucleosomePositions[[i]])
            points(nucl, rep(y_pos, length(nucl)), ylim = c(y_min, y_max),
                    xlim = c(x_min, x_max), col = posColors[i], pch = 19)
        }
    } else {
        points(start(nucleosomePositions), rep(-(step),
                length(nucleosomePositions)),
                ylim = c(y_min, y_max), xlim = c(x_min, x_max),
                col = posColors[1], pch = 19)
    }

    # Add legend
    legend("top", posNames, fill = posColors, bty = "n", horiz = TRUE)
}

#' @title Split a \code{GRanges} containing reads in a list of smaller
#' segments for the \code{rjmcmc} function.
#'
#' @description Split a \code{GRanges} of reads (as example, the reads from
#' a chromosome) in a \code{list} of smaller \code{GRanges} sot that the
#' \code{rjmcmc} function can be run on each segments.
#'
#' @param reads a \code{GRanges}, the reads that need to be segmented.
#'
#' @param zeta a positive \code{integer} or \code{numeric}, the length
#' of the nucleosomes. Default: 147.
#'
#' @param delta a positive \code{integer} or \code{numeric}, the accepted
#' range of overlapping section between segments. The overlapping section
#' being \code{zeta} + \code{delta}.
#'
#' @param maxLength a positive \code{integer} or \code{numeric}, the
#' length of each segment.
#'
#' @return a \code{GRangesList} containing all the segments.
#'
#' @examples
#'
#' ## Load synthetic dataset of reads
#' data(syntheticNucleosomeReads)
#'
#' ## Use dataset of reads to create GRanges object
#' sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
#'     ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
#'     end = syntheticNucleosomeReads$dataIP$end),
#'     strand = syntheticNucleosomeReads$dataIP$strand)
#'
#' # Segmentation of the reads
#' segmentation(reads = sampleGRanges, zeta = 147, delta = 50,
#' maxLength = 1000)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges splitAsList
#' @export
segmentation <- function(reads, zeta = 147, delta, maxLength) {

    validateSegmentationParameters(reads, zeta, delta, maxLength)

    # Set min and max positions
    posMin <- min(start(reads))
    posMax <- max(end(reads))

    # Create segments
    starts = seq(posMin, posMax, by = (maxLength - (zeta + delta)))
    subject = GRanges(seqlevels(reads), IRanges(starts, width=maxLength))
    hits = findOverlaps(reads, subject)
    splitAsList(reads[queryHits(hits)], subjectHits(hits))
}

#' @title Nucleosome positioning mapping on a large segment, up to a chromosome
#'
#' @description Use of a fully Bayesian hierarchical model for chromosome-wide
#' profiling of nucleosome positions based on high-throughput short-read
#' data (MNase-Seq data). Beware that for a genome-wide profiling, each
#' chromosome must be treated separatly. This function is optimized to run
#' on an entire chromosome.
#'
#' The function will process by splittingg the \code{GRanges} of reads
#' (as example, the reads from a chromosome) in a \code{list} of smaller
#' \code{GRanges} segments that can be run by the
#' \code{rjmcmc} function. All those steps are done automatically.
#'
#' @param reads a \code{GRanges}, the forward and reverse
#' reads that need to be segmented.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
#' Default: \code{NULL}.
#'
#' @param zeta a positive \code{integer} or \code{numeric}, the length
#' of the nucleosomes. Default: 147.
#'
#' @param delta a positive \code{integer} or \code{numeric}, the accepted
#' range of overlapping section between segments. The overlapping section
#' being \code{zeta} + \code{delta}.
#'
#' @param maxLength a positive \code{integer} or \code{numeric}, the
#' length of each segment.
#'
#' @param nbrIterations a positive \code{integer} or \code{numeric}, the
#' number of iterations. Non-integer values of
#' \code{nbrIterations} will be casted to \code{integer} and truncated towards
#' zero.
#'
#' @param kMax a positive \code{integer} or \code{numeric}, the maximum number
#' of degrees of freedom per region. Non-integer values
#' of \code{kMax} will be casted to \code{integer} and truncated towards zero.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution. Default: 3.
#'
#' @param minInterval a \code{numeric}, the minimum distance between two
#' nucleosomes.
#'
#' @param maxInterval a \code{numeric}, the maximum distance between two
#' nucleosomes.
#'
#' @param minReads a positive \code{integer} or \code{numeric}, the minimum
#' number of reads in a potential canditate region. Non-integer values
#' of \code{minReads} will be casted to \code{integer} and truncated towards
#' zero. Default: 5.
#'
#' @param adaptIterationsToReads a \code{logical} indicating if the number
#' of iterations must be modified in function of the number of reads.
#' Default: \code{TRUE}.
#'
#' @param vSeed a \code{integer}. A seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: -1.
#'
#' @param nbCores a positive \code{integer}, the number
#' of cores used to run in parallel. Default: 1.
#'
#' @param dirOut a \code{character} string. The name of the directory
#' where 2 directories are created (if they don't already exists).
#' The directory "dirOut/results" contents the rjmcmc results for each segment.
#' The directory "dirOut/done" contents file a log file for each segment in
#' RData format. If the log file for a segment is in the directory,
#' the program considers that it is has been processed and run the next
#' segment. Default: "out".
#'
#' @param saveSEG a \code{logical}. When \code{TRUE}, a RDS file containing
#' the segments generated by  \code{\link{segmentation}} function is
#' saved in directory named from paramter \code{dirOut}.
#' Default: \code{FALSE}.
#'
#' @param saveAsRDS a \code{logical}. When \code{TRUE}, a RDS file containing
#' the complete output of the \code{rjmcmc} function is created.
#' Default: \code{FALSE}.
#'
#' @return a \code{list} of class
#' "rjmcmcNucleosomesBeforeAndAfterPostTreatment" containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{GRanges} containing the positions of the nucleosomes.
#'     \item kPost a \code{integer}, the number of nucleosomes after
#' post-treatment and '*' as strand. The \code{seqnames} of the \code{GRanges}
#' correspond to the \code{seqName} input value. \code{NA} when no nucleosome
#' is detected.
#'     \item muPost a \code{GRanges} containing the positions of the
#' nucleosomes after post-treament and '*' as strand. The \code{seqnames}
#' of the \code{GRanges} correspond to the \code{seqName} input value.
#' \code{NA} when no nucleosome is detected.
#' }
#'
#' @examples
#'
#' ## Load synthetic dataset of reads
#' data(syntheticNucleosomeReads)
#'
#' ## Use dataset of reads to create GRanges object
#' sampleGRanges <- GRanges(syntheticNucleosomeReads$dataIP)
#'
#' ## Run nucleosome detection on the entire sample
#' \dontrun{result <- rjmcmcCHR(reads = sampleGRanges, zeta = 147, delta=50,
#' maxLength=1200, nbrIterations = 1000, lambda = 3, kMax = 30,
#' minInterval = 146, maxInterval = 292, minReads = 5, vSeed = 10113,
#' nbCores = 2, saveAsRDS = FALSE)}
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom GenomicRanges strand
#' @export
rjmcmcCHR <- function(reads, seqName = NULL, zeta = 147,
                        delta, maxLength,
                        nbrIterations, kMax, lambda = 3,
                        minInterval, maxInterval, minReads = 5,
                        adaptIterationsToReads = TRUE, vSeed = -1,
                        nbCores = 1, dirOut="out",
                        saveAsRDS = FALSE, saveSEG = TRUE){

    if(!dir.exists(dirOut)){
        dir.create(dirOut)
    }

    dirDone <- paste0(dirOut, "/done")
    dirResults <- paste0(dirOut, "/results")

    if(!dir.exists(dirResults)){
        dir.create(dirResults)
    }

    if(!dir.exists(dirDone)){
        dir.create(dirDone)
    }

    ## Only keep reads associated to the specified chromosome
    if (!is.null(seqName)) {
        reads <- reads[seqnames(reads) == seqName]
    }

    ## Split reads into segments
    seg <- segmentation(reads, zeta, delta, maxLength)

    if(saveSEG){
        options(digits.secs = 2)
        file_name <- gsub(Sys.time(), pattern = "[:. ]", replacement = "_",
                            perl = TRUE)
        saveRDS(object = seg,
                file = paste0(dirOut,"/seg_", file_name, ".RDS"))
    }

    nbSeg <- length(seg)

    param <- SnowParam(workers = nbCores, stop.on.error = TRUE)

    ## Run RJMCMC on each segment
    a <- bplapply(1:nbSeg, FUN = runCHR, seg, niter = nbrIterations,
                    kmax = kMax, lambda = lambda, ecartmin = minInterval,
                    ecartmax = maxInterval, minReads = minReads,
                    adaptNbrIterations = adaptIterationsToReads,
                    dirOut = dirOut,
                    vSeed = vSeed, saveAsRDS = saveAsRDS, BPPARAM = param)

    ## Merge results from all segments
    results <- mergeAllRDSFilesFromDirectory(dirResults)

    ## Post-treatment of the nucleosomes
    resultPostTreatement <- postTreatment(reads = reads,
                                seqName = seqName, results,
                                chrLength=max(start(reads), end(reads)) + 1000)

    results$muPost <- resultPostTreatement

    results$kPost <- length(resultPostTreatement)

    class(results)<-"rjmcmcNucleosomesBeforeAndAfterPostTreatment"

    return(results)
}
