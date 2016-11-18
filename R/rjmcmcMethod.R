#' @title Nucleosome positioning mapping on a segment
#'
#' @description Use of a fully Bayesian hierarchical model for chromosome-wide
#' profiling of nucleosome positions based on high-throughput short-read
#' data (MNase-Seq data). Beware that for a genome-wide profiling, each
#' chromosome must be treated separatly. This function is optimized to run
#' on segments that are smaller sections of the chromosome.
#'
#' @param forwardandReverseReads a \code{GRanges} containing forward and
#' reverse reads.
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
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes. \code{NA} when no nucleosome is
#' detected.
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
#' result <- rjmcmc(forwardandReverseReads = reads_demo_01,
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
rjmcmc <- function(forwardandReverseReads,
                    nbrIterations, kMax, lambda = 3,
                    minInterval, maxInterval, minReads = 5,
                    adaptIterationsToReads = TRUE, vSeed = -1,
                    saveAsRDS = FALSE) {

    # Get call information
    cl <- match.call()

    # Parameters validation
    validateRJMCMCParameters(forwardandReverseReads = forwardandReverseReads,
                            nbrIterations = nbrIterations,
                            kMax = kMax,
                            lambda = lambda,
                            minInterval = minInterval,
                            maxInterval = maxInterval,
                            minReads = minReads,
                            adaptIterationsToReads = adaptIterationsToReads,
                            vSeed = vSeed)

    resultRJMCMC <- NULL

    if (length(forwardandReverseReads) > 0) {

        startPosForwardReads <- start(forwardandReverseReads[
                                        strand(forwardandReverseReads) == "+"])

        startPosReverseReads <- end(forwardandReverseReads[
                            strand(forwardandReverseReads) == "-"])

        # Find nucleosome positions
        if(length(startPosForwardReads) > 0 & length(startPosReverseReads) > 0){
            resultRJMCMC <- rjmcmcNucleo(startPosForwardReads,
                                        startPosReverseReads,
                                        nbrIterations, kMax, lambda,
                                        minInterval, maxInterval, minReads,
                                        adaptIterationsToReads, vSeed)
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
        ## Set values when no nucleosome can be found
        # Find k value with the maximum of iterations
        iterPerK <- data.frame(k = resultRJMCMC$k, it = resultRJMCMC$it)
        sumIterPerK <- aggregate(it ~ k, data = iterPerK, sum)
        maxRow <- which.max( sumIterPerK[,"it"])
        k <- sumIterPerK$k[maxRow]
        # Find mu values associated to the k value
        mu <- resultRJMCMC$muHat[k,][1:k]
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
#'     \item mu a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
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
#'     \item mu a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
#' }
#'
#' @examples
#'
#' ## Use RDS files present in the RJMCMC package
#' RDSFiles <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#' full.names = TRUE, pattern = "*rds")
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
#' @param forwardandReverseReads a \code{GRanges} containing forward and
#' reverse reads. The \code{GRanges} should contain at least one read.
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
#' indicating the lenght of the current chromosome. The length of the
#' chromosome is used to ensure that the consensus positions are all
#' located inside the chromosome.
#'
#' @return a \code{array} of \code{numeric}, the updated values of the
#' nucleosome positions. When no nucleosome is present, \code{NULL} is
#' returned.
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_02)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(forwardandReverseReads = reads_demo_02,
#'             nbrIterations = 1000, lambda = 2, kMax = 30,
#'             minInterval = 146, maxInterval = 490, minReads = 3, vSeed = 11)
#'
#' ## Before post-treatment
#' result
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' postResult <- postTreatment(forwardandReverseReads = reads_demo_02,
#'             result, 100, 73500)
#'
#' ## After post-treatment
#' postResult
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
postTreatment <- function(forwardandReverseReads,
                            resultRJMCMC, extendingSize = 74L, chrLength) {

    ## Validate parameters
    validatePrepMergeParameters(forwardandReverseReads,
                                        resultRJMCMC, extendingSize, chrLength)

    ## Run post merging function and return results
    return(postMerge(forwardandReverseReads,
                resultRJMCMC, extendingSize, chrLength))
}


#' @title Generate a graph of nucleosome positions with read coverage
#'
#' @description Generate a graph for
#' a \code{list} or a \code{vector} of nucleosome positions. In presence of
#' only one prediction (with multiples nucleosome positions), a \code{vector}
#' is used. In presence of more thant one predictions (as example, before and
#' after post-treatment or results from different software), a \code{list} with
#' one entry per prediction is used. All predictions must have been obtained
#' using the same reads.
#'
#' @param nucleosomePositions a \code{list} or a \code{vector} of
#' \code{numeric}, the nucleosome positions for one or
#' multiples predictions are obtained using the same reads. In presence of
#' only one prediction (with multiples nucleosome positions), a \code{vector}
#' is used. In presence of more thant one predictions (as example, before and
#' after post-treatment or results from different software), a \code{list} with
#' one entry per prediction is used.
#'
#' @param reads an \code{IRanges} containing all the reads.
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
#' data(reads_demo_01)
#'
#' result <- rjmcmc(forwardandReverseReads = reads_demo_01,
#'             nbrIterations = 4000, lambda = 2, kMax = 30,
#'             minInterval = 146, maxInterval = 292, minReads = 5,
#'             vSeed = 10213)
#'
#' reads <-IRanges(start = start(reads_demo_01),
#'             end = end(reads_demo_01))
#'
#' ## Create graph using the synthetic map
#' plotNucleosomes(nucleosomePositions = result$mu, reads = reads)
#'
#' @author Astrid Deschenes
#' @importFrom IRanges coverage
#' @importFrom graphics plot lines abline points legend polygon
#' @importFrom grDevices rainbow
#' @importFrom BiocGenerics start end
#' @export
plotNucleosomes <- function(nucleosomePositions, reads, xlab = "position",
                                ylab = "coverage", names=NULL) {

    validatePlotNucleosomesParameters(nucleosomePositions, reads, xlab,
                                        ylab, names)

    ## Set variables differently if vector or list
    if (!is.atomic(nucleosomePositions)) {
        nbrItems <-length(nucleosomePositions)
        posColors <- c(rainbow(nbrItems), "gray")
        if (is.null(names)) {
            extraNames <- names(nucleosomePositions)
        } else {
            extraNames <- names
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

    ## Set Y axis maximum range
    y_max <- max(coverage(reads), na.rm = TRUE) + 10

    ## Step in between each result, when more than one result
    step = ceiling(y_max / 80)

    ## Always set Y axis minimum to zero
    y_min <- -1 - (step*nbrItems)

    ## Set X axis minimum ans maximum
    x_min <- min(c(unlist(nucleosomePositions), start(reads), end(reads)))
    x_min <- floor(x_min)
    x_max <- max(c(unlist(nucleosomePositions), start(reads), end(reads)))
    x_max <- ceiling(x_max)

    # Plot coverage
    coverage <- c(0, as.integer(coverage(reads)), 0)
    position <- c(0, 1:(length(coverage) - 1))
    plot(coverage(reads), type = "l", col = "gray",
            ylim = c(y_min, y_max), xlim = c(x_min, x_max), xlab = xlab,
            ylab = ylab)
    polygon(c(x_min, position, 0), c(0, coverage, 0), col="gray",
            border = "gray", ylim = c(y_min, y_max), xlim = c(x_min, x_max))

    # Plot nucleosome positions
    if (nbrItems > 1) {
        for (i in 1:nbrItems) {
            y_pos = (-(step)) * i
            nucl <- nucleosomePositions[[i]]
            points(nucl, rep(y_pos, length(nucl)), ylim = c(y_min, y_max),
                    xlim = c(x_min, x_max), col = posColors[i], pch = 19)
        }
    } else {
        points(nucleosomePositions, rep(-(step), length(nucleosomePositions)),
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
#' @param dataIP a \code{GRanges}, the reads that need to be segmented.
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
#' @return a \code{list} of \code{GRanges}, the list of segments.
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
#' segmentation(sampleGRanges, zeta = 147, delta = 50, maxLength = 1000)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @export
segmentation <- function(dataIP, zeta = 147, delta, maxLength) {

    validateSegmentationParameters(dataIP, zeta, delta, maxLength)

    # Set min and max position
    posMin <- min(start(dataIP))
    posMax <- max(end(dataIP))

    # Segment GRanges
    lapply(seq(posMin, posMax, by = (maxLength - (zeta + delta))),
            function(x, dataIP, zeta, delta, maxLength){
                dataIP[start(dataIP) >= x & start(dataIP) <= (x + maxLength)]
            },
            dataIP=dataIP, zeta=zeta, delta=delta, maxL = maxLength)
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
#' @param forwardandReverseReads a \code{GRanges}, the forward and reverse
#' reads that need to be segmented.
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
#'     \item mu a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
#'     \item muPost a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes after post-treament.
#' }
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
#' ## Run nucleosome detection on the entire sample
#' \dontrun{result <- rjmcmcCHR(forwardandReverseReads=sampleGRanges,
#'              zeta = 147, delta=50, maxLength=1200,
#'              nbrIterations = 1000, lambda = 3, kMax = 30,
#'              minInterval = 146, maxInterval = 292, minReads = 5,
#'              vSeed = 10113, nbCores = 2, saveAsRDS = FALSE)}
#'
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom GenomicRanges strand
#' @export
rjmcmcCHR <- function(forwardandReverseReads, zeta = 147, delta, maxLength,
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

    seg <- segmentation(forwardandReverseReads, zeta, delta, maxLength)

    if(saveSEG){
        options(digits.secs = 2)
        file_name <- gsub(Sys.time(), pattern = "[:. ]", replacement = "_",
                            perl = TRUE)
        saveRDS(object = seg,
                file = paste0(dirOut,"/seg_", file_name, ".RDS"))
    }

    nbSeg <- length(seg)

    param <- SnowParam(workers = nbCores, stop.on.error = TRUE)

    a <- bplapply(1:nbSeg, FUN = runCHR, seg, niter = nbrIterations,
                    kmax = kMax, lambda = lambda, ecartmin = minInterval,
                    ecartmax = maxInterval, minReads = minReads,
                    adaptNbrIterations = adaptIterationsToReads,
                    dirOut = dirOut,
                    vSeed = vSeed, saveAsRDS = saveAsRDS, BPPARAM = param)

    results <- mergeAllRDSFilesFromDirectory(dirResults)

    resultPostTreatement <- postTreatment(forwardandReverseReads =
                                                forwardandReverseReads,
                                results,
                                chrLength=max(start(forwardandReverseReads),
                                        end(forwardandReverseReads)) + 1000)

    results$muPost <- resultPostTreatement

    results$kPost <- length(resultPostTreatement)

    class(results)<-"rjmcmcNucleosomesBeforeAndAfterPostTreatment"

    return(results)
}
