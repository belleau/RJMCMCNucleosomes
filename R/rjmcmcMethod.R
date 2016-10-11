#' @title Nucleosome positioning mapping
#'
#' @description Use of a fully Bayesian hierarchical model for chromosome-wide
#' profiling of nucleosome positions based on high-throughput short-read
#' data (MNase-Seq data). Beware that for a genome-wide profiling, each
#' chromosome must be treated separatly.
#'
#' @param startPosForwardReads a \code{vector} of \code{numeric}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of \code{numeric}, the
#' start position of all the reverse reads. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param nbrIterations a positive \code{integer} or \code{numeric}, the
#' number of iterations. Non-integer values of
#' \code{nbrIterations} will be casted to \code{integer} and truncated towards
#' zero.
#'
#' @param kMax a positive \code{integer} or \code{numeric}, the maximum number
#' of nucleosomes per region. Non-integer values
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
#' data(reads_demo)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 1000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5, vSeed = 10113)
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
#' @author Rawane Samb, Pascal Belleau, Astrid Deschênes
#' @importFrom stats aggregate
#' @export
rjmcmc <- function(startPosForwardReads, startPosReverseReads,
                    nbrIterations, kMax, lambda = 3,
                    minInterval, maxInterval, minReads = 5,
                    adaptIterationsToReads = TRUE, vSeed = -1) {

    # Get call information
    cl <- match.call()

    # Parameters validation
    validateRJMCMCParameters(startPosForwardReads = startPosForwardReads,
                            startPosReverseReads = startPosReverseReads,
                            nbrIterations = nbrIterations,
                            kMax = kMax,
                            lambda = lambda,
                            minInterval = minInterval,
                            maxInterval = maxInterval,
                            minReads = minReads,
                            adaptIterationsToReads = adaptIterationsToReads,
                            vSeed = vSeed)

    # Find nucleosome positions
    resultRJMCMC <- rjmcmcNucleo(startPosForwardReads, startPosReverseReads,
                                 nbrIterations, kMax, lambda,
                                 minInterval, maxInterval, minReads,
                                 adaptIterationsToReads, vSeed)

    if (is.null(resultRJMCMC)) {
        ## Set values when no nucleosome can be found
        k = 0
        mu = NA
        k_max = NA
    } else {
        ## Set values when no nucleosome can be found
        # Find k value with the maximum of iterations
        iterPerK <- data.frame(k=resultRJMCMC$k, it=resultRJMCMC$it)
        sumIterPerK <- aggregate(it ~ k, data=iterPerK, sum)
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
#' TODO : A faire pour nouveau format
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
#'     \item sigmaf a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta a \code{vector} of \code{numeric} of length
#' \code{k}, the distance between the maxima of the forward and
#' reverse reads position densities for each nucleosome.
#' }
#'
#' @examples
#'
#' ## Use a directory present in the RJMCMC package
#' directoryWithRDSFiles <- system.file("extdata", package = "RJMCMCNucleosomes")
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
#' @author Pascal Belleau, Astrid Deschênes
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
#'     \item sigmaf a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta a \code{vector} of \code{numeric} of length
#' \code{k}, the distance between the maxima of the forward and
#' reverse reads position densities for each nucleosome.
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
#' @author Pascal Belleau, Astrid Deschênes
#' @export
mergeRDSFiles <- function(RDSFiles) {

    ## Validate parameters
    validateRDSFilesParameters(RDSFiles)

    ## Return merge information provided by each file
    return(mergeAllRDSFiles(RDSFiles))
}


#' @title A post treatment function to merge closely positioned nucleosomes,
#' from the same chromosome, identified by the \code{\link{rjmcmc}} function.
#'
#' @description A helper function which merges closely positioned nucleosomes
#' to rectify the over splitting and provide a more conservative approach.
#' Beware that each chromosome must be treated separatly.
#'
#' @param startPosForwardReads a \code{vector} of \code{numeric}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of \code{numeric}, the
#' start position of all the reverse reads. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param resultRJMCMC an object of \code{class}
#' "rjmcmcNucleosomes" or "rjmcmcNucleosomesMerge", the information
#'  about nucleosome positioning for an entire chromosome or a region that must
#'  be treated as one unit.
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
#' nucleosome positions.
#'
#' @examples
#'
#' ## TODO : A faire
#'
#' ## Fix seed
#' set.seed(1132)
#'
#' ## Loading dataset
#' data(reads_demo)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 1000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5)
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' ##postResult <- postTreatment(startPosForwardReads = reads_demo$readsForward,
#' ##         startPosReverseReads = reads_demo$readsReverse, result, 74, 73500)
#'
#' ##postResult
#'
#' @author Pascal Belleau, Astrid Deschênes
#' @export
postTreatment <- function(startPosForwardReads, startPosReverseReads,
                           resultRJMCMC, extendingSize = 74L, chrLength) {

    ## Validate parameters
    validatePrepMergeParameters(startPosForwardReads, startPosReverseReads,
                                        resultRJMCMC, extendingSize, chrLength)

    ## Run post merging function and return results
    return(postMerge(startPosForwardReads, startPosReverseReads,
              resultRJMCMC, extendingSize, chrLength))
}


#' @title Generate a graph of nucleosome position
#'
#' @description Generate a graph for
#' a list of nucleosomes
#'
#' @param nucleosomesPosition a list of \code{numeric}, the nucleosome positions
#'
#' @param reads an \code{IRanges} containing all the reads
#'
#' @param xlab a \code{character} string containing the label of the x-axis
#'
#' @param ylab a \code{character} string containing the label of the y-axis
#'
#' @return a graph containing the nucleosome position and the read coverage
#'
#' @examples
#'
#' ## Generate a synthetic sample with 10 well-positioned nucleosomes, 2 fuzzy
#' ## nucleosomes and 2 deleted nucleosomes using a Student distribution
#' ## with a variance of 10 for the well-positioned nucleosomes,
#' ## a variance of 20 for the fuzzy nucleosomes
#' library(nucleoSim)
#' nucleosomeSample <- syntheticNucReadsFromDist(wp.num=10, wp.del=2,
#' wp.var=10, fuz.num=2, fuz.var=20, max.cover=100, dist="Student",
#' nuc.len=147, len.var=12, read.len=45, lin.len=20, rnd.seed=155, offset=100)
#'
#' dataIP <-nucleosomeSample$dataIP
#'
#' forwardReads <- dataIP[dataIP$strand == "+",]$start
#' reverseReads <- dataIP[dataIP$strand == "-",]$end
#'
#' result <- rjmcmc(startPosForwardReads = forwardReads,
#'          startPosReverseReads = reverseReads,
#'          nbrIterations = 4000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5, vSeed = 10213)
#'
#' reads <-IRanges(start = dataIP$start, end=dataIP$end)
#'
#' ## Create graph using the synthetic map
#' plotNucleosomes(nucleosomesPosition = result$mu, reads = reads)
#'
#' @author Astrid Deschenes
#' @importFrom IRanges coverage
#' @importFrom graphics plot lines abline points legend polygon
#' @importFrom stats start end
#' @export
plotNucleosomes <- function(nucleosomesPosition, reads, xlab="position",
                                ylab="coverage") {


    ## Set Y axis maximum range
    y_max <- max(coverage(reads), na.rm = TRUE) + 10

    ## Always set Y axis minimum to zero
    y_min <- 0

    ## Set X axis minimum ans maximum
    x_min <- min(nucleosomesPosition, start(reads), end(reads))
    x_min <- floor(x_min)
    x_max <- max(nucleosomesPosition, start(reads), end(reads))
    x_max <- ceiling(x_max)

    # Plot coverage
    coverage <- c(0, as.integer(coverage(reads)), 0)
    position <- c(x_min, 1:(length(coverage) - 1))
    plot(position, coverage, type = "l", col = "gray",
         ylim = c(y_min, y_max), xlim = c(x_min, x_max), xlab=xlab, ylab=ylab)
    polygon(c(x_min, position, 0), c(0, coverage, 0), col="gray", border = "gray")

    # Plot nucleosome positions
    points(nucleosomesPosition, rep(0, length(nucleosomesPosition)),
            col = "forestgreen",  pch = 19)

    # Add legend
    legend("top", c("Nucleosome", "Coverage"),
               fill = c("forestgreen", "gray"), bty = "n",
               horiz = TRUE)
}

