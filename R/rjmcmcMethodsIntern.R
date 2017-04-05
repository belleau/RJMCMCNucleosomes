#' @title Interface for the RJMCMC nucleosome mapping method in C++
#'
#' @description Function that calls the core of the nucleosome positioning
#' mapping function that is implemented in C++.
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
#' @return a \code{list} containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item k_max a \code{integer}, the maximum number of nucleosomes
#' obtained during the iteration process.
#'     \item it a \code{vector} of \code{integer} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item nbState a \code{integer}, the number of changes of state.
#'     \item mu a \code{matrix} of \code{numeric} with \code{k_max} columns
#' and \code{nbState} row containing, in each row, the \code{mu} values
#' associated the the state identified by the row number.
#'     \item muHat a \code{matrix} of \code{numeric} with \code{k_max} columns
#' and \code{k_max} rows containing, in each row, the mean \code{mu} values
#' associated the number of nucleosomes detected. The row number
#' corresponds to the number of nucleosomes detected.
#'     \item nbK a \code{vector} of length \code{k_max} containing
#' \code{integer}, the number of iterations
#' which detected a specific number of nucleosomes. The position in the vector
#' correspond to the number of nucleosomes.
#' }
#'
#' @examples
#'
## Loading dataset
#' data(reads_demo_01)
#'
#' forward <- start(reads_demo_01[strand(reads_demo_01) == "+"])
#' reverse <- end(reads_demo_01[strand(reads_demo_01) == "-"])
#'
#' ## Run nucleosome positioning
#' result <- RJMCMCNucleosomes:::rjmcmcNucleo(
#'             startPosForwardReads = forward,
#'             startPosReverseReads = reverse,
#'             nbrIterations = 1000, lambda = 2, kMax = 30,
#'             minInterval = 146, maxInterval = 292, minReads = 5,
#'             adaptIterationsToReads = TRUE, vSeed = -1)
#'
#' ## Print the final estimation of the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
#'
rjmcmcNucleo <- function(startPosForwardReads,
                                startPosReverseReads,
                                nbrIterations, kMax, lambda,
                                minInterval, maxInterval,
                                minReads = 5L,
                                adaptIterationsToReads = TRUE, vSeed = -1) {
    # Call C++ function
    .Call('C_RJMCMCNucleosomes_rjmcmcNucleo',
            PACKAGE = 'RJMCMCNucleosomes',
            startPosForwardReads, startPosReverseReads,
            nbrIterations, kMax, lambda, minInterval,
            maxInterval, minReads, adaptIterationsToReads, vSeed)
}


#' @title Parameters validation for the \code{\link{rjmcmc}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{rjmcmc}} function.
#'
#' @param reads a \code{GRanges} containing all forward
#' and reverse reads.The start positions of both reads are going to be used
#' for the analysis. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
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
#' of the Poisson distribution.
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
#' zero.
#'
#' @param adaptIterationsToReads a \code{logical} indicating if the number
#' of iterations must be modified in function of the number of reads.
#'
#' @param vSeed a \code{integer}. A seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
#'     ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
#'     strand = Rle(strand(c("-", "+", "-", "+", "-")), c(1, 2, 2, 3, 2)))
#'
#' ## The function returns 0 when all paramaters are valid
#' RJMCMCNucleosomes:::validateRJMCMCParameters(reads = reads,
#' seqName = "chr1", nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1,
#' minInterval = 100, maxInterval = 200, adaptIterationsToReads = TRUE,
#' vSeed = 100)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{RJMCMCNucleosomes:::validateRJMCMCParameters(
#' reads = NA, seqName = "chr1",
#' nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1, minInterval = 100,
#' maxInterval = 200, adaptIterationsToReads = TRUE, vSeed = -1)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber runValue
#' @keywords internal
validateRJMCMCParameters <- function(reads,seqName, nbrIterations,
                                        kMax, lambda,
                                        minInterval, maxInterval, minReads,
                                        adaptIterationsToReads, vSeed) {

    ## Validate that the reads is a GRanges
    if (!(class(reads) == "GRanges" )) {
        stop(paste0("reads must be a GRanges"))
    }

    if (is.null(seqName) &&
        (length(runValue(seqnames(reads))) > 1)) {
        stop(paste0("seqName must be the name of one of the chromosomes ",
                "present in the GRanges"))
    }

    if (!is.null(seqName) && !is.character(seqName)) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## When seqName present, it needs to be the name of one of the
    ## chromosomes present in the GRanges
    if (!is.null(seqName) && is.character(seqName) && !(seqName %in%
                            runValue(seqnames(reads)))) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## Validate the nbrIterations parameter
    if (!(isSingleInteger(nbrIterations) || isSingleNumber(nbrIterations)) ||
            as.integer(nbrIterations) < 1) {
        stop("nbrIterations must be a positive integer or numeric")
    }

    ## Validate the kMax parameter
    if (!(isSingleInteger(kMax) || isSingleNumber(kMax)) ||
            as.integer(kMax) < 1) {
        stop("kMax must be a positive integer or numeric")
    }

    ## Validate the minReads parameter
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
            as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate the lambda parameter
    if (!isSingleNumber(lambda) || lambda <= 0) {
        stop("lambda must be a positive numeric")
    }

    ## Validate that adaptIterationsToReads is a logical
    if (!is.logical(adaptIterationsToReads)) {
        stop("adaptIterationsToReads must be a logical.")
    }

    ## Validate that vSeed is a numeric value
    if (!isSingleNumber(vSeed)) {
        stop("vSeed must be a numeric value.")
    }

    return(0)
}


#' @title Merge nucleosome information
#'
#' @description Merge nucleosome information present in multiple RDS files.
#'
#' @param arrayOfFiles a \code{array}, the name of each file that must be
#' used to merge nucleosome information.
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
#' ## Loading two files containing nucleosomes informations for two sections of
#' ## the same chromosome
#' file_1 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#'                 pattern = "RJMCMC_seg_01.RDS",
#'                 full.names = TRUE)
#'
#' file_2 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#'                 pattern = "RJMCMC_seg_02.RDS",
#'                 full.names = TRUE)
#'
#' ## Merging nucleosomes informations from the two files
#' result <- RJMCMCNucleosomes:::mergeAllRDSFiles(c(file_1, file_2))
#'
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom GenomicRanges GRangesList
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal
#'
mergeAllRDSFiles <- function(arrayOfFiles) {

    ## Create list that will contain data from all files
    result        <- list()
    result$k      <- 0
    mu <- GRangesList()

    ## Extract information from each file
    for (fileName in arrayOfFiles) {
        data <- readRDS(file = fileName)
        ## Only use data from rjmcmcNucleosomes or rjmcmcNucleosomesMerge class
        if ((is(data, "rjmcmcNucleosomesMerge") ||
                        is(data, "rjmcmcNucleosomes")) &
                        class(data$mu) == "GRanges") {
            result$k      <- result$k + data$k
            if(class(data$mu)=="GRanges"){
                mu[[fileName]] <- data$mu
            }
        }
    }

    ## Ensure that all values are ordered in ascending order of mu
    result$mu     <- sort(unlist(mu, use.names=FALSE))

    ## Assign class type to list
    class(result)<-"rjmcmcNucleosomesMerge"

    return(result)
}


#' @title Parameters validation for the \code{\link{mergeRDSFiles}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{mergeRDSFiles}} function.
#'
#' @param RDSFiles a \code{array}, the names of all RDS used to merge
#' nucleosome information. The files must contain R object of \code{class}
#' "rjmcmcNucleosomes" or "rjmcmcNucleosomesMerge".
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Loading a file
#' file_test <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#' pattern = "RJMCMC_seg_02.RDS", full.names = TRUE)
#'
#' ## Testing using a real file
#' RJMCMCNucleosomes:::validateRDSFilesParameters(file_test)
#'
#' @author Astrid Deschenes
#' @keywords internal
#'
validateRDSFilesParameters <- function(RDSFiles) {

    ## Validate the RDSFiles parameters
    if (!is.character(RDSFiles) || (length(RDSFiles) == 0)) {
        stop("RDSFiles must be a list of valid RDS files")
    }

    ## Validate that all files exist
    for (fileName in RDSFiles) {
        if (!file.exists(fileName)) {
            stop("The file \'", fileName, "\' does not exist.")
        }
    }

    return(0)
}


#' @title Parameters validation for the
#' \code{\link{mergeAllRDSFilesFromDirectory}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{mergeAllRDSFilesFromDirectory}} function.
#'
#' @param directory a \code{character}, the
#' name of the directory (relative or absolute path) containing RDS files.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load an existing directory
#' directory <- system.file("extdata", package = "RJMCMCNucleosomes")
#'
#' ## Testing using a real directory
#' RJMCMCNucleosomes:::validateDirectoryParameters(directory)
#'
#' @author Astrid Deschenes
#' @keywords internal
#'
validateDirectoryParameters <- function(directory) {

    ## Validate that the directory exist
    if (!file.exists(directory)) {
        stop("The directory \'", directory, "\' does not exist.")
    }

    return(0)
}


#' @title Parameters validation for the \code{\link{postMerge}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{postMerge}} function.
#'
#' @param reads a \code{GRanges} containing all forward
#' and reverse reads.The start positions of both reads are going to be used
#' for the analysis. Beware that the start position of
#' a reverse read is always higher that the end positition. The \code{GRanges}
#' should at least contain one read.
#'
#' @param seqName a \code{character} string containing the label of the
#' chromosome, present in the \code{GRanges} object, that will be used. The
#' \code{NULL} value is accepted when only one seqname is
#' present in the \code{GRanges}; the only seqname present will be used.
#'
#' @param resultRJMCMC an object of \code{class}
#' "rjmcmcNucleosomes" or "rjmcmcNucleosomesMerge" that contain information
#'  about nucleosome positioning for an entire chromosome.
#'
#' @param extendingSize a positive \code{numeric} or a positive \code{integer}
#' indicating the size of the consensus region used to group closeley
#' positioned nucleosomes.The minimum size of the consensus region is equal to
#' twice the value of the \code{extendingSize} parameter. The numeric will
#' be treated as an integer.
#'
#' @param chrLength a positive \code{numeric} or a positive \code{integer}
#' indicating the lenght of the current chromosome. The length of the
#' chromosome is used to ensure that the consensus positions are all
#' located inside the chromosome.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset containing forward and reverse reads
#' data(reads_demo_01)
#'
#' ## Load dataset containing nucleosome information
#' file_002 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#' pattern = "RJMCMC_seg_02.RDS", full.names = TRUE)
#' nucleosome_info <- readRDS(file_002)
#'
#' ## The function returns 0 when all parameters are valid
#' RJMCMCNucleosomes:::validatePrepMergeParameters(reads = reads_demo_01,
#' seqName = "chr_SYNTHETIC", resultRJMCMC = nucleosome_info,
#' extendingSize = 74, chrLength = 10000000)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{RJMCMCNucleosomes:::validatePrepMergeParameters(
#' reads = c(72400, 72431, 72428, 72429, 72426),
#' resultRJMCMC = NA, extendingSize = 74, chrLength = 10000000)}
#'
#' @author Astrid Deschenes
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
#'
validatePrepMergeParameters <- function(reads, seqName,
                                            resultRJMCMC, extendingSize,
                                            chrLength) {

    ## Validate that the reads is a GRanges
    if (!(class(reads) == "GRanges" )) {
        stop(paste0("reads must be a GRanges"))
    }

    ## Validate that the reads is not empty
    if (length(reads) == 0 ) {
        stop(paste0("reads must be a non-empty GRanges"))
    }

    ## Validate that when seqName is NULL, only one seqname
    ## is present in GRanges
    if (is.null(seqName) &&
        (length(runValue(seqnames(reads))) > 1)) {
        stop(paste0("seqName must be the name of one of the chromosomes ",
                    "present in the GRanges"))
    }

    ## Validate that when not NULL, seqName is a string
    if (!is.null(seqName) && !is.character(seqName)) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## When seqName present, it needs to be the name of one of the
    ## chromosomes present in the GRanges
    if (!is.null(seqName) && is.character(seqName) && !(seqName %in%
                runValue(seqnames(reads)))) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## Validate the resultRJMCMC parameter
    if (!is(resultRJMCMC, "rjmcmcNucleosomesMerge") &&
                !is(resultRJMCMC, "rjmcmcNucleosomes")) {
        stop(paste0("resultRJMCMC must be an object of class",
                    "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'."))
    }

    ## Validate the extendingSize parameter
    if (!(isSingleInteger(extendingSize) || isSingleNumber(extendingSize)) ||
            as.integer(extendingSize) < 1) {
        stop("extendingSize must be a positive integer or numeric")
    }

    ## Validate the chrLength parameter
    if (!(isSingleInteger(chrLength) || isSingleNumber(chrLength)) ||
            as.integer(chrLength) < 1) {
        stop("chrLength must be a positive integer or numeric")
    }

    return(0)
}


#' @title Parameters validation for the \code{\link{plotNucleosomes}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{plotNucleosomes}} function.
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
#'
#' @param xlab a \code{character} string containing the label of the x-axis.
#'
#' @param ylab a \code{character} string containing the label of the y-axis.
#'
#' @param names a \code{vector} of a \code{character} string containing the
#' label of each prediction set. The \code{vector} must be the same length of
#' the \code{nucleosomePositions} \code{list} or 1 in presence of a
#' \code{vector}.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load GRanges dataset
#' data(reads_demo_01)
#'
#' ## Load RJMCMC result
#' data(RJMCMC_result)
#'
#' ## The function returns 0 when all parameters are valid
#' RJMCMCNucleosomes:::validatePlotNucleosomesParameters(nucleosomePositions =
#' RJMCMC_result$mu, reads = reads_demo_01, seqName = "chr_SYNTHETIC",
#' xlab = "position", ylab = "coverage", names = c("test"))
#'
#' ## The function raises an error when at least one paramater is not valid
#' #\dontrun{RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
#' #nucleosomePositions = c("hi"), reads = reads,
#' #xlab = "position", ylab = "coverage", names = c("test"))}
#'
#' #\dontrun{RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
#' #nucleosomePositions = RJMCMC_result$mu, reads = reads_demo_01,
#' #seqName = "chr_SYNTHETIC", xlab = "position", ylab = "coverage",
#' #names = c("test_one", "test_false"))}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
#'
validatePlotNucleosomesParameters <- function(nucleosomePositions,
                    reads, seqName, xlab, ylab, names) {

    ## Validate that nucleosomePositions is a vector
    if(!(class(nucleosomePositions) == "GRanges" ||
            class(nucleosomePositions) == "GRangesList")){
        stop("nucleosomePositions must be a \'GRanges\' or a \'GRangesList\'")
    }

    ## Validate that reads is a GRanges object
    if (!(class(reads) == "GRanges")) {
        stop("reads must be an object of class \'GRanges\'")
    }

    ## Validate that the reads is not empty
    if (length(reads) == 0 ) {
        stop(paste0("reads must be a non-empty GRanges"))
    }

    ## Validate that when seqName is NULL, only one seqname
    ## must be present in GRanges
    if (is.null(seqName) &&
        (length(runValue(seqnames(reads))) > 1)) {
        stop(paste0("seqName must be the name of one of the chromosomes ",
                    "present in the GRanges when more than one chromosome is ",
                    "present"))
    }

    ## When seqName present, it needs to be the name of one of the
    ## chromosomes present in the GRanges
    if (!is.null(seqName) && is.character(seqName) && !(seqName %in%
                    runValue(seqnames(reads)))) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## Validate that when not NULL, seqName is a string
    if (!is.null(seqName) && !is.character(seqName)) {
        stop(paste0("seqName must be a character string corresponding to the ",
                    "name of one of the chromosomes present in the GRanges"))
    }

    ## Validate that xlab is a IRanges object
    if (!is.character(xlab)) {
        stop("xlab must be a character string")
    }

    ## Validate that ylab is a IRanges object
    if (!is.character(ylab)) {
        stop("ylab must be a character string")
    }

    ## Validate that names is a vector of character strings with
    ## length corresponding to nucleosomePositions
    if (!is.null(names)) {
        if (!is.vector(names)) {
            stop("names must be a vector or a list of character strings")
        }
        if (!class(nucleosomePositions) == "GRangesList") {
            if(!(length(names) == 1) || !is.character(names)) {
                stop("names must be a vector of one character string")
            }
        } else {
            if(!(length(names) == length(nucleosomePositions))) {
                stop(paste0("names must be a vector containing the same ",
                    "number of character string as the number of entries ",
                    "in nucleosomesPositions list"))
            }
        }
    }

    return(0)
}

#' @title Parameters validation for the \code{\link{segmentation}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{segmentation}} function.
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
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load synthetic dataset of reads
#' data(syntheticNucleosomeReads)
#'
#' ## Use dataset of reads to create GRanges object
#' sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
#' ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
#' end = syntheticNucleosomeReads$dataIP$end),
#' strand = syntheticNucleosomeReads$dataIP$strand)
#'
#' ## The function returns 0 when all parameters are valid
#' RJMCMCNucleosomes:::validateSegmentationParameters(reads = sampleGRanges,
#' zeta = 147, delta = 30, maxLength = 12000)
#'
#' ## The function raises an error when at least one paramater is not valid
#' #\dontrun{RJMCMCNucleosomes:::validateSegmentationParameters(
#' #reads = c(100), zeta = 147, delta = 30, maxLength = 12000)}
#'
#' #\dontrun{RJMCMCNucleosomes:::validateSegmentationParameters(
#' #reads = sampleGRanges, zeta = "hi", delta = 30, maxLength = 12000)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
#'
validateSegmentationParameters <- function(reads, zeta = 147, delta,
                                            maxLength) {
    # Validate that dataIP is a GRanges
    if(!is(reads,"GRanges"))
    {
        stop("reads must be \'GRanges\' object.")
    }

    ## Validate the zeta parameter
    if (!(isSingleInteger(zeta) || isSingleNumber(zeta)) ||
        as.integer(zeta) < 1) {
        stop("zeta must be a positive integer or numeric")
    }

    ## Validate the delta parameter
    if (!(isSingleInteger(delta) || isSingleNumber(delta)) ||
        as.integer(delta) < 1) {
        stop("delta must be a positive integer or numeric")
    }

    ## Validate the maxLength parameter
    if (!(isSingleInteger(maxLength) || isSingleNumber(maxLength)) ||
        as.integer(maxLength) < 1) {
        stop("maxLength must be a positive integer or numeric")
    }

    return(0)
}

#' @title A internal post treatment function to merge closely positioned
#' nucleosomes, from the same chromosome,
#' identified by the \code{\link{rjmcmc}} function
#'
#' @description A internal helper function which merges closely positioned
#' nucleosomes to rectify the over splitting and provide a more conservative
#' approach. Beware that each chromosome must be treated separatly.
#'
#' The function uses the Bioconductor \code{package} \code{consensusSeeker} to
#' group closely positioned nucleosomes.
#'
#' @param reads a \code{GRanges} containing all forward
#' and reverse reads.The start positions of both reads are going to be used
#' for the analysis. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param resultRJMCMC an object of class 'rjmcmcNucleosomes' or
#' 'rjmcmcNucleosomesMerge' containing informations about nucleosomes.
#'
#' @param extendingSize a positive \code{numeric} or a positive \code{integer}
#' indicating the size of the consensus region used to group closeley
#' positioned nucleosomes.The minimum size of the consensus region is equal to
#' twice the value of the \code{extendingSize} parameter. The numeric will
#' be treated as an integer.
#'
#' @param chrLength a positive \code{numeric} or a positive \code{integer}
#' indicating the lenght of the current chromosome. The length of the
#' chromosome is used to ensure that the consensus positions are all
#' located inside the chromosome.
#'
#' @param minReads a positive \code{integer} or \code{numeric}, the minimum
#' number of reads in a potential canditate region. Non-integer values
#' of \code{minReads} will be casted to \code{integer} and truncated towards
#' zero. Default: 5.
#'
#' @return a \code{array} of \code{numeric}, the updated values of the
#' nucleosome positions. When no nucleosome is present, \code{NULL} is
#' returned.
#'
#' @examples
#'
#' ## Loading dataset
#' data(RJMCMC_result)
#' data(reads_demo_02)
#'
#' ## Results before post-treatment
#' RJMCMC_result$mu
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' postResult <- RJMCMCNucleosomes:::postMerge(reads = reads_demo_02,
#' resultRJMCMC = RJMCMC_result, extendingSize = 80, chrLength = 73500)
#'
#' ## Results after post-treatment
#' postResult
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom consensusSeekeR findConsensusPeakRegions
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqnames genome
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics sapply
#' @keywords internal
#'
postMerge <- function(reads, resultRJMCMC, extendingSize,
                        chrLength, minReads = 5, seqName = NULL)
{
    ## Only keep reads associated to the specified chromosome

    if (!is.null(seqName)) {
        reads <- reads[seqnames(reads) == seqName]
    }else{
        seqName=as.character(unique(seqnames(reads)))
    }
    genomeCur <- as.character(genome(reads)[which(names(genome(reads))
                                                    == seqName)])
    ## Prepare information about reads
    segReads <- list(yF = numeric(), yR = numeric())
    segReads$yF <- start(reads[strand(reads) == "+" ])
    segReads$yR <- end(reads[strand(reads)  == "-" ])

    ## Prepare Seqinfo object using chromosome length
    seqinfo <- Seqinfo(c(seqName), c(chrLength), FALSE, genomeCur)

    finalResult <- NULL

    ## Prepare a GRanges using nucleosome positions
    ## Only apply merge when at least one nucleosome is present
    if (!is.null(resultRJMCMC$mu) && class(resultRJMCMC$mu) == "GRanges"
        && length(resultRJMCMC$mu) > 0) {
        nbMu <- length(resultRJMCMC$mu)
        rjmcmc_peak <- resultRJMCMC$mu
        nbPeaks <- length(rjmcmc_peak)
        names(rjmcmc_peak) <- rep("RJMCMC", nbPeaks)
        rjmcmc_peak$name   <- paste0("RJMCMC_", 1:nbPeaks)

        ## Find nucleosomes present in same regions
        result <- findConsensusPeakRegions(peaks = c(rjmcmc_peak),
                                        chrInfo = seqinfo,
                                        extendingSize = extendingSize,
                                        expandToFitPeakRegion = FALSE,
                                        shrinkToFitPeakRegion = FALSE,
                                        minNbrExp = 1)

        overlapsPeak <- findOverlaps(query = result$consensusRanges,
                                    subject =  rjmcmc_peak)

        allOverlap <- queryHits(overlapsPeak)
        uniqueOverlap <- unique(allOverlap)
        nbOverlap <- length(uniqueOverlap)

        ## Interval used to check for reads
        maxLimit <- 74 + extendingSize
        minLimit <- 74 - extendingSize

        result <- sapply(uniqueOverlap, FUN = function(x) {
            current <- subjectHits(overlapsPeak)[allOverlap == x]
            if(length(current) > 1) {
                ## When more than one nucleosome present, use mean position
                valCentral <- mean(start(resultRJMCMC$mu)[current])
                # - (74 + extendingSize)
                a <- min(start(resultRJMCMC$mu)[current])
                # + (74 - extendingSize)
                b <- max(start(resultRJMCMC$mu)[current])
                ## A miminum number of reads is needed closed to the position
                if(length(segReads$yF[segReads$yF >= (a - maxLimit) &
                    segReads$yF <= (b - minLimit)]) >= minReads &
                    length(segReads$yR[segReads$yR >= (a + minLimit) &
                    segReads$yR <= (b + maxLimit)]) >= minReads) {
                    ## Calculate the new position of the nucleosome
                    newMu <- (mean(segReads$yF[segReads$yF >= (a - maxLimit) &
                                segReads$yF <= (b - minLimit)]) +
                            (mean(segReads$yR[segReads$yR >= (a + minLimit) &
                                    segReads$yR <= (b + maxLimit)]) -
                                mean(segReads$yF[segReads$yF >= (a - maxLimit) &
                            segReads$yF <= (b - minLimit) ]))/2)
                    return(newMu)
                } else {
                    ## Nucleosomes that do not respect condition are flushed
                    return(NA)
                }
            } else {
                ## Return the current mu
                return(start(resultRJMCMC$mu)[current])
            }
        }) # end sapply

        if (!all(is.na(result))) {
            finalResult <- GRanges(seqnames = rep(seqName,
                                        length(as.numeric(
                                                result[!is.na(result)]
                                                ))),
                                    ranges=IRanges(
                                start=round(as.numeric(result[!is.na(result)])),
                                end=round(as.numeric(result[!is.na(result)]))),
                            strand=rep("*", length(
                                as.numeric(result[!is.na(result)]))))
        }
    }

    return(finalResult)
}



#' @title Run \code{\link{rjmcmc}} on multiples segments and merge results.
#'
#' @description Run \code{\link{rjmcmc}} on a segment that is contained
#' in a \code{list} of segments. Files generated by the function are all saved
#' in a directory specified by user. A RData log file is created when a
#' segment has been run while the result is saved in a RDS file.
#'
#' If the same output directory is used more than once, the
#' \code{\link{rjmcmc}} won't be called for segments that have the à
#' corresponding RData log file.
#'
#' @param p a \code{integer}, the position of the segment to treat in a
#' \code{list} of \code{GRanges}.
#'
#' @param seg a \code{list} a \code{GRanges} containing the segments
#' to be process.
#'
#' @param niter a positive \code{integer} or \code{numeric}, the
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
#' @param adaptNbrIterations a \code{logical} indicating if the number
#' of iterations must be modified in function of the number of reads.
#'
#' @param vSeed a \code{integer}. A seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: -1.
#'
#' @param saveAsRDS a \code{logical}. When \code{TRUE}, a RDS file containing
#' the complete output of the c++ rjmcmc() function is created.
#' Default : \code{FALSE}.
#'
#' @param maxLength a positive \code{integer} or \code{numeric}, the
#' length of each segment.
#'
#' @return 0.
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
#' seg <- segmentation(sampleGRanges, zeta = 147, delta = 50, maxLength = 1000)
#' \dontrun{
#' dir.create("out")
#' dir.create("out/done")
#' dir.create("out/results")
#'
#' runCHR(p=1, seg=seg, niter=1000, kmax=330,lambda=3,
#'             ecartmin=147, ecartmax=297, minReads=5)}
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom GenomicRanges strand
#' @keywords internal
#'
runCHR <- function(p, seg, niter, kmax, lambda,
                    ecartmin, ecartmax, minReads, adaptNbrIterations,
                    vSeed=-1, saveAsRDS=FALSE, dirOut="out")
{
    dirDone <- paste0(dirOut, "/done")
    dirResults <- paste0(dirOut, "/results")
    validateDirectoryParameters(dirDone)
    validateDirectoryParameters(dirResults)

    nameDone <- paste0(dirDone,"/rjmcm_seg_", p, ".RData")

    if (!file.exists(nameDone) ) {
        nameRDS  <- paste0(dirResults,"/rjmcmc_seg_", p, ".rds")
        print(paste0("Doing: ", nameRDS))
        listeSeg <- rjmcmc(reads = seg[[p]], nbrIterations = niter,
                            kMax = kmax, lambda=lambda, minInterval = ecartmin,
                            maxInterval = ecartmax, minReads = minReads,
                            vSeed = vSeed, saveAsRDS = saveAsRDS,
                            adaptIterationsToReads = adaptNbrIterations)

        if(all(is.na(listeSeg)))
        {
            print(paste0("Seg NA ", nameRDS))
        }
        saveRDS(listeSeg, nameRDS)

        print(paste0("Done: ", nameRDS))
        save(nameRDS, file=nameDone)
    }

    return(0)
}


