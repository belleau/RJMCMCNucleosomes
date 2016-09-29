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
#' of nucleosomes.
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
#' \item \code{kMax} a \code{integer}, the maximum number of nucleosomes
#' obtained during the iteration process.
#' \item \code{sigmaf} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#' \item \code{sigmar} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#' \item \code{delta} a \code{vector} of \code{numeric} of length
#' \code{k}, the distance between the maxima of the forward and reverse reads
#' position densities for each nucleosome.
#' \item \code{df} a \code{vector} of \code{numeric} of length
#' \code{k}, the degrees of freedom for each nucleosome.
#' \item \code{w} a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' \item \code{qmu} a \code{matrix} of \code{numerical} with a number of rows
#' of \code{k}, the 2.5\% and 97.5\% quantiles of each \code{mu}.
#' \item \code{qsigmaf} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the variance of the
#' forward reads for each nucleosome.
#' \item \code{qsigmar} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the variance the
#' reverse reads for each nucleosome.
#' \item \code{qdelta} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the distance between
#' the maxima of the forward and reverse reads
#' position densities for each nucleosome.
#' \item \code{qdf} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the degrees of freedom
#' for each nucleosome.
#' \item \code{qw} a \code{matrix} of \code{numerical} with a number of rows
#' of \code{k}, the 2.5\% and 97.5\% quantiles of the weight for each
#' nucleosome.
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
#' ## TODO : a faire avec nouveau format
#'
#' ## Print the final estimation of the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' ## Print the maximum number of nucleosomes obtained during the iteration
#' ## process
#' result$kMax
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
#     nf              <- length(startPosForwardReads)
#     nr              <- length(startPosReverseReads)
#     nbrReads        <- nf + nr
#     if (adaptIterationsToReads) {
#         nbrIterations <- ifelse(nbrReads <= 12, 1000, nbrIterations)
#     }

    resultRJMCMC <- rjmcmcNucleo(startPosForwardReads, startPosReverseReads,
                                 nbrIterations, kMax, lambda,
                                 minInterval, maxInterval, minReads,
                                 adaptIterationsToReads, vSeed)

    # Find k value with the maximum of iterations
    iterPerK <- data.frame(k=resultRJMCMC$K, it=resultRJMCMC$it)
    sumIterPerK <- aggregate(it ~ k, data=iterPerK, sum)
    maxRow <- which.max( sumIterPerK[,"it"])
    k <- sumIterPerK$k[maxRow]

    # Find mu values associated to the k value
    mu <- resultRJMCMC$muHat[k,]

    result <- list(
        call = cl,
        k = k,
        mu = mu[1:k],
        k_max = resultRJMCMC$KMax
    )

    class(result)<-"rjmcmcNucleosomes"

#     result <- list(
#         call    = cl,
#         K       = k,
#         k       = km,
#         mu      = mu_hat,
#         sigmaf  = sigmaf_hat,
#         sigmar  = sigmar_hat,
#         delta   = delta_hat,
#         df      = df_hat,
#         w       = w_hat,
#         qmu     = qmu,
#         qsigmaf = qsigmaf,
#         qsigmar = qsigmar,
#         qdelta  = qdelta,
#         qdf     = qdf,
#         qw      = qw
#     )

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
#' directoryWithRDSFiles <- system.file("extdata", package = "RJMCMC")
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


#' @title Merge nucleosome information for selected RDS files.
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
#' RDSFiles <- dir(system.file("extdata", package = "RJMCMC"),
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
#' @author Pascal Belleau, Astrid Deschenes
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

