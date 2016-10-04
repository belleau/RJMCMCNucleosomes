#' @title Interface for the RJMCMC nucleosome mapping method in C++
#'
#' @description Function that calls the core of the nucleosome positioning
#' mapping function that is implemented in C++.
#'
#' @param k a positive \code{integer}, the number of nucleosomes.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
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
#' @return TODO
#'
#' @examples
#'
## Loading dataset
#' data(reads_demo)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmcNucleo(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 1000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5,
#'          adaptIterationsToReads = TRUE, vSeed = -1)
#'
#' ## Print the final estimation of the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' @author Pascal Belleau
#' @importFrom Rcpp evalCpp
#' @useDynLib RJMCMCNucleosomes, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @keywords internal
#'
rjmcmcNucleo <- function(startPosForwardReads,
                              startPosReverseReads,
                              nbrIterations, kMax, lambda,
                              minInterval, maxInterval,
                              minReads = 5L,
                              adaptIterationsToReads = TRUE, vSeed = -1) {
    .Call('RJMCMCNucleosomes_rjmcmcNucleo',
          PACKAGE = 'RJMCMCNucleosomes',
          startPosForwardReads, startPosReverseReads,
          nbrIterations, kMax, lambda, minInterval,
          maxInterval, minReads, adaptIterationsToReads, vSeed)
}


#' @title Birth Submove Probability
#'
#' @description Calculation of the birth submove probability of adding a new
#' nucleosome using a truncated Poisson distribution.
#'
#' @param k a positive \code{integer}, the number of nucleosomes.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
#'
#' @param kMax a positive \code{integer}, the maximum number of nucleosomes
#' authorized. When \code{k} is equal or superior to \code{kMax}, the
#' returned value is \code{0}. Default: \code{30}.
#'
#' @return a \code{numeric} value. The value \code{0} when \code{k} is equal
#' or superior to \code{kMax} or when \code{k} is equal to \code{1}.
#'
#' @examples
#'
#' ## Return the birth submove probability
#' RJMCMCNucleosomes:::Bk(k = 14L, lambda = 1L, kMax = 30L)
#'
#' ## Zero is returned when k = 1
#' RJMCMCNucleosomes:::Bk(k = 1L, lambda = 3L, kMax = 30L)
#'
#' ## Zero is returned when k is superior to kMax
#' RJMCMCNucleosomes:::Bk(k = 31L, lambda = 2L, kMax = 30L)
#'
#' @author Rawane Samb
#' @importFrom stats dpois
#' @keywords internal
Bk <- function(k, lambda, kMax = 30) {
    ifelse((k >= kMax), 0,
            0.5 * min(1, lambda/ (k+1))) # min(1, dpois(k + 1, lambda) / dpois(k, lambda)))
}


#' @title Random deviate from a truncated normal distribution
#'
#' @description Generate a random deviate value from a normal distribution.
#' The returned value is included inside a specified range
#' ]\code{minValue},\code{maxValue}[
#' specified by user. The mean and variance of the normal distribution is
#' also specified by user.
#'
#' @param mu a \code{numeric} value, the mean of the normal distribution.
#'
#' @param sigma a non-negative \code{numeric}, the variance of the normal
#' distribution.
#'
#' @param minValue a \code{numeric} value, the inferior boundary of the
#' range in which the output value must be located. The returned value has to
#' be superior to \code{minValue}.
#'
#' @param maxValue a \code{numeric} value, the superior boundary of the range
#' in which the output value must be located. The returned value has to be
#' inferior to \code{maxValue}.
#'
#' @return a \code{numeric} which is superior to \code{minValue} and inferior
#' to \code{maxValue}.
#'
#' @examples
#'
#' ## Set seed to replicate results
#' set.seed(3333)
#'
#' ## Return a value, between 1000 and 3000, generated froma a normal
#' ## distribution with a mean of 2000 and a variance of 30.
#' RJMCMCNucleosomes:::tnormale(2000, 30, 1000, 30000)
#'
#' @author Rawane Samb
#' @importFrom stats rnorm
#' @keywords internal
tnormale <- function(mu, sigma, minValue, maxValue) {
    repeat {
        y <- rnorm(1, mu, sd = sqrt(sigma))
        if (y > minValue & y < maxValue) break()
    }
    return(y)
}


#' @title Student Mixture Model
#'
#' @description Generation a value from a Student Mixture distribution.
#'
#' @param i an \code{integer}, a count parameter.
#'
#' @param k a positive \code{integer} value, the number of nucleosomes in the
#' analyzed region.
#'
#' @param weight a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each
#' nucleosome. The sum of all \code{weight} values must be equal to \code{1}.
#' The length of \code{weight} must be equal to \code{k}.
#'
#' @param mu a \code{vector} of positive \code{integer} of length \code{k},
#' the positions of all the nucleosomes in the analyzed region. The length
#' of \code{weight} must be equal to \code{k}.
#'
#' @param sigma a \code{vector} of \code{numeric} of length \code{k}, the
#' variance for each nucleosome. The length of \code{sigma} must be equal
#' to \code{k}.
#'
#' @param dfr a \code{vector} of \code{numeric} of length \code{k}, the degrees
#' of freedom for each nucleosome. The length of \code{dfr} must
#' be equal to \code{k}.
#'
#' @return a \code{numerical}, the value generated from a Student Mixture
#' distribution.
#'
#' @examples
#'
#' ## Return a value generated from a student mixture
#' RJMCMCNucleosomes:::student.mixture(i = 1L, k = 4L, weight = c(0.1, 0.3, 0.34, 0.26),
#' mu = c(12L, 15L, 25L, 44L), sigma = c(4, 7, 6, 5), dfr = c(5L, 3L, 12L, 4L))
#'
#' @importFrom stats runif rt
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
student.mixture <- function(i, k, weight, mu, sigma, dfr) {
    # Adding zero to the weight vector and calculating the cumulative sums
    sumWeight <- cumsum(c(0, weight))

    u <- runif(1, 0, 1)

    # Get the maximal position where the sum of weight is inferior to u
    position <- max(which(sumWeight < u))

    return(mu[position] + sqrt(sigma[position]) * rt(1, dfr[position]))
}


#' @title Normal Mixture Model
#'
#' @description Generation a value from a Normal Mixture distribution.
#'
#' @param i a \code{integer}, a count parameter.
#'
#' @param k a positive \code{integer} value, the number of nucleosomes in the
#' analyzed region.
#'
#' @param weight a \code{vector} of length \code{k}, the weight for each
#' nucleosome. The sum of all \code{weight} values must be equal to \code{1}.
#' The length of \code{weight} must be equal to \code{k}.
#'
#' @param mu a \code{vector} of positive \code{integer} of length \code{k},
#' the positions of all the nucleosomes in the analyzed region. The length
#' of \code{weight} must be equal to \code{k}.
#'
#' @param sigma a \code{vector} of length \code{k}, the variance for each
#' nucleosome. The length of \code{sigma} must be equal to \code{k}.
#'
#' @return a \code{numerical}, the value generated from a Normal Mixture
#' distribution.
#'
#' @examples
#'
#' ## Return a value generated from a normal mixture
#' RJMCMCNucleosomes:::normal.mixture(i = 1L, k = 4L, weight = c(0.2, 0.3, 0.24, 0.26),
#' mu = c(12L, 15L, 25L, 44L), sigma = c(4, 7, 6, 5))
#'
#' @importFrom stats runif
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
normal.mixture <- function(i, k, weight, mu, sigma) {
    # Adding zero to the weight vector and calculating the cumulative sums
    sumWeight <- cumsum(c(0, weight))

    u <- runif(1, 0, 1)

    # Get the maximal position where the sum of weight is inferior to u
    position <- max(which(sumWeight < u))

    return(rnorm(1, mu[position], sd = sqrt(sigma[position])))
}


#' @title Prior density of \eqn{mu}
#'
#' @description Computes the prior density of \eqn{mu} conditionally to
#' the number of nucleosomes.
#'
#' For more information on the calculation of the prior density of \eqn{mu},
#' see Proposotion 1 and equation (11) of the cited article.
#'
#' @param mu a \code{vector} of positive \code{integer} containing the
#' positions of all nucleosomes.
#'
#' @param readPositions a \code{vector} of positive \code{integer}
#' corresponding to the
#' positions of all reads, including forward and reverse strands. The
#' values insinde \code{readPositions} must be sorted.
#'
#' @return  a \code{numeric}, the exact prior density of \code{mu} given the
#' number of nucleosomes.
#'
#' @references Samb R., Khadraoui K., Belleau P., Deschenes A., Lakhal L.
#' and Droit A. Using
#' informative Multinomial-Dirichlet prior in a t-mixture with
#' reversible jump estimation of nucleosome positions for genome-wide
#' profiling. Statistical Applications in Genetics and Molecular Biology.
#' Volume 14, Issue 6, Pages 517-532, ISSN (Online) 1544-6115,
#' ISSN (Print) 2194-6302, DOI: 10.1515/sagmb-2014-0098, December 2015
#'
#' @examples
#'
#' ## Sorted vector of read positions, including forward and reverse
#' readPositions <- c(9909L, 9928L, 9935L, 26603L, 26616L, 26632L, 26636L,
#' 26640L, 44900L, 44902L, 44909L,  44910L, 44910L, 44918L,
#' 44924L, 44931L, 44935L, 44942L, 44946L)
#'
#' ## Position of the group of nucleosomes
#' mu <- c(10000L, 26700L, 45000L)
#'
#' ## Calculation of the exact prior density of mu
#' density <- RJMCMCNucleosomes:::priorMuDensity(mu, readPositions)
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
priorMuDensity <- function(mu, readPositions) {
    ## Get the number of nucleosomes
    k <- length(mu)

    ## Create a matrix used in the calculation of the priors
    basicMatrix <- matrix(0L, nrow = k, ncol = k)
    for (i in 1:k) {
        basicMatrix[i, i] <- 1L
    }
    if (k > 1) {
        for (i in 2:k) {
            basicMatrix[i , i - 1] <- -1L
        }
    }
    omega <- t(basicMatrix) %*% basicMatrix

    ## Calculating the range (R)
    R <- max(readPositions) - min(readPositions)

    ## Calculating the mean (E)
    E <- (max(readPositions) + min(readPositions))/2
    tau <- 1/R^2
    M <- rep(E, k)
    const <- (pi/(2*tau))^{-k/2}

    ## The calculation of the prior
    ## Equation 11 in the cited article
    return(const * exp(-(tau/2) * (t(mu - M) %*% omega %*% (mu - M))))
}


#' @title Element with the hightest number of occurences
#'
#' @description Returned the \code{integer} with the highest number of
#' occurences in a \code{vector}.
#' When more than one \code{integer} have the highest number of occurences,
#' \code{NA} is returned.
#'
#' @param sample a \code{numeric} \code{vector} (of positive \code{integer}
#' values). If the elements of \code{sample} are \code{numeric} but not
#' \code{integer}, the elements are truncated by \code{as.integer}.
#'
#' @return  a \code{integer} with the highest number of occurences or
#' \code{NA} when more than one \code{integer} have the highest number
#' of occurences.
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
#' @examples
#'
#' ## Return the element with the hightest number of occurence
#' data01 <- c(1L, 2L, 5L, 10L, 5L, 10L, 5L)
#' RJMCMCNucleosomes:::elementWithHighestMode(data01)
#'
#' data02 <- c(3L, 6L, 4L, 3L, 6L)
#' RJMCMCNucleosomes:::elementWithHighestMode(data02)
#'
elementWithHighestMode <- function(sample) {
    tabsample <- tabulate(sample)
    maxOccurence <- tabsample == max(tabsample)
    ifelse(sum(maxOccurence) == 1, which(maxOccurence), NA)
}


#' @title Parameters validation for the \code{\link{rjmcmc}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{rjmcmc}} function.
#'
#' @param startPosForwardReads a \code{vector} of positive \code{integer}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of positive \code{integer}, the
#' positions of all the reverse reads. Beware that the start position of
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
#' ## The function returns 0 when all paramaters are valid
#' RJMCMCNucleosomes:::validateRJMCMCParameters(startPosForwardReads = c(72400,
#' 72431, 72428, 72429, 72426), startPosReverseReads = c(72520, 72523, 72521,
#' 72533, 72511), nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1,
#' minInterval = 100, maxInterval = 200, adaptIterationsToReads = TRUE,
#' vSeed = 100)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{RJMCMCNucleosomes:::validateRJMCMCParameters(startPosForwardReads =
#' c(72400, 72431, 72428, 72429, 72426), startPosReverseReads = NA,
#' nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1, minInterval = 100,
#' maxInterval = 200, adaptIterationsToReads = TRUE, vSeed = -1)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRJMCMCParameters <- function(startPosForwardReads, startPosReverseReads,
                                    nbrIterations, kMax, lambda,
                                    minInterval, maxInterval, minReads,
                                    adaptIterationsToReads, vSeed) {
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

    ## Validate that the startPosForwardReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosForwardReads) || !is.numeric(startPosForwardReads)
        || length(startPosForwardReads) < 1)
    {
        stop(paste0("startPosForwardReads must be a non-empty vector of ",
                    "numeric values."))
    }

    ## Validate that the startPosReverseReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosReverseReads) || !is.numeric(startPosReverseReads)
        || length(startPosReverseReads) < 1)
    {
        stop(paste0("startPosReverseReads must be a non-empty vector of ",
                    "numeric values."))
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
#' TODO: to modify when new format will be set.
#'
#' @param arrayOfFiles a \code{array}, the name of each file that must be
#' used to merge nucleosome information.
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
#' ## Loading two files containing nucleosomes informations for two sections of
#' ## the same chromosome
#' file_100 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#'                 pattern = "yeastRes_Chr1_Seg_100.rds",
#'                 full.names = TRUE)
#'
#' file_101 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#'                 pattern = "yeastRes_Chr1_Seg_101.rds",
#'                 full.names = TRUE)
#'
#' ## Merging nucleosomes informations from the two files
#' result <- RJMCMCNucleosomes:::mergeAllRDSFiles(c(file_100, file_101))
#'
#' @importFrom methods is
#' @author Pascal Belleau, Astrid Deschênes
#' @keywords internal
#'
mergeAllRDSFiles <- function(arrayOfFiles) {

    ## Create list that will contain data from all files
    result        <- list()
    result$k      <- 0
    result$mu     <- array(dim = c(0))
    result$sigmaf <- array(dim = c(0))
    result$sigmar <- array(dim = c(0))
    result$delta  <- array(dim = c(0))
    result$df     <- array(dim = c(0))

    ## Extract information from each file
    for (fileName in arrayOfFiles) {
        data <- readRDS(file = fileName)
        ## Only use data from rjmcmcNucleosomes or rjmcmcNucleosomesMerge class
        if (is(data, "rjmcmcNucleosomesMerge") ||
                is(data, "rjmcmcNucleosomes")) {
            result$k      <- result$k + data$k
            result$mu     <- c(result$mu, data$mu)
            result$sigmaf <- c(result$sigmaf, data$sigmaf)
            result$sigmar <- c(result$sigmar, data$sigmar)
            result$delta  <- c(result$delta, data$delta)
            result$df     <- c(result$df, data$df)
        }
    }

    ## Ensure that all values are ordered in ascending order of mu
    newOrder      <- order(result$mu)
    result$mu     <- result$mu[newOrder]
    result$sigmaf <- result$sigmaf[newOrder]
    result$sigmar <- result$sigmar[newOrder]
    result$delta  <- result$delta[newOrder]
    result$df     <- result$df[newOrder]

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
#' pattern = "yeastRes_Chr1_Seg_002.rds", full.names = TRUE)
#'
#' ## Testing using a real file
#' RJMCMCNucleosomes:::validateRDSFilesParameters(c(file_test))
#'
#' @author Astrid Deschenes
#' @keywords internal
#'
validateRDSFilesParameters <- function(RDSFiles) {

    ## Validate the RDSFiles parameters
    if (is.null(RDSFiles) || is.na(RDSFiles)) {
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
#' @param startPosForwardReads a \code{vector} of positive \code{integer}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of positive \code{integer}, the
#' positions of all the reverse reads. Beware that the start position of
#' a reverse read is always higher that the end position.
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
#' data(reads_demo)
#'
#' ## Load dataset containing nucleosome information
#' file_002 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
#' pattern = "yeastRes_Chr1_Seg_002.rds", full.names = TRUE)
#' nucleosome_info <- readRDS(file_002)
#'
#' ## The function returns 0 when all parameters are valid
#' RJMCMCNucleosomes:::validatePrepMergeParameters(startPosForwardReads =
#' reads_demo$readsForward, startPosReverseReads = reads_demo$readsReverse,
#' resultRJMCMC = nucleosome_info, extendingSize = 74, chrLength = 10000000)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{RJMCMCNucleosomes:::validatePrepMergeParameters(startPosForwardReads = c(72400,
#' 72431, 72428, 72429, 72426), startPosReverseReads = c(72522, 72531, 72528,
#' 72559, 72546), resultRJMCMC = NA, extendingSize = 74, chrLength = 10000000)}
#'
#' @author Astrid Deschenes
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
#'
validatePrepMergeParameters <- function(startPosForwardReads,
                                            startPosReverseReads,
                                            resultRJMCMC, extendingSize,
                                            chrLength) {

    ## Validate that the startPosForwardReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosForwardReads) || !is.numeric(startPosForwardReads)
        || length(startPosForwardReads) < 1)
    {
        stop(paste0("startPosForwardReads must be a non-empty vector of ",
                    "numeric values."))
    }

    ## Validate that the startPosReverseReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosReverseReads) || !is.numeric(startPosReverseReads)
        || length(startPosReverseReads) < 1)
    {
        stop(paste0("startPosReverseReads must be a non-empty vector of ",
                    "numeric values."))
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
#' @param startPosFrowardReads a \code{vector} of \code{numeric}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of \code{numeric}, the
#' start position of all the reverse reads. Beware that the start position of
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
#'
#' @return a \code{array} of \code{numeric}, the updated values of the
#' nucleosome positions.
#'
#' @examples
#'
#' ## Loading dataset
#' data(RJMCMC_result)
#' data(reads_demo)
#'
#' ## Results before post-treatment
#' RJMCMC_result$mu
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' postResult <- RJMCMCNucleosomes:::postMerge(startPosForwardReads =
#' reads_demo$readsForward, startPosReverseReads = reads_demo$readsReverse,
#' resultRJMCMC = RJMCMC_result, extendingSize = 80, chrLength = 73500)
#'
#' ## Results after post-treatment
#' postResult
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom consensusSeekeR findConsensusPeakRegions
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqnames
#' @importFrom S4Vectors queryHits subjectHits
#' @keywords internal
#'
postMerge <- function(startPosForwardReads, startPosReverseReads,
                                resultRJMCMC, extendingSize, chrLength, minReads = 5)
{
    ## Prepare information about reads
    segReads <- list(yF = numeric(), yR = numeric())
    segReads$yF <- startPosForwardReads
    segReads$yR <- startPosReverseReads

    ## Prepare Seqinfo object using chromosome length
    seqinfo <- Seqinfo(c("chrI"), c(chrLength), FALSE, "mock1")

    ## Prepare first GRanges using nucleosome positions
    nbMu <- length(resultRJMCMC$mu)
    rjmcmc_peak <- GRanges(seqnames = rep('chrI', nbMu),
                        IRanges(resultRJMCMC$mu, resultRJMCMC$mu),
                        seqinfo = seqinfo)
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

    ## Treat each overlapping region separatly
    newMu <- array(dim = nbOverlap)
    cpt <- 1L
    for(position in 1:nbOverlap){
        ## Extract nucleosomes present in the current overlapping region
        current <- subjectHits(overlapsPeak)[queryHits(overlapsPeak) ==
                                                    uniqueOverlap[position]]
        if(length(current) > 1) {
            ## When more than one nucleosome present, use mean position
            valCentral <- mean(resultRJMCMC$mu[current])
            a <- min(resultRJMCMC$mu[current]) # - (74 + extendingSize)
            b <- max(resultRJMCMC$mu[current]) # + (74 - extendingSize)

            if(length(segReads$yF[segReads$yF >= (a - maxLimit) &
                                    segReads$yF <= (b - minLimit)]) >= minReads &
                length(segReads$yR[segReads$yR >= (a + minLimit) &
                                    segReads$yR <= (b + maxLimit)]) >= minReads) {
                ## Calculate the new position of the nucleosome
                newMu[cpt] <- (mean(segReads$yF[segReads$yF >= (a - maxLimit) &
                                    segReads$yF <= (b - minLimit)]) +
                            (mean(segReads$yR[segReads$yR >= (a + minLimit) &
                                    segReads$yR <= (b + maxLimit)]) -
                            mean(segReads$yF[segReads$yF >= (a - maxLimit) &
                                    segReads$yF <= (b - minLimit) ]))/2)
                cpt <- cpt + 1L
            }

            ## Nucleosomes not respecting the condition are flushed
        } else {
            newMu[cpt] <- resultRJMCMC$mu[current]
            cpt <- cpt + 1L
        }
    }

    return(as.numeric(newMu[!is.na(newMu)]))
}
