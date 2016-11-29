#' @title Formated output of predicted nucleosomes
#'
#' @description Generated a formated output of a list marked as
#' an \code{rjmcmcNucleosomesBeforeAndAfterPostTreatment} class
#'
#' @method print rjmcmcNucleosomesBeforeAndAfterPostTreatment
#'
#' @param x the output object from \code{rjmcmcCHR}
#' function to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @return an object of class
#' \code{rjmcmcNucleosomesBeforeAndAfterPostTreatment}
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
#' ## Print result
#' \dontrun{print(result)}
#'
#' @author Astrid Deschenes
#' @export
print.rjmcmcNucleosomesBeforeAndAfterPostTreatment <- function(x, ...) {
    # Print title before printing the content
    cat(paste0("\nRJMCMCNucleosomes - Predicted nucleosomes Before and ",
                "After Post-Treatment\n"))
    cat("BEFORE POST-TREATMENT\n")
    cat("Number of nucleosomes:\n")
    print(x$k, ...)
    cat("\nNucleosomes positions:\n")
    print(x$mu, ...)
    cat("\nAFTER POST-TREATMENT\n")
    cat("Number of nucleosomes:\n")
    print(x$kPost, ...)
    cat("\nNucleosomes positions:\n")
    print(x$muPost, ...)
    invisible(x)
}
