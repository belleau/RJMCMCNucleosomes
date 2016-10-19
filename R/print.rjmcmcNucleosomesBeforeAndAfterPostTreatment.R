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
#' @examples
#'
#' ## TODO
#'
#' @author Astrid Deschenes
#' @export
print.rjmcmcNucleosomesBeforeAndAfterPostTreatment <- function(x, ...) {
    # Print title before printing the content
    cat(paste0("\nRJMCMCNucleosomes - Predicted nucleosomes Before and ",
               "After Post-Treatment\n"))
    cat("BEFORE POST-TREATMENT\n")
    cat("\nNumber of nucleosomes:\n")
    print(x$k, ...)
    cat("\nNucleosomes positions:\n")
    print(x$mu, ...)
    cat("AFTER POST-TREATMENT\n")
    cat("\nNumber of nucleosomes:\n")
    print(x$kPost, ...)
    cat("\nNucleosomes positions:\n")
    print(x$muPost, ...)
    invisible(x)
}
