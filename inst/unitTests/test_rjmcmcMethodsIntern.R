###################################################
# Created by Astrid Deschenes
# 2015-06-12
###################################################

###################################################
## Test the rjmcmcMethodsIntern.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMCNucleosomes" )
}

### }}}


file_002 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                pattern = "yeastRes_Chr1_Seg_002.rds",
                full.names = TRUE)

data_002 <- readRDS(file_002)

data(RJMCMC_result)
data(reads_demo)


#########################################################
## elementWithHighestMode() function
#########################################################

test.elementWithHighestMode_results_with_different_vectors <- function() {
    set.seed(101)
    obs <- lapply(list(test1 = c(1, 2, 3, 3, 21, 22),
                    test2 = c(2, 4, 2, 4, 4, 0),
                    test3 = c(12, 13, 13, 12, 1),
                    test4 = c(23, 22, 21, 20, 22, 22, 21, 21, 22)),
                  FUN = function(x) {RJMCMCNucleosomes:::elementWithHighestMode(x)})
    exp <- list(test1 = 3L, test2 = 4L, test3 = NA, test4 = 22L)
    message <- paste0(" elementWithHighestMode_results_with_different_vectors() ",
                    "- Not all tested vectors generated",
                    " expected values.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## priorMuDensity() function
#########################################################

# test.priorMuDensity_with_1000_reads <- function() {
#     set.seed(101)
#     k <- 3L
#     nbrReads <- 500L
#     mu <- c(10000L, 26700L, 45000L)
#     sigma <- rep(400L, k)
#     delta  <- rep(147, k)
#     weight <- c(0.3, 0.2, 0.5)
#     readsForward <- sapply(1:nbrReads, RJMCMCNucleosomes:::normal.mixture, k = k,
#                         w = weight, mu = mu - delta/2, sigma = sigma)
#     readsReverse <- sapply(1:nbrReads, RJMCMCNucleosomes:::normal.mixture, k = k,
#                         w = weight, mu = mu - delta/2, sigma = sigma)
#     reads <- sort(c(readsForward, readsReverse))
#     obs <- RJMCMCNucleosomes:::priorMuDensity(mu, reads)
#     exp <- 8.0883349761e-15
#     message <- paste0(" priorMuDensity_with_1000_reads() ",
#                     "- The result is not the expected value.")
#     checkEqualsNumeric(obs, exp, msg = message)
# }

test.priorMuDensity_results_with_various_values_of_mu <- function() {
    set.seed(101)
    obs <- mapply(list(A=c(10200, 10300), B=c(10108, 10206, 10222),
                    C=c(10333, 10455, 10899)),
                    FUN = function(x) { RJMCMCNucleosomes:::priorMuDensity(x,
                                            10000:11000) })
    exp <- c(6.0557145969e-7, 4.6807063829e-10, 4.5053092518e-10)
    message <- paste0(" priorMuDensity_results_with_various_values_of_mu() ",
                    "- Not all tested data with various mu generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}


#########################################################
## validatePrepMergeParameters() function
#########################################################

## Test the result when startPosForwardReads is NA
test.validatePrepMergeParameters_startPosForwardReads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
                        startPosForwardReads = NA,
                        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                        72429.24, 72426.08),
                        resultRJMCMC = data_002, extendingSize = 11,
                        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                    "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_NA() ",
                      "- NA for startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is empty
test.validatePrepMergeParameters_startPosForwardReads_empty <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_empty() ",
                      "- empty startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is not number
test.validatePrepMergeParameters_startPosForwardReads_not_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c("A", "B"),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_not_number() ",
                      "- not number startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is NA
test.validatePrepMergeParameters_startPosReverseReads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = NA,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_NA() ",
                      "- NA for startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is empty
test.validatePrepMergeParameters_startPosReverseReads_empty <- function() {
    seqinfo <- GenomeInfoDb::Seqinfo(c("chr1"), c(1000000), NA, "mock1")
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_empty() ",
                      "- empty startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is not number
test.validatePrepMergeParameters_startPosReverseReads_not_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c("A", "B"),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_not_number() ",
                      "- not number startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is NA
test.validatePrepMergeParameters_resultRJMCMC_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = NA, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_not_number() ",
                      "- NA resultRJMCMC did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is a number
test.validatePrepMergeParameters_resultRJMCMC_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = 33, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_resultRJMCMC_number() ",
                      "- number resultRJMCMC did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is a string
test.validatePrepMergeParameters_nbBase_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = "ALLO",
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_number() ",
                      "- string nbBase did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is an array
test.validatePrepMergeParameters_nbBase_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = c(10, 11),
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_string() ",
                      "- array nbBase did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is a string
test.validatePrepMergeParameters_chrLength_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = "5000"), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- string chrLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is an array
test.validatePrepMergeParameters_chrLength_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = c(100, 200)), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- array chrLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validatePrepMergeParameters_all_valid <- function() {
    obs <- RJMCMCNucleosomes:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = 200000)
    exp <- 0
    message <- paste0(" test.validatePrepMergeParameters_all_valid() ",
                      "- All valid parameters did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateRJMCMCParameters() function
#########################################################

## Test the result when nbrIterations is NA
test.validateRJMCMCParameters_nbrIterations_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                    72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                    72429.24, 72426.08),
        nbrIterations = NA,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_NA() ",
                      "- NA for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is zero
test.validateRJMCMCParameters_nbrIterations_zero <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 0,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Zero for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is negative
test.validateRJMCMCParameters_nbrIterations_negative <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = -1,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Negative value for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is NA
test.validateRJMCMCParameters_kMax_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = NA, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_NA() ",
                      "- NA value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is zero
test.validateRJMCMCParameters_kMax_zero <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters (
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 0, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_zero() ",
                      "- Zero value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is negative
test.validateRJMCMCParameters_kMax_negative <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters (
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = -1, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_negative() ",
                      "- Negative value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is NA
test.validateRJMCMCParameters_minReads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = NA, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is zero
test.validateRJMCMCParameters_minReads_zero <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 0, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is negative
test.validateRJMCMCParameters_minReads_negative <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 3, lambda = 1, minReads = -1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is NA
test.validateRJMCMCParameters_lambda_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = NA, minReads = 2, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is zero
test.validateRJMCMCParameters_lambda_zero <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 0, minReads = 3, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is negative
test.validateRJMCMCParameters_lambda_negative <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = -1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is NA
test.validateRJMCMCParameters_startPosForwardReads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = NA,
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosForwardReads_NA() ",
                      "- NA value for startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is empty array
test.validateRJMCMCParameters_startPosForwardReads_empty <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosForwardReads_empty() ",
                        "- Empty array for startPosForwardReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is NA
test.validateRJMCMCParameters_startPosReverseReads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = NA,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosReverseReads_NA() ",
                        "- NA value for startPosReverseReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is empty array
test.validateRJMCMCParameters_startPosReverseReads_empty_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosReverseReads_empty_array() ",
                        "- Empty array for startPosReverseReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when adaptIterationsToReads is string
test.validateRJMCMCParameters_adaptIterationsToReads_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = "allo"), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_string() ",
                        "- String for adaptIterationsToReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when adaptIterationsToReads is number
test.validateRJMCMCParameters_adaptIterationsToReads_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = 33), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_number() ",
                        "- Number value for adaptIterationsToReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when vSeed is not a number
test.validateRJMCMCParameters_vSeed_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE, vSeed = "Allo"), error=conditionMessage)
    exp <- "vSeed must be a numeric value."
    message <- paste0(" test.validateRJMCMCParameters_vSeed_number() ",
                      "- String value for vSeed did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validateRJMCMCParameters_all_valid <- function() {
    obs <- RJMCMCNucleosomes:::validateRJMCMCParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = TRUE, vSeed = 1002)
    exp <- 0
    message <- paste0(" test.validateParameters_all_valid() ",
                      "- All valid parameters did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

#########################################################
## validateRDSFilesParameters() function
#########################################################

## Test the result when RDSFiles is NA
test.validateRDSFilesParameters_RDSFiles_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRDSFilesParameters(
        RDSFiles = NA), error=conditionMessage)
    exp <- "RDSFiles must be a list of valid RDS files"
    message <- paste0(" test.validateRDSFilesParameters_RDSFiles_NA() ",
                        "- NA for RDSFiles did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when RDSFiles is empty array
test.validateRDSFilesParameters_RDSFiles_empty_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRDSFilesParameters(
        RDSFiles = c()), error=conditionMessage)
    exp <- "RDSFiles must be a list of valid RDS files"
    message <- paste0(" test.validateRDSFilesParameters_RDSFiles_empty_array() ",
                        "- Empty array for RDSFiles did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## postMerge() function
#########################################################


## Test the result of postMerge() function
test.postMerge_good_01 <- function() {
    obs <- RJMCMCNucleosomes:::postMerge(startPosForwardReads = reads_demo$readsForward,
                                startPosReverseReads = reads_demo$readsReverse,
                                resultRJMCMC = RJMCMC_result,
                                extendingSize = 10, chrLength = 80000)
    exp <- c(72434.766272478853, 72544.048047704578, 73146.590899701128)
    message <- paste0(" test.postMerge_good_01() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}


## Test the result of postMerge() function
test.postMerge_good_02 <- function() {
    obs <- RJMCMCNucleosomes:::postMerge(startPosForwardReads = reads_demo$readsForward,
                              startPosReverseReads = reads_demo$readsReverse,
                              resultRJMCMC = RJMCMC_result,
                              extendingSize = 100, chrLength = 80000)
    exp <- c(72452.452375092398, 73146.590899701128)
    message <- paste0(" test.postMerge_good_02() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}
