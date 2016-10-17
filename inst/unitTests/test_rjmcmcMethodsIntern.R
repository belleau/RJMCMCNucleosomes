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
                pattern = "newSeg_2.rds",
                full.names = TRUE)

data_002 <- readRDS(file_002)

data(RJMCMC_result)
data(reads_demo)
data(syntheticNucleosomeReads)

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
## validatePlotNucleosomesParameters() function
#########################################################

## Test the result when nucleosomePositions is NA
test.validatePlotNucleosomesParameters_nucleosomePositions_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = NA, reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions can only contain numerical values"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_NA() ",
                      "- NA for nucleosomePositions did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomesParameters is empty vector
test.validatePlotNucleosomesParameters_nucleosomePositions_empty_vector <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(), reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions must be a \'list\' or a \'vector\' of numeric values"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_empty_vector() ",
                      "- Empty vector for nucleosomePositions did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomesParameters is not numerical
test.validatePlotNucleosomesParameters_nucleosomePositions_not_numerical <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c("hi", "test"), reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions can only contain numerical values"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_not_numerical() ",
                      "- Not numeric for nucleosomePositions did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomesParameters is not a list of numerical
test.validatePlotNucleosomesParameters_nucleosomePositions_list_not_numerical <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(a=c(1, 3), b=c("hi", "test")), reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions can only contain numerical values"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_list_not_numerical() ",
                      "- Not list of numeric for nucleosomePositions did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is null
test.validatePlotNucleosomesParameters_reads_null <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(1001, 1003), reads = NULL, xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "reads must be an object of class \'IRanges\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_reads_null() ",
                      "- Not list of numeric for reads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is not IRanges
test.validatePlotNucleosomesParameters_reads_not_IRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(1001, 1003), reads = c(100,200), xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "reads must be an object of class \'IRanges\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_reads_not_IRanges() ",
                      "- Not IRanges for reads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when xlab is not a character string
test.validatePlotNucleosomesParameters_xlab_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(1001, 1003), reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = 33,
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "xlab must be a character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_xlab_not_string() ",
                      "- Not character string for xlab did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when ylab is not a character string
test.validatePlotNucleosomesParameters_ylab_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(1001, 1003), reads = IRanges(start=c(950, 969),
        end=c(1020, 1022)), xlab = "x",
        ylab = c(44,33), names=c("test")),
        error=conditionMessage)
    exp <- "ylab must be a character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_xylab_not_string() ",
                      "- Not character string for xlab did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names is not a character string
test.validatePlotNucleosomesParameters_names_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(1001, 1003), reads = IRanges(start=c(950, 969),
                                                             end=c(1020, 1022)), xlab = "x",
        ylab = "y", names=c(33)),
        error=conditionMessage)
    exp <- "names must be a vector of one character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_string() ",
                      "- Not character string for names did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names length does not fit nucleosomePositions entries
test.validatePlotNucleosomesParameters_names_not_good_length_01 <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = list(a=c(1001), b=c(1003)),
        reads = IRanges(start=c(950, 969), end=c(1020, 1022)), xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "names must be a vector containing the same number of character string as the number of entries in nucleosomesPositions list"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_good_length_01
                      () ",
                      "- Not good length for names did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names length does not fit nucleosomePositions entries
test.validatePlotNucleosomesParameters_names_not_good_length_02 <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = list(a=c(1001)),
        reads = IRanges(start=c(950, 969), end=c(1020, 1022)), xlab = "x",
        ylab = "y", names=c("test", "test02")),
        error=conditionMessage)
    exp <- "names must be a vector containing the same number of character string as the number of entries in nucleosomesPositions list"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_good_length_02() ",
                      "- Not good length for names did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test tha all valid parameters return zero
test.validatePlotNucleosomesParameters_all_good  <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = list(a=c(1001)),
        reads = IRanges(start=c(950, 969), end=c(1020, 1022)), xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- 0
    message <- paste0(" test.validatePlotNucleosomesParameters_all_good() ",
                      "- All good parameters did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateSegmentationParameters() function
#########################################################

## Test the result when dataIP is NA
test.validateSegmentationParameters_dataIP_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = NA, zeta = 147, delta = 12, maxLength = 20000),
        error=conditionMessage)
    exp <- "dataIP must be \'GRanges\' object."
    message <- paste0(" test.validateSegmentationParameters_dataIP_NA() ",
                      "- NA for dataIP did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when dataIP is not GRanges
test.validateSegmentationParameters_dataIP_not_GRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = c(1, 3, 2), zeta = 147, delta = 12, maxLength = 20000),
        error=conditionMessage)
    exp <- "dataIP must be \'GRanges\' object."
    message <- paste0(" test.validateSegmentationParameters_dataIP_not_GRanges() ",
                      "- Not GRanges for dataIP did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zeta a vector of numeric
test.validateSegmentationParameters_zeta_vector <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
        ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
        end = syntheticNucleosomeReads$dataIP$end),
        strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = c(147, 12), delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_vector() ",
                      "- Vector for zeta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zeta is zero
test.validateSegmentationParameters_zeta_zero <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 0, delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_zero() ",
                      "- Zero value for zeta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zeta is negative
test.validateSegmentationParameters_zeta_negative <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = -1, delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_negative() ",
                      "- Negative value for zeta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta a vector of numeric
test.validateSegmentationParameters_delta_vector <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = c(11, 21),
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_vector() ",
                      "- Vector for delta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta is zero
test.validateSegmentationParameters_delta_zero <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = 0,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_zero() ",
                      "- Zero value for delta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta is negative
test.validateSegmentationParameters_delta_negative <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = -1,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_negative() ",
                      "- Negative value for delta did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength a vector of numeric
test.validateSegmentationParameters_maxLength_vector <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = 12,
        maxLength = c(10, 20)),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_vector() ",
                      "- Vector for maxLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength is zero
test.validateSegmentationParameters_maxLength_zero <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = 12,
        maxLength = 0),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_zero() ",
                      "- Zero value for maxLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength is negative
test.validateSegmentationParameters_maxLength_negative <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = 12,
        maxLength = -1),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_negative() ",
                      "- Negative value for maxLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test when all parameters are valids
test.validateSegmentationParameters_all_valid <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        dataIP = sampleGRanges, zeta = 147, delta = 12,
        maxLength = 2000),
        error=conditionMessage)
    exp <- 0
    message <- paste0(" test.validateSegmentationParameters_all_valid() ",
                      "- All valid parameters did not  ",
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
