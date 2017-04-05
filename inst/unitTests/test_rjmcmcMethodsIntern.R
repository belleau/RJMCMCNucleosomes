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
                pattern = "RJMCMC_seg_02.RDS",
                full.names = TRUE)

data_002 <- readRDS(file_002)

data(RJMCMC_result)
data(reads_demo_01)
data(reads_demo_02)
data(syntheticNucleosomeReads)

#########################################################
## validatePrepMergeParameters() function
#########################################################

## Test the result when reads is NA
test.validatePrepMergeParameters_reads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
                        reads = NA, seqName = NULL,
                        resultRJMCMC = data_002, extendingSize = 11,
                        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("reads must be a GRanges")
    message <- paste0(" test.validatePrepMergeParameters_reads_NA() ",
                      "- NA for reads did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is empty GRanges
test.validatePrepMergeParameters_reads_empty <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = GRanges(), seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("reads must be a non-empty GRanges")
    message <- paste0(" test.validatePrepMergeParameters_reads_empty() ",
                      "- reads as empty GRanges did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is not GRanges
test.validatePrepMergeParameters_reads_not_GRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = c("A", "B"), seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("reads must be a GRanges")
    message <- paste0(" test.validatePrepMergeParameters_reads_not_GRanges() ",
                      "- not GRanges reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName is NULL with GRanges with multiples chromosomes
test.validatePrepMergeParameters_seqName_NULL_GRanges_complex <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("seqName must be the name of one of the chromosomes ",
                  "present in the GRanges")
    message <- paste0(" test.validatePrepMergeParameters_seqName_NULL_GRanges_complex() ",
                      "- Complex GRanges for reads ",
                      "and NULL for seqName did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not a string
test.validatePrepMergeParameters_seqName_not_string <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads, seqName = 333,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validatePrepMergeParameters_seqName_not_string() ",
                      "- seqName not a string ",
                      "did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not in GRanges
test.validatePrepMergeParameters_seqName_not_in_GRanges <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads, seqName = "chr3",
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validatePrepMergeParameters_seqName_not_in_GRanges() ",
                      "- seqName not in GRanges ",
                      "did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is NA
test.validatePrepMergeParameters_resultRJMCMC_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = NA, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_not_number() ",
                      "- NA resultRJMCMC did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is a number
test.validatePrepMergeParameters_resultRJMCMC_number <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = 33, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_resultRJMCMC_number() ",
                      "- number resultRJMCMC did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is a string
test.validatePrepMergeParameters_nbBase_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = "ALLO",
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_number() ",
                      "- string nbBase did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is an array
test.validatePrepMergeParameters_nbBase_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = c(10, 11),
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_string() ",
                      "- array nbBase did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is a string
test.validatePrepMergeParameters_chrLength_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = "5000"), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- string chrLength did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is an array
test.validatePrepMergeParameters_chrLength_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = c(100, 200)), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- array chrLength did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validatePrepMergeParameters_all_valid <- function() {
    obs <- RJMCMCNucleosomes:::validatePrepMergeParameters(
        reads = reads_demo_02, seqName = NULL,
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = 200000)
    exp <- 0
    message <- paste0(" test.validatePrepMergeParameters_all_valid() ",
                      "- All valid parameters did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateRJMCMCParameters() function
#########################################################

## Test the result when nbrIterations is NA
test.validateRJMCMCParameters_nbrIterations_NA <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                    ranges = IRanges(101:110, end = 111:120,
                                        names = head(letters, 10)),
                    strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                        c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL,
        nbrIterations = NA,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_NA() ",
                      "- NA for nbrIterations did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is zero
test.validateRJMCMCParameters_nbrIterations_zero <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL,
        nbrIterations = 0,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Zero for nbrIterations did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is negative
test.validateRJMCMCParameters_nbrIterations_negative <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = -1,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Negative value for nbrIterations did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is NA
test.validateRJMCMCParameters_kMax_NA <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = NA, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_NA() ",
                      "- NA value for kMax did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is zero
test.validateRJMCMCParameters_kMax_zero <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters (
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 0, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_zero() ",
                      "- Zero value for kMax did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is negative
test.validateRJMCMCParameters_kMax_negative <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters (
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = -1, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_negative() ",
                      "- Negative value for kMax did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is NA
test.validateRJMCMCParameters_minReads_NA <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = NA, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for minReads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is zero
test.validateRJMCMCParameters_minReads_zero <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 0, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for minReads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is negative
test.validateRJMCMCParameters_minReads_negative <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 3, lambda = 1, minReads = -1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for minReads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is NA
test.validateRJMCMCParameters_lambda_NA <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = NA, minReads = 2, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for lambda did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is zero
test.validateRJMCMCParameters_lambda_zero <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL,
        nbrIterations = 2,
        kMax = 10, lambda = 0, minReads = 3, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for lambda did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is negative
test.validateRJMCMCParameters_lambda_negative <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = -1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for lambda did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is NA
test.validateRJMCMCParameters_reads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = NA, nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("reads must be a GRanges")
    message <- paste0(" test.validateRJMCMCParameters_reads_NA() ",
                      "- NA value for reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is empty
test.validateRJMCMCParameters_reads_empty <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = GRanges(), seqName = NULL,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292, vSeed = -1,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- 0
    message <- paste0(" test.validateRJMCMCParameters_reads_empty() ",
                        "- Empty GRanges for reads did not ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName is NULL with GRanges with multiples chromosomes
test.validateRJMCMCParameters_seqName_NULL_GRanges_complex <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL, nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292, vSeed = -1,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("seqName must be the name of one of the chromosomes ",
                  "present in the GRanges")
    message <- paste0(" test.validateRJMCMCParameters_seqName_NULL_GRanges_complex() ",
                      "- Complex GRanges for forwardandReverseReads ",
                      "and NULL for seqName did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not in GRanges
test.validateRJMCMCParameters_seqName_not_string <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = 333,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292, vSeed = -1,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validateRJMCMCParameters_seqName_not_string() ",
                      "- seqName not a string ",
                      "did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not in GRanges
test.validateRJMCMCParameters_seqName_not_in_GRanges <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr3",
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292, vSeed = -1,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validateRJMCMCParameters_seqName_not_in_GRanges() ",
                      "- seqName not in GRanges ",
                      "did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}


## Test the result when adaptIterationsToReads is string
test.validateRJMCMCParameters_adaptIterationsToReads_string <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = "allo"), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_string() ",
                        "- String for adaptIterationsToReads did not ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when adaptIterationsToReads is number
test.validateRJMCMCParameters_adaptIterationsToReads_number <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = NULL,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = 33), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_number() ",
                        "- Number value for adaptIterationsToReads did not ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when vSeed is not a number
test.validateRJMCMCParameters_vSeed_number <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE, vSeed = "Allo"), error=conditionMessage)
    exp <- "vSeed must be a numeric value."
    message <- paste0(" test.validateRJMCMCParameters_vSeed_number() ",
                      "- String value for vSeed did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validateRJMCMCParameters_all_valid <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- RJMCMCNucleosomes:::validateRJMCMCParameters(
        reads = reads, seqName = "chr1",
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = TRUE, vSeed = 1002)
    exp <- 0
    message <- paste0(" test.validateParameters_all_valid() ",
                      "- All valid parameters did not ",
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
                        "- NA for RDSFiles did not ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when RDSFiles is empty array
test.validateRDSFilesParameters_RDSFiles_empty_array <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateRDSFilesParameters(
        RDSFiles = c()), error=conditionMessage)
    exp <- "RDSFiles must be a list of valid RDS files"
    message <- paste0(" test.validateRDSFilesParameters_RDSFiles_empty_array() ",
                        "- Empty array for RDSFiles did not ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validatePlotNucleosomesParameters() function
#########################################################

## Test the result when nucleosomePositions is NA
test.validatePlotNucleosomesParameters_nucleosomePositions_NA <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = NA, reads  = reads,
        seqName = "chr1",
        xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions must be a \'GRanges\' or a \'GRangesList\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_NA() ",
                      "- NA for nucleosomePositions did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomesParameters is empty vector
test.validatePlotNucleosomesParameters_nucleosomePositions_empty_vector <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(), reads = reads,
        seqName = "chr1",
        xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions must be a \'GRanges\' or a \'GRangesList\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_empty_vector() ",
                      "- Empty vector for nucleosomePositions did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomePositions is not GRanges
test.validatePlotNucleosomesParameters_nucleosomePositions_string <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c("hi", "test"), reads = reads,
        seqName = "chr1", xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "nucleosomePositions must be a \'GRanges\' or a \'GRangesList\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_string() ",
                      "- String for nucleosomePositions did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nucleosomesParameters is not a GRanges
test.validatePlotNucleosomesParameters_nucleosomePositions_not_GRanges <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = c(a=c(1, 3), b=c("hi", "test")),
        reads = reads, seqName = "chr1", xlab = "x", ylab = "y",
        names=c("test")), error=conditionMessage)
    exp <- "nucleosomePositions must be a \'GRanges\' or a \'GRangesList\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_nucleosomePositions_not_GRanges() ",
                      "- Not GRanges for nucleosomePositions did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is null
test.validatePlotNucleosomesParameters_reads_null <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu, reads = NULL,
        seqName = NULL, xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "reads must be an object of class \'GRanges\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_reads_null() ",
                      "- NULL reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is empty
test.validatePlotNucleosomesParameters_reads_empty_GRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu, reads = GRanges(),
        seqName = NULL, xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "reads must be a non-empty GRanges"
    message <- paste0(" test.validatePlotNucleosomesParameters_reads_empty_GRanges() ",
                      "- Empty reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName is null and GRanges complexe
test.validatePlotNucleosomesParameters_seqName_null_with_complexe_GRanges <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(8,2)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu, reads = reads,
        seqName = NULL, xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- paste0("seqName must be the name of one of the chromosomes ",
        "present in the GRanges when more than one chromosome is present")
    message <- paste0(" test.validatePlotNucleosomesParameters_seqName_null_with_complexe_GRanges() ",
                      "- Null seqName with complexe GRanges did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


## Test the result when reads is not GRanges
test.validatePlotNucleosomesParameters_reads_not_GRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu, reads = c(1,2),
        seqName = NULL, xlab = "x",
        ylab = "y", names=c("test")), error=conditionMessage)
    exp <- "reads must be an object of class \'GRanges\'"
    message <- paste0(" test.validatePlotNucleosomesParameters_forwardandReverseReads_not_GRanges() ",
                      "- Not GRanges for reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not string
test.validatePlotNucleosomesParameters_seqName_not_string <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions =  data_002$mu,
        reads = reads, seqName = 44, xlab = "x", ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validatePlotNucleosomesParameters_seqName_not_string() ",
                      "- SeqName not a string did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when seqName not in GRanges
test.validatePlotNucleosomesParameters_seqName_not_in_GRanges <- function() {
    reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
                     ranges = IRanges(101:110, end = 111:120,
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")),
                                  c(1, 2, 2, 3, 2)))
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions =  data_002$mu,
        reads = reads, seqName = "chr1_test", xlab = "x", ylab = "y",
        names=c("test")), error=conditionMessage)
    exp <- paste0("seqName must be a character string corresponding to the ",
                  "name of one of the chromosomes present in the GRanges")
    message <- paste0(" test.validatePlotNucleosomesParameters_seqName_not_in_GRanges() ",
                      "- SeqName not in GRanges did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when xlab is not a character string
test.validatePlotNucleosomesParameters_xlab_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu, reads = reads_demo_01,
        seqName = "chr_SYNTHETIC", xlab = 33,
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "xlab must be a character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_xlab_not_string() ",
                      "- Not character string for xlab did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when ylab is not a character string
test.validatePlotNucleosomesParameters_ylab_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu,
        reads = reads_demo_01, seqName = "chr_SYNTHETIC", xlab = "x",
        ylab = c(44,33), names=c("test")),
        error=conditionMessage)
    exp <- "ylab must be a character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_xylab_not_string() ",
                      "- Not character string for xlab did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names is not a character string
test.validatePlotNucleosomesParameters_names_not_string <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = data_002$mu,
        reads = reads_demo_01, seqName = "chr_SYNTHETIC", xlab = "x",
        ylab = "y", names=c(33)),
        error=conditionMessage)
    exp <- "names must be a vector of one character string"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_string() ",
                      "- Not character string for names did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names length does not fit nucleosomePositions entries
test.validatePlotNucleosomesParameters_names_not_good_length_01 <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = GRangesList(a=data_002$mu, b=data_002$mu),
        reads = reads_demo_01, seqName = NULL, xlab = "x",
        ylab = "y", names=c("test")),
        error=conditionMessage)
    exp <- "names must be a vector containing the same number of character string as the number of entries in nucleosomesPositions list"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_good_length_01
                      () ",
                      "- Not good length for names did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when names length does not fit nucleosomePositions entries
test.validatePlotNucleosomesParameters_names_not_good_length_02 <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = GRangesList(a=data_002$mu),
        reads = reads_demo_01, seqName = NULL, xlab = "x",
        ylab = "y", names=c("test", "test02")), error=conditionMessage)
    exp <- "names must be a vector containing the same number of character string as the number of entries in nucleosomesPositions list"
    message <- paste0(" test.validatePlotNucleosomesParameters_names_not_good_length_02() ",
                      "- Not good length for names did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test tha all valid parameters return zero
test.validatePlotNucleosomesParameters_all_good  <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validatePlotNucleosomesParameters(
        nucleosomePositions = GRangesList(a=data_002$mu), seqName = "chr_SYNTHETIC",
        reads = reads_demo_01, xlab = "x",
        ylab = "y", names=c("test")), error=conditionMessage)
    exp <- 0
    message <- paste0(" test.validatePlotNucleosomesParameters_all_good() ",
                      "- All good parameters did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateSegmentationParameters() function
#########################################################

## Test the result when reads is NA
test.validateSegmentationParameters_reads_NA <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = NA, zeta = 147, delta = 12, maxLength = 20000),
        error=conditionMessage)
    exp <- "reads must be \'GRanges\' object."
    message <- paste0(" test.validateSegmentationParameters_reads_NA() ",
                      "- NA for reads did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when reads is not GRanges
test.validateSegmentationParameters_reads_not_GRanges <- function() {
    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = c(1, 3, 2), zeta = 147, delta = 12, maxLength = 20000),
        error=conditionMessage)
    exp <- "reads must be \'GRanges\' object."
    message <- paste0(" test.validateSegmentationParameters_reads_not_GRanges() ",
                      "- Not GRanges for reads did not ",
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
        reads = sampleGRanges, zeta = c(147, 12), delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_vector() ",
                      "- Vector for zeta did not ",
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
        reads = sampleGRanges, zeta = 0, delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_zero() ",
                      "- Zero value for zeta did not ",
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
        reads = sampleGRanges, zeta = -1, delta = 12,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "zeta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_zeta_negative() ",
                      "- Negative value for zeta did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta a vector of numeric
test.validateSegmentationParameters_delta_vector <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = c(11, 21),
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_vector() ",
                      "- Vector for delta did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta is zero
test.validateSegmentationParameters_delta_zero <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = 0,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_zero() ",
                      "- Zero value for delta did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when delta is negative
test.validateSegmentationParameters_delta_negative <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = -1,
        maxLength = 20000),
        error=conditionMessage)
    exp <- "delta must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_delta_negative() ",
                      "- Negative value for delta did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength a vector of numeric
test.validateSegmentationParameters_maxLength_vector <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                            ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                            strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = 12,
        maxLength = c(10, 20)),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_vector() ",
                      "- Vector for maxLength did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength is zero
test.validateSegmentationParameters_maxLength_zero <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = 12,
        maxLength = 0),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_zero() ",
                      "- Zero value for maxLength did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when maxLength is negative
test.validateSegmentationParameters_maxLength_negative <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = 12,
        maxLength = -1),
        error=conditionMessage)
    exp <- "maxLength must be a positive integer or numeric"
    message <- paste0(" test.validateSegmentationParameters_maxLength_negative() ",
                      "- Negative value for maxLength did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test when all parameters are valids
test.validateSegmentationParameters_all_valid <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- tryCatch(RJMCMCNucleosomes:::validateSegmentationParameters(
        reads = sampleGRanges, zeta = 147, delta = 12,
        maxLength = 2000),
        error=conditionMessage)
    exp <- 0
    message <- paste0(" test.validateSegmentationParameters_all_valid() ",
                      "- All valid parameters did not ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

#########################################################
## postMerge() function
#########################################################


## Test the result of postMerge() function
test.postMerge_good_01 <- function() {

    obs <- RJMCMCNucleosomes:::postMerge(reads = reads_demo_02,
                                resultRJMCMC = RJMCMC_result,
                                extendingSize = 10, chrLength = 80000)
    listMu <- c(10072, 10241,
                10574, 10665,
                10744)
    exp <- GRanges(seqnames = rep("chr_SYNTHETIC",length(listMu)),
                      ranges=IRanges(start=listMu,end=listMu),
                      strand=rep("*", length(listMu)))

    message <- paste0(" test.postMerge_good_01() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result of postMerge() function when mu is set to NA
test.postMerge_mu_NA_returned <- function() {

    testRJMCMC <- RJMCMC_result
    testRJMCMC$mu <- NA

    ## Reads are not located where the nucleosomes are
    obs <- RJMCMCNucleosomes:::postMerge(reads = reads_demo_01,
                                        resultRJMCMC = testRJMCMC,
                                        extendingSize = 100, chrLength = 80000)
    exp <- NULL
    message <- paste0(" test.postMerge_mu_NA_returned() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result of postMerge() function when empty mu is passed
test.postMerge_mu_empty_returned <- function() {

    testRJMCMC <- RJMCMC_result
    testRJMCMC$mu <- c()

    ## Reads are not located where the nucleosomes are
    obs <- RJMCMCNucleosomes:::postMerge(reads = reads_demo_01,
                                        resultRJMCMC = testRJMCMC,
                                        extendingSize = 100, chrLength = 80000)
    exp <- NULL
    message <- paste0(" test.postMerge_mu_NA_returned() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}

