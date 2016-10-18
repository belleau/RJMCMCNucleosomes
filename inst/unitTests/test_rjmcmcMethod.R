###################################################
# Created by Astrid Deschenes
# 2015-06-30
###################################################

###################################################
## Test the rjmcmcMethod.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMCNucleosomes" )
}

### }}}

data(reads_demo)
data(reads_demo_02)
data(RJMCMC_result)
data(syntheticNucleosomeReads)

DIRECTORY <- system.file("extdata", package = "RJMCMCNucleosomes")

file_001 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "newSeg_1.rds",
                        full.names = TRUE)

file_002 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "newSeg_2.rds",
                        full.names = TRUE)

file_003 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "newSeg_3.rds",
                        full.names = TRUE)


###########################################################
## RJMCMC() function
###########################################################

test.rjmcmc_one_read_forward_and_one_read_reverse <- function() {

    obs <- rjmcmc(startPosForwardReads = c(1),
                  startPosReverseReads = c(2),
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10,
                  vSeed = 2211)

    exp.k           <- 1
    exp.k_max       <- 1
    exp.mu          <- c(1.017962247133255)

    message     <- paste0(" test.rjmcmc_one_read_forward_and_one_read_reverse() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
}


test.rjmcmc_good_result_01 <- function() {

    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                        startPosReverseReads = reads_demo$readsReverse,
                        nbrIterations = 100, lambda = 2, kMax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5,
                        vSeed = 1001)

    exp.k           <- 4
    exp.k_max       <- 5
    exp.mu          <- c(72393.0264332128, 72457.6656495946, 72555.4429585791, 72990.7142984954)

    message     <- paste0(" rjmcmc_good_result_01() ",
                       "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
}

test.rjmcmc_good_result_02 <- function() {

    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                    startPosReverseReads = reads_demo$readsReverse,
                    nbrIterations = 200, lambda = 3, kMax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5,
                    vSeed = 201)

    exp.k           <- 4
    exp.k_max       <- 4
    exp.mu          <- c(72325.7014073403, 72357.7522020026, 72678.9998565154, 72875.6569724952)

    message     <- paste0(" rjmcmc_good_result_02() ",
                      "- RJMCMC did not generated expected values")


    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
}

test.rjmcmc_good_result_03 <- function() {

    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                  startPosReverseReads = reads_demo$readsReverse,
                  nbrIterations = 110, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 335,
                  vSeed = 2011)

    exp.k           <- 4
    exp.k_max       <- 4
    exp.mu          <- c(72564.6044994216, 72863.5808919129, 72972.7389043033, 73207.8950726675)

    message     <- paste0(" rjmcmc_good_result_03() ",
                           "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
}

test.rjmcmc_good_result_04 <- function() {

    obs <- rjmcmc(startPosForwardReads = reads_demo_02$readsForward,
                  startPosReverseReads = reads_demo_02$readsReverse,
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10,
                  vSeed = 2211)

    exp.k           <- 5
    exp.k_max       <- 6
    exp.mu          <- c(18839.4528131919, 19201.6625115471, 19343.2693614808, 19415.4983810908, 19515.0250764107)

    message     <- paste0(" test.rjmcmc_good_result_04() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
}




###########################################################
## mergeAllRDSFilesFromDirectory() function
###########################################################


test.mergeAllRDSFilesFromDirectory_notExisting <- function() {
    dir_01 <- "/toto1/toto2/toto3/toto4/toto5/"
    dir_02 <- "/toto5/toto4/toto3/toto2/toto1/"

    dir <- NULL
    if (!file.exists(dir_01)) {
        dir <- dir_01
    } else {
        if (!file.exists(dir_02)) {
            dir <- dir_02
        }
    }

    if (!is.null(dir)) {
        obs <- tryCatch(mergeAllRDSFilesFromDirectory(dir),
                    error=conditionMessage)
        exp <- paste0("The directory \'", dir,
                  "\' does not exist.")
        message <- paste0(" test.mergeResultFilesInDirectory_notExisting() ",
                      "- A not existing directory did not generated ",
                      "expected message.")
        checkEquals(obs, exp, msg = message)
    }
}

test.mergeAllRDSFilesFromDirectory_good <- function() {

    obs <- mergeAllRDSFilesFromDirectory(DIRECTORY)
    exp <- list()
    exp$k <- 16
    exp$mu <- c(10092.474777629515302, 10242.340786347993344,
                10410.315021756090573, 10546.628207912892321,
                11134.263941022001745, 11244.139414670189581,
                11380.302471258171863, 11412.657313642583176,
                11578.516646129490255, 11868.408425173787691,
                12058.054137626086231, 12235.422415610730241,
                12276.903444548192056, 12412.604063330700228,
                12473.443670263422973, 12585.175197400152683)

    class(exp) <- "rjmcmcNucleosomesMerge"

    message <- paste0(" test.mergeAllRDSFilesFromDirectory_good() ",
                      "- The mergeAllRDSFilesFromDirectory() did not generated ",
                      "expected output.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## mergeRDSFiles() function
###########################################################

test.mergeRDSFiles_notExisting <- function() {
    file_01 <- "/toto1/toto2/toto3/toto4/toto5/improbable_file_01335320111.RDS"
    file_02 <- "/toto5/toto4/toto3/toto2/toto1/improbable_file_01335320111.RDS"

    fileName <- array(dim = c(0))
    if (!file.exists(file_01)) {
        fileName <- c(file_01)
    }
    if (!file.exists(file_02)) {
        fileName <- c(fileName, file_02)
    }

    if (length(fileName) > 0) {
        obs <- tryCatch(mergeRDSFiles(fileName),
                        error=conditionMessage)
        exp <- paste0("The file \'", fileName[1],
                      "\' does not exist.")
        message <- paste0(" test.mergeRDSFiles_notExisting() ",
                          "- A not existing file did not generated ",
                          "expected message.")
        checkEquals(obs, exp, msg = message)
    }
}

test.mergeRDSFiles_good <- function() {

    files <- c(file_001, file_002, file_003)

    obs <- mergeRDSFiles(files)
    exp <- list()
    exp$k <- 16
    exp$mu <- c(10092.474777629515302, 10242.340786347993344,
                10410.315021756090573, 10546.628207912892321,
                11134.263941022001745, 11244.139414670189581,
                11380.302471258171863, 11412.657313642583176,
                11578.516646129490255, 11868.408425173787691,
                12058.054137626086231, 12235.422415610730241,
                12276.903444548192056, 12412.604063330700228,
                12473.443670263422973, 12585.175197400152683)

    class(exp) <- "rjmcmcNucleosomesMerge"

    message <- paste0(" test.mergeRDSFiles_good() ",
                      "- The mergeRDSFiles() did not generated ",
                      "expected output.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## postTreatment() function
###########################################################

test.postTreatment_good_01 <- function() {

    obs <- postTreatment(startPosForwardReads = reads_demo$readsForward,
                              startPosReverseReads = reads_demo$readsReverse,
                              resultRJMCMC = RJMCMC_result,
                              extendingSize = 20,
                              chrLength = 80000)

    exp <- c(72434.76627247885335236788,
                72544.04804770457849372178,
                73146.59089970112836454064)

    message <- paste0(" test.postTreatment_good_01() ",
                      "- posTreatment() did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}

test.postTreatment_good_02 <- function() {

    obs <- postTreatment(startPosForwardReads = reads_demo$readsForward,
                        startPosReverseReads = reads_demo$readsReverse,
                        resultRJMCMC = RJMCMC_result,
                        extendingSize = 200,
                        chrLength = 80000)

    exp <- c(72533.80877122777746990323,
             73146.59089970112836454064)

    message <- paste0(" test.postTreatment_good_02() ",
                      "- posTreatment() did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## segmentation() function
###########################################################

test.segmentation_good_01 <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                    ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                    end = syntheticNucleosomeReads$dataIP$end),
                    strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- segmentation(sampleGRanges, zeta =  147, delta = 20, maxLength = 20000)

    message <- paste0(" test.segmentation_good_01() ",
                      "- segmentation() did not generated expected result.")

    exp.len = 3
    exp.01.len = 10504
    exp.02.len = 11818
    exp.03.len = 9686

    checkTrue(is.list(obs), ms = message)
    checkEquals(length(obs), exp.len, ms = message)
    checkEquals(length(obs[[1]]), exp.01.len, ms = message)
    checkTrue(is(obs[[1]],"GRanges"), ms = message)
    checkEquals(length(obs[[2]]), exp.02.len, ms = message)
    checkTrue(is(obs[[2]],"GRanges"), ms = message)
    checkEquals(length(obs[[3]]), exp.03.len, ms = message)
    checkTrue(is(obs[[3]],"GRanges"), ms = message)
}

test.segmentation_good_02  <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                             ranges = IRanges(start = syntheticNucleosomeReads$dataIP$start,
                                              end = syntheticNucleosomeReads$dataIP$end),
                             strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- segmentation(sampleGRanges, zeta =  142, delta = 40, maxLength = 15000)

    message <- paste0(" test.segmentation_good_02() ",
                      "- segmentation() did not generated expected result.")

    exp.len = 4
    exp.01.len = 7972
    exp.02.len = 8496
    exp.03.len = 9362
    exp.04.len = 6390

    checkTrue(is.list(obs), ms = message)
    checkEquals(length(obs), exp.len, ms = message)
    checkEquals(length(obs[[1]]), exp.01.len, ms = message)
    checkTrue(is(obs[[1]],"GRanges"), ms = message)
    checkEquals(length(obs[[2]]), exp.02.len, ms = message)
    checkTrue(is(obs[[2]],"GRanges"), ms = message)
    checkEquals(length(obs[[3]]), exp.03.len, ms = message)
    checkTrue(is(obs[[3]],"GRanges"), ms = message)
}

