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
    library( "GenomicRanges")
}

### }}}

data(reads_demo_01)
data(reads_demo_02)
data(RJMCMC_result)
data(syntheticNucleosomeReads)

DIRECTORY <- system.file("extdata", package = "RJMCMCNucleosomes")

file_001 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "RJMCMC_seg_01.RDS",
                        full.names = TRUE)

file_002 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "RJMCMC_seg_02.RDS",
                        full.names = TRUE)

file_003 <- dir(system.file("extdata", package = "RJMCMCNucleosomes"),
                        pattern = "RJMCMC_seg_03.RDS",
                        full.names = TRUE)


###########################################################
## RJMCMC() function
###########################################################

test.rjmcmc_one_read_forward_and_one_read_reverse <- function() {

    reads <- GRanges(seqnames = Rle(c("chr1"), c(2)),
                            ranges = IRanges(101:102, end = 111:112,
                                        names = head(letters, 2)),
                            strand = Rle(strand(c("-", "+")), c(1, 1)))
    obs <- rjmcmc(reads = reads, nbrIterations = 210, lambda = 3, kMax = 30,
                    minInterval = 100, maxInterval = 200, minReads = 10,
                    vSeed = 2211)

    exp.k           <- 1
    exp.k_max       <- 1
    exp.mu          <- c(102)
    exp.strand      <- c('*')
    exp.seqnames    <- c("chr1")

    message     <- paste0(" test.rjmcmc_one_read_forward_and_one_read_reverse() ",
                            "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(strand(obs$mu), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
}


test.rjmcmc_good_result_01 <- function() {

    obs <- rjmcmc(reads = reads_demo_01,
                        nbrIterations = 100, lambda = 2, kMax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5,
                        vSeed = 1001)

    exp.k           <- 2
    exp.k_max       <- 2
    exp.mu          <- as.integer(c(10080, 10146))
    exp.strand      <- c('*', '*')
    exp.seqnames    <- rep("chr_SYNTHETIC", 2)

    message     <- paste0(" rjmcmc_good_result_01() ",
                        "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(strand(obs$mu), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
}

test.rjmcmc_good_result_02 <- function() {

    obs <- rjmcmc(reads  = reads_demo_01,
                    nbrIterations = 200, lambda = 3, kMax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5,
                    vSeed = 201)

    exp.k           <- 1
    exp.k_max       <- 3
    exp.mu          <- c(10056)
    exp.strand      <- c('*')
    exp.seqnames    <- c("chr_SYNTHETIC")

    message     <- paste0(" rjmcmc_good_result_02() ",
                        "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(strand(obs$mu), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
}

test.rjmcmc_good_result_03 <- function() {

    obs <- rjmcmc(reads = reads_demo_01,
                  nbrIterations = 110, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 335,
                  vSeed = 2011)

    exp.k           <- 2
    exp.k_max       <- 4
    exp.mu          <- c(10058, 10453)
    exp.strand      <- c('*', '*')
    exp.seqnames    <- rep("chr_SYNTHETIC", 2)

    message     <- paste0(" rjmcmc_good_result_03() ",
                           "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(strand(obs$mu), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
}

test.rjmcmc_good_result_04 <- function() {

    obs <- rjmcmc(reads = reads_demo_01,
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10,
                  vSeed = 2211)

    exp.k           <- 1
    exp.k_max       <- 2
    exp.mu          <- c(10077)
    exp.strand      <- c('*')
    exp.seqnames    <- c("chr_SYNTHETIC")

    message     <- paste0(" test.rjmcmc_good_result_04() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(strand(obs$mu), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
}

test.rjmcmc_good_result_with_GRanges_with_multiple_names  <- function() {

    reads <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(5, 5)),
                     ranges = IRanges(start=c(11:15, 1106:1110), end = c(21:25, 1116:1120),
                                      names = head(letters, 10)),
                     strand = Rle(strand(c("-", "+", "-", "+", "-")), c(2, 2, 2, 2, 2)))
    obs <- rjmcmc(reads = reads, seqName = "chr2",
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10,
                  vSeed = 2211)

    exp.k           <- 1
    exp.k_max       <- 1
    exp.mu          <- c(1107)
    exp.strand      <- c('*')
    exp.seqnames    <- c("chr2")

    message     <- paste0(" test.rjmcmc_good_result_with_GRanges_with_multiple_names() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$k_max, exp.k_max, msg = message)
    checkTrue(class(obs$mu) == "GRanges", msg = message)
    checkEqualsNumeric(start(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(end(obs$mu), exp.mu, msg = message)
    checkEqualsNumeric(as.vector(strand(obs$mu)), exp.strand, msg = message)
    checkEqualsNumeric(seqnames(obs$mu), exp.seqnames, msg = message)
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
    exp$k <- 11
    listMu <- c(10077, 10236,
                10406, 10571,
                10744, 10842,
                10846, 10896,
                10906, 11410,
                11580)
    exp$mu <- GRanges(seqnames = rep("chr_SYNTHETIC",length(listMu)),
                      ranges=IRanges(start=listMu,end=listMu),
                      strand=rep("*", length(listMu)))

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
    exp$k <- 11
    listMu <- c(10077, 10236,
                10406, 10571,
                10744, 10842,
                10846, 10896,
                10906, 11410,
                11580)
    exp$mu <- GRanges(seqnames = rep("chr_SYNTHETIC",length(listMu)),
                      ranges=IRanges(start=listMu,end=listMu),
                      strand=rep("*", length(listMu)))

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

    obs <- postTreatment(reads = reads_demo_02,
                              resultRJMCMC = RJMCMC_result,
                              extendingSize = 20,
                              chrLength = 80000)

    listMu <- c(10072, 10241,
             10574,10665,
             10744)

    exp <- GRanges(seqnames = rep("chr_SYNTHETIC",length(listMu)),
                ranges=IRanges(start=listMu,end=listMu),
                strand=rep("*", length(listMu)))

    message <- paste0(" test.postTreatment_good_01() ",
                      "- posTreatment() did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}


test.postTreatment_good_02 <- function() {

    reads <- reads_demo_02
    seqlevels(reads)<- c(seqlevels(reads), "chr2")
    seqnames(reads)[301:450]<-Rle(values = "chr2", lengths = c(150))

    obs <- postTreatment(reads = reads, seqName = "chr_SYNTHETIC",
                         resultRJMCMC = RJMCMC_result,
                         extendingSize = 20,
                         chrLength = 80000)

    listMu <- c(10072, 10241, 10574, 10744)

    exp <- GRanges(seqnames = rep("chr_SYNTHETIC",length(listMu)),
                   ranges=IRanges(start=listMu,end=listMu),
                   strand=rep("*", length(listMu)))

    message <- paste0(" test.postTreatment_good_02() ",
                      "- posTreatment() did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}

###########################################################
## segmentation() function
###########################################################

test.segmentation_good_01 <- function() {

    sampleGRanges <- GRanges(seqnames = syntheticNucleosomeReads$dataIP$chr,
                    ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                    strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- segmentation(sampleGRanges, zeta =  147, delta = 20, maxLength = 20000)

    message <- paste0(" test.segmentation_good_01() ",
                      "- segmentation() did not generated expected result.")

    exp.len = 3
    exp.01.len = 10504
    exp.02.len = 11873
    exp.03.len = 9686

    checkTrue(class(obs)=="GRangesList", ms = message)
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
                    ranges = IRanges(
                                start = syntheticNucleosomeReads$dataIP$start,
                                end = syntheticNucleosomeReads$dataIP$end),
                            strand = syntheticNucleosomeReads$dataIP$strand)

    obs <- segmentation(sampleGRanges, zeta =  142, delta = 40,
                            maxLength = 15000)

    message <- paste0(" test.segmentation_good_02() ",
                        "- segmentation() did not generated expected result.")

    exp.len = 4
    exp.01.len = 7972
    exp.02.len = 8572
    exp.03.len = 9362
    exp.04.len = 6390

    checkTrue(class(obs)=="GRangesList", ms = message)
    checkEquals(length(obs), exp.len, ms = message)
    checkEquals(length(obs[[1]]), exp.01.len, ms = message)
    checkTrue(is(obs[[1]],"GRanges"), ms = message)
    checkEquals(length(obs[[2]]), exp.02.len, ms = message)
    checkTrue(is(obs[[2]],"GRanges"), ms = message)
    checkEquals(length(obs[[3]]), exp.03.len, ms = message)
    checkTrue(is(obs[[3]],"GRanges"), ms = message)
}


###########################################################
## rjmcmcCHR() function
###########################################################

test.rjmcmcCHR_good_01 <- function() {

        temp_dir <- "test_rjmcmcCHR_good_01"
        tryCatch({
            reads <- GRanges(syntheticNucleosomeReads$dataIP[1:500,])

            obs <- rjmcmcCHR(reads = reads,
                        zeta = 147, delta = 50, maxLength = 1200,
                        dirOut = temp_dir,
                        nbrIterations = 1000, lambda = 3, kMax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5,
                        vSeed = 10113, nbCores = 2, saveAsRDS = FALSE)
            message <- paste0(" test.rjmcmcCHR_good_01() ",
                        "- rjmcmcCHR() did not generated expected result.")
            exp.k <- 2
            exp.kPost <- 1
            exp.mu <- c(1082, 1193)
            exp.muPost <- c(1188)
            exp.muStrand <- c('*', '*')
            exp.muPost <- c(1188)
            exp.muPostStrand <- c('*')

            checkTrue(is.list(obs), ms = message)
            checkEquals(obs$k, exp.k, ms = message)
            checkEquals(obs$kPost, exp.kPost, ms = message)
            checkEquals(start(obs$mu), exp.mu, ms = message)
            checkEquals(start(obs$muPost), exp.muPost, ms = message)
            checkEquals(end(obs$mu), exp.mu, ms = message)
            checkEquals(as.vector(strand(obs$mu)), exp.muStrand, ms = message)
            checkEquals(start(obs$muPost), exp.muPost, ms = message)
            checkEquals(end(obs$muPost), exp.muPost, ms = message)
            checkEquals(as.vector(strand(obs$muPost)), exp.muPostStrand, ms = message)

        }, finally = {
            if (dir.exists(temp_dir)) {
                unlink(temp_dir, recursive = TRUE, force = FALSE)
            }
        }
    )
    ## Double check
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

test.rjmcmcCHR_good_02 <- function() {

    temp_dir <- "test_rjmcmcCHR_good_02"
    tryCatch({
        reads <- GRanges(syntheticNucleosomeReads$dataIP[1:500,])
        seqlevels(reads)<- c(seqlevels(reads), "chr2")
        seqnames(reads)[400:500]<-Rle(values = "chr2", lengths = c(101))
        obs <- rjmcmcCHR(reads = reads, seqName = "chr_SYNTHETIC",
                         zeta = 147, delta = 50, maxLength = 1200,
                         dirOut = temp_dir,
                         nbrIterations = 1000, lambda = 3, kMax = 30,
                         minInterval = 146, maxInterval = 292, minReads = 5,
                         vSeed = 10113, nbCores = 2, saveAsRDS = FALSE)
        message <- paste0(" test.rjmcmcCHR_good_02() ",
                          "- rjmcmcCHR() did not generated expected result.")
        exp.k <- 2
        exp.kPost <- 2
        exp.mu <- c(1075, 1250)
        exp.muPost <- c(1075, 1250)
        exp.muStrand <- c('*', '*')
        exp.muPost <- c(1075, 1250)
        exp.muPostStrand <- c('*', '*')

        checkTrue(is.list(obs), ms = message)
        checkEquals(obs$k, exp.k, ms = message)
        checkEquals(obs$kPost, exp.kPost, ms = message)
        checkEquals(start(obs$mu), exp.mu, ms = message)
        checkEquals(start(obs$muPost), exp.muPost, ms = message)
        checkEquals(end(obs$mu), exp.mu, ms = message)
        checkEquals(as.vector(strand(obs$mu)), exp.muStrand, ms = message)
        checkEquals(start(obs$muPost), exp.muPost, ms = message)
        checkEquals(end(obs$muPost), exp.muPost, ms = message)
        checkEquals(as.vector(strand(obs$muPost)), exp.muPostStrand, ms = message)

    }, finally = {
        if (dir.exists(temp_dir)) {
            unlink(temp_dir, recursive = TRUE, force = FALSE)
        }
    }
    )
    ## Double check
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}
