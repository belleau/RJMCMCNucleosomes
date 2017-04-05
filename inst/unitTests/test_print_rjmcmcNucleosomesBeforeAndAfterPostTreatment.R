###################################################
# Created by Astrid Deschenes
# 2015-03-08
###################################################

########################################################################
## Test the print.rjmcmcNucleosomesBeforeAndAfterPostTreatment function
########################################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMCNucleosomes" )
}

### }}}

data(syntheticNucleosomeReads)

####################################################################
## print.rjmcmcNucleosomesBeforeAndAfterPostTreatment() function
####################################################################

test.print_rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value <- function() {

    sampleGRanges <- GRanges(syntheticNucleosomeReads$dataIP)

    result <- rjmcmcCHR(reads = sampleGRanges, zeta = 147, delta=50,
                    maxLength=1200, nbrIterations = 50, lambda = 3, kMax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5, vSeed = 222,
                    nbCores = 1, saveAsRDS = FALSE)

    resultPrint <- print(result)

    message <- paste0("test.print_rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value() ",
                           "- print method did not returned expected value")

    checkEquals(result, resultPrint, msg = message)
}
