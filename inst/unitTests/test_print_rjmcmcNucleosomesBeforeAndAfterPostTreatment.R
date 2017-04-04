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

data("reads_demo_02")

####################################################################
## print.rjmcmcNucleosomesBeforeAndAfterPostTreatment() function
####################################################################

test.print_rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value <- function() {

    ## Nucleosome detection
    result <- rjmcmc(reads = reads_demo_02,
                    seqName = "chr_SYNTHETIC", nbrIterations = 1000,
                    lambda = 2, kMax = 30, minInterval = 146,
                    maxInterval = 490, minReads = 3, vSeed = 11)

    ##Post-treatment function which merged closely positioned nucleosomes
    postResult <- postTreatment(reads = reads_demo_02,
                            seqName = "chr_SYNTHETIC", result, 100, 73500)

    result <- print(postResult)

    message <- paste0("test.print_rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value() ",
                           "- print method did not returned expected value")

    checkEquals(postResult, result, msg = message)
}
