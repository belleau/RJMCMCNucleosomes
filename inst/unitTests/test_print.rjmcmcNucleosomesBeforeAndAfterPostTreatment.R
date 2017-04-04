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

data("reads_demo_01")

####################################################################
## print.rjmcmcNucleosomesBeforeAndAfterPostTreatment() function
####################################################################

test.print_rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value <- function() {

    ##directoryWithRDSFiles <- system.file("extdata", package = "RJMCMCNucleosomes")

    ##resultInit <- mergeAllRDSFilesFromDirectory(directoryWithRDSFiles)

    ##  forward <- start(reads_demo_01[strand(reads_demo_01) == "+"])
    ## reverse <- end(reads_demo_01[strand(reads_demo_01) == "-"])
    ## postResult <- postTreatment(startPosForwardReads = reads_demo$readsForward,
    ##                             startPosReverseReads = reads_demo$readsReverse,
    ##                             resultInit, 74, 73500)

    ##result <- print(postResult)

    ## message <- paste0(" test.print.rjmcmcNucleosomesBeforeAndAfterPostTreatment_returned_value() ",
    ##                   "- print method did not returned expected value")

    ##checkEquals(postResult, result, msg = message)

    ## TODO
    return(TRUE)
}
