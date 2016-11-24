## ----loadingPackage, warning=FALSE, message=FALSE------------------------
library(RJMCMCNucleosomes)

## ----createSample, collapse=TRUE, message=FALSE--------------------------
## Load nucleoSim package
library(nucleoSim)

val.num       <- 50     ### Number of well-positioned nucleosomes
val.del       <- 10     ### Number of well-positioned nucleosomes to delete
val.var       <- 30     ### variance associated to well-positioned nucleosomes
val.fuz       <- 10     ### Number of fuzzy nucleosomes
val.fuz.var   <- 50     ### variance associated to fuzzy nucleosomes
val.max.cover <- 70     ### Maximum coverage for one nucleosome
val.nuc.len   <- 147    ### Distance between nucleosomes
val.len.var   <- 10     ### Variance associated to the length of the reads
val.lin.len   <- 20     ### The length of the DNA linker
val.rnd.seed  <- 100    ### Set seed when result needs to be reproducible
val.offset    <- 10000  ### The number of bases used to offset 
                        ### all nucleosomes and reads

## Create sample using a Normal distribution
sample <- nucleoSim::syntheticNucReadsFromDist(wp.num=val.num,
                                    wp.del=val.del, wp.var=val.var,
                                    fuz.num=val.del, fuz.var=val.fuz.var,
                                    max.cover=val.max.cover, 
                                    nuc.len=val.nuc.len,
                                    len.var=val.len.var, 
                                    lin.len=val.lin.len,
                                    rnd.seed=val.rnd.seed,
                                    distr="Normal", offset=val.offset)

## Create visual representation of the synthetic nucleosome sample
plot(sample)

## ----segment01, warning=FALSE, collapse=TRUE, message=FALSE--------------
## Load needed packages
library(GenomicRanges)

## Transform sample dataset into GRanges object
sampleGRanges <- GRanges(seqnames = sample$dataIP$chr, 
                          ranges = IRanges(start = sample$dataIP$start, 
                                           end = sample$dataIP$end), 
                          strand = sample$dataIP$strand)

## Segment sample into candidate regions
sampleSegmented <- segmentation(dataIP = sampleGRanges, zeta = 147, 
                                delta = 40, maxLength = 1000)

## Number of segments created
length(sampleSegmented)

## ----runRJMCMC01 , warning=FALSE, collapse=TRUE--------------------------
## Extract the first segment 
segment01 <- sampleSegmented[[1]]

## Run RJMCMC analysis
## A higher number of iterations is recommanded for real analysis
resultSegment01 <- rjmcmc(forwardandReverseReads  = segment01,
                            nbrIterations = 100000, lambda = 3, kMax = 30,
                            minInterval = 100, maxInterval = 200, minReads = 5,
                            vSeed = 1921)

## Print the predicted nucleosomes for the first segment
resultSegment01

## ----regroup01, warning=FALSE, collapse=TRUE, message=FALSE--------------
## The directory containing the results of all segments
## On RDS file has been created for each segment
directory <- system.file("extdata", "demo_vignette", 
                            package = "RJMCMCNucleosomes")

## Merging predicted nucleosomes from all segments
resultsAllMerged <- mergeAllRDSFilesFromDirectory(directory)

resultsAllMerged

## ----postProcess01, collapse=TRUE, message=FALSE-------------------------
## Split reads from the initial sample data into forward and reverse subsets
reads <- GRanges(sample$dataIP)

## Number of nucleosomes before the post-treatment
resultsAllMerged$k

## Use the post-treatment function
resultsPostTreatment <- postTreatment(forwardandReverseReads = reads,
                            resultRJMCMC = resultsAllMerged,
                            extendingSize = 15,
                            chrLength = max(start(reads), end(reads)) + 1000)


## Number of nucleosomes after the post-treatment
length(resultsPostTreatment)

## ----visualisation, collapse=TRUE, message=FALSE, fig.height=6.5, fig.width=8----
## Extract reads to create an IRange
reads <-IRanges(start = sample$dataIP$start, end=sample$dataIP$end)

## Merge predictions from before post-treatment and after post-treatment in 
## a list so that both results will be shown in graph
resultsBeforeAndAfter <- list(Sample = c(sample$wp$nucleopos, 
                                            sample$fuz$nucleopos),
                              BeforePostTreatment = resultsAllMerged$mu,    
                              AfterPostTreatment = resultsPostTreatment)

## Create graph using nucleosome positions and reads
## The plot will shows : 
## 1. nucleosomes from the sample
## 2. nucleosomes detected by rjmcmc() function
## 3. nucleosomes obtained after post-treament 
plotNucleosomes(resultsBeforeAndAfter, reads, 
                    names=c("Sample", "RJMCMC", "After Post-Treatment"))

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

