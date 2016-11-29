#' RJMCMCNucleosomes: Bayesian hierarchical model for genome-wide
#' nucleosome positioning with high-throughput short-read data (MNase-Seq)
#'
#' This package does nucleosome positioning using informative
#' Multinomial-Dirichlet prior in a t-mixture with reversible jump
#' estimation of nucleosome positions for genome-wide profiling.
#'
#' @docType package
#'
#' @name RJMCMCNucleosomes-package
#'
#' @aliases RJMCMCNucleosomes-package RJMCMCNucleosomes
#'
#' @author  Pascal Belleau,
#' Rawane Samb,
#' Astrid Deschênes,
#' Khader Khadraoui,
#' Lajmi Lakhal and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} { for profiling of nucleosome positions for a
#'     segment}
#'     \item \code{\link{rjmcmcCHR}} { for profiling of nucleosome positions
#'     for a large region. The function will take care of spliting and
#'     merging.}
#'     \item \code{\link{segmentation}} { for spliting a \code{GRanges}
#'     containing reads in a list of smaller segments for
#'     the \code{rjmcmc} function.}
#'     \item \code{\link{postTreatment}} { for merging closely positioned
#'     nucleosomes}
#'     \item \code{\link{mergeRDSFiles}} { for merging nucleosome information
#'     from selected RDS files.}
#'     \item \code{\link{plotNucleosomes}} { for generating a graph containing
#'     the nucleosome positions and the read coverage.}
#' }
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib RJMCMCNucleosomes
#' @keywords package
NULL

#' Forward reads and reverse reads in \code{GRanges} format
#' (for demo purpose).
#'
#' A group of forward and reverse reads, in a \code{GRanges}, that can be
#' used to test the \code{rjmcmc} function.
#'
#' @name reads_demo_01
#'
#' @docType data
#'
#' @aliases reads_demo_01
#'
#' @format A \code{GRanges} containing forward and reverse reads.
#'
#' @return A \code{GRanges} containing forward and reverse reads.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
#' }
#'
#' @usage data(reads_demo_01)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_01)
#'
#' ## Nucleosome positioning
#' rjmcmc(reads = reads_demo_01, nbrIterations = 100, lambda = 3, kMax = 30,
#'             minInterval = 146, maxInterval = 292, minReads = 5)
#'
NULL

#' Forward reads and reverse reads in \code{GRanges} format
#' (for demo purpose).
#'
#' A group of forward and reverse reads that can be used to test the
#' \code{rjmcmc} function.
#'
#' @name reads_demo_02
#'
#' @docType data
#'
#' @aliases reads_demo_02
#'
#' @format A \code{GRanges} containing forward and reverse reads.
#'
#' @return A \code{GRanges} containing forward and reverse reads.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
#'     \item \code{\link{rjmcmcCHR}} { for profiling of nucleosome positions
#'     for a large region. The function will take care of spliting and
#'     merging.}
#'     \item \code{\link{segmentation}} { for spliting a \code{GRanges}
#'     containing reads in a list of smaller segments for
#'     the \code{rjmcmc} function.}
#'     \item \code{\link{postTreatment}} { for merging closely positioned
#'     nucleosomes}
#'     \item \code{\link{mergeRDSFiles}} { for merging nucleosome information
#'     from selected RDS files.}
#'     \item \code{\link{plotNucleosomes}} { for generating a graph containing
#'     the nucleosome positions and the read coverage.}
#' }
#'
#' @usage data(reads_demo_02)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_02)
#'
#' ## Nucleosome positioning
#' ## Since there is only one chromosome present in reads_demo_02, the name
#' ## of the chromosome does not need to be specified
#' rjmcmc(reads = reads_demo_02, nbrIterations = 150, lambda = 3, kMax = 30,
#'             minInterval = 144, maxInterval = 290, minReads = 6)
#'
NULL

#' Nucleosomes obtained by running RJMCMC function using reads
#' from reads_demo_02 dataset (for demo purpose).
#'
#' A \code{list} of \code{class}
#' "rjmcmcNucleosomes" which contains the information about the
#' detected nucleosomes.
#'
#' @name RJMCMC_result
#'
#' @docType data
#'
#' @aliases RJMCMC_result
#'
#' @format A \code{list} of \code{class} "rjmcmcNucleosomes" containing:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{k} a \code{integer}, the final estimation of the number
#' of nucleosomes. \code{0} when no nucleosome is detected.
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes. \code{NA} when no nucleosome is
#' detected.
#' \item \code{k_max} a \code{integer}, the maximum number of nucleosomes
#' obtained during the iteration process. \code{NA} when no nucleosome is
#' detected.
#' }
#'
#' @return A \code{list} of \code{class} "rjmcmcNucleosomes" containing:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{k} a \code{integer}, the final estimation of the number
#' of nucleosomes. \code{0} when no nucleosome is detected.
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes. \code{NA} when no nucleosome is
#' detected.
#' \item \code{k_max} a \code{integer}, the maximum number of nucleosomes
#' obtained during the iteration process. \code{NA} when no nucleosome is
#' detected.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
#'     \item \code{\link{rjmcmcCHR}} { for profiling of nucleosome positions
#'     for a large region. The function will take care of spliting and
#'     merging.}
#'     \item \code{\link{segmentation}} { for spliting a \code{GRanges}
#'     containing reads in a list of smaller segments for
#'     the \code{rjmcmc} function.}
#'     \item \code{\link{postTreatment}} { for merging closely positioned
#'     nucleosomes}
#'     \item \code{\link{mergeRDSFiles}} { for merging nucleosome information
#'     from selected RDS files.}
#'     \item \code{\link{plotNucleosomes}} { for generating a graph containing
#'     the nucleosome positions and the read coverage.}
#' }
#'
#' @usage data(RJMCMC_result)
#' @docType data
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(RJMCMC_result)
#' data(reads_demo_02)
#'
#' ## Results before post-treatment
#' RJMCMC_result$mu
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' postResult <- postTreatment(reads = reads_demo_02,
#'     extendingSize = 60, chrLength = 100000, resultRJMCMC = RJMCMC_result)
#'
#' ## Results after post-treatment
#' postResult
#'
NULL


#' Simulated dataset of reads generated by \code{nucleoSim} package
#' (for demo purpose).
#'
#' A \code{list} of \code{class}
#' "syntheticNucReads" which contains the information about synthetic reads
#' related to nucleosomes. The datset has been created using a total of 300
#' well-positioned nucleosomes, 30 fuzzy nucleosomes with variance of reads
#' following a Normal distribution.
#'
#' @name syntheticNucleosomeReads
#'
#' @docType data
#'
#' @aliases syntheticNucleosomeReads
#'
#' @format A \code{list} containing:
#' \itemize{
#'     \item \code{call} the called that generated the dataset.
#'     \item \code{dataIP} a \code{data.frame} with the chromosome name, the
#'     starting and ending positions and the direction of all forward and
#'     reverse reads for all well-positioned and fuzzy nucleosomes. Paired-end
#'     reads are identified with an unique id.
#'     \item \code{wp} a \code{data.frame} with the positions of all the
#'     well-positioned nucleosomes, as well as the number of paired-reads
#'     associated to each one.
#'     \item \code{fuz} a \code{data.frame} with the positions of all the
#'     fuzzy nucleosomes, as well as the number of paired-reads associated
#'     to each one.
#'     \item \code{paired} a \code{data.frame} with the starting and ending
#'     positions of the reads used to generate the paired-end reads.
#'     Paired-end reads are identified with an unique id.
#' }
#'
#' @return A \code{list} containing:
#' \itemize{
#'     \item \code{call} the called that generated the dataset.
#'     \item \code{dataIP} a \code{data.frame} with the chromosome name, the
#'     starting and ending positions and the direction of all forward and
#'     reverse reads for all well-positioned and fuzzy nucleosomes. Paired-end
#'     reads are identified with an unique id.
#'     \item \code{wp} a \code{data.frame} with the positions of all the
#'     well-positioned nucleosomes, as well as the number of paired-reads
#'     associated to each one.
#'     \item \code{fuz} a \code{data.frame} with the positions of all the
#'     fuzzy nucleosomes, as well as the number of paired-reads associated
#'     to each one.
#'     \item \code{paired} a \code{data.frame} with the starting and ending
#'     positions of the reads used to generate the paired-end reads.
#'     Paired-end reads are identified with an unique id.
#' }
#'
#' @usage data(syntheticNucleosomeReads)
#'
#' @keywords datasets
#'
NULL
