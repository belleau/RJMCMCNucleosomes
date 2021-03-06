% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rjmcmcMethodsIntern.R
\name{validateRJMCMCParameters}
\alias{validateRJMCMCParameters}
\title{Parameters validation for the \code{\link{rjmcmc}}
function}
\usage{
validateRJMCMCParameters(reads, seqName, nbrIterations, kMax, lambda,
  minInterval, maxInterval, minReads, adaptIterationsToReads, vSeed)
}
\arguments{
\item{reads}{a \code{GRanges} containing all forward
and reverse reads.The start positions of both reads are going to be used
for the analysis. Beware that the start position of
a reverse read is always higher that the end positition.}

\item{seqName}{a \code{character} string containing the label of the
chromosome, present in the \code{GRanges} object, that will be used. The
\code{NULL} value is accepted when only one seqname is
present in the \code{GRanges}; the only seqname present will be used.}

\item{nbrIterations}{a positive \code{integer} or \code{numeric}, the
number of iterations. Non-integer values of
\code{nbrIterations} will be casted to \code{integer} and truncated towards
zero.}

\item{kMax}{a positive \code{integer} or \code{numeric}, the maximum number
of nucleosomes per region. Non-integer values
of \code{kMax} will be casted to \code{integer} and truncated towards zero.}

\item{lambda}{a positive \code{numeric}, the theorical mean
of the Poisson distribution.}

\item{minInterval}{a \code{numeric}, the minimum distance between two
nucleosomes.}

\item{maxInterval}{a \code{numeric}, the maximum distance between two
nucleosomes.}

\item{minReads}{a positive \code{integer} or \code{numeric}, the minimum
number of reads in a potential canditate region. Non-integer values
of \code{minReads} will be casted to \code{integer} and truncated towards
zero.}

\item{adaptIterationsToReads}{a \code{logical} indicating if the number
of iterations must be modified in function of the number of reads.}

\item{vSeed}{a \code{integer}. A seed used when reproducible results are
needed. When a value inferior or equal to zero is given, a random integer
is used.}
}
\value{
\code{0} indicating that all parameters validations have been
successful.
}
\description{
Validation of all parameters needed by the public
\code{\link{rjmcmc}} function.
}
\examples{

reads <- GRanges(seqnames = Rle(c("chr1"), c(10)),
    ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "-", "+", "-")), c(1, 2, 2, 3, 2)))

## The function returns 0 when all paramaters are valid
RJMCMCNucleosomes:::validateRJMCMCParameters(reads = reads,
seqName = "chr1", nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1,
minInterval = 100, maxInterval = 200, adaptIterationsToReads = TRUE,
vSeed = 100)

## The function raises an error when at least one paramater is not valid
\dontrun{RJMCMCNucleosomes:::validateRJMCMCParameters(
reads = NA, seqName = "chr1",
nbrIterations = 2, kMax = 10, lambda = 1, minReads = 1, minInterval = 100,
maxInterval = 200, adaptIterationsToReads = TRUE, vSeed = -1)}

}
\author{
Astrid Deschenes
}
\keyword{internal}
