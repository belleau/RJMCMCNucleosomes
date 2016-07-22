#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
//

// [[Rcpp::export]]
List rjmcmcNucleosomes(SEXP startPosForwardReads, SEXP startPosReverseReads, long nbrIterations, int kMax, int lambda, int minInterval, int maxInterval, int minReads = 5, bool adaptIterationsToReads = true) {
    Rcpp::IntegerVector startFReads(startPosForwardReads);
    Rcpp::IntegerVector startRReads(startPosReverseReads);

    int nf, nr;
    long tot;
    regionState currentState;
    currentState.insert(nf);
    nf = startFReads.size();
    nr = startRReads.size();
    tot = nbrIterations + kMax;
    List nbSeq = List::create(Rcpp::Named("nf") = nf, Rcpp::Named("nr") = nr, Rcpp::Named("tot") = tot, Rcpp::Named("tot") = currentState.empty());
    return nbSeq;
}
