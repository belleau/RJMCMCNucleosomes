#include <Rcpp.h>
#include <iostream>
//#include "SpaceState.h"
//#include "SpaceNucleosome.h"
#include "PartitionAll.h"
#include "SpaceDirichlet.h"

//typedef space_process::SpaceState regionState;

using namespace Rcpp;
using namespace std;
//using namespace space_process;

// Below is a simple example of exporting a C++ function to R. You can
//

// [[Rcpp::export]]
List rjmcmcNucleo(SEXP startPosForwardReads, SEXP startPosReverseReads,
                        long nbrIterations, int kMax, int lambda,
                        int minInterval, int maxInterval, int minReads = 5,
                        bool adaptIterationsToReads = true) {
    IntegerVector startFReads(startPosForwardReads); // *startFReads = new IntegerVector(startPosForwardReads);
    IntegerVector startRReads(startPosReverseReads); // *startRReads = new IntegerVector(startPosReverseReads);

    int nf=1, nr;
    long tot;
    //startFReads[1] = 1;
    //space_process::SpaceState currentState(startFReads, startRReads, 147);
    //space_process::SpaceNucleosome currentState(startFReads, startRReads, 147);
    space_process::PartitionAll bla(startFReads, startRReads, 147);

    space_process::SpaceDirichlet currentState(startFReads, startRReads, 147, bla);
    //cout << "Aye " << startFReads[1] << "\n";
    bla.initMu(currentState.newMu(), 3);
    double mu = currentState.newMu();
    //cout << " Mu " << mu << "\n";
    //currentState.insert(nf);
    nf = startFReads.size();
    //nf = currentState.getP();
    nr = startRReads.size();
    tot = nbrIterations + kMax;
    List nbSeq = List::create(Rcpp::Named("nf") = nf, Rcpp::Named("nr") = nr, Rcpp::Named("tot") = tot);
    return nbSeq;
}
