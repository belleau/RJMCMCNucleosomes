#include <Rcpp.h>
#include <iostream>
//#include "SpaceState.h"
#include "SpaceNucleosomeD.h"
#include "PartitionAll.h"

#include "NucleoDirichletPA.h"
#include "Factory.h"
#include "bla1.h"

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
    std::vector<double> fReads = Rcpp::as<std::vector<double> >(startFReads);
    std::vector<double> rReads = Rcpp::as<std::vector<double> >(startRReads);

    int nf=1, nr;
    long tot;
    //startFReads[1] = 1;



    //space_process::SpaceState currentState(startFReads, startRReads, 147);
    //space_process::SpaceNucleosome currentState(startFReads, startRReads, 147);
    //space_process::PartitionAll bla(startFReads, startRReads, 147);

    //space_process::Factory<space_process::SpaceNucleosomeD, space_process::PartitionAll> truc;
    // space_process::NucleoDirichletPA

    //space_process::PartitionAll<space_process::NucleoDirichletPA> bla(fReads, rReads, 147);
    //cout << "Bla " << fReads[0] << "\n";
    space_process::PartitionAll<space_process::NucleoDirichletPA> currentState(fReads, rReads, 147);

    currentState.insertD(10010,3);
    currentState.insertD(10020,3);
    currentState.insertD(10030,3);
    currentState.insertD(10040,3);
    currentState.insertD(10050,3);
    currentState.insertD(10060,3);
    currentState.evalPriorMuDensity();

    //space_process::SpaceDirichlet<space_process::PartitionAll<space_process::NucleoDirichletPA> > currentState(fReads, rReads, 147);
    //cout << "Bla " << fReads[0] << "\n";
    //currentState.initMu( 3);
    //cout << "Aye " << startFReads[1] << "\n";
    //bla.initMu(currentState.newMu(), 3);
    //double mu = currentState.newMu();
    //cout << " Mu " << mu << "\n";
    //currentState.insert(nf);
    nf = startFReads.size();
    //nf = currentState.getP();
    nr = startRReads.size();
    tot = nbrIterations + kMax;
    List nbSeq = List::create(Rcpp::Named("nf") = nf, Rcpp::Named("nr") = nr, Rcpp::Named("tot") = tot);
    return nbSeq;
}
