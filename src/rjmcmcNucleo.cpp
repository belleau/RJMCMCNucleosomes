#include <Rcpp.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <gsl/gsl_randist.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "PartitionAll.h"
#include "SimulationNucleoD.h"


#include "SegmentSeq.h"

//typedef space_process::SpaceState regionState;

using namespace Rcpp;
using namespace std;
using namespace space_process;

// Below is a simple example of exporting a C++ function to R. You can
//

// [[Rcpp::export]]
List rjmcmcNucleo(SEXP startPosForwardReads, SEXP startPosReverseReads,
                        long nbrIterations, int kMax, int lambda,
                        int minInterval, int maxInterval, int minReads = 5,
                        bool adaptIterationsToReads = true, int vSeed = -1) {
	/*********************************************************************
	 * init R var
	 *********************************************************************/
	NumericVector startFReads(startPosForwardReads); // *startFReads = new IntegerVector(startPosForwardReads);
    NumericVector startRReads(startPosReverseReads); // *startRReads = new IntegerVector(startPosReverseReads);

    std::vector<double> fReads = Rcpp::as<std::vector<double> >(startFReads);
    std::vector<double> rReads = Rcpp::as<std::vector<double> >(startRReads);

    long nr = 0; // Number of read
    //long tot = 0;
    //int kMaxO = 0; // the biggest K
    //int mO = 2;
    //kMax = 30;


    /*********************************************************************
	 * init random number generator
	 *********************************************************************/
    const gsl_rng_type * T;
    gsl_rng *rng;
	//long seed;

	T = gsl_rng_default;

	rng = gsl_rng_alloc (T);     // pick random number generator
	if(vSeed <= 0){
		vSeed = time (NULL) * getpid();
	}
	gsl_rng_set (rng, vSeed);

	/*********************************************************************
	 * Space Nucleosome and Segment
	 *********************************************************************/

	std::vector< PartitionAll<NucleoDirichletPA> *> res; // vector of space nucleosome accepted

	SegmentSeq seg(fReads, rReads, 147);                 // Reads of the segment



	nr = seg.sizeFReads() + seg.sizeRReads();
    if(adaptIterationsToReads){
    	if(nr <= 12){
    		nbrIterations = 1000;
    	}
    }

    List resO;

	SimulationNucleoD<PartitionAll<NucleoDirichletPA> > pourv(seg, rng, kMax, nbrIterations);
	if(pourv.initMu(lambda)){
		pourv.simulate();
		pourv.statSim();
		resO = pourv.simRapport();
	}
	else{
		resO = R_NilValue;
	}

    return resO;
}


static const R_CallMethodDef callMethods[] = {
    {".rjmcmcNucleo", (DL_FUNC) &rjmcmcNucleo, 10},
    {NULL, NULL, 0}
};

void R_init_rjmcmcNucleo(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
