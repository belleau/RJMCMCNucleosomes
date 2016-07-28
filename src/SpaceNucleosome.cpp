/*
 * SpaceNucleosome.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: belleau
 */

#include "SpaceNucleosome.h"

using namespace std;
using namespace space_process;
using namespace Rcpp;

SpaceNucleosome::SpaceNucleosome(Rcpp::IntegerVector const  &fReads,
		Rcpp::IntegerVector const &rReads, int zeta):
		SpaceState(fReads, rReads, zeta){
}

SpaceNucleosome::~SpaceNucleosome() {
	// TODO Auto-generated destructor stub
}

