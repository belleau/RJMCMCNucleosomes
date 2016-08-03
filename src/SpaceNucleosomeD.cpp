/*
 * SpaceNucleosome.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: belleau
 */

#include "SpaceNucleosomeD.h"

using namespace std;
using namespace space_process;
using namespace Rcpp;

SpaceNucleosomeD::SpaceNucleosomeD(Rcpp::IntegerVector const  &fReads,
		Rcpp::IntegerVector const &rReads, int zeta){
}

int SpaceNucleosomeD::size() {
	return(d_nucleosomes.size());
}

bool SpaceNucleosomeD::empty() {
	return(d_nucleosomes.empty());
}

void SpaceNucleosomeD::insert(NucleoDirichlet &u) {
	d_nucleosomes.push_back(u);
}

SpaceNucleosomeD::~SpaceNucleosomeD() {
	// TODO Auto-generated destructor stub
}

