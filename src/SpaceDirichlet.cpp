/*
 * SpaceDirichlet.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: belleau
 */

#include "SpaceDirichlet.h"
#include <gsl/gsl_randist.h>

namespace space_process {

SpaceDirichlet::SpaceDirichlet(Rcpp::IntegerVector const  &fReads,
		Rcpp::IntegerVector const &rReads, int zeta, SpaceNucleosomeD nucleosomes):SpaceState(fReads, rReads, zeta), d_nucleosomes(nucleosomes) {
	// TODO Auto-generated constructor stub

}

SpaceDirichlet::~SpaceDirichlet() {

}

double SpaceDirichlet::newMu(double minPos, double maxPos) {
	return(gsl_ran_flat(d_rng, minPos, maxPos));
}
double SpaceDirichlet::newMu() {
	return(gsl_ran_flat(d_rng, (double)d_minPos, (double)d_maxPos));
}

} /* namespace space_process */
