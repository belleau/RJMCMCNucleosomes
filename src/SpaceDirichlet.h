/*
 * SpaceDirichlet.h
 *
 *  Created on: Jul 28, 2016
 *      Author: belleau
 */

#ifndef SPACEDIRICHLET_H_
#define SPACEDIRICHLET_H_

#include "SpaceState.h"
#include "SpaceNucleosomeD.h"

namespace space_process {

class SpaceDirichlet: public SpaceState {
	SpaceNucleosomeD d_nucleosomes;
public:
	SpaceDirichlet(Rcpp::IntegerVector const  &fReads,
			Rcpp::IntegerVector const &rReads, int zeta, SpaceNucleosomeD nucleosomes);
	virtual ~SpaceDirichlet();

	double newMu(double minPos, double maxPos);
	double newMu();
};

} /* namespace space_process */

#endif /* SPACEDIRICHLET_H_ */
