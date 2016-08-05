/*
 * SpaceDirichlet.h
 *
 *  Created on: Jul 28, 2016
 *      Author: belleau
 */

#ifndef SPACEDIRICHLET_H_
#define SPACEDIRICHLET_H_

#include <gsl/gsl_randist.h>
#include "SpaceState.h"
#include "SpaceNucleosomeD.h"

namespace space_process {

template<typename SpaceModif>  /***** BEWARE SpaceModif
									Must inherit
									from SpaceNucleosomeD *****/
class SpaceDirichlet: public SpaceState<SpaceModif> {

public:
	SpaceDirichlet(std::vector<double> const &fReads,
			std::vector<double> const &rReads, int zeta)
	:SpaceState<SpaceModif>(fReads, rReads, zeta){
	};

	virtual ~SpaceDirichlet(){};

	double newMu(double minPos, double maxPos){
		return(gsl_ran_flat(this->rng(), minPos, maxPos));
	};
	double newMu(){
		return(gsl_ran_flat(this->rng(), this->minPos(), this->maxPos()));
	};
};

} /* namespace space_process */

#endif /* SPACEDIRICHLET_H_ */
