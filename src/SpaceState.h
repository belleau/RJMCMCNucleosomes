/*
 * SpaceState.h
 *
 *  Created on: Jul 22, 2016
 *      Author: Pascal Belleau
 */



#ifndef SPACESTATE_H_
#define SPACESTATE_H_

#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <unistd.h>

#include "SpaceNucleosome.h"

namespace space_process{

template<typename SpaceModif>  /***** BEWARE SpaceModif Must inherit from SpaceNucleosome *****/
class SpaceState {
	gsl_rng *d_rng;  // random number generator
	protected:
		SpaceModif d_nucleosomes;
		double d_minPos, d_maxPos;

		std::vector<double> const &d_startFReads;
		std::vector<double> const &d_startRReads; /* vector of
													 reads start
													 position
													 ( foward and
													 reverse) */
		const long d_sizeFReads, d_sizeRReads;
		const int d_zeta;

	public:
		SpaceState(std::vector<double> const &fReads, std::vector<double> const &rReads, int zeta)
		:d_startFReads(fReads), d_startRReads(rReads), d_zeta(zeta),
		 d_sizeFReads(fReads.size()), d_sizeRReads(rReads.size()),
		 d_nucleosomes(fReads, rReads, zeta, sizeFReads(), sizeRReads()){

			d_minPos = std::min(*(std::min_element(d_startFReads.begin(), d_startFReads.end())), *(std::min_element(d_startRReads.begin(), d_startRReads.end())));
			d_maxPos = std::max(*(std::max_element(d_startFReads.begin(), d_startFReads.end())), *(std::max_element(d_startRReads.begin(), d_startRReads.end())));

			const gsl_rng_type * T;
			long seed;

			T = gsl_rng_default;

			d_rng = gsl_rng_alloc (T);     // pick random number generator
			seed = time (NULL) * getpid();
			gsl_rng_set (d_rng, seed);                  // set seed
		}
;
		virtual ~SpaceState(){};

		long sizeFReads(){
			return(d_sizeFReads);
		};
		long sizeRReads(){
			return(d_sizeFReads);
		};
		double minPos(){
			return(d_minPos);
		}
		double maxPos(){
			return(d_maxPos);
		}
	protected:
		gsl_rng * rng(){
			return(d_rng);
		};
	};

}; /* namespace space_process */

#endif /* SPACESTATE_H_ */
