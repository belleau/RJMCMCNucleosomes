/*
 * SpaceNucleosome.h
 *
 * Operation nucleosome from the space
 * Now only assign reads to nucleosomes
 *
 *  Created on: Jul 26, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOME_H_
#define SPACENUCLEOSOME_H_

#include "SpaceState.h"

namespace space_process{

	class SpaceNucleosome: virtual public SpaceState {

	public:
		SpaceNucleosome(Rcpp::IntegerVector const  &fReads,
				Rcpp::IntegerVector const &rReads, int zeta);
		virtual ~SpaceNucleosome();
	};
}

#endif /* SPACENUCLEOSOME_H_ */
