/*
 * SpaceNucleosomeD.h
 *
 * Operation nucleosome from the space
 * Now only assign reads to nucleosomes
 *
 *  Created on: Jul 26, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOMED_H_
#define SPACENUCLEOSOMED_H_

#include <Rcpp.h>
#include "SpaceNucleosome.h"
#include "NucleoDirichlet.h"

namespace space_process{

	template<typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
	class SpaceNucleosomeD: public SpaceNucleosome<NucleoD>{

	public:

		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta){
		};

		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta, long sizeFReads, long sizeRReads)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta, sizeFReads, sizeRReads){
		};

		virtual ~SpaceNucleosomeD(){};


	};
}

#endif /* SPACENUCLEOSOMED_H_ */
