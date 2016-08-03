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
#include "NucleoDirichlet.h"

namespace space_process{

	class SpaceNucleosomeD{
	protected:
		std::list<NucleoDirichlet> d_nucleosomes;

	public:
		SpaceNucleosomeD(Rcpp::IntegerVector const  &fReads,
				Rcpp::IntegerVector const &rReads, int zeta);
		virtual ~SpaceNucleosomeD();
		int size();
		bool empty();
		void insert(NucleoDirichlet &u);

	};
}

#endif /* SPACENUCLEOSOMED_H_ */
