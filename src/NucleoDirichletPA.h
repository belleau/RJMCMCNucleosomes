/*
 * NucleoDirichletPA.h
 *
 *  Created on: Jul 29, 2016
 *      Author: belleau
 */

#ifndef NUCLEODIRICHLETPA_H_
#define NUCLEODIRICHLETPA_H_

#include "NucleoDirichlet.h"

namespace space_process {

class NucleoDirichletPA: public NucleoDirichlet {
	public:
		NucleoDirichletPA(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng);
		virtual ~NucleoDirichletPA();
		double testT();
		void testFRStart();

};

} /* namespace space_process */

#endif /* NUCLEODIRICHLETPA_H_ */
