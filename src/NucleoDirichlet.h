/*
 * NucleoDirichlet.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEODIRICHLET_H_
#define NUCLEODIRICHLET_H_

#include "Nucleosome.h"

namespace space_process {

class NucleoDirichlet: public Nucleosome {
	int d_df;
	double d_bF, d_bR; /* Kbr without divided by w */
	double delta;
public:
	NucleoDirichlet(double mu);
	virtual ~NucleoDirichlet();
	double testT();
//	void setDf(int df);

//	void setBf();
};

} /* namespace space_process */

#endif /* NUCLEODIRICHLET_H_ */
