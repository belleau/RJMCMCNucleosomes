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
	double d_bF, d_bR; /* Kbr divided by w */
	double d_avg;
	double d_delta;
public:
	NucleoDirichlet(double mu, int df);
	virtual ~NucleoDirichlet();
	double testT();

	void setAvg(double avg);
	double avg();

	void setDelta(double delta);
	double delta();

	void setDf(int df);
	int df();

	void setBF(double bF);
	double bF();
	void setBR(double bR);
	double bR();

	void evalSigmaF();
	void evalSigmaR();
//	void setBf();
};

} /* namespace space_process */

#endif /* NUCLEODIRICHLET_H_ */
