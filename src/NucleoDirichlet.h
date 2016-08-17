/*
 * NucleoDirichlet.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEODIRICHLET_H_
#define NUCLEODIRICHLET_H_


#include "Nucleosome.h"
#include <iostream>
#include <gsl/gsl_randist.h>

namespace space_process {

class NucleoDirichlet: public Nucleosome {
	typedef typename std::vector<double> containerD;
	typedef typename containerD::const_iterator iteratorD;
	int d_df;
	containerD d_bF, d_bR; /* Kbr divided by w */
	double d_delta;
public:
	NucleoDirichlet(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng);
	NucleoDirichlet(double mu, SegmentSeq const &segSeq, gsl_rng *rng);
	virtual ~NucleoDirichlet();
	double testT();


	void setDelta(double delta);
	double delta();

	void setDf(int df);
	int df();
	/*
	void setBF(double bF);
	double bF();
	void setBR(double bR);
	double bR();
*/

	void evalSigmaF();
	void evalSigmaR();

	void evalDelta();
	void evalBF();
	iteratorD bFBegin() const;
	iteratorD bFEnd() const;

	void evalBR();
	iteratorD bRBegin() const;
	iteratorD bREnd() const;

//	void setBf();
};

} /* namespace space_process */

#endif /* NUCLEODIRICHLET_H_ */
