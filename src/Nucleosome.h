/*
 * Nucleosome.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

#include <Rcpp.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "SegmentSeq.h"

namespace space_process {

class Nucleosome {

	SegmentSeq const &d_segSeq;

	std::vector<double>::iterator d_startF, d_endF;
	std::vector<double>::iterator d_startR, d_endR;
	long d_sizeF, d_sizeR;

	double d_mu;
	double d_sigmaF, d_sigmaR;

	gsl_rng *d_rng;
public:
	Nucleosome(double pos, SegmentSeq const &segSeq, gsl_rng *rng);
	virtual ~Nucleosome();
	void setStartF(std::vector<double>::iterator startF);
	void setEndF(std::vector<double>::iterator endF);
	void setStartR(std::vector<double>::iterator startR);
	void setEndR(std::vector<double>::iterator endR);
	void setSizeF(int);
	void setSizeR(int);
	int sizeF();
	int sizeR();
	void setFStartPos(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, int n);
	void setRStartPos(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, int n);
	double mu();

	double varRead(std::vector<double>::iterator start, std::vector<double>::iterator end, int n);

	void evalSigmaF();
	void evalSigmaR();
	void setSigmaF(double sigmaF);
	double sigmaF();
	void setSigmaR(double sigmaR);
	double sigmaR();
	double zeta();
	double deltaMin();
	double deltaMax();
protected:
	std::vector<double>::iterator startF() ;
	std::vector<double>::iterator endF();
	std::vector<double>::iterator startR();
	std::vector<double>::iterator endR();

	gsl_rng * rng(){
		return(d_rng);
	};
};

}/* namespace space_process */

#endif /* NUCLEOSOME_H_ */
