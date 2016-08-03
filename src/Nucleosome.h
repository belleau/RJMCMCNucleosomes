/*
 * Nucleosome.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

#include <Rcpp.h>

namespace space_process {

class Nucleosome {
	std::vector<double>::iterator d_startF, d_endF;
	std::vector<double>::iterator d_startR, d_endR;
	int d_sizeF, d_sizeR;
	double d_mu;
	double d_sigmaF, d_sigmaR;


public:
	Nucleosome(double pos);
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

	double varRead(std::vector<double>::iterator start, std::vector<double>::iterator end, int n);

	void evalSigmaF();
	void evalSigmaR();
	void setSigmaF(double sigmaF);
	double sigmaF();
	void setSigmaR(double sigmaR);
	double sigmaR();
protected:
	std::vector<double>::iterator startF() ;
	std::vector<double>::iterator endF();
	std::vector<double>::iterator startR();
	std::vector<double>::iterator endR();

};

}/* namespace space_process */

#endif /* NUCLEOSOME_H_ */
