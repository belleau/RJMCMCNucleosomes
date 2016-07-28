/*
 * Nucleosome.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

class Nucleosome {
	double mu;
	long *d_fStartPos, *d_rStartPos;
	double sigmaf, sigmar;
public:
	Nucleosome(double pos);
	virtual ~Nucleosome();

};

#endif /* NUCLEOSOME_H_ */
