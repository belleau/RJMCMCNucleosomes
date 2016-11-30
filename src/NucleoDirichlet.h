/*
 * NucleoDirichlet.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef NUCLEODIRICHLET_H_
#define NUCLEODIRICHLET_H_


#include "Nucleosome.h"
#include <Rcpp.h>
#include <gsl/gsl_randist.h>

namespace space_process {

class NucleoDirichlet: public Nucleosome {
    int d_df;
    std::vector<double> d_bF, d_bR; /* Kbr divided by w */
    double d_delta;
    double d_tP;
public:
    NucleoDirichlet(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng);
    NucleoDirichlet(double mu, SegmentSeq const &segSeq, gsl_rng *rng);
    virtual ~NucleoDirichlet();
    double testT();

    double tP();
    void setDelta(double delta);
    double delta();

    void setDf(int df);
    int df();

    void evalSigmaF();
    void evalSigmaR();

    void evalDelta();
    void evalBF();

    void setBF(std::vector<double> &bF);

    std::vector<double> &bF();

    std::vector<double>::const_iterator bFBegin() const;
    std::vector<double>::const_iterator bFEnd() const;

    void evalBR();
    void setBR(std::vector<double> &bR);
    std::vector<double> &bR();
    std::vector<double>::const_iterator bRBegin() const;
    std::vector<double>::const_iterator bREnd() const;

//	void setBf();
};

} /* namespace space_process */

#endif /* NUCLEODIRICHLET_H_ */
