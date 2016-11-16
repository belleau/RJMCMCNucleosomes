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

    std::vector<double>::const_iterator d_startF, d_endF;
    std::vector<double>::const_iterator d_startR, d_endR;
    long d_sizeF, d_sizeR, d_dimN;

    double d_mu;
    double d_sigmaF, d_sigmaR;

    gsl_rng *d_rng;
public:
    Nucleosome(double pos, SegmentSeq const &segSeq, gsl_rng *rng);
    virtual ~Nucleosome();
    void setStartF(std::vector<double>::const_iterator startF);
    void setEndF(std::vector<double>::const_iterator endF);
    void setStartR(std::vector<double>::const_iterator startR);
    void setEndR(std::vector<double>::const_iterator endR);
    void setSizeF(int);
    void setSizeR(int);
    void setDimN(long);
    int sizeF();
    int sizeR();
    long dimN();
    void setFStartPos(std::vector<double>::const_iterator fStart, std::vector<double>::const_iterator fEnd, int n);
    void setRStartPos(std::vector<double>::const_iterator fStart, std::vector<double>::const_iterator fEnd, int n);
    double mu();

    double varRead(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end, int n);

    void evalSigmaF();
    void evalSigmaR();
    void setSigmaF(double sigmaF);
    double sigmaF();
    void setSigmaR(double sigmaR);
    double sigmaR();
    double zeta();
    double deltaMin();
    double deltaMax();
    std::vector<double>::const_iterator startF() ;
    std::vector<double>::const_iterator endF();


protected:

    std::vector<double>::const_iterator startR();
    std::vector<double>::const_iterator endR();

    std::vector<double>::const_iterator beginFR();
    std::vector<double>::const_iterator endFR();
    long sizeFR();

    std::vector<double>::const_iterator beginRR();
    std::vector<double>::const_iterator endRR();
    long sizeRR();

    gsl_rng * rng(){
        return(d_rng);
    };

};

}/* namespace space_process */

#endif /* NUCLEOSOME_H_ */
