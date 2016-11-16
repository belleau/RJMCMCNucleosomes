/*
 * NucleoDirichletPA.cpp
 *
 *  Created on: Jul 29, 2016
 *      Author: belleau
 */

#include "NucleoDirichletPA.h"


namespace space_process {

NucleoDirichletPA::NucleoDirichletPA(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng):
    NucleoDirichlet(mu, df, segSeq, rng){
    setSizeF(-1);
    setSizeR(-1);

}

NucleoDirichletPA::NucleoDirichletPA(double mu, SegmentSeq const &segSeq, gsl_rng *rng):
    NucleoDirichlet(mu, segSeq, rng){
    setSizeF(-1);
    setSizeR(-1);

}


NucleoDirichletPA::~NucleoDirichletPA() {
}

double NucleoDirichletPA::aF(){
    return(d_aF);
}

void NucleoDirichletPA::setAF(double aF){
    d_aF = aF;
}

double NucleoDirichletPA::aR(){
    return(d_aR);
}
void NucleoDirichletPA::setAR(double aR){
    d_aR = aR;
}

void NucleoDirichletPA::setAvg(double avg){
    d_avg = avg;
}

double NucleoDirichletPA::avg(){
    return(d_avg);
}


} /* namespace space_process */
