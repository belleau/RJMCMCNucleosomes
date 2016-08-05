/*
 * NucleoDirichlet.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */
#include "NucleoDirichlet.h"
#include <iostream>
#include <gsl/gsl_randist.h>
using namespace std;

namespace space_process {

NucleoDirichlet::NucleoDirichlet(double mu, int df, SegmentSeq &segSeq):
	Nucleosome(mu, segSeq){
	d_df = df;
}

NucleoDirichlet::~NucleoDirichlet() {
	// TODO Auto-generated destructor stub
}

void NucleoDirichlet::setAvg(double avg){
	d_avg = avg;
}

double NucleoDirichlet::NucleoDirichlet::avg(){
	return(d_avg);
}

void NucleoDirichlet::setDelta(double delta){
	d_delta = delta;
}

double NucleoDirichlet::delta(){
	return(d_delta);
}

void NucleoDirichlet::setDf(int df){
	d_df = df;
}

int NucleoDirichlet::df(){
	return(d_df);
}
/*
void NucleoDirichlet::setBF(double bF){
	d_bF = bF;
}

double NucleoDirichlet::bF(){
	return(d_bF);
}

void NucleoDirichlet::setBR(double bR){
	d_bR = bR;
}

double NucleoDirichlet::bR(){
	return(d_bR);
}
*/
void NucleoDirichlet::evalSigmaF(){
	setSigmaF(-1.0);
	if(d_df > 0){
		setSigmaF(varRead(startF(), endF(), sizeF()) * (d_df - 2) / d_df);
	}
}

void NucleoDirichlet::evalSigmaR(){
	setSigmaR(-1.0);
	if(d_df > 0){
		setSigmaR(varRead(startR(), endR(), sizeR()) * (d_df - 2) / d_df);
	}
}

void NucleoDirichlet::evalBF(std::vector<double> const &fReads){
	// dt((paramValues$startPSF - muValue[m] + deltaValue[m]/2)/sqrt(sigmafValue[m]), dfValue[m]) / sqrt(sigmafValue[m]

	// std::transform(myv1.begin(), myv1.end(), myv1.begin(),
    //std::bind1st(std::multiplies<T>(),3));


	//d_bF
	double tmp = mu() + delta();

	for(vector<double>::const_iterator it = fReads.begin(); it != fReads.end(); it++){

	}
	//double bF = gsl_ran_tdist_pdf(x, df());

}

void NucleoDirichlet::evalBR(){
	 // dt((paramValues$startPSF - muValue[m] + deltaValue[m]/2)/sqrt(sigmafValue[m]), dfValue[m]) / sqrt(sigmafValue[m]

	//double x =
	//double bF = gsl_ran_tdist_pdf(x, df());

}


} /* namespace space_process */
