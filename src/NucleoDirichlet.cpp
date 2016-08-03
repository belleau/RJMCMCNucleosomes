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

NucleoDirichlet::NucleoDirichlet(double mu, int df):
	Nucleosome(mu){
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


} /* namespace space_process */
