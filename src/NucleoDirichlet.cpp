/*
 * NucleoDirichlet.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */
#include "NucleoDirichlet.h"
using namespace std;

namespace space_process {

typedef typename std::vector<double> containerD;
typedef typename containerD::const_iterator iteratorD;

NucleoDirichlet::NucleoDirichlet(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng):
	Nucleosome(mu, segSeq, rng), d_df(df){
}

NucleoDirichlet::NucleoDirichlet(double mu, SegmentSeq const &segSeq, gsl_rng *rng):
	Nucleosome(mu, segSeq, rng){
}

NucleoDirichlet::~NucleoDirichlet() {
	// TODO Auto-generated destructor stub
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
	if(d_df > 2){
		setSigmaF(varRead(startF(), endF(), sizeF()) * (d_df - 2) / d_df);
	}
}

void NucleoDirichlet::evalSigmaR(){
	setSigmaR(-1.0);
	if(d_df > 2){
		setSigmaR(varRead(startR(), endR(), sizeR()) * (d_df - 2) / d_df);
	}
}

void NucleoDirichlet::evalDelta(){
	/* compute truncated normal */
	if(sigmaF() > 0 && sigmaR() > 0){
		double sigma = sqrt( 1 /(1/sigmaF() + 1/sigmaR()));
		double tmpDelta = 147;
		do{

			tmpDelta = zeta() + gsl_ran_gaussian(rng(), sigma);

		}while(tmpDelta > deltaMin() &&  tmpDelta < deltaMax());
		setDelta(tmpDelta);
	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}

}

void NucleoDirichlet::evalBF(){

	/* reinit d_bF */
	d_bF.clear();
	d_bF.resize(sizeFR());
	if(sigmaF() > 0.000001)
	{
		double tmp = mu() + delta()/2;
		double sdF = sqrt(sigmaF());
		long i = 0;
		for(vector<double>::const_iterator it = beginFR(); it != endFR(); it++){
			//d_bF.push_back(gsl_ran_tdist_pdf(((*it - tmp ) / sdF), df()));
			d_bF[i++] = gsl_ran_tdist_pdf(((*it - tmp ) / sdF), df()) / sdF;
		}
	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}
}

iteratorD NucleoDirichlet::bFBegin() const{
	return(d_bF.begin());
}

iteratorD NucleoDirichlet::bFEnd() const{
	return(d_bF.end());
}

iteratorD NucleoDirichlet::bRBegin() const{
	return(d_bR.begin());
}

iteratorD NucleoDirichlet::bREnd() const{
	return(d_bR.end());
}


void NucleoDirichlet::evalBR(){

	/* reinit d_bF */
	d_bR.clear();
	d_bR.resize(sizeRR());

	if(sigmaR() > 0.000001){

		double tmp = mu() + delta()/2;
		double sdR = sqrt(sigmaR());
		long i = 0;
		for(vector<double>::const_iterator it = beginRR(); it != endRR(); it++){
			//d_bF.push_back(gsl_ran_tdist_pdf(((*it - tmp ) / sdF), df()));
			d_bR[i++] = gsl_ran_tdist_pdf(((*it - tmp ) / sdR), df()) / sdR;
		}
	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}
}



} /* namespace space_process */
