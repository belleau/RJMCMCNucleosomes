/*
 * NucleoDirichlet.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */
#include "NucleoDirichlet.h"
using namespace std;

namespace space_process {

//typedef typename std::vector<double> containerD;
//typedef typename containerD::const_iterator iteratorD;

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

		//double sigma = sqrt( 1 /(1/sigmaF() + 1/sigmaR()));
		double sigma = sqrt( 1 /(sigmaF() + sigmaR()));
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
		double tmp = mu() - delta()/2.0;
		double sdF = sqrt(sigmaF());
		long i = 0;

		for(vector<double>::const_iterator it = beginFR(); it != endFR(); it++){ //run on all seq

			d_bF[i++] = gsl_ran_tdist_pdf(((*it - tmp ) / sdF), df()) / sdF;
		}

	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}
}

double NucleoDirichlet::tP(){
	return(d_tP);
}

void NucleoDirichlet::testF(){
	std::cout << "tmp " << mu() + delta()/2.0 << " delta " << delta() << "\n";
	/*cout << "Read ";
	for(vector<double>::const_iterator it = startF(); it != endF(); it++){ //run on all seq
		cout << ", " << (*it);
	}
	cout << "\n";*/
}

void NucleoDirichlet::evalBF1(){

	/* reinit d_bF */
	d_bF.clear();
	d_bF.resize(sizeFR());
	d_tP = 0;
	if(sigmaF() > 0.000001)
	{
		double tmp = mu() - delta()/2.0;
		double sdF = sqrt(sigmaF());
		long i = 0;
		//if(1/sdF > 0.5)
		//std::cout << "tmp " << 1/sdF << " df " << df() << "\n";
		//std::cout << "tmp ";
		//d_tP = mu(); // + delta()/2;

		d_tP = sdF;
		//std::cout << " sigmafBF " << sigmaF() << "\n";
		for(vector<double>::const_iterator it = beginFR(); it != endFR(); it++){ //run on all seq
			//std::cout << " t " << ((*it - tmp ) / sdF);
			//d_tP += ((*it - tmp ) / sdF);
			d_bF[i++] = gsl_ran_tdist_pdf(((*it - tmp ) / sdF), df()) / sdF;
		}
		//d_tP /= i;
		//std::cout << "\n";
		//cout << "Sigma " << sigmaF() << " " << sdF << "\n";
	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}
}


void NucleoDirichlet::setBF(std::vector<double> &bF){
	d_bF = bF;
}

std::vector<double> & NucleoDirichlet::bF(){
	return(d_bF);
}

std::vector<double>::const_iterator NucleoDirichlet::bFBegin() const{
	return(d_bF.begin());
}

std::vector<double>::const_iterator NucleoDirichlet::bFEnd() const{
	return(d_bF.end());
}

std::vector<double>::const_iterator NucleoDirichlet::bRBegin() const{
	return(d_bR.begin());
}

std::vector<double>::const_iterator NucleoDirichlet::bREnd() const{
	return(d_bR.end());
}


void NucleoDirichlet::evalBR(){

	/* reinit d_bF */
	d_bR.clear();
	d_bR.resize(sizeRR());

	if(sigmaR() > 0.000001){

		double tmp = mu() + delta()/2.0;
		double sdR = sqrt(sigmaR());
		long i = 0;
		for(vector<double>::const_iterator it = beginRR(); it != endRR(); it++){

			d_bR[i++] = gsl_ran_tdist_pdf(((*it - tmp ) / sdR), df()) / sdR;
		}
	}
	else{
		cerr << "sigmaF or sigmaR not bigger than 0\n";
		//exit(1);
	}
}

void NucleoDirichlet::setBR(std::vector<double> &bR){
	d_bR = bR;
}

std::vector<double> &NucleoDirichlet::bR(){
	return(d_bR);
}

} /* namespace space_process */
