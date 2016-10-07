/*
 * SimulationNucleoD.h
 *
 *  Created on: Sep 27, 2016
 *      Author: belleau
 */

#ifndef SIMULATIONNUCLEOD_H_
#define SIMULATIONNUCLEOD_H_


#include "SimulationNucleo.h"
#include <math.h>


namespace space_process {

template< typename NucleoSpace>
class SimulationNucleoD: public SimulationNucleo<NucleoSpace> {
	typedef std::vector<NucleoSpace *> NucleoSim;
	typedef typename NucleoSim::iterator itState;

	int d_kMax;
	Rcpp::List d_resStat;


public:
	SimulationNucleoD(SegmentSeq const &segSeq,
			gsl_rng * rng, int kMax, long nbIteration=10000);
	virtual ~SimulationNucleoD();

	bool initMu(int lambda, int df = 3);
	int kMax();
	bool sampler();
	void simulate();
	double computeRho();
	void initResStat();

	void statSim(){
		int i = 0;

		//initResStat();

		Rcpp::NumericVector listK = Rcpp::NumericVector( Rcpp::Dimension(this->sizeState()));
		Rcpp::NumericMatrix mu = Rcpp::NumericMatrix( Rcpp::Dimension(this->sizeState(), this->kMaxS()));
		Rcpp::IntegerVector listIt = Rcpp::IntegerVector( Rcpp::Dimension(this->sizeState()));
		Rcpp::IntegerVector nbK = Rcpp::IntegerVector(this->kMaxS());
		Rcpp::NumericVector muHat = Rcpp::NumericVector(Rcpp::Dimension(this->kMaxS(), this->kMaxS()));

		for(itState it = this->beginState(); it != this->endState();it++){
			listK[i] = (*it)->valK();
			listIt[i] = (*it)->iteration();
			nbK[(*it)->valK()-1] += (*it)->iteration();

			std::vector<double> tmp = (*it)->mu();

			for(int j = 0; j < this->kMaxS(); j++){
			   if(j < listK[i]){
				   mu[i  + j * this->sizeState()] = tmp[j];
				   muHat[((*it)->valK()-1) + j * this->kMaxS()] += (*it)->iteration() * tmp[j];
			   }
			   else{
				   mu[i + j * this->sizeState()] = 0;
			   }
		   }
			i++;
		}

		for(int j = 0; j < this->kMaxS(); j++){
		   for(int l = 0; l < this->kMaxS(); l++){
			   if(nbK[j] > 0)
				   muHat[j + l * this->kMaxS()] /= nbK[j];
		   }
	   }

		d_resStat = Rcpp::List::create( Rcpp::Named("k") = listK
				, Rcpp::Named("k_max") = this->kMaxS(), Rcpp::Named("it") = listIt
				, Rcpp::Named("nbState") = this->sizeState(), Rcpp::Named("mu") = mu
				, Rcpp::Named("muHat") = muHat
				, Rcpp::Named("nbK") = nbK);
	};

	Rcpp::List simRapport();
};

#include "SimulationNucleoD.tpp"
} /* namespace space_process */

#endif /* SIMULATIONNUCLEOD_H_ */
