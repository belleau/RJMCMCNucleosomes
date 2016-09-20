/*
 * SpaceNucleosome.h
 *
 *  Created on: Aug 5, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOME_H_
#define SPACENUCLEOSOME_H_

#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <list>
#include <math.h>
#include <time.h>

#include "Nucleosome.h"
#include "SegmentSeq.h"

namespace space_process {

template<typename NucleoClass>    /***** BEWARE NucleoClass Must inherit from Nucleosome *****/
class SpaceNucleosome {
    typedef std::list<NucleoClass*> containerNucleo;
    typedef typename containerNucleo::iterator itNucleo;
    typedef std::vector<NucleoClass*> vecNucleo;
	typedef typename vecNucleo::iterator itVNucleo;
	typedef std::vector<itNucleo> vecItNucleo;
	typedef typename vecItNucleo::iterator itItVNucleo;

    SegmentSeq const &d_segSeq;
	containerNucleo d_nucleosomes; /* List of nucleosomes */

	vecNucleo d_modNucleo;
	vecItNucleo d_addNucleo;

	int d_valK;
	int d_Max;
	gsl_rng *d_rng;  // random number generator
	long d_iteration;

public:


	SpaceNucleosome(SegmentSeq const &segSeq):
		d_segSeq(segSeq), d_valK(0), d_iteration(0){
		setRNG();
	};// Finir la cascade

	SpaceNucleosome(SegmentSeq const &segSeq, int seed):
			d_segSeq(segSeq), d_valK(0), d_iteration(0){
			setRNG(seed);
		};// Finir la cascade

	SpaceNucleosome(SegmentSeq const &segSeq, gsl_rng * rng):
		d_segSeq(segSeq), d_rng(rng), d_valK(0), d_iteration(0){

	}
	//const gsl_rng_type * T
/*	SpaceNucleosome(std::vector<double> const  &fReads,
			std::vector<double> const &rReads, int zeta)
		:d_segSeq(fReads, rReads, zeta), d_valK(0){
		setDefault();
	};

	SpaceNucleosome(std::vector<double> const  &fReads,
			std::vector<double> const &rReads, int zeta,
			long sizeFReads, long sizeRReads)
		:d_segSeq(fReads, rReads, zeta, sizeFReads, sizeRReads)
		, d_valK(0){

		setDefault();
	};*/


	virtual ~SpaceNucleosome(){};

	int size(){
		return(d_nucleosomes.size());
	};

	bool empty(){
		return(d_nucleosomes.empty());
	};

	void pushNucleo(NucleoClass *u){
		d_nucleosomes.push_back(u);
		d_valK++;
	};

	void insertNucleo(itNucleo it, NucleoClass *u){
		d_nucleosomes.insert(it, u);
		d_valK++;
	};
/*	void setDeltaMin(int deltaMin){
		d_segSeq.setDeltaMin(deltaMin);
	};
	void setDeltaMax(int deltaMax){
		d_segSeq.setDeltaMax(deltaMax);
	};*/

	void setRng(gsl_rng *rng){
		d_rng = rng;
	};

	int valK(){
		return(d_valK);
	};

	double minPos(){
		return(d_segSeq.minPos());
	};

	double maxPos(){
		return(d_segSeq.maxPos());
	};

	long sizeFReads(){
		return(d_segSeq.sizeFReads());
	}

	long sizeRReads(){
		return(d_segSeq.sizeRReads());
	}

	void displayMu(){
		std::cout << "Mu";
		for(itNucleo it = d_nucleosomes.begin() ; it != d_nucleosomes.end(); it++){
			std::cout << " " << (*it)->mu();
			std::cout << " : " << (*it)->avg();
		}
		std::cout << "\n";
	}

	Rcpp::NumericVector mu(){
		Rcpp::NumericVector mu = Rcpp::NumericVector(valK());
		int i = 0;
		for(itNucleo it = d_nucleosomes.begin() ; it != d_nucleosomes.end(); it++){
			mu[i++] = (*it)->mu();
		}
		return(mu);
	}

	void eraseNucleo(itNucleo it){ //itNucleo it
		d_nucleosomes.erase(it);
		d_valK--;
		//d_nucleosomes.erase(d_nucleosomes.begin());
	}

	void addIteration(){
		d_iteration++;
	};

	long iteration(){
		return(d_iteration);
	}

protected:
	gsl_rng * rng(){
		return(d_rng);
	};

	SegmentSeq const &segSeq(){
		return(d_segSeq);
	};

	itNucleo nucleoBegin(){
		return(d_nucleosomes.begin());
	};

	itNucleo nucleoEnd(){
		return(d_nucleosomes.end());
	};
	itNucleo nucleosomes(itNucleo itPos, int start, int pos){
		int i = start;
		do{
			itPos++;
		}while(pos > i++ && itPos != d_nucleosomes.end());
		if(pos < i)
		{
			itPos--;
		}
		return(itPos);
	};



	void pushModNucleo(NucleoClass *u){
		d_modNucleo.push_back(u);
	};

	void pushAddNucleo(itNucleo &u){
		d_addNucleo.push_back(u);
	};

	void resetNucleo(){
		for(itNucleo it = d_nucleosomes.begin(); it != d_nucleosomes.end();it++){
			if(*it != NULL){
				delete *it;
				*it = NULL;
			}
		}
		d_nucleosomes.clear();
	};

	void resetMod(){

		d_modNucleo.clear();
	};
	void displayMod(){

		std::cout << "Mod ";
		for(itVNucleo it = d_modNucleo.begin(); it != d_modNucleo.end();it++){
			std::cout << " M " << (*it)->mu();
		}
		std::cout << "\n";
		d_modNucleo.clear();
	};

	void resetAdd(){

		for(itItVNucleo it = d_addNucleo.begin(); it != d_addNucleo.end();it++){
			if(**it != NULL){
				delete **it;
				**it = NULL;
			}
		}
		d_addNucleo.clear();
	};

	void clearAdd(){
		d_addNucleo.clear();
	}

	containerNucleo &nucleosomes(){
		return(d_nucleosomes);
	};

	void setNucleosomes(containerNucleo &nucleosomes){

		d_nucleosomes = nucleosomes;
	};

	void setValK(int k){
		d_valK= k;
	}
private:

	void setRNG(){
		const gsl_rng_type * T;
		long seed;

		T = gsl_rng_default;

		d_rng = gsl_rng_alloc (T);     // pick random number generator
		seed = time (NULL) * getpid();
		gsl_rng_set (d_rng, seed);                  // set seed
	};

	void setRNG(int seed){
		const gsl_rng_type * T;

		T = gsl_rng_default;

		d_rng = gsl_rng_alloc (T);     // pick random number generator

		gsl_rng_set (d_rng, seed);                  // set seed
	};


/*	long sizeFReads(){
		return(d_sizeFReads);
	};
	long sizeRReads(){
		return(d_sizeRReads);
	};
protected:
	std::vector<double> const  startFReads(){
		return(d_startFReads);
	};

	std::vector<double> const  startRReads(){
		return(d_startRReads);
	};*/

}; /* Class SpaceNucleosome */

} /* namespace space_process */

#endif /* SPACENUCLEOSOME_H_ */
