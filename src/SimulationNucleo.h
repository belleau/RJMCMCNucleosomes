/*
 * SimulationNucleo.h
 *
 *  Created on: Sep 27, 2016
 *      Author: belleau
 */

#ifndef SIMULATIONNUCLEO_H_
#define SIMULATIONNUCLEO_H_

#include <Rcpp.h>
#include <gsl/gsl_randist.h>


namespace space_process {

template< typename NucleoSpace>
class SimulationNucleo{
	typedef std::vector<NucleoSpace *> NucleoSim;
	typedef typename NucleoSim::iterator itState;

	NucleoSim d_results;
	NucleoSpace *d_currentState, *d_mod;

	gsl_rng *d_rng;              // random number generator

	long d_nbIterations;         // Number of Iterations

	long d_kMaxS;                  // max of K in all simulation
	double d_rhoP1;
	SegmentSeq const &d_segSeq;


public:
	SimulationNucleo(SegmentSeq const &segSeq,
			gsl_rng * rng, long nbIteration=10000);

	virtual ~SimulationNucleo();



	int sizeState();
	void pushState();


	long nbIterations();

	long kMaxS();
	void setKMaxS(long kMaxS);

	double rhoP1();


	//virtual bool sampler();

protected:
	void setCurrentState(NucleoSpace *currentState);

	inline gsl_rng * rng(){
		return(d_rng);
	};

	SegmentSeq const &segSeq();

	void currentClone();

	inline NucleoSpace * currentState(){
		return(d_currentState);
	};

	inline NucleoSpace * mod(){
		return(d_mod);
	};

	void setRhoP1(double rhoP1);

	void acceptMod();

	itState beginState(){
		return(d_results.begin());
	}

	itState endState(){
		return(d_results.end());
	}

}; /* class SimulationNucleo */

#include "SimulationNucleo.tpp"

} /* namespace space_process */


#endif /* SIMULATIONNUCLEO_H_ */
