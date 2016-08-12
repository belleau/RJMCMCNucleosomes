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

#include "Nucleosome.h"
#include "SegmentSeq.h"

namespace space_process {

template<typename NucleoClass>    /***** BEWARE NucleoClass Must inherit from Nucleosome *****/
class SpaceNucleosome {
    typedef std::list<NucleoClass*> containerNucleo;
    typedef typename containerNucleo::const_iterator itNucleo;

	//SegmentSeq const d_segSeq; /* sera bientot une reference quand plusieurs SpaceNucleosome */
    SegmentSeq const &d_segSeq;
	containerNucleo d_nucleosomes; /* List of nucleosomes */
	int d_valK;
	gsl_rng *d_rng;  // random number generator

public:


	SpaceNucleosome(SegmentSeq const &segSeq):
		d_segSeq(segSeq){
		setRNG();
	};// Finir la cascade

	SpaceNucleosome(SegmentSeq const &segSeq, int seed):
			d_segSeq(segSeq){
			setRNG(seed);
		};// Finir la cascade

	SpaceNucleosome(SegmentSeq const &segSeq, gsl_rng * rng):
		d_segSeq(segSeq), d_rng(rng){

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


	void insert(NucleoClass *u){
		d_nucleosomes.push_back(u);
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
