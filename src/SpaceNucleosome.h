/*
 * SpaceNucleosome.h
 *
 *  Created on: Aug 5, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOME_H_
#define SPACENUCLEOSOME_H_

#include "Nucleosome.h"
#include "SegmentSeq.h"

namespace space_process {

template<typename NucleoClass>    /***** BEWARE NucleoClass Must inherit from Nucleosome *****/
class SpaceNucleosome {

	SegmentSeq const d_segSeq;


	std::list<NucleoClass> d_nucleosomes;
public:
	SpaceNucleosome(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta)
		:d_segSeq(fReads, rReads, zeta){
	};

	SpaceNucleosome(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta, long sizeFReads, long sizeRReads)
		:d_segSeq(fReads, rReads, zeta, sizeFReads, sizeRReads){
	};
	virtual ~SpaceNucleosome(){};

	int size(){
		return(d_nucleosomes.size());
	};

	bool empty(){
		return(d_nucleosomes.empty());
	};

	void insert(NucleoClass &u){
		d_nucleosomes.push_back(u);
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
