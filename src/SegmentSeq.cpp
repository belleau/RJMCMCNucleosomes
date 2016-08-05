/*
 * SegmentSeq.cpp
 *
 *  Created on: Aug 5, 2016
 *      Author: belleau
 */

#include "SegmentSeq.h"

namespace space_process {


SegmentSeq::SegmentSeq(std::vector<double> const &fReads,
			std::vector<double> const &rReads, int zeta)
	:d_startFReads(fReads), d_startRReads(rReads), d_zeta(zeta),
	 d_sizeFReads(fReads.size()), d_sizeRReads(rReads.size()){

	setMinMax();
}

SegmentSeq::SegmentSeq(std::vector<double> const &fReads,
			std::vector<double> const &rReads, int zeta,
			long sizeFReads, long sizeRReads)
	:d_startFReads(fReads), d_startRReads(rReads), d_zeta(zeta),
	 d_sizeFReads(sizeFReads), d_sizeRReads(sizeRReads){

	setMinMax();
}

SegmentSeq::~SegmentSeq() {
	// TODO Auto-generated destructor stub
}

void SegmentSeq::setMinMax(){
	d_minPos = std::min(*(std::min_element(d_startFReads.begin(),
					d_startFReads.end())),
			*(std::min_element(d_startRReads.begin(),
					d_startRReads.end())));

	d_maxPos = std::max(*(std::max_element(d_startFReads.begin(),
					d_startFReads.end())),
			*(std::max_element(d_startRReads.begin(),
					d_startRReads.end())));
}

long SegmentSeq::sizeFReads(){
	return(d_sizeFReads);
};

long SegmentSeq::sizeRReads(){
	return(d_sizeFReads);
};

double SegmentSeq::minPos(){
	return(d_minPos);
}

double SegmentSeq::maxPos(){
	return(d_maxPos);
}

int SegmentSeq::zeta(){
	return(d_zeta);
}

} /* namespace space_process */
