/*
 * SegmentSeq.h
 *
 *  Created on: Aug 5, 2016
 *      Author: belleau
 */

#ifndef SEGMENTSEQ_H_
#define SEGMENTSEQ_H_
#include <Rcpp.h>

namespace space_process {

class SegmentSeq {
	const long d_sizeFReads, d_sizeRReads;

	double d_minPos, d_maxPos;
	const int d_zeta;

	std::vector<double> const &d_startFReads;
	std::vector<double> const &d_startRReads; /* vector of
												 reads start
												 position
												 ( foward and
												 reverse) */
public:
	SegmentSeq(std::vector<double> const &fReads,
			std::vector<double> const &rReads, int zeta);
	SegmentSeq(std::vector<double> const &fReads,
				std::vector<double> const &rReads, int zeta,
				long sizeFReads, long sizeRReads);
	virtual ~SegmentSeq();

	long sizeFReads();
	long sizeRReads();
	double minPos();
	double maxPos();
	int zeta();
private:
	void setMinMax();
};

} /* namespace space_process */

#endif /* SEGMENTSEQ_H_ */
