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
    int d_deltaMin, d_deltaMax;
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

    long sizeFReads() const;
    long sizeRReads() const;
    double minPos() const;
    double maxPos() const;
    int zeta() const;
    void setDeltaMin(int deltaMin);
    void setDeltaMax(int deltaMax);
    int deltaMin() const;
    int deltaMax() const;
    std::vector<double>::const_iterator beginFR() const;
    std::vector<double>::const_iterator endFR() const;
    std::vector<double>::const_iterator beginRR() const;
    std::vector<double>::const_iterator endRR() const;

private:
    void setMinMax();
    void setDefault();
};

} /* namespace space_process */

#endif /* SEGMENTSEQ_H_ */
