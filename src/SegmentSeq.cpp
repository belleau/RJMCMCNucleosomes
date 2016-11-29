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
    :d_sizeFReads(fReads.size()), d_sizeRReads(rReads.size()),
        d_zeta(zeta), d_startFReads(fReads), d_startRReads(rReads)
     {

    setDefault();
}

SegmentSeq::SegmentSeq(std::vector<double> const &fReads,
            std::vector<double> const &rReads, int zeta,
            long sizeFReads, long sizeRReads)
   :d_sizeFReads(sizeFReads), d_sizeRReads(sizeRReads),
        d_zeta(zeta), d_startFReads(fReads), d_startRReads(rReads){

    setDefault();
}

SegmentSeq::~SegmentSeq() {
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

void SegmentSeq::setDefault(){
    setMinMax();
    d_deltaMin = d_zeta - 5;
    d_deltaMax = d_zeta + 5;
}

long SegmentSeq::sizeFReads() const{
    return(d_sizeFReads);
};

long SegmentSeq::sizeRReads() const{
    return(d_sizeRReads);
};

double SegmentSeq::minPos() const{
    return(d_minPos);
}

double SegmentSeq::maxPos() const{
    return(d_maxPos);
}

int SegmentSeq::zeta() const{
    return(d_zeta);
}

void SegmentSeq::setDeltaMin(int deltaMin){
    d_deltaMin = deltaMin;
}
void SegmentSeq:: setDeltaMax(int deltaMax){
    d_deltaMax = deltaMax;
}
int SegmentSeq::deltaMin() const{
    return(d_deltaMax);
}
int SegmentSeq::deltaMax() const{
    return(d_deltaMin);
}
std::vector<double>::const_iterator SegmentSeq::beginFR() const{
    return(d_startFReads.begin());
}

std::vector<double>::const_iterator SegmentSeq::endFR() const{
    return(d_startFReads.end());
}

std::vector<double>::const_iterator SegmentSeq::beginRR() const{
    return(d_startRReads.begin());
}

std::vector<double>::const_iterator SegmentSeq::endRR() const{
    return(d_startRReads.end());
}


} /* namespace space_process */
