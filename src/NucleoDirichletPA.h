/*
 * NucleoDirichletPA.h
 *
 *  Created on: Jul 29, 2016
 *      Author: belleau
 */

#ifndef NUCLEODIRICHLETPA_H_
#define NUCLEODIRICHLETPA_H_

#include "NucleoDirichlet.h"

namespace space_process {

class NucleoDirichletPA: public NucleoDirichlet {
    double d_aF, d_aR;
    double d_avg;
    public:
        NucleoDirichletPA(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng);
        NucleoDirichletPA(double mu, SegmentSeq const &segSeq, gsl_rng *rng);
        virtual ~NucleoDirichletPA();

        double aF();
        void setAF(double aF);
        double aR();
        void setAR(double aR);

        void setAvg(double avg);
        double avg();

        double testT();
        void testFRStart();

};

} /* namespace space_process */

#endif /* NUCLEODIRICHLETPA_H_ */
