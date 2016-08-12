/*
 * PartitionAll.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef PARTITIONALL_H_
#define PARTITIONALL_H_
#include <iostream>
#include "SpaceNucleosomeD.h"
#include "NucleoDirichletPA.h"

namespace space_process {

template <typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
class PartitionAll: public SpaceNucleosomeD<NucleoD> {
public:

    typedef PartitionAll<NucleoD> NucleoSpace;

    typedef typename std::vector<double> containerD;
    typedef typename containerD::const_iterator iteratorD;
private:
    std::vector<double> *d_y;
    long d_ySize;
	//std::vector<double> &d_y;
public:

	PartitionAll(SegmentSeq const &segSeq)
		:SpaceNucleosomeD<NucleoD>(segSeq), d_y(new containerD){ //
        //d_y = new std::vector<double>;
        //d_y->push_back(1.0);

		d_y->insert(yBegin(), segSeq.beginFR(), segSeq.endFR());
		d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

		std::sort(d_y->begin(),d_y->end());
		d_ySize = (*d_y).size();
	};


	PartitionAll(SegmentSeq const &segSeq, int seed, containerD *y, long ySize)
		:SpaceNucleosomeD<NucleoD>(segSeq, seed), d_y(y){ //

	};

	PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, containerD *y, long ySize)
		:SpaceNucleosomeD<NucleoD>(segSeq, rng), d_y(y), d_ySize(ySize){ //

	};
	/*
	PartitionAll(std::vector<double> const  &fReads,
			std::vector<double> const &rReads, int zeta)
	:SpaceNucleosomeD<NucleoD>(fReads, rReads, zeta) {

		y = fReads;
	    y.insert(y.end(),rReads.begin(),rReads.end());
	    sort(y.begin(),y.end());
	};
	PartitionAll(std::vector<double> const  &fReads,
			std::vector<double> const &rReads,
			int zeta, long sizeFReads, long sizeRReads)
	:SpaceNucleosomeD<NucleoD>(fReads, rReads, zeta, sizeFReads, sizeRReads) {

		y = fReads;
	    y.insert(y.end(),rReads.begin(),rReads.end());
	    sort(y.begin(),y.end());
	};
	*/
	virtual ~PartitionAll(){};




	void initMu(int df){
		if(this->empty()){
			if(!(yEmpty()))
			{
				int cpt = 0;
				bool flag = true;
				do{
					double mu= gsl_ran_flat(this->rng(), this->minPos(), this->maxPos());
					NucleoD *u = new NucleoD(mu, df, this->segSeq(), this->rng());
					cpt++;
					(*u).setAvg(accumulate( yBegin(), yEnd(), 0.0)/ySize());

					//long t = setFoward(y[0], *u.avg(), *u);
					if(setFoward(y(0), (*u).avg(), *u) > 1)
					{
						iteratorD last = yEnd();
						/* end = (*(--last) + 1) to include the last read in sigmaR */
						if(setReverse((*u).avg(), (*(--last) + 1), *u) > 1)
						{
							flag = false;
							(*u).evalSigmaF();
							(*u).evalSigmaR();
							(*u).evalDelta();
							(*u).evalBF();
							(*u).evalBR();
							//w.push_back(1);
						}
					}

				}while(flag && cpt == 1000);
				if(flag){
					std::cerr << "Problem with the number of reads to initialise mu\n";
					exit(1);
				}
			}
			else{
				std::cerr << "No reads \n";
				exit(1);
			}
		}
	};

	long setFoward(double start, double end, NucleoD &u){
		long l = 0;
		iteratorD fStart = yBegin();
		iteratorD fEnd =  yEnd();
		int cpt = getLimit(start,end, fStart, fEnd, l);
		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};

	long setFoward(iteratorD fStart, iteratorD fEnd, double start, double end, NucleoD &u){
		long l = 0;
		int cpt = getLimit(start,end, fStart, fEnd, l);
		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};
	long setReverse(double start, double end, NucleoD &u){
		long l = 0;
		iteratorD rStart = yBegin();
		iteratorD rEnd =  yEnd();
		int cpt = getLimit(start,end, rStart, rEnd, l);
		u.setRStartPos(rStart, rEnd, cpt);
		return(l);
	};
	long setReverse(std::vector<double>::iterator rStart, std::vector<double>::iterator rEnd, double start, double end, NucleoD &u){
		long l = 0;
		int cpt = getLimit(start,end, rStart, rEnd, l);
		u.setRStartPos(rStart, rEnd, cpt);
		return(l);
	};

	NucleoSpace * clone(){

		NucleoSpace *a = new NucleoSpace(this->segSeq(), this->rng(), d_y, ySize());

		a->setValK(this->valK());
		a->setNucleosomes(this->nucleosomes());

		a->setMeanRead(this->meanRead());
		a->setR2(this->r2());
		a->setCMuDensity(this->cMuDensity());
        return(a);
	};

	void reset(){
		delete d_y;
	}

private:
	int getLimit(double start, double end, iteratorD &startIt, iteratorD &endIt, long &l){
		iteratorD it=yBegin();
		bool flag=1;
		int cpt = 0;
		double pr = -1.0;
		l = 0;
		while(flag && it!=yEnd()){
			if(*it >= start){
				startIt = it;

				while(*it < end  && it != yEnd()){
					if(pr < (*it + 0.0001)){
						l++;
					}
					pr = *it;
					it++;
					cpt++;
				}
				endIt = it;
				flag=0;
			}else{
				it++;
			}
		}

		return(cpt);
	};

	inline bool yEmpty(){
	    	return((*d_y).empty());
	};

    inline bool ySize(){
        	return(d_ySize);
	};

	inline iteratorD yBegin(){
		return((*d_y).begin());
	};

	inline iteratorD yEnd(){
		return((*d_y).end());
	};

	inline double y(int i){
		return((*d_y)[i]);
	};


};

} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
