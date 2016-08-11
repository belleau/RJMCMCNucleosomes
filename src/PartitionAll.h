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
	std::vector<double> y;
public:
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
	virtual ~PartitionAll(){};

	void initMu(int df){
		if(this->empty()){
			if(!(y.empty()))
			{
				int cpt = 0;
				bool flag = true;
				do{
					double mu= gsl_ran_flat(this->rng(), this->minPos(), this->maxPos());
					NucleoD u(mu, df, this->segSeq(), this->rng());
					cpt++;
					u.setAvg(accumulate( y.begin(), y.end(), 0.0)/y.size());

					long t = setFoward(y[0], u.avg(), u);
					if(setFoward(y[0], u.avg(), u) > 1)
					{
						std::vector<double>::iterator last = y.end();
						/* end = (*(--last) + 1) to include the last read in sigmaR */
						if(setReverse(u.avg(), (*(--last) + 1), u) > 1)
						{
							flag = false;
							u.evalSigmaF();
							u.evalSigmaR();
							u.evalDelta();
							u.evalBF();
							u.evalBR();
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
		std::vector<double>::iterator fStart = y.begin();
		std::vector<double>::iterator fEnd =  y.end();
		int cpt = getLimit(start,end, fStart, fEnd, l);
		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};

	long setFoward(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, double start, double end, NucleoD &u){
		long l = 0;
		int cpt = getLimit(start,end, fStart, fEnd, l);
		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};
	long setReverse(double start, double end, NucleoD &u){
		long l = 0;
		std::vector<double>::iterator rStart = y.begin();
		std::vector<double>::iterator rEnd =  y.end();
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

private:
	int getLimit(double start, double end, std::vector<double>::iterator &startIt, std::vector<double>::iterator &endIt, long &l){
		std::vector<double>::iterator it=y.begin();
		bool flag=1;
		int cpt = 0;
		double pr = -1.0;
		l = 0;
		while(flag && it!=y.end()){
			if(*it >= start){
				startIt = it;

				while(*it < end  && it != y.end()){
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

};

} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
