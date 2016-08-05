/*
 * PartitionAll.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef PARTITIONALL_H_
#define PARTITIONALL_H_

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
	/*	int cpt = 0;
		for(vector<double>::iterator it=y.begin(); it!=y.end(); ++it){
		    cpt++;
			cout << "Val " << *it << " cpt " << cpt << "\n";
		}
		cout <<"\n\n";
	*/
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

	void initMu(double Mu, int df){
		if(this->empty()){
			if(!(y.empty()))
			{
				NucleoD u(Mu, df, this->startFReads(), this->startRReads(), this->sizeFReads(), this->sizeRReads());

				u.setAvg(accumulate( y.begin(), y.end(), 0.0)/y.size());

				setFoward(y[0], u.avg(), u);
				std::vector<double>::iterator last = y.end();
				setReverse(u.avg(), *(--last), u);
				u.evalSigmaF();
				u.evalSigmaR();
			}
			else{
				std::cerr << "No reads \n";
				exit(1);
			}
		}
	};
	void setFoward(double start, double end, NucleoD &u){

		std::vector<double>::iterator fStart = y.begin();
		std::vector<double>::iterator fEnd =  y.end();
		int cpt = getLimit(start,end, fStart, fEnd);
		u.setFStartPos(fStart, fEnd, cpt);
	};
	void setFoward(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, double start, double end, NucleoD &u){

		int cpt = getLimit(start,end, fStart, fEnd);
		u.setFStartPos(fStart, fEnd, cpt);
	};
	void setReverse(double start, double end, NucleoD &u){

		std::vector<double>::iterator rStart = y.begin();
		std::vector<double>::iterator rEnd =  y.end();
		int cpt = getLimit(start,end, rStart, rEnd);
		u.setRStartPos(rStart, rEnd, cpt);
	};
	void setReverse(std::vector<double>::iterator rStart, std::vector<double>::iterator rEnd, double start, double end, NucleoD &u){

		int cpt = getLimit(start,end, rStart, rEnd);
		u.setRStartPos(rStart, rEnd, cpt);
	};

private:
	int getLimit(double start, double end, std::vector<double>::iterator &startIt, std::vector<double>::iterator &endIt){
		std::vector<double>::iterator it=y.begin();
		bool flag=1;
		int cpt = 0;
		while(flag && it!=y.end()){
			if(*it >= start){
				startIt = it;
				while(*it <= end  && it != y.end()){
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
