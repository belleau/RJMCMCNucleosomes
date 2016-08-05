/*
 * PartitionAll.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#include "PartitionAll.h"
#include "NucleoDirichletPA.h"

using namespace std;
using namespace space_process;
using namespace Rcpp;


namespace space_process {
/*
PartitionAll::PartitionAll(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta):SpaceNucleosomeD(fReads, rReads, zeta) {

	y = fReads;
    y.insert(y.end(),rReads.begin(),rReads.end());
    sort(y.begin(),y.end());
	int cpt = 0;
	for(vector<double>::iterator it=y.begin(); it!=y.end(); ++it){
	    cpt++;
		cout << "Val " << *it << " cpt " << cpt << "\n";
	}
	cout <<"\n\n";

}


PartitionAll::PartitionAll(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta, long sizeFReads, long sizeRReads):SpaceNucleosomeD(fReads, rReads, zeta, sizeFReads, sizeRReads) {

	y = fReads;
    y.insert(y.end(),rReads.begin(),rReads.end());
    sort(y.begin(),y.end());
}

int PartitionAll::getLimit(double start, double end, vector<double>::iterator &startIt, vector<double>::iterator &endIt){
	vector<double>::iterator it=y.begin();
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
}

void PartitionAll::setFoward(double start, double end, NucleoDirichletPA &u){

	vector<double>::iterator fStart = y.begin();
	vector<double>::iterator fEnd =  y.end();
	int cpt = getLimit(start,end, fStart, fEnd);
	u.setFStartPos(fStart, fEnd, cpt);
}

void PartitionAll::setFoward(vector<double>::iterator fStart, vector<double>::iterator fEnd, double start, double end, NucleoDirichletPA &u){

	int cpt = getLimit(start,end, fStart, fEnd);
	u.setFStartPos(fStart, fEnd, cpt);
}

void PartitionAll::setReverse(double start, double end, NucleoDirichletPA &u){

	vector<double>::iterator rStart = y.begin();
	vector<double>::iterator rEnd =  y.end();
	int cpt = getLimit(start,end, rStart, rEnd);
	u.setRStartPos(rStart, rEnd, cpt);
}

void PartitionAll::setReverse(vector<double>::iterator rStart, vector<double>::iterator rEnd, double start, double end, NucleoDirichletPA &u){

	int cpt = getLimit(start,end, rStart, rEnd);
	u.setRStartPos(rStart, rEnd, cpt);
}


void PartitionAll::initMu(double Mu, int df){
	if(empty()){
		if(!(y.empty()))
		{
			NucleoDirichletPA u(Mu, df, d_startFReads, d_startRReads, d_sizeFReads, d_sizeRReads);

			u.setAvg(accumulate( y.begin(), y.end(), 0.0)/y.size());

			setFoward(y[0], u.avg(), u);
			vector<double>::iterator last = y.end();
			setReverse(u.avg(), *(--last), u);
			u.evalSigmaF();
			u.evalSigmaR();
		}
		else{
			cerr << "No reads \n";
			exit(1);
		}
	}
}


PartitionAll::~PartitionAll() {
	// TODO Auto-generated destructor stub
}
*/
} /* namespace space_process */
