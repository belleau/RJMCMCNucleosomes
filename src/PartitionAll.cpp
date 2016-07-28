/*
 * PartitionAll.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#include "PartitionAll.h"
#include "NucleoDirichlet.h"

using namespace std;
using namespace space_process;
using namespace Rcpp;


namespace space_process {

PartitionAll::PartitionAll(Rcpp::IntegerVector const  &fReads,
		Rcpp::IntegerVector const &rReads, int zeta):SpaceNucleosome(fReads, rReads, zeta), SpaceState(fReads, rReads, zeta) {
	// TODO Auto-generated constructor stub
	y = as<vector<double> >(d_startFReads);
	std::vector<double> x = as<vector<double> >(d_startRReads);
	y.insert(y.end(),x.begin(),x.end());
	sort(y.begin(),y.end());

	/*for(vector<double>::iterator it=y.begin(); it!=y.end(); ++it){
	    	cout << "Val " << *it << "\n";
	}*/
	NucleoDirichlet bla(32);
	bla.testT();
}

PartitionAll::~PartitionAll() {
	// TODO Auto-generated destructor stub
}

} /* namespace space_process */
