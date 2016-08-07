/*
 * NucleoDirichletPA.cpp
 *
 *  Created on: Jul 29, 2016
 *      Author: belleau
 */

#include "NucleoDirichletPA.h"
#include <iostream>


using namespace std;

namespace space_process {

NucleoDirichletPA::NucleoDirichletPA(double mu, int df, SegmentSeq const &segSeq, gsl_rng *rng):
	NucleoDirichlet(mu, df, segSeq, rng){
	setSizeF(-1);
	setSizeR(-1);

}

NucleoDirichletPA::~NucleoDirichletPA() {
	// TODO Auto-generated destructor stub
}

double NucleoDirichletPA::testT(){
	// TODO Auto-generated constructor stub
	double x = gsl_ran_tdist_pdf(0.4,3.0);
	cout << " T " << x << "\n";
	return(x);
}


void NucleoDirichletPA::testFRStart()
{
	cout << "n " << sizeF() << "\n";
	int cpt = 0;
	/*for(vector<double>::iterator it = d_startF; it != d_endF; it++){
		cout << "F " << *it << " c " << ++cpt << "\n";
	}
	cout << "n " << d_sizeF << "\n";
	for(vector<double>::iterator it = d_startR; it != d_endR; it++){
			cout << "R " << *it << " c " << ++cpt << "\n";
	}*/

	cout << "sigmaF " << sigmaF() << "\n";
}


} /* namespace space_process */
