/*
 * NucleoDirichlet.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */
#include "NucleoDirichlet.h"
#include <iostream>
#include <gsl/gsl_randist.h>

namespace space_process {

NucleoDirichlet::NucleoDirichlet(double mu):
	Nucleosome(mu){
	// TODO Auto-generated constructor stub

}

NucleoDirichlet::~NucleoDirichlet() {
	// TODO Auto-generated destructor stub
}

double NucleoDirichlet::testT(){
	// TODO Auto-generated constructor stub
	double x = gsl_ran_tdist_pdf(0.4,3.0);
	std::cout << " T " << x << "\n";
	return(x);
}


} /* namespace space_process */
