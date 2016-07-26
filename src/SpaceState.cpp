/*
 * SpaceState.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: belleau
 */

#include <Rcpp.h>

#include "SpaceState.h"

using namespace std;
using namespace space_process;

SpaceState::SpaceState() {
	// TODO Auto-generated constructor stub
	p = 1;
}

SpaceState::~SpaceState() {
	// TODO Auto-generated destructor stub
}

int SpaceState::getP(){
	return(p);
}
