/*
 * SpaceState.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: belleau
 */

#include <Rcpp.h>
#include <iostream>
#include "SpaceState.h"

using namespace std;
using namespace space_process;
using namespace Rcpp;

/* les valeurs des startPos des reads passer peuvent etre modifier mais sont dans
 * des constantes qui ne le peuvent pas
 */
SpaceState::SpaceState(IntegerVector const &fReads, IntegerVector const &rReads, int zeta)
:d_startFReads(fReads), d_startRReads(rReads), d_zeta(zeta)
{

	minPos = min(long(min(d_startFReads)), long(min(d_startRReads)));
	maxPos = max(long(max(d_startFReads)), long(max(d_startRReads)));

	cout << "zeta " << d_zeta << "\n";
	cout <<"F "<< min(d_startFReads) << "\n";
	cout <<"R "<< d_startRReads.size() << "\n";
}

SpaceState::~SpaceState() {
	// TODO Auto-generated destructor stub
}


