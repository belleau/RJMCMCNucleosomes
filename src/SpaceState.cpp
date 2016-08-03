/*
 * SpaceState.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: belleau
 */

#include <Rcpp.h>
#include <iostream>
#include "SpaceState.h"
#include <gsl/gsl_randist.h>
#include <unistd.h>

using namespace std;
using namespace space_process;
using namespace Rcpp;

/* les valeurs des startPos des reads passer peuvent etre modifier mais sont dans
 * des constantes qui ne le peuvent pas
 */
SpaceState::SpaceState(IntegerVector const &fReads, IntegerVector const &rReads, int zeta)
:d_startFReads(fReads), d_startRReads(rReads), d_zeta(zeta)
{

	d_minPos = min(long(min(d_startFReads)), long(min(d_startRReads)));
	d_maxPos = max(long(max(d_startFReads)), long(max(d_startRReads)));

	const gsl_rng_type * T;
	long seed;

	T = gsl_rng_default;

	d_rng = gsl_rng_alloc (T);     // pick random number generator
	seed = time (NULL) * getpid();
	gsl_rng_set (d_rng, seed);                  // set seed

/*
	cout << "zeta " << d_zeta << "\n";
	cout <<"F "<< min(d_startFReads) << "\n";
	cout <<"R "<< d_startRReads.size() << "\n";
	*/
}

SpaceState::~SpaceState() {
	// TODO Auto-generated destructor stub
}


