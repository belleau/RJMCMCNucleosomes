/*
 * SpaceState.h
 *
 *  Created on: Jul 22, 2016
 *      Author: Pascal Belleau
 */



#ifndef SPACESTATE_H_
#define SPACESTATE_H_

#include <Rcpp.h>
#include <gsl/gsl_randist.h>

namespace space_process{

	class SpaceState {
	protected:
		long d_minPos, d_maxPos;
		gsl_rng *d_rng;  // random number generator
		Rcpp::IntegerVector const d_startFReads;
		Rcpp::IntegerVector const d_startRReads; /* vector of
																 reads start
																 position
																 ( foward and
																 reverse) */
		const int d_zeta;

	public:
		SpaceState(Rcpp::IntegerVector const  &fReads, Rcpp::IntegerVector const &rReads, int zeta);
		virtual ~SpaceState();

	};

}; /* namespace space_process */

#endif /* SPACESTATE_H_ */
