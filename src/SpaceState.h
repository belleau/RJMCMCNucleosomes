/*
 * SpaceState.h
 *
 *  Created on: Jul 22, 2016
 *      Author: Pascal Belleau
 */



#ifndef SPACESTATE_H_
#define SPACESTATE_H_

#include <Rcpp.h>

namespace space_process{

	class SpaceState {
	protected:
		long minPos, maxPos;
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
