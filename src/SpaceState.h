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
		int p;
		//Rcpp::IntegerVector startFReads, startRReads;

	public:
		SpaceState();//{ p = 1; }
		virtual ~SpaceState(); //{}
		int getP(); // {return(p);}
	};

}; /* namespace space_process */

#endif /* SPACESTATE_H_ */
