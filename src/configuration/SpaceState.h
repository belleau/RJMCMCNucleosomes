/*
 * SpaceState.h
 *
 *  Created on: Jul 22, 2016
 *      Author: belleau
 */

#ifndef SPACESTATE_H_
#define SPACESTATE_H_

namespace space_process{

class SpaceState {
	std::list<&Nucleosome> d_currentMu
public:
	SpaceState();
	virtual ~SpaceState();
};

}; /* namespace space_process */

#endif /* SPACESTATE_H_ */
