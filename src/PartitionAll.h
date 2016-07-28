/*
 * PartitionAll.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef PARTITIONALL_H_
#define PARTITIONALL_H_

#include "SpaceNucleosome.h"

namespace space_process {

class PartitionAll: public SpaceNucleosome {
	std::vector<double> y;
public:
	PartitionAll(Rcpp::IntegerVector const  &fReads,
			Rcpp::IntegerVector const &rReads, int zeta);
	virtual ~PartitionAll();
};

} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
