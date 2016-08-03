/*
 * PartitionAll.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef PARTITIONALL_H_
#define PARTITIONALL_H_

#include "SpaceNucleosomeD.h"
#include "NucleoDirichletPA.h"

namespace space_process {

class PartitionAll: public SpaceNucleosomeD {
	std::vector<double> y;
public:
	PartitionAll(Rcpp::IntegerVector const  &fReads,
			Rcpp::IntegerVector const &rReads, int zeta);
	virtual ~PartitionAll();
	void initMu(double Mu, int df);
	void setFoward(double start, double end, NucleoDirichletPA &u);
	void setFoward(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, double start, double end, NucleoDirichletPA &u);
	void setReverse(double start, double end, NucleoDirichletPA &u);
	void setReverse(std::vector<double>::iterator rStart, std::vector<double>::iterator rEnd, double start, double end, NucleoDirichletPA &u);
private:
	int getLimit(double start, double end, std::vector<double>::iterator &startIt, std::vector<double>::iterator &endIt);

};

} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
