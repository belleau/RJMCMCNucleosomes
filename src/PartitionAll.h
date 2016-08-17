/*
 * PartitionAll.h
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#ifndef PARTITIONALL_H_
#define PARTITIONALL_H_
#include <iostream>
#include "SpaceNucleosomeD.h"
#include "NucleoDirichletPA.h"

namespace space_process {

template <typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
class PartitionAll: public SpaceNucleosomeD<NucleoD> {

    typedef PartitionAll<NucleoD> NucleoSpace;
    typedef std::list<NucleoD*> containerNucleo;
    typedef typename containerNucleo::iterator iteratorNucleo;

    typedef typename std::vector<double> containerD;
    typedef typename containerD::const_iterator iteratorD;

    std::vector<double> *d_y;
    long d_ySize;
	//std::vector<double> &d_y;
public:

	PartitionAll(SegmentSeq const &segSeq)
		:SpaceNucleosomeD<NucleoD>(segSeq), d_y(new containerD){ //
        //d_y = new std::vector<double>;
        //d_y->push_back(1.0);

		d_y->insert(yBegin(), segSeq.beginFR(), segSeq.endFR());
		d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

		std::sort(d_y->begin(),d_y->end());
		d_ySize = (*d_y).size();
	};

	PartitionAll(SegmentSeq const &segSeq, int seed)
		:SpaceNucleosomeD<NucleoD>(segSeq,seed), d_y(new containerD){ //
		//d_y = new std::vector<double>;
		//d_y->push_back(1.0);

		d_y->insert(yBegin(), segSeq.beginFR(), segSeq.endFR());
		d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

		std::sort(d_y->begin(),d_y->end());
		d_ySize = (*d_y).size();
	};

	PartitionAll(SegmentSeq const &segSeq, int seed, containerD *y, long ySize)
		:SpaceNucleosomeD<NucleoD>(segSeq, seed), d_y(y){ //

	};

	PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, containerD *y, long ySize)
		:SpaceNucleosomeD<NucleoD>(segSeq, rng), d_y(y), d_ySize(ySize){ //

	};
	/*
	PartitionAll(std::vector<double> const  &fReads,
			std::vector<double> const &rReads, int zeta)
	:SpaceNucleosomeD<NucleoD>(fReads, rReads, zeta) {

		y = fReads;
	    y.insert(y.end(),rReads.begin(),rReads.end());
	    sort(y.begin(),y.end());
	};
	PartitionAll(std::vector<double> const  &fReads,
			std::vector<double> const &rReads,
			int zeta, long sizeFReads, long sizeRReads)
	:SpaceNucleosomeD<NucleoD>(fReads, rReads, zeta, sizeFReads, sizeRReads) {

		y = fReads;
	    y.insert(y.end(),rReads.begin(),rReads.end());
	    sort(y.begin(),y.end());
	};
	*/
	virtual ~PartitionAll(){};




	bool initMu(int df){

		bool flag = true;
		NucleoD *u;
		if(this->empty()){
			if(!(yEmpty()))
			{
				int cpt = 0;

				do{
					double mu= gsl_ran_flat(this->rng(), this->minPos(), this->maxPos());

					u = new NucleoD(mu, df, this->segSeq(), this->rng());
					cpt++;
					double end = this->maxPos();
					flag = setNucleoD(u, y(0), end); //
                    if(flag) // Not enough read foward or reverse
                	{
                    	this->pushNucleo(u);
						this->evalW();
            		}
                    else{
                    	delete u;
                    }

				}while(!(flag) && cpt == 1000);
				if(!(flag)){
					std::cerr << "Problem with the number of reads to initialise mu\n";
					//exit(1);
				}
			}
			else{
				std::cerr << "No reads \n";
				//exit(1);
			}
		}

		return(flag);
	};

/*
	long setFoward(double start, double end, NucleoD &u){
		long l = 0;
		iteratorD fStart = yBegin();
		iteratorD fEnd =  yEnd();
		int cpt = getLimit(start,end, fStart, fEnd, l);

		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};*/


	long setFoward(iteratorD fStart, iteratorD fEnd, double start, double end, NucleoD &u){
		long l = 0;
		int cpt = getLimit(start,end, fStart, fEnd, l);

		u.setFStartPos(fStart, fEnd, cpt);
		return(l);
	};

	/*
	long setReverse(double start, double end, NucleoD &u){
		long l = 0;
		iteratorD rStart = yBegin();
		iteratorD rEnd =  yEnd();
		int cpt = getLimit(start,end, rStart, rEnd, l);

		u.setRStartPos(rStart, rEnd, cpt);
		return(l);
	};*/

	long setReverse(iteratorD rStart, iteratorD rEnd, double start, double end, NucleoD &u){

		long l = 0;
		int cpt = getLimit(start,end, rStart, rEnd, l);

		u.setRStartPos(rStart, rEnd, cpt);
		return(l);
	};



	NucleoSpace * clone(){

		NucleoSpace *a = new NucleoSpace(this->segSeq(), this->rng(), d_y, ySize());

		a->setValK(this->valK());
		a->setNucleosomes(this->nucleosomes());
		a->setLambda(this->lambda());
		a->setMeanRead(this->meanRead());
		a->setR2(this->r2());
		a->setCMuDensity(this->cMuDensity());

        return(a);
	};

	bool setNucleoD1(NucleoD *u, double aF, double aR){

		bool flag = false;
		long l;
		iteratorD startIt, endIt;
		int dimNucleo = getLimit( aF, aR, startIt, endIt, l);

		if(l > 1){ /* More than one distinct reads between aF and aR */

			(*u).setAvg(accumulate( startIt, endIt, 0.0)/dimNucleo);

			if(setFoward(startIt, endIt, aF, (*u).avg(), *u) > 1){

				if(setReverse(startIt, endIt, (*u).avg(), aR, *u) > 1){

					flag = true;
					(*u).evalSigmaF();
					(*u).evalSigmaR();
					(*u).evalDelta();
					(*u).evalBF();
					(*u).evalBR();
					(*u).setAF(aF);
					(*u).setAR(aR);
				}

			}
		}
		return(flag);
	}
	bool setNucleoD(NucleoD *u, double aF, double aR){

			bool flag = false;
			long l;
			iteratorD startIt, endIt;
			int dimNucleo = getLimit( aF, aR, startIt, endIt, l);

			if(l > 1){

				(*u).setAvg(accumulate( startIt, endIt, 0.0)/dimNucleo);

				if(setFoward(startIt, endIt, aF, (*u).avg(), *u) > 1){

					if(setReverse(startIt, endIt, (*u).avg(), aR, *u) > 1){

						(*u).evalSigmaF();
						(*u).evalSigmaR();
						if((*u).sigmaF() > 0.000001 && (*u).sigmaR() > 0.000001){
							(*u).evalDelta();
							(*u).evalBF();
							(*u).evalBR();
							(*u).setAF(aF);
							(*u).setAR(aR);
							flag = true;

						}
					}

				}
			}
			return(flag);
		}

	void birth(){
		int cpt = 0;
		bool flag = false;
		iteratorNucleo it1, it2;
		NucleoD *uBef, *uBirth, *uNext;
		double muBef, muBirth, muNext;
		int k = this->valK();
		int i = 0;
		try{
			do{
				uBef = NULL;
				uBirth = NULL;
				uNext = NULL;

				flag = false;
				cpt++;
				 /* gsl_ran_flat (const gsl_rng * r, double a, double b) */

				i = (int) gsl_ran_flat (this->rng(), 0, k+1); // nucleo between i-1 (minPos()) and i (maxPos())

				//std::cout << "I " << i << " k " << k <<  "\n";


				//long vMin, vMax;

				double startBirth;
				double aFBirth, aRBirth;



				if(i > 0){
					it1 = this->nucleosomes(this->nucleoBegin(), 0, i-1); // go to the position i-1

					muBef = (*it1)->mu();

					if(i < k){
						it2 = this->nucleosomes(it1, i-1, i);
						muNext = (*it2)->mu();
					}
					else{
						muNext = (*it1)->aR(); // maxPos
					}
				}
				else{  // i == 0
					it2 = this->nucleoBegin();
					muBef = (*it2)->aF();  // minPos
					muNext = (*it2)->mu();
				}


				muBirth = gsl_ran_flat(this->rng(), muBef, muNext);// jusqu'ici ok
				//std::cout << "Start " << this->minPos() << " End " << this->maxPos() << "\n";
				//std::cout << "Bef " << muBef << " Birth " << muBirth << " Next " << muNext << "\n";


				if(i > 0) // Modify nucleosome i-1
				{
					aFBirth = gsl_ran_flat(this->rng(), muBef, muBirth);

					uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1
					if(!(setNucleoD(uBef, (*it1)->aF(), aFBirth))){
						flag = true;
					}
					if(!(flag) && i < k){ // Modify nucleosome i+1

						aRBirth = gsl_ran_flat(this->rng(), muBirth, muNext);
						uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New i+1
						if(!(setNucleoD(uNext, aRBirth, (*it2)->aR()))){
							flag = true;
						}
					}
					else{
						aRBirth = this->maxPos() + 1;
					}
					if(!(flag)){

						int dF = (int) gsl_ran_flat (this->rng(), 3, 31);

						uBirth = new NucleoD(muBirth, dF, this->segSeq(), this->rng());

						if(!(setNucleoD(uBirth, aFBirth, aRBirth))){
							flag = true;
						}
					}
				}
				else{ // i == 0

					aFBirth = muBef;
					aRBirth = gsl_ran_flat(this->rng(), muBirth, muNext);
					uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New mu 0
					if(!(setNucleoD(uNext, aRBirth, (*it2)->aR()))){
						flag = true;
					}
					if(!(flag)){

						int df = (int) gsl_ran_flat (this->rng(), 3, 31);
						uBirth = new NucleoD(muBirth, df, this->segSeq(), this->rng());

						if(!(setNucleoD(uBirth, aFBirth, aRBirth))){
							flag = true;
						}
					}
				} // end i == 0
				if(flag){
					delete uBef;
					uBef = NULL;
					delete uNext;
					uNext = NULL;
					delete uBirth;
					uBirth = NULL;
				}

			}while(flag && cpt == 1000);

			if(!(flag))
			{

				if(i > 0){

					this->pushModNucleo(*it1);
					*it1 = uBef;
					this->pushAddNucleo(it1);

				}


				if(i < k){
					this->pushModNucleo(*it2);
					*it2 = uNext;
					this->insertNucleo(it2, uBirth);
					this->pushAddNucleo(it2);
					this->pushAddNucleo(--it2);
				}
				else{
					this->insertNucleo(this->nucleoEnd(), uBirth);
					this->pushAddNucleo(--(this->nucleoEnd()));
				}

				//delete uBirth;
			}
		}
		catch(std::bad_alloc&) {
		    std::cout << "Memory problem\n";
		    std::cerr << "Memory problem\n";
		}
	}

	void reset(){
		delete d_y;
		d_y = NULL;
		this->resetNucleo();
		/* delete les nucleosomes */

	}
	void accept(){
		this->resetMod();
	}
	void reject(){
		delete d_y;
		d_y = NULL;
		this->resetAdd();
	}

private:
	int getLimit(double start, double end, iteratorD &startIt, iteratorD &endIt, long &l){

		iteratorD it=yBegin();
		bool flag=1;
		int cpt = 0;
		double pr = -1.0;
		l = 0;
		while(flag && it!=yEnd()){
			if(*it >= start){
				startIt = it;

				while(*it < end  && it != yEnd()){
					if(pr < (*it + 0.000001)){
						l++;
					}
					pr = *it;
					it++;
					cpt++;
				}
				endIt = it;
				flag=0;
			}else{
				it++;
			}
		}

		return(cpt);
	};

	inline bool yEmpty(){
	    	return((*d_y).empty());
	};

    inline bool ySize(){
        	return(d_ySize);
	};

	inline iteratorD yBegin(){
		return((*d_y).begin());
	};

	inline iteratorD yEnd(){
		return((*d_y).end());
	};

	inline double y(int i){
		return((*d_y)[i]);
	};


};

} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
