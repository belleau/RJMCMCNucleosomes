#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "PartitionAll.h"
#include "SimulationNucleoD.h"


#include "SegmentSeq.h"

//typedef space_process::SpaceState regionState;

using namespace Rcpp;
using namespace std;
using namespace space_process;

// Below is a simple example of exporting a C++ function to R. You can
//

// [[Rcpp::export]]
List rjmcmcNucleo(SEXP startPosForwardReads, SEXP startPosReverseReads,
                        long nbrIterations, int kMax, int lambda,
                        int minInterval, int maxInterval, int minReads = 5,
                        bool adaptIterationsToReads = true, int vSeed = -1) {
	/*********************************************************************
	 * init R var
	 *********************************************************************/
	NumericVector startFReads(startPosForwardReads); // *startFReads = new IntegerVector(startPosForwardReads);
    NumericVector startRReads(startPosReverseReads); // *startRReads = new IntegerVector(startPosReverseReads);

    std::vector<double> fReads = Rcpp::as<std::vector<double> >(startFReads);
    std::vector<double> rReads = Rcpp::as<std::vector<double> >(startRReads);

    long nr = 0; // Number of read
    long tot = 0;
    int kMaxO = 0; // the biggest K
    int mO = 1;
    //kMax = 30;

    /*********************************************************************
	 * Debug var
	 *********************************************************************/

    bool dispRho = false;           // display rho for debug
    bool displayInt = false;
    int typeMv = 0;                 // move type for debug
    double nbType[5] = {0,0,0,0,0}; // number of type move for debug


    /*********************************************************************
	 * init random number generator
	 *********************************************************************/
    const gsl_rng_type * T;
    gsl_rng *rng;
	//long seed;

	T = gsl_rng_default;

	rng = gsl_rng_alloc (T);     // pick random number generator
	if(vSeed <= 0){
		vSeed = time (NULL) * getpid();
	}
	gsl_rng_set (rng, vSeed);

	/*********************************************************************
	 * Space Nucleosome and Segment
	 *********************************************************************/

	std::vector< PartitionAll<NucleoDirichletPA> *> res; // vector of space nucleosome accepted

	SegmentSeq seg(fReads, rReads, 147);                 // Reads of the segment



	nr = seg.sizeFReads() + seg.sizeRReads();
    if(adaptIterationsToReads){
    	if(nr <= 12){
    		nbrIterations = 1000;
    	}
    }

    List resO;

    if(mO > 1){

		SimulationNucleoD<PartitionAll<NucleoDirichletPA> > pourv(seg, rng, kMax, nbrIterations);
		if(pourv.initMu()){
			pourv.simulate();
			pourv.statSim();
			resO = pourv.simRapport();
		}
		else{
			resO = R_NilValue;
		}
    }
    else{
    /*********************************************************************
	 * Init Space Nucleosome
	 *********************************************************************/

		PartitionAll<NucleoDirichletPA> *currentState = new PartitionAll<NucleoDirichletPA>(seg, rng, kMax);
		(*currentState).initMu1( 3);// test si ok
		(*currentState).prepSpace();
		(*currentState).addIteration();
		res.push_back(currentState);



		for(long i = 0; i< nbrIterations;i++){

			double pt;
			bool flag = false; // Generate a new valide (with enought read in the partition) move

			if(displayInt && i%10000 == 0){
				cout << i << "\n";
				(*currentState).displayMu();
				cout << "K " << (*currentState).valK() << "\n";
			}


			PartitionAll<NucleoDirichletPA> *mod = (*currentState).clone();
			double rho = 1;
			pt = gsl_ran_flat (rng, 0, 1);
			if((*currentState).valK() > 1){

				//std::cout << "dK " << (*currentState).dK() << " BK " << (*currentState).bK() << " u " << pt << "\n";
				if(pt > (*currentState).dK()){
					if(pt <= ((*currentState).dK() + (*currentState).bK()) ){
						flag = (*mod).birthR();
						if(flag){
							(*mod).prepSpace();
							rho = (*mod).rhoP2() / (*currentState).bK();
							rho *= (*mod).qalloc();

							if(dispRho)
							cout << "Birth " << rho << "\n";
							typeMv = 4;
						}
					}
					else{
						flag = (*mod).mhR();
						if(flag){
							(*mod).prepSpace();
							rho = 1.0;
							//rho = (*mod).rhoP2() / (*currentState).bK();

							if(dispRho)
							cout << "MH " << rho << "\n";
							typeMv = 5;
						}
					}
				}
				else{

					flag = (*mod).death();
					if(flag){
						(*mod).prepSpace();
						rho = (*mod).bK() / (*currentState).rhoP2();
						rho /= (*mod).qalloc();
						if(dispRho)
						cout << "Death " << rho << "\n";
						typeMv = 3;
					}
				}
			}
			else
			{
				if(pt <= 0.5){
					flag = (*mod).birthR();
					if(flag){
						(*mod).prepSpace();
						rho = (*mod).rhoP2() / (*currentState).bK(); // << (*mod).rhoP2()
						rho *= (*mod).qalloc();
						//cout << "v0 1 " << (*mod).kD() << " " << (*currentState).kD() << "\n";
						if(dispRho)
						cout << "Birth1 "  << rho << "\n";
						typeMv = 1;
					}
				}
				else{
					flag = (*mod).mhR();
					if(flag){
						(*mod).prepSpace();
						rho = 1.0;
						//rho = (*mod).rhoP2() / (*currentState).bK();

						if(dispRho)
						cout << "MH1 " << rho << "\n";
						typeMv = 2;
					}
				}
			}

			if(flag){
				nbType[typeMv-1]++;
				if(dispRho)
				cout << "v0 " << (*mod).kD() << " " << (*currentState).kD() << "\n";

				double tmp = exp(((*mod).kD() - (*currentState).kD()));

				rho *= tmp;
				tmp = ((*mod).priorMuDensity() / (*currentState).priorMuDensity());

				rho *= tmp;


				tmp = ((*mod).multinomial() / (*currentState).multinomial());

				rho *= tmp;

				rho = std::min(1.0, rho);

				pt = gsl_ran_flat (rng, 0, 1);

				if(rho >= pt){
						currentState = mod;
						(*currentState).addIteration();
						res.push_back(currentState);
						kMaxO = max(kMaxO, (*currentState).valK());
				}
				else{
					(*currentState).addIteration();
					(*mod).reject();
					delete mod;
				}
			}
			else{
				(*currentState).addIteration();
				delete mod;
			}

		} // End for iteration


		tot = res.size();


		Rcpp::NumericVector listK = Rcpp::NumericVector( Rcpp::Dimension(tot));
		Rcpp::NumericMatrix mu = Rcpp::NumericMatrix( Rcpp::Dimension(tot, kMaxO));
		Rcpp::IntegerVector listIt = Rcpp::IntegerVector( Rcpp::Dimension(tot));
		Rcpp::IntegerVector nbK = Rcpp::IntegerVector(kMaxO);
		Rcpp::NumericVector muHat = Rcpp::NumericVector(Rcpp::Dimension(kMaxO, kMaxO));

	   int i = 0;
	   for(std::vector< PartitionAll<NucleoDirichletPA> *>::iterator it = res.begin(); it != res.end();it++){

		   listK[i] = (*it)->valK();
		   listIt[i] = (*it)->iteration();
		   nbK[(*it)->valK()-1]+= (*it)->iteration();

		   //Rcpp::NumericVector tmp = (*it)->mu();
		   std::vector<double> tmp = (*it)->mu();

		   for(int j = 0; j < kMaxO; j++){
			   if(j < listK[i]){
				   mu[i  + j * tot] = tmp[j];
				   muHat[((*it)->valK()-1) + j * kMaxO] += (*it)->iteration() * tmp[j];

			   }
			   else{
				   mu[i + j * tot] = 0;
			   }
		   }
		   i++;
	   }


	   int posM = 0;
	   int curM = -1;
	   for(int j = 0; j < kMaxO; j++){
		   if(nbK[j] > curM){
			   posM = j + 1;
		   }
		   for(int l = 0; l < kMaxO; l++){
			   if(nbK[j] > 0)
				   muHat[j + l * kMaxO] /= nbK[j];
		   }
	   }



	   resO = List::create( Rcpp::Named("k") = listK
						, Rcpp::Named("k_max") = kMaxO, Rcpp::Named("it") = listIt
						, Rcpp::Named("tot") = tot, Rcpp::Named("mu") = mu
						, Rcpp::Named("muHat") = muHat
						, Rcpp::Named("nbK") = nbK); //(*currentState).mu()

		for(std::vector< PartitionAll<NucleoDirichletPA> *>::iterator it = res.begin(); it != res.end();it++){
			if(it == res.begin()){
				(*it)->delCurrent();
			}
			else{
				(*it)->delMod();
			}
			delete (*it);
		}
    }
    //return nbSeq;
    return resO;
}
