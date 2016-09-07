#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <math.h>

#include "PartitionAll.h"

//#include "NucleoDirichletPA.h"
//#include "Factory.h"
//#include "bla2.h"
//#include "bla1.h"

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
                        bool adaptIterationsToReads = true) {
    NumericVector startFReads(startPosForwardReads); // *startFReads = new IntegerVector(startPosForwardReads);
    NumericVector startRReads(startPosReverseReads); // *startRReads = new IntegerVector(startPosReverseReads);
    std::vector<double> fReads = Rcpp::as<std::vector<double> >(startFReads);
    std::vector<double> rReads = Rcpp::as<std::vector<double> >(startRReads);

    int nf=1, nr;
    long tot;

    const gsl_rng_type * T;
    gsl_rng *rng;
	long seed;

	T = gsl_rng_default;

	rng = gsl_rng_alloc (T);     // pick random number generator
	seed = kMax; //time (NULL) * getpid();
	gsl_rng_set (rng, seed);




	std::vector< PartitionAll<NucleoDirichletPA> *> res;

    SegmentSeq seg(fReads, rReads, 147);
    PartitionAll<NucleoDirichletPA> *currentState = new PartitionAll<NucleoDirichletPA>(seg, rng);
    (*currentState).initMu1( 3);// test si ok
    (*currentState).prepSpace();
    (*currentState).addIteration();
   res.push_back(currentState);

   bool dispRho = adaptIterationsToReads;
   double bla = 0;
   int cptBla = 0;
   int test3 = 0;
   double nbType[5] = {0,0,0,0,0};
   bool vStop = true;
   int kMaxO = (*currentState).valK();
   //int *k = new int[nbrIterations];
   //double **mu = new double *[nbrIterations];

   for(long i = 0; i< nbrIterations;i++){
	   double valPourv = 0;
	   double pt;
	   bool flag = false;

	  /*if(i%10000 == 0)
	   {
		   cout << i << "\n";
		   (*currentState).displayMu();
		   //if(dispRho)
		   cout << "K " << (*currentState).valK() << "\n";
	   }*/


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
						valPourv = (*mod).rhoP2();
						rho *= (*mod).qalloc();
						//bla += (*mod).kD();
						//bla += (*mod).tB();
						//cptBla++;
						//cout << "v0 " << (*mod).kD() << " " << (*currentState).kD() << "\n";
						if(dispRho)
						cout << "Birth " << rho << "\n";
						test3 = 4;

					}
					else{
						cout << "Aye\n";
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
						test3 = 5;
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
					test3 = 3;
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
					test3 = 1;
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
					test3 = 2;
				}
    		}
    	}

    	if(flag){
    		nbType[test3-1]++;
    		if(dispRho)
    		cout << "v0 " << (*mod).kD() << " " << (*currentState).kD() << "\n";
    		//valPourv = rho;
    		double tmp = exp(((*mod).kD() - (*currentState).kD()));
    		//cout << " v1 " << tmp;
    		rho *= tmp;
    		tmp = ((*mod).priorMuDensity() / (*currentState).priorMuDensity());
    		//cout << " v2 " << tmp;
			rho *= tmp;


    		tmp = ((*mod).multinomial() / (*currentState).multinomial());
    		//cout << " v3 " << tmp;
			rho *= tmp;

    		rho = std::min(1.0, rho);
    		//cout << " Rho " << rho << "\n";
			pt = gsl_ran_flat (rng, 0, 1);
			if(rho >= pt){
				/*std::cout << "type " << test3 << " i " << (*mod).tB();
								cout << " k " << (*mod).valK() << " kd " << (*mod).kD();
								cout << " " << nbType[0] << " " << nbType[1] << " " << nbType[2] << " " << nbType[3] << " " << nbType[4];
								cout << " rho " << rho << " pt " << pt;
								cout << "\n";*/
					// clean currentState
					/*if(vStop){*/
					//(*currentState).delCurrent();
					//delete currentState;
					//(*mod).accept();
					currentState = mod;
					(*currentState).addIteration();
					res.push_back(currentState);
					kMaxO = max(kMaxO, (*currentState).valK());
					//(*currentState).displayMu();

				/*}
				if((*currentState).valK() == 3)
					vStop = false;*/
				/*std::cout << "type " << test3 << " i " << (*currentState).tB();
				cout << " k " << (*currentState).valK() << " kd " << (*currentState).kD();
				cout << " " << nbType[0] << " " << nbType[1] << " " << nbType[2] << " " << nbType[3] << " " << nbType[4];
				cout << "\n";*/
				//fill(&nbType[0], &nbType[5], 0);

				/*if(test3){
					if((*currentState).valK() == 4){
						std::cout << "test3 " << (i+1) << "\n";
						//(*currentState).displayMu();
						test3 = false;
					}
				}*/
			}
			else{
				/*if(test3 == 4){
					std::cout << "Rtype " << test3 << " i " << (*mod).tB();
					cout << " k " << (*mod).valK() << " kd " << (*mod).kD();
					cout << " kdC " << (*currentState).kD();
					cout << " rho " << rho << " pt " << pt << " pourv " << valPourv;
					//cout << " " << nbType[0] << " " << nbType[1] << " " << nbType[2] << " " << nbType[3] << " " << nbType[4];
					cout << "\n";
				}*/
				(*currentState).addIteration();
				(*mod).reject();
				delete mod;
			}
    	}
    	else{
    		(*currentState).addIteration();
			delete mod;
    	}

    }

    /*cout << "\n";
   	(*currentState).displayMu();

    cout << "K " << (*currentState).valK() << "\n";
    //cout << "Ok " << bla << " cpt " << cptBla << " m " << bla / cptBla << "\n";*/
    tot = res.size();
    //cout << "Res " << tot << "\n";

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

	   Rcpp::NumericVector tmp = (*it)->mu();

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
	   //cout << " " <<  (*it)->iteration() << ":" << (*it)->valK();

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

    //nf = startFReads.size();

    //nf = currentState.getP();
    //nr = startRReads.size();


    List nbSeq = List::create( Rcpp::Named("K") = listK
    				, Rcpp::Named("KMax") = kMaxO, Rcpp::Named("it") = listIt
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

    //(*currentState).delCurrent();
    //delete currentState;
    return nbSeq;
}
