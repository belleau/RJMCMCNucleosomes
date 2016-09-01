#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <math.h>

//#include "SpaceState.h"
#include "SpaceNucleosomeD.h"
#include "PartitionAll.h"

#include "NucleoDirichletPA.h"
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
    //startFReads[1] = 1;


    /*space_process::bla2<int> a;
    a.insert(1);
    a.insert(2);
    a.insert(3);
    cout << "Size a " << a.size() << "\n";

    space_process::bla2<int> b(a);
    b.insert(5);

    cout << "Size a " << a.size() << "\n";
    cout << "Size b " << b.size() << "\n";
    space_process::bla2<int> *c = a.clone();

    (*c).insert(6);
    (*c).insert(6);
    (*c).insert(6);

    cout << "Size a " << a.size() << "\n";
    cout << "Size b " << b.size() << "\n";
    cout << "Size c " << (*c).size() << "\n";*/

    //space_process::SpaceState currentState(startFReads, startRReads, 147);
    //space_process::SpaceNucleosome currentState(startFReads, startRReads, 147);
    //space_process::PartitionAll bla(startFReads, startRReads, 147);

    //space_process::Factory<space_process::SpaceNucleosomeD, space_process::PartitionAll> truc;
    // space_process::NucleoDirichletPA

    //space_process::PartitionAll<space_process::NucleoDirichletPA> bla(fReads, rReads, 147);
    //cout << "Bla " << fReads[0] << "\n";

    /*space_process::PartitionAll<space_process::NucleoDirichletPA> currentState(fReads, rReads, 147);

    currentState.insertD(10010,3);
    currentState.insertD(10020,3);
    currentState.insertD(10030,3);
    currentState.insertD(10040,3);
    currentState.insertD(10050,3);
    currentState.insertD(10060,3);
    currentState.evalPriorMuDensity();*/
/*
    bla1<int> o;
    o.pourv();
*/
    const gsl_rng_type * T;
    gsl_rng *rng;
	long seed;

	T = gsl_rng_default;

	rng = gsl_rng_alloc (T);     // pick random number generator
	seed = kMax; //time (NULL) * getpid();
	gsl_rng_set (rng, seed);



	/*PartitionAll<NucleoDirichletPA> *mod = (*currentState).clone();
	    (*mod).birth();
		(*mod).accept();
		(*mod).birth();
		(*mod).accept();
		(*mod).birth();
		(*mod).accept();
		(*mod).birth();
		(*mod).accept();
		(*mod).birth();
		(*mod).accept();
		(*mod).birth();
		(*mod).accept();
		(*mod).birth();
		//(*mod).accept();
		delete currentState;
		(*mod).accept();
		currentState = mod;
		mod = (*currentState).clone();
		cout << "New\n";
		(*mod). displayMu();
		(*mod).death();
		(*mod).accept();
		(*mod).reset();*/


    SegmentSeq seg(fReads, rReads, 147);
    PartitionAll<NucleoDirichletPA> *currentState = new PartitionAll<NucleoDirichletPA>(seg, rng);
    (*currentState).initMu1( 3);// test si ok
    (*currentState).prepSpace();



   bool dispRho = adaptIterationsToReads;
   double bla = 0;
   int cptBla = 0;
   int test3 = 0;
   double nbType[5] = {0,0,0,0,0};
   bool vStop = true;
   for(long i = 0; i< nbrIterations;i++){
	   double valPourv = 0;
	   double pt;
	   bool flag = false;

	  if(i%10000 == 0)
	   {
		   cout << i << "\n";
		   (*currentState).displayMu();
		   if(dispRho)
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
					(*currentState).delCurrent();
					delete currentState;
					(*mod).accept();
					currentState = mod;
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
				(*mod).reject();
				delete mod;
			}
    	}
    	else{
			delete mod;
    	}

    }

    /*cout << "\n";
   	(*currentState).displayMu();

    cout << "K " << (*currentState).valK() << "\n";*/
    //cout << "Ok " << bla << " cpt " << cptBla << " m " << bla / cptBla << "\n";
    //(*mod).reset();
    //delete mod;
    /*currentState.insertD(10010,3);
	currentState.insertD(10020,3);
	currentState.insertD(10030,3);
	currentState.insertD(10040,3);
	currentState.insertD(10050,3);
	currentState.insertD(10060,3);
	currentState.evalPriorMuDensity();

	PartitionAll<NucleoDirichletPA> *mod = currentState.clone();

	cout << "valKc " << currentState.valK() << "\n";
	cout << "valKm " << mod->valK() << "\n";
	cout << "sizem " << mod->size() << "\n";
	mod->insertD(10070,3);
	mod->insertD(10080,3);
	cout << "valKc " << currentState.valK() << "\n";
	cout << "valKm " << mod->valK() << "\n";
	cout << "sizem " << mod->size() << "\n";
	mod->evalW();
	mod->evalKdDim();
	mod->birth();*/
	//mod->w();
    //space_process::SpaceDirichlet<space_process::PartitionAll<space_process::NucleoDirichletPA> > currentState(fReads, rReads, 147);
    //cout << "Bla " << fReads[0] << "\n";
    //cout << "initMu " << currentState.initMu( 3) << "\n";
    //cout << "Aye " << startFReads[1] << "\n";
    //bla.initMu(currentState.newMu(), 3);
    //double mu = currentState.newMu();
    //cout << " Mu " << mu << "\n";
    //currentState.insert(nf);
    nf = startFReads.size();
    //nf = currentState.getP();
    nr = startRReads.size();
    tot = nbrIterations + kMax;

    List nbSeq = List::create(Rcpp::Named("nf") = nf, Rcpp::Named("nr") = nr, Rcpp::Named("tot") = tot, Rcpp::Named("mu") = (*currentState).mu());
    return nbSeq;
}
