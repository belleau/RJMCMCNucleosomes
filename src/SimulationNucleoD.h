/*
 * SimulationNucleoD.h
 *
 *  Created on: Sep 27, 2016
 *      Author: belleau
 */

#ifndef SIMULATIONNUCLEOD_H_
#define SIMULATIONNUCLEOD_H_


#include "SimulationNucleo.h"
#include <math.h>


namespace space_process {

    template< typename NucleoSpace>
    class SimulationNucleoD: public SimulationNucleo<NucleoSpace> {
        typedef std::vector<NucleoSpace *> NucleoSim;
        typedef typename NucleoSim::iterator itState;

        int d_kMax;
        Rcpp::List d_resStat;


    public:
        SimulationNucleoD(SegmentSeq const &segSeq,
                gsl_rng * rng, int kMax, long nbIteration=10000);
        virtual ~SimulationNucleoD();

        bool initMu(int lambda, int df = 3);
        int kMax();
        bool sampler();
        void simulate();
        double computeRho();
        void initResStat();

        void statSim();

        Rcpp::List simRapport();
    };

/*******************************************************************
 * Implementaton
 *******************************************************************/


    template< typename NucleoSpace>
    SimulationNucleoD<NucleoSpace>::SimulationNucleoD(SegmentSeq const &segSeq,
            gsl_rng * rng, int kMax, long nbIteration):
        SimulationNucleo<NucleoSpace>(segSeq, rng, nbIteration),
        d_kMax(kMax), d_resStat(R_NilValue) {
    }

    template< typename NucleoSpace>
    SimulationNucleoD<NucleoSpace>::~SimulationNucleoD() {

    }

    template< typename NucleoSpace>
    bool SimulationNucleoD<NucleoSpace>::initMu(int lambda, int df){
        bool flag= true;
        try{
            NucleoSpace *tmp = new NucleoSpace(this->segSeq(), this->rng(), kMax());
            (*tmp).setLambda(lambda);
            flag = (*tmp).initMu1(df);
            if(flag){
                (*tmp).prepSpace();
                (*tmp).addIteration();
                this->setCurrentState(tmp);
                this->pushState();
                this->setKMaxS(1);
            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
        return(flag);
    }

    template< typename NucleoSpace>
    int SimulationNucleoD<NucleoSpace>::kMax(){
        return(d_kMax);
    }

    template< typename NucleoSpace>
    bool SimulationNucleoD<NucleoSpace>::sampler(){
        bool flag = false;
        double u = gsl_ran_flat (this->rng(), 0, 1);
        double rhoP1 = 1.0;
        //int typeMv = 0;

        this->currentClone();   // init mod

        if((*(this->currentState())).valK() > 1){

            if(u > (*(this->currentState())).dK()){
                if(u <= ((*(this->currentState())).dK() + (*(this->currentState())).bK()) ){

                    flag = (*(this->mod())).birthR();

                    if(flag){
                        (*(this->mod())).prepSpace();
                        rhoP1 = (*(this->mod())).rhoP2() / (*(this->currentState())).bK();
                        rhoP1 *= (*(this->mod())).qalloc();
                        //typeMv = 4;
                    }
                }
                else{
                    flag = (*(this->mod())).mhR();
                    if(flag){
                        (*(this->mod())).prepSpace();
                        /* rhoP1 = 1.0; in the case mh */
                        //typeMv = 5;
                    }
                }
            }
            else{
                flag = (*(this->mod())).death();
                if(flag){
                    (*(this->mod())).prepSpace();
                    rhoP1 = (*(this->mod())).bK() / (*(this->currentState())).rhoP2();
                    rhoP1 /= (*(this->mod())).qalloc();
                    //typeMv = 3;
                }
            } /* End death */
        } /* End  case K > 1 */
        else{  /* case K == 1 */
            if(u <= 0.5){
                flag = (*(this->mod())).birthR();
                if(flag){
                    (*(this->mod())).prepSpace();
                    rhoP1 = (*(this->mod())).rhoP2() / (*(this->currentState())).bK(); // << (*mod).rhoP2()
                    rhoP1 *= (*(this->mod())).qalloc();
                    //typeMv = 1;
                }
            }
            else{
                flag = (*(this->mod())).mhR();
                if(flag){
                    (*(this->mod())).prepSpace();
                    /* rhoP1 = 1.0; in the case mh */
                    //typeMv = 2;
                }
            }
        }
        this->setRhoP1(rhoP1);
        //cout << " typeMv " << typeMv << "\n";
        return(flag);
    }

    template< typename NucleoSpace>
    double SimulationNucleoD<NucleoSpace>::computeRho(){

        double rho = this->rhoP1();
        rho *= exp(((*(this->mod())).kD() - (*(this->currentState())).kD()));
        rho *= ((*(this->mod())).priorMuDensity() / (*(this->currentState())).priorMuDensity());
        rho *= ((*(this->mod())).multinomial() / (*(this->currentState())).multinomial());
        rho = std::min(1.0, rho);
        return(rho);

    }


    template< typename NucleoSpace>
    void SimulationNucleoD<NucleoSpace>::simulate(){
        for(long i = 0; i< this->nbIterations();i++){

            bool flag = false; // Generate a new valide (with enought read in the partition) move
            double rho = 1;
            double u = 0;

            flag = sampler();
            if(flag){
                rho = computeRho();

                u = gsl_ran_flat (this->rng(), 0, 1);

                //cout << " u " << u << " Rho " << rho << "\n";
                if(rho > u){                   /* Accept mod */
                    this->acceptMod();
                    //(*(this->currentState())).addIteration();

                    this->pushState();

                    this->setKMaxS( fmax( this->kMaxS(), (*(this->currentState())).valK()));
                }
                else{                          /* Reject mod */
                    (*(this->currentState())).addIteration();
                    (*(this->mod())).reject(); // delete nucleosome add in mod
                    delete this->mod();
                }
            }
            else{               /* Enable to birth, mh or death in this iteration*/

                (*(this->currentState())).addIteration();
                delete this->mod();
            }
        }
    }


    template< typename NucleoSpace>
    Rcpp::List SimulationNucleoD<NucleoSpace>::simRapport(){
        return(d_resStat);

    }

    template< typename NucleoSpace>
    void SimulationNucleoD<NucleoSpace>::statSim(){

        int i = 0;

        //initResStat();

        Rcpp::NumericVector listK = Rcpp::NumericVector( Rcpp::Dimension(this->sizeState()));
        Rcpp::NumericMatrix mu = Rcpp::NumericMatrix( Rcpp::Dimension(this->sizeState(), this->kMaxS()));
        Rcpp::IntegerVector listIt = Rcpp::IntegerVector( Rcpp::Dimension(this->sizeState()));
        Rcpp::IntegerVector nbK = Rcpp::IntegerVector(this->kMaxS());
        Rcpp::NumericVector muHat = Rcpp::NumericVector(Rcpp::Dimension(this->kMaxS(), this->kMaxS()));

        for(itState it = this->beginState(); it != this->endState();it++){
            listK[i] = (*it)->valK();
            listIt[i] = (*it)->iteration();
            nbK[(*it)->valK()-1] += (*it)->iteration();

            std::vector<double> tmp = (*it)->mu();

            for(int j = 0; j < this->kMaxS(); j++){
                if(j < listK[i]){
                    mu[i  + j * this->sizeState()] = tmp[j];
                    muHat[((*it)->valK()-1) + j * this->kMaxS()] += (*it)->iteration() * tmp[j];
                }
                else{
                    mu[i + j * this->sizeState()] = 0;
                }
            }
            i++;
        }

        for(int j = 0; j < this->kMaxS(); j++){
            for(int l = 0; l < this->kMaxS(); l++){
                if(nbK[j] > 0)
                muHat[j + l * this->kMaxS()] /= nbK[j];
            }
        }

        d_resStat = Rcpp::List::create( Rcpp::Named("k") = listK
                , Rcpp::Named("k_max") = this->kMaxS(), Rcpp::Named("it") = listIt
                , Rcpp::Named("nbState") = this->sizeState(), Rcpp::Named("mu") = mu
                , Rcpp::Named("muHat") = muHat
                , Rcpp::Named("nbK") = nbK);
    }

} /* namespace space_process */

#endif /* SIMULATIONNUCLEOD_H_ */
