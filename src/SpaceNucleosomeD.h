/*
 * SpaceNucleosomeD.h
 *
 * Operation nucleosome from the space
 * Now only assign reads to nucleosomes
 *
 *  Created on: Jul 26, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOMED_H_
#define SPACENUCLEOSOMED_H_

#include <Rcpp.h>
#include <math.h>
//#include <gsl/gsl_blas.h>
#include "SpaceNucleosome.h"
#include "NucleoDirichlet.h"

namespace space_process{

    template<typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
    class SpaceNucleosomeD: public SpaceNucleosome<NucleoD>{
        typedef std::list<NucleoD*> containerNucleo;
        typedef typename containerNucleo::const_iterator itNucleo;
        double *d_w;
        double d_kD;
        double d_priorMuDensity;
        double d_multinomial;
        double d_qalloc;


        int d_lambda;
        unsigned int *d_dim;
        int d_c;
        int d_dfMax;

        double d_meanRead, d_r2, d_cMuDensity;
        double d_tB;

    public:
        SpaceNucleosomeD(SegmentSeq const &segSeq);

        SpaceNucleosomeD(SegmentSeq const &segSeq, int seed);

        SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng);

        SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng, int dfMax);

        virtual ~SpaceNucleosomeD(){};

        double tB();

        void setTB(double tB);

        inline int lambda(){
            return(d_lambda);
        };

        void setLambda(int l);

        int dfMax();

        double qalloc();

        void setQalloc(double qalloc);

        double meanRead();

        double r2();

        void evalW();

        double cMuDensity();

        void setPriorMuDensity(double priorMuDensity);

        double priorMuDensity();

        void insertD(double mu, int df);

        void evalPriorMuDensity();

        void evalDim();

        void evalKdDim();

        void evalMultinomial();

        double multinomial();

        double dK();

        double bK();

        double kD();

        double rhoP1();

        double rhoP2();

    protected:
        void setMeanRead(double meanRead);

        void setR2(double r2);

        void setCMuDensity(double cMuDensity);

        void delCurrentD();

    private:
        void setDefault();

    };

/*******************************************************************
 * Implementaton
 *******************************************************************/

    //template <typename NucleoD>
    //SpaceNucleosomeD<NucleoD>::

    template <typename NucleoD>
    SpaceNucleosomeD<NucleoD>::SpaceNucleosomeD(SegmentSeq const &segSeq)
        :SpaceNucleosome<NucleoD>(segSeq), d_w(NULL), d_dim(NULL),
            d_dfMax(30){
        setDefault();
    }

    template <typename NucleoD>
    SpaceNucleosomeD<NucleoD>::SpaceNucleosomeD(SegmentSeq const &segSeq, int seed)
        :SpaceNucleosome<NucleoD>(segSeq, seed), d_w(NULL),
            d_dim(NULL), d_dfMax(30){
        setDefault();
    }

    template <typename NucleoD>
    SpaceNucleosomeD<NucleoD>::SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng)
        :SpaceNucleosome<NucleoD>(segSeq, rng), d_w(NULL), d_dim(NULL),
            d_dfMax(30){
        setDefault();
    }

    template <typename NucleoD>
    SpaceNucleosomeD<NucleoD>::SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng, int dfMax)
        :SpaceNucleosome<NucleoD>(segSeq, rng), d_w(NULL), d_dim(NULL),
            d_dfMax(dfMax){
        setDefault();
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::tB(){
        return(d_tB);
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setTB(double tB){
        d_tB = tB;
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setLambda(int l){
        d_lambda = l;
    }

    template <typename NucleoD>
    int SpaceNucleosomeD<NucleoD>::dfMax(){
        return(d_dfMax);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::qalloc(){
        return(d_qalloc);
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setQalloc(double qalloc){
        d_qalloc = qalloc;
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::meanRead(){
        return(d_meanRead);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::r2(){
        return(d_r2);
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::evalW(){
        /* gsl_ran_dirichlet (const gsl_rng * r, size_t K, const double alpha[], double theta[]) */
        int k = this->valK();

        try{
            double *alpha = new double[k];

            std::fill_n(alpha, k, 1.0);

            d_w = new double[k];
            gsl_ran_dirichlet (this->rng(), k, alpha, d_w);

            delete[] alpha;
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }

    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::cMuDensity(){
        return(d_cMuDensity);
    };

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setPriorMuDensity(double priorMuDensity){
        d_priorMuDensity = priorMuDensity;
    };

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::priorMuDensity(){
        return(d_priorMuDensity);
    };

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::insertD(double mu, int df){
        try{
            NucleoD *u = new NucleoD(mu, df, this->segSeq(), this->rng());
            this->insert(u);
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }

    };

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::evalPriorMuDensity(){

        double m = meanRead(); /* Mean of the read*/
        double result = 0;


        itNucleo nucleoIt = this->nucleoBegin();
        result = 0;
        if(this->valK() == 1){
            result = 2 * pow((**nucleoIt).mu()  - m, 2);
        }
        else{
            if(this->valK() == 2){
                double v [2]= {(**nucleoIt++).mu()  - m, (**nucleoIt).mu() - m};
                result = (2 * v[0] -  v[1]) * v[0] + (v[1] - v[0]) *v[1];
            }
            if(this->valK() > 2)
            {
                /* matrix multiplication t(mu) Omega mu */
                double v [3] = {(**nucleoIt++).mu() - m, 0, 0};
                v[1] = (**nucleoIt++).mu()  - m;
                result = (2 * v[0] -  v[1]) * v[0];

                int i = 2;
                do{
                    v[i%3] =  (**nucleoIt++).mu() - m;

                    result += (2 * v[(i-1)%3] - v[(i)%3] - v[(i-2)%3]) * v[(i-1)%3];
                    i++;
                }while(nucleoIt != this->nucleoEnd());
                // result du dernier


                result +=  (v[(i-1)%3] - v[(i-2)%3]) * v[(i-1)%3];

            }
        }

        double tmpC = pow(cMuDensity(), -1 * this->valK() / 2.0);


        setPriorMuDensity(tmpC * exp(- result/ (2.0 * r2()) ));
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::evalDim(){

        try{
            d_dim = new unsigned int[this->valK()];

            int n = 0;
            for(itNucleo it = this->nucleoBegin(); it != this->nucleoEnd(); it++)
            {

                d_dim[n++] = (*it)->dimN();
            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::evalKdDim(){

        try{
            d_dim = new unsigned int[this->valK()];
            int i;
            int s = this->sizeFReads() + this->sizeRReads();

            double *yRead = new double[s];

            std::fill_n(yRead, s, 0.0);

            int n = 0;
            for(itNucleo it = this->nucleoBegin(); it != this->nucleoEnd(); it++)
            {
                i = 0;
                d_dim[n] = (*it)->dimN();//(*it)->sizeF() + (*it)->sizeR();
                for(std::vector<double>::const_iterator  itF = (*it)->bFBegin(); itF != (*it)->bFEnd(); itF++){
                    yRead[i++] += d_w[n] * (*itF);
                }
                for(std::vector<double>::const_iterator  itR = (*it)->bRBegin(); itR != (*it)->bREnd(); itR++){
                    yRead[i++] += d_w[n] * (*itR);
                }

                n++;

            }

            d_kD = 0;

            for(int j = 0; j < s; j++){
                d_kD += log(yRead[j]);
            }


            delete[] yRead;
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::evalMultinomial(){
        d_multinomial = gsl_ran_multinomial_pdf (this->valK(), d_w, d_dim);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::multinomial(){
        return(d_multinomial);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::dK(){
        double d = 0;
        if(this->valK() > 1){
            d = 0.5 * std::min(1.0, this->valK() / ((double) lambda()));
        }
        return(d);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::bK(){
        double d = 0;
        d = 0.5 * std::min(1.0,  lambda() / (double)(this->valK() + 1.0) );
        return(d);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::kD(){
        return(d_kD);
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::rhoP1(){
        return(priorMuDensity() * multinomial());
    }

    template <typename NucleoD>
    double SpaceNucleosomeD<NucleoD>::rhoP2(){
        return(dK() * (lambda() / ((double) this->valK())) );
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setMeanRead(double meanRead){
        d_meanRead = meanRead;
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setR2(double r2){
        d_r2 = r2;
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setCMuDensity(double cMuDensity){
        d_cMuDensity = cMuDensity;
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::delCurrentD(){
        delete[] d_dim;
        d_dim = NULL;
        delete[] d_w;
        d_w = NULL;
    }

    template <typename NucleoD>
    void SpaceNucleosomeD<NucleoD>::setDefault(){
        d_r2 = pow((this->maxPos() - this->minPos()),2);
        d_meanRead = (this->maxPos() + this->minPos()) / 2.0;
        d_cMuDensity = M_PI * d_r2 / 2.0;
        d_lambda = 3;
    }


}

#endif /* SPACENUCLEOSOMED_H_ */
