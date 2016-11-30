/*
 * SpaceNucleosome.h
 *
 *  Created on: Aug 5, 2016
 *      Author: belleau
 */

#ifndef SPACENUCLEOSOME_H_
#define SPACENUCLEOSOME_H_

#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <list>
#include <math.h>
#include <time.h>

#include "Nucleosome.h"
#include "SegmentSeq.h"

namespace space_process {

template<typename NucleoClass>    /***** BEWARE NucleoClass Must inherit from Nucleosome *****/
class SpaceNucleosome {
    typedef std::list<NucleoClass*> containerNucleo;
    typedef typename containerNucleo::iterator itNucleo;
    typedef std::vector<itNucleo> vecItNucleo;
    typedef typename vecItNucleo::iterator itItVNucleo;

    typedef std::vector<NucleoClass*> vecNucleo;
    typedef typename vecNucleo::iterator itVNucleo;

    vecNucleo d_modNucleo;
    vecItNucleo d_addNucleo;

    SegmentSeq const &d_segSeq;
    containerNucleo d_nucleosomes; /* List of nucleosomes */



    int d_valK;            // Number of nucleosome
    gsl_rng *d_rng;       // random number generator
    long d_iteration;     // Number of iteration not accept

public:
    SpaceNucleosome(SegmentSeq const &segSeq);

    SpaceNucleosome(SegmentSeq const &segSeq, int seed);

    SpaceNucleosome(SegmentSeq const &segSeq, gsl_rng * rng);


    virtual ~SpaceNucleosome(){};

    int size();

    bool empty();

    void pushNucleo(NucleoClass *u);

    void insertNucleo(itNucleo it, NucleoClass *u);

    void setRng(gsl_rng *rng);

    int valK();

    double minPos();

    double maxPos();

    long sizeFReads();

    long sizeRReads();

    void displayMu();

    std::vector<double> mu();

    void eraseNucleo(itNucleo it);

    void addIteration();

    long iteration();

protected:
    gsl_rng * rng();

    SegmentSeq const &segSeq();

    void pushModNucleo(NucleoClass *u);

    void pushAddNucleo(itNucleo &u);

    void resetNucleo();

    void resetMod();

    void resetAdd();

    void clearAdd();

    containerNucleo &nucleosomes();

    void setNucleosomes(containerNucleo &nucleosomes);

    void setValK(int k);
private:

    void setRNG();

    void setRNG(int seed);


protected:
    itNucleo nucleoBegin(){
        return(d_nucleosomes.begin());
    };

    itNucleo nucleoEnd(){
        return(d_nucleosomes.end());
    };
    itNucleo nucleosomes(itNucleo itPos, int start, int pos){
        int i = start;
        do{
            itPos++;
        }while(pos > i++ && itPos != d_nucleosomes.end());
        if(pos < i)
        {
            itPos--;
        }
        return(itPos);
    };

}; /* Class SpaceNucleosome */

/*******************************************************************
 * Implementaton
 *******************************************************************/

    //template<typename NucleoClass>
    //SpaceNucleosome<NucleoClass>::

    template<typename NucleoClass>
    SpaceNucleosome<NucleoClass>::SpaceNucleosome(SegmentSeq const &segSeq):
        d_segSeq(segSeq), d_valK(0), d_iteration(0){
        setRNG();
    }

    template<typename NucleoClass>
    SpaceNucleosome<NucleoClass>::SpaceNucleosome(SegmentSeq const &segSeq, int seed):
        d_segSeq(segSeq), d_valK(0), d_iteration(0){
        setRNG(seed);
    }

    template<typename NucleoClass>
    SpaceNucleosome<NucleoClass>::SpaceNucleosome(SegmentSeq const &segSeq, gsl_rng * rng):
        d_segSeq(segSeq), d_valK(0), d_rng(rng), d_iteration(0){
    }

    template<typename NucleoClass>
    int SpaceNucleosome<NucleoClass>::size(){
        return(d_nucleosomes.size());
    }

    template<typename NucleoClass>
    bool SpaceNucleosome<NucleoClass>::empty(){
        return(d_nucleosomes.empty());
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::pushNucleo(NucleoClass *u){
        d_nucleosomes.push_back(u);
        d_valK++;
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::insertNucleo(itNucleo it, NucleoClass *u){
        d_nucleosomes.insert(it, u);
        d_valK++;
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::setRng(gsl_rng *rng){
        d_rng = rng;
    }

    template<typename NucleoClass>
    int SpaceNucleosome<NucleoClass>::valK(){
        return(d_valK);
    }

    template<typename NucleoClass>
    double SpaceNucleosome<NucleoClass>::minPos(){
        return(d_segSeq.minPos());
    }

    template<typename NucleoClass>
    double SpaceNucleosome<NucleoClass>::maxPos(){
        return(d_segSeq.maxPos());
    }

    template<typename NucleoClass>
    long SpaceNucleosome<NucleoClass>::sizeFReads(){
        return(d_segSeq.sizeFReads());
    }

    template<typename NucleoClass>
    long SpaceNucleosome<NucleoClass>::sizeRReads(){
        return(d_segSeq.sizeRReads());
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::displayMu(){
        Rcpp::Rcout << "Mu";
        for(itNucleo it = d_nucleosomes.begin() ; it != d_nucleosomes.end(); it++){
            Rcpp::Rcout << " " << (*it)->mu();
            Rcpp::Rcout << " : " << (*it)->avg();
        }
        Rcpp::Rcout << "\n";
    }

    template<typename NucleoClass>
    std::vector<double> SpaceNucleosome<NucleoClass>::mu(){
        //Rcpp::NumericVector mu = Rcpp::NumericVector(valK());
        std::vector<double> mu(valK(), 0.0);
        int i = 0;
        for(itNucleo it = d_nucleosomes.begin() ; it != d_nucleosomes.end(); it++){
            mu[i++] = (*it)->mu();
        }
        return(mu);
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::eraseNucleo(itNucleo it){ //itNucleo it
        d_nucleosomes.erase(it);
        d_valK--;
        //d_nucleosomes.erase(d_nucleosomes.begin());
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::addIteration(){
        d_iteration++;
    }

    template<typename NucleoClass>
    long SpaceNucleosome<NucleoClass>::iteration(){
        return(d_iteration);
    }

    template<typename NucleoClass>
    gsl_rng * SpaceNucleosome<NucleoClass>::rng(){
        return(d_rng);
    }

    template<typename NucleoClass>
    SegmentSeq const &SpaceNucleosome<NucleoClass>::segSeq(){
        return(d_segSeq);
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::pushModNucleo(NucleoClass *u){
        d_modNucleo.push_back(u);
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::pushAddNucleo(itNucleo &u){
        d_addNucleo.push_back(u);
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::resetNucleo(){
        for(itNucleo it = d_nucleosomes.begin(); it != d_nucleosomes.end();it++){
            if(*it != NULL){
                delete *it;
                *it = NULL;
            }
        }
        d_nucleosomes.clear();
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::resetMod(){
        d_modNucleo.clear();
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::resetAdd(){

        for(itItVNucleo it = d_addNucleo.begin(); it != d_addNucleo.end();it++){
            if(**it != NULL){
                delete **it;
                **it = NULL;
            }
        }
        d_addNucleo.clear();
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::clearAdd(){
        d_addNucleo.clear();
    }

    template<typename NucleoClass>
    std::list<NucleoClass*> &SpaceNucleosome<NucleoClass>::nucleosomes(){
        return(d_nucleosomes);
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::setNucleosomes(containerNucleo &nucleosomes){

        d_nucleosomes = nucleosomes;
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::setValK(int k){
        d_valK= k;
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::setRNG(){
        const gsl_rng_type * T;
        long seed;

        T = gsl_rng_default;

        d_rng = gsl_rng_alloc (T);     // pick random number generator
        seed = time (NULL) * getpid();
        gsl_rng_set (d_rng, seed);                  // set seed
    }

    template<typename NucleoClass>
    void SpaceNucleosome<NucleoClass>::setRNG(int seed){
        const gsl_rng_type * T;

        T = gsl_rng_default;

        d_rng = gsl_rng_alloc (T);     // pick random number generator

        gsl_rng_set (d_rng, seed);                  // set seed
    }


} /* namespace space_process */

#endif /* SPACENUCLEOSOME_H_ */
