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

template <typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
class PartitionAll: public SpaceNucleosomeD<NucleoD> {
    typedef PartitionAll<NucleoD> NucleoSpace;
    typedef std::list<NucleoD*> containerNucleo;
    typedef typename containerNucleo::iterator iteratorNucleo;

    std::vector<double> *d_y;
    long d_ySize;

public:

    PartitionAll(SegmentSeq const &segSeq);

    PartitionAll(SegmentSeq const &segSeq, int seed);

    PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng);

    PartitionAll(SegmentSeq const &segSeq, int seed
            , std::vector<double> *y, long ySize);

    PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng
            , std::vector<double> *y, long ySize);

    PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, int dfMax);

    PartitionAll(SegmentSeq const &segSeq, int seed
            , std::vector<double> *y, long ySize, int dfMax);

    PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng
            , std::vector<double> *y, long ySize, int dfMax);

    virtual ~PartitionAll(){};

    bool initMu(int df);

    bool initMu1(int df);

    long setFoward(std::vector<double>::const_iterator fStart
            , std::vector<double>::const_iterator fEnd
            , double start, double end, NucleoD &u);


    long setReverse(std::vector<double>::const_iterator rStart
            , std::vector<double>::const_iterator rEnd
            , double start, double end, NucleoD &u);

    NucleoSpace * clone();

    bool setNucleoD(NucleoD *u, double aF, double aR);

    bool setNucleoDR(NucleoD *u, double aF, double aR, NucleoD *old);

    bool birth();

    bool birthR();

    bool death();

    bool deathR();

    bool mh();

    bool mhR();

    void prepSpace();

    void delCurrent();

    void delMod();

    void reset();

    void reject();

private:
    int getLimit(double start, double end
            , std::vector<double>::const_iterator &startIt
            , std::vector<double>::const_iterator &endIt
            , long &l, bool excEnd = false);

    inline bool yEmpty(){
            return((*d_y).empty());
    };

    inline bool ySize(){
            return(d_ySize);
    };

    inline std::vector<double>::iterator yBegin(){
        return((*d_y).begin());
    };

    inline std::vector<double>::iterator yEnd(){
        return((*d_y).end());
    };

    inline double y(int i){
        return((*d_y)[i]);
    };


};

/*******************************************************************
 * Implementaton
 *******************************************************************/

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq)
        :SpaceNucleosomeD<NucleoD>(segSeq), d_y(new std::vector<double>){

        d_y->insert(yEnd(), segSeq.beginFR(), segSeq.endFR());
        d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

        std::sort(d_y->begin(),d_y->end());
        d_ySize = (*d_y).size();
    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, int seed)
        :SpaceNucleosomeD<NucleoD>(segSeq,seed), d_y(new std::vector<double>){
        d_y->insert(yEnd(), segSeq.beginFR(), segSeq.endFR());
        d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

        std::sort(d_y->begin(),d_y->end());
        d_ySize = (*d_y).size();
    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng)
        :SpaceNucleosomeD<NucleoD>(segSeq,rng), d_y(new std::vector<double>){
        d_y->insert(yEnd(), segSeq.beginFR(), segSeq.endFR());
        d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

        std::sort(d_y->begin(),d_y->end());
        d_ySize = (*d_y).size();
    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, int seed, std::vector<double> *y, long ySize)
        :SpaceNucleosomeD<NucleoD>(segSeq, seed), d_y(y){ //

    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, std::vector<double> *y, long ySize)
        :SpaceNucleosomeD<NucleoD>(segSeq, rng), d_y(y), d_ySize(ySize){

    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, int dfMax)
        :SpaceNucleosomeD<NucleoD>(segSeq,rng, dfMax), d_y(new std::vector<double>){
        d_y->insert(yEnd(), segSeq.beginFR(), segSeq.endFR());
        d_y->insert(yEnd(), segSeq.beginRR(), segSeq.endRR());

        std::sort(d_y->begin(),d_y->end());
        d_ySize = (*d_y).size();
    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, int seed, std::vector<double> *y, long ySize, int dfMax)
        :SpaceNucleosomeD<NucleoD>(segSeq, seed, dfMax), d_y(y){
    }

    template <typename NucleoD>
    PartitionAll<NucleoD>::PartitionAll(SegmentSeq const &segSeq, gsl_rng * rng, std::vector<double> *y, long ySize, int dfMax)
        :SpaceNucleosomeD<NucleoD>(segSeq, rng, dfMax), d_y(y), d_ySize(ySize){
    }


    template <typename NucleoD>
    bool PartitionAll<NucleoD>::initMu(int df){

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
                    }
                    else{
                        delete u;
                    }

                }while(!(flag) && cpt == 1000);
                if(!(flag)){
                    Rcpp::stop("Problem with the number of reads to initialise mu\n");
                }
            }
            else{
                Rcpp::stop("No reads\n");
            }
        }

        return(flag);
    }
    template <typename NucleoD>
    bool PartitionAll<NucleoD>::initMu1(int df){

        bool flag = true;
        NucleoD *u;
        if(this->empty()){
            if(!(yEmpty()))
            {
                double mu= gsl_ran_flat(this->rng(), this->minPos(), this->maxPos());

                u = new NucleoD(mu, df, this->segSeq(), this->rng());

                double start = this->minPos();
                double end = this->maxPos();
                u->setSigmaF(1);
                u->setSigmaR(1);
                u->setDf(df);
                u->setDelta(gsl_ran_flat(this->rng(), 0, 2 * (mu - start)));
                setFoward(yBegin(), yEnd(), start, mu, *u);
                setReverse(yBegin(), yEnd(), mu, end, *u);
                u->evalBF();
                u->evalBR();
                u->setAF(start);
                u->setAR(end);
                this->pushNucleo(u);

            }
            else{
                Rcpp::stop("No reads\n");
            }
        }

        return(flag);
    }

    template <typename NucleoD>
    long PartitionAll<NucleoD>::setFoward(std::vector<double>::const_iterator fStart
                        , std::vector<double>::const_iterator fEnd
                        , double start, double end, NucleoD &u){
        long l = 0;
        int cpt = getLimit(start,end, fStart, fEnd, l);

        u.setFStartPos(fStart, fEnd, cpt);
        return(l);
    }

    template <typename NucleoD>
    long PartitionAll<NucleoD>::setReverse(std::vector<double>::const_iterator rStart
            , std::vector<double>::const_iterator rEnd
            , double start, double end, NucleoD &u){

        long l = 0;
        int cpt = getLimit(start,end, rStart, rEnd, l);

        u.setRStartPos(rStart, rEnd, cpt);
        return(l);
    }

    template <typename NucleoD>
    PartitionAll<NucleoD> * PartitionAll<NucleoD>::clone(){

        NucleoSpace *a = new NucleoSpace(this->segSeq()
                                , this->rng(), d_y, ySize()
                                , this->dfMax());

        a->setValK(this->valK());
        a->setNucleosomes(this->nucleosomes());
        a->setLambda(this->lambda());
        a->setMeanRead(this->meanRead());
        a->setR2(this->r2());
        a->setCMuDensity(this->cMuDensity());

        return(a);
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::setNucleoD(NucleoD *u, double aF, double aR){

        bool flag = false;
        long l;
        std::vector<double>::const_iterator startIt, endIt;
        int dimNucleo = getLimit( aF, aR, startIt, endIt, l, true);
        (*u).setDimN(dimNucleo);
        if(l > 1){

            (*u).setAvg(accumulate( startIt, endIt, 0.0)/ ((double) dimNucleo));

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

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::setNucleoDR(NucleoD *u, double aF
                                        , double aR, NucleoD *old){

        bool flag = true;
        long l;
        std::vector<double>::const_iterator startIt, endIt;
        int dimNucleo = getLimit( aF, aR, startIt, endIt, l, true);
        (*u).setDimN(dimNucleo);
        if(l > 1){

            (*u).setAvg(accumulate( startIt, endIt, 0.0)/((double) dimNucleo));

            setFoward(startIt, endIt, aF, (*u).avg(), *u);

            setReverse(startIt, endIt, (*u).avg(), aR, *u);

            (*u).setSigmaF((*old).sigmaF());
            (*u).setSigmaR((*old).sigmaR());
            (*u).setDelta((*old).delta());
            (*u).setBF( (*old).bF() );
            (*u).setBR( (*old).bR() );
            (*u).setAF(aF);
            (*u).setAR(aR);


        }
        return(flag);
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::birth(){
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

                i = (int) gsl_ran_flat (this->rng(), 0, k+1); // select nucleo between k+1 for (maxPos())

                //double startBirth;
                double aFBirth, aRBirth;

                if(i > 0){
                    it1 = this->nucleosomes(this->nucleoBegin(), 0, i-1); // go to the position i-1

                    muBef = (*it1)->mu();

                    if(i < k){
                        //it2 = this->nucleosomes(it1, i-1, i);
                        it2 = it1;
                        it2++;
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

                muBirth = gsl_ran_flat(this->rng(), muBef, muNext); // New mu



                if(i > 0) // Modify nucleosome i-1
                {
                    aFBirth = gsl_ran_flat(this->rng(), muBef, muBirth);

                    uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1

                    // Eval the var interne of uBef
                    // and check if enough reads between (*it1)->aF(), aFBirth
                    flag = !(setNucleoD(uBef, (*it1)->aF(), aFBirth));


                    if(!(flag) && i < k){ // Modify nucleosome i+1

                        aRBirth = gsl_ran_flat(this->rng(), muBirth, muNext);
                        uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New i+1
                        flag = !(setNucleoD(uNext, aRBirth, (*it2)->aR()));

                    }
                    else{
                        aRBirth = this->maxPos() + 1;
                    }
                    if(!(flag)){

                        int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);

                        uBirth = new NucleoD(muBirth, df, this->segSeq(), this->rng());
                        flag = !(setNucleoD(uBirth, aFBirth, aRBirth));

                    }
                }
                else{ // i == 0

                    aFBirth = muBef;
                    aRBirth = gsl_ran_flat(this->rng(), muBirth, muNext);
                    uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New mu 0
                    flag = !(setNucleoD(uNext, aRBirth, (*it2)->aR()));

                    if(!(flag)){

                        int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);
                        uBirth = new NucleoD(muBirth, df, this->segSeq(), this->rng());
                        flag = !(setNucleoD(uBirth, aFBirth, aRBirth));

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
                this->setQalloc(muNext - muBef);
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
            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
        return(!(flag));
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::birthR(){
        int cpt = 0;
        bool flag = false;
        iteratorNucleo it1, it2;
        NucleoD *uBef, *uBirth, *uNext;
        double muBef, muBirth, muNext;
        int k = this->valK();
        int i = 0;
        try{
            do{
                double aFBirth, aRBirth;
                uBef = NULL;
                uBirth = NULL;
                uNext = NULL;

                flag = false;
                cpt++;

                i = (int) gsl_ran_flat (this->rng(), 0, k); // select nucleo
                                    //no birth between mu(k) and (maxPos())

                //double startBirth;



                this->setTB(i);
                if(i > 0){
                    it1 = this->nucleosomes(this->nucleoBegin(), 0, i-1); // go to the position i-1
                    muBef = (*it1)->mu();

                    it2 = it1;
                    it2++;   // go to i

                    muNext = (*it2)->mu();

                }
                else{
                    muBef = this->minPos();
                    it2 = this->nucleoBegin();
                    muNext = (*it2)->mu();

                }

                muBirth = gsl_ran_flat(this->rng(), muBef, muNext); // New mu


                aRBirth = gsl_ran_flat(this->rng(), muBirth, muNext);

                if(i == 0){
                    aFBirth = this->minPos();
                }
                else{
                    aFBirth = gsl_ran_flat(this->rng(), muBef, muBirth);
                }

                if(i > 0) // Modify nucleosome i-1
                {


                    uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1

                    // Eval the var interne of uBef
                    // and check if enough reads between (*it1)->aF(), aFBirth

                    setNucleoDR(uBef, (*it1)->aF(), aFBirth, *it1);

                    //if(!(flag)){ // Modify nucleosome i

                    uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New i+1
                    //flag = !(setNucleoDR(uNext, aRBirth, (*it2)->aR(), *it2));
                    setNucleoDR(uNext, aRBirth, (*it2)->aR(), *it2);
                    //}

                    //if(!(flag)){

                    int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);

                    uBirth = new NucleoD(muBirth, df, this->segSeq(), this->rng());

                    flag = !(setNucleoD(uBirth, aFBirth, aRBirth));

                }
                else{ // i == 0

                    uNext = new NucleoD(muNext, (*it2)->df(), this->segSeq(), this->rng()); // New mu 0

                    setNucleoDR(uNext, aRBirth, (*it2)->aR(), *it2);



                    int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);
                    uBirth = new NucleoD(muBirth, df, this->segSeq(), this->rng());

                    flag = !(setNucleoD(uBirth, aFBirth, aRBirth));

                } // end i == 0
                if(flag){
                    delete uBef;
                    uBef = NULL;
                    delete uNext;
                    uNext = NULL;
                    delete uBirth;
                    uBirth = NULL;
                }

            }while(flag && cpt < 1000);

            if(!(flag))
            {
                this->setQalloc(muNext - muBef);
                if(i > 0){

                    this->pushModNucleo(*it1);
                    *it1 = uBef;
                    this->pushAddNucleo(it1);

                }



                this->pushModNucleo(*it2);
                *it2 = uNext;
                this->insertNucleo(it2, uBirth);


                this->pushAddNucleo(it2);
                this->pushAddNucleo(--it2);


            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
        return(!(flag));
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::death(){
        int cpt = 0;
        bool flag = false;
        iteratorNucleo it1, it2, it3;
        NucleoD *uBef, *uNext;
        double muBef, muNext, vNext;
        int k = this->valK();
        int i = 0;

        if(k > 1){ // Don't remove the last one
            try{
                do{
                    uBef = NULL;
                    uNext = NULL;
                    double a = this->maxPos() + 1;
                    flag = false;
                    cpt++;

                    /* select the nucleosome to remove */
                    i = (int) gsl_ran_flat (this->rng(), 0, k);

                    /* get the nucleosome */
                    it2 = this->nucleosomes(this->nucleoBegin(), 0, i);

                    if( i > 0){ // get nucleo before
                        it1 = it2; // go to the position i-1
                        it1--;
                        muBef = (*it1)->mu();
                        //a = (*it1)->aR();
                    }
                    else{ /* case i == 0 */
                        muBef = this->minPos();
                    }

                    if( i < (k-1) ){ // get nucleo after
                        it3 = it2;
                        it3++; // go to the position i+1
                        muNext = (*it3)->mu();
                        vNext = muNext;
                    }
                    else{ /* case i is the last nucleo */

                        vNext = this->maxPos();
                        muNext = vNext + 1; /* to include the last read*/
                    }

                    if(i > 0){

                        uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1

                        if(i < (k-1)){
                            a = gsl_ran_flat(this->rng(), muBef, muNext);
                        } // by default a = maxPos() + 1

                        /* Init uBef and validate the number reads between
                         * a and (*it1)->aR() */
                        flag = !(setNucleoD(uBef, (*it1)->aF(), a));
                    }
                    else{
                        a = this->minPos();
                    }

                    /* here a is the aR of uBef and aF of uNext */
                    if(!(flag) && i < (k-1)){ // Modify nucleosome i+1

                        uNext = new  NucleoD(muNext, (*it3)->df()
                                        , this->segSeq(), this->rng()); // New i+1
                        /* Init uBef and validate the number reads between
                         * a and (*it3)->aR() */
                        flag = !(setNucleoD(uNext, a, (*it3)->aR()));

                    }

                    if(flag){
                        delete uBef;
                        uBef = NULL;
                        delete uNext;
                        uNext = NULL;
                    }

                }while(flag && cpt == 1000);

                if(!(flag))
                {
                    this->setQalloc(vNext - muBef);

                    if(i > 0){ // Update i-1

                        this->pushModNucleo(*it1);
                        *it1 = uBef;
                        this->pushAddNucleo(it1);

                    }

                    if(i < (k-1)){ // Update i+1

                        this->pushModNucleo(*it3);
                        *it3 = uNext;
                        this->pushAddNucleo(it3);
                    }

                    /* Remove i */
                    this->pushModNucleo(*it2);
                    this->eraseNucleo(it2);

                }
            }
            catch(std::bad_alloc&) {
                Rcpp::stop("Memory problem\n");
            }
        }
        else{
            flag = true; // k <= 1
        }
        return(!(flag));
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::deathR(){
        int cpt = 0;
        bool flag = false;
        iteratorNucleo it1, it2, it3;
        NucleoD *uBef, *uNext;
        double muBef, muNext;
        int k = this->valK();
        int i = 0;

        if(k > 1){ // Don't remove the last one
            try{
                do{
                    uBef = NULL;
                    uNext = NULL;
                    double a = this->maxPos() + 1;
                    flag = false;
                    cpt++;
                    i = (int) gsl_ran_flat (this->rng(), 0, k);
                    this->setTB(i);
                    it2 = this->nucleosomes(this->nucleoBegin(), 0, i);

                    if( i > 0){ // Update nucleo before
                        it1 = it2; // go to the position i-1
                        it1--;
                        muBef = (*it1)->mu();


                    }
                    else{
                        muBef = this->minPos();
                    }

                    if( i < (k-1) ){ // Update nucleo after
                        it3 = it2;
                        it3++; // go to the position i+1
                        muNext = (*it3)->mu();

                    }
                    else{
                        muNext = this->maxPos();
                    }


                    if(i > 0){

                        a = (*it1)->aR();

                    }
                    else{
                        a = this->minPos();
                    }
                    if(i < (k-1)){ // Modify nucleosome i+1

                        uNext = new  NucleoD(muNext, (*it3)->df(), this->segSeq(), this->rng()); // New i+1
                        // Enough reads between a and (*it3)->aR()
                        flag = !(setNucleoDR(uNext, a, (*it3)->aR(), *it3));


                    }
                    else{ // i == k-1

                        uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng());

                        // Enough reads between a and (*it3)->aR()
                        //flag = !(setNucleoDR(uBef, (*it1)->aF(), this->maxPos(), *it1));
                        setNucleoDR(uBef, (*it1)->aF(), this->maxPos(), *it1);
                    }
                    if(flag){
                        delete uBef;
                        uBef = NULL;
                        delete uNext;
                        uNext = NULL;
                    }

                }while(flag && cpt < 1000);

                if(!(flag))
                {
                    this->setQalloc(muNext - muBef);


                    if(i < (k-1)){

                        this->pushModNucleo(*it3);
                        *it3 = uNext;
                        this->pushAddNucleo(it3);
                    }
                    else{
                        this->pushModNucleo(*it1);
                        *it1 = uBef;
                        this->pushAddNucleo(it1);
                    }

                    this->pushModNucleo(*it2);

                    this->eraseNucleo(it2);

                }
            }
            catch(std::bad_alloc&) {
                Rcpp::stop("Memory problem\n");
            }
        }
        else{
            flag = true; // k <= 1
        }
        return(!(flag));
    }


    template <typename NucleoD>
    bool PartitionAll<NucleoD>::mh(){
        int cpt = 0;
        bool flag = false;
        iteratorNucleo it1, it2, it3;
        NucleoD *uBef, *uMH, *uNext;
        double muBef, muMH, muNext, vNext;
        double aF, aR;
        int k = this->valK();
        int i = 0;


        try{
            do{
                uBef = NULL;
                uMH = NULL;
                uNext = NULL;

                flag = false;
                cpt++;
                i = (int) gsl_ran_flat (this->rng(), 0, k);
                it2 = this->nucleosomes(this->nucleoBegin(), 0, i);

                if( i > 0){ // Update nucleo before
                    it1 = it2; // go to the position i-1
                    it1--;
                    muBef = (*it1)->mu();

                }
                else{
                    muBef = this->minPos();
                }

                if( i < (k-1) ){ // Update nucleo after
                    it3 = it2;
                    it3++; // go to the position i+1
                    muNext = (*it3)->mu();
                    vNext = muNext;
                }
                else{

                    vNext = this->maxPos();
                    muNext = vNext + 1;
                }

                muMH = gsl_ran_flat(this->rng(), muBef, vNext); // New mu
                int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);

                uMH = new  NucleoD(muMH, df, this->segSeq(), this->rng());

                aF = gsl_ran_flat(this->rng(), muBef, muMH);
                aR = gsl_ran_flat(this->rng(), muMH, muNext);

                flag = !(setNucleoD(uMH, aF, aR));

                if(!(flag)){
                    if(i > 0){

                        uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1

                        // Enough reads between a and (*it1)->aR()
                        flag = !(setNucleoD(uBef, (*it1)->aF(), aF));
                    }
                    if(!(flag) && i < (k-1)){ // Modify nucleosome i+1

                        uNext = new  NucleoD(muNext, (*it3)->df(), this->segSeq(), this->rng()); // New i+1
                        // Enough reads between a and (*it3)->aR()
                        flag = !(setNucleoD(uNext, aR, (*it3)->aR()));
                    }
                }
                if(flag){
                    delete uBef;
                    uBef = NULL;
                    delete uNext;
                    uNext = NULL;
                }

            }while(flag && cpt == 1000);

            if(!(flag))
            {
                this->setQalloc(vNext - muBef);
                this->pushModNucleo(*it2);
                *it2 = uMH;
                this->pushAddNucleo(it2);

                if(i > 0){

                    this->pushModNucleo(*it1);
                    *it1 = uBef;
                    this->pushAddNucleo(it1);

                }

                if(i < (k-1)){

                    this->pushModNucleo(*it3);
                    *it3 = uNext;
                    this->pushAddNucleo(it3);
                }

            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
        return(!(flag));
    }

    template <typename NucleoD>
    bool PartitionAll<NucleoD>::mhR(){
        int cpt = 0;
        bool flag = false;
        iteratorNucleo it1, it2, it3;
        NucleoD *uBef, *uMH, *uNext;
        double muBef, muMH, muNext, vNext;
        double aF, aR;
        int k = this->valK();
        int i = 0;

        try{
            do{
                uBef = NULL;
                uMH = NULL;
                uNext = NULL;

                flag = false;
                cpt++;
                if(k == 1){
                    i = 0;
                    it2 = this->nucleoBegin();
                }
                else{
                    i = (int) gsl_ran_flat (this->rng(), 1, k);
                    it2 = this->nucleosomes(this->nucleoBegin(), 0, i);
                }
                this->setTB(i);


                // go to the position i-1 exists
                // because  i > 0
                if(i>0){
                    it1 = it2;
                    it1--;
                    muBef = (*it1)->mu();
                }
                else
                {
                    muBef = this->minPos();
                }

                if( i < (k-1) ){ // Update nucleo after
                    it3 = it2;
                    it3++; // go to the position i+1
                    muNext = (*it3)->mu();
                    vNext = muNext;
                }
                else{

                    vNext = this->maxPos();
                    muNext = vNext + 1;
                }

                muMH = gsl_ran_flat(this->rng(), muBef, vNext); // New mu
                int df = (int) gsl_ran_flat (this->rng(), 3, this->dfMax() + 1);

                uMH = new  NucleoD(muMH, df, this->segSeq(), this->rng());

                aF = gsl_ran_flat(this->rng(), muBef, muMH);
                aR = gsl_ran_flat(this->rng(), muMH, muNext);

                flag = !(setNucleoD(uMH, aF, aR));

                if(!(flag)){
                    if(i > 0){
                        uBef = new  NucleoD(muBef, (*it1)->df(), this->segSeq(), this->rng()); // New i-1

                        // Enough reads between a and (*it1)->aR()
                        //flag = !(setNucleoDR(uBef, (*it1)->aF(), aF, *it1));

                        setNucleoDR(uBef, (*it1)->aF(), aF, *it1);
                    }

                    //if(!(flag) && i < (k-1)){ // Modify nucleosome i+1
                    if(i < (k-1)){ // Modify nucleosome i+1

                        uNext = new  NucleoD(muNext, (*it3)->df(), this->segSeq(), this->rng()); // New i+1
                        // Enough reads between a and (*it3)->aR()
                        //flag = !(setNucleoDR(uNext, aR, (*it3)->aR(), *it3));

                        setNucleoDR(uNext, aR, (*it3)->aR(), *it3);
                    }
                }
                if(flag){
                    delete uBef;
                    uBef = NULL;
                    delete uNext;
                    uNext = NULL;
                }

            }while(flag && cpt < 1000);

            if(!(flag))
            {
                this->setQalloc(vNext - muBef);

                this->pushModNucleo(*it2);
                *it2 = uMH;
                this->pushAddNucleo(it2);

                if(i > 0){

                    this->pushModNucleo(*it1);
                    *it1 = uBef;
                    this->pushAddNucleo(it1);

                }

                if(i < (k-1)){

                    this->pushModNucleo(*it3);
                    *it3 = uNext;
                    this->pushAddNucleo(it3);
                }

            }
        }
        catch(std::bad_alloc&) {
            Rcpp::stop("Memory problem\n");
        }
        return(!(flag));
    }

    template <typename NucleoD>
    void PartitionAll<NucleoD>::prepSpace(){
        this->evalW();
        this->evalKdDim();
        this->evalPriorMuDensity();
        this->evalMultinomial();

    }

    template <typename NucleoD>
    void PartitionAll<NucleoD>::delCurrent(){

        this->delCurrentD();
        this->resetAdd();
        this->resetMod();
        //this->resetNucleo();
        reset();
    }

    template <typename NucleoD>
    void PartitionAll<NucleoD>::delMod(){
        this->delCurrentD();
        this->resetAdd();
        this->resetMod();
        //reset();
    }

    template <typename NucleoD>
    void PartitionAll<NucleoD>::reset(){
        delete d_y;
        d_y = NULL;
        this->resetNucleo();
        /* delete les nucleosomes */

    }

    template <typename NucleoD>
    void PartitionAll<NucleoD>::reject(){
        this->resetAdd();
    }

    template <typename NucleoD>
    int PartitionAll<NucleoD>::getLimit(double start, double end
            , std::vector<double>::const_iterator &startIt
            , std::vector<double>::const_iterator &endIt
            , long &l, bool excEnd){

        std::vector<double>::const_iterator it=yBegin();
        bool flag=1;
        int cpt = 0;
        double pr = -1.0;
        double e = -0.000001;
        if(excEnd)
            e *= -1;

        l = 0;
        while(flag && it!=yEnd()){
            if(*it >= start){
                startIt = it;

                while(end - *it   > e  && it != yEnd()){
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
    }


} /* namespace space_process */

#endif /* PARTITIONALL_H_ */
