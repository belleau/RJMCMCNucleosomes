/*
 * Nucleosome.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#include "Nucleosome.h"
using namespace std;

namespace space_process {


Nucleosome::Nucleosome(double pos, SegmentSeq const &segSeq, gsl_rng *rng)
    :d_segSeq(segSeq){
    d_mu = pos;
    d_rng = rng;
}

Nucleosome::~Nucleosome() {
}

void Nucleosome::setFStartPos(std::vector<double>::const_iterator fStart, std::vector<double>::const_iterator fEnd, int n){
    d_startF = fStart;
    d_endF = fEnd;
    d_sizeF = n;
}

void Nucleosome::setRStartPos(std::vector<double>::const_iterator rStart, std::vector<double>::const_iterator rEnd, int n){
    d_startR = rStart;
    d_endR = rEnd;
    d_sizeR = n;
}

void Nucleosome::setStartF(std::vector<double>::const_iterator startF){
    d_startF = startF;
}

vector<double>::const_iterator Nucleosome::startF(){
    return(d_startF);
}

void Nucleosome::setEndF(std::vector<double>::const_iterator endF){
    d_endF = endF;
}

std::vector<double>::const_iterator Nucleosome::endF(){
    return(d_endF);
}

void Nucleosome::setStartR(std::vector<double>::const_iterator startR){
    d_startR = startR;
}

std::vector<double>::const_iterator Nucleosome::startR(){
    return(d_startR);
}

void Nucleosome::setEndR(std::vector<double>::const_iterator endR){
    d_endR = endR;
}

std::vector<double>::const_iterator Nucleosome::endR(){
    return(d_endR);
}

void Nucleosome::setDimN(long dimN){
    d_dimN = dimN;
}

long Nucleosome::dimN(){
    return(d_dimN);
}


void Nucleosome::setSizeF(int sizeF){
    d_sizeF = sizeF;
}

int Nucleosome::sizeF(){
    return(d_sizeF);
}

void Nucleosome::setSizeR(int sizeR){
    d_sizeR = sizeR;
}

int Nucleosome::sizeR(){
    return(d_sizeR);
}

double Nucleosome::mu(){
    return(d_mu);
}

void Nucleosome::setSigmaF(double sigmaF){
    d_sigmaF = sigmaF;
}

double Nucleosome::sigmaF(){
    return(d_sigmaF);
}

void Nucleosome::setSigmaR(double sigmaR){
    d_sigmaR = sigmaR;
}

double Nucleosome::sigmaR(){
    return(d_sigmaR);
}

double Nucleosome::deltaMin(){
    return(double(d_segSeq.deltaMin()));
}
double Nucleosome::deltaMax(){
    return(double(d_segSeq.deltaMax()));
}

double Nucleosome::zeta(){
    return(double(d_segSeq.zeta()));
}

vector<double>::const_iterator Nucleosome::beginFR(){
    return(d_segSeq.beginFR());
}

vector<double>::const_iterator Nucleosome::endFR(){
    return(d_segSeq.endFR());
}

long Nucleosome::sizeFR(){
    return(d_segSeq.sizeFReads());
}

vector<double>::const_iterator Nucleosome::beginRR(){
    return(d_segSeq.beginRR());
}

vector<double>::const_iterator Nucleosome::endRR(){
    return(d_segSeq.endRR());
}

long Nucleosome::sizeRR(){
    return(d_segSeq.sizeRReads());
}

double Nucleosome::varRead(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end, int n){
    double var = -1.0;
    int pv = 0;
    if(n>0){
        double avg = accumulate(start, end, 0.0) / n;

        double sq_sum = 0.0;
        //std::vector<double>::const_iterator tmp = end;
        //tmp--;

        for(std::vector<double>::const_iterator it = start; it != end;it++){
            sq_sum += (*it - avg) * (*it - avg);
            pv++;
        }
        var = sq_sum / (n - 1);
    }

    return(var);
}

void Nucleosome::evalSigmaF(){
    d_sigmaF = varRead(startF(), endF(), sizeF());
}

void Nucleosome::evalSigmaR(){
    d_sigmaR = varRead(startR(), endR(), sizeR());
}


} /* namespace space_process */
