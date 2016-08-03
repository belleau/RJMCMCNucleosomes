/*
 * Nucleosome.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: belleau
 */

#include "Nucleosome.h"
using namespace std;

namespace space_process {


Nucleosome::Nucleosome(double pos) {
	// TODO Auto-generated constructor stub
	d_mu = pos;
}

Nucleosome::~Nucleosome() {
	// TODO Auto-generated destructor stub
}

void Nucleosome::setFStartPos(std::vector<double>::iterator fStart, std::vector<double>::iterator fEnd, int n){
	d_startF = fStart;
	d_endF = fEnd;
	d_sizeF = n;
}

void Nucleosome::setRStartPos(std::vector<double>::iterator rStart, std::vector<double>::iterator rEnd, int n){
	d_startR = rStart;
	d_endR = rEnd;
	d_sizeR = n;
}

void Nucleosome::setStartF(std::vector<double>::iterator startF){
	d_startF = startF;
}

std::vector<double>::iterator Nucleosome::startF(){
	return(d_startF);
}

void Nucleosome::setEndF(std::vector<double>::iterator endF){
	d_endF = endF;
}

std::vector<double>::iterator Nucleosome::endF(){
	return(d_endF);
}

void Nucleosome::setStartR(std::vector<double>::iterator startR){
	d_startR = startR;
}

std::vector<double>::iterator Nucleosome::startR(){
	return(d_startR);
}

void Nucleosome::setEndR(std::vector<double>::iterator endR){
	d_endR = endR;
}

std::vector<double>::iterator Nucleosome::endR(){
	return(d_endR);
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

double Nucleosome::varRead(std::vector<double>::iterator start, std::vector<double>::iterator end, int n){
	double var = -1.0;
	if(n>0){
		double avg = accumulate(start, end, 0.0) / n;
		double sq_sum = 0.0;

        for(std::vector<double>::iterator it = start; it != end;it++){
            sq_sum = (*it - avg) * (*it - avg);
        }
		var = sq_sum / n;
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
