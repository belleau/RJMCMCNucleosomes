/*
 * SimulationNucleo.tpp
 *
 *  Created on: Sep 27, 2016
 *      Author: belleau
 */


//using namespace Rcpp;
//using namespace std;


template< typename NucleoSpace>
SimulationNucleo<NucleoSpace>::SimulationNucleo(SegmentSeq const &segSeq,
		gsl_rng * rng, long nbIteration):
		d_segSeq(segSeq), d_rng(rng),
		d_nbIterations(nbIteration),
		d_currentState(NULL), d_kMaxS(0){

}

//template< typename NucleoSpace> SimulationNucleo<NucleoSpace>::

template< typename NucleoSpace>
SimulationNucleo<NucleoSpace>::~SimulationNucleo(){

}

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::setCurrentState(NucleoSpace *currentState){
	d_currentState = currentState;
}


template< typename NucleoSpace>
SegmentSeq const &SimulationNucleo<NucleoSpace>::segSeq(){
	return(d_segSeq);
};

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::pushState(){
	d_results.push_back(d_currentState);
}

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::setRhoP1(double rhoP1){
	d_rhoP1 = rhoP1;
}

template< typename NucleoSpace>
long SimulationNucleo<NucleoSpace>::nbIterations(){
	return(d_nbIterations);
}

template< typename NucleoSpace>
long SimulationNucleo<NucleoSpace>::kMaxS(){
	return(d_kMaxS);
}

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::setKMaxS(long kMaxS){
	d_kMaxS = kMaxS;
}

template< typename NucleoSpace>
int SimulationNucleo<NucleoSpace>::sizeState(){
	return(d_results.size());
}

template< typename NucleoSpace>
double SimulationNucleo<NucleoSpace>::rhoP1(){
	return(d_rhoP1);
}

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::currentClone(){
	d_mod = (*currentState()).clone();
}

template< typename NucleoSpace>
void SimulationNucleo<NucleoSpace>::acceptMod(){
	d_currentState = d_mod;
	(*currentState()).addIteration();
}


