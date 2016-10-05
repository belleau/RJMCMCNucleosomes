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
//		typedef typename std::vector<double> containerD;
//		typedef typename containerD::const_iterator iteratorD;
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
		SpaceNucleosomeD(SegmentSeq const &segSeq)
			:SpaceNucleosome<NucleoD>(segSeq), d_w(NULL), d_dim(NULL),
			 d_dfMax(30){
			setDefault();
		};

		SpaceNucleosomeD(SegmentSeq const &segSeq, int seed)
			:SpaceNucleosome<NucleoD>(segSeq, seed), d_w(NULL),
			 d_dim(NULL), d_dfMax(30){
			setDefault();
		};

		SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng)
			:SpaceNucleosome<NucleoD>(segSeq, rng), d_w(NULL), d_dim(NULL),
			 d_dfMax(30){
			setDefault();
		};

		SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng, int dfMax)
			:SpaceNucleosome<NucleoD>(segSeq, rng), d_w(NULL), d_dim(NULL),
			 d_dfMax(dfMax){
			setDefault();
		};

/*
		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta){
			setDefault();
		};

		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta, long sizeFReads, long sizeRReads)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta, sizeFReads, sizeRReads){
			setDefault();
		};
*/
		virtual ~SpaceNucleosomeD(){}; //delete[] d_w;
		double tB(){
			return(d_tB);
		};
		void setTB(double tB){
			d_tB = tB;
		}

		int lambda(){
			return(d_lambda);
		};

		void setLambda(int l){
			d_lambda = l;
		};

		int dfMax(){
			return(d_dfMax);
		};

		double qalloc(){
			return(d_qalloc);
		};

		void setQalloc(double qalloc){

			d_qalloc = qalloc;
		};

		double meanRead(){
			return(d_meanRead);
		};

		double r2()
		{
			return(d_r2);
		};

		void evalW(){
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
				std::cout << "Memory problem\n";
				std::cerr << "Memory problem\n";
			}

		};


		double cMuDensity(){
			return(d_cMuDensity);
		};

		void setPriorMuDensity(double priorMuDensity){
			d_priorMuDensity = priorMuDensity;
		};

		double priorMuDensity(){
			return(d_priorMuDensity);
		};

		void insertD(double mu, int df){
			NucleoD *u = new NucleoD(mu, df, this->segSeq(), this->rng());
			this->insert(u);
		};

		void evalPriorMuDensity(){
			/* Matrix  omega
			 * R <- max(readPositions) - min(readPositions)
			 * tau <- 1/R^2
			 * E <- (max(readPositions) + min(readPositions))/2
			 * M <- rep(E, k)
			 * const <- (pi/(2*tau))^{-k/2}
			 *
			 * equation 11
			 * const * exp(-(tau/2) * (t(mu - M) %*% omega %*% (mu - M)))
			 *
			 * */
			//std::cout.precision(17);
			double m = meanRead(); /* Mean of the read*/
			double result = 0;
			int k = this->valK();

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
					double v [3] = {(**nucleoIt++).mu() - m, (**nucleoIt++).mu()  - m, 0};
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
		};

		void evalDim(){

			try{
				d_dim = new unsigned int[this->valK()];

				int n = 0;
				for(itNucleo it = this->nucleoBegin(); it != this->nucleoEnd(); it++)
				{

					d_dim[n++] = (*it)->dimN();
				}
			}
			catch(std::bad_alloc&) {
				std::cout << "Memory problem\n";
				std::cerr << "Memory problem\n";
			}
		}

		void evalKdDim(){
			/*std::cout << "sizeF " << this->sizeFReads() << "\n";
			std::cout << "sizeR " << this->sizeRReads() << "\n";*/

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
				std::cout << "Memory problem\n";
				std::cerr << "Memory problem\n";
			}
		};

		void evalKdDim1(){
					/*std::cout << "sizeF " << this->sizeFReads() << "\n";
					std::cout << "sizeR " << this->sizeRReads() << "\n";*/
                    //d_tB = 0;
					try{
						d_dim = new unsigned int[this->valK()];
						int i;
						int s = this->sizeFReads() + this->sizeRReads();

						double *yRead = new double[s];
						double *pourv = new double[s];

						std::fill_n(yRead, s, 0.0);
						std::fill_n(pourv, s, 0.0);

						int n = 0;
						for(itNucleo it = this->nucleoBegin(); it != this->nucleoEnd(); it++)
						{
							i = 0;
							d_dim[n] = (*it)->dimN();//(*it)->sizeF() + (*it)->sizeR();
							for(std::vector<double>::const_iterator  itF = (*it)->bFBegin(); itF != (*it)->bFEnd(); itF++){
								pourv[i] += (*itF);
								yRead[i] += d_w[n] * (*itF);
								i++;
							}
							for(std::vector<double>::const_iterator  itR = (*it)->bRBegin(); itR != (*it)->bREnd(); itR++){
								pourv[i] += (*itR);
								yRead[i] += d_w[n] * (*itR);
								i++;
							}

							n++;

						}

						d_kD = 0;
						int bla = 0;
						for(int j = 0; j < s; j++){
							if(j%2 == 0){
								//d_tB += pourv[j];
								//d_tB += log(yRead[j]);
								bla++;
								/*if(bla == 5){
									std::cout << "yf5 " << log(yRead[j]) << "\n";
								}*/
							}
							d_kD += log(yRead[j]);
						}
						//d_tB /= bla;

						delete[] yRead;
						delete[] pourv;
					}
					catch(std::bad_alloc&) {
						std::cout << "Memory problem\n";
						std::cerr << "Memory problem\n";
					}
				};

		void evalMultinomial(){
			d_multinomial = gsl_ran_multinomial_pdf (this->valK(), d_w, d_dim);
		};

		double multinomial(){
			return(d_multinomial);
		};

		double dK(){
			double d = 0;
			if(this->valK() > 1){
				d = 0.5 * std::min(1.0, this->valK() / ((double) d_lambda));
			}
			return(d);
		};

		double bK(){
			double d = 0;
			d = 0.5 * std::min(1.0,  d_lambda / (double)(this->valK() + 1.0) );
			return(d);
		};
		double kD(){
			return(d_kD);
		}
		double rhoP1(){

			return(priorMuDensity() * multinomial());
		};

		double rhoP2(){
			//std::cout << "dK " << dK() << " " << lambda() << " " << this->valK() << " " << qalloc() << "\n";
			return(dK() * (lambda() / ((double) this->valK())) );
		};

	protected:
		void setMeanRead(double meanRead){
			d_meanRead = meanRead;
		};

		void setR2(double r2){
			d_r2 = r2;
		};

		void setCMuDensity(double cMuDensity){
			d_cMuDensity = cMuDensity;
		};

		void delCurrentD(){
			delete[] d_dim;
			d_dim = NULL;
			delete[] d_w;
			d_w = NULL;
		};

	private:
		void setDefault(){
			d_r2 = pow((this->maxPos() - this->minPos()),2);
			d_meanRead = (this->maxPos() + this->minPos()) / 2.0;
			d_cMuDensity = M_PI * d_r2 / 2.0;
			d_lambda = 3;
		};

	};
}

#endif /* SPACENUCLEOSOMED_H_ */
