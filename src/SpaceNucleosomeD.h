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
#include <gsl/gsl_blas.h>
#include "SpaceNucleosome.h"
#include "NucleoDirichlet.h"

namespace space_process{

	template<typename NucleoD>    /***** BEWARE NucleoD Must inherit from SpaceNucleosomeD *****/
	class SpaceNucleosomeD: public SpaceNucleosome<NucleoD>{
		typedef std::list<NucleoD> containerNucleo;
		typedef typename containerNucleo::iterator itNucleo;
		std::vector<double> d_w;
		double d_kD;
		double d_priorMuDensity;

		double d_meanRead, d_r2, d_cMuDensity;


	public:

		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta){
			setDefault();
		};

		SpaceNucleosomeD(std::vector<double> const  &fReads, std::vector<double> const &rReads, int zeta, long sizeFReads, long sizeRReads)
			:SpaceNucleosome<NucleoD>(fReads, rReads, zeta, sizeFReads, sizeRReads){
			setDefault();
		};

		virtual ~SpaceNucleosomeD(){};

		double meanRead(){
			return(d_meanRead);
		};

		double r2()
		{
			return(d_r2);
		};

		double cMuDensity(){
			return(d_cMuDensity);
		};

		void setPriorMuDensity(double priorMuDensity){
			d_priorMuDensity = priorMuDensity;
		}

		void insertD(double mu, int df){
			NucleoD u(mu, df, this->segSeq(), this->rng());
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




            //nucleoIt = this->nucleoBegin();

			itNucleo nucleoIt = this->nucleoBegin();
			result = 0;
			if(this->valK() == 1){
				result = 2 * pow((*nucleoIt).mu()  - m, 2);
			}
			else{
				if(this->valK() == 2){
					double v [2]= {(*nucleoIt++).mu()  - m, (*nucleoIt).mu() - m};
					result = (2 * v[0] -  v[1]) * v[0] + (v[1] - v[0]) *v[1];
				}
				if(this->valK() > 2)
				{
					double v [3] = {(*nucleoIt++).mu() - m, (*nucleoIt++).mu()  - m, 0};
					result = (2 * v[0] -  v[1]) * v[0];

					int i = 2;
					do{
						v[i%3] =  (*nucleoIt++).mu() - m;

						result += (2 * v[(i-1)%3] - v[(i)%3] - v[(i-2)%3]) * v[(i-1)%3];// verifier le modulo
						i++;
					}while(nucleoIt != this->nucleoEnd());
					// result du dernier


					result +=  (v[(i-1)%3] - v[(i-2)%3]) * v[(i-1)%3];

				}
			}

			double tmpC = pow(cMuDensity(), -1 * this->valK() / 2);


			setPriorMuDensity(tmpC * exp(- result/ (2 * r2()) ));
		};

		/* *******************************************
		 *
		 * Version gsl of evalPriorMuDensity
		 * evalGSLPriorMuDensity
		 * Unit: milliseconds
		 * min       lq     mean   median       uq      max neval
		 * 53.57089 55.60392 58.43092 57.84131 60.53918 73.7484   200
		 * evalPriorMuDensity
		 * 2.370297 2.539221 2.754932 2.667075 2.85793 4.749791   200
		 *
		 ********************************************  */
		double evalGSLPriorMuDensity(){
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


			/* Version gsl slower */
			for(int a = 0; a < 100000; a++){
			itNucleo nucleoIt = this->nucleoBegin();


			double *omegaD = new double[k * k]();
			double *vD = new double[k]();



			for(int i = 0; i < k; i++){
				omegaD[i * k + i] = 2.0;
				vD[i] = (*nucleoIt++).mu() - m;
				//gsl_matrix_set(omega, i, i, 2.0);
				//gsl_matrix_set(v, 1, i, (*nucleoIt).mu() - m);
				//nucleoIt++;
				if(i > 0){
					omegaD[i * k +i-1] = -1.0;
					omegaD[(i-1) * k + i] = -1.0;
					//gsl_matrix_set(omega, i, i-1, -1.0);

					//gsl_matrix_set(omega, i-1, i, -1.0);
				}
			}
			omegaD[(k-1) * k + (k-1)] = 1.0;

			double *cD = new double[k]();
			double *dD = new double[1]();
			gsl_matrix_view omega = gsl_matrix_view_array(omegaD, k, k);
			gsl_matrix_view v = gsl_matrix_view_array(vD, k, 1);
			gsl_matrix_view c = gsl_matrix_view_array(cD, 1, k);
			gsl_matrix_view d = gsl_matrix_view_array(dD, 1, 1);
			gsl_blas_dgemm (CblasTrans, CblasNoTrans,
										  1.0, &v.matrix, &omega.matrix,
										  0.0, &c.matrix);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
													  1.0, &c.matrix, &v.matrix,
													  0.0, &d.matrix);

			delete[] omegaD;
			delete[] vD;
			delete[] cD;
			delete[] dD;
			}


			//double tmpC = pow(cMuDensity(), -1 * this->valK() / 2);
			//tmpC * exp(- result/ (2 * r2()) )
			return(result);
		};

	private:
		void setDefault(){
			d_r2 = pow((this->maxPos() - this->minPos()),2);
			d_meanRead = (this->maxPos() + this->minPos()) / 2;
			d_cMuDensity = M_PI * d_r2 / 2;
		}

	};
}

#endif /* SPACENUCLEOSOMED_H_ */
