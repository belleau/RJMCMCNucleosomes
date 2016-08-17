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
		typedef std::list<NucleoD*> containerNucleo;
		typedef typename containerNucleo::const_iterator itNucleo;
		typedef typename std::vector<double> containerD;
		typedef typename containerD::const_iterator iteratorD;
		double *d_w;
		double d_kD;
		double d_priorMuDensity;
		double d_multinomial;

		int d_lambda;
		int *d_dim;

		double d_meanRead, d_r2, d_cMuDensity;


	public:
		SpaceNucleosomeD(SegmentSeq const &segSeq)
			:SpaceNucleosome<NucleoD>(segSeq){
			setDefault();
		};

		SpaceNucleosomeD(SegmentSeq const &segSeq, int seed)
			:SpaceNucleosome<NucleoD>(segSeq, seed){
		};

		SpaceNucleosomeD(SegmentSeq const &segSeq, gsl_rng * rng)
			:SpaceNucleosome<NucleoD>(segSeq, rng){

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

		int lambda(){
			return(d_lambda);
		};

		void setLambda(int l){
			d_lambda = l;
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

			double *alpha = new double[k];
			//memset(alpha, 1.0, k);
			std::fill_n(alpha, k, 1.0);
			double *theta = new double[k];
			gsl_ran_dirichlet (this->rng(), k, alpha, theta);

			d_w = theta;
			delete[] alpha;
		};


		double cMuDensity(){
			return(d_cMuDensity);
		};

		void setPriorMuDensity(double priorMuDensity){
			d_priorMuDensity = priorMuDensity;
		}

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

			double tmpC = pow(cMuDensity(), -1 * this->valK() / 2);


			setPriorMuDensity(tmpC * exp(- result/ (2 * r2()) ));
		};


		void evalKdDim(){
			/*std::cout << "sizeF " << this->sizeFReads() << "\n";
			std::cout << "sizeR " << this->sizeRReads() << "\n";*/
			d_dim = new int[this->valK()];
			int i;
			int s = this->sizeFReads() + this->sizeRReads();
			double *yRead = new double[s];
			//double *yR = new double[this->sizeRReads()];
			std::fill_n(yRead, s, 0.0);
			//std::fill_n(yR, this->sizeRReads(), 0.0);
			int n = 0;
			for(itNucleo it = this->nucleoBegin(); it != this->nucleoEnd(); it++)
			{
				i = 0;
				d_dim[n] = (*it)->sizeF() + (*it)->sizeR();
				for(iteratorD  itF = (*it)->bFBegin(); itF != (*it)->bFEnd(); itF++){
					yRead[i++] += d_w[n] * (*itF);
				}
				for(iteratorD  itF = (*it)->bRBegin(); itF != (*it)->bREnd(); itF++){
					yRead[i++] += d_w[n++] * (*itF);
				}

			}
			d_kD = 0;
			for(int j = 0; j < s; j++){
				d_kD += log(yRead[j]);
			}

            //delete[] yF;
            delete[] yRead;
		}

		void evalMultinomial(){
			/* double gsl_ran_multinomial_pdf (size_t K, const double p[], const unsigned int n[])*/
			d_multinomial = gsl_ran_multinomial_pdf (this->valK(), d_w, d_dim);
		}

		double dK(){
			double d = 0;
			if(this->valK() > 1){
				d = 0.5 * std::min(1, this->valK() / d_lambda);
			}
			return(d);
		}
		double bK(){
			double d = 0;
			if(this->valK() > 1){
				d = 0.5 * std::min(1, d_lambda / (this->valK() + 1) );
			}
			return(d);
		}

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
	private:
		void setDefault(){
			d_r2 = pow((this->maxPos() - this->minPos()),2);
			d_meanRead = (this->maxPos() + this->minPos()) / 2;
			d_cMuDensity = M_PI * d_r2 / 2;
			d_lambda = 3;
		}

	};
}

#endif /* SPACENUCLEOSOMED_H_ */
