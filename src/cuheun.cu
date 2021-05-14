#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include "../inc/cumalloc.h"
#include "../inc/cudefine.h"
#include "../inc/cuthermal.h"
#include "../inc/params1.h"
#include "../inc/NeighbourList.h"
#include "../inc/mathfuncs.h"

namespace cuheun {

	double *Sdashnx, *Sdashny, *Sdashnz;
	double *DelSx, *DelSy, *DelSz;
	double *Htx, *Hty, *Htz;  

	__constant__ double c_dzp;
	__constant__ double c_dtau; 
	__constant__ double c_hdtau; 
	__constant__ double c_lambda;
	__constant__ double c_lambdap;
	__constant__ int c_Nq;
	__constant__ double c_angle;

	void allocate_heun_consts(){
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_lambda), &params::lambda, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_lambdap), &params::lambdaPrime, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dtau), &params::dtau, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dzp), &params::d_z_prime, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_hdtau), &params::half_dtau, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_Nq), &params::Nq, sizeof(int)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_angle), &params::angle, sizeof(double)));
	}



	__global__ void cuRotfun(int N, double *dSx1d, double  *dSy1d, double  *dSz1d){
	
		const int i = blockDim.x*blockIdx.x + threadIdx.x;
	
       		if (i < N){

		double vec[3];
		if (( i % c_Nq == 0) || (i % c_Nq == 3) || (i % c_Nq == 5) || (i % c_Nq == 6)) {

                        vec[0] = dSx1d[i];
                        vec[1] = dSy1d[i] * cos(c_angle) - dSz1d[i] * sin(c_angle);
                        vec[2] = dSy1d[i] * sin(c_angle) + dSz1d[i] * cos(c_angle);


			dSx1d[i] = vec[0];
                        dSy1d[i] = vec[1];
                        dSz1d[i] = vec[2];


		}
		}
	}

	__global__ void cuHeun1(int N, double *Thermal_Fluct, float *gvalsx1, float *gvalsy1, float *gvalsz1, int *dx_adj1, int *dadjncy1, double *Htx, double *Hty, double *Htz, double *dSx1d, double *dSy1d, double *dSz1d, double *dJx, double *dJy, double *dJz, double *Hapx, double *Hapy, double *Hapz, double *DelSx,  double *DelSy, double *DelSz, double *Sdashnx, double *Sdashny, double *Sdashnz){

		const int a = blockDim.x*blockIdx.x + threadIdx.x;

		if (a < N){
			Htx[a] = static_cast<double>(gvalsx1[a]) * Thermal_Fluct[a];
			Hty[a] = static_cast<double>(gvalsy1[a]) * Thermal_Fluct[a];
			Htz[a] = static_cast<double>(gvalsz1[a]) * Thermal_Fluct[a];

			double Han[3] = {0.0, 0.0, c_dzp * dSz1d[a]};
			double Hex[3] = {0.0, 0.0, 0.0};
			int counting = dx_adj1[a];

			for (int b = dx_adj1[a]; b < dx_adj1[a+1]; b++){
				Hex[0] += dJx[counting] * (dSx1d[dadjncy1[b]]);
				Hex[1] += dJy[counting] * (dSy1d[dadjncy1[b]]);
				Hex[2] += dJz[counting] * (dSz1d[dadjncy1[b]]);
				counting++;
			}

			double Hnew[3];
			Hnew[0] = Htx[a] + Hapx[a] + Han[0] + Hex[0];
			Hnew[1] = Hty[a] + Hapy[a] + Han[1] + Hex[1];
			Hnew[2] = Htz[a] + Hapz[a] + Han[2] + Hex[2];

			double SxH[3];
			SxH[0] = dSy1d[a] * Hnew[2] - dSz1d[a] * Hnew[1];
			SxH[1] = dSz1d[a] * Hnew[0] - dSx1d[a] * Hnew[2];
			SxH[2] = dSx1d[a] * Hnew[1] - dSy1d[a] * Hnew[0];

			double SxSxH[3];
			SxSxH[0] = dSy1d[a] * SxH[2] - dSz1d[a] * SxH[1];
			SxSxH[1] = dSz1d[a] * SxH[0] - dSx1d[a] * SxH[2];
			SxSxH[2] = dSx1d[a] * SxH[1] - dSy1d[a] * SxH[0];

			DelSx[a] = - 1 * c_lambdap * (SxH[0] + c_lambda * SxSxH[0]);
			DelSy[a] = - 1 * c_lambdap * (SxH[1] + c_lambda * SxSxH[1]);
			DelSz[a] = - 1 * c_lambdap * (SxH[2] + c_lambda * SxSxH[2]);

			double Sdash[3];
			Sdash[0] = dSx1d[a] + (DelSx[a] * c_dtau);
			Sdash[1] = dSy1d[a] + (DelSy[a] * c_dtau);
			Sdash[2] = dSz1d[a] + (DelSz[a] * c_dtau);

			double cinvmag = 1 / sqrt(Sdash[0] * Sdash[0] + Sdash[1] * Sdash[1] + Sdash[2] * Sdash[2]);

			Sdashnx[a] = cinvmag * Sdash[0];
			Sdashny[a] = cinvmag * Sdash[1];
			Sdashnz[a] = cinvmag * Sdash[2];
		} 
	}

	__global__ void cuHeun2(int N, int *dx_adj1, int *dadjncy1, double *Htx, double *Hty, double *Htz, double *dSx1d, double *dSy1d, double *dSz1d, double *dJx, double *dJy, double *dJz, double *Hapx, double *Hapy, double *Hapz, double *DelSx,  double *DelSy, double *DelSz, double *Sdashnx, double *Sdashny, double *Sdashnz){

		const int a = blockDim.x * blockIdx.x + threadIdx.x;

		if (a < N){
			double Han_dash[3] = {0.0, 0.0, c_dzp * Sdashnz[a]};
			double Hex_dash[3] = {0.0, 0.0, 0.0};
			int counting = dx_adj1[a];

			// Exchange interaction prime
			for (int b = dx_adj1[a]; b < dx_adj1[a+1]; b++){

				Hex_dash[0] += dJx[counting] * (Sdashnx[dadjncy1[b]]);
				Hex_dash[1] += dJy[counting] * (Sdashny[dadjncy1[b]]);
				Hex_dash[2] += dJz[counting] * (Sdashnz[dadjncy1[b]]);
				counting++;
			}

			double Hnew_dash[3];
			Hnew_dash[0] = Htx[a] + Hapx[a] + Han_dash[0] + Hex_dash[0];
			Hnew_dash[1] = Hty[a] + Hapy[a] + Han_dash[1] + Hex_dash[1];
			Hnew_dash[2] = Htz[a] + Hapz[a] + Han_dash[2] + Hex_dash[2];

			double SxHd[3];
			SxHd[0] = Sdashny[a] * Hnew_dash[2] - Sdashnz[a] * Hnew_dash[1];
			SxHd[1] = Sdashnz[a] * Hnew_dash[0] - Sdashnx[a] * Hnew_dash[2];
			SxHd[2] = Sdashnx[a] * Hnew_dash[1] - Sdashny[a] * Hnew_dash[0];

			double SxSxHd[3];
			SxSxHd[0] = Sdashny[a] * SxHd[2] - Sdashnz[a] * SxHd[1];
			SxSxHd[1] = Sdashnz[a] * SxHd[0] - Sdashnx[a] * SxHd[2];
			SxSxHd[2] = Sdashnx[a] * SxHd[1] - Sdashny[a] * SxHd[0];

			double DelSd[3];
			DelSd[0] = -1 * c_lambdap * (SxHd[0] + c_lambda * SxSxHd[0]);
			DelSd[1] = -1 * c_lambdap * (SxHd[1] + c_lambda * SxSxHd[1]);
			DelSd[2] = -1 * c_lambdap * (SxHd[2] + c_lambda * SxSxHd[2]);

			double Sn[3];
			Sn[0] = dSx1d[a] + c_hdtau * (DelSx[a] + DelSd[0]);
			Sn[1] = dSy1d[a] + c_hdtau * (DelSy[a] + DelSd[1]);
			Sn[2] = dSz1d[a] + c_hdtau * (DelSz[a] + DelSd[2]);

			double cinvmag1 = 1 / sqrt(Sn[0] * Sn[0] + Sn[1] * Sn[1] + Sn[2] * Sn[2]);
			dSx1d[a] = cinvmag1 * Sn[0];
			dSy1d[a] = cinvmag1 * Sn[1];
			dSz1d[a] = cinvmag1 * Sn[2];
		}
	}

	void testing(){

		Array<double> testingx;
		Array<double> testingy, testingz;
		double sumx = 0;
		double sumy = 0;
		double sumz = 0;
		testingx.resize(params::Nspins);
		testingy.resize(params::Nspins);
		testingz.resize(params::Nspins);

		CUDA_CALL(cudaMemcpy(testingx.ptr(), Sdashnx, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(testingy.ptr(), Sdashny, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(testingz.ptr(), Sdashnz, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));

		for (int i = 0; i < params::Nspins; i++){
			sumx += testingx[i];
			sumy += testingy[i];
			sumz += testingz[i];
		}
		std::cout << sumx << " " << sumy << " " << sumz << std::endl;	
	}

}
