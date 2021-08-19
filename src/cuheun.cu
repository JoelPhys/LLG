#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include "../inc/cumalloc.h"
#include "../inc/cudefine.h"
#include "../inc/cuthermal.h"
#include "../inc/config.h"
#include "../inc/NeighbourList.h"
#include "../inc/mathfuncs.h"

namespace cuheun {

	double *Sdashnx, *Sdashny, *Sdashnz;
	double *DelSx, *DelSy, *DelSz;
	double *Htx, *Hty, *Htz;  

	__constant__ double c_dxup;	
	__constant__ double c_dyup;
	__constant__ double c_dzup;
	__constant__ double c_dxcp;
	__constant__ double c_dycp;
	__constant__ double c_dzcp;
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
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dxup), &params::dxup, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dyup), &params::dyup, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dzup), &params::dzup, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dxcp), &params::dxcp, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dycp), &params::dycp, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_dzcp), &params::dzcp, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_hdtau), &params::half_dtau, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_Nq), &params::Nq, sizeof(int)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_angle), &params::angle, sizeof(double)));
	}

	__global__ void cuFixSpins1(int N, int *dlw, int *drw,  double *dsurfx, double *dsurfy, double *dsurfz, double *dSx1d, double *dSy1d, double *dSz1d){

		const int a = blockDim.x * blockIdx.x + threadIdx.x;

		if (a < N){

			dSx1d[dlw[a]] =  dsurfx[a];
			dSy1d[dlw[a]] =  dsurfy[a];
			dSz1d[dlw[a]] =  dsurfz[a];

			// dSx1d[drw[a]] = -1.0;
			// dSy1d[drw[a]] =  0.0;
			// dSz1d[drw[a]] =  0.0;
		}

	}

	__global__ void cuFixSpins2(int N, int *dlw, int *drw, double *dsurfx, double *dsurfy, double *dsurfz, double *Sdashnx, double *Sdashny, double *Sdashnz){

		const int a = blockDim.x * blockIdx.x + threadIdx.x;

		if (a < N){

			Sdashnx[dlw[a]] =  dsurfx[a];
			Sdashny[dlw[a]] =  dsurfy[a];
			Sdashnz[dlw[a]] =  dsurfz[a];

			// Sdashnx[drw[a]] = -1.0;
			// Sdashny[drw[a]] =  0.0;
			// Sdashnz[drw[a]] =  0.0;
		}

	}

	__global__ void cuRotfun(int N, double *dSx1d, double  *dSy1d, double  *dSz1d){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){

			double vec[3];
			if (( i % c_Nq == 1) || (i % c_Nq == 3)) {

				vec[0] = dSx1d[i] * cos(c_angle) - dSy1d[i] * sin(c_angle);
				vec[1] = dSx1d[i] * sin(c_angle) + dSy1d[i] * cos(c_angle);
				vec[2] = dSz1d[i];

				dSx1d[i] = vec[0];
				dSy1d[i] = vec[1];
				dSz1d[i] = vec[2];

			}
		}
	}

	__global__ void cuHeun1(int N,  double time, double *Thermal_Fluct, float *gvalsx1, float *gvalsy1, float *gvalsz1, int *dx_adj1, int *dadjncy1, double *Htx, double *Hty, double *Htz, double *dSx1d, double *dSy1d, double *dSz1d, double *dJx, double *dJy, double *dJz, double *Hapx, double *Hapy, double *Hapz, double *DelSx,  double *DelSy, double *DelSz, double *Sdashnx, double *Sdashny, double *Sdashnz){

		const int a = blockDim.x*blockIdx.x + threadIdx.x;

		if (a < N){
			Htx[a] = static_cast<double>(gvalsx1[a]) * Thermal_Fluct[a];
			Hty[a] = static_cast<double>(gvalsy1[a]) * Thermal_Fluct[a];
			Htz[a] = static_cast<double>(gvalsz1[a]) * Thermal_Fluct[a];

			double Huni[3];
			Huni[0] = c_dxup * dSx1d[a]; 
			Huni[1] = c_dyup * dSy1d[a]; 
			Huni[2] = c_dzup * dSz1d[a]; 

			double Hcub[3];
			Hcub[0] = c_dxcp * dSx1d[a] * dSx1d[a] * dSx1d[a]; 
			Hcub[1] = c_dycp * dSy1d[a] * dSy1d[a] * dSy1d[a]; 
			Hcub[2] = c_dzcp * dSz1d[a] * dSz1d[a] * dSz1d[a];

			double Hex[3] = {0.0, 0.0, 0.0};
			int counting = dx_adj1[a];

			for (int b = dx_adj1[a]; b < dx_adj1[a+1]; b++){
				Hex[0] += dJx[counting] * (dSx1d[dadjncy1[b]]);
				Hex[1] += dJy[counting] * (dSy1d[dadjncy1[b]]);
				Hex[2] += dJz[counting] * (dSz1d[dadjncy1[b]]);
				counting++;
			}

			double Hnew[3];
			Hnew[0] = Htx[a] + Hapx[a] + Huni[0] + Hcub[0] + Hex[0];
			Hnew[1] = Hty[a] + Hapy[a] + Huni[1] + Hcub[1] + Hex[1];
			Hnew[2] = Htz[a] + Hapz[a] + Huni[2] + Hcub[2] + Hex[2];

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

	__global__ void cuHeun2(int N, double time, int *dx_adj1, int *dadjncy1, double *Htx, double *Hty, double *Htz, double *dSx1d, double *dSy1d, double *dSz1d, double *dJx, double *dJy, double *dJz, double *Hapx, double *Hapy, double *Hapz, double *DelSx,  double *DelSy, double *DelSz, double *Sdashnx, double *Sdashny, double *Sdashnz){

		const int a = blockDim.x * blockIdx.x + threadIdx.x;

		if (a < N){

			double Huni_dash[3];
			Huni_dash[0] = c_dxup * Sdashnx[a];
			Huni_dash[1] = c_dyup * Sdashny[a];
			Huni_dash[2]=  c_dzup * Sdashnz[a];

			double Hcub_dash[3];
			Hcub_dash[0] = c_dxcp * Sdashnx[a] * Sdashnx[a] * Sdashnx[a];
			Hcub_dash[1] = c_dycp * Sdashny[a] * Sdashny[a] * Sdashny[a];
			Hcub_dash[2]=  c_dzcp * Sdashnz[a] * Sdashnz[a] * Sdashnz[a];

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
			Hnew_dash[0] = Htx[a] + Hapx[a] + Huni_dash[0] + Hcub_dash[0] + Hex_dash[0];
			Hnew_dash[1] = Hty[a] + Hapy[a] + Huni_dash[1] + Hcub_dash[1] + Hex_dash[1];
			Hnew_dash[2] = Htz[a] + Hapz[a] + Huni_dash[2] + Hcub_dash[2] + Hex_dash[2];

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

			// if (a == 0){
			// 	double vec[3];
			// 	vec[0] = cos(0.00001 * time);
			// 	vec[1] = sin(0.00001 * time);
			// 	vec[2] = 0.0;

			// 	dSx1d[a] = vec[0];
			// 	dSy1d[a] = vec[1];
			// 	dSz1d[a] = vec[2];
			// }
		}
	}

	// void testing(int i){

	// 	Array<double> testingx;
	// 	Array<double> testingy, testingz;
	// 	testingx.resize(params::Nspins);
	// 	// testingy.resize(params::Nspins);
	// 	// testingz.resize(params::Nspins);

	// 	CUDA_CALL(cudaMemcpy(testingx.ptr(), Te, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));
	// 	// CUDA_CALL(cudaMemcpy(testingy.ptr(), Hapy, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));
	// 	// CUDA_CALL(cudaMemcpy(testingz.ptr(), Hapz, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));

    //     // std::cout << "STEP: " << i << std::endl;
    //     for (int a = 0; a < 1; a++){
	// 	    std::cout << testingx(a) << " ";	
    //     }
    //     // std::cout << std::endl;
	// }




}
