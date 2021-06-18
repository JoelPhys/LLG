#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include "../inc/array.h"
#include "../inc/NeighbourList.h"
#include "../inc/params1.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cudefine.h"
#include "../inc/cufields.h"

namespace cufields {

	__global__ void square_pulse(int N, double time, double start_time, double end_time, double height, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){

			if ((time >= start_time) && (time < end_time)){
				Hapx[i] = height;
				Hapy[i] = 0.0;
				Hapz[i] = 0.0;  
			}
			else {
				Hapx[i] = 0.0;
				Hapy[i] = 0.0;
				Hapz[i] = 0.0; 
			}

		}

	}

	__global__ void gaussian_pulse(int N, double time, double height, double std_dev, double centre_pos,  double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;



		if (i < N){

			double gauss;
			gauss = height * exp(-1 * (((i - centre_pos) * (i - centre_pos))/(2 * std_dev * std_dev)));

			Hapx[i] = gauss;
			Hapy[i] = 0.0;
			Hapz[i] = 0.0; 

		}


	}

	void testing(){

		Array<double> testingx;
		Array<double> testingy, testingz;
		testingx.resize(params::Nspins);
		testingy.resize(params::Nspins);
		testingz.resize(params::Nspins);

		CUDA_CALL(cudaMemcpy(testingx.ptr(), cuglob::Hapx, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(testingy.ptr(), cuglob::Hapy, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(testingz.ptr(), cuglob::Hapz, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));

		std::cout << testingx(0) << " " << testingy(0) << " " << testingz(0) << std::endl;	
	}

}
