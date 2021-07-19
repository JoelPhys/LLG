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

		double gauss;
		gauss = height * exp(-1 * (((time - centre_pos) * (time - centre_pos))/(2 * std_dev * std_dev)));

		if (i < N){

			Hapx[i] = gauss;
			Hapy[i] = 0.0;
			Hapz[i] = 0.0; 

		}


	}

		__global__ void multi_cycle_pulse(int N, double time, double height, double std_dev, double centre_pos, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * exp(-1 * (((time - centre_pos) * (time - centre_pos))/(2 * std_dev * std_dev))) * sin(2*M_PI*freq*(time - centre_pos));

		if (i < N){

			Hapx[i] = gauss;
			Hapy[i] = 0.0;
			Hapz[i] = 0.0; 

		}


	}

}
