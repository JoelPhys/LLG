// cpp header files
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>

// my header files
#include "../inc/config.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/cuheun.h"
#include "../inc/cumalloc.h"
#include "../inc/cudefine.h"
#include "../inc/cufields.h"
#include "../inc/cuthermal.h"
#include "../inc/neighbourlist.h"

namespace cufields {

	__global__ void uniform(int N, double x, double y, double z, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			Hapx[i] = x;
			Hapy[i] = y;
			Hapz[i] = z;  
		}

	}

	__global__ void uniform_staggered(int nsites, int *dsublat_sites, int N, double x, double y, double z, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			
			int sublatsites = dsublat_sites[i % nsites];
			
			if (sublatsites == 0){
				Hapx[i] = x;
				Hapy[i] = y;
				Hapz[i] = z; 
			} 
			else if (sublatsites == 1){
				Hapx[i] = -1*x;
				Hapy[i] = -1*y;
				Hapz[i] = -1*z; 
			} 
		}

	}

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

		__global__ void square_pulse_staggered(int nsites, int *dsublat_sites, int N, double time, double start_time, double end_time, double height, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){

			if ((time >= start_time) && (time < end_time)){
					//if (( i % 4 == 0) || (i % 4 == 2)) {
				int sublatsites = dsublat_sites[i % nsites];
			
				if (sublatsites == 0){
						Hapx[i] = 0.0;
						Hapy[i] = 0.0;
						Hapz[i] = height;  
					}
				else if (sublatsites == 1){
						Hapx[i] = 0.0;
						Hapy[i] = 0.0;
						Hapz[i] = -1.0 * height;  
					}
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
		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev)));

		if (i < N){

			Hapx[i] = 0.0;
			Hapy[i] = gauss;
			Hapz[i] = 0.0; 

		}


	}

		__global__ void gaussian_pulse_staggered(int nsites, int *dsublat_sites, int N, double time, double height, double std_dev, double centre_pos,  double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev)));

		if (i < N){
			
			int sublatsites = dsublat_sites[i % nsites]; 
			
			if (sublatsites == 0) {
				Hapx[i] = 0.0;
				Hapy[i] = gauss;
				Hapz[i] = 0.0;  
			}
			else if (sublatsites == 1) {
				Hapx[i] = 0.0;
				Hapy[i] = -1.0 * gauss;
				Hapz[i] = 0.0;  
			}
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

	__global__ void multi_cycle_pulse_staggered(int nsites, int *dsublat_sites, int N, double time, double height, double std_dev, double centre_pos, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev))) * sin(2.0*M_PI*freq*(time - centre_pos));

		if (i < N){
 
			int sublatsites = dsublat_sites[i % nsites];
			
			if (sublatsites == 0){
				Hapx[i] = cos(0.25*M_PI)*gauss;
				Hapy[i] = cos(0.25*M_PI)*gauss;
				Hapz[i] = 0.0;  
			}
			else if (sublatsites == 1){
				Hapx[i] = -1.0 * cos(0.25*M_PI)*gauss;
				Hapy[i] = -1.0 * cos(0.25*M_PI)*gauss;
				Hapz[i] = 0.0;  
			}
		}


	}

	__global__ void sine_pulse(int N, double time, double height, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * sin(2.0*M_PI*freq*time);

		if (i < N){
			Hapx[i] = 0.0;
			Hapy[i] = gauss;
			Hapz[i] = 0.0;  
		}


	}

	__global__ void sine_pulse_linear(int N, double time, double height, double freq, double kpoint, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			double gauss1, gauss2;
			double sitepos = i; //  Nq*((i/(Nq*Lz*Ly)) % Lx)+qlayer;
			gauss1 = height * sin(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
			gauss2 = height * cos(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
			Hapx[i] = gauss1;
			Hapy[i] = 0.0;
			Hapz[i] = 0.0;  
		}


	}
	
	__global__ void sine_pulse_circular(int N, double time, double height, double freq, double kpoint, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			double gauss1, gauss2;
			double sitepos = i; //  Nq*((i/(Nq*Lz*Ly)) % Lx)+qlayer;
			gauss1 = height * sin(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
			gauss2 = height * cos(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
			Hapx[i] = gauss2;
			Hapy[i] = gauss1;
			Hapz[i] = 0.0;  
		}


	}


	__global__ void sine_pulse_staggered(int nsites, int *dsublat_sites, int N, double time, double height, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * sin(2.0*M_PI*freq*time);

		if (i < N){
			int sublatsites = dsublat_sites[i % nsites];

			if (sublatsites == 0){
				Hapx[i] = 0.0;
				Hapy[i] = gauss;
				Hapz[i] = 0.0;  
			}
			else if (sublatsites == 1){
				Hapx[i] = 0.0;
				Hapy[i] = -1.0 * gauss;
				Hapz[i] = 0.0;  
			}
		}


	}

	void testing(int i){

		Array<double> testingx;
		Array<double> testingy, testingz;
		testingx.resize(params::Nspins);
		// testingy.resize(params::Nspins);
		// testingz.resize(params::Nspins);

		CUDA_CALL(cudaMemcpy(testingx.ptr(), cuglob::Hapy, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		// CUDA_CALL(cudaMemcpy(testingy.ptr(), Hapy, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));
		// CUDA_CALL(cudaMemcpy(testingz.ptr(), Hapz, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));

	    std::cout << i << " ";
	    for (int a = 0; a < 4; a++){
		    std::cout << testingx(a) << " ";	
	    }
	    std::cout << std::endl;
	}

}
