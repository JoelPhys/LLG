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
#include "../inc/fields.h"
#include "../inc/cumalloc.h"
#include "../inc/cudefine.h"
#include "../inc/cufields.h"
#include "../inc/cuthermal.h"
#include "../inc/neighbourlist.h"

namespace cufields {


	// device constant variables
	__constant__ double c_direc_mag;
	__constant__ int c_direcx;
	__constant__ int c_direcy;
	__constant__ int c_direcz;
	int* d_sublat_stag;

	void allocate_field_variables(){

		int test1,test2,test3;
		test1=fields::direc[0];
		test2=fields::direc[1];
		test3=fields::direc[2];
		
		// variables in constant memory
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_direc_mag), &fields::direc_mag, sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_direcx), &test1, sizeof(int)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_direcy), &test2, sizeof(int)));
		CUDA_CALL(cudaMemcpyToSymbol(*(&c_direcz), &test3, sizeof(int)));
		

		// variables in global memory
		CUDA_CALL(cudaMalloc((void**)&d_sublat_stag, sizeof(int)*params::Nsublat));
		CUDA_CALL(cudaMemset(d_sublat_stag, 0.0, sizeof(int)*params::Nsublat));
		CUDA_CALL(cudaMemcpy(d_sublat_stag,  &fields::sublat_stag[0], sizeof(int)*params::Nsublat, cudaMemcpyHostToDevice));

	}


	//// pointer function for cuda fields: input is effective field and sublattice, time, 
	//__device__ void (*cuda_field_function_pointer)(double *, double * double *);

	//// function for each field
	//__device__ void uniform(double *eff_field, double *sublat, double temp){
	//	eff_fieldx = Bx;
	//	eff_fieldy = By;
	//	eff_fieldz = z;	
	//}	
	//__device__ void uniform_staggered(){


	//	if (i < N){
	//		
	//		int sublatsites = dsublat_sites[i % nsites];
	//		
	//		if (sublatsites == 0){
	//			Hapx[i] = x;
	//			Hapy[i] = y;
	//			Hapz[i] = z; 
	//		} 
	//		else if (sublatsites == 1){
	//			Hapx[i] = -1*x;
	//			Hapy[i] = -1*y;
	//			Hapz[i] = -1*z; 
	//		} 
	//	}

	//}

	__global__ void uniform(int nsites, int *dsublat_sites, int *d_sublat_stag, int N, double x, double y, double z, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			int sublatsites = dsublat_sites[i % nsites];
			Hapx[i] = d_sublat_stag[sublatsites]*x;
			Hapy[i] = d_sublat_stag[sublatsites]*y;
			Hapz[i] = d_sublat_stag[sublatsites]*z;  
		}

	}

	__global__ void square_pulse(int nsites, int *dsublat_sites, int *d_sublat_stag, int N, double time, double start_time, double end_time, double height, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){

			if ((time >= start_time) && (time < end_time)){

				int sublatsites = dsublat_sites[i % nsites];

				Hapx[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcx*height;
				Hapy[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcy*height;
				Hapz[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcz*height;
			}
			else {
				Hapx[i] = 0.0;
				Hapy[i] = 0.0;
				Hapz[i] = 0.0; 
			}

		}

	}

	__global__ void gaussian_pulse(int nsites, int *dsublat_sites, int *d_sublat_stag, int N, double time, double height, double std_dev, double centre_pos,  double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev)));

		if (i < N){

			int sublatsites = dsublat_sites[i % nsites];

			Hapx[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcx*gauss;
			Hapy[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcy*gauss;
			Hapz[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcz*gauss;

		}


	}

	__global__ void multi_cycle_pulse(int nsites, int *dsublat_sites, int *d_sublat_stag, int N, double time, double height, double std_dev, double centre_pos, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev))) * sin(2.0*M_PI*freq*(time - centre_pos));

		if (i < N){

			int sublatsites = dsublat_sites[i % nsites];

			Hapx[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcx*gauss;
			Hapy[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcy*gauss;
			Hapz[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcz*gauss;

		}


	}


	// TODO: fix this section, not sure how to code up circular and linear fields at present
	__global__ void sine_pulse_linear(int N, double time, double height, double freq, double kpoint, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			double gauss1;
			double sitepos = i; //  Nq*((i/(Nq*Lz*Ly)) % Lx)+qlayer;
			gauss1 = height * sin(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
			Hapx[i] = gauss1;
			Hapy[i] = 0.0;
			Hapz[i] = 0.0;  
		}


	}

	// TODO: fix this section, not sure how to code up circular and linear fields at present
	__global__ void sine_pulse_circular(int N, double time, double height, double freq, double kpoint, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		if (i < N){
			double gauss1 = 0;
			double gauss2 = 0;
			double sitepos = i; //  Nq*((i/(Nq*Lz*Ly)) % Lx)+qlayer;
			double kstep;
			for (int k = 0; k < N; k++){
				kstep = static_cast<double>(k)/static_cast<double>(N);	
				gauss1 = height * sin(kstep * M_PI * sitepos + 2.0*M_PI*freq*time);
				gauss2 = height * cos(kstep * M_PI * sitepos + 2.0*M_PI*freq*time);
			}
			Hapx[i] = gauss2;
			Hapy[i] = gauss1;
			Hapz[i] = 0.0;  
		}

		//if (i < N){
		//	double gauss1, gauss2;
		//	double sitepos = i; //  Nq*((i/(Nq*Lz*Ly)) % Lx)+qlayer;
		//	gauss1 = height * sin(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
		//	gauss2 = height * cos(kpoint * M_PI * sitepos + 2.0*M_PI*freq*time);
		//	Hapx[i] = gauss2;
		//	Hapy[i] = gauss1;
		//	Hapz[i] = 0.0;  
		//}


	}


	__global__ void sine_pulse(int nsites, int *dsublat_sites, int *d_sublat_stag,  int N, double time, double height, double freq, double *Hapx, double *Hapy, double *Hapz){

		const int i = blockDim.x*blockIdx.x + threadIdx.x;

		double gauss;
		gauss = height * sin(2.0*M_PI*freq*time);

		if (i < N){

			int sublatsites = dsublat_sites[i % nsites];

			Hapx[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcx*gauss;
			Hapy[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcy*gauss;
			Hapz[i] = d_sublat_stag[sublatsites]*c_direc_mag*c_direcz*gauss;

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
