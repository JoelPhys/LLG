#ifndef __CUTHERMAL_H__
#define __CUTHERMAL_H__

#include <curand.h>
#include <cuda.h>
#include <cuda_runtime.h>

namespace cuthermal {

	extern float *gvalsx, *gvalsy, *gvalsz;
	extern double *Te, *Tp, *P_it;
	extern double *dtfa;
	extern int *dzlayer;

	void init_cuthermal(double equilibrium_temp);
	void curand_generator();
	void gen_thermal_noise();
	void destroy_generator();
	__global__ void ttm(double , int , double *, double *, double *);
	__global__ void ttf(double , int , double *, double *, int *);

	// debugging
	void testing(int i);
}

#endif
