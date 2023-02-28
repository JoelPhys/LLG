#ifndef __CUFIELDS_H__
#define __CUFIELDS_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>


namespace cufields {

    extern double start_time;
    extern double end_time;
    extern double height;
	extern int *d_sublat_stag;

	void allocate_field_variables();
    extern __global__ void uniform(int, int *, int *, int, double, double, double, double *, double *, double *);
    extern __global__ void square_pulse(int, int*, int*, int, double, double, double, double, double *, double *, double *);
	extern __global__ void gaussian_pulse(int, int *, int *, int, double, double, double, double,  double *, double *, double *);
    extern __global__ void multi_cycle_pulse(int, int*, int*, int , double, double, double, double,  double, double *, double *, double *);
    extern __global__ void sine_pulse(int , int *, int *, int, double , double , double , double *, double *, double *);
	extern __global__ void sine_pulse_circular(int , double , double , double , double , double *, double *, double *, int, int *,int *);
	extern __global__ void sine_pulse_linear(int , double , double , double , double , double *, double *, double *, int, int *,int *);

    void testing(int i);
}

#endif
