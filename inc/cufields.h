#ifndef __CUFIELDS_H__
#define __CUFIELDS_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>


namespace cufields {

    extern double start_time;
    extern double end_time;
    extern double height;

    extern __global__ void uniform(int, double, double, double, double *, double *, double *);
    extern __global__ void uniform_staggered(int, int *, int, double, double, double, double *, double *, double *);
    extern __global__ void square_pulse(int, double, double, double, double, double *, double *, double *);
    extern __global__ void square_pulse_staggered(int, int *, int, double, double, double, double, double *, double *, double *);
    extern __global__ void gaussian_pulse(int, double, double, double, double,  double *, double *, double *);
    extern __global__ void gaussian_pulse_staggered(int, int *, int, double, double, double, double,  double *, double *, double *);
    extern __global__ void multi_cycle_pulse(int, double, double, double, double,  double, double *, double *, double *);
    extern __global__ void multi_cycle_pulse_staggered(int, int *, int, double, double, double, double,  double, double *, double *, double *);
    extern __global__ void sine_pulse(int , double , double , double , double *, double *, double *);
    extern __global__ void sine_pulse_staggered(int, int *, int , double , double , double , double *, double *, double *);

    void testing(int i);
}

#endif
