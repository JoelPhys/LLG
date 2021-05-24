#ifndef __CUFIELDS_H__
#define __CUFIELDS_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>


namespace cufields {

    extern double start_time;
    extern double end_time;
    extern double height;

    extern __global__ void square_pulse(int, double, double, double, double, double *, double *, double *);
    void testing();

}

#endif