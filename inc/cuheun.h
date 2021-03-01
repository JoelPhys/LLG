#ifndef __CUHEUN_H__
#define __CUHEUN_H__

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include "../inc/params1.h"
#include "../inc/cumalloc.h"


namespace cuheun {

	extern double *Sdashnx, *Sdashny, *Sdashnz;
	extern double *DelSx, *DelSy, *DelSz;
	extern double *Htx, *Hty, *Htz;

	void allocate_heun_consts();
	extern __global__ void cuHeun1(int, double *, float *, float *, float *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,  double *, double *, double *, double *, double *);
	extern __global__ void cuHeun2(int , int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,  double *, double *, double *, double *, double *);
	
	//For Debug
	void testing();
}

#endif
