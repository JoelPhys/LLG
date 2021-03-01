#ifndef __CUMALLOC_H__
#define __CUMALLOC_H__

#include "../inc/array.h"

namespace cuglob {

	extern double *dSx1d, *dSy1d, *dSz1d;
	extern double *Hapx, *Hapy, *Hapz;
	extern double *dJx, *dJy, *dJz;
	extern int *dx_adj, *dadjncy;
	extern double *dtfa;
	extern Array<double> pJx, pJy, pJz;
	extern Array<int> px_adj, padjncy;

	void device_info();
	void allocate_heun_memory();
	void copy_field_to_device();
	void allocate_Jij_memory();
	void copy_jij_to_device();
	void allocate_device_consts();
	void copy_spins_to_device();
	void copy_spins_to_host();
	void clear_memory();
	void copy_thermal_to_device(double Thermal_Fluct);


	inline void check_cuda_errors(const char *filename, const int line_number)
	{
#ifdef DEBUG
		cudaThreadSynchronize();
		cudaError_t error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
			exit(-1);
		}
#endif
	};


}

#endif
