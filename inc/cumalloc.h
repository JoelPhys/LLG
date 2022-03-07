#ifndef __CUMALLOC_H__
#define __CUMALLOC_H__

#include "../inc/array.h"

namespace cuglob {

	extern double *dSx1d, *dSy1d, *dSz1d;
	extern double *Hapx, *Hapy, *Hapz;
	extern double *dJx, *dJy, *dJz;
	extern int *dx_adj, *dadjncy;
	extern int *dlw, *drw;
	extern Array<double> pJx, pJy, pJz;
	extern Array<int> px_adj, padjncy;
	extern int *dsimspin;

	//testing for hedgehog;
	extern double *dsurfx, *dsurfy, *dsurfz;

	// Damping
	extern double *c_lambda, *c_lambdap;

	//testing
	extern double *dJx_new, *dJy_new, *dJz_new;
	extern int *djind;

	void device_info();
	void allocate_heun_memory();
	void allocate_Jij_memory();
	void allocate_device_consts();
	void copy_temp_to_device(double equilibium_temp);
	void copy_field_to_device();
	void copy_spins_to_device();
	void copy_spins_to_host();
	void copy_field_to_host();
	void copy_damp_to_device();
	void copy_jij_to_device();
	void copy_dw_to_device();
	void clear_memory();


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
