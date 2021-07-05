#ifndef _CUFUNCS_H_
#define _CUFUNCS_H_

#include <cuda.h> 
#include <curand.h>
#include <cuda_runtime.h>

namespace cufuncs {

	extern int threadsperbock;
	extern int bpg;
	void cuRotation();
	void cuDomainWall();
	void cuSquarePulse(double time);
	void cuGaussPulse(double time);
	void cuMultiPulse(double time);
	void init_device_vars();
	void integration(double time);
}
#endif
