#ifndef _CUINTEGRATE_H_
#define _CUINTEGRATE_H_

#include <cuda.h> 
#include <curand.h>
#include <cuda_runtime.h>

namespace cuint {

	extern int threadsperbock;
	extern int bpg;
	void cuRotation();
	void init_device_vars();
	void integration();
}
#endif
