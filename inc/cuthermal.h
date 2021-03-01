#ifndef __CUTHERMAL_H__
#define __CUTHERMAL_H__

#include <curand.h>
#include <cuda.h>
#include <cuda_runtime.h>

namespace cuthermal {

	extern float *gvalsx, *gvalsy, *gvalsz;

	void curand_generator();
	void gen_thermal_noise();
	void destroy_generator();
}

#endif
