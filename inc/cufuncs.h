#ifndef _CUFUNCS_H_
#define _CUFUNCS_H_

#include <cuda.h> 
#include <curand.h>
#include <cuda_runtime.h>
#include <sstream>

namespace cufuncs {

	extern int threadsperbock;
	extern int bpg;
	void cuRotation();
	void cuDomainWall();
	void cuSquarePulse(double time, double start_time, double end_time, double height);
	void cuGaussPulse(double time);
	void cuMultiPulse(double time);
	void cuTemperature(std::string type, double time, double ttm_start);
	void init_device_vars();
	void integration(double time);
}
#endif
