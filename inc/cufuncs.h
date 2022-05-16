#ifndef _CUFUNCS_H_
#define _CUFUNCS_H_


#include <sstream>

namespace cufuncs {

	void cuRotation();
	void cuDomainWall();
	void cuFields(std::string type, double time, double start_time, double end_time, double height);
	void cuTemperature(std::string type, double time, double ttm_start);
	void init_device_vars();
	void integration(double time);
	void initial_energy(double time);
	void integrationDW(double time);
}
#endif
