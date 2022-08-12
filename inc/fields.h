#ifndef _FIELDS_H_
#define _FIELDS_H_

// cpp header files
#include <string>

// my header files
#include "../inc/array.h"

namespace fields {

	extern Array<double> H_appx;
	extern Array<double> H_appy;
	extern Array<double> H_appz; 

	extern double start_time;
	extern double end_time;
	extern double height;
	extern double std_dev;
	extern double centre_pos;
	extern double freq;
	extern double kpoint;
	extern std::string type;
	extern double cuniform[3];

	void readfields();
	void square_pulse(double time);
	void gaussian_pulse(double time);
	void multi_cycle_pulse(double time);
	void calculate(double time);
}
#endif
