#ifndef _FIELDS_H_
#define _FIELDS_H_

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

	void readfields();

}
#endif
