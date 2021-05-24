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

	void readfields();

}
#endif
