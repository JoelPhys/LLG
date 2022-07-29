
#ifndef _THERMAL_H_
#define _THERMAL_H_

// cpp header files
#include <vector>

namespace thermal {

    // Temperature
    extern std::string temptype;
    extern double ttm_start;
    extern double temp_gradient;


	
	// Define global variables
    extern std::vector<double> ttime;     // time
    extern std::vector<double> tz;        // z plane
    extern std::vector<double> te;        // electron temp
    extern std::vector<double> pp;        // pump power
    extern std::vector<double> tp;        // phonon temp

    extern int nsteps;
    extern int nz;
	
	void initthermal(double temp);
	void ttm(double time);
	void ttmtofile();
    void closettmfile();
	//void ReadThermalFile();
}
#endif
