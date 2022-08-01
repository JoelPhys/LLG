
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
	
	// Two-temperature model global variables
	extern double gamma_e;           
    extern double Cp;                
    extern double kappa_0;           
    extern double delta;             
    extern double Gep;               
    extern double P_0;               
    extern double t0;                
    extern double tau;               
	extern double oneOvrdzdz;
    extern double oneOvr2dz;
	
	// functions
	void initthermal(double temp);
	void ttm(double time);
	void ttmtofile(double time);
    void closettmfile();
	//void ReadThermalFile();
}
#endif
