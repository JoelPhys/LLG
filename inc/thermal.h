
#ifndef _THERMAL_H_
#define _THERMAL_H_

// cpp header files
#include <vector>

namespace thermal {

    // Define global variables
    extern std::vector<double> ttime;     // time
    extern std::vector<double> tz;        // z plane
    extern std::vector<double> te;        // electron temp
    extern std::vector<double> pp;        // pump power
    extern std::vector<double> tp;        // phonon temp

    extern int nsteps;
    extern int nz;

    void ReadThermalFile();

}
#endif