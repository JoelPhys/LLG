
#ifndef __THERMAL_H__
#define __THERMAL_H__

#include <fstream>
#include <vector>
#include <iostream>
#include "../inc/config.h"
#include "../inc/libconfig.h++"
#include "../inc/util.h"


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