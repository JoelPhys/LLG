#ifndef _neighbourlist_H_
#define _neighbourlist_H_

    #include <iostream>
    #include <sstream>
    #include <cmath>
    #include <fstream>
    #include <vector>
    #include <iomanip>
    #include <random>
    #include <algorithm>
    #include "config.h"
    #include "array.h"
    #include "array.h"
    #include "array2d.h"
    #include "array3d.h"

    namespace neigh {

        // Globals
        extern std::vector<unsigned int> adjncy;
        extern std::vector<double> Jijy_prime;
        extern std::vector<double> Jijz_prime;
        extern std::vector<double> Jijx_prime;
        extern std::vector<unsigned int> x_adj;

        // Functions
        void init();
        void ReadFile();
        void InteractionMatrix();
        void Heun(double Thermal_Fluct);

    }

#endif
