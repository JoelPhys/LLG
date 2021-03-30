#ifndef __NEIGHBOURLIST_H__
#define __NEIGHBOURLIST_H__

    #include <iostream>
    #include <sstream>
    #include <cmath>
    #include <fstream>
    #include <vector>
    #include <iomanip>
    #include <random>
    #include <algorithm>
    // #include "params.h"
    #include "array.h"
    #include "array2d.h"
    #include "array3d.h"

    namespace neigh {

        // Globals ================================================================= //
        extern Array2D<double> H_thermal;
        extern Array2D<double> Delta_S;
        extern Array<double> Sx1d;
        extern Array<double> Sy1d;
        extern Array<double> Sz1d;
        extern Array<double> S_dash_normedx1d;
        extern Array<double> S_dash_normedy1d;
        extern Array<double> S_dash_normedz1d;
        extern std::vector<unsigned int> adjncy;
        extern std::vector<double> Jijy_prime;
        extern std::vector<double> Jijz_prime;
        extern std::vector<double> Jijx_prime;
        extern std::vector<unsigned int> x_adj;
        // ======================================================================== //

        void IntialisePointersNL();
        void ReadFile();
        void InteractionMatrix();
        void Heun(double Thermal_Fluct);

    }

#endif
