
#ifndef _PARAMS1_H_
#define _PARAMS1_H_

    #include <vector>
    #include <string>
    #include <cmath>
    #include <string>
    #include "array2d.h"
    #include "array.h"
    #include "libconfig.h++"


    namespace params {

        extern double k_B, mu_b, gamma;

        extern double dt, Nt, dtau, half_dtau;   
        extern int relaxtime;

        extern double lambda, lambdaPrime, mu_s, INVmu_s, d_z, thermal_const, d_z_prime;

        extern Array<double> H_appx;
        extern Array<double> H_appy;
        extern Array<double> H_appz;

        extern int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat, NmomentsSubLat;
        extern int Idx, Idy, Idz; // For integer lattice
        extern double a1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;

        extern int xdimS, ydimS, zdimS, start;
        extern double dt_spinwaves;
        
        
        // Specifies how the sublattices will be sorted in output
        extern std::string afmflag;
        
        //output file location
        extern std::string filepath;

        // Jij SETTINGS
        extern std::string Jij_filename;
        extern std::string Jij_units;
        extern std::string changesign;
        extern bool JijCutoff, Jijhalf;
        extern double Jij_min;
        extern int ibtoq;

        //atom sites
        extern double Plat[3][3];
        
        extern std::vector< std::vector<double> > sites;
        extern std::vector< std::vector<double> > Isites; // For integer lattice
        extern std::vector< std::vector<double> > initm; // For integer lattice

        extern double PlatINV[3][3];

        extern libconfig::Config cfg;

        void intitialiseConfig(const char* cfg_filename);
        void readparams();

    }

#endif
