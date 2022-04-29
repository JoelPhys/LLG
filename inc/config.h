
#ifndef _PARAMS1_H_
#define _PARAMS1_H_

    // cpp header files
    #include <vector>
    #include <string>

    // my header files
    #include "libconfig.h++"


    namespace params {

        extern std::string simtype;

        extern double k_B, mu_b, gamma;

        extern double dt, Nt, dtau, half_dtau;   
        extern int relaxtime, outputstep;

        extern std::vector<double> lambda;
        extern std::vector<double> lambdaPrime;
	    extern std::vector<double> thermal_const;
        extern std::vector<double> mu_s;
        extern std::vector<double> INVmu_s;

        // Uniaxial Anisotropy
        extern double dxu, dyu, dzu;
        extern double dxup, dyup, dzup;

        // Uniaxial Anisotropy
        extern double dzc, dzcp;

        extern int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat;
        extern int Idx, Idy, Idz; // For integer lattice
        extern double a1, b1, c1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;
        extern std::vector<int> sublat_sites;
        extern std::vector<int> NmomentsSubLat;

        extern int xdimS, ydimS, zdimS, start;
        extern double dt_spinwaves;
        extern double sg_spinwaves;

        //Rotation angle
        extern double angle;	
        
        // Specifies how the sublattices will be sorted in output
        extern std::string afmflag;

        // Input file format
        extern std::string format;
        
        //output file location
        extern std::string filepath;
        extern std::string filepath_sw;

        // Jij SETTINGS
        extern std::string Jij_filename;
        extern std::string Jij_units;
        extern std::string changesign;
        extern bool JijCutoff, Jijhalf;
        extern double Jij_min;
        extern int ibtoq;

        //atom sites
        extern double Plat[3][3];

        //Boundary Conditions
        extern std::string xbound;
        extern std::string ybound;
        extern std::string zbound;
        
        extern std::vector< std::vector<double> > sites;
        extern std::vector< std::vector<double> > Isites; // For integer lattice
        extern std::vector< std::vector<double> > initm; // For integer lattice

        // Temperature
        extern std::string temptype;
        extern double ttm_start;
        extern double temp_gradient;

        // Output Lattice
        extern bool OutputLattice;		
        extern int OutputLatticeStep;	
        extern std::string OutputLatticeFilepath;		

        extern bool OutputToTerminal;

        extern double PlatINV[3][3];

        extern libconfig::Config cfg;

        void banner();
        void intitialiseConfig(const char* cfg_filename);
        void readparams();

    }

#endif
