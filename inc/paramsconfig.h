#ifndef __PARAMS_H__
#define __PARAMS_H__

    #include <vector>
    #include <string>
    #include <cmath>
    #include <libconfig.h++>


    namespace config {

        // PHYSICAL CONSTANTS
        const double k_B = 1.3807e-23;
        const double mu_b = 9.2740e-24;
        const double gamma = 1.76e11;

        // SIMULATION PARAMETERS
        const double dt = 1e-16; 
        const double Nt = 350000; 
        const double dtau = gamma * dt;
        const double half_dtau = dtau * 0.5;    

        // MATERIAL CONSTANTS
        const double lambda = 1;
        const double lambdaPrime = 1 / (1+(lambda*lambda));
        const double mu_s = 1.5 * mu_b;
        const double INVmu_s = 1 / mu_s;
        const double d_z = 0;
        const double thermal_const = sqrt( (2 * lambda * k_B)  / (mu_s * dtau) );
        const double d_z_prime = 2 * ( d_z / mu_s );


        // EXTERNAL FIELD
        const double H_app[3] = {0,0,0};

        // SYSTEM DIMENSIONS
        const int Lx = 10;
        const int Ly = 20;
        const int Lz = 40;
        const int Nq = 1;
        const double a1 = 0.287e-9;
        const int ax = 2;
        const int ay = 2;
        const int az = 2;
        const int Nspins = Nq*Lx*Ly*Lz;
        const int Nmoments = (1*Lx*Ly*Lz); 
        const double NsitesINV_S = 1/(Lx*Ly*Lz);
        const double xdim = ax*Lx;
        const double ydim = ay*Ly;
        const double zdim = az*Lz;
        const double NsitesINV = 1/(xdim*ydim*zdim);
        const int zdimC = zdim/2+1;


        // playing with spinwaves
        int xdimS = Lx;
        int ydimS = Ly;
        int zdimS = Lz/2 + 1;
        double dt_spinwaves = 1e-14;
        int start = 50000;

        // Jij SETTINGS
        std::string Jij_filename = "../SimpleCrystal_3D/SC_test1.dat";
        bool Jijhalf = false;
        std::string Jij_units = "J";
        bool JijCutoff = false;
        const double Jij_min = 0.01;
        const int ibtoq = -1;
        bool changesign = false;

        //atom sites
        double Plat[3][3] = {{1.0, 0.0, 0.0},
                            { 0.0, 1.0, 0.0},
                            { 0.0, 0.0, 1.0}};
        
        double sites[1][3] = {0.0, 0.0, 0.0};

        double PlatINV[3][3];

    }


#endif