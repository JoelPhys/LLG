#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sstream>
#include "../inc/mathfuncs.h"
#include "../inc/array2d.h"
#include "../inc/params1.h"
#include "../inc/error.h"

namespace params {

    libconfig::Config cfg;

    double k_B, mu_b, gamma;
    double dt, Nt, dtau, half_dtau;   
    double lambda, lambdaPrime, mu_s, INVmu_s, d_z, thermal_const, d_z_prime;
    int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat, NmomentsSubLat;
    int Idx, Idy, Idz; // For integer lattice
    double a1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;
    int xdimS, ydimS, zdimS, start;
    double dt_spinwaves;

    int relaxtime;

    // Jij SETTINGS
    std::string Jij_filename;
    std::string Jij_units;
    bool JijCutoff, changesign, Jijhalf;
    double Jij_min;
    int ibtoq;

    // Applied field settings
    std::string Type;
    bool Hfield;
    Array2D<double> H_app;

    // Lattive Vectors
    double Plat[3][3];
    double PlatINV[3][3];

    std::vector< std::vector<double> > sites;
    std::vector< std::vector<double> > Isites;

    //Intialise Config File ====================================================================================================================================//
    void intitialiseConfig(const char* cfg_filename){

        if(!cfg_filename)
        {
            error::Info(__FILE__,__LINE__);
            error::Message("You must give a config file, exiting");
        }
        try
        {
            cfg.readFile(cfg_filename);
        }
        catch(const libconfig::FileIOException &fioex)
        {
            error::Info(__FILE__,__LINE__);
            error::Message("I/O error while reading config file");
        }
        catch(const libconfig::ParseException &pex)
        {
            error::Info(__FILE__,__LINE__);
            std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
            exit(EXIT_FAILURE);
        }
        cfg.setAutoConvert(true);
    }
    //===========================================================================================================================================================//

    // Read Parameters ==========================================================================================================================================//
    void readparams(){

        k_B = cfg.lookup("PhysicalConsts.BoltzmannConstant");
        mu_b = cfg.lookup("PhysicalConsts.BohrMagneton");
        gamma = cfg.lookup("PhysicalConsts.GyromagneticRatio");

        dt = cfg.lookup("Time.SizeOfStep");

        Nt = cfg.lookup("Time.NumberOfSteps");
        relaxtime = cfg.lookup("Time.RelaxationTime");
        dtau = gamma * dt;
        half_dtau = 0.5 * dtau;   

        lambda = cfg.lookup("MaterialConsts.lambda");
        mu_s = cfg.lookup("MaterialConsts.mu_s");
        d_z = cfg.lookup("MaterialConsts.d_z");
        a1 = cfg.lookup("MaterialConsts.a");
        mu_s *= mu_b;
        INVmu_s = 1 / mu_s;
        thermal_const = sqrt( (2 * lambda * k_B)  / (mu_s * dtau) );
        d_z_prime = 2 * ( d_z / mu_s );
        lambdaPrime = 1 / (1+(lambda*lambda));

        // system dimensions
        Lx = cfg.lookup("Geom.UnitCellsInX");
        Ly = cfg.lookup("Geom.UnitCellsInY");
        Lz = cfg.lookup("Geom.UnitCellsInZ");
        Nq = cfg.lookup("Geom.NumberOfSites");
        Nsublat = cfg.lookup("Geom.NumberOfSublat");
        ax = 2;
        ay = 2;
        az = 2;

        Nspins = Nq*Lx*Ly*Lz;
        Nmoments = (Nq*Lx*Ly*Lz); 
        NmomentsSubLat = (Nsublat*Lx*Ly*Lz);
        NsitesINV_S = 1/(Lx*Ly*Lz);
        xdim = ax*Lx;
        ydim = ay*Ly;
        zdim = az*Lz;
        NsitesINV = 1/(xdim*ydim*zdim);
        zdimC = zdim/2+1;

        //Read Site positions ==============================================================================
        sites.resize(Nq);
        for (int s = 0; s < Nq; s++){
            sites[s].resize(3);

            std::stringstream sstr;
            sstr << "Site" << s;
            std::string str = sstr.str();

            libconfig::Setting& setting = cfg.lookup("Sites");  
            sites[s][0] = setting[str.c_str()][0];
            sites[s][1] = setting[str.c_str()][1];
            sites[s][2] = setting[str.c_str()][2];
        }
        //=======================================================================================================
        
        // Read integer site position ===========================================================================
        Isites.resize(Nq);
        std::cout << "Isites = " << std::endl;
        for (int s = 0; s < Nq; s++){
            Isites[s].resize(3);

            std::stringstream sstr;
            sstr << "ISite" << s;
            std::string str = sstr.str();

            libconfig::Setting& setting = cfg.lookup("IntegerSites");  
            Isites[s][0] = setting[str.c_str()][0];
            Isites[s][1] = setting[str.c_str()][1];
            Isites[s][2] = setting[str.c_str()][2];
        }
        //Read integer lattice spacing
        Idx = cfg.lookup("IntegerSites.Idx");
        Idy = cfg.lookup("IntegerSites.Idy");
        Idz = cfg.lookup("IntegerSites.Idz");

        //=======================================================================================================

        //Read Lattice Vectors
        std::cout << "Lattice Vectors = " << std::endl;
        for (int v = 0; v < 3; v++){

            std::stringstream sstr1;
            sstr1 << "LatVec" << v;
            std::string str1 = sstr1.str();

            libconfig::Setting& setting = cfg.lookup("LatticeVectors");   
            Plat[v][0] = setting[str1.c_str()][0];
            Plat[v][1] = setting[str1.c_str()][1];
            Plat[v][2] = setting[str1.c_str()][2];
            std::cout << Plat[v][0] << " " << Plat[v][1] << " " << Plat[v][2] << std::endl;
        }


       
        //Read external field
        H_app.resize(Nmoments,3);
        H_app.IFill(0);
        Type = cfg.lookup("ExternalField.Type").c_str(); 
        libconfig::Setting& setting1 = cfg.lookup("ExternalField");

        if (Type == "Uniform") {     
            for (int a = 0; a < Nmoments; a++){     
                H_app(a,0) = setting1["Field"][0]; 
                H_app(a,1) = setting1["Field"][1]; 
                H_app(a,2) = setting1["Field"][2];
            }    
        }
        else if (Type == "Split") {   
            for (int a = 0; a < Nmoments; a++){     
                if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 2)) {
                    H_app(a,0) = setting1["Field"][0]; 
                    H_app(a,1) = setting1["Field"][1]; 
                    H_app(a,2) = setting1["Field"][2];
                }
                else if ((modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
                    H_app(a,0) = -1.0 * static_cast<double>(setting1["Field"][0]); 
                    H_app(a,1) = -1.0 * static_cast<double>(setting1["Field"][1]); 
                    H_app(a,2) = -1.0 * static_cast<double>(setting1["Field"][2]);
                }
                else std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
            }
        }
        else std::cout << "WARNING: Unknown Field Type." << std::endl;

        start = cfg.lookup("Spinwaves.StartTime");
        dt_spinwaves = cfg.lookup("Spinwaves.TimeStep");
        Jij_filename = cfg.lookup("Exchange.InputFile").c_str();        
        Jij_units = cfg.lookup("Exchange.Units").c_str();   
        JijCutoff = cfg.lookup("Exchange.Cutoff");    
        changesign = cfg.lookup("Exchange.ChangeSign");  
        Jijhalf = cfg.lookup("Exchange.Double_Jij");    
        Jij_min = cfg.lookup("Exchange.CutoffEnergy");    
        ibtoq = cfg.lookup("Exchange.ibtoq");    

    }
    //========================================================================================================================================//


}
