#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <time.h>
#include "../inc/fftw3.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/mathfuncs.h"
#include "../inc/params1.h"
#include "../inc/geom.h"
#include "../inc/NeighbourList.h"
#include "../inc/util.h"
// #include "inc/FFT.h"
#include "../inc/spinwaves.h"

#define IMAG 1
#define REAL 0

int main(int argc, char* argv[]){

    // functions ============================================================================================== //
    params::intitialiseConfig(argv[1]); 
    params::readparams();
    geom::CreateLattice();
    geom::CountDistinct();
    geom::CreateIntLattice();
    neigh::ReadFile();
    neigh::InteractionMatrixJerome();
    neigh::IntialisePointersNL();
    util::InitUtil();
    IdentityMatrix();
    spinwaves::initialiseFFT();
    // ======================================================================================================== //

    // ======= Temperature ==================================================================================== //
    const double Temp = (atof(argv[2]));
    const double thermal_fluct = params::thermal_const * sqrt(Temp);
    std::cout << "Temperature = " << Temp << std::endl;
    // ========================================================================================================= //

    // ======== Set clock ====================================================================================== //
    clock_t begin, end;
    double time_spent;
    begin = clock();
    // ========================================================================================================= //
    
    if ((std::string(argv[3]) == "1") || (std::string(argv[3]) == "3")){
        // ======= Initiliase Spin Position ======================================================================== //
        Array4D<double> Sx4, Sy4, Sz4;
        Sx4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
        Sy4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
        Sz4.resize(params::Lx, params::Ly, params::Lz, params::Nq);

        Sx4.IFill(0);
        Sy4.IFill(0);
        Sz4.IFill(0);


        int count1d = 0;

        for (int x = 0; x < params::Lx; x++){
            for (int y = 0; y < params::Ly; y++){
                for (int z = 0; z < params::Lz; z++){

                    Sx4(x,y,z,0) = 1;
                    Sx4(x,y,z,1) = -1;
                    Sx4(x,y,z,2) = 1;
                    Sx4(x,y,z,3) = -1;

                    neigh::Sx1d(count1d + 0) = Sx4(x,y,z,0);
                    neigh::Sx1d(count1d + 1) = Sx4(x,y,z,1);
                    neigh::Sx1d(count1d + 2) = Sx4(x,y,z,2);
                    neigh::Sx1d(count1d + 3) = Sx4(x,y,z,3);

                    count1d += params::Nq; 
                }
            }
        }
        // ================================================================================================== //
    }
    if (std::string(argv[3]) == "1"){
         // Relax to equilibrium magnetisation (10 ps) ======================================================= //
        std::stringstream sstr_eq;
        sstr_eq << "output/magnetisation/equilibrium_T_" << Temp << ".txt";
        std::ofstream equilfile;
        equilfile.open(sstr_eq.str());


        for (int i = 0; i < params::relaxtime; i++){
            neigh::Heun(thermal_fluct);
        }

        for (int j = 0; j < params::Nmoments; j++){
            equilfile << neigh::Sx1d[j] << " " << neigh::Sy1d[j] << " " << neigh::Sz1d[j] << std::endl;
        }

        equilfile << std::flush;
        equilfile.close();

        end = clock();
        std::cout << std::setprecision(10) << "Equilibration time = " << (double)(end - begin) / CLOCKS_PER_SEC << " seconds" << std::endl; 
        // ================================================================================================== //
    }
    if (std::string(argv[3]) == "2"){

        // Read in equilibrium spin values ================================================================== //
        std::stringstream sstr_eq;
        sstr_eq << "output/magnetisation/equilibrium_T_" << Temp << ".txt";
        std::ifstream equilibrationfile(sstr_eq.str());
        
        if (!equilibrationfile){
            std::cout << "ERROR: Could not open equilibrium file" << std::endl;
            exit(0);
        }

        double sx, sy, sz;
        int i = 0;
        while (equilibrationfile >> sx >> sy >> sz){
            neigh::Sx1d[i] = sx;
            neigh::Sy1d[i] = sy;
            neigh::Sz1d[i] = sz;
            i++;
        }
        equilibrationfile.close();
        // ================================================================================================== //
    }
    if ((std::string(argv[3]) == "2") || (std::string(argv[3]) == "3")){
        // Initialise some variables ======================================================================== // 
        double t = 0;
        double tau = 0;

        int c;
        c = params::dt_spinwaves / params::dt;
 
        util::InitOutputFile(Temp);

        // ========== LOOP THROUGH TIMESTEPS ================================================================ //
        for (int i = 0; i < params::Nt; i++){
            
            t = t + params::dt;
            tau = tau + params::dtau;

            util::ResetMag();
            neigh::Heun(thermal_fluct);
            util::SortSublat();
            util::MagLength();
            // Rotation();
            util::OutputMagToTerm(i);
            util::OutputMagToFile(i);
            util::SumMag(i);

            // SPINWAVES ===================================================================================== //
            // flip a random spin for spinwaves
            if (i >= params::start) {
                neigh::Sx1d(0) = 0;
                neigh::Sy1d(0) = 0;
                neigh::Sz1d(0) = 1;
            }

            if ((i >= params::start) && (i % c == 0)){
                spinwaves::file_spnwvs << spinwaves::icount * params::dt_spinwaves << "\t";
                spinwaves::FFTspace();      
            }
            // ================================================================================================ //
        }
        // ==================================================================================================== //
    
        // Carry out time FFT once simulation is complete
        spinwaves::FFTtime();

        // output sum of magnetisation
        util::OutputSumMag();
        
        end = clock();
        std::cout << std::setprecision(10) << "Simulation Time = " << (double)(end - begin) / CLOCKS_PER_SEC << std::endl; 

        // CLOSE FILE
        util::CloseMagFile();
    }

    return 0;
}
