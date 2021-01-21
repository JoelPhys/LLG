#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <time.h>
#include <fftw3.h>
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/mathfuncs.h"
#include "../inc/params1.h"
#include "../inc/geom.h"
#include "../inc/NeighbourList.h"
// #include "inc/FFT.h"
#include "../inc/spinwaves.h"

#define IMAG 1
#define REAL 0

int main(int argc, char* argv[]){

    // functions ============================================================================================== //
    // params::intitialiseConfig("/Users/Hirst/Documents/PhD/LLG_code/Mn2Au/config_files/af_gga.cfg");
    params::intitialiseConfig(argv[1]); 
    params::readparams();
    // spinwaves::initialiseFFT();
    geom::CreateLattice();
    neigh::ReadFile();
    neigh::InteractionMatrixJerome();
    neigh::IntialisePointersNL();
    // ======================================================================================================== //

    // ======= Temperature ==================================================================================== //
    const double Temp = (atof(argv[2]));  // REMOVE THIS ASAP
    const double thermal_fluct = params::thermal_const * sqrt(Temp);
    std::cout << "Temperature = " << Temp << std::endl;
    // ========================================================================================================= //
    

    // ======= Initiliase Spin Position ======================================================================== //
    Array4D<double> Sx4, Sy4, Sz4;
    Sx4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
    Sy4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
    Sz4.resize(params::Lx, params::Ly, params::Lz, params::Nq);

    Sz4.IFill(0);
    Sx4.IFill(0);
    Sy4.IFill(0);

    int count1d = 0;

    for (int x = 0; x < params::Lx; x++){
        for (int y = 0; y < params::Ly; y++){
            for (int z = 0; z < params::Lz; z++){

                Sz4(x,y,z,0) = 1;
                // Sz4(x,y,z,1) = -1;
                // Sz4(x,y,z,2) = 1;
                // Sz4(x,y,z,3) = -1;

                neigh::Sz1d(count1d + 0) = Sz4(x,y,z,0);
                // neigh::Sz1d(count1d + 1) = Sz4(x,y,z,1);
                // neigh::Sz1d(count1d + 2) = Sz4(x,y,z,2);
                // neigh::Sz1d(count1d + 3) = Sz4(x,y,z,3);

                count1d += params::Nq; 
            }
        }
    }
        std::cout << __LINE__ << std::endl;


    double M[params::Nsublat][3];
    double Mmag[params::Nsublat];
    double MdivMs[params::Nsublat];
    int c;
    c = params::dt_spinwaves / params::dt;
    // ================================================================================================== //
        
    // ======= OPEN OUTPUT FILE ======================================================================== //
    std::stringstream sstr;
    sstr << "output/magnetisation/mag_tsteps_" << params::Nt << "_T_" << Temp << ".txt";
    std::ofstream myfile;
    myfile.open(sstr.str());
    // ================================================================================================= //

    // ========== Time + Temp Variables ================================================================= //
    double t = 0;
    double tau = 0;

    clock_t begin, end;
    double time_spent;
    begin = clock();
    // ================================================================================================== //

    std::cout << c << std::endl;

    // ========== LOOP THROUGH TIMESTEPS ================================================================ //
    for (int i = 0; i < params::Nt; i++){
        
        t = t + params::dt;
        tau = tau + params::dtau;

        // Reset magnetisation after each timestep ======================================================= //
        for (int ii = 0; ii < params::Nsublat; ii++){ 
            M[ii][0] = 0;
            M[ii][1] = 0;
            M[ii][2] = 0;
        }
        // =============================================================================================== //

        // INTERGRATE ==================================================================================== //
        neigh::Heun(thermal_fluct);
        // =============================================================================================== //

        // SPINWAVES ===================================================================================== //
        // flip a random spin for spinwaves
        // if (i >= params::start) {
        //     neigh::Sx1d(0) = 1;
        //     neigh::Sy1d(0) = 0;
        //     neigh::Sz1d(0) = 0;
        // }

        // if ((i >= params::start) && (i % c == 0)){
        //     spinwaves::file_spnwvs << spinwaves::icount * params::dt_spinwaves << "\t";
        //     spinwaves::FFTspace();      
        // }
        // ==========================================================================================

        for (int a = 0; a < params::Nspins; a++){     
            if (modfunc(params::Nq,a) == 0){
                M[0][0] += neigh::Sx1d(a);
                M[0][1] += neigh::Sy1d(a);
                M[0][2] += neigh::Sz1d(a); 
            }
            else if (modfunc(params::Nq,a) == 1) {
                M[1][0] += neigh::Sx1d(a);
                M[1][1] += neigh::Sy1d(a);
                M[1][2] += neigh::Sz1d(a); 
            }
            else if (modfunc(params::Nq,a) == 2) {
                M[0][0] += neigh::Sx1d(a);
                M[0][1] += neigh::Sy1d(a);
                M[0][2] += neigh::Sz1d(a);
            }
            else if (modfunc(params::Nq,a) == 3) {
                M[1][0] += neigh::Sx1d(a);
                M[1][1] += neigh::Sy1d(a);
                M[1][2] += neigh::Sz1d(a); 
            }
            else std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
        }

        for (int l = 0; l < params::Nsublat; l++){
            Mmag[l] = sqrt(M[l][0] * M[l][0] + M[l][1] * M[l][1] + M[l][2] * M[l][2]);
            MdivMs[l] = Mmag[l] / (params::NmomentsSubLat);
        }
        
        // Output to file
        // myfile << i << " ";
        // myfile << M[0] / params::Nmoments << " " << M[1] / params::Nmoments << " " << M[2] / params::Nmoments << " " << MdivMs << "\n"; 

        std::cout << i << " ";

        for (int l = 0; l < params::Nsublat; l++){
            for (int m = 0; m < 3; m++){
                std::cout << M[l][m] / params::NmomentsSubLat << "\t"; 
            }
            std::cout << MdivMs[l] << "\t";
        }
        std::cout << "\n";

    }
    // ==================================================================================================== //
    
    // Carry out time FFT once simulation is complete
        // spinwaves::FFTtime();

    end = clock();
    std::cout << std::setprecision(10) << "Time in seconds = " << (double)(end - begin) / CLOCKS_PER_SEC << std::endl; 

    // CLOSE FILE
    myfile << std::flush;
    myfile.close();

    return 0;
}
