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
// #include "inc/FFT.h"
#include "../inc/spinwaves.h"

#define IMAG 1
#define REAL 0

int main(int argc, char* argv[]){

    // functions ============================================================================================== //
    params::intitialiseConfig(argv[1]); 
    params::readparams();
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

    double M[params::Nsublat][3];
    double Mmag[params::Nsublat];
    double MdivMs[params::Nsublat];
    double MdivMsSum[params::Nsublat];
    double Msum[params::Nsublat][3];
    double MsumSQR[params::Nsublat][3];
    int isum = 0;
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
        std::cout << i << " ";

        for (int l = 0; l < params::Nsublat; l++){
            for (int m = 0; m < 3; m++){
                std::cout << M[l][m] / params::NmomentsSubLat << "\t"; 
            }
            std::cout << MdivMs[l] << "\t";
        }
        std::cout << "\n";

        // Sum magnetisation and average over temperature
        if (i > (params::Nt / 10) - 1) {
            for (int l = 0; l < params::Nsublat; l++){
                for (int m = 0; m < 3; m++){
                    Msum[l][m] += M[l][m] / params::NmomentsSubLat;
                    MsumSQR[l][m] += (M[l][m] / params::NmomentsSubLat) * (M[l][m] / params::NmomentsSubLat);
                }
                MdivMsSum[l] += MdivMs[l];
            }
            isum++;
        }

    }
    // ==================================================================================================== //

    std::cout << "For averaging: " << std::endl;

    for (int l = 0; l < params::Nsublat; l++){
            for (int m = 0; m < 3; m++){
                    std::cout << Msum[l][m] <<  "\t";
            }
    }

    std::cout << MdivMsSum[0] << "\t" << MdivMsSum[1];
    std::cout << "\n";

    for (int l = 0; l < params::Nsublat; l++){
            for (int m = 0; m < 3; m++){
                    std::cout << MsumSQR[l][m] <<  "\t";
            }
    }
    std::cout << "\n";
    std::cout << isum << std::endl;
 
    // Carry out time FFT once simulation is complete
    // spinwaves::FFTtime();

    end = clock();
    std::cout << std::setprecision(10) << "Time in seconds = " << (double)(end - begin) / CLOCKS_PER_SEC << std::endl; 

    // CLOSE FILE
    myfile << std::flush;
    myfile.close();

    return 0;
}
