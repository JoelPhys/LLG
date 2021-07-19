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
#include "../inc/fields.h"
#include "../inc/geom.h"
#include "../inc/NeighbourList.h"
#include "../inc/util.h"
// #include "inc/FFT.h"
#include "../inc/spinwaves.h"

//Cuda Header files
#ifdef CUDA
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include "../inc/cuthermal.h"
#include "../inc/cuheun.h"
#include "../inc/cufields.h"
#include "../inc/cumalloc.h"
#include "../inc/cudefine.h"
#include "../inc/cufuncs.h"
#endif

#define IMAG 1
#define REAL 0



int main(int argc, char* argv[]){

	// functions ============================================================================================== //
	params::intitialiseConfig(argv[1]); 
	params::readparams();
	fields::readfields();
	geom::CreateLattice();
	geom::CountDistinct();
	geom::CreateIntLattice();
	// geom::InitDomainWall();
	neigh::ReadFile();
	neigh::InteractionMatrix();
	neigh::IntialisePointersNL();
	util::InitUtil();
	IdentityMatrix();
	// spinwaves::initialiseFFT();
	#define TERM params::term;
	// ======================================================================================================== //

	// ======= Temperature ==================================================================================== //
	const double Temp = (atof(argv[2]));
	const double thermal_fluct = params::thermal_const * sqrt(Temp);
	std::cout << "Temperature = " << Temp << "(K)" << std::endl;
	std::cout << "Thermal Fluct = " << thermal_fluct << std::endl;
	// ========================================================================================================= //

	// ======== Set clock ====================================================================================== //
	clock_t begin, end;
	double time_spent;
	begin = clock();
	// ========================================================================================================= //

	if ((std::string(argv[3]) == "1") || (std::string(argv[3]) == "3")){

		// ======= Initiliase Spin Position ======================================================================== //
		geom::InitSpins();


		#ifdef CUDA
		std::cout << "CUDA Simulation" << std::endl;
		cuglob::device_info();
		cuglob::allocate_heun_memory();
		cuheun::allocate_heun_consts();
		cuthermal::init_cuthermal(Temp);
		cuglob::copy_temp_to_device(Temp);
		cuglob::copy_spins_to_device();
		cuglob::copy_field_to_device();
		// cuglob::copy_dw_to_device();
		cuglob::copy_jij_to_device();
		cufuncs::init_device_vars();
		cuthermal::curand_generator();
		#endif            
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

		util::InitMagFile(Temp);
		// util::InitDWFile(Temp);
		
		// ========== LOOP THROUGH TIMESTEPS ================================================================ //
		for (int i = 0; i < params::Nt; i++){

			#ifdef CUDA
			// if (i ==  params::Nt / 2) {
			// 	std::cout << "Rotation matrix applied with angle " << params::angle << " (rad) at time t = " << std::scientific << i * params::dt << " (s)" << std::endl;
			// 	cufuncs::cuRotation();
			// }
			#endif

			if (i % params::outputstep == 0){
				#ifdef CUDA
				cuglob::copy_spins_to_host();
				#endif	
				util::ResetMag();
				util::SortSublat();
				util::MagLength();
				cuthermal::testing(static_cast<double>(i));

				if (params::OutputToTerminal == true){
					util::OutputMagToTerm(i);
				}
				util::OutputMagToFile(i);
				
				
				// util::OutputDWtoFile(i);
				// if ((i >= params::start)){
			    // 		spinwaves::file_spnwvs << spinwaves::icount * params::dt_spinwaves << "\t";
			    // 		spinwaves::FFTspace();      
				// }
			}

			t = t + params::dt;
			tau = tau + params::dtau;

			#ifdef CUDA
			cufuncs::cuTemperature(params::temptype, static_cast<double>(i) * params::dt, params::ttm_start);
			// cufuncs::cuMultiPulse(static_cast<double>(i));
			cuthermal::gen_thermal_noise();
			cufuncs::integration(static_cast<double>(i));
			// cufuncs::cuDomainWall();
			#else
			neigh::Heun(thermal_fluct);
			#endif
			
		}
		// ==================================================================================================== //

		// Carry out time FFT once simulation is complete
		// spinwaves::FFTtime();

		// output sum of magnetisation
		// util::OutputSumMag();

		end = clock();
		std::cout << std::setprecision(10) << "Simulation Time = " << (double)(end - begin) / CLOCKS_PER_SEC << std::endl; 

		// CLOSE FILE
		util::CloseMagFile();

		//Deallocate Device memory
		#ifdef CUDA 
		cuthermal::destroy_generator();
		#endif
	}

	return 0;
}
