// cpp header files
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <time.h>
#include <cstring>

// my header files
#include "../inc/util.h"
#include "../inc/geom.h"
#include "../inc/heun.h"
#include "../inc/spins.h"
#include "../inc/fftw3.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/fields.h"
#include "../inc/thermal.h"
#include "../inc/defects.h"
#include "../inc/defines.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/spinwaves.h"
#include "../inc/mathfuncs.h"
#include "../inc/neighbourlist.h"

//openMP Header Files
#ifdef MP
#include <omp.h>
#include "../inc/mpheun.h"
#endif

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

int main(int argc, char* argv[]){

	
	// functions ============================================================================================== //
	params::banner();
	params::intitialiseConfig(argv[1]); 
	params::readparams();
	spins::init();
	fields::readfields();
	geom::CreateLattice();
	geom::CountDistinct();
	geom::CreateIntLattice();
	defects::init();
	util::init();
	#ifdef MP
		mpheun::init();
	#else
		heun::init();
	#endif
	
	if (params::simtype == "Spinwaves"){
		spinwaves::init();
	}
	// ======================================================================================================== //

	// ======= Initiliase Spin Position ======================================================================== //
	spins::populate();
	std::string arg4, arg3 = argv[3];
	if ((argc==4) && (arg3 == "2")){
		std::cout << "ERROR: No external spin file provided. Exiting." << std::endl;
		exit(0);
	}
	else if (argv[4]!=NULL){
		DEBUGGER;
		arg4 = argv[4];
		util::readexternalspins(arg4);
		DEBUGGER;
	}
	defects::populate();
	spins::randomise();
	
	// testing for hedgehog
	if (params::simtype == "DW"){
		geom::initdw();
	}

	// ======= Temperature ==================================================================================== //
	TITLE("TEMPERATURE")
	const double Temp = (atof(argv[2]));

	// check if argument is a number
	char *endptr;
	double d = strtod(argv[2], &endptr);
	int ok = endptr == argv[2] + strlen(argv[2]);
	if (!ok) {
		std::cout << "ERROR: Temperature is not a number. \nExiting." << std::endl;
		exit(0);
	}
	
	// Print temperature
	INFO_OUT("Initial Temperature: ", Temp << " [K]");
	thermal::initthermal(Temp);
	// ========================================================================================================= //
	

	neigh::ReadFile();
	neigh::InteractionMatrix();
	#ifdef MP
		std::cout << "MMMMMMMMMMMMPPPPPPPPPPPPPPP" << std::endl;
		mpheun::init();
	#else
		heun::init();
	#endif

	#ifdef CUDA
	std::cout << "CUDA Simulation" << std::endl;
	cuglob::device_info();
	cuglob::allocate_heun_memory();
	cufields::allocate_field_variables();
	cuheun::allocate_heun_consts();
	cuthermal::init_cuthermal(Temp);
	cuglob::copy_temp_to_device(Temp);
	cuglob::copy_spins_to_device();
	cuglob::copy_field_to_device();
	cuglob::copy_damp_to_device();
	cuglob::copy_energy_to_device();
	if (params::simtype == "DW"){
		cuglob::copy_dw_to_device();
	}
	cuglob::copy_jij_to_device();
	cufuncs::init_device_vars();
	cuthermal::curand_generator();
	#endif    


	// Initialise some variables ======================================================================== // 
	double t = 0;
	double tau = 0;

	int c;
	c = spinwaves::dt_spinwaves / params::dt;

	util::InitMagFile(Temp);
	util::InitFldFile(Temp);

	if (params::simtype == "DW"){		
	util::InitDWFile(Temp);
	}
	TITLE("SIMULATION STARTING");
	util::startclock();

	// ========== LOOP THROUGH TIMESTEPS ================================================================ //
	for (int i = 0; i < params::Nt; i++){

        if (i == 0) {
            std::cout << "Rotation matrix applied with angle " << params::angle << " (rad) at time t = " << std::scientific << i * params::dt << " (s)" << std::endl;
        #ifdef CUDA
            cufuncs::cuRotation();
        #elif MP
			mpheun::rotation();
		#else
            heun::rotation();
        #endif
		} 


		//if (i == 5000){
		//	for (int l = 0; l < params::Nq; l++){
		//		params::lambda[l] = 0.0;
		//		params::thermal_const[l] = 0.0;
		//		params::lambdaPrime[l] = 1 / (1+(params::lambda[l]*params::lambda[l]));
		//		cuglob::copy_damp_to_device();
		//		cuglob::copy_temp_to_device(Temp);
		//	}
		//}
	
		// Output Lattice subroutine	
		//if ((i % params::OutputLatticeStep == 0) && (i >= params::OutputLatticeStart) ){
		//	if (params::OutputLattice == true){
		//		util::OutputLatticetoFile(Temp);
		//	}
		//}

		// Magnetisation subroutine
		if (i % params::outputstep == 0){
			#ifdef CUDA
			cuglob::copy_spins_to_host();
			cuglob::copy_field_to_host();
			cuglob::copy_energy_to_host();
			#endif	
			util::ResetMag();
			util::SortSublat();
			util::MagLength();
			if (params::OutputToTerminal == true){
				util::OutputMagToTerm(i);
			}
			//std::cout << i*params::dt << " " << params::mu_s[0]*0.5*7.246e22*heun::spin_temp<< "\n";
			util::OutputMagToFile(i);
			util::OutputFldToFile(i);
			if (params::simtype == "DW"){		
				util::OutputDWtoFile(i);
			}

			// If using two-temperature model - output temperature profile to file
			if (thermal::temptype == "ttm"){
				thermal::ttmtofile(static_cast<double>(i) * params::dt);
			}
		}

		// Output Lattice subroutine	
		if ((i % params::OutputLatticeStep == 0) && (i >= params::OutputLatticeStart) ){
			if (params::OutputLattice == true){
				util::OutputLatticetoFile(Temp);
			}
		}

		// Spinwaves subroutine
		if (params::simtype == "Spinwaves"){	
			if ((i >= spinwaves::start)){
				if (i % spinwaves::int_dt_spinwaves == 0){
					if (i % params::outputstep == 0){
						spinwaves::file_spnwvs << spinwaves::icount * spinwaves::dt_spinwaves << "\t";
						spinwaves::FFTspace();      
					}
					else {
						#ifdef CUDA
						cuglob::copy_spins_to_host();
						//cuglob::copy_field_to_host();
						//cuglob::copy_energy_to_host();
						#endif
						spinwaves::file_spnwvs << spinwaves::icount * spinwaves::dt_spinwaves << "\t";
						spinwaves::FFTspace();      
					}
				}
			}
		}

		// Increment time
		t = t + params::dt;
		tau = tau + params::dtau;


		// Calculate temperature and fields on cpu. This is done even when integrating on GPU for outputting to file
		thermal::cputemperature(static_cast<double>(i) * params::dt);
		fields::calculate(static_cast<double>(i) * params::dt);	

		#ifdef CUDA
			cufuncs::cuTemperature(thermal::temptype, static_cast<double>(i) * params::dt, thermal::ttm_start);
			cufuncs::cuFields(fields::type, static_cast<double>(i)  * params::dt, fields::start_time, fields::end_time, fields::height);
			cuthermal::gen_thermal_noise();
			cufuncs::integration(static_cast<double>(i));
			// cufields::testing(i);
		#elif MP
			mpheun::integration();
		#else	
			heun::integration();
		#endif


	}
	// ==================================================================================================== //

	if (params::OutputLattice == true){
		util::OutputLatticetoFile(Temp);
	}

	// Carry out time FFT once simulation is complete
	if (params::simtype == "Spinwaves"){
		spinwaves::FFTtime();
	}
						
	TITLE("FINISHED");
	util::endclock();
	util::CloseMagFile();

	if (thermal::temptype == "ttm"){
		thermal::closettmfile();
	}
	return 0;
}
