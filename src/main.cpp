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
#include "../inc/defects.h"
#include "../inc/defines.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/spinwaves.h"
#include "../inc/mathfuncs.h"
#include "../inc/neighbourlist.h"

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
	neigh::ReadFile();
	neigh::InteractionMatrix();
	heun::init();
	util::init();

	if (params::simtype == "spinwaves"){
		spinwaves::init();
	}
	// ======================================================================================================== //

	// ======= Temperature ==================================================================================== //
	const double Temp = (atof(argv[2]));
	const double thermal_fluct = params::thermal_const * sqrt(Temp);
	INFO_OUT("Temperature: ", Temp << "(K)");
	INFO_OUT("Thermal Fluct: ", thermal_fluct);
	// ========================================================================================================= //

	// ======= Initiliase Spin Position ======================================================================== //
	spins::populate();
	defects::populate();
	util::readexternalspins(argv[3]);

	// testing for hedgehog
	if (params::simtype == "DW"){
		geom::initdw();
	}

	#ifdef CUDA
	std::cout << "CUDA Simulation" << std::endl;
	cuglob::device_info();
	cuglob::allocate_heun_memory();
	cuheun::allocate_heun_consts();
	cuthermal::init_cuthermal(Temp);
	cuglob::copy_temp_to_device(Temp);
	cuglob::copy_spins_to_device();
	cuglob::copy_field_to_device();
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
	c = params::dt_spinwaves / params::dt;

	util::InitMagFile(Temp);

	if (params::simtype == "DW"){		
	util::InitDWFile(Temp);
	}

	TITLE("SIMULATION STARTING");
	util::startclock();
		
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

			if (params::OutputToTerminal == true){
				util::OutputMagToTerm(i);
			}
			util::OutputMagToFile(i);
			
			if (params::simtype == "DW"){		
				util::OutputDWtoFile(i);
			}
			if (params::simtype == "spinwaves"){
				if ((i >= params::start)){
					spinwaves::file_spnwvs << spinwaves::icount * params::dt_spinwaves << "\t";
					spinwaves::FFTspace();      
				}
			}
		}

		t = t + params::dt;
		tau = tau + params::dtau;

		#ifdef CUDA
			cufuncs::cuTemperature(params::temptype, static_cast<double>(i) * params::dt, params::ttm_start);
			cufuncs::cuFields(fields::type, static_cast<double>(i)  * params::dt, fields::start_time, fields::end_time, fields::height);
			cuthermal::gen_thermal_noise();
			cufuncs::integration(static_cast<double>(i));
		#else	
			heun::integration(thermal_fluct);
		#endif
		
	}
	// ==================================================================================================== //

	util::OutputLatticetoTerm();

	// Carry out time FFT once simulation is complete
	if (params::simtype == "spinwaves"){
	spinwaves::FFTtime();
	}
						
	TITLE("FINISHED");
	util::endclock();
	util::CloseMagFile();


	return 0;
}
