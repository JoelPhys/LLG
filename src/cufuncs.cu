// cpp header files
#include <cuda.h>
#include <sstream>
#include <iomanip>
#include <curand.h>
#include <iostream>
#include <cuda_runtime.h>

// my header files
#include "../inc/geom.h"
#include "../inc/neighbourlist.h"
#include "../inc/fields.h"
#include "../inc/config.h"
#include "../inc/cuheun.h"
#include "../inc/cufuncs.h"
#include "../inc/cufields.h"
#include "../inc/cumalloc.h"
#include "../inc/cuthermal.h"



namespace cufuncs {

	int threadsperblock;
	int bpg;
	int nspinsdw;

	void init_device_vars(){
		threadsperblock = 256;
        bpg = (params::Nspins + threadsperblock - 1) / threadsperblock;
		nspinsdw = params::Ly*params::Lz*params::Nq;

		//testing for hedgehog
		nspinsdw = (2*(params::Lx*params::Ly) + 2*(params::Lx*(params::Ly-2)) + 2*((params::Lx-2)*(params::Ly-2)))*params::Nq;
	}

	void cuDomainWall(){
		cuheun::cuFixSpins1<<<bpg,threadsperblock>>>(nspinsdw, cuglob::dlw, cuglob::drw, cuglob::dsurfx, cuglob::dsurfy, cuglob::dsurfz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d);
	}

	void cuFields(std::string type, double time, double start_time, double end_time, double height){
		if (type == "Uniform"){
			cufields::uniform<<<bpg,threadsperblock>>>(params::Nspins, fields::cuniform[0], fields::cuniform[1], fields::cuniform[2], cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Uniform_Staggered"){
			cufields::uniform_staggered<<<bpg,threadsperblock>>>(params::Nspins, fields::cuniform[0], fields::cuniform[1], fields::cuniform[2], cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Square_Pulse"){
			cufields::square_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, start_time, end_time, height, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Gaussian_Pulse"){
			cufields::gaussian_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Multi_Cycle_Pulse"){
			cufields::multi_cycle_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, fields::freq, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Square_Pulse_Staggered"){
			cufields::square_pulse_staggered<<<bpg,threadsperblock>>>(params::Nspins, time, start_time, end_time, height, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Gaussian_Pulse_Staggered"){
			cufields::gaussian_pulse_staggered<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Multi_Cycle_Pulse_Staggered"){
			cufields::multi_cycle_pulse_staggered<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, fields::freq, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else if (type == "Sine_Pulse"){
			cufields::sine_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::freq, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
		}
		else {
			std::cout << "ERROR: Unkown field type: " << type << std::endl;
			exit(0);		
		}
		
	}
	

	void cuTemperature(std::string type, double time, double ttm_start){
		if (type == "ttm"){
			if (time < ttm_start){

				cuthermal::ttf<<<bpg,threadsperblock>>>(time, params::Nspins, cuthermal::dconst, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
			}
			else {
				// std::cout << std::scientific << time << " " << ttm_start << std::endl;
				cuthermal::ttm<<<bpg,threadsperblock>>>(time - ttm_start, params::Lz, cuthermal::Te, cuthermal::Tp, cuthermal::P_it);
				cuthermal::ttf<<<bpg,threadsperblock>>>(time - ttm_start, params::Nspins, cuthermal::dconst, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
			}
		}
		else if (type == "constant"){
			cuthermal::ttf<<<bpg,threadsperblock>>>(time - ttm_start, params::Nspins, cuthermal::dconst, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
		}
		else if (type == "uniform_gradient"){
			cuthermal::ttfg<<<bpg,threadsperblock>>>(time - ttm_start, params::Nspins, cuthermal::dconst, cuthermal::dtfa, cuthermal::Te, cuthermal::dxlayer, params::temp_gradient);
		}
	}
	
	void cuRotation(){
		cuheun::cuRotfun<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d); 
	}

	void integration(double time){
		cuheun::cuHeun1<<<bpg,threadsperblock>>>(cuglob::djind, neigh::nsimspin, time, cuglob::dsimspin, cuglob::c_lambda, cuglob::c_lambdap, cuthermal::dtfa, cuthermal::gvalsx, cuthermal::gvalsy, cuthermal::gvalsz, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx_new, cuglob::dJy_new, cuglob::dJz_new, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx,  cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
		cuheun::cuHeun2<<<bpg,threadsperblock>>>(cuglob::djind, neigh::nsimspin, time, cuglob::dsimspin, cuglob::c_lambda, cuglob::c_lambdap, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx_new, cuglob::dJy_new, cuglob::dJz_new, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx, cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
	}

}
