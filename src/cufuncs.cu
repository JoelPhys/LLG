#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <sstream>
#include "../inc/params1.h"
#include "../inc/cumalloc.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cufuncs.h"
#include "../inc/cufields.h"
#include "../inc/fields.h"
#include "../inc/geom.h"
#include <iostream>

namespace cufuncs {

	int threadsperblock;
	int bpg;
	int nspinsdw;

	void init_device_vars(){
		threadsperblock = 256;
        	bpg = (params::Nspins + threadsperblock - 1) / threadsperblock;
		nspinsdw = params::Ly*params::Lz*params::Nq;
	}

	void cuDomainWall(){
		cuheun::cuFixSpins1<<<bpg,threadsperblock>>>(nspinsdw, cuglob::dlw, cuglob::drw, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d);
	}

	void cuSquarePulse(double time){
		cufields::square_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::start_time, fields::end_time, fields::height, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
	}

	void cuGaussPulse(double time){
		cufields::gaussian_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
	}

	void cuMultiPulse(double time){
		cufields::multi_cycle_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, fields::freq, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
	}

	void cuTemperature(std::string type, double time, double ttm_start){
		if (type == "ttm"){
			if (time < ttm_start){

				cuthermal::ttf<<<bpg,threadsperblock>>>(time, params::Nspins, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
			}
			else {
				cuthermal::ttm<<<bpg,threadsperblock>>>(time - ttm_start, params::Lz, cuthermal::Te, cuthermal::Tp, cuthermal::P_it);
				cuthermal::ttf<<<bpg,threadsperblock>>>(time - ttm_start, params::Nspins, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
			}
		}
		else if (type == "constant"){
			cuthermal::ttf<<<bpg,threadsperblock>>>(time - ttm_start, params::Nspins, cuthermal::dtfa, cuthermal::Te, cuthermal::dzlayer);
		}
	}
	
	void cuRotation(){
		cuheun::cuRotfun<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d); 
	}

	void integration(double time){
		cuheun::cuHeun1<<<bpg,threadsperblock>>>(params::Nspins, time, cuthermal::dtfa, cuthermal::gvalsx, cuthermal::gvalsy, cuthermal::gvalsz, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx,  cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
		// cuheun::cuFixSpins2<<<bpg,threadsperblock>>>(nspinsdw, cuglob::dlw, cuglob::drw, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
		cuheun::cuHeun2<<<bpg,threadsperblock>>>(params::Nspins, time, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx, cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
	}
}
