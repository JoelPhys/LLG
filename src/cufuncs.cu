#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include "../inc/params1.h"
#include "../inc/cumalloc.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cufuncs.h"
#include "../inc/cufields.h"
#include "../inc/fields.h"
#include <iostream>

namespace cufuncs {

	int threadsperblock;
	int bpg;

	void init_device_vars(){
		threadsperblock = 256;
        bpg = (params::Nspins + threadsperblock - 1) / threadsperblock;
	}

	void cuSquarePulse(double time){
		cufields::square_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::start_time, fields::end_time, fields::height, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
	}

	void cuGaussPulse(double time){
		cufields::gaussian_pulse<<<bpg,threadsperblock>>>(params::Nspins, time, fields::height, fields::std_dev, fields::centre_pos, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz);
	}


	void cuRotation(){
		cuheun::cuRotfun<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d); 
	}

	void integration(double time){
		cuheun::cuHeun1<<<bpg,threadsperblock>>>(params::Nspins, time, cuglob::dtfa, cuthermal::gvalsx, cuthermal::gvalsy, cuthermal::gvalsz, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx,  cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
		cuheun::cuHeun2<<<bpg,threadsperblock>>>(params::Nspins, time, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx, cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
	}
}
