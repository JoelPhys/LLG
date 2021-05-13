#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include "../inc/params1.h"
#include "../inc/cumalloc.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cuintegrate.h"
#include <iostream>
namespace cuint {

	int threadsperblock;
	int bpg;

	void init_device_vars(){
		threadsperblock = 512;
        	bpg = (params::Nspins + threadsperblock - 1) / threadsperblock;
	}

	void cuRotation(){
		cuheun::cuRotfun<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d); 
	}

	void integration(){
		cuheun::cuHeun1<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dtfa, cuthermal::gvalsx, cuthermal::gvalsy, cuthermal::gvalsz, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx,  cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
		cuheun::cuHeun2<<<bpg,threadsperblock>>>(params::Nspins, cuglob::dx_adj, cuglob::dadjncy, cuheun::Htx, cuheun::Hty, cuheun::Htz, cuglob::dSx1d, cuglob::dSy1d, cuglob::dSz1d, cuglob::dJx, cuglob::dJy, cuglob::dJz, cuglob::Hapx, cuglob::Hapy, cuglob::Hapz, cuheun::DelSx, cuheun::DelSy, cuheun::DelSz, cuheun::Sdashnx, cuheun::Sdashny, cuheun::Sdashnz);
	}
}
