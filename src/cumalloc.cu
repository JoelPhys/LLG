#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include "../inc/array.h"
#include "../inc/neighbourlist.h"
#include "../inc/config.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cudefine.h"
#include "../inc/fields.h"
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/defines.h"
#include <cstring> 


namespace cuglob {

	double *dSx1d, *dSy1d, *dSz1d;
	double *Hapx, *Hapy, *Hapz;
	double *dJx, *dJy, *dJz;
	int *dlw, *drw;
	int *dx_adj, *dadjncy;
	Array<double> pJx, pJy, pJz;
	Array<int> px_adj, padjncy;

	//testing for hedgehog
	double *dsurfx, *dsurfy, *dsurfz;

	int device = 0;

	void device_info(){
		cudaGetDevice(&device);
		struct cudaDeviceProp properties;
		cudaGetDeviceProperties(&properties, device);
		TITLE("CUDA DEVICE PROPERTIES");
		std::cout.width(75); std::cout << std::left << "Device name:"; std::cout << properties.name << std::endl;
		std::cout.width(75); std::cout << std::left << "Memory Clock Rate (KHz):"; std::cout << properties.memoryClockRate << std::endl;
    	std::cout.width(75); std::cout << std::left << "Memory Bus Width (bits):"; std::cout << properties.memoryBusWidth << std::endl;
    	std::cout.width(75); std::cout << std::left << "Peak Memory Bandwidth (GB/s):"; std::cout << 2.0*properties.memoryClockRate*(properties.memoryBusWidth/8)/1.0e6 << std::endl;
		std::cout.width(75); std::cout << std::left << "multiprocessors:"; std::cout << properties.multiProcessorCount << std::endl;
		std::cout.width(75); std::cout << std::left << "max threads per processor:"; std::cout << properties.maxThreadsPerMultiProcessor << std::endl;
		std::cout.width(75); std::cout << std::left << "max threads per block:"; std::cout << properties.maxThreadsPerBlock << std::endl;	
	}

	void clear_memory(){
		CUDA_CALL(cudaFree(dSx1d));
		CUDA_CALL(cudaFree(dSy1d));
		CUDA_CALL(cudaFree(dSz1d));
		CUDA_CALL(cudaFree(cuheun::Sdashnx));
		CUDA_CALL(cudaFree(cuheun::Sdashny));
		CUDA_CALL(cudaFree(cuheun::Sdashnz));
		CUDA_CALL(cudaFree(cuheun::DelSx));
		CUDA_CALL(cudaFree(cuheun::DelSy));
		CUDA_CALL(cudaFree(cuheun::DelSz));
		CUDA_CALL(cudaFree(cuheun::Htx));
		CUDA_CALL(cudaFree(cuheun::Hty));
		CUDA_CALL(cudaFree(cuheun::Htz));
		CUDA_CALL(cudaFree(Hapx));
		CUDA_CALL(cudaFree(Hapy));
		CUDA_CALL(cudaFree(Hapz));
		CUDA_CALL(cudaFree(dJx));
		CUDA_CALL(cudaFree(dJy));
		CUDA_CALL(cudaFree(dJz));
		CUDA_CALL(cudaFree(dx_adj));
		CUDA_CALL(cudaFree(dadjncy));
		CUDA_CALL(cudaFree(cuthermal::gvalsx));
		CUDA_CALL(cudaFree(cuthermal::gvalsy));
		CUDA_CALL(cudaFree(cuthermal::gvalsz));
		CUDA_CALL(cudaFree(cuthermal::Te));
		CUDA_CALL(cudaFree(cuthermal::P_it));
		CUDA_CALL(cudaFree(cuthermal::Tp));
		CUDA_CALL(cudaFree(dlw));
		CUDA_CALL(cudaFree(drw));
		INFO_OUT("memory deallocated on GPU device: ", "success" << std::endl);
	}


	void allocate_heun_memory(){

		//Device spin variables
		CUDA_CALL(cudaMalloc((void**)&dSx1d, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&dSy1d, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&dSz1d, sizeof(double)*params::Nspins));

		//Device Heun variables
		CUDA_CALL(cudaMalloc((void**)&cuheun::Sdashnx, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Sdashnx, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::Sdashny, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Sdashny, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::Sdashnz, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Sdashnz, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::DelSx, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::DelSx, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::DelSy, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::DelSy, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::DelSz, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::DelSz, 0.0, sizeof(double) * params::Nspins));

		// Stochastic Magnetic Field
		CUDA_CALL(cudaMalloc((void**)&cuheun::Htx, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Htx, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::Hty, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Hty, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuheun::Htz, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(cuheun::Htz, 0.0, sizeof(double) * params::Nspins));

		//external field
		CUDA_CALL(cudaMalloc((void**)&Hapx, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(Hapx, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&Hapy, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(Hapy, 0.0, sizeof(double) * params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&Hapz, sizeof(double)*params::Nspins));
		CUDA_CALL(cudaMemset(Hapz, 0.0, sizeof(double) * params::Nspins));

		// Random number arrays
		CUDA_CALL(cudaMalloc((void**)&cuthermal::gvalsx, sizeof(float)*params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::gvalsy, sizeof(float)*params::Nspins));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::gvalsz, sizeof(float)*params::Nspins));
		CUDA_CALL(cudaMemset(cuthermal::gvalsx, 0.0, sizeof(float) * params::Nspins));
		CUDA_CALL(cudaMemset(cuthermal::gvalsy, 0.0, sizeof(float) * params::Nspins));
		CUDA_CALL(cudaMemset(cuthermal::gvalsz, 0.0, sizeof(float) * params::Nspins));

		// Jij matrices
		CUDA_CALL(cudaMalloc((void**)&dJx, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dJy, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dJz, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dx_adj, sizeof(int)*neigh::x_adj.size()));
		CUDA_CALL(cudaMalloc((void**)&dadjncy, sizeof(int)*neigh::adjncy.size()));

		// Domain Wall arrays
		CUDA_CALL(cudaMalloc((void**)&dlw, sizeof(int)*geom::lw.size()));
		CUDA_CALL(cudaMalloc((void**)&drw, sizeof(int)*geom::rw.size()));
		CUDA_CALL(cudaMemset(dlw, 0.0, sizeof(int) * geom::lw.size()));
		CUDA_CALL(cudaMemset(drw, 0.0, sizeof(int) * geom::rw.size()));

		//testing for hedgehog
		CUDA_CALL(cudaMalloc((void**)&dsurfx, sizeof(double)*geom::surfx.size()));
		CUDA_CALL(cudaMalloc((void**)&dsurfy, sizeof(double)*geom::surfy.size()));
		CUDA_CALL(cudaMalloc((void**)&dsurfz, sizeof(double)*geom::surfz.size()));
		CUDA_CALL(cudaMemset(dsurfx, 0.0, sizeof(double) * geom::surfx.size()));
		CUDA_CALL(cudaMemset(dsurfy, 0.0, sizeof(double) * geom::surfy.size()));
		CUDA_CALL(cudaMemset(dsurfz, 0.0, sizeof(double) * geom::surfz.size()));

		// Temperature arrays
		CUDA_CALL(cudaMalloc((void**)&cuthermal::Te, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::Tp, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::P_it, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::dzlayer, sizeof(int)*geom::zlayer.size()));

		CUDA_CALL(cudaMemset(cuthermal::Te, 0.0, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMemset(cuthermal::Tp, 0.0, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMemset(cuthermal::P_it, 0.0, sizeof(double)*params::Lz));
		CUDA_CALL(cudaMalloc((void**)&cuthermal::dtfa, sizeof(double)* params::Nspins));
		CUDA_CALL(cudaMemset(cuthermal::dtfa, 0.0, sizeof(double) * params::Nspins)); 


		// clear memory at exit of program
		atexit(clear_memory);

	}

	void copy_temp_to_device(double equilibium_temp){
		Array<double> rTe; 
		Array<double> rTp;

		rTe.resize(params::Lz);
		rTp.resize(params::Lz);
		for (int i = 0; i < params::Lz; i++){
			rTe[i] = equilibium_temp;
			rTp[i] = equilibium_temp;
		}

		CUDA_CALL(cudaMemcpy(cuthermal::Te, rTe.ptr(), sizeof(double) * params::Lz, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(cuthermal::Tp, rTp.ptr(), sizeof(double) * params::Lz, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(geom::zlayer.ptr(), cuthermal::dzlayer, sizeof(int) * params::Nspins, cudaMemcpyDeviceToHost));

	}

	void copy_field_to_device(){
		CUDA_CALL(cudaMemcpy(Hapx, fields::H_appx.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(Hapy, fields::H_appy.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(Hapz, fields::H_appz.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
	}	

	void copy_jij_to_device(){

		//resize 1d arrays
		pJx.resize(neigh::Jijx_prime.size());
		pJy.resize(neigh::Jijy_prime.size());
		pJz.resize(neigh::Jijz_prime.size());	

		px_adj.resize(neigh::x_adj.size());
		padjncy.resize(neigh::adjncy.size());
		
		for (int a = 0; a < neigh::Jijx_prime.size(); a++){
			pJx(a) = neigh::Jijx_prime[a];
			pJy(a) = neigh::Jijy_prime[a];
			pJz(a) = neigh::Jijz_prime[a];
		}

		for (int a = 0; a < neigh::x_adj.size(); a++){
			px_adj(a) = neigh::x_adj[a];
		}

		for (int a = 0; a < neigh::adjncy.size(); a++){
			padjncy(a) = neigh::adjncy[a];
		}		

		CUDA_CALL(cudaMemcpy(dJx, pJx.ptr(), sizeof(double) * neigh::Jijx_prime.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dJy, pJy.ptr(), sizeof(double) * neigh::Jijy_prime.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dJz, pJz.ptr(), sizeof(double) * neigh::Jijz_prime.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dx_adj, px_adj.ptr(), sizeof(int) * (neigh::x_adj.size()), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dadjncy, padjncy.ptr(), sizeof(int) * (neigh::adjncy.size()), cudaMemcpyHostToDevice));

	}

	void copy_spins_to_device(){
		CUDA_CALL(cudaMemcpy(dSx1d, spins::sx1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dSy1d, spins::sy1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dSz1d, spins::sz1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
	}

	void copy_dw_to_device(){
		CUDA_CALL(cudaMemcpy(dlw, geom::lw.ptr(), sizeof(int) * geom::lw.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(drw, geom::rw.ptr(), sizeof(int) * geom::rw.size(), cudaMemcpyHostToDevice));

		//testing for hedgehog
		CUDA_CALL(cudaMemcpy(dsurfx, geom::surfx.ptr(), sizeof(double) * geom::surfx.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dsurfy, geom::surfy.ptr(), sizeof(double) * geom::surfy.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dsurfz, geom::surfz.ptr(), sizeof(double) * geom::surfz.size(), cudaMemcpyHostToDevice));
	}

	void copy_spins_to_host(){
		CUDA_CALL(cudaMemcpy(spins::sx1d.ptr(), dSx1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(spins::sy1d.ptr(), dSy1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(spins::sz1d.ptr(), dSz1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
	}



}
