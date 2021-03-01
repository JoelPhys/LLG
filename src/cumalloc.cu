#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include "../inc/array.h"
#include "../inc/NeighbourList.h"
#include "../inc/params1.h"
#include "../inc/cuheun.h"
#include "../inc/cuthermal.h"
#include "../inc/cudefine.h"

namespace cuglob {

	double *dSx1d, *dSy1d, *dSz1d;
	double *Hapx, *Hapy, *Hapz;
	double *dJx, *dJy, *dJz;
	int *dx_adj, *dadjncy;
	double *dtfa;
	Array<double> pJx, pJy, pJz;
	Array<int> px_adj, padjncy;

	int device = 0;

	void device_info(){
		cudaGetDevice(&device);
		struct cudaDeviceProp properties;
		cudaGetDeviceProperties(&properties, device);
		std::cout << "using " << properties.multiProcessorCount << " multiprocessors" << std::endl;
		std::cout << "max threads per processor " << properties.maxThreadsPerMultiProcessor << std::endl;
		std::cout << "max threads per block " << properties.maxThreadsPerBlock << std::endl;	
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
		std::cout << "memory deallocated on device" << std::endl;
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


		//thermal array
		CUDA_CALL(cudaMalloc((void**)&dtfa, sizeof(double)* params::Nspins));
		CUDA_CALL(cudaMemset(dtfa, 0.0, sizeof(double) * params::Nspins));


		// Jij matrices
		CUDA_CALL(cudaMalloc((void**)&dJx, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dJy, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dJz, sizeof(double)*neigh::Jijx_prime.size()));
		CUDA_CALL(cudaMalloc((void**)&dx_adj, sizeof(int)*neigh::x_adj.size()));
		CUDA_CALL(cudaMalloc((void**)&dadjncy, sizeof(int)*neigh::adjncy.size()));

		// clear memory at exit of program
		atexit(clear_memory);

	}

	void copy_field_to_device(){
		CUDA_CALL(cudaMemcpy(Hapx, params::H_appx.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(Hapy, params::H_appy.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(Hapz, params::H_appz.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
	}

	void copy_thermal_to_device(double Thermal_Fluct){
		Array<double> tfa;
		tfa.resize(params::Nspins);
		for (int a = 0; a < params::Nspins; a++){
			tfa(a) = Thermal_Fluct;
		}

		CUDA_CALL(cudaMemcpy(dtfa, tfa.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));    
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
		CUDA_CALL(cudaMemcpy(dSx1d, neigh::Sx1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dSy1d, neigh::Sy1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(dSz1d, neigh::Sz1d.ptr(), sizeof(double) * params::Nspins, cudaMemcpyHostToDevice));
	}

	void copy_spins_to_host(){
		CUDA_CALL(cudaMemcpy(neigh::Sx1d.ptr(), dSx1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(neigh::Sy1d.ptr(), dSy1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		CUDA_CALL(cudaMemcpy(neigh::Sz1d.ptr(), dSz1d, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
	}



}
