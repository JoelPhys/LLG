// cpp header files
#include <ctime>
#include <cuda.h>
#include <curand.h>
#include <iostream>
#include <curand_kernel.h>

// my header files
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/thermal.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/cudefine.h"

namespace cuthermal {


    // two temperature model constants
    __constant__ double c_gamma_e;
	__constant__ double c_Cp;
	__constant__ double c_kappa_0;
	__constant__ double c_delta;
	__constant__ double c_Gep;
	__constant__ double c_P_0;
	__constant__ double c_t0;
	__constant__ double c_tau;
	__constant__ double c_Nz;
	__constant__ double c_dz;
    __constant__ double c_dt;
    __constant__ double c_oneOvrdzdz;
    __constant__ double c_oneOvr2dz;
    // __constant__ double c_thermal_const;
    __constant__ int c_Nq;

    // Thermal Constant
	double *dconst;

    // thermal gradient constants
    __constant__ double c_grad;

    double *Te, *Tp, *P_it;
	double *dtfa;
    int *dxlayer, *dylayer, *dzlayer;

    // stochastic noise variables
    float *gvalsx, *gvalsy, *gvalsz;
    curandGenerator_t gen;


    void init_cuthermal(double equilibrium_temp){


		std::cout << thermal::gamma_e << std::endl;           
		std::cout << thermal::Cp << std::endl;                
		std::cout << thermal::kappa_0 << std::endl;           
		std::cout << thermal::delta << std::endl;             
		std::cout << thermal::Gep << std::endl;               
		std::cout << thermal::P_0 << std::endl;               
		std::cout << thermal::t0 << std::endl;                
		std::cout << thermal::tau << std::endl;               
		std::cout << thermal::oneOvrdzdz << std::endl; 
		std::cout << thermal::oneOvr2dz << std::endl; 


        // Constants two temperature model
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_gamma_e), &thermal::gamma_e, sizeof(double)));                        
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Cp), &thermal::Cp, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_kappa_0), &thermal::kappa_0, sizeof(double)));                        
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_delta), &thermal::delta, sizeof(double)));                            
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Gep), &thermal::Gep, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_P_0), &thermal::P_0, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_t0), &thermal::t0, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_tau), &thermal::tau, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Nz), &params::Lz, sizeof(int)));                                     
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_dz), &params::c1, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_dt), &params::dt, sizeof(double)));                          
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Nq), &params::Nq, sizeof(int)));
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_oneOvrdzdz), &thermal::oneOvrdzdz, sizeof(double)));
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_oneOvr2dz), &thermal::oneOvr2dz, sizeof(double)));

        // Constants for thermal gradient
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_grad), &thermal::temp_gradient, sizeof(double)));
    }

    void destroy_generator(){
        curandDestroyGenerator(gen);
        INFO_OUT("generator destroyed: ", "success");
    }

    void curand_generator(){
	std::time_t result = std::time(nullptr);
	int seed = static_cast<int>(result);
    INFO_OUT("time since epoch:",result << " [s]");
	CURAND_CALL(curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_MTGP32));
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, seed));
	INFO_OUT("Curand Seed = ",seed);

    atexit(destroy_generator);
    }

    void gen_thermal_noise(){
        CURAND_CALL(curandGenerateNormal(gen, gvalsx, params::Nspins, 0.0, 1.0));
        CURAND_CALL(curandGenerateNormal(gen, gvalsy, params::Nspins, 0.0, 1.0));
        CURAND_CALL(curandGenerateNormal(gen, gvalsz, params::Nspins, 0.0, 1.0));
    }
    
    __global__ void ttm(double time, int Nz, double *Te, double *Tp, double *P_it)
    {
        const int i = blockDim.x * blockIdx.x + threadIdx.x; 
        
        if (i < Nz){

            double Tep1;
            double Tpp1;
            //double z;

            // if (i == 0)
            // {
                P_it[i]=c_P_0*exp(-((time-c_t0)/c_tau)*((time-c_t0)/c_tau));
                Tep1=Te[0] + (c_dt/(c_gamma_e*Te[0]))*(c_Gep*(Tp[0]-Te[0]) + P_it[0] + c_kappa_0*( (Te[0]/Tp[0]) * 2.0*(Te[1]-Te[0])*c_oneOvrdzdz));
                Tpp1=Tp[0]+(c_dt*c_Gep/c_Cp)*(Te[0]-Tp[0]);
            // }
            // if (i == Nz-1)
            // {
            //     z=static_cast<double>(Nz-1)*c_dz;
            //     P_it[Nz-1]=c_P_0*exp(-((time-c_t0)/c_tau)*((time-c_t0)/c_tau))*exp(-z/c_delta);
            //     Tep1=Te[Nz-1]+(c_dt/(c_gamma_e*Te[Nz-1]))*(c_Gep*(Tp[Nz-1]-Te[Nz-1])+P_it[Nz-1]+c_kappa_0*( (Te[Nz-1]/Tp[Nz-1]) * 2.0*(Te[Nz-2]-Te[Nz-1])*c_oneOvrdzdz));
            //     Tpp1=Tp[Nz-1]+(c_dt*c_Gep/c_Cp)*(Te[Nz-1]-Tp[Nz-1]);
            // }
            // if ((1 <= i) && (i < Nz-1))
            // {
            //     z=static_cast<double>(i)*c_dz;
            //     P_it[i]=c_P_0*exp(-((time-c_t0)/c_tau)*((time-c_t0)/c_tau))*exp(-z/c_delta);
            //     Tep1=Te[i] + (c_dt/(c_gamma_e*Te[i]))*(c_Gep*(Tp[i]-Te[i]) + P_it[i]+c_kappa_0*( (Te[i]/Tp[i]) * (Te[i+1]-2.0*Te[i]+Te[i-1])*c_oneOvrdzdz+(Tp[i]*((Te[i+1]-Te[i-1])*c_oneOvr2dz) - Te[i]*(Tp[i+1]-Tp[i-1])*c_oneOvr2dz)/(Tp[i]*Tp[i])*((Te[i+1]-Te[i-1])*c_oneOvr2dz)));
            //     Tpp1=Tp[i]+(c_dt*c_Gep/c_Cp)*(Te[i]-Tp[i]);
            // }

            //update the values of Te[i] and Tp[i]
            Te[i]=Tep1;
            Tp[i]=Tpp1;
        }
    }

    __global__ void ttf(double time, int N, double *dconst, double *dtfa, double *Te, int *dzlayer)
    {

        const int a = blockDim.x * blockIdx.x + threadIdx.x; 
        int siteincell = a % c_Nq;

        if (a < N){
            dtfa[a] = dconst[siteincell] * sqrt(Te[dzlayer[a]]);
        }
    }

    __global__ void ttfg(double time, int N, double *dconst, double *dtfa, double *Te, int *dlayer,  double grad)
    {

        const int a = blockDim.x * blockIdx.x + threadIdx.x; 
        int siteincell = a % c_Nq;

        if (a < N){

            dtfa[a] = dconst[siteincell] * sqrt(dlayer[a] * grad);
        }
    }


    void testing(int i){

		Array<double> testingx;
		testingx.resize(params::Lz);
		CUDA_CALL(cudaMemcpy(testingx.ptr(), Te, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));
	    std::cout << i << " ";
        std::cout << testingx(0) << " " << testingx(1) << std::endl;	
	}


    
}
