// cpp header files
#include <ctime>
#include <cuda.h>
#include <curand.h>
#include <iostream>
#include <curand_kernel.h>

// my header files
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/cudefine.h"

namespace cuthermal {



    //two temperature model variables
    double gamma_e;
    double Cp;
    double kappa_0;
    double delta;
    double Gep;
    double P_0;
    double t0;
    double tau;
    int    Nz;
    double dz;
    double dt;
    double Tinit;
    double oneOvrdzdz;
    double oneOvr2dz;

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
    __constant__ double c_Tinit;
    __constant__ double c_oneOvrdzdz;
    __constant__ double c_oneOvr2dz;
    __constant__ double c_thermal_const;

    // thermal gradient constants
    __constant__ double c_grad;

    double *Te, *Tp, *P_it;
	double *dtfa;
    int *dxlayer, *dylayer, *dzlayer;

    // stochastic noise variables
    float *gvalsx, *gvalsy, *gvalsz;
    curandGenerator_t gen;


    void init_cuthermal(double equilibrium_temp){

        //two temperature model variables
        gamma_e=125.0;                //gamma_e defines the electron specific heat through, C_e = gamma_e * T_e. [J/m^3/K^2]
        Cp=3e6;                       //Specific heat of phonons. [J/m^3/K]
        kappa_0=11.0;                 //kappa_0 defines the thermal heat conductivity (kappa) through, kappa = kappa_0 * T_e/T_p [J/m/K/s]
        delta=20.0e-9;                //Penetration depth of laser. [m]
        Gep=10e17;                    //Electron-phonon coupling [ J/m^3/s/K]
        P_0=2.0e21;                   //Pump fluence prefactor, P_0. P(z,t)=P_0*exp(-((t-t0)/tau)**2)*exp(-z/delta) [ J/m^3/s]
        t0=100e-15;                   //Pump temporal offset [s]
        tau=200e-15;                   //Pump temporal full width half max [s]
        Nz=100;                       //number of unit cells in z-direction (assumed uniform heating perpendicular [unit cells in z]
        dz=0.3e-9;                    //lattice constant (or difference between planes) [m]
        dt=1e-16;                     //Timestep [s]
        Tinit=equilibrium_temp;
        oneOvrdzdz=1./(dz*dz);
        oneOvr2dz=1./(2.0*dz);

        // Constants two temperature model
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_gamma_e), &gamma_e, sizeof(double)));                        
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Cp), &Cp, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_kappa_0), &kappa_0, sizeof(double)));                        
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_delta), &delta, sizeof(double)));                            
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Gep), &Gep, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_P_0), &P_0, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_t0), &t0, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_tau), &tau, sizeof(double)));                                
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Nz), &Nz, sizeof(int)));                                     
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_dz), &dz, sizeof(double)));                                  
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_dt), &params::dt, sizeof(double)));                          
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_thermal_const), &params::thermal_const, sizeof(double)));
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_Tinit), &Tinit, sizeof(double)));
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_oneOvrdzdz), &oneOvrdzdz, sizeof(double)));
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_oneOvr2dz), &oneOvr2dz, sizeof(double)));

        // Constants for thermal gradient
        CUDA_CALL(cudaMemcpyToSymbol(*(&c_grad), &params::temp_gradient, sizeof(double)));
    }

    void destroy_generator(){
        curandDestroyGenerator(gen);
        INFO_OUT("generator destroyed: ", "success");
    }

    void curand_generator(){
	std::time_t result = std::time(nullptr);
	int seed = static_cast<int>(result);
    std::cout << "time since epoch  = " << result << " (s)" << std::endl;
	CURAND_CALL(curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_MTGP32));
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, seed));
	std::cout << "Curand Seed = " << seed << std::endl;

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
            double z;

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

    __global__ void ttf(double time, int N, double *dtfa, double *Te, int *dzlayer)
    {

        const int a = blockDim.x * blockIdx.x + threadIdx.x; 

        if (a < N){
            dtfa[a] = c_thermal_const * sqrt(Te[dzlayer[a]]);
        }
    }

    __global__ void ttfg(double time, int N, double *dtfa, double *Te, int *dlayer,  double grad)
    {

        const int a = blockDim.x * blockIdx.x + threadIdx.x; 

        if (a < N){
            dtfa[a] = c_thermal_const * sqrt(dlayer[a] * grad);
        }
    }


    void testing(int i){

		Array<double> testingx;
		Array<double> testingy, testingz;
		testingx.resize(params::Nspins);
		testingy.resize(params::Lx);
		// testingz.resize(params::Nspins);

		CUDA_CALL(cudaMemcpy(testingx.ptr(), dtfa, sizeof(double) * params::Nspins, cudaMemcpyDeviceToHost));
		// CUDA_CALL(cudaMemcpy(testingy.ptr(), Te, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));
		// CUDA_CALL(cudaMemcpy(testingz.ptr(), Hapz, sizeof(double) * params::Lz, cudaMemcpyDeviceToHost));

	    std::cout << i << " ";
	    for (int a = 0; a < params::Lx; a++){
		    std::cout << std::scientific << a * params::a1 << " ";
            std::cout << testingx(geom::LatCount(a,2,2,0)) << std::endl;	
            // std::cout << testingy(a) << std::endl;	
	    }
	    std::cout << std::endl;
	}


    
}
