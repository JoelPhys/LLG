#include <curand.h>
#include <cuda.h>
#include <curand_kernel.h>
#include "../inc/cudefine.h"
#include "../inc/params1.h"
#include <iostream>
#include <ctime>

namespace cuthermal {

    float *gvalsx, *gvalsy, *gvalsz;
    curandGenerator_t gen;

    void curand_generator(){
	std::time_t result = std::time(nullptr);
	int seed = static_cast<int>(result);
        std::cout << "time since epoch  = " << result << " (s)" << std::endl;
	CURAND_CALL(curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_MTGP32));
        CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, seed));
	std::cout << "Curand Seed = " << seed << std::endl;
    }

    void gen_thermal_noise(){
        CURAND_CALL(curandGenerateNormal(gen, gvalsx, params::Nspins, 0.0, 1.0));
        CURAND_CALL(curandGenerateNormal(gen, gvalsy, params::Nspins, 0.0, 1.0));
        CURAND_CALL(curandGenerateNormal(gen, gvalsz, params::Nspins, 0.0, 1.0));
    }

    void destroy_generator(){
        curandDestroyGenerator(gen);
	std::cout << "Curand Generator Destroyed" << std::endl;
    }


}
