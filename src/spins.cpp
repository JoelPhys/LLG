#include <iostream>

#include "../inc/mathfuncs.h"
#include "../inc/array.h"
#include "../inc/spins.h"
#include "../inc/config.h"

namespace spins {

    // Global Variables
    Array<double> sx1d;
    Array<double> sy1d;
    Array<double> sz1d;

    // Testing for hedgehog
    Array<double> surfx;
    Array<double> surfy;
    Array<double> surfz;

    Array<int> lw;
    Array<int> rw;
    Array<int> zlayer;
    

    void init(){

        sx1d.resize(params::Nspins);
        sy1d.resize(params::Nspins);
        sz1d.resize(params::Nspins);
        sx1d.IFill(0);
        sy1d.IFill(0);
        sz1d.IFill(0);

    }

    void populate(){


		int count1d = 0;

		if (params::afmflag != "NiO"){
            for (int x = 0; x < params::Lx; x++){
                for (int y = 0; y < params::Ly; y++){
                    for (int z = 0; z < params::Lz; z++){
                        for (int q = 0; q < params::Nq; q++){

                            // testing domain walls
                            // if (x < params::Lx/2){
                            //     Sx4(x,y,z,q) = params::initm[q][0];
                            //     Sy4(x,y,z,q) = params::initm[q][1];
                            //     Sz4(x,y,z,q) = params::initm[q][2];                       
                            // }
                            // else if (x > params::Lx/2){
                            //     Sx4(x,y,z,q) = -1 * params::initm[q][0];
                            //     Sy4(x,y,z,q) = -1 * params::initm[q][1];
                            //     Sz4(x,y,z,q) = -1 * params::initm[q][2];                         
                            // }
                            // else if (x == params::Lx/2){
                            //     Sx4(x,y,z,q) = 0.0;
                            //     Sz4(x,y,z,q) = 1.0;                        
                            // }

                            sx1d(count1d + q) = params::initm[q][0];
                            sy1d(count1d + q) = params::initm[q][1];
                            sz1d(count1d + q) = params::initm[q][2];

                        }	
                        count1d += params::Nq; 
                    }
                }
            }
		}
		else {
            for (int a = 0; a < params::Nspins; a++){
                if (a / params::Nq % 2 == 0){
                    if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
                        sz1d(a) = 1.0; 
                    }
                    else if (modfunc(params::Nq,a) == 2) {
                        sz1d(a) = -1.0;
                    }
                    else {
                        std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
                        exit(0);
                    }
                }
                else if (a / params::Nq % 2 == 1) {
                    if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
                        sz1d(a) = -1.0; 
                    }
                    else if (modfunc(params::Nq,a) == 2) {
                        sz1d(a) = 1.0;
                    }
                    else {
                        std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
                        exit(0);
                    }
                }
            }
        }
    }


}