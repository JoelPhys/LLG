//cpp header files
#include <iostream>
#include <cmath>
#include <random>
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/spins.h"
#include "../inc/config.h"
#include "../inc/mathfuncs.h"

namespace spins {

    // Global Variables
    Array<double> sx1d;
    Array<double> sy1d;
    Array<double> sz1d;

    // Total Energy
    Array<double> Ex;
    Array<double> Ey;
    Array<double> Ez;

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


        Ex.resize(params::Nspins);
        Ey.resize(params::Nspins);
        Ez.resize(params::Nspins);
        Ex.IFill(0);
        Ey.IFill(0);
        Ez.IFill(0);
    }

    void populate(){


		int count1d = 0;
        if (params::afmflag == "test"){

            //Normal Distribution for Stochastic noise
            std::normal_distribution<double> distribution(0.0,1.0);
            std::random_device device;
            std::mt19937 generator(device());
    
            for (int i = 0; i < params::Nspins; i++){
                double v1=0,v2=0,s=2.0,ss=0.0;
                while(s>1.0)
                {
                    v1=2.0*distribution(generator)-1.0;
                    v2=2.0*distribution(generator)-1.0;
                    s=v1*v1+v2*v2;
                }
                ss=sqrt(1.0-s);
                sx1d(i)=2.0*v1*ss;
                sy1d(i)=2.0*v2*ss;
                sz1d(i)=1.0-2.0*s;
            }
        }
		else if (params::afmflag != "NiO"){
            for (int x = 0; x < params::Lx; x++){
                for (int y = 0; y < params::Ly; y++){
                    for (int z = 0; z < params::Lz; z++){
                        for (int q = 0; q < params::Nq; q++){

                            // testing domain walls
                            //if (x < params::Lx/2){
                            //    sx1d(geom::LatCount(x,y,z,q)) = params::initm[q][0];
                            //    sy1d(geom::LatCount(x,y,z,q)) = params::initm[q][1];
                            //    sz1d(geom::LatCount(x,y,z,q)) = params::initm[q][2];                       
                            //}
                            //else if (x > params::Lx/2){
                            //    sx1d(geom::LatCount(x,y,z,q)) = -1 * params::initm[q][0];
                            //    sy1d(geom::LatCount(x,y,z,q)) = -1 * params::initm[q][1];
                            //    sz1d(geom::LatCount(x,y,z,q)) = -1 * params::initm[q][2];                         
                            //}
                            //else if (x == params::Lx/2){
                            //    sx1d(geom::LatCount(x,y,z,q)) = 0.5;
                            //    sz1d(geom::LatCount(x,y,z,q)) = 0.5;                        
                            //}

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
