

#include <cmath>
#include <vector>
#include <random>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/array.h"
#include "../inc/fields.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/thermal.h"
#include "../inc/mathfuncs.h"
#include "../inc/libconfig.h++"
#include "../inc/neighbourlist.h"

namespace heun {

    // ============================== INTIALISE HEUN ARRAYS ======================= //

    Array2D<double> Delta_S;
    Array2D<double> H_thermal;
    Array<double> S_dash_normedx1d, S_dash_normedy1d, S_dash_normedz1d;

    // function variables
    int a;
    int counting;
	int siteincell;
    double S_new[3];
    double invmag, invmag1;
    double guassian_vals[3];   
    double H_uni[3], H_uni_dash[3];
    double H_cub[3], H_cub_dash[3];
    double H_new[3], H_new_dash[3];
    double H_exch[3], H_exch_dash[3];
    double S_dash[3],  Delta_S_dash[3];
    double ScrossP[3], CrossP1[3], CrossP2[3], CrossP3[3], CrossP4[3];

    // Random number generation for stochastic noise
    std::normal_distribution<double> distribution(0.0,1.0);
    std::random_device device;
    std::mt19937 generator(device());

    // ============================================================================= //

    void init(){
      
        
        
        Delta_S.resize(params::Nspins, 3);
        H_thermal.resize(params::Nspins, 3);
		Delta_S.IFill(0);
		H_thermal.IFill(0);         
   
   		S_dash_normedx1d.resize(params::Nspins);
        S_dash_normedy1d.resize(params::Nspins);
        S_dash_normedz1d.resize(params::Nspins);
        S_dash_normedx1d.IFill(0);
        S_dash_normedy1d.IFill(0);
        S_dash_normedz1d.IFill(0);   

    }

    void integration(){

        for (int c = 0; c < neigh::simspin.size(); c++){
			
            a = neigh::simspin[c];
		
    		siteincell = a % params::Nq;
			
			for (int k = 0; k < 3; k++){
                guassian_vals[k] = distribution(generator);
            }

			H_thermal(a,0) = guassian_vals[0] * params::thermal_const[siteincell] * sqrt(thermal::Te[geom::zlayer[a]]);
            H_thermal(a,1) = guassian_vals[1] * params::thermal_const[siteincell] * sqrt(thermal::Te[geom::zlayer[a]]);
            H_thermal(a,2) = guassian_vals[2] * params::thermal_const[siteincell] * sqrt(thermal::Te[geom::zlayer[a]]);
		
			// Uniaxial anisotropy
            H_uni[0] = params::dxup[siteincell] * spins::sx1d(a);
            H_uni[1] = params::dyup[siteincell] * spins::sy1d(a);
            H_uni[2] = params::dzup[siteincell] * spins::sz1d(a);

            // Cubic Anisotropy
            //H_cub[0] = params::dzcp * spins::sx1d(a) * spins::sx1d(a) * spins::sx1d(a);
            //H_cub[1] = params::dzcp * spins::sy1d(a) * spins::sy1d(a) * spins::sy1d(a);
            //H_cub[2] = params::dzcp * spins::sz1d(a) * spins::sz1d(a) * spins::sz1d(a);
			//H_cub[0] = 2.0 *  0.1787343119 * spins::sx1d(a) * spins::sy1d(a) * spins::sy1d(a); // 2 * 0.04  meV / mu_s 
			//H_cub[1] = 2.0 *  0.1787343119 * spins::sy1d(a) * spins::sx1d(a) * spins::sx1d(a); // 2 * 0.04  meV / mu_s 
			//H_cub[2] = 4.0 * -0.0670253669 * spins::sz1d(a) * spins::sz1d(a) * spins::sz1d(a); // 4 * 0.015 meV / mu_s
			H_cub[0] = 0.0;
			H_cub[1] = 0.0;
			H_cub[2] = 0.0;
            
			// Exchange interaction
            H_exch[0] = 0;
            H_exch[1] = 0;
            H_exch[2] = 0;

            // counting = neigh::x_adj[a];

            for (int b = neigh::x_adj[c]; b < neigh::x_adj[c+1]; b++){

                // H_exch[0] += neigh::Jijx_prime[counting] * (spins::sx1d(neigh::adjncy[b]));
                // H_exch[1] += neigh::Jijy_prime[counting] * (spins::sy1d(neigh::adjncy[b]));
                // H_exch[2] += neigh::Jijz_prime[counting] * (spins::sz1d(neigh::adjncy[b]));
                H_exch[0] += neigh::Jijx[neigh::jind[b]] * (spins::sx1d(neigh::adjncy[b]));
                H_exch[1] += neigh::Jijy[neigh::jind[b]] * (spins::sy1d(neigh::adjncy[b]));
                H_exch[2] += neigh::Jijz[neigh::jind[b]] * (spins::sz1d(neigh::adjncy[b]));
                // counting++;
            }
            
			H_new[0] = H_thermal(a,0) + fields::H_appx(a) + H_uni[0] + H_cub[0] + H_exch[0];
            H_new[1] = H_thermal(a,1) + fields::H_appy(a) + H_uni[1] + H_cub[1] + H_exch[1];
            H_new[2] = H_thermal(a,2) + fields::H_appz(a) + H_uni[2] + H_cub[2] + H_exch[2];

            ScrossP[0] = spins::sx1d(a);
            ScrossP[1] = spins::sy1d(a);
            ScrossP[2] = spins::sz1d(a);

            CrossP(ScrossP, H_new, CrossP1);
            CrossP(ScrossP, CrossP1, CrossP2);

            Delta_S(a,0) = -params::lambdaPrime[siteincell] * (CrossP1[0] + params::lambda[siteincell]* CrossP2[0]);
            Delta_S(a,1) = -params::lambdaPrime[siteincell] * (CrossP1[1] + params::lambda[siteincell]* CrossP2[1]);
            Delta_S(a,2) = -params::lambdaPrime[siteincell] * (CrossP1[2] + params::lambda[siteincell]* CrossP2[2]);

            S_dash[0] = spins::sx1d(a) + (Delta_S(a,0) * params::dtau);
            S_dash[1] = spins::sy1d(a) + (Delta_S(a,1) * params::dtau);
            S_dash[2] = spins::sz1d(a) + (Delta_S(a,2) * params::dtau);

            invmag = 1 / sqrt(S_dash[0] * S_dash[0] + S_dash[1] * S_dash[1] + S_dash[2] * S_dash[2]);
            S_dash_normedx1d(a) = invmag * S_dash[0];
            S_dash_normedy1d(a) = invmag * S_dash[1];
            S_dash_normedz1d(a) = invmag * S_dash[2];   
        }

        for (int c = 0; c < neigh::simspin.size(); c++){

            a = neigh::simspin[c];

			siteincell = a % params::Nq;

            //  Uniaxial Anisototropy
            H_uni_dash[0]= params::dxup[siteincell] * S_dash_normedx1d(a);
            H_uni_dash[1]= params::dyup[siteincell] * S_dash_normedy1d(a);
            H_uni_dash[2]= params::dzup[siteincell] * S_dash_normedz1d(a);

            // Cubic Ansisotropy
            //H_cub_dash[0]= params::dzcp * S_dash_normedx1d(a) * S_dash_normedx1d(a) * S_dash_normedx1d(a);
            //H_cub_dash[1]= params::dzcp * S_dash_normedy1d(a) * S_dash_normedy1d(a) * S_dash_normedy1d(a);
            //H_cub_dash[2]= params::dzcp * S_dash_normedz1d(a) * S_dash_normedz1d(a) * S_dash_normedz1d(a);
			//H_cub_dash[0] = 2.0 *  0.1787343119 * S_dash_normedx1d(a) * S_dash_normedy1d(a) * S_dash_normedy1d(a); // 2 * 0.04  meV / mu_s 
			//H_cub_dash[1] = 2.0 *  0.1787343119 * S_dash_normedy1d(a) * S_dash_normedx1d(a) * S_dash_normedx1d(a); // 2 * 0.04  meV / mu_s 
			//H_cub_dash[2] = 4.0 * -0.0670253669 * S_dash_normedz1d(a) * S_dash_normedz1d(a) * S_dash_normedz1d(a); // 4 * 0.015 meV / mu_s
			H_cub_dash[0] = 0.0;
			H_cub_dash[1] = 0.0;
			H_cub_dash[2] = 0.0;
            
			H_exch_dash[0] = 0;
            H_exch_dash[1] = 0;
            H_exch_dash[2] = 0;

            // Exchange interaction prime
            // counting = neigh::x_adj[a];
            for (int b = neigh::x_adj[c]; b < neigh::x_adj[c+1]; b++){
                // H_exch_dash[0] += neigh::Jijx_prime[counting] * (S_dash_normedx1d(neigh::adjncy[b]));
                // H_exch_dash[1] += neigh::Jijy_prime[counting] * (S_dash_normedy1d(neigh::adjncy[b]));
                // H_exch_dash[2] += neigh::Jijz_prime[counting] * (S_dash_normedz1d(neigh::adjncy[b]));
                H_exch_dash[0] += neigh::Jijx[neigh::jind[b]] * (S_dash_normedx1d(neigh::adjncy[b]));
                H_exch_dash[1] += neigh::Jijy[neigh::jind[b]] * (S_dash_normedy1d(neigh::adjncy[b]));
                H_exch_dash[2] += neigh::Jijz[neigh::jind[b]] * (S_dash_normedz1d(neigh::adjncy[b]));
                // counting++;
            }

            H_new_dash[0] = H_thermal(a,0) + fields::H_appx(a) + H_uni_dash[0] + H_cub_dash[0] + H_exch_dash[0];
            H_new_dash[1] = H_thermal(a,1) + fields::H_appy(a) + H_uni_dash[1] + H_cub_dash[1] + H_exch_dash[1];
            H_new_dash[2] = H_thermal(a,2) + fields::H_appz(a) + H_uni_dash[2] + H_cub_dash[2] + H_exch_dash[2];

            // Calculate Corrector and Normalise

            ScrossP[0] = S_dash_normedx1d(a);
            ScrossP[1] = S_dash_normedy1d(a);
            ScrossP[2] = S_dash_normedz1d(a);

            CrossP(ScrossP, H_new_dash, CrossP3);
            CrossP(ScrossP, CrossP3, CrossP4);

            Delta_S_dash[0] = -params::lambdaPrime[siteincell] * (CrossP3[0] + params::lambda[siteincell]* CrossP4[0]);
            Delta_S_dash[1] = -params::lambdaPrime[siteincell] * (CrossP3[1] + params::lambda[siteincell]* CrossP4[1]);
            Delta_S_dash[2] = -params::lambdaPrime[siteincell] * (CrossP3[2] + params::lambda[siteincell]* CrossP4[2]); 

            S_new[0] = spins::sx1d(a) + params::half_dtau * (Delta_S(a,0) + Delta_S_dash[0]);
            S_new[1] = spins::sy1d(a) + params::half_dtau * (Delta_S(a,1) + Delta_S_dash[1]);
            S_new[2] = spins::sz1d(a) + params::half_dtau * (Delta_S(a,2) + Delta_S_dash[2]);

            invmag1 = 1 / sqrt(S_new[0] * S_new[0] + S_new[1] * S_new[1] + S_new[2] * S_new[2]);
            spins::sx1d(a) = invmag1 * S_new[0];
            spins::sy1d(a) = invmag1 * S_new[1];
            spins::sz1d(a) = invmag1 * S_new[2];
        }
    }
}
