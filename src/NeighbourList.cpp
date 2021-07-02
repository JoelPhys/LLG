
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include "../inc/mathfuncs.h"
#include "../inc/fields.h"
#include "../inc/params1.h"
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/libconfig.h++"
#include "../inc/NeighbourList.h"

namespace neigh {

    // ==== INTIALISE HEUN ARRAYS ======================= //

    Array2D<double> H_thermal;
    Array2D<double> Delta_S;
    Array<double> Sx1d;
    Array<double> Sy1d;
    Array<double> Sz1d;
    Array<double> S_dash_normedx1d;
    Array<double> S_dash_normedy1d;
    Array<double> S_dash_normedz1d;


    double guassian_vals[3];
    double ani1[3];
    double H_uni[3];
    double H_uni_dash[3];
    double H_cub[3];
    double H_cub_dash[3];
    double H_new[3];
    double H_new_dash[3];
    double S_dash[3];
    double Delta_S_dash[3];
    double S_new[3];
    double CrossP1[3];
    double CrossP2[3];
    double CrossP3[3];
    double CrossP4[3];
    double invmag;
    double invmag1;
    double ScrossP[3];
    double H_exch[3];
    double H_exch_dash[3];
    int counting;
    std::normal_distribution<double> distribution(0.0,1.0);
    std::random_device device;
    std::mt19937 generator(device());
    // ================================================ //

    int length;
    std::vector<unsigned int> adjncy;
    std::vector<double> Jijy_prime;
    std::vector<double> Jijz_prime;
    std::vector<double> Jijx_prime;
    std::vector<unsigned int> x_adj;
    std::vector<double> Nx;
    std::vector<double> Ny;
    std::vector<double> Nz;
    std::vector<int> ib;
    std::vector<int> jb;

    std::vector<double> Jijx;
    std::vector<double> Jijy;
    std::vector<double> Jijz;
    

    void IntialisePointersNL(){
        H_thermal.resize(params::Nspins, 3);
        Delta_S.resize(params::Nspins, 3);
        H_thermal.IFill(0);
        Delta_S.IFill(0);

        Sx1d.resize(params::Nspins);
        Sy1d.resize(params::Nspins);
        Sz1d.resize(params::Nspins);
        S_dash_normedx1d.resize(params::Nspins);
        S_dash_normedy1d.resize(params::Nspins);
        S_dash_normedz1d.resize(params::Nspins);

        Sx1d.IFill(0);
        Sy1d.IFill(0);
        Sz1d.IFill(0);
        S_dash_normedx1d.IFill(0);
        S_dash_normedy1d.IFill(0);
        S_dash_normedz1d.IFill(0);     
    }

    void ReadFile(){

        std::ifstream input(params::Jij_filename);
        if (!input){
            std::cout << "ERROR: Could not open file" << std::endl;
            exit(0);
        }

        double sum1;

        if (params::format == "Jerome"){
            double a, b, c, d, e, f;
            double count = 0;

            while (input >> a >> b >> c >> d >> e >> f)
            {
                if (params::JijCutoff == false) {
                    ib.push_back(a);
                    jb.push_back(b);
                    Nx.push_back(c);
                    Ny.push_back(d);
                    Nz.push_back(e);
                    Jijx.push_back(f);
                    Jijy.push_back(f);
                    Jijz.push_back(f);
                }
                else if (params::JijCutoff == true) {
                    if (std::abs(f) < params::Jij_min) { 

                    }
                    else {
                        ib.push_back(a);
                        jb.push_back(b);
                        Nx.push_back(c);
                        Ny.push_back(d);
                        Nz.push_back(e);
                        Jijx.push_back(f);
                        Jijy.push_back(f);
                        Jijz.push_back(f);
                        count++;
                    }
                }
                else {
                    std::cout << "WARNING: Unasigned Jij cutoff flag" << std::endl;
                }
            }
            std::cout << "Jij input file has been read" << std::endl;
        }
        else if (params::format == "diag"){
            double a, b, c, d, e, f, g, h;
            double count = 0;

            while (input >> a >> b >> c >> d >> e >> f >> g >> h)
            {
                if (params::JijCutoff == false) {
                    ib.push_back(a);
                    jb.push_back(b);
                    Nx.push_back(c);
                    Ny.push_back(d);
                    Nz.push_back(e);
                    Jijx.push_back(f);
                    Jijy.push_back(g);
                    Jijz.push_back(h);
                }
                else if (params::JijCutoff == true) {
                    if (std::abs(f) < params::Jij_min) { 

                    }
                    else {
                        ib.push_back(a);
                        jb.push_back(b);
                        Nx.push_back(c);
                        Ny.push_back(d);
                        Nz.push_back(e);
                        Jijx.push_back(f);
                        Jijy.push_back(g);
                        Jijz.push_back(h);
                        count++;
                    }
                }
                else {
                    std::cout << "WARNING: Unasigned Jij cutoff flag" << std::endl;
                }
            }
            std::cout << "Jij input file has been read" << std::endl;
        }
	else {
		std::cout << "ERROR: Unknown exchange file type " << std::endl;
		exit(0);
	}

        length = Jijx.size();
        
        if (params::Jijhalf == true) {
            for (int i = 0; i < length; i++){
                Jijx[i] = Jijx[i] * 2;
                Jijy[i] = Jijy[i] * 2;
                Jijz[i] = Jijz[i] * 2;
            }
            std::cout << "Jij values have been doubled" << std::endl;
        }
        else if (params::Jijhalf == false) {
            std::cout << "Jij values have NOT been doubled" << std::endl;

        }
        else {
            std::cout << "WARNING: Unasigned Jij double flag" << std::endl;
        }

        for (int i = 0; i < length; i++){
            sum1 += std::abs(Jijx[i]);
        }

        std::cout << "Sum of Jij = " << sum1 << std::endl;
        input.close();
    }

    void InteractionMatrix() {

        int NxP[length];
        int NyP[length];
        int NzP[length];
        double vecX, vecY, vecZ;

        if (params::changesign == "Y"){
            for (int i = 0; i < length; i++){
                if (ib[i] == jb[i]){
                    Jijx[i] = Jijx[i];
                    Jijy[i] = Jijy[i];
                    Jijz[i] = Jijz[i];
                }
                else if (((ib[i] == 3) && (jb[i] == 5)) || ((ib[i] == 5) && (jb[i] == 3))) {
                    Jijx[i] = Jijx[i];
                    Jijy[i] = Jijy[i];
                    Jijz[i] = Jijz[i];
                }
                else if (((ib[i] == 4) && (jb[i] == 6)) || ((ib[i] == 6) && (jb[i] == 4))) {
                    Jijx[i] = Jijx[i];
                    Jijy[i] = Jijy[i];
                    Jijz[i] = Jijz[i];
                }
                else  {
                    Jijx[i] = -1 * Jijx[i];
                    Jijy[i] = -1 * Jijy[i];
                    Jijz[i] = -1 * Jijz[i];
                }
            }
        }
        else if (params::changesign == "N"){
        }
        else if (params::changesign == "ferromagnetic"){
            for (int i = 0; i < length; i++){
                Jijx[i] = std::abs(Jijx[i]);
                Jijy[i] = std::abs(Jijy[i]);
                Jijz[i] = std::abs(Jijz[i]);
            }
        }
        else {
            std::cout << "WARNING: unassigned changesign flag." << std::endl;
        }

        // Find inverse of lattice vectors
        Inverse3x3(params::Plat, params::PlatINV);
	
        // Convert to unit cell vector
        for (int i = 0; i < length; i++){
            
            vecX = Nx[i] + params::sites[ib[i]+params::ibtoq][0] - params::sites[jb[i]+params::ibtoq][0];
            vecY = Ny[i] + params::sites[ib[i]+params::ibtoq][1] - params::sites[jb[i]+params::ibtoq][1];
            vecZ = Nz[i] + params::sites[ib[i]+params::ibtoq][2] - params::sites[jb[i]+params::ibtoq][2];

            NxP[i] = nearbyint((params::PlatINV[0][0] * vecX) + (params::PlatINV[0][1] * vecY) + (params::PlatINV[0][2] * vecZ));
            NyP[i] = nearbyint((params::PlatINV[1][0] * vecX) + (params::PlatINV[1][1] * vecY) + (params::PlatINV[1][2] * vecZ));
            NzP[i] = nearbyint((params::PlatINV[2][0] * vecX) + (params::PlatINV[2][1] * vecY) + (params::PlatINV[2][2] * vecZ));

        }
        
	int xval;
        int yval;
        int zval;
        int qval;
        int adjcounter = 0;

        x_adj.push_back(0);

        // ========== Neighbour List =========== //
        if (params::Jij_units == "mRy") {
            for (int x = 0; x < params::Lx; ++x){                // Depth
                for (int y = 0; y < params::Ly; ++y){            // Row
                    for (int z = 0; z < params::Lz; ++z){        // Column
                        for (int q = 0; q < params::Nq; ++q){    // Unit Cell
                            for (int i = 0; i < length; i++){

                                if (q == ib[i]+params::ibtoq){

                                    if (params::xbound == "periodic"){
                                        xval = modfunc(params::Lx, NxP[i] + x);
                                    }
                                    else if (params::xbound == "open"){
                                        xval = NxP[i] + x;
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown xbound " << std::endl;
                                        exit(0);
                                    }
                                    if (params::ybound == "periodic"){
                                        yval = modfunc(params::Ly, NyP[i] + y);
                                    }
                                    else if (params::ybound == "open"){
                                        yval = NyP[i] + y; 
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown ybound " << std::endl;
                                        exit(0);
                                    }
                                    if (params::zbound == "periodic"){
                                        zval = modfunc(params::Lz, NzP[i] + z);
                                    }
                                    else if (params::zbound == "open"){
                                        zval = NzP[i] + z; 
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown zbound " << std::endl;
                                        exit(0);
                                    }

                                    qval = jb[i]+params::ibtoq;

                                    if (((xval >= 0) && (yval >= 0) && (zval >= 0)) && ((xval < params::Lx) && (yval < params::Ly) && (zval < params::Lz))) {
                                        adjncy.push_back(geom::LatCount(xval, yval, zval, qval));
                                        adjcounter++;

                                        Jijx_prime.push_back( ( Jijx[i]  * 2.179872e-21) / params::mu_s);
                                        Jijy_prime.push_back( ( Jijy[i]  * 2.179872e-21) / params::mu_s);
                                        Jijz_prime.push_back( ( Jijz[i]  * 2.179872e-21) / params::mu_s);
                                    }
                                }
                            }
                            x_adj.push_back(adjcounter);
                        }
                    }
                }
            }  
        }
        else if (params::Jij_units == "J") {
            for (int x = 0; x < params::Lx; ++x){                // Depth
                for (int y = 0; y < params::Ly; ++y){            // Row
                    for (int z = 0; z < params::Lz; ++z){        // Column
                        for (int q = 0; q < params::Nq; ++q){    // Unit Cell
                            for (int i = 0; i < length; i++){

                                if (q == ib[i]+params::ibtoq){

                                    if (params::xbound == "periodic"){
                                        xval = modfunc(params::Lx, NxP[i] + x);
                                    }
                                    else if (params::xbound == "open"){
                                        xval = NxP[i] + x;
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown xbound " << std::endl;
                                        exit(0);
                                    }
                                    if (params::ybound == "periodic"){
                                        yval = modfunc(params::Ly, NyP[i] + y);
                                    }
                                    else if (params::ybound == "open"){
                                        yval = NyP[i] + y; 
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown ybound " << std::endl;
                                        exit(0);
                                    }
                                    if (params::zbound == "periodic"){
                                        zval = modfunc(params::Lz, NzP[i] + z);
                                    }
                                    else if (params::zbound == "open"){
                                        zval = NzP[i] + z; 
                                    }
                                    else {
                                        std::cout << "ERROR: Unknown zbound " << std::endl;
                                        exit(0);
                                    }

                                    qval = jb[i]+params::ibtoq;

                                    if (((xval >= 0) && (yval >= 0) && (zval >= 0)) && ((xval < params::Lx) && (yval < params::Ly) && (zval < params::Lz))) {
                                        adjncy.push_back(geom::LatCount(xval, yval, zval, qval));
                                        adjcounter++;

                                        Jijx_prime.push_back( ( Jijx[i] ) / params::mu_s);
                                        Jijy_prime.push_back( ( Jijy[i] ) / params::mu_s);
                                        Jijz_prime.push_back( ( Jijz[i] ) / params::mu_s);
                                    }
                                }
                            }
                            x_adj.push_back(adjcounter);
                        }
                    }
                }
            }  
        }
        else {
            std::cout << "ERROR: Unknown Jij units" << std::endl;
            exit(0);
        }

        std::cout << x_adj[0] << " " << x_adj[1] << " " << x_adj[2] <<  std::endl;

        double Jijsize = 0;
        for (int i = 0; i < adjncy.size() / (x_adj.size()-1); i++){
            Jijsize += (Jijx[i]);
        }

        std::cout << "length of x_adj = " << x_adj.size() << std::endl;
        std::cout << "length of adjncy = " << adjncy.size() << std::endl;
        std::cout << "average number of neighbours = " << adjncy.size() / (x_adj.size()-1) << std::endl;
        std::cout << "length of Jij = " << Jijz_prime.size() << std::endl;
        std::cout << "sum of neighbouring Jijs = " << Jijsize << " (J)" << std::endl;
    }

    void Heun(double Thermal_Fluct){

        for (int a = 0; a < params::Nspins; a++){

            for (int k = 0; k < 3; k++){
                guassian_vals[k] = distribution(generator);
            }

            H_thermal(a,0) = guassian_vals[0] * Thermal_Fluct;
            H_thermal(a,1) = guassian_vals[1] * Thermal_Fluct;
            H_thermal(a,2) = guassian_vals[2] * Thermal_Fluct;

            // Uniaxial anisotropy
            H_uni[0] = params::dxup * Sx1d(a);
            H_uni[1] = params::dyup * Sy1d(a);
            H_uni[2] = params::dzup * Sz1d(a);

            // Cubic Anisotropy
            H_cub[0] = params::dxcp * Sx1d(a) * Sx1d(a) * Sx1d(a);
            H_cub[1] = params::dycp * Sy1d(a) * Sy1d(a) * Sy1d(a);
            H_cub[2] = params::dzcp * Sz1d(a) * Sz1d(a) * Sz1d(a);

            // Exchange interaction
            H_exch[0] = 0;
            H_exch[1] = 0;
            H_exch[2] = 0;

            counting = x_adj[a];

            for (int b = x_adj[a]; b < x_adj[a+1]; b++){
                H_exch[0] += Jijx_prime[counting] * (Sx1d(adjncy[b]));
                H_exch[1] += Jijy_prime[counting] * (Sy1d(adjncy[b]));
                H_exch[2] += Jijz_prime[counting] * (Sz1d(adjncy[b]));
                counting++;
            }

            H_new[0] = H_thermal(a,0) + fields::H_appx(a) + H_uni[0] + H_cub[0] + H_exch[0];
            H_new[1] = H_thermal(a,1) + fields::H_appy(a) + H_uni[1] + H_cub[1] + H_exch[1];
            H_new[2] = H_thermal(a,2) + fields::H_appz(a) + H_uni[2] + H_cub[2] + H_exch[2];

            ScrossP[0] = Sx1d(a);
            ScrossP[1] = Sy1d(a);
            ScrossP[2] = Sz1d(a);

            CrossP(ScrossP, H_new, CrossP1);
            CrossP(ScrossP, CrossP1, CrossP2);

            Delta_S(a,0) = -params::lambdaPrime * (CrossP1[0] + params::lambda * CrossP2[0]);
            Delta_S(a,1) = -params::lambdaPrime * (CrossP1[1] + params::lambda * CrossP2[1]);
            Delta_S(a,2) = -params::lambdaPrime * (CrossP1[2] + params::lambda * CrossP2[2]);

            S_dash[0] = Sx1d(a) + (Delta_S(a,0) * params::dtau);
            S_dash[1] = Sy1d(a) + (Delta_S(a,1) * params::dtau);
            S_dash[2] = Sz1d(a) + (Delta_S(a,2) * params::dtau);

            invmag = 1 / sqrt(S_dash[0] * S_dash[0] + S_dash[1] * S_dash[1] + S_dash[2] * S_dash[2]);
            S_dash_normedx1d(a) = invmag * S_dash[0];
            S_dash_normedy1d(a) = invmag * S_dash[1];
            S_dash_normedz1d(a) = invmag * S_dash[2];   
        }

        for (int a = 0; a < params::Nspins; a++){
            
            //  Uniaxial Anisototropy
            H_uni_dash[0]= params::dxup * S_dash_normedx1d(a);
            H_uni_dash[1]= params::dyup * S_dash_normedy1d(a);
            H_uni_dash[2]= params::dzup * S_dash_normedz1d(a);

            // Cubic Ansisotropy
            H_cub_dash[0]= params::dxcp * S_dash_normedx1d(a) * S_dash_normedx1d(a) * S_dash_normedx1d(a);
            H_cub_dash[1]= params::dycp * S_dash_normedy1d(a) * S_dash_normedy1d(a) * S_dash_normedy1d(a);
            H_cub_dash[2]= params::dzcp * S_dash_normedz1d(a) * S_dash_normedz1d(a) * S_dash_normedz1d(a);

            H_exch_dash[0] = 0;
            H_exch_dash[1] = 0;
            H_exch_dash[2] = 0;

            // Exchange interaction prime
            counting = x_adj[a];
            for (int b = x_adj[a]; b < x_adj[a+1]; b++){
                H_exch_dash[0] += Jijx_prime[counting] * (S_dash_normedx1d(adjncy[b]));
                H_exch_dash[1] += Jijy_prime[counting] * (S_dash_normedy1d(adjncy[b]));
                H_exch_dash[2] += Jijz_prime[counting] * (S_dash_normedz1d(adjncy[b]));
                counting++;
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

            Delta_S_dash[0] = -params::lambdaPrime * (CrossP3[0] + params::lambda * CrossP4[0]);
            Delta_S_dash[1] = -params::lambdaPrime * (CrossP3[1] + params::lambda * CrossP4[1]);
            Delta_S_dash[2] = -params::lambdaPrime * (CrossP3[2] + params::lambda * CrossP4[2]); 

            S_new[0] = Sx1d(a) + params::half_dtau * (Delta_S(a,0) + Delta_S_dash[0]);
            S_new[1] = Sy1d(a) + params::half_dtau * (Delta_S(a,1) + Delta_S_dash[1]);
            S_new[2] = Sz1d(a) + params::half_dtau * (Delta_S(a,2) + Delta_S_dash[2]);

            invmag1 = 1 / sqrt(S_new[0] * S_new[0] + S_new[1] * S_new[1] + S_new[2] * S_new[2]);
            Sx1d(a) = invmag1 * S_new[0];
            Sy1d(a) = invmag1 * S_new[1];
            Sz1d(a) = invmag1 * S_new[2];
        }
        // double sumx, sumy, sumz;
        // for (int i = 0; i < params::Nspins; i++){
        //     sumx += H_thermal(i,0);
        //     sumy += H_thermal(i,1);
        //     sumz += H_thermal(i,2);
        // }
        // std::cout << sumx << " " << sumy << " " << sumz << std::endl;
    }

}
