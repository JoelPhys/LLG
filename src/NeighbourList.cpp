
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <cstring>
#include "../inc/defines.h"
#include "../inc/spins.h"
#include "../inc/mathfuncs.h"
#include "../inc/fields.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/libconfig.h++"
#include "../inc/NeighbourList.h"

namespace neigh {

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

    void ReadFile(){

        TITLE("READING IN EXCHANGE FILE");

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
            std::cout << "Jij input file has been read in successfully" << std::endl;
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

        std::cout.width(75); std::cout << std::left << "length of x_adj:"; std::cout << x_adj.size() << std::endl;
        std::cout.width(75); std::cout << std::left << "length of adjncy:"; std::cout << adjncy.size() << std::endl;
        std::cout.width(75); std::cout << std::left << "average number of neighbours:"; std::cout << adjncy.size() / (x_adj.size()-1) << std::endl;
        std::cout.width(75); std::cout << std::left << "length of Jij:"; std::cout << Jijz_prime.size() << std::endl;
        std::cout.width(75); std::cout << std::left << "sum of neighbouring Jijs:"; std::cout << Jijsize << " (J)" << std::endl;
    }

    

}
