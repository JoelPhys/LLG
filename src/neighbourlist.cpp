// cpp header files
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <cstring>

// my header files
#include "../inc/defines.h"
#include "../inc/spins.h"
#include "../inc/mathfuncs.h"
#include "../inc/fields.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/error.h"
#include "../inc/defects.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/libconfig.h++"
#include "../inc/neighbourlist.h"

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

    //testing
    std::vector<int> jind;

    std::vector<double> Jijx;
    std::vector<double> Jijy;
    std::vector<double> Jijz;
    
    std::vector<int> simspin;
    int nsimspin;

    void ReadFile(){

        TITLE("EXCHANGE INFORMATION");
		INFO_OUT("Exchange filename: ", params::Jij_filename);        
		INFO_OUT("Exhchange Cutoff:", params::JijCutoff);
		INFO_OUT("Exhchange Energy Minimum:", params::Jij_min << " (" << params::Jij_units << ")");
		if (params::Jijhalf == true) {INFO_OUT("Have Jij values been doubled", "Yes");}
		else if (params::Jijhalf == false) {INFO_OUT("Have Jij values been doubled", "No");}


        std::ifstream input(params::Jij_filename);
        if (!input){
            std::cout << "ERROR: Could not open file" << std::endl;
            exit(0);
        }

        double sum1;

		std::string line;
		getline(input, line);         //read the first line of your file to string
		std::stringstream s;
		s << line;                   //send the line to the stringstream object...
		
		int how_many_cols = 0;    
		double value;
		
		while(s >> value) how_many_cols++;  //while there's something in the line, increase the number of columns
	
		// reset input back to first line
        input.clear();
        input.seekg(0);


        if (params::format == "Jerome"){
            double a, b, c, d, e, f;
            double count = 0;

			if (how_many_cols != 6){
				std::cout << "ERROR: Incorrect number of cols in Jij file. Exiting" << std::endl;
				exit(0);
			}

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
            INFO_OUT("Has the Jij input file has been read in successfully:","Yes");
        }
        else if (params::format == "diag"){
            double a, b, c, d, e, f, g, h;
            double count = 0;

			if (how_many_cols != 8){
				std::cout << "ERROR: Incorrect number of cols in Jij file. Exiting" << std::endl;
				exit(0);
			}

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
            INFO_OUT("Have Jij values been doubled","Yes");
        }
        else if (params::Jijhalf == false) {
            INFO_OUT("Have Jij values been doubled","No");
        }
        else {
            std::cout << "WARNING: Unasigned Jij double flag" << std::endl;
        }
    }

    void InteractionMatrix() {

        double NxP[length];
        double NyP[length];
        double NzP[length];
        double vecX, vecY, vecZ;
        
        // variables to check if Jij vector sits within integer unit cell
        double near[3];
        double close_to_int[3];

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

            if (params::Jij_units == "mRy") {
                Jijx[i] *= (2.179872e-21) / params::mu_s[ib[i]+params::ibtoq];
                Jijy[i] *= (2.179872e-21) / params::mu_s[ib[i]+params::ibtoq];
                Jijz[i] *= (2.179872e-21) / params::mu_s[ib[i]+params::ibtoq];
            }
            else if (params::Jij_units == "J") {
                Jijx[i] *= (1.0) / params::mu_s[ib[i]+params::ibtoq];
                Jijy[i] *= (1.0) / params::mu_s[ib[i]+params::ibtoq];
                Jijz[i] *= (1.0) / params::mu_s[ib[i]+params::ibtoq];
            }   

            vecX = Nx[i] + params::sites[ib[i]+params::ibtoq][0] - params::sites[jb[i]+params::ibtoq][0];
            vecY = Ny[i] + params::sites[ib[i]+params::ibtoq][1] - params::sites[jb[i]+params::ibtoq][1];
            vecZ = Nz[i] + params::sites[ib[i]+params::ibtoq][2] - params::sites[jb[i]+params::ibtoq][2];

            NxP[i] = (params::PlatINV[0][0] * vecX) + (params::PlatINV[0][1] * vecY) + (params::PlatINV[0][2] * vecZ);
            NyP[i] = (params::PlatINV[1][0] * vecX) + (params::PlatINV[1][1] * vecY) + (params::PlatINV[1][2] * vecZ);
            NzP[i] = (params::PlatINV[2][0] * vecX) + (params::PlatINV[2][1] * vecY) + (params::PlatINV[2][2] * vecZ);

            //NxP[i] = nearbyint((params::PlatINV[0][0] * vecX) + (params::PlatINV[0][1] * vecY) + (params::PlatINV[0][2] * vecZ));
            //NyP[i] = nearbyint((params::PlatINV[1][0] * vecX) + (params::PlatINV[1][1] * vecY) + (params::PlatINV[1][2] * vecZ));
            //NzP[i] = nearbyint((params::PlatINV[2][0] * vecX) + (params::PlatINV[2][1] * vecY) + (params::PlatINV[2][2] * vecZ));

            //check that the exchange is an integer value of unit cells. If not, there's likely an issue with the lattice vectors


            near[0] = nearbyint(NxP[i]);
            near[1] = nearbyint(NyP[i]);
            near[2] = nearbyint(NzP[i]);
            close_to_int[0] = std::abs(near[0]-NxP[i]);
            close_to_int[1] = std::abs(near[1]-NyP[i]);
            close_to_int[2] = std::abs(near[2]-NzP[i]);

            if ((close_to_int[0] > 0.001) || (close_to_int[1] > 0.001) || (close_to_int[2] > 0.001)){
                std::cout << "ERROR:  Unable to map Jij vector to unit cell." << std::endl;
				std::cout << "Line in Jij file: " << i+1 << std::endl;
				std::cout << "Vectors: " << vecX << " " << vecY << " " << vecZ << std::endl;
				std::cout << "Jij unit cell position: " << NxP[i] << " " << NyP[i] << " " << NzP[i] << std::endl;
                std::cout << "Jij values: " << Jijx[i] << " " << Jijy[i] << " " << Jijz[i] << std::endl;
                std::cout << "Exiting. ";
                std::cout << std::endl;
                exit(0);
            }

            NxP[i] = nearbyint(NxP[i]);
            NyP[i] = nearbyint(NyP[i]);
            NzP[i] = nearbyint(NzP[i]);

        }
        
	    int xval;
        int yval;
        int zval;
        int qval;
        int adjcounter = 0;
        int defectcounter = 0;

        bool xm,xp,ym,yp,zm,zp,bound;
        

        x_adj.push_back(0);

        // ========== Neighbour List =========== //
        if (params::Jij_units == "mRy") {
            for (int x = 0; x < params::Lx; ++x){                // Depth
                for (int y = 0; y < params::Ly; ++y){            // Row
                    for (int z = 0; z < params::Lz; ++z){        // Column
                        for (int q = 0; q < params::Nq; ++q){    // Unit Cell

                            // IRREGULAR SHAPE FOR HEDGEHOG SIMULATIONS - NEED A BETTER LONG TERM FIX FOR THIE
                            // if ((x != 0) && (x != params::Lx-1)){
                            //     xp = spins::sx1d(geom::LatCount(x+1,y  ,z  ,q)) == 0.0;
                            //     xm = spins::sx1d(geom::LatCount(x-1,y  ,z  ,q)) == 0.0;
                            // }
                            // if ((y != 0) && (y != params::Ly-1)){
                            //     yp = spins::sx1d(geom::LatCount(x  ,y+1,z  ,q)) == 0.0;
                            //     ym = spins::sx1d(geom::LatCount(x  ,y-1,z  ,q)) == 0.0; 
                            // }      
                            // if ((z != 0) && (z != params::Lz-1)){
                            //     zp = spins::sx1d(geom::LatCount(x  ,y  ,z+1,q)) == 0.0;
                            //     zm = spins::sx1d(geom::LatCount(x  ,y  ,z-1,q)) == 0.0;
                            // }
                            // bound = (xp || xm || yp || ym || zp || zm) == 1.0;

                            if ((params::xbound == "fixed") && ((x == 0) || (x == params::Lx-1)));
                            else if ((params::ybound == "fixed") && ((y == 0) || (y == params::Ly-1)));
                            else if ((params::zbound == "fixed") && ((z == 0) || (z == params::Lz-1)));
                            else if ((params::xbound == "fixed") && (params::ybound == "fixed") && (params::zbound == "fixed") && (bound == 1.0));
                            else if ((defects::list.size() != 0) && (geom::LatCount(x,y,z,q) == defects::list[defectcounter])){
                                defectcounter++;
                            }
                            else {
                                simspin.push_back(geom::LatCount(x,y,z,q));

                                for (int i = 0; i < length; i++){

                                    if (q == ib[i]+params::ibtoq){

                                        if (params::xbound == "periodic"){
                                            xval = modfunc(params::Lx, static_cast<int>(NxP[i]) + x);
                                        }
                                        else if ((params::xbound == "open") || (params::xbound == "fixed")) {
                                            xval = static_cast<int>(NxP[i]) + x;
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown xbound " << std::endl;
                                            exit(0);
                                        }
                                        if (params::ybound == "periodic"){
                                            yval = modfunc(params::Ly, static_cast<int>(NyP[i]) + y);
                                        }
                                        else if ((params::ybound == "open") || (params::ybound == "fixed")){
                                            yval = static_cast<int>(NyP[i]) + y; 
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown ybound " << std::endl;
                                            exit(0);
                                        }
                                        if (params::zbound == "periodic"){
                                            zval = modfunc(params::Lz, static_cast<int>(NzP[i]) + z);
                                        }
                                        else if ((params::zbound == "open") || (params::zbound == "fixed")){
                                            zval = static_cast<int>(NzP[i]) + z; 
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown zbound " << std::endl;
                                            exit(0);
                                        }

                                        qval = jb[i]+params::ibtoq;

                                        if (((xval >= 0) && (yval >= 0) && (zval >= 0)) && ((xval < params::Lx) && (yval < params::Ly) && (zval < params::Lz))) {
                                            adjncy.push_back(geom::LatCount(xval, yval, zval, qval));
                                            adjcounter++;


                                            jind.push_back(i);
                                            Jijx_prime.push_back( Jijx[i] /** 2.179872e-21) / params::mu_s*/);
                                            Jijy_prime.push_back( Jijy[i] /** 2.179872e-21) / params::mu_s*/);
                                            Jijz_prime.push_back( Jijz[i] /** 2.179872e-21) / params::mu_s*/);
                                        }
                                    }
                                }
                                x_adj.push_back(adjcounter);
                            }
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

                            // IRREGULAR SHAPE FOR HEDGEHOG SIMULATIONS - NEED A BETTER LONG TERM FIX FOR THIS
                            if ((x != 0) && (x != params::Lx-1)){
                                xp = spins::sx1d(geom::LatCount(x+1,y  ,z  ,q)) == 0.0;
                                xm = spins::sx1d(geom::LatCount(x-1,y  ,z  ,q)) == 0.0;
                            }
                            if ((y != 0) && (y != params::Ly-1)){
                                yp = spins::sx1d(geom::LatCount(x  ,y+1,z  ,q)) == 0.0;
                                ym = spins::sx1d(geom::LatCount(x  ,y-1,z  ,q)) == 0.0; 
                            }
                            if ((z != 0) && (z != params::Lz-1)){
                                zp = spins::sx1d(geom::LatCount(x  ,y  ,z+1,q)) == 0.0;
                                zm = spins::sx1d(geom::LatCount(x  ,y  ,z-1,q)) == 0.0;
                            }
                            bound = (xp || xm || yp || ym || zp || zm) == 1.0;

                            // std::cout << x << " " << y << " " << z << " " << q << "\n";


                            if ((params::xbound == "fixed") && ((x == 0) || (x == params::Lx-1)));
                            else if ((params::ybound == "fixed") && ((y == 0) || (y == params::Ly-1)));
                            else if ((params::zbound == "fixed")&& ((z == 0) || (z == params::Lz-1)));
                            else if ((params::xbound == "fixed") && (params::ybound == "fixed") && (params::zbound == "fixed") && (bound == 1.0));
                            else if ((defects::list.size() != 0) && (geom::LatCount(x,y,z,q) == defects::list[defectcounter])){
                                defectcounter++;
                            }
                            else {      
                                simspin.push_back(geom::LatCount(x,y,z,q));

                                for (int i = 0; i < length; i++){

                                    if (q == ib[i]+params::ibtoq){

                                        if (params::xbound == "periodic"){
                                            xval = modfunc(params::Lx, static_cast<int>(NxP[i]) + x);
                                        }
                                        else if ((params::xbound == "open") || (params::xbound == "fixed")) {
                                            xval = static_cast<int>(NxP[i]) + x;
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown xbound " << std::endl;
                                            exit(0);
                                        }
                                        if (params::ybound == "periodic"){
                                            yval = modfunc(params::Ly, static_cast<int>(NyP[i]) + y);
                                        }
                                        else if ((params::ybound == "open") || (params::ybound == "fixed")){
                                            yval = static_cast<int>(NyP[i]) + y; 
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown ybound " << std::endl;
                                            exit(0);
                                        }
                                        if (params::zbound == "periodic"){
                                            zval = modfunc(params::Lz, static_cast<int>(NzP[i]) + z);
                                        }
                                        else if ((params::zbound == "open") || (params::zbound == "fixed")){
                                            zval = static_cast<int>(NzP[i]) + z; 
                                        }
                                        else {
                                            std::cout << "ERROR: Unknown zbound " << std::endl;
                                            exit(0);
                                        }

                                        qval = jb[i]+params::ibtoq;

                                        if (((xval >= 0) && (yval >= 0) && (zval >= 0)) && ((xval < params::Lx) && (yval < params::Ly) && (zval < params::Lz))) {
                                            adjncy.push_back(geom::LatCount(xval, yval, zval, qval));
                                            adjcounter++;

                                            jind.push_back(i);
                                            Jijx_prime.push_back( (Jijx[i]) /*/ params::mu_s*/);
                                            Jijy_prime.push_back( (Jijy[i]) /*/ params::mu_s*/);
                                            Jijz_prime.push_back( (Jijz[i]) /*/ params::mu_s*/);
                                        }
                                    }
                                }
                                x_adj.push_back(adjcounter);
                            }
                        }
                    }
                }
            }  
        }
        else {
            std::cout << "ERROR: Unknown Jij units" << std::endl;
            exit(0);
        }

        // number of spins being simulated
        nsimspin = simspin.size();

        double Jijsize = 0;
        for (int i = 0; i < adjncy.size() / (x_adj.size()-1); i++){
            Jijsize += (Jijx[i]);
        }

		// Check x_adj is the same length as the number of spins
		if (x_adj.size() != params::Nspins+1 && ((params::xbound != "fixed") || (params::xbound != "fixed") || (params::xbound != "fixed"))){
			std::cout << "ERROR encountered at file " << __FILE__ << " at line " << __LINE__ << ": x_adj is not the same size as the number of atoms. Exiting." << std::endl;
			std::cout << "x_adj size: " << x_adj.size() << "\n";
		  	std::cout << "number of atoms: " << params::Nspins << std::endl;	
			exit(0);
		}

        INFO_OUT("length of x_adj:", x_adj.size());
        INFO_OUT("length of adjncy:", adjncy.size());
        INFO_OUT("average number of neighbours:", adjncy.size() / (x_adj.size()-1));
        INFO_OUT("length of Jij:", Jijz_prime.size());
        INFO_OUT("sum of neighbouring Jijs:", Jijsize << " (J)");
    }

    

}
