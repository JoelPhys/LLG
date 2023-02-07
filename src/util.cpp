// cpp header files
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

// my header files
#include "../inc/util.h"
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/array.h"
#include "../inc/fields.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/array2d.h"
#include "../inc/mathfuncs.h"
#include "../inc/neighbourlist.h"


namespace util {


	double E;
	Array<double> Mt;
	Array<double> Mmag;
	Array<double> MdivMs;
	Array<double> MdivMsSum;
	Array<double> sumx, sumy, sumz;

	Array2D<double> M;
	Array2D<double> Msum;
	Array2D<double> MsumSQR;

	int lc;

	int isum = 0;
	int sublatindex;

	clock_t begin, end;

	std::ofstream magfile;
	std::ofstream dwfile;
	std::ofstream fldfile;
	std::ofstream latfile;
	int lati;

	// for averaging over lattice
	double avgx,avgy,avgz;

	void init(){
		Mt.resize(3);
		M.resize(params::Nsublat,3);
		Mmag.resize(params::Nsublat);
		MdivMs.resize(params::Nsublat);
		MdivMsSum.resize(params::Nsublat);
		Msum.resize(params::Nsublat,3);
		MsumSQR.resize(params::Nsublat,3);
		sumx.resize(params::Lx);
		sumy.resize(params::Lx);
		sumz.resize(params::Lx);
	}

	void InitMagFile(double temp){
		std::stringstream sstr;
		// sstr << params::filepath << "mag_tsteps_" << params::Nt << "_T_" << temp << ".dat";
		sstr << params::filepath << "mag_tsteps_" << params::Nt << "_T_" << std::setw(4) << std::setfill('0') << temp << ".out";
		magfile.open(sstr.str());
	}

	void InitDWFile(double temp){
		std::stringstream sstr;
		sstr << params::filepath << "dw_T_" << std::setw(4) << std::setfill('0') << temp << ".out";
		dwfile.open(sstr.str());
	}

	void InitFldFile(double temp){
		std::stringstream sstr;
		sstr << params::filepath << "field_" << std::setw(4) << std::setfill('0') << temp << ".out";
		fldfile.open(sstr.str());
	}

	void ResetMag(){

		// Total magnetisation
		Mt(0) = 0.0;
		Mt(1) = 0.0;
		Mt(2) = 0.0;

		// Sublattice magnetisation
		for (int ii = 0; ii < params::Nsublat; ii++){ 
			M(ii,0) = 0;
			M(ii,1) = 0;
			M(ii,2) = 0;
		}
	}

	//sort sites into sublattice
	void SortSublat(){

		E = 0.0;
		
		for (int a = 0; a < params::Nspins; a++){
			
			// Total Magnetisation
			Mt(0) += spins::sx1d(a);			
			Mt(1) += spins::sy1d(a);
			Mt(2) += spins::sz1d(a);

			sublatindex = params::sublat_sites[a % params::Nq];
			M(sublatindex,0) += spins::sx1d(a);
			M(sublatindex,1) += spins::sy1d(a);
			M(sublatindex,2) += spins::sz1d(a);

			E += spins::Ex(a);

		}

		// if (params::afmflag == "SC"){
		// 	for (int a = 0; a < params::Nspins; a++){   

		// 		// Total Magnetisation
		// 		Mt(0) += spins::sx1d(a);			
		// 		Mt(1) += spins::sy1d(a);
		// 		Mt(2) += spins::sz1d(a);

		// 		if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 3) || (modfunc(params::Nq,a) == 5) || (modfunc(params::Nq,a) == 6)) {
		// 			M(0,0) += spins::sx1d(a);
		// 			M(0,1) += spins::sy1d(a);
		// 			M(0,2) += spins::sz1d(a); 
		// 		}
		// 		else {
		// 			M(1,0) += spins::sx1d(a);
		// 			M(1,1) += spins::sy1d(a);
		// 			M(1,2) += spins::sz1d(a); 
		// 		}
		// 	}
		// }
		// else if (params::afmflag == "N"){
		// 	for (int a = 0; a < params::Nspins; a++){   
		// 		M(0,0) += spins::sx1d(a);
		// 		M(0,1) += spins::sy1d(a);
		// 		M(0,2) += spins::sz1d(a); 
		// 	}
		// }
		// else if (params::afmflag == "CuMnAs"){
		// 	for (int a = 0; a < params::Nspins; a++){   

		// 		// Total Magnetisation
		// 		Mt(0) += spins::sx1d(a);			
		// 		Mt(1) += spins::sy1d(a);
		// 		Mt(2) += spins::sz1d(a);

		// 		if (modfunc(params::Nq,a) == 0){
		// 			M(0,0) += spins::sx1d(a);
		// 			M(0,1) += spins::sy1d(a);
		// 			M(0,2) += spins::sz1d(a); 
		// 		}
		// 		else if (modfunc(params::Nq,a) == 1) {
		// 			M(1,0) += spins::sx1d(a);
		// 			M(1,1) += spins::sy1d(a);
		// 			M(1,2) += spins::sz1d(a); 
		// 		}
		// 		else {
		// 			std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
		// 			exit(0);
		// 		}
		// 	}
		// }
		// else if (params::afmflag == "Mn2Au"){
		// 	for (int a = 0; a < params::Nspins; a++){   

		// 		// Total Magnetisation
		// 		Mt(0) += spins::sx1d(a);			
		// 		Mt(1) += spins::sy1d(a);
		// 		Mt(2) += spins::sz1d(a);

		// 		if (modfunc(params::Nq,a) == 0){
		// 			M(0,0) += spins::sx1d(a);
		// 			M(0,1) += spins::sy1d(a);
		// 			M(0,2) += spins::sz1d(a); 
		// 		}
		// 		else if (modfunc(params::Nq,a) == 1) {
		// 			M(1,0) += spins::sx1d(a);
		// 			M(1,1) += spins::sy1d(a);
		// 			M(1,2) += spins::sz1d(a); 
		// 		}
		// 		else if (modfunc(params::Nq,a) == 2) {
		// 			M(0,0) += spins::sx1d(a);
		// 			M(0,1) += spins::sy1d(a);
		// 			M(0,2) += spins::sz1d(a);
		// 		}
		// 		else if (modfunc(params::Nq,a) == 3) {
		// 			M(1,0) += spins::sx1d(a);
		// 			M(1,1) += spins::sy1d(a);
		// 			M(1,2) += spins::sz1d(a); 
		// 		}
		// 		else {
		// 			std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
		// 			exit(0);
		// 		}
		// 	}
		// }
		// else if (params::afmflag == "NiO"){
		// 	for (int a = 0; a < params::Nspins; a++){   

		// 		// Total Magnetisation
		// 		Mt(0) += spins::sx1d(a);			
		// 		Mt(1) += spins::sy1d(a);
		// 		Mt(2) += spins::sz1d(a);

		// 		if (a / params::Nq % 2 == 0){
		// 			if (modfunc(params::Nq,a) == 0){
		// 				M(1,0) += spins::sx1d(a);
		// 				M(1,1) += spins::sy1d(a);
		// 				M(1,2) += spins::sz1d(a); 
		// 			}
		// 			else if (modfunc(params::Nq,a) == 1) {
		// 				M(1,0) += spins::sx1d(a);
		// 				M(1,1) += spins::sy1d(a);
		// 				M(1,2) += spins::sz1d(a); 
		// 			}
		// 			else if (modfunc(params::Nq,a) == 2) {
		// 				M(0,0) += spins::sx1d(a);
		// 				M(0,1) += spins::sy1d(a);
		// 				M(0,2) += spins::sz1d(a);
		// 			}
		// 			else if (modfunc(params::Nq,a) == 3) {
		// 				M(1,0) += spins::sx1d(a);
		// 				M(1,1) += spins::sy1d(a);
		// 				M(1,2) += spins::sz1d(a); 
		// 			}
		// 			else {
		// 				std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
		// 				exit(0);
		// 			}
		// 		}
		// 		else if (a / params::Nq % 2 == 1) {
		// 			if (modfunc(params::Nq,a) == 0){
		// 				M(0,0) += spins::sx1d(a);
		// 				M(0,1) += spins::sy1d(a);
		// 				M(0,2) += spins::sz1d(a); 
		// 			}
		// 			else if (modfunc(params::Nq,a) == 1) {
		// 				M(0,0) += spins::sx1d(a);
		// 				M(0,1) += spins::sy1d(a);
		// 				M(0,2) += spins::sz1d(a); 
		// 			}
		// 			else if (modfunc(params::Nq,a) == 2) {
		// 				M(1,0) += spins::sx1d(a);
		// 				M(1,1) += spins::sy1d(a);
		// 				M(1,2) += spins::sz1d(a);
		// 			}
		// 			else if (modfunc(params::Nq,a) == 3) {
		// 				M(0,0) += spins::sx1d(a);
		// 				M(0,1) += spins::sy1d(a);
		// 				M(0,2) += spins::sz1d(a); 
		// 			}
		// 			else {
		// 				std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
		// 			}	
		// 		}
		// 	}
		// }
		// else {
		// 	std::cout << "ERROR: Unassigned afmflag" << std::endl;
		// 	exit(0);
		// }
	}

	void MagLength(){

		for (int l = 0; l < params::Nsublat; l++){
			Mmag(l) = sqrt(M(l,0) * M(l,0) + M(l,1) * M(l,1) + M(l,2) * M(l,2));
			MdivMs(l) = Mmag(l) / (params::NmomentsSubLat[l]);
		}
	}

	void SumMag(int i){
		if (i > (params::Nt / 10) - 1) {
			for (int l = 0; l < params::Nsublat; l++){
				for (int m = 0; m < 3; m++){
					Msum(l,m) += M(l,m) / params::NmomentsSubLat[l];
					MsumSQR(l,m) += (M(l,m) / params::NmomentsSubLat[l]) * (M(l,m) / params::NmomentsSubLat[l]);
				}
				MdivMsSum[l] += MdivMs[l];
			}
			isum++;
		}
	}

	void OutputSumMag(){
		std::cout << "For averaging: " << std::endl;

		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				std::cout << Msum(l,m) <<  "\t";
			}
		}

		std::cout << MdivMsSum[0] << "\t" << MdivMsSum[1];
		std::cout << "\n";

		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				std::cout << MsumSQR(l,m) <<  "\t";
			}
		}
		std::cout << "\n";
		std::cout << isum << std::endl;
	}

	void OutputMagToTerm(int i){
		std::cout << std::scientific << i << " " << i * params::dt  <<  " | ";
		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				std::cout << std::fixed << std::setprecision(6) << M(l,m) / params::NmomentsSubLat[l] << "\t"; 
			}
			std::cout << std::fixed << std::setprecision(6) << MdivMs(l) << "\t | \t";
		}

		//output Neel Vector
		// for (int m = 0; m < 3; m++){
		// 	std::cout  <<(M(0,m) / params::NmomentsSubLat) + (M(1,m) / params::NmomentsSubLat) << "\t";
		// }

		// output Total magnetisation
		if (params::afmflag != "N"){
			std::cout << Mt(0) / params::Nspins << "\t" << Mt(1)  / params::Nspins<< "\t" << Mt(2) / params::Nspins;
		}

		//output total energy
		//std::cout << "\t | \t" << E;
		
		std::cout << "\n";
	}

	void OutputMagToFile(int i){

		magfile << i * params::dt << " ";

		//Output sublattice magnetisation
		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				magfile << M(l,m) / params::NmomentsSubLat[l] << "\t"; 
			}
			magfile << MdivMs(l) << "\t";
		}

		//output Neel Vector
		// for (int m = 0; m < 3; m++){
		// 	magfile << (M(0,m) / params::NmomentsSubLat) + (M(1,m) / params::NmomentsSubLat) << "\t";
		// }

		// output Total magnetisation
		if (params::afmflag != "N"){
			magfile << Mt(0)  / params::Nspins << "\t" << Mt(1)  / params::Nspins<< "\t" << Mt(2) / params::Nspins;
		}

		//output total energy
		magfile << "\t" << E;
		
		magfile << "\n";
	}

	void CloseMagFile(){
		magfile << std::flush;
		fldfile << std::flush;
		magfile.close();
		fldfile.close();
	}
	

	void OutputLatticetoFile(double temp){

		std::stringstream sstr;
		sstr << params::OutputLatticeFilepath << "temp_" << std::setw(4) << std::setfill('0') << temp << "_file_" << std::setw(4) << std::setfill('0') << lati << ".lat";
		latfile.open(sstr.str());
		
		if (params::OutputLatticeAverageOver == "false"){
			for (int i = 0; i < params::Lx; i++){
				for (int j = 0; j < params::Ly; j++){
					for (int k = 0; k < params::Lz; k++){
						for (int q = 0; q < params::Nq; q++){

						
							

							latfile << i << " " << j << " " << k << " " << q << " ";
							//latfile << geom::latticeX(i,j,k,q) << " ";
							//latfile << geom::latticeY(i,j,k,q) << " ";
							//latfile << geom::latticeZ(i,j,k,q) << " ";
							latfile << spins::sx1d(geom::LatCount(i,j,k,q)) << " ";
							latfile << spins::sy1d(geom::LatCount(i,j,k,q)) << " ";
							latfile << spins::sz1d(geom::LatCount(i,j,k,q)) << "\n";
						
						
						}
					}
				}
			}
		}
		else if (params::OutputLatticeAverageOver == "x"){
			OutputLatticeAverageOverX();
		}	
		else if (params::OutputLatticeAverageOver == "y"){
			OutputLatticeAverageOverY();
		}	
		else if (params::OutputLatticeAverageOver == "z"){
			OutputLatticeAverageOverZ();
		}
		else if (params::OutputLatticeAverageOver == "q"){
			OutputLatticeAverageOverQ();
		}
        else {
            std::cout << "ERROR: Unknown Util.OutputLatticeAverageOver specified in config file. \n";
            std::cout << "Util.OutputLatticeAverageOver = " << params::OutputLatticeAverageOver << std::endl;
            std::cout << "Exiting." << std::endl;
            exit(0);
        }
	
		// Reset averaging over components
		avgx = 0.0;	avgy = 0.0;	avgz = 0.0;

		lati++;
		latfile << std::flush;
		latfile.close(); 

	}

	void OutputDWtoFile(int i){
		for (int h = 0; h < params::Lx; h++){
			sumx(h) = 0.0;
			sumy(h) = 0.0;
			sumz(h) = 0.0;
			for (int j =0; j < params::Ly; j++){
				for (int k =0; k < params::Lz; k++){
					lc = geom::LatCount(h,j,k,0);
					sumx(h) += spins::sx1d(lc);
					sumy(h) += spins::sy1d(lc);
					sumz(h) += spins::sz1d(lc);
				}
			}
			sumx(h) /= params::Ly*params::Lz;
			sumy(h) /= params::Ly*params::Lz;
			sumz(h) /= params::Ly*params::Lz;
			dwfile << sumx(h) << " " << sumy(h) << " " << sumz(h) << " ";
		}
		dwfile << "\n";
	}

	void OutputFldToFile(int i){
	
		// output time
		fldfile << static_cast<double>(i) * params::dt << " ";
		
		// Output Field for each sublattice
		for (int i = 0; i < params::uniquesublat.size(); i++){
			fldfile << fields::H_appx[params::uniquesublat[i]] << " ";
			fldfile << fields::H_appy[params::uniquesublat[i]] << " ";
			fldfile << fields::H_appz[params::uniquesublat[i]] << " ";
		}

		fldfile << "\n";
	}

	void startclock(){
		double time_spent;
		begin = clock();
	}

	void endclock(){
		end = clock();
		double endtime = (double)(end - begin) / CLOCKS_PER_SEC;
		INFO_OUT("Simulation Time: ", std::setprecision(10) << endtime << " seconds"); 
	}

	void readexternalspins(std::string cmdline){

		// spicify filename
		std::stringstream sstr_eq;
		sstr_eq << cmdline;
		// std::cout << sstr_eq.str() << std::endl;
		std::ifstream equilibrationfile(sstr_eq.str());

		INFO_OUT("External Spin file:", cmdline);

		// Check if equilibrium file could be opened
		if (!equilibrationfile){
			std::cout << "ERROR: Could not open equilibrium file: " << sstr_eq.str() << std::endl;
			exit(0);
		}

		double sx, sy, sz;
		int posx, posy, posz, posq;
		int count = 0;

		// Check number of cols in file to see if input row is per atom OR per unit cell
		std::string line;
		getline(equilibrationfile,line);
		std::stringstream s;
		s << line;
		int numcols = 0;
		double value;
		while(s >> value) numcols++;

		// reset file to first line
		equilibrationfile.clear();
		equilibrationfile.seekg(0);

		// Loop through file
		if (numcols == 6){
			
			// output number of cols
			std::cout << "Equilibration file is per UNIT CELL" << std::endl;
			
			// loop through file
			while (equilibrationfile >> posx >> posy >> posz >> sx >> sy >> sz){
				
				for (int q = 0; q < params::Nq; q++){
					spins::sx1d(geom::LatCount(posx,posy,posz,q)) = sx;
					spins::sy1d(geom::LatCount(posx,posy,posz,q)) = sy;
					spins::sz1d(geom::LatCount(posx,posy,posz,q)) = sz;
				}
				count++;
			}
			// check number of lines is equal to number of unit cells
			if (count != params::Lx*params::Ly*params::Lz){
				std::cout << "ERROR: Wrong number of spins in file. Exiting." << std::endl;
				INFO_OUT("Number of lines in Equilibration File:", count);
				INFO_OUT("Number of unit cells:", params::Lx*params::Ly*params::Lz);
				exit(0);
			}	
		}
		else if (numcols == 7){

			// Output number of cols
			std::cout << "Equilibration file is per ATOM" << std::endl;
			
			// loop through file
			while (equilibrationfile >> posx >> posy >> posz >> posq >> sx >> sy >> sz){
				
				spins::sx1d(geom::LatCount(posx,posy,posz,posq)) = sx;
				spins::sy1d(geom::LatCount(posx,posy,posz,posq)) = sy;
				spins::sz1d(geom::LatCount(posx,posy,posz,posq)) = sz;
				count++;
			}
			
			// Check number of lines is equal to number of atoms
			if (count != params::Lx*params::Ly*params::Lz*params::Nq){
				std::cout << "ERROR: Wrong number of spins in file. Exiting." << std::endl;
				exit(0);
			}	
		}
		else {
			std::cout << "Unknown Number of cols in equilibration file." << std::endl;
			INFO_OUT("Number of cols:",numcols);
			exit(0);
		}

		
		

		// close file
		equilibrationfile.close();

	}


	void OutputLatticeAverageOverX(){

		for (int j = 0; j < params::Ly; j++){
			for (int k = 0; k < params::Lz; k++){
				for (int q = 0; q < params::Nq; q++){				
					for (int i = 0; i < params::Lx; i++){
						avgx +=	spins::sx1d(geom::LatCount(i,j,k,q));
						avgy +=	spins::sy1d(geom::LatCount(i,j,k,q));
						avgz +=	spins::sz1d(geom::LatCount(i,j,k,q));
					}
					
					// divive by number of layers along direction
					latfile << j << " " << k << " " << q << " ";
					latfile << avgx/params::Lx << " ";
					latfile << avgy/params::Lx << " ";
					latfile << avgz/params::Lx << "\n";
					
					// reset averages	
					avgx = 0.0; avgy = 0; avgz = 0;
				}
			}
		}
	}

	void OutputLatticeAverageOverY(){

		for (int i = 0; i < params::Lx; i++){
			for (int k = 0; k < params::Lz; k++){
				for (int q = 0; q < params::Nq; q++){
					for (int j = 0; j < params::Ly; j++){
							avgx +=	spins::sx1d(geom::LatCount(i,j,k,q));
							avgy +=	spins::sy1d(geom::LatCount(i,j,k,q));
							avgz +=	spins::sz1d(geom::LatCount(i,j,k,q));
					}

					// divive by number of layers along direction
					latfile << i << " " << k << " " << q << " ";
					latfile << avgx/params::Ly << " ";
					latfile << avgy/params::Ly << " ";
					latfile << avgz/params::Ly << "\n";

					// reset averages	
					avgx = 0.0; avgy = 0; avgz = 0;
				
				}
			}
		}
	}

	void OutputLatticeAverageOverZ(){

		for (int i = 0; i < params::Lx; i++){
			for (int j = 0; j < params::Ly; j++){
				for (int q = 0; q < params::Nq; q++){
					for (int k = 0; k < params::Lz; k++){
							avgx +=	spins::sx1d(geom::LatCount(i,j,k,q));
							avgy +=	spins::sy1d(geom::LatCount(i,j,k,q));
							avgz +=	spins::sz1d(geom::LatCount(i,j,k,q));
					}
					
					// divive by number of layers along direction
					latfile << i << " " << j << " " << q << " ";
					latfile << avgx/params::Lz << " ";
					latfile << avgy/params::Lz << " ";
					latfile << avgz/params::Lz << "\n";

					// reset averages	
					avgx = 0.0; avgy = 0; avgz = 0;
				}
			}
		}
	}
	
	void OutputLatticeAverageOverQ(){

		for (int i = 0; i < params::Lx; i++){
			for (int j = 0; j < params::Ly; j++){
				for (int k = 0; k < params::Lz; k++){
					for (int q = 0; q < params::Nq; q++){
							avgx +=	spins::sx1d(geom::LatCount(i,j,k,q));
							avgy +=	spins::sy1d(geom::LatCount(i,j,k,q));
							avgz +=	spins::sz1d(geom::LatCount(i,j,k,q));
					}
					latfile << i << " " << j << " " << k << " ";
					latfile << avgx/params::Nq << " ";
					latfile << avgy/params::Nq << " ";
					latfile << avgz/params::Nq << "\n";

					// reset averages	
					avgx = 0.0; avgy = 0; avgz = 0;
				}
			}
		}
	}
}
