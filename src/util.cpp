#include <cmath>
#include <fstream>
#include "../inc/util.h"
#include "../inc/params1.h"
#include "../inc/NeighbourList.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/mathfuncs.h"
#include "../inc/geom.h"



namespace util {


	Array2D<double> M;
	Array<double> Mt;
	Array<double> Mmag;
	Array<double> MdivMs;
	Array<double> MdivMsSum;
	Array2D<double> Msum;
	Array2D<double> MsumSQR;
	Array<double> sumx, sumy, sumz;
	int lc;

	int isum = 0;

	std::ofstream magfile;
	std::ofstream dwfile;

	void InitUtil(){
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
		sstr << params::filepath << "mag_tsteps_" << params::Nt << "_T_" << temp << ".dat";
		magfile.open(sstr.str());
	}

	void InitDWFile(double temp){
		std::stringstream sstr;
		sstr << params::filepath << "dw_T_" << temp << ".dat";
		dwfile.open(sstr.str());
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


		if (params::afmflag == "SC"){
			for (int a = 0; a < params::Nspins; a++){   

				// Total Magnetisation
				Mt(0) += neigh::Sx1d(a);			
				Mt(1) += neigh::Sy1d(a);
				Mt(2) += neigh::Sz1d(a);

				if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 3) || (modfunc(params::Nq,a) == 5) || (modfunc(params::Nq,a) == 6)) {
					M(0,0) += neigh::Sx1d(a);
					M(0,1) += neigh::Sy1d(a);
					M(0,2) += neigh::Sz1d(a); 
				}
				else {
					M(1,0) += neigh::Sx1d(a);
					M(1,1) += neigh::Sy1d(a);
					M(1,2) += neigh::Sz1d(a); 
				}
			}
		}
		else if (params::afmflag == "N"){
			for (int a = 0; a < params::Nspins; a++){   
				M(0,0) += neigh::Sx1d(a);
				M(0,1) += neigh::Sy1d(a);
				M(0,2) += neigh::Sz1d(a); 
			}
		}
		else if (params::afmflag == "Mn2Au"){
			for (int a = 0; a < params::Nspins; a++){   

				// Total Magnetisation
				Mt(0) += neigh::Sx1d(a);			
				Mt(1) += neigh::Sy1d(a);
				Mt(2) += neigh::Sz1d(a);

				if (modfunc(params::Nq,a) == 0){
					M(0,0) += neigh::Sx1d(a);
					M(0,1) += neigh::Sy1d(a);
					M(0,2) += neigh::Sz1d(a); 
				}
				else if (modfunc(params::Nq,a) == 1) {
					M(1,0) += neigh::Sx1d(a);
					M(1,1) += neigh::Sy1d(a);
					M(1,2) += neigh::Sz1d(a); 
				}
				else if (modfunc(params::Nq,a) == 2) {
					M(0,0) += neigh::Sx1d(a);
					M(0,1) += neigh::Sy1d(a);
					M(0,2) += neigh::Sz1d(a);
				}
				else if (modfunc(params::Nq,a) == 3) {
					M(1,0) += neigh::Sx1d(a);
					M(1,1) += neigh::Sy1d(a);
					M(1,2) += neigh::Sz1d(a); 
				}
				else {
					std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
					exit(0);
				}
			}
		}
		else if (params::afmflag == "NiO"){
			for (int a = 0; a < params::Nspins; a++){   

				// Total Magnetisation
				Mt(0) += neigh::Sx1d(a);			
				Mt(1) += neigh::Sy1d(a);
				Mt(2) += neigh::Sz1d(a);

				if (a / params::Nq % 2 == 0){
					if (modfunc(params::Nq,a) == 0){
						M(1,0) += neigh::Sx1d(a);
						M(1,1) += neigh::Sy1d(a);
						M(1,2) += neigh::Sz1d(a); 
					}
					else if (modfunc(params::Nq,a) == 1) {
						M(1,0) += neigh::Sx1d(a);
						M(1,1) += neigh::Sy1d(a);
						M(1,2) += neigh::Sz1d(a); 
					}
					else if (modfunc(params::Nq,a) == 2) {
						M(0,0) += neigh::Sx1d(a);
						M(0,1) += neigh::Sy1d(a);
						M(0,2) += neigh::Sz1d(a);
					}
					else if (modfunc(params::Nq,a) == 3) {
						M(1,0) += neigh::Sx1d(a);
						M(1,1) += neigh::Sy1d(a);
						M(1,2) += neigh::Sz1d(a); 
					}
					else {
						std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
						exit(0);
					}
				}
				else if (a / params::Nq % 2 == 1) {
					if (modfunc(params::Nq,a) == 0){
						M(0,0) += neigh::Sx1d(a);
						M(0,1) += neigh::Sy1d(a);
						M(0,2) += neigh::Sz1d(a); 
					}
					else if (modfunc(params::Nq,a) == 1) {
						M(0,0) += neigh::Sx1d(a);
						M(0,1) += neigh::Sy1d(a);
						M(0,2) += neigh::Sz1d(a); 
					}
					else if (modfunc(params::Nq,a) == 2) {
						M(1,0) += neigh::Sx1d(a);
						M(1,1) += neigh::Sy1d(a);
						M(1,2) += neigh::Sz1d(a);
					}
					else if (modfunc(params::Nq,a) == 3) {
						M(0,0) += neigh::Sx1d(a);
						M(0,1) += neigh::Sy1d(a);
						M(0,2) += neigh::Sz1d(a); 
					}
					else {
						std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
					}	
				}
			}
		}
		else {
			std::cerr << "ERROR: Unassigned afmflag" << std::endl;
			exit(0);
		}
	}

	void MagLength(){

		for (int l = 0; l < params::Nsublat; l++){
			Mmag(l) = sqrt(M(l,0) * M(l,0) + M(l,1) * M(l,1) + M(l,2) * M(l,2));
			MdivMs(l) = Mmag(l) / (params::NmomentsSubLat);
		}
	}

	void SumMag(int i){
		if (i > (params::Nt / 10) - 1) {
			for (int l = 0; l < params::Nsublat; l++){
				for (int m = 0; m < 3; m++){
					Msum(l,m) += M(l,m) / params::NmomentsSubLat;
					MsumSQR(l,m) += (M(l,m) / params::NmomentsSubLat) * (M(l,m) / params::NmomentsSubLat);
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
		std::cout << std::fixed << i <<  " | ";
		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				std::cout << std::fixed << std::setprecision(6) << M(l,m) / params::NmomentsSubLat << "\t"; 
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
		std::cout << "\n";
	}

	void OutputMagToFile(int i){

		magfile << i << " ";

		//Output sublattice magnetisation
		for (int l = 0; l < params::Nsublat; l++){
			for (int m = 0; m < 3; m++){
				magfile << M(l,m) / params::NmomentsSubLat << "\t"; 
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
		magfile << "\n";
	}

	void CloseMagFile(){
		magfile << std::flush;
		magfile.close();
	}

	void OutputDWtoFile(int i){
		for (int h = 0; h < params::Lx; h++){
			sumx(h) = 0.0;
			sumy(h) = 0.0;
			sumz(h) = 0.0;
			for (int j =0; j < params::Ly; j++){
				for (int k =0; k < params::Lz; k++){
					lc = geom::LatCount(h,j,k,2);
					sumx(h) += neigh::Sx1d(lc);
					sumy(h) += neigh::Sy1d(lc);
					sumz(h) += neigh::Sz1d(lc);
				}
			}
			sumx(h) /= params::Ly*params::Lz;
			sumy(h) /= params::Ly*params::Lz;
			sumz(h) /= params::Ly*params::Lz;
			dwfile << sumx(h) << " " << sumy(h) << " " << sumz(h) << " ";
		}
		dwfile << "\n";
	}



}
