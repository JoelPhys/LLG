#include <cmath>
#include <fstream>
#include "../inc/util.h"
#include "../inc/params1.h"
#include "../inc/NeighbourList.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/mathfuncs.h"



namespace util {


	Array2D<double> M;
	Array<double> Mmag;
	Array<double> MdivMs;
	Array<double> MdivMsSum;
	Array2D<double> Msum;
	Array2D<double> MsumSQR;
	int isum = 0;

	std::ofstream magfile;

	void InitUtil(){
		M.resize(params::Nsublat,3);
		Mmag.resize(params::Nsublat);
		MdivMs.resize(params::Nsublat);
		MdivMsSum.resize(params::Nsublat);
		Msum.resize(params::Nsublat,3);
		MsumSQR.resize(params::Nsublat,3);
	}

	void InitOutputFile(double temp){
		std::stringstream sstr;
		sstr << "/home/b6033256/Mn2Au/output/mag_tsteps_" << params::Nt << "_T_" << temp << ".txt";
		magfile.open(sstr.str());
	}

	void ResetMag(){
		for (int ii = 0; ii < params::Nsublat; ii++){ 
			M(ii,0) = 0;
			M(ii,1) = 0;
			M(ii,2) = 0;
		}
	}

	//sort sites into sublattice
	void SortSublat(){

		for (int a = 0; a < params::Nspins; a++){   
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

		// output Neel Vector
		for (int m = 0; m < 3; m++){
			std::cout  <<(M(0,m) / params::NmomentsSubLat) + (M(1,m) / params::NmomentsSubLat) << "\t";
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
		for (int m = 0; m < 3; m++){
			magfile << (M(0,m) / params::NmomentsSubLat) + (M(1,m) / params::NmomentsSubLat) << "\t";
		}
		magfile << "\n";
	}

	void CloseMagFile(){
		magfile << std::flush;
		magfile.close();
	}


}
