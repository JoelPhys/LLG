// cpp header files
#include <cmath>
#include <vector>
#include <random>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>

// my header files
#include "../inc/geom.h"
#include "../inc/fftw3.h"
#include "../inc/spins.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/mathfuncs.h"
#include "../inc/spinwaves.h"
#include "../inc/libconfig.h++"
#include "../inc/neighbourlist.h"

namespace spinwaves {

	// Macros for FFT
	#define REAL 0
	#define IMAG 1

	// initialise arrays
	Array3D<double> stcf;
	Array3D<fftw_complex> FFTstcf;
	Array2D<fftw_complex> FFTstcfarray;
	Array<fftw_complex> FFTstcfT;
	Array<fftw_complex> stcfT;

	//Testing for Mn2Au
	Array4D<fftw_complex> FFTstcfarray_test;

	// FFTW plans
	fftw_plan Slat;
	fftw_plan Stime;

	std::ofstream file_spnwvs;

	int Npoints;
	int icount;

	// variables for Hamming Windowing function
	double windowing = 0;
	double hammA = 0.54;
	double hammB = 0.46;

	double norm;

	int lval, mval, nval;

	// Kpath variables from config giles
	std::vector<double> kpathx;
	std::vector<double> kpathy;
	std::vector<double> kpathz;
	int kpath_length;

	// testing with 1d
	Array<double> stcf1d;
	Array<fftw_complex> FFTstcf1d;
	fftw_plan Slat1d;

	void init(){

	
		// Read kpath from config file
		libconfig::Setting& setting = params::cfg.lookup("Spinwaves");

		// Check if kpaths are the same lengths
		if ((setting["kpathx"].getLength() == setting["kpathy"].getLength()) && (setting["kpathy"].getLength() == setting["kpathx"].getLength())){

			kpath_length = setting["kpathx"].getLength();
			
			// Add kpath from config file to cfg array
			for (int v = 0; v < kpath_length; v++){
				kpathx.push_back(setting["kpathx"][v]);
				kpathy.push_back(setting["kpathy"][v]);
				kpathz.push_back(setting["kpathz"][v]);
				v++;
			}
		}
		else {
			std::cout << "ERROR: kpaths are different lengths. Exiting." << std::endl;
			exit(0); 
		}

		// Print kpath to log 


		// calculate number of points in time array
		Npoints  = ceil((params::Nt - params::start) / (params::dt_spinwaves / params::dt));
		INFO_OUT("number of timepoints for spinwaves = ", Npoints);
		
		
		icount = 0;
		norm = 1 / (geom::Ix * geom::Iy * geom::Iz * Npoints);

		stcf.resize(geom::Ix,geom::Iy,geom::Iz);
		FFTstcf.resize(geom::Ix,geom::Iy,geom::IzC);
		stcfT.resize(Npoints);
		FFTstcfT.resize(Npoints);

		// testing for 1d
		stcf1d.resize(params::Lx);
		FFTstcf1d.resize(params::Lx/2+1);
		FFTstcfarray.resize(Npoints, params::Lx/2+1);

		//TESTING WITH 2D outputs
		FFTstcfarray_test.resize(Npoints, geom::Ix, geom::Iy, geom::IzC);

		// FFT plans
		fftw_set_timelimit(60);
		Slat = fftw_plan_dft_r2c_3d(geom::Ix,geom::Iy,geom::Iz, stcf.ptr(), FFTstcf.ptr(), FFTW_MEASURE);
		fftw_set_timelimit(60);
		Stime = fftw_plan_dft_1d(Npoints, FFTstcfT.ptr(), stcfT.ptr(), FFTW_FORWARD, FFTW_MEASURE);

		
		stcfT.IFill(0);
		FFTstcf.IFill(0);
		FFTstcfT.IFill(0);
		stcf.IFill(0);
		FFTstcfarray.IFill(0);
		
		//TESTING MN2AU
		fftw_set_timelimit(60);
		Slat1d = fftw_plan_dft_r2c_1d(params::Lx, stcf1d.ptr(), FFTstcf1d.ptr(), FFTW_MEASURE);
		stcf1d.IFill(0);
		FFTstcf1d.IFill(0);

		// intialise output file 
		std::stringstream spnwvs;
		spnwvs << "output/spinwaves/spinwaves.txt";
		file_spnwvs.open(spnwvs.str());
	}

	void FFTspace(){

		// // convert back to 3D array from 1D list 
		for (int l = 0; l < geom::Ix; l += params::Idx){
			for (int m = 0; m < geom::Iy; m += params::Idy){
				for (int n = 0; n < geom::Iz; n += params::Idz){
					for (int q = 0; q < params::Nq; q++){
						lval = l + params::Isites[q][0];
						mval = m + params::Isites[q][1];
						nval = n + params::Isites[q][2];
						stcf(lval,mval,nval) = spins::sx1d(geom::Scount(lval,mval,nval)) * spins::sx1d(geom::Scount(lval,mval,nval)) + spins::sz1d(geom::Scount(lval,mval,nval)) * spins::sz1d(geom::Scount(lval,mval,nval));
					}
				}
			}
		}

		// compute fft
		fftw_execute(Slat);

		// //TESTING WITH MN2AU DELETE IF NOT USING
		// for (int l = 0; l < params::Lx; l++){
		//     stcf1d(l) = spins::sx1d(geom::LatCount(l,1,1,0)) * spins::sx1d(geom::LatCount(l,1,1,0)) + spins::sz1d(geom::LatCount(l,1,1,0)) * spins::sz1d(geom::LatCount(l,1,1,0));
		// }

		// testing for Mn2Au
		// fftw_execute(Slat1d);

		//windowing function
		windowing = hammA - hammB * cos((2 * M_PI * icount) / (Npoints  - 1));     

		// for (int j = 0; j < params::Lx/2+1; j++){

		// //     if (j < kpts){
		// //         FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf(j,0,0)[REAL];
		// //         FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf(j,0,0)[IMAG];
		// //         file_spnwvs << FFTstcf(j,0,0)[REAL] << " " << FFTstcf(j,0,0)[IMAG] << "\t";

		//     // TESTING MN2AU
		//     FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf1d(j)[REAL];
		//     FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf1d(j)[IMAG];
		//     // file_spnwvs << FFTstcf1d(j)[REAL] << " " << FFTstcf1d(j)[IMAG] << "\t";
		// }

		for (int j = 0; j < geom::Ix; j++){
			for (int l = 0; l < geom::Iy; l++){
				for (int k = 0; k < geom::IzC; k++){

					FFTstcfarray_test(icount,j,l,k)[REAL] = windowing * FFTstcf(j,l,k)[REAL];
					FFTstcfarray_test(icount,j,l,k)[IMAG] = windowing * FFTstcf(j,l,k)[IMAG];
				}
			}
		}

		file_spnwvs << "\n";
		icount++;

	}

	void FFTtime(){

		file_spnwvs.close();

 		double j1, kpoint;
		double os1, os2, os[icount/2], freq;
		double freqstep = (1/params::dt_spinwaves) / icount; 

		// Create output file for peaks
		std::stringstream peakstring;
		peakstring << params::filepath_sw << "peaks.txt";
		std::ofstream peaksout;
		peaksout.open(peakstring.str());
		peaksout << std::setprecision(10);

		// loop over all elememts in k array
		double kpts = geom::Ix; 

		// double kpathx[7] = {0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.423940489};
		// double kpathy[7] = {0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0};
		// double kpathz[7] = {0.0, 0.0, 0.5, 0.5, 0.0, 1.0, 1.0};

		// double kpathx[13] = {0, 14, 14, 14,  0,  0, 11, 17,  0, 14, 17, 11,  0};
		// double kpathy[13] = {0, 14, 14,  0,  0,  0,  0,  0,  0, 14, 11, 11,  0};
		// double kpathz[13] = {0,  0, 14, 14,  0, 29, 30,  0,  0,  0,  0, 30, 30};

		// double kpathx[13] = {0, 15, 15, 15,  0,  0, 12, 18,  0, 15, 17, 12,  0};
		// double kpathy[13] = {0, 15, 15,  0,  0,  0,  0,  0,  0, 15, 11, 12,  0};
		// double kpathz[13] = {0,  0, 15, 15,  0, 30, 30,  0,  0,  0,  0, 30, 30};

		// double kpathx[8] = {0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0};
		// double kpathy[8] = {0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0};
		// double kpathz[8] = {0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};

		// double kpathx[4] = {0.0, 0.0, 0.5, 0.0};
		// double kpathy[4] = {0.0, 0.5, 0.5, 0.0};
		// double kpathz[4] = {0.0, 0.0, 0.0, 0.0};

		int ratio[3];

		int from[3], to[3], in[3];

		int counter = 0;
		double xout;

		// TESTING  ==============================================================================================================

		for (int p = 0; p < kpath_length-1; p++){
			
			// x component
			if (kpathx[p+1] > kpathx[p]){
				if (p == 0) {from[0] = static_cast<int>(kpathx[p] * params::Lx);}
				if (p >= 1) {from[0] = static_cast<int>(kpathx[p] * params::Lx)+1;}
				to[0] = static_cast<int>(kpathx[p+1] * params::Lx);
				in[0] = 1; //static_cast<int>(std::abs(kpathx[p+1] - kpathx[p])/(kpathx[p+1] - kpathx[p]));
			}
			else if (kpathx[p+1] == kpathx[p]){
				from[0] = static_cast<int>(kpathx[p] * params::Lx);
				to[0] = static_cast<int>(kpathx[p] * params::Lx)+1;
				in[0] = 0;
				if (from[0] != 0){
					from[0] += -1;
					to[0] += -1;
				}
			}
			else {
				from[0] = static_cast<int>(kpathx[p] * params::Lx)-2;
				to[0] = static_cast<int>(kpathx[p+1] * params::Lx)-1;
				in[0] = -1; //static_cast<int>(std::abs(kpathx[p+1] - kpathx[p])/(kpathx[p+1] - kpathx[p]));
			}

			// y component
			if (kpathy[p+1] > kpathy[p]){
				if (p == 0) {from[1] = static_cast<int>(kpathy[p] * params::Ly);}
				if (p >= 1) {from[1] = static_cast<int>(kpathy[p] * params::Ly)+1;}
				to[1] = static_cast<int>(kpathy[p+1] * geom::Iy);
				in[1] = 1; //static_cast<int>(std::abs(kpathy[p+1] - kpathy[p])/(kpathy[p+1] - kpathy[p]));
			}
			else if (kpathy[p+1] == kpathy[p]){
				from[1] = static_cast<int>(kpathy[p] * params::Ly);
				to[1] = static_cast<int>(kpathy[p] * params::Ly)+1;
				in[1] = 0;
				if (from[1] != 0){
					from[1] += -1;
					to[1] += -1;
				}
			}
			else {
				from[1] = static_cast<int>(kpathy[p] * params::Ly)-2;
				to[1] = static_cast<int>(kpathy[p+1] * params::Ly)-1;
				in[1] = -1; //static_cast<int>(std::abs(kpathy[p+1] - kpathy[p])/(kpathy[p+1] - kpathy[p]));
			}

			// z component
			if (kpathz[p+1] > kpathz[p]){
				if (p == 0) {from[2] = static_cast<int>(kpathz[p] * params::Lz);}
				if (p >= 1) {from[2] = static_cast<int>(kpathz[p] * params::Lz)+1;}
				to[2] = static_cast<int>(kpathz[p+1] * params::Lz);
				in[2] = 1; //static_cast<int>(std::abs(kpathy[p+1] - kpathy[p])/(kpathy[p+1] - kpathy[p]));
			}
			else if (kpathz[p+1] == kpathz[p]){
				from[2] = static_cast<int>(kpathz[p] * params::Lz);
				to[2] = static_cast<int>(kpathz[p] * params::Lz)+1;
				in[2] = 0;
				if (from[2] != 0){
					from[2] += -1;
					to[2] += -1;
				}
			}
			else {
				from[2] = static_cast<int>(kpathz[p] * params::Lz)-2;
				to[2] = static_cast<int>(kpathz[p+1] * params::Lz)-1;
				in[2] = -1; //static_cast<int>(std::abs(kpathy[p+1] - kpathy[p])/(kpathy[p+1] - kpathy[p]));
			}
			// DEBUGGER;
			// if (kpathx[p+1] > kpathx[p]){
			// 	from[0] = kpathx[p];
			// 	to[0] = kpathx[p+1];
			// 	in[0] = 1;
			// }
			// else if (kpathx[p+1] == kpathx[p]){
			// 	from[0] = kpathx[p];
			// 	to[0] = kpathx[p]+1;
			// 	in[0] = 0;
			// }
			// else {
			// 	from[0] = kpathx[p];
			// 	to[0] = kpathx[p+1];
			// 	in[0] = -1;
			// }

			// // y component
			// if (kpathy[p+1] > kpathy[p]){
			// 	from[1] = kpathy[p];
			// 	to[1] = kpathy[p+1];
			// 	in[1] = 1;
			// }
			// else if (kpathy[p+1] == kpathy[p]){
			// 	from[1] = kpathy[p];
			// 	to[1] = kpathy[p]+1;
			// 	in[1] = 0;
			// }
			// else {
			// 	from[1] = kpathy[p];
			// 	to[1] = kpathy[p+1];
			// 	in[1] = -1;
			// }

			// // z component
			// if (kpathz[p+1] > kpathz[p]){
			// 	from[2] = kpathz[p];
			// 	to[2] = kpathz[p+1];
			// 	in[2] = 1;
			// 	if (p == 10){
			// 		in[2] = 5;
			// 	}
			// }
			// else if (kpathz[p+1] == kpathz[p]){
			// 	from[2] = kpathz[p];
			// 	to[2] = kpathz[p]+1;	
			// 	in[2] = 0;
			// }
			// else {
			// 	from[2] = kpathz[p];
			// 	to[2] = kpathz[p+1];
			// 	in[2] = -1;
			// 	if (p == 6){
			// 		in[2] = -5;
			// 	}
			// }


			// std::cout << "from = " << from[0] << " " << from[1] << " " << from[2] << std::endl;
			// std::cout << "to = " << to[0] << " " << to[1] << " " << to[2] << std::endl;
			// std::cout << "in = " << in[0] << " " << in[1] << " " << in[2] << std::endl;

			int a = from[0];
			int b = from[1];
			int c = from[2];


			int largestdiff = std::abs(from[0]- to[0]);

			for(int i = 1; i < 3; ++i){
				// Change < to > if you want to find the smallest element
				if (largestdiff < std::abs(from[i]- to[i])){
					largestdiff = std::abs(from[i]- to[i]);
				}
			}

			// if (p == 10) largestdiff = 7;
			// if (p == 6) largestdiff = 7;

			for (int x = 0; x < largestdiff; x++){
				
				for (int np = 0; np < icount; np++){
					FFTstcfT(np)[REAL] = FFTstcfarray_test(np,a,b,c)[REAL];
					FFTstcfT(np)[IMAG] = FFTstcfarray_test(np,a,b,c)[IMAG];
				}
				// for (int np = 0; np < icount; np++){
				// 	FFTstcfT(np)[REAL] = FFTstcfarray(np,a)[REAL];
				// 	FFTstcfT(np)[IMAG] = FFTstcfarray(np,a)[IMAG];
				// }
				std::cout << " TEST " << /*largestdiff << " " << p << " " <<*/ a << " " << b << " " << c << std::endl;


				if (a != to[0]){
					a += in[0];
				}
				if (b != to[1]){
					b += in[1];
				}
				if (c != to[2]){
					c += in[2];
				}

				fftw_execute(Stime);

				// Create output files for k vectors
				std::stringstream sstr2;
				sstr2 << params::filepath_sw;
				sstr2 << "kx" << std::setw(4) << std::setfill('0') << counter;
				// sstr2 << "ky" << std::setw(4) << std::setfill('0') << b;
				// sstr2 << "kz" << std::setw(4) << std::setfill('0') << c;
				sstr2 << ".txt";

				std::ofstream kzout;
				kzout.open(sstr2.str());
				kzout << std::setprecision(10);

				for (int j = 0; j < icount/2; j++){


					// output one side spectrum;
					j1 = icount-j-1;     
					os1 = stcfT(j)[REAL] * stcfT(j)[REAL] + stcfT(j)[IMAG] * stcfT(j)[IMAG];
					os2 = stcfT(j1)[REAL] * stcfT(j1)[REAL] + stcfT(j1)[IMAG] * stcfT(j1)[IMAG];

					os[j] = os1 + os2;
				}

				double largest = os[0];
				double index;

				// Find largest value in k_z array
				for (int jj = 1; jj < icount/2; jj++){
					if (largest < os[jj]){
						largest = os[jj];
						index = jj;
					}
				}

				// output peaks as frequency values to a file for comparison against LSWT
				double ai = in[0]/params::a1;
				double bi = in[1]/params::b1;
				double ci = in[2]/params::c1;
				double mag = sqrt(ai*ai + bi*bi + ci*ci);
				xout += mag;

				peaksout << xout << " " << index * ( 1 / (params::dt_spinwaves * Npoints)) << "\n";
				counter++;

				for (int kk = 0; kk < icount/2; kk++){
					os[kk] /= largest;
					freq = kk * freqstep;
					// kzout << os[kk] << "\n";
				}

				int n = icount/2;
				int s = icount/2 + icount/2 - 1;
				int mean = icount/4;
				double sg = params::sg_spinwaves;
				double sum = 0;
				double g[n];
				double w[s];

				for (int i = 0; i < n; i++){     
					g[i] = exp(-1 * ((i - mean) * (i - mean)) / (2 * sg * sg));
					sum += g[i];
				}

				//Ensure sum of gaussian is 0
				for (int i = 0; i < n; i++){
					g[i] /= sum;         
				}        

				// Convolution
				for (int k = 0; k < s; k++){
					w[k] = 0;

					for (int i = 0; i < s; i++){

						//w[k] += in[i] * g[i];
						if (k-i >= 0 && k-i < n && i < n){
							w[k] += os[i] * g[k-i];
						}
					}

				}

				largest = w[0];
				for (int jj = 0; jj < s; jj++){
					if (largest < w[jj]){
						largest = w[jj];
					}
				}
				for (int kk = 0; kk < s; kk++){
					w[kk] /= largest;
				}

				for (int kk = 0; kk < n; kk++){
					if (kk - mean >= 0){
						kzout << w[kk] << "\n";
					}
				}

				kzout << std::flush;
				kzout.close(); 

			}
		}
			
		// ==================================================================================================================================

			fftw_destroy_plan(Slat);
			fftw_destroy_plan(Stime);

			peaksout << std::flush;
			peaksout.close();

		}

	}
