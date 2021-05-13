#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include "../inc/fftw3.h"
// #include "params.h"
#include "../inc/params1.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/mathfuncs.h"
#include "../inc/params.h"
#include "../inc/geom.h"
#include "../inc/spinwaves.h"
#include "../inc/NeighbourList.h"

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
    // Array<double> stcf_new;
    // Array<fftw_complex> FFTstcf_new;
    // fftw_plan Slat_new;

    fftw_plan Slat;
    fftw_plan Stime;

    std::ofstream file_spnwvs;

    double windowing = 0;

    int Npoints;
    int icount;

    double hammA = 0.54;
    double hammB = 0.46;

    double norm;

    int lval, mval, nval;

    void initialiseFFT(){
        Npoints  = ceil((params::Nt - params::start) / (params::dt_spinwaves / params::dt));
        norm = 1 / (geom::Ix * geom::Iy * geom::Iz * Npoints);

        std::cout << "Npoints = " << Npoints << std::endl;
        icount = 0;


        //testing for Mn2Au - Let's try going along z axis using q = 1 and q = 2.
        // stcf_new.resize(params::Lz * 2);
        // FFTstcf_new.resize(params::Lz + 1); // Array size is x/2 + 1 following FFT
        // Slat_new = fftw_plan_dft_r2c_1d(params::Lz * 2, stcf_new.ptr(), FFTstcf_new.ptr(), FFTW_MEASURE);


        stcf.resize(geom::Ix,geom::Iy,geom::Iz);
        FFTstcf.resize(geom::Ix,geom::Iy,geom::IzC);
        stcfT.resize(Npoints);
        FFTstcfT.resize(Npoints);
        // FFTstcfarray.resize(params::Nt,3*geom::Ix);
        FFTstcfarray.resize(params::Nt,geom::IzC);

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
                        stcf(lval,mval,nval) = neigh::Sx1d(geom::Scount(lval,mval,nval)) * neigh::Sx1d(geom::Scount(lval,mval,nval)) + neigh::Sy1d(geom::Scount(lval,mval,nval)) * neigh::Sy1d(geom::Scount(lval,mval,nval));
                    }
                }
            }
        }

        // testing for Mn2Au
        // int cnt = 0;
        // double x2;
        // double y2;
        // for (int a = 0; a < params::Nspins; a++){
        //     if (a < params::Lz*4 && (( a % 4 == 1) || (a % 4 == 2))){
        //         x2 = neigh::Sx1d(a) * neigh::Sx1d(a);
        //         y2 = neigh::Sy1d(a) * neigh::Sy1d(a);
        //         stcf_new(cnt) = x2 + y2;
        //         cnt++;
        //     }
        // }
        // fftw_execute(Slat_new);

        // compute fft
        fftw_execute(Slat);

        //windowing function
        windowing = hammA - hammB * cos((2 * M_PI * icount) / (Npoints  - 1));

        //lets try and program a k path
        double kpath[5][3] = {{0,0,0},{1,0,0},{1,1,0},{0,0,0},{1,1,1}};

        // for loop for difference between the k points
       // for (int i =0; i < 5; i++){
       //     int dx = kpath[i+1][0] - kpath[i][0];
       //     int dy = kpath[i+1][1] - kpath[i][1];
       //     int dz = kpath[i+1][2] - kpath[i][2];
       // }

        double kpts = geom::IzC;      

        for (int j = 0; j < kpts; j++){

            if (j < geom::IzC){
                FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf(0,0,j)[REAL];
                FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf(0,0,j)[IMAG];
                file_spnwvs << FFTstcf(0,0,j)[REAL] << " " << FFTstcf(0,0,j)[IMAG] << "\t";
            }
            // else if (geom::Ix <= j && j < 2*geom::Ix){
            //     FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf(j-geom::Ix,geom::Iy-1,0)[REAL];
            //     FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf(j-geom::Ix,geom::Iy-1,0)[IMAG];
            //     file_spnwvs << FFTstcf(j-geom::Ix,geom::Iy-1,0)[REAL] << " " << FFTstcf(j-geom::Ix,geom::Iy-1,0)[IMAG] << "\t";
            // }
            // else if (2*geom::Ix <= j && j < 3*geom::Ix){
            //     FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf(3*geom::Ix-1-j,3*geom::Iy-1-j,0)[REAL];
            //     FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf(3*geom::Ix-1-j,3*geom::Iy-1-j,0)[IMAG];
            //     file_spnwvs << FFTstcf(3*geom::Ix-1-j,3*geom::Iy-1-j,0)[REAL] << " " << FFTstcf(3*geom::Ix-1-j,3*geom::Iy-1-j,0)[IMAG] << "\t";
            // }
            // else {
            //     std::cout << "ERROR" << std::endl;
            // }
        }

        // testing for Mn2Au
        // double kpts = params::Lz+1;      

        // for (int j = 0; j < kpts; j++){
        //         FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf_new(j)[REAL];
        //         FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf_new(j)[IMAG];
        //         //file_spnwvs << FFTstcf(0,0,j)[REAL] << " " << FFTstcf(0,0,j)[IMAG] << "\t";
        // }

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
        peakstring << "output/spinwaves/peaks.txt";
        std::ofstream peaksout;
        peaksout.open(peakstring.str());
        peaksout << std::setprecision(10);

        // loop over all elememts in k array
        double kpts = geom::IzC; 

        for (int z = 0; z < kpts; z++){

            for (int j = 0; j < icount; j++){
                FFTstcfT(j)[REAL] = FFTstcfarray(j,z)[REAL];
                FFTstcfT(j)[IMAG] = FFTstcfarray(j,z)[IMAG];
            }

            fftw_execute(Stime);

            for (int j = 0; j < icount/2; j++){


                // output one side spectrum;
                j1 = icount-j-1;     
                os1 = stcfT(j)[REAL] * stcfT(j)[REAL] + stcfT(j)[IMAG] * stcfT(j)[IMAG];
                os2 = stcfT(j1)[REAL] * stcfT(j1)[REAL] + stcfT(j1)[IMAG] * stcfT(j1)[IMAG];

                os[j] = os1 + os2;
            }

            // convert index to k space
            kpoint = z * ( ( (2 * M_PI)/ ( params::a1)) / ( geom::IzC ));

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
            peaksout << kpoint << " " << index * ( 1 / (params::dt_spinwaves * Npoints)) << "\n";

            // Create output files for k vectors
            std::stringstream sstr2;
            sstr2 << "output/spinwaves/kz";
            sstr2 << std::setw(3) << std::setfill('0') << z;
            sstr2 << ".txt";

            std::ofstream kzout;
            kzout.open(sstr2.str());
            kzout << std::setprecision(10);

            for (int kk = 0; kk < icount/2; kk++){
                os[kk] /= largest;
                freq = kk * freqstep;
                // kzout << os[kk] << "\n";
            }

            int n = icount/2;
            int s = icount/2 + icount/2 - 1;
            int mean = icount/4;
            double sg = 10;
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
            for (int jj = 1; jj < s; jj++){
                if (largest < w[jj]){
                    largest = w[jj];
                }
            }
            for (int kk = 0; kk < s; kk++){
                w[kk] /= largest;
            }

        // for (int kk = 0; kk < n; kk++){
        //     kzout << os[kk] << " ";
        //     kzout << w[kk + mean] << "\n";
        // }

        for (int kk = 0; kk < n; kk++){
            if (kk - mean >= 0){
                kzout << w[kk] << std::endl;
            }
	    }

            kzout << std::flush;
            kzout.close();


        }

        fftw_destroy_plan(Slat);
        fftw_destroy_plan(Stime);

        peaksout << std::flush;
        peaksout.close();

    }

}
