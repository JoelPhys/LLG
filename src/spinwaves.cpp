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
    Array3D<fftw_complex> FFTstcfarray_test;

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

        stcf.resize(geom::Ix,geom::Iy,geom::Iz);
        FFTstcf.resize(geom::Ix,geom::Iy,geom::IzC);
        stcfT.resize(Npoints);
        FFTstcfT.resize(Npoints);
        FFTstcfarray.resize(Npoints, geom::Ix);

        //TESTING WITH 2D outputs
        FFTstcfarray_test.resize(Npoints, geom::Ix, geom::IzC);

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
        // stcf1d.IFill(0);
        // FFTstcf1d.IFill(0);

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
                        stcf(lval,mval,nval) = neigh::Sy1d(geom::Scount(lval,mval,nval)) * neigh::Sy1d(geom::Scount(lval,mval,nval)) + neigh::Sz1d(geom::Scount(lval,mval,nval)) * neigh::Sz1d(geom::Scount(lval,mval,nval));
                    }
                }
            }
        }

        // compute fft
        fftw_execute(Slat);

        // //TESTING WITH MN2AU DELETE IF NOT USING
        // for (int l = 0; l < params::Lx; l++){
        //     stcf1d(l) = neigh::Sy1d(geom::LatCount(l,1,1,2)) * neigh::Sy1d(geom::LatCount(l,1,1,2)) + neigh::Sz1d(geom::LatCount(l,1,1,2)) * neigh::Sz1d(geom::LatCount(l,1,1,2));
        // }

        // // testing for Mn2Au
        // fftw_execute(Slat1d);

        //windowing function
        windowing = hammA - hammB * cos((2 * M_PI * icount) / (Npoints  - 1));

        double kpts = geom::Ix;      

        // for (int j = 0; j < geom::Ix; j++){

        //     if (j < kpts){
        //         FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf(j,0,0)[REAL];
        //         FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf(j,0,0)[IMAG];
        //         file_spnwvs << FFTstcf(j,0,0)[REAL] << " " << FFTstcf(j,0,0)[IMAG] << "\t";

        //         //TESTING MN2AU
        //         // FFTstcfarray(icount,j)[REAL] = windowing * FFTstcf1d(j)[REAL];
        //         // FFTstcfarray(icount,j)[IMAG] = windowing * FFTstcf1d(j)[IMAG];
        //         // file_spnwvs << FFTstcf1d(j)[REAL] << " " << FFTstcf1d(j)[IMAG] << "\t";

        //     }
        // }

        for (int j = 0; j < geom::Ix; j++){
            for (int k = 0; k < geom::IzC; k++){

                FFTstcfarray_test(icount,j,k)[REAL] = windowing * FFTstcf(j,0,k)[REAL];
                FFTstcfarray_test(icount,j,k)[IMAG] = windowing * FFTstcf(j,0,k)[IMAG];

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
        peakstring << "output/spinwaves/peaks.txt";
        std::ofstream peaksout;
        peaksout.open(peakstring.str());
        peaksout << std::setprecision(10);

        // loop over all elememts in k array
        double kpts = geom::Ix; 


        // TESTING  ==============================================================================================================
        for (int j = 0; j < geom::Ix; j++){
            for (int k = 0; k < geom::IzC; k++){
                
                for (int np = 0; np < icount; np++){
                    FFTstcfT(np)[REAL] = FFTstcfarray_test(np,j,k)[REAL];
                    FFTstcfT(np)[IMAG] = FFTstcfarray_test(np,j,k)[IMAG];
                }

                fftw_execute(Stime);

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
                peaksout << index * ( 1 / (params::dt_spinwaves * Npoints)) << " ";
                
            }
            peaksout << "\n";
        }
        // ==================================================================================================================================

        // for (int z = 0; z < kpts; z++){


        //     // Create output files for k vectors
        //     std::stringstream sstr2;
        //     sstr2 << "output/spinwaves/kz";
        //     sstr2 << std::setw(4) << std::setfill('0') << z;
        //     sstr2 << ".txt";

        //     std::ofstream kzout;
        //     kzout.open(sstr2.str());
        //     kzout << std::setprecision(10);

        //     // convert index to k space
        //     kpoint = z * ( ( (2 * M_PI)/ ( params::a1)) / ( geom::IzC ));


        //     for (int j = 0; j < icount; j++){
        //         FFTstcfT(j)[REAL] = FFTstcfarray(j,z)[REAL];
        //         FFTstcfT(j)[IMAG] = FFTstcfarray(j,z)[IMAG];
        //     }

        //     fftw_execute(Stime);

        //     for (int j = 0; j < icount/2; j++){


        //         // output one side spectrum;
        //         j1 = icount-j-1;     
        //         os1 = stcfT(j)[REAL] * stcfT(j)[REAL] + stcfT(j)[IMAG] * stcfT(j)[IMAG];
        //         os2 = stcfT(j1)[REAL] * stcfT(j1)[REAL] + stcfT(j1)[IMAG] * stcfT(j1)[IMAG];

        //         os[j] = os1 + os2;
        //     }

        //     double largest = os[0];
        //     double index;

        //     // Find largest value in k_z array
        //     for (int jj = 1; jj < icount/2; jj++){
        //         if (largest < os[jj]){
        //             largest = os[jj];
        //             index = jj;
        //         }
        //     }

        //     // output peaks as frequency values to a file for comparison against LSWT
        //     peaksout << kpoint << " " << index * ( 1 / (params::dt_spinwaves * Npoints)) << "\n";

        //     for (int kk = 0; kk < icount/2; kk++){
        //         os[kk] /= largest;
        //         freq = kk * freqstep;
        //         // kzout << os[kk] << "\n";
        //     }

        //     int n = icount/2;
        //     int s = icount/2 + icount/2 - 1;
        //     int mean = icount/4;
        //     double sg = 10;
        //     double sum = 0;
        //     double g[n];
        //     double w[s];
            
        //     for (int i = 0; i < n; i++){	
        //         g[i] = exp(-1 * ((i - mean) * (i - mean)) / (2 * sg * sg));
        //         sum += g[i];
        //     }

        //     //Ensure sum of gaussian is 0
        //     for (int i = 0; i < n; i++){
        //         g[i] /= sum;		
        //     }	

        //     // Convolution
        //     for (int k = 0; k < s; k++){
        //         w[k] = 0;

        //         for (int i = 0; i < s; i++){
                
        //             //w[k] += in[i] * g[i];
        //             if (k-i >= 0 && k-i < n && i < n){
        //                 w[k] += os[i] * g[k-i];
        //             }
        //         }

        //     }

        //     largest = w[0];
        //     for (int jj = 1; jj < s; jj++){
        //         if (largest < w[jj]){
        //             largest = w[jj];
        //         }
        //     }
        //     for (int kk = 0; kk < s; kk++){
        //         w[kk] /= largest;
        //     }

        //     for (int kk = 0; kk < n; kk++){
        //         if (kk - mean >= 0){
        //             kzout << w[kk] << "\n";
        //         }
        //     }

        //         kzout << std::flush;
        //         kzout.close();


        // }

        fftw_destroy_plan(Slat);
        fftw_destroy_plan(Stime);

        peaksout << std::flush;
        peaksout.close();

    }

}
