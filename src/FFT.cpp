#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include "fftw3.h"
#include "../inc/FFT.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"
#include "../inc/mathfuncs.h"

namespace FFT {

    // Macros for FFT
    #define REAL 0
    #define IMAG 1

    Array2D<double> H_thermal;
    Array2D<double> Delta_S;
    Array3D<double> Sx, Sy, Sz;
    Array3D<double> Sxd, Syd, Szd;
    Array3D<double> Hxi, Hyi, Hzi;
    Array3D<double> Jxx, Jyy, Jzz; 
    Array3D<fftw_complex> FFTJxx, FFTJyy, FFTJzz;
    Array3D<fftw_complex> Hxk, Hyk, Hzk;
    Array3D<fftw_complex> FFTSx, FFTSy, FFTSz;
    Array3D<fftw_complex> FFTSx1, FFTSy1, FFTSz1;
    fftw_plan planFFTSx;
    fftw_plan planFFTSy; 
    fftw_plan planFFTSz;
    fftw_plan planFFTSx1;
    fftw_plan planFFTSy1;
    fftw_plan planFFTSz1;
    fftw_plan planFFTHx;
    fftw_plan planFFTHy;
    fftw_plan planFFTHz;

    void FFT_InteractionMatrix(){


            std::ifstream input("FePt_Jij.dat");
            double a, b, c, d, e, f;
            std::vector<double> Nx;
            std::vector<double> Ny;
            std::vector<double> Nz;
            std::vector<double> jxx;
            std::vector<double> jyy;
            std::vector<double> jzz;

            while (input >> a >> b >> c >> d >> e >> f)
            {
                Nx.push_back(a);
                Ny.push_back(b);
                Nz.push_back(c);
                jxx.push_back(d);
                jyy.push_back(e);
                jzz.push_back(f);
                
            }
            input.close();

            Jxx.resize(geom::Ix, geom::Iy, geom::Iz);
            Jyy.resize(geom::Ix, geom::Iy, geom::Iz);
            Jzz.resize(geom::Ix, geom::Iy, geom::Iz);
            FFTJxx.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTJyy.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTJzz.resize(geom::Ix, geom::Iy, geom::IzC);

            fftw_plan planFFTJxx = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Jxx.ptr(), FFTJxx.ptr(), FFTW_MEASURE);
            fftw_plan planFFTJyy = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Jyy.ptr(), FFTJyy.ptr(), FFTW_MEASURE);
            fftw_plan planFFTJzz = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Jzz.ptr(), FFTJzz.ptr(), FFTW_MEASURE);

            Jxx.IFill(0);
            Jyy.IFill(0);
            Jzz.IFill(0); 
            FFTJxx.IFill(0);
            FFTJyy.IFill(0);
            FFTJzz.IFill(0);

            for (int i = 0; i < jxx.size(); i++){
                Jxx(modfunc(geom::Ix, Nx[i]), modfunc(geom::Iy, Ny[i]), modfunc(geom::Iz, Nz[i])) = jxx[i];
                Jyy(modfunc(geom::Ix, Nx[i]), modfunc(geom::Iy, Ny[i]), modfunc(geom::Iz, Nz[i])) = jyy[i];
                Jzz(modfunc(geom::Ix, Nx[i]), modfunc(geom::Iy, Ny[i]), modfunc(geom::Iz, Nz[i])) = jzz[i];
            } 

            fftw_execute(planFFTJxx);
            fftw_execute(planFFTJyy);
            fftw_execute(planFFTJzz);
            fftw_destroy_plan(planFFTJxx);
            fftw_destroy_plan(planFFTJyy);
            fftw_destroy_plan(planFFTJzz);
        }

        void FFT_initialise(){
            Sx.resize(geom::Ix, geom::Iy, geom::Iz);
            Sy.resize(geom::Ix, geom::Iy, geom::Iz);
            Sz.resize(geom::Ix, geom::Iy, geom::Iz);
            Sxd.resize(geom::Ix, geom::Iy, geom::Iz);
            Syd.resize(geom::Ix, geom::Iy, geom::Iz);
            Szd.resize(geom::Ix, geom::Iy, geom::Iz);  
            FFTSx.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTSy.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTSz.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTSx1.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTSy1.resize(geom::Ix, geom::Iy, geom::IzC);
            FFTSz1.resize(geom::Ix, geom::Iy, geom::IzC);
            Hxk.resize(geom::Ix, geom::Iy, geom::IzC);
            Hyk.resize(geom::Ix, geom::Iy, geom::IzC);
            Hzk.resize(geom::Ix, geom::Iy, geom::IzC);
            Hxi.resize(geom::Ix, geom::Iy, geom::Iz);
            Hyi.resize(geom::Ix, geom::Iy, geom::Iz);
            Hzi.resize(geom::Ix, geom::Iy, geom::Iz);

            planFFTSx = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Sx.ptr(), FFTSx.ptr(), FFTW_MEASURE);
            planFFTSy = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Sy.ptr(), FFTSy.ptr(), FFTW_MEASURE);
            planFFTSz = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Sz.ptr(), FFTSz.ptr(), FFTW_MEASURE);
            planFFTSx1 = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Sxd.ptr(), FFTSx1.ptr(), FFTW_MEASURE);
            planFFTSy1 = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Syd.ptr(), FFTSy1.ptr(), FFTW_MEASURE);
            planFFTSz1 = fftw_plan_dft_r2c_3d(geom::Ix, geom::Iy, geom::Iz, Szd.ptr(), FFTSz1.ptr(), FFTW_MEASURE);
            planFFTHx = fftw_plan_dft_c2r_3d(geom::Ix, geom::Iy, geom::Iz, Hxk.ptr(), Hxi.ptr(), FFTW_MEASURE);
            planFFTHy = fftw_plan_dft_c2r_3d(geom::Ix, geom::Iy, geom::Iz, Hyk.ptr(), Hyi.ptr(), FFTW_MEASURE);
            planFFTHz = fftw_plan_dft_c2r_3d(geom::Ix, geom::Iy, geom::Iz, Hzk.ptr(), Hzi.ptr(), FFTW_MEASURE);

            Sxd.IFill(0);
            Syd.IFill(0);
            Szd.IFill(0);
            FFTSx.IFill(0);
            FFTSy.IFill(0);
            FFTSz.IFill(0);
            FFTSx1.IFill(0);
            FFTSy1.IFill(0);
            FFTSz1.IFill(0);
            Hxi.IFill(0);
            Hyi.IFill(0);
            Hzi.IFill(0);
            Hxk.IFill(0);
            Hyk.IFill(0);
            Hzk.IFill(0);
        }

        void fft1(){
            fftw_execute(planFFTSx);
            fftw_execute(planFFTSy);
            fftw_execute(planFFTSz);

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdimC; z++){
                        Hxk(x,y,z)[REAL] = (FFTJxx(x,y,z)[REAL] * FFTSx(x,y,z)[REAL]) - (FFTJxx(x,y,z)[IMAG] * FFTSx(x,y,z)[IMAG]);
                        Hyk(x,y,z)[REAL] = (FFTJyy(x,y,z)[REAL] * FFTSy(x,y,z)[REAL]) - (FFTJyy(x,y,z)[IMAG] * FFTSy(x,y,z)[IMAG]);
                        Hzk(x,y,z)[REAL] = (FFTJzz(x,y,z)[REAL] * FFTSz(x,y,z)[REAL]) - (FFTJzz(x,y,z)[IMAG] * FFTSz(x,y,z)[IMAG]);
                        Hxk(x,y,z)[IMAG] = (FFTJxx(x,y,z)[REAL] * FFTSx(x,y,z)[IMAG]) + (FFTJxx(x,y,z)[IMAG] * FFTSx(x,y,z)[REAL]);
                        Hyk(x,y,z)[IMAG] = (FFTJyy(x,y,z)[REAL] * FFTSy(x,y,z)[IMAG]) + (FFTJyy(x,y,z)[IMAG] * FFTSy(x,y,z)[REAL]);
                        Hzk(x,y,z)[IMAG] = (FFTJzz(x,y,z)[REAL] * FFTSz(x,y,z)[IMAG]) + (FFTJzz(x,y,z)[IMAG] * FFTSz(x,y,z)[REAL]);
                    }
                } 
            }

            fftw_execute(planFFTHx);
            fftw_execute(planFFTHy);
            fftw_execute(planFFTHz);

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        Hxi(x,y,z) *= params::NsitesINV;
                        Hyi(x,y,z) *= params::NsitesINV;
                        Hzi(x,y,z) *= params::NsitesINV;  
                    }
                }
            }

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        if ((geom::Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            Sx1d[geom::Scount(x,y,z)] = Sx(x,y,z);
                            Sy1d[geom::Scount(x,y,z)] = Sy(x,y,z);
                            Sz1d[geom::Scount(x,y,z)] = Sz(x,y,z);
                            Hx1d[geom::Scount(x,y,z)] = Hxi(x,y,z);
                            Hy1d[geom::Scount(x,y,z)] = Hyi(x,y,z);
                            Hz1d[geom::Scount(x,y,z)] = Hzi(x,y,z);
                        }
                    }
                }
            }
        }

        void fft2(){
            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        if ((geom::Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            S_dash_normedx(x,y,z) = S_dash_normedx1d[geom::Scount(x,y,z)];
                            S_dash_normedy(x,y,z) = S_dash_normedy1d[geom::Scount(x,y,z)];
                            S_dash_normedz(x,y,z) = S_dash_normedz1d[geom::Scount(x,y,z)];
                        }
                    }
                }
            }
            
            fftw_execute(planFFTSx1);
            fftw_execute(planFFTSy1);
            fftw_execute(planFFTSz1);

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdimC; z++){
                        Hxk(x,y,z)[REAL] = (FFTJxx(x,y,z)[REAL] * FFTSx1(x,y,z)[REAL]) - (FFTJxx(x,y,z)[IMAG] * FFTSx1(x,y,z)[IMAG]);
                        Hyk(x,y,z)[REAL] = (FFTJyy(x,y,z)[REAL] * FFTSy1(x,y,z)[REAL]) - (FFTJyy(x,y,z)[IMAG] * FFTSy1(x,y,z)[IMAG]);
                        Hzk(x,y,z)[REAL] = (FFTJzz(x,y,z)[REAL] * FFTSz1(x,y,z)[REAL]) - (FFTJzz(x,y,z)[IMAG] * FFTSz1(x,y,z)[IMAG]);
                        Hxk(x,y,z)[IMAG] = (FFTJxx(x,y,z)[REAL] * FFTSx1(x,y,z)[IMAG]) + (FFTJxx(x,y,z)[IMAG] * FFTSx1(x,y,z)[REAL]);
                        Hyk(x,y,z)[IMAG] = (FFTJyy(x,y,z)[REAL] * FFTSy1(x,y,z)[IMAG]) + (FFTJyy(x,y,z)[IMAG] * FFTSy1(x,y,z)[REAL]);
                        Hzk(x,y,z)[IMAG] = (FFTJzz(x,y,z)[REAL] * FFTSz1(x,y,z)[IMAG]) + (FFTJzz(x,y,z)[IMAG] * FFTSz1(x,y,z)[REAL]);
                    }
                } 
            }

            fftw_execute(planFFTHx);
            fftw_execute(planFFTHy);
            fftw_execute(planFFTHz);

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        Hxi(x,y,z) *= params::NsitesINV;
                        Hyi(x,y,z) *= params::NsitesINV;
                        Hzi(x,y,z) *= params::NsitesINV;

                        if ((geom::Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            Hx1d[geom::Scount(x,y,z)] = Hxi(x,y,z);
                            Hy1d[geom::Scount(x,y,z)] = Hyi(x,y,z);
                            Hz1d[geom::Scount(x,y,z)] = Hzi(x,y,z);
                        }
                    }
                }
            }
        }

        void FFT_destroyplans(){
            fftw_destroy_plan(planFFTHx);
            fftw_destroy_plan(planFFTHy);
            fftw_destroy_plan(planFFTHz);
            fftw_destroy_plan(planFFTSx);
            fftw_destroy_plan(planFFTSy);
            fftw_destroy_plan(planFFTSz);
            fftw_destroy_plan(planFFTSx1);
            fftw_destroy_plan(planFFTSy1);
            fftw_destroy_plan(planFFTSz1);
        }

}