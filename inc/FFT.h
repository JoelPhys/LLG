#ifndef __FFT_H__
#define __FFT_H__

    #include <iostream>
    #include <sstream>
    #include <cmath>
    #include <fstream>
    #include <vector>
    #include <iomanip>
    #include <random>
    #include <algorithm>
    #include "fftw3.h"
    // #include "params.h"
    #include "params1.h"
    #include "array2d.h"
    #include "array3d.h"
    #include "array4d.h"
    #include "mathfuncs.h"

    namespace FFT {

        // Macros for FFT
        #define REAL 0
        #define IMAG 1


        Array2D<double> H_thermal;
        Array2D<double> Delta_S;
        Array3D<double> Sx, Sy, Sz;
        Array3D<double> S_dash_normedx, S_dash_normedy, S_dash_normedz;
        Array3D<double> Hxi, Hyi, Hzi;
        Array3D<double> Jxx, Jyy, Jzz; 
        Array3D<fftw_complex> FFTJxx, FFTJyy, FFTJzz;
        Array3D<fftw_complex> Hxk, Hyk, Hzk;
        Array3D<fftw_complex> FFTSx, FFTSy, FFTSz;
        Array3D<fftw_complex> FFTSx1, FFTSy1, FFTSz1;
        Array3D<int> Scount;
        fftw_plan planFFTSx;
        fftw_plan planFFTSy; 
        fftw_plan planFFTSz;
        fftw_plan planFFTSx1;
        fftw_plan planFFTSy1;
        fftw_plan planFFTSz1;
        fftw_plan planFFTHx;
        fftw_plan planFFTHy;
        fftw_plan planFFTHz;

        double Sx1d[params::Nspins];
        double Sy1d[params::Nspins];
        double Sz1d[params::Nspins];
        double Hx1d[params::Nspins];
        double Hy1d[params::Nspins];
        double Hz1d[params::Nspins];
        double S_dash_normedx1d[params::Nspins];
        double S_dash_normedy1d[params::Nspins];
        double S_dash_normedz1d[params::Nspins];
        double H_ani[3];
        double H_ani_dash[3];
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

        double guassian_vals[3];
        std::normal_distribution<double> distribution(0.0,1.0);
        std::random_device device;
        std::mt19937 generator(device());

        void IntialisePointersFFT(){
            H_thermal.resize(params::Nspins, 3);
            Delta_S.resize(params::Nspins, 3);
            H_thermal.IFill(0);
            Delta_S.IFill(0);
        }


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

            Jxx.resize(params::xdim, params::ydim, params::zdim);
            Jyy.resize(params::xdim, params::ydim, params::zdim);
            Jzz.resize(params::xdim, params::ydim, params::zdim);
            FFTJxx.resize(params::xdim, params::ydim, params::zdimC);
            FFTJyy.resize(params::xdim, params::ydim, params::zdimC);
            FFTJzz.resize(params::xdim, params::ydim, params::zdimC);

            fftw_plan planFFTJxx = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Jxx.ptr(), FFTJxx.ptr(), FFTW_MEASURE);
            fftw_plan planFFTJyy = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Jyy.ptr(), FFTJyy.ptr(), FFTW_MEASURE);
            fftw_plan planFFTJzz = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Jzz.ptr(), FFTJzz.ptr(), FFTW_MEASURE);

            Jxx.IFill(0);
            Jyy.IFill(0);
            Jzz.IFill(0);
            FFTJxx.IFill(0);
            FFTJyy.IFill(0);
            FFTJzz.IFill(0);

            for (int i = 0; i < jxx.size(); i++){
                Jxx(modfunc(params::xdim, Nx[i]), modfunc(params::ydim, Ny[i]), modfunc(params::zdim, Nz[i])) = jxx[i];
                Jyy(modfunc(params::xdim, Nx[i]), modfunc(params::ydim, Ny[i]), modfunc(params::zdim, Nz[i])) = jyy[i];
                Jzz(modfunc(params::xdim, Nx[i]), modfunc(params::ydim, Ny[i]), modfunc(params::zdim, Nz[i])) = jzz[i];
            } 

            fftw_execute(planFFTJxx);
            fftw_execute(planFFTJyy);
            fftw_execute(planFFTJzz);
            fftw_destroy_plan(planFFTJxx);
            fftw_destroy_plan(planFFTJyy);
            fftw_destroy_plan(planFFTJzz);
        }

        void FFT_initialise(){

            Sx.resize(params::xdim, params::ydim, params::zdim);
            Sy.resize(params::xdim, params::ydim, params::zdim);
            Sz.resize(params::xdim, params::ydim, params::zdim);
            Scount.resize(params::xdim, params::ydim, params::zdim);
            S_dash_normedx.resize(params::xdim, params::ydim, params::zdim);
            S_dash_normedy.resize(params::xdim, params::ydim, params::zdim);
            S_dash_normedz.resize(params::xdim, params::ydim, params::zdim);  
            FFTSx.resize(params::xdim, params::ydim, params::zdimC);
            FFTSy.resize(params::xdim, params::ydim, params::zdimC);
            FFTSz.resize(params::xdim, params::ydim, params::zdimC);
            FFTSx1.resize(params::xdim, params::ydim, params::zdimC);
            FFTSy1.resize(params::xdim, params::ydim, params::zdimC);
            FFTSz1.resize(params::xdim, params::ydim, params::zdimC);
            Hxk.resize(params::xdim, params::ydim, params::zdimC);
            Hyk.resize(params::xdim, params::ydim, params::zdimC);
            Hzk.resize(params::xdim, params::ydim, params::zdimC);
            Hxi.resize(params::xdim, params::ydim, params::zdim);
            Hyi.resize(params::xdim, params::ydim, params::zdim);
            Hzi.resize(params::xdim, params::ydim, params::zdim);

            planFFTSx = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Sx.ptr(), FFTSx.ptr(), FFTW_MEASURE);
            planFFTSy = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Sy.ptr(), FFTSy.ptr(), FFTW_MEASURE);
            planFFTSz = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, Sz.ptr(), FFTSz.ptr(), FFTW_MEASURE);
            planFFTSx1 = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, S_dash_normedx.ptr(), FFTSx1.ptr(), FFTW_MEASURE);
            planFFTSy1 = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, S_dash_normedy.ptr(), FFTSy1.ptr(), FFTW_MEASURE);
            planFFTSz1 = fftw_plan_dft_r2c_3d(params::xdim, params::ydim, params::zdim, S_dash_normedz.ptr(), FFTSz1.ptr(), FFTW_MEASURE);
            planFFTHx = fftw_plan_dft_c2r_3d(params::xdim, params::ydim, params::zdim, Hxk.ptr(), Hxi.ptr(), FFTW_MEASURE);
            planFFTHy = fftw_plan_dft_c2r_3d(params::xdim, params::ydim, params::zdim, Hyk.ptr(), Hyi.ptr(), FFTW_MEASURE);
            planFFTHz = fftw_plan_dft_c2r_3d(params::xdim, params::ydim, params::zdim, Hzk.ptr(), Hzi.ptr(), FFTW_MEASURE);

            S_dash_normedx.IFill(0);
            S_dash_normedy.IFill(0);
            S_dash_normedz.IFill(0);
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
            Scount.IFill(0);
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
                        if ((Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            Sx1d[Scount(x,y,z)] = Sx(x,y,z);
                            Sy1d[Scount(x,y,z)] = Sy(x,y,z);
                            Sz1d[Scount(x,y,z)] = Sz(x,y,z);
                            Hx1d[Scount(x,y,z)] = Hxi(x,y,z);
                            Hy1d[Scount(x,y,z)] = Hyi(x,y,z);
                            Hz1d[Scount(x,y,z)] = Hzi(x,y,z);
                        }
                    }
                }
            }
        }

        void FFT_Heun1(double Thermal_Fluct){
            for (int a = 0; a < params::Nspins; a++){

                for (int g = 0; g < 3; g++){
                    guassian_vals[g] = distribution(generator);
                }

                H_thermal(a,0)  = guassian_vals[0] * Thermal_Fluct;
                H_thermal(a,1)  = guassian_vals[1] * Thermal_Fluct;
                H_thermal(a,2)  = guassian_vals[2] * Thermal_Fluct;

                // Uniaxial anisotropy in Z axis
                H_ani[0]=0;
                H_ani[1]=0;
                H_ani[2]=params::d_z_prime * Sz1d[a];

                H_new[0] = H_thermal(a,0) + params::H_app[0] + H_ani[0] + (Hx1d[a] * params::INVmu_s);
                H_new[1] = H_thermal(a,1) + params::H_app[1] + H_ani[1] + (Hy1d[a] * params::INVmu_s);
                H_new[2] = H_thermal(a,2) + params::H_app[2] + H_ani[2] + (Hz1d[a] * params::INVmu_s);

                ScrossP[0] = Sx1d[a];
                ScrossP[1] = Sy1d[a];
                ScrossP[2] = Sz1d[a];

                CrossP(ScrossP, H_new, CrossP1);
                CrossP(ScrossP, CrossP1, CrossP2);

                Delta_S(a,0)  = -params::lambdaPrime * (CrossP1[0] + params::lambda * CrossP2[0]);
                Delta_S(a,1)  = -params::lambdaPrime * (CrossP1[1] + params::lambda * CrossP2[1]);
                Delta_S(a,2)  = -params::lambdaPrime * (CrossP1[2] + params::lambda * CrossP2[2]);

                S_dash[0] = Sx1d[a] + (Delta_S(a,0)  * params::dtau);
                S_dash[1] = Sy1d[a] + (Delta_S(a,1)  * params::dtau);
                S_dash[2] = Sz1d[a] + (Delta_S(a,2)  * params::dtau);

                invmag = 1 / sqrt(S_dash[0] * S_dash[0] + S_dash[1] * S_dash[1] + S_dash[2] * S_dash[2]);
                S_dash_normedx1d[a] = invmag * S_dash[0];
                S_dash_normedy1d[a] = invmag * S_dash[1];
                S_dash_normedz1d[a] = invmag * S_dash[2];
            }
        }

        void FFT_Heun2(double Thermal_Fluct){
            for (int a = 0; a < params::Nspins; a++){

                H_ani_dash[0]=0;
                H_ani_dash[1]=0;
                H_ani_dash[2]=params::d_z_prime * S_dash_normedz1d[a];

                H_new_dash[0] = H_thermal(a,0)  + params::H_app[0] + H_ani_dash[0] + (Hx1d[a] * params::INVmu_s);
                H_new_dash[1] = H_thermal(a,1)  + params::H_app[1] + H_ani_dash[1] + (Hy1d[a] * params::INVmu_s);
                H_new_dash[2] = H_thermal(a,2)  + params::H_app[2] + H_ani_dash[2] + (Hz1d[a] * params::INVmu_s);

                // Calculate Corrector and Normalise
                ScrossP[0] = S_dash_normedx1d[a];
                ScrossP[1] = S_dash_normedy1d[a];
                ScrossP[2] = S_dash_normedz1d[a];

                CrossP(ScrossP, H_new_dash, CrossP3);
                CrossP(ScrossP, CrossP3, CrossP4);

                Delta_S_dash[0] = -params::lambdaPrime * (CrossP3[0] + params::lambda * CrossP4[0]);
                Delta_S_dash[1] = -params::lambdaPrime * (CrossP3[1] + params::lambda * CrossP4[1]);
                Delta_S_dash[2] = -params::lambdaPrime * (CrossP3[2] + params::lambda * CrossP4[2]); 

                S_new[0] = Sx1d[a] + params::half_dtau * (Delta_S(a,0)  + Delta_S_dash[0]);
                S_new[1] = Sy1d[a] + params::half_dtau * (Delta_S(a,1)  + Delta_S_dash[1]);
                S_new[2] = Sz1d[a] + params::half_dtau * (Delta_S(a,2)  + Delta_S_dash[2]);

                invmag1 = 1 / sqrt(S_new[0] * S_new[0] + S_new[1] * S_new[1] + S_new[2] * S_new[2]);
                Sx1d[a] = invmag1 * S_new[0];
                Sy1d[a] = invmag1 * S_new[1];
                Sz1d[a] = invmag1 * S_new[2];                
            }

            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        if ((Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            Sx(x,y,z) = Sx1d[Scount(x,y,z)];
                            Sy(x,y,z) = Sy1d[Scount(x,y,z)];
                            Sz(x,y,z) = Sz1d[Scount(x,y,z)];
                        }
                    }
                }
            }

        }

        void fft2(){
            for (int x = 0; x < params::xdim; x++){
                for (int y = 0; y < params::ydim; y++){
                    for (int z = 0; z < params::zdim; z++){
                        if ((Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            S_dash_normedx(x,y,z) = S_dash_normedx1d[Scount(x,y,z)];
                            S_dash_normedy(x,y,z) = S_dash_normedy1d[Scount(x,y,z)];
                            S_dash_normedz(x,y,z) = S_dash_normedz1d[Scount(x,y,z)];
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

                        if ((Scount(x,y,z) == 0) && (x > 0 || y > 0 || z > 0 )){

                        }
                        else {
                            Hx1d[Scount(x,y,z)] = Hxi(x,y,z);
                            Hy1d[Scount(x,y,z)] = Hyi(x,y,z);
                            Hz1d[Scount(x,y,z)] = Hzi(x,y,z);
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

#endif