#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <time.h>
#include <fftw3.h>
#include "inc/array2d.h"
#include "inc/array3d.h"
#include "inc/array4d.h"
#include "inc/mathfuncs1.h"

// Macros for FFT
#define REAL 0
#define IMAG 1

const int Lx = 20;
const int Ly = 20;
const int Lz = 20;
const int Nq = 1;
const int ax = 2;
const int ay = 2;
const int az = 2;
const int Nspins = Nq*Lx*Ly*Lz;
const double xdim = ax*Lx;
const double ydim = ay*Ly;
const double zdim = az*Lz;
const double NsitesINV = 1/(xdim*ydim*zdim);
const int zdimC = zdim/2+1;

// GLOBAL CONSTANT
const double dt = 1E-16; 
const double Nt = 10000; 
const double lambda = 1; 
const double lambdaPrime = 1 / (1+(lambda*lambda));
const double k_B = 1.3807e-23; 
const double mu_b = 9.2740e-24;
const double gamma1 = 1.76E11;
const double dtau = gamma1 * dt;
const double half_dtau = dtau * 0.5;

// MATERIAL CONSTANTS
const double d_z = 0; // Anisotropy Energy
const double mu_s = 1.5 * mu_b; // Atomic Spin Moment
const double INVmu_s = 1 / mu_s;
const double d_z_prime = 2 * ( d_z / mu_s );
const double thermal_const = sqrt( (2 * lambda * k_B)  / (mu_s * dtau) );
double H_app[3] = {0,0,0};

//atom sites
double Plat[3][3] = {{ 1.000000, 0.000000,  0.000000},
                     { 0.000000, 1.000000,  0.000000},
                     { 0.000000, 0.000000,  1.000000}};
 

double sites[1][3] = {{ 0.00000,  0.00000,  0.00000}};


double PlatINV[3][3];

// ==== INTIALISE HEUN ARRAYS ======================= //
double guassian_vals[3];
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
double M[3];
double H_exch[3];
double H_exch_dash[3];
int counting;
std::normal_distribution<double> distribution(0.0,1.0);
std::random_device device;
std::mt19937 generator(device());
// ================================================ //

Array2D<double> H_thermal;
Array2D<double> Delta_S;
Array3D<double> S_dash_normedx, S_dash_normedy, S_dash_normedz;
Array3D<double> Sx, Sy, Sz;
Array3D<double> Hxi, Hyi, Hzi;
Array3D<double> Jxx, Jyy, Jzz; 
Array3D<fftw_complex> FFTJxx, FFTJyy, FFTJzz;
Array3D<fftw_complex> Hxk, Hyk, Hzk;
Array3D<fftw_complex> FFTSx, FFTSy, FFTSz;
Array3D<fftw_complex> FFTSx1, FFTSy1, FFTSz1;
Array3D<int> Scount;
fftw_plan planFFTSx;
fftw_plan planFFTSy; 
fftw_plan planFFTSz;;
fftw_plan planFFTSx1;
fftw_plan planFFTSy1;
fftw_plan planFFTSz1;
fftw_plan planFFTHx;
fftw_plan planFFTHy;
fftw_plan planFFTHz;

double Sx1d[Nq*Lx*Ly*Lz];
double Sy1d[Nq*Lx*Ly*Lz];
double Sz1d[Nq*Lx*Ly*Lz];
double Hx1d[Nq*Lx*Ly*Lz];
double Hy1d[Nq*Lx*Ly*Lz];
double Hz1d[Nq*Lx*Ly*Lz];
double S_dash_normedx1d[Nq*Lx*Ly*Lz];
double S_dash_normedy1d[Nq*Lx*Ly*Lz];
double S_dash_normedz1d[Nq*Lx*Ly*Lz];

Array4D<double> latticeX;
Array4D<double> latticeY;
Array4D<double> latticeZ;
Array4D<double> LatCount;

// NEIGHBOUR LIST //
std::vector<int> adjncy;
std::vector<double> Jijy_prime;
std::vector<double> Jijz_prime;
std::vector<double> Jijx_prime;
std::vector<int> x_adj;

void CreateLattice(){
    int counter = 0;
    int Cx = 0;
    int Cy = 0;
    int Cz = 0;

    latticeX.resize(Lx,Ly,Lz,Nq);
    latticeY.resize(Lx,Ly,Lz,Nq);
    latticeZ.resize(Lx,Ly,Lz,Nq);
    LatCount.resize(Lx,Ly,Lz,Nq);

    latticeX.IFill(0);
    latticeY.IFill(0);
    latticeZ.IFill(0);
    LatCount.IFill(0);

    for (int x = 0; x < Lx; ++x){ 
        Cy = 0;          
        for (int y = 0; y < Ly; ++y){     
            Cz = 0;     
            for (int z = 0; z < Lz; z++){   
                for (int q = 0; q < Nq; q++){   
                    latticeX(x,y,z,q) = sites[q][0] + Plat[0][0]*Cx + Plat[1][0]*Cy + Plat[2][0]*Cz;
                    latticeY(x,y,z,q) = sites[q][1] + Plat[0][1]*Cx + Plat[1][1]*Cy + Plat[2][1]*Cz;
                    latticeZ(x,y,z,q) = sites[q][2] + Plat[0][2]*Cx + Plat[1][2]*Cy + Plat[2][2]*Cz;
                    LatCount(x,y,z,q) = counter;
                    counter++;
                }
                Cz++;
            }
            Cy++;
        }
        Cx++;
    }
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

    Jxx.resize(xdim, ydim, zdim);
    Jyy.resize(xdim, ydim, zdim);
    Jzz.resize(xdim, ydim, zdim);
    FFTJxx.resize(xdim, ydim, zdimC);
    FFTJyy.resize(xdim, ydim, zdimC);
    FFTJzz.resize(xdim, ydim, zdimC);

    fftw_plan planFFTJxx = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Jxx.ptr(), FFTJxx.ptr(), FFTW_MEASURE);
    fftw_plan planFFTJyy = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Jyy.ptr(), FFTJyy.ptr(), FFTW_MEASURE);
    fftw_plan planFFTJzz = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Jzz.ptr(), FFTJzz.ptr(), FFTW_MEASURE);

    Jxx.IFill(0);
    Jyy.IFill(0);
    Jzz.IFill(0);
    FFTJxx.IFill(0);
    FFTJyy.IFill(0);
    FFTJzz.IFill(0);

    for (int i = 0; i < jxx.size(); i++){
        Jxx(modfunc(xdim, Nx[i]), modfunc(ydim, Ny[i]), modfunc(zdim, Nz[i])) = jxx[i];
        Jyy(modfunc(xdim, Nx[i]), modfunc(ydim, Ny[i]), modfunc(zdim, Nz[i])) = jyy[i];
        Jzz(modfunc(xdim, Nx[i]), modfunc(ydim, Ny[i]), modfunc(zdim, Nz[i])) = jzz[i];
    } 

    fftw_execute(planFFTJxx);
    fftw_execute(planFFTJyy);
    fftw_execute(planFFTJzz);
    fftw_destroy_plan(planFFTJxx);
    fftw_destroy_plan(planFFTJyy);
    fftw_destroy_plan(planFFTJzz);
}

void FFT_initialise(){

    Sx.resize(xdim, ydim, zdim);
    Sy.resize(xdim, ydim, zdim);
    Sz.resize(xdim, ydim, zdim);
    Scount.resize(xdim, ydim, xdim);
    S_dash_normedx.resize(xdim, ydim, zdim);
    S_dash_normedy.resize(xdim, ydim, zdim);
    S_dash_normedz.resize(xdim, ydim, zdim);  
    FFTSx.resize(xdim, ydim, zdimC);
    FFTSy.resize(xdim, ydim, zdimC);
    FFTSz.resize(xdim, ydim, zdimC);
    FFTSx1.resize(xdim, ydim, zdimC);
    FFTSy1.resize(xdim, ydim, zdimC);
    FFTSz1.resize(xdim, ydim, zdimC);
    Hxk.resize(xdim, ydim, zdimC);
    Hyk.resize(xdim, ydim, zdimC);
    Hzk.resize(xdim, ydim, zdimC);
    Hxi.resize(xdim, ydim, zdim);
    Hyi.resize(xdim, ydim, zdim);
    Hzi.resize(xdim, ydim, zdim);

    planFFTSx = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Sx.ptr(), FFTSx.ptr(), FFTW_MEASURE);
    planFFTSy = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Sy.ptr(), FFTSy.ptr(), FFTW_MEASURE);
    planFFTSz = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, Sz.ptr(), FFTSz.ptr(), FFTW_MEASURE);
    planFFTSx1 = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, S_dash_normedx.ptr(), FFTSx1.ptr(), FFTW_MEASURE);
    planFFTSy1 = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, S_dash_normedy.ptr(), FFTSy1.ptr(), FFTW_MEASURE);
    planFFTSz1 = fftw_plan_dft_r2c_3d(xdim, ydim, zdim, S_dash_normedz.ptr(), FFTSz1.ptr(), FFTW_MEASURE);
    planFFTHx = fftw_plan_dft_c2r_3d(xdim, ydim, zdim, Hxk.ptr(), Hxi.ptr(), FFTW_MEASURE);
    planFFTHy = fftw_plan_dft_c2r_3d(xdim, ydim, zdim, Hyk.ptr(), Hyi.ptr(), FFTW_MEASURE);
    planFFTHz = fftw_plan_dft_c2r_3d(xdim, ydim, zdim, Hzk.ptr(), Hzi.ptr(), FFTW_MEASURE);

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

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdimC; z++){
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

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdim; z++){
                Hxi(x,y,z) *= NsitesINV;
                Hyi(x,y,z) *= NsitesINV;
                Hzi(x,y,z) *= NsitesINV;  
            }
        }
    }

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdim; z++){
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
    for (int a = 0; a < Nspins; a++){

        for (int g = 0; g < 3; g++){
            guassian_vals[g] = distribution(generator);
        }

        H_thermal(a,0)  = guassian_vals[0] * Thermal_Fluct;
        H_thermal(a,1)  = guassian_vals[1] * Thermal_Fluct;
        H_thermal(a,2)  = guassian_vals[2] * Thermal_Fluct;

        // Uniaxial anisotropy in Z axis
        H_ani[0]=0;
        H_ani[1]=0;
        H_ani[2]=d_z_prime * Sz1d[a];

        H_new[0] = H_thermal(a,0) + H_app[0] + H_ani[0] + (Hx1d[a] * INVmu_s);
        H_new[1] = H_thermal(a,1) + H_app[1] + H_ani[1] + (Hy1d[a] * INVmu_s);
        H_new[2] = H_thermal(a,2) + H_app[2] + H_ani[2] + (Hz1d[a] * INVmu_s);

        ScrossP[0] = Sx1d[a];
        ScrossP[1] = Sy1d[a];
        ScrossP[2] = Sz1d[a];

        CrossP(ScrossP, H_new, CrossP1);
        CrossP(ScrossP, CrossP1, CrossP2);

        Delta_S(a,0)  = -lambdaPrime * (CrossP1[0] + lambda * CrossP2[0]);
        Delta_S(a,1)  = -lambdaPrime * (CrossP1[1] + lambda * CrossP2[1]);
        Delta_S(a,2)  = -lambdaPrime * (CrossP1[2] + lambda * CrossP2[2]);

        S_dash[0] = Sx1d[a] + (Delta_S(a,0)  * dtau);
        S_dash[1] = Sy1d[a] + (Delta_S(a,1)  * dtau);
        S_dash[2] = Sz1d[a] + (Delta_S(a,2)  * dtau);

        invmag = 1 / sqrt(S_dash[0] * S_dash[0] + S_dash[1] * S_dash[1] + S_dash[2] * S_dash[2]);
        S_dash_normedx1d[a] = invmag * S_dash[0];
        S_dash_normedy1d[a] = invmag * S_dash[1];
        S_dash_normedz1d[a] = invmag * S_dash[2];
    }
}

void FFT_Heun2(double Thermal_Fluct){
    for (int a = 0; a < Nspins; a++){

        H_ani_dash[0]=0;
        H_ani_dash[1]=0;
        H_ani_dash[2]=d_z_prime * S_dash_normedz1d[a];

        H_new_dash[0] = H_thermal(a,0)  + H_app[0] + H_ani_dash[0] + (Hx1d[a] * INVmu_s);
        H_new_dash[1] = H_thermal(a,1)  + H_app[1] + H_ani_dash[1] + (Hy1d[a] * INVmu_s);
        H_new_dash[2] = H_thermal(a,2)  + H_app[2] + H_ani_dash[2] + (Hz1d[a] * INVmu_s);

        // Calculate Corrector and Normalise
        ScrossP[0] = S_dash_normedx1d[a];
        ScrossP[1] = S_dash_normedy1d[a];
        ScrossP[2] = S_dash_normedz1d[a];

        CrossP(ScrossP, H_new_dash, CrossP3);
        CrossP(ScrossP, CrossP3, CrossP4);

        Delta_S_dash[0] = -lambdaPrime * (CrossP3[0] + lambda * CrossP4[0]);
        Delta_S_dash[1] = -lambdaPrime * (CrossP3[1] + lambda * CrossP4[1]);
        Delta_S_dash[2] = -lambdaPrime * (CrossP3[2] + lambda * CrossP4[2]); 

        S_new[0] = Sx1d[a] + half_dtau * (Delta_S(a,0)  + Delta_S_dash[0]);
        S_new[1] = Sy1d[a] + half_dtau * (Delta_S(a,1)  + Delta_S_dash[1]);
        S_new[2] = Sz1d[a] + half_dtau * (Delta_S(a,2)  + Delta_S_dash[2]);

        invmag1 = 1 / sqrt(S_new[0] * S_new[0] + S_new[1] * S_new[1] + S_new[2] * S_new[2]);
        Sx1d[a] = invmag1 * S_new[0];
        Sy1d[a] = invmag1 * S_new[1];
        Sz1d[a] = invmag1 * S_new[2];

        M[0] += Sx1d[a];
        M[1] += Sy1d[a];
        M[2] += Sz1d[a]; 
        
    }

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdim; z++){
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
    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdim; z++){
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

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdimC; z++){
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

    for (int x = 0; x < xdim; x++){
        for (int y = 0; y < ydim; y++){
            for (int z = 0; z < zdim; z++){
                Hxi(x,y,z) *= NsitesINV;
                Hyi(x,y,z) *= NsitesINV;
                Hzi(x,y,z) *= NsitesINV;

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

void CSR_InteractionMatrix(){

    std::ifstream input("/Users/Hirst/Documents/PhD/LLG_code/SimpleCrystal_3D/input/SC_test1.dat");
    double a, b, c, d, e, f;
    std::vector<double> Nx;
    std::vector<double> Ny;
    std::vector<double> Nz;
    std::vector<double> Jij;
    std::vector<int> ib;
    std::vector<int> jb;
    int zdimC = static_cast<double>(zdim)/2 + 1;
    double sum;


    Inverse3x3(Plat, PlatINV);

    while (input >> a >> b >> c >> d >> e >> f)
    {
        // if (std::abs(f) < 0.01) { // Jij less than 1% of max value

        // }
        // else {
            ib.push_back(a);
            jb.push_back(b);
            Nx.push_back(c);
            Ny.push_back(d);
            Nz.push_back(e);
            Jij.push_back(f);
            sum += std::abs(f);
        // }
    }

    std::cout << "Sum of Jij = " << sum << std::endl;
    input.close();

    int length = Jij.size();
    int NxP[length];
    int NyP[length];
    int NzP[length];

    // Convert to unit cell vector
    for (int i = 0; i < length; i++){

        NxP[i] = ( Nx[i] + sites[ib[i]-3][0] - sites[jb[i]-3][0] ) * PlatINV[0][0] + ( Ny[i] + sites[ib[i]-3][1] - sites[jb[i]-3][1] ) * PlatINV[1][0] + ( Nz[i] + sites[ib[i]-3][2] - sites[jb[i]-3][2] ) * PlatINV[2][0]; 
        NyP[i] = ( Nx[i] + sites[ib[i]-3][0] - sites[jb[i]-3][0] ) * PlatINV[0][1] + ( Ny[i] + sites[ib[i]-3][1] - sites[jb[i]-3][1] ) * PlatINV[1][1] + ( Nz[i] + sites[ib[i]-3][2] - sites[jb[i]-3][2] ) * PlatINV[2][1]; 
        NzP[i] = ( Nx[i] + sites[ib[i]-3][0] - sites[jb[i]-3][0] ) * PlatINV[0][2] + ( Ny[i] + sites[ib[i]-3][1] - sites[jb[i]-3][1] ) * PlatINV[1][2] + ( Nz[i] + sites[ib[i]-3][2] - sites[jb[i]-3][2] ) * PlatINV[2][2]; 

    }

    int xval;
    int yval;
    int zval;
    int qval;
    int adjcounter = 0;

    x_adj.push_back(0);

    // ========== Neighbour List =========== //

    for (int x = 0; x < Lx; ++x){                // Depth
        for (int y = 0; y < Ly; ++y){            // Row
            for (int z = 0; z < Lz; ++z){        // Column
                for (int q = 0; q < Nq; ++q){    // Unit Cell
                    for (int i = 0; i < length; i++){

                        if (q == ib[i]-1){
                            xval = modfunc(Lx, NxP[i] + x);
                            yval = modfunc(Ly, NyP[i] + y); 
                            zval = modfunc(Lz, NzP[i] + z);
                            qval = jb[i]-1;
                            adjncy.push_back(LatCount(xval, yval, zval, qval));
                            adjcounter++;

                            // std::cout << xval<< " " << yval << " " << zval << " " << Jij[i] << std::endl;
                            // std::cout << latticeX(xval,yval,zval,qval) << "\t" << latticeY(xval,yval,zval,qval) << "\t" << latticeZ(xval,yval,zval,qval) << "\t" << Jij[i] << std::endl;

                            // Jijx_prime.push_back( ( Jij[i] * 2.179872e-21) / mu_s);
                            // Jijy_prime.push_back( ( Jij[i] * 2.179872e-21) / mu_s);
                            // Jijz_prime.push_back( ( Jij[i] * 2.179872e-21) / mu_s);

                            Jijx_prime.push_back( ( Jij[i]) / mu_s);
                            Jijy_prime.push_back( ( Jij[i]) / mu_s);
                            Jijz_prime.push_back( ( Jij[i]) / mu_s);

                        }
                    }
                    x_adj.push_back(adjcounter);
                }
            }
        }
    }    
}

void CSR_Heun(double Thermal_Fluct){

        for (int a = 0; a < Nspins; a++){

        for (int k = 0; k < 3; k++){
            guassian_vals[k] = distribution(generator);
        }

        H_thermal(a,0) = guassian_vals[0] * Thermal_Fluct;
        H_thermal(a,1) = guassian_vals[1] * Thermal_Fluct;
        H_thermal(a,2) = guassian_vals[2] * Thermal_Fluct;

        // Uniaxial anisotropy in Z axis
        H_ani[0]=0;
        H_ani[1]=0;
        H_ani[2]=d_z_prime * Sz1d[a];

        // Exchange interaction
        H_exch[0] = 0;
        H_exch[1] = 0;
        H_exch[2] = 0;
        counting = 0;

        for (int b = x_adj[a]; b < x_adj[a+1]; b++){
            H_exch[0] += Jijx_prime[counting] * (Sx1d[adjncy[b]]);
            H_exch[1] += Jijy_prime[counting] * (Sy1d[adjncy[b]]);
            H_exch[2] += Jijz_prime[counting] * (Sz1d[adjncy[b]]);
            counting++;
        }

        H_new[0] = H_thermal(a,0) + H_app[0] + H_ani[0] + H_exch[0];
        H_new[1] = H_thermal(a,1) + H_app[1] + H_ani[1] + H_exch[1];
        H_new[2] = H_thermal(a,2) + H_app[2] + H_ani[2] + H_exch[2];

        ScrossP[0] = Sx1d[a];
        ScrossP[1] = Sy1d[a];
        ScrossP[2] = Sz1d[a];

        CrossP(ScrossP, H_new, CrossP1);
        CrossP(ScrossP, CrossP1, CrossP2);

        Delta_S(a,0) = -lambdaPrime * (CrossP1[0] + lambda * CrossP2[0]);
        Delta_S(a,1) = -lambdaPrime * (CrossP1[1] + lambda * CrossP2[1]);
        Delta_S(a,2) = -lambdaPrime * (CrossP1[2] + lambda * CrossP2[2]);

        S_dash[0] = Sx1d[a] + (Delta_S(a,0) * dtau);
        S_dash[1] = Sy1d[a] + (Delta_S(a,1) * dtau);
        S_dash[2] = Sz1d[a] + (Delta_S(a,2) * dtau);

        invmag = 1 / sqrt(S_dash[0] * S_dash[0] + S_dash[1] * S_dash[1] + S_dash[2] * S_dash[2]);
        S_dash_normedx1d[a] = invmag * S_dash[0];
        S_dash_normedy1d[a] = invmag * S_dash[1];
        S_dash_normedz1d[a] = invmag * S_dash[2];   
    }

    for (int a = 0; a < Nspins; a++){

        H_ani_dash[0]=0;
        H_ani_dash[1]=0;
        H_ani_dash[2]=d_z_prime * S_dash_normedz1d[a];

        H_exch_dash[0] = 0;
        H_exch_dash[1] = 0;
        H_exch_dash[2] = 0;
        counting = 0;

        // Exchange interaction prime
        for (int b = x_adj[a]; b < x_adj[a+1]; b++){

            H_exch_dash[0] += Jijx_prime[counting] * (S_dash_normedx1d[adjncy[b]]);
            H_exch_dash[1] += Jijy_prime[counting] * (S_dash_normedy1d[adjncy[b]]);
            H_exch_dash[2] += Jijz_prime[counting] * (S_dash_normedz1d[adjncy[b]]);
            counting++;
        }

        H_new_dash[0] = H_thermal(a,0) + H_app[0] + H_ani_dash[0] + H_exch_dash[0];
        H_new_dash[1] = H_thermal(a,1) + H_app[1] + H_ani_dash[1] + H_exch_dash[1];
        H_new_dash[2] = H_thermal(a,2) + H_app[2] + H_ani_dash[2] + H_exch_dash[2];

        // Calculate Corrector and Normalise

        ScrossP[0] = S_dash_normedx1d[a];
        ScrossP[1] = S_dash_normedy1d[a];
        ScrossP[2] = S_dash_normedz1d[a];

        CrossP(ScrossP, H_new_dash, CrossP3);
        CrossP(ScrossP, CrossP3, CrossP4);

        Delta_S_dash[0] = -lambdaPrime * (CrossP3[0] + lambda * CrossP4[0]);
        Delta_S_dash[1] = -lambdaPrime * (CrossP3[1] + lambda * CrossP4[1]);
        Delta_S_dash[2] = -lambdaPrime * (CrossP3[2] + lambda * CrossP4[2]); 

        S_new[0] = Sx1d[a] + half_dtau * (Delta_S(a,0) + Delta_S_dash[0]);
        S_new[1] = Sy1d[a] + half_dtau * (Delta_S(a,1) + Delta_S_dash[1]);
        S_new[2] = Sz1d[a] + half_dtau * (Delta_S(a,2) + Delta_S_dash[2]);

        invmag1 = 1 / sqrt(S_new[0] * S_new[0] + S_new[1] * S_new[1] + S_new[2] * S_new[2]);
        Sx1d[a] = invmag1 * S_new[0];
        Sy1d[a] = invmag1 * S_new[1];
        Sz1d[a] = invmag1 * S_new[2];

        M[0] += Sx1d[a];
        M[1] += Sy1d[a];
        M[2] += Sz1d[a]; 
    }

}


int main(int argc, char* argv[]){

    std::cout << __LINE__ << std::endl;

    CreateLattice();
    CSR_InteractionMatrix();

    // ======= Temperature ==================================================================================== //
    const double Temp = atof(argv[1]);
    const double thermal_fluct = thermal_const * sqrt(Temp);
    // ========================================================================================================= //

    H_thermal.resize(Nspins, 3);
    Delta_S.resize(Nspins, 3);
    H_thermal.IFill(0);
    Delta_S.IFill(0);

    for (int a = 0; a < Nspins; a++){ 
        Sx1d[a] = 0.0;
        Sy1d[a] = 0.0;
        Sz1d[a] = 1.0;
    }
    // ================================================================================================== //
        
    // ======= OPEN OUTPUT FILE ======================================================================== //
    std::stringstream sstr;
    sstr << "mag_tsteps_" << Nt << "_T_" << Temp << ".txt";
    std::ofstream myfile;
    myfile.open(sstr.str());
    // ================================================================================================= //

    // ========== Time + Temp Variables ================================================================= //
    double t = 0;
    double tau = 0;

    clock_t begin, end;
    double time_spent;
    begin = clock();
    // ================================================================================================== //


    // ========== LOOP THROUGH TIMESTEPS ================================================================ //
    for (int i = 0; i < Nt; i++){
        
        t = t + dt;
        tau = tau + dtau;

        M[0] = 0;
        M[1] = 0;
        M[2] = 0;

        CSR_Heun(thermal_fluct);


        double Mmag = sqrt(M[0] * M[0] + M[1] * M[1] + M[2] * M[2]);
        double MdivMs = Mmag / (Nspins);
 
        // Output to file
        std::cout << i << " ";
        std::cout << M[0] / Nspins << " " << M[1] / Nspins << " " << M[2] / Nspins << " " << MdivMs << "\n";  
    }
    // ==================================================================================================== //

    myfile << std::flush;

    end = clock();
    std::cout << std::setprecision(10) << "Time in seconds = " << (double)(end - begin) / CLOCKS_PER_SEC << std::endl;


    // Destroy Plans
    // FFTdestroyplans();

    // CLOSE FILE
    myfile.close();

    return 0;
}