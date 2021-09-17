// cpp header files
#include <cmath>
#include <iostream>
#include <fftw3.h>

// my header files
#include "../inc/array3d.h"

// Prevent division by 0
double eps = 1e-18;

// demag tensor
Array3D<double> N11; 
Array3D<double> N22; 
Array3D<double> N33; 
Array3D<double> N12; 
Array3D<double> N13; 
Array3D<double> N23; 

// demag tensor in fourier space
Array3D<fftw_complex> FFTN11;
Array3D<fftw_complex> FFTN22;
Array3D<fftw_complex> FFTN33;
Array3D<fftw_complex> FFTN12;
Array3D<fftw_complex> FFTN13;
Array3D<fftw_complex> FFTN23;


int cellx = 5;
int celly = 5;
int cellz = 5;

double dr1 = 5e-9;
double dr2 = 5e-9;
double dr3 = 5e-9;


// fft plans
fftw_plan N11plan;
fftw_plan N22plan;
fftw_plan N33plan;
fftw_plan N12plan;
fftw_plan N13plan;
fftw_plan N23plan;

// Modulo function
int modfunc(int L, int x){
	int y = ((x % L) + L) % L;
	return y;
}

void plans(){

    FFTN11.resize(2*cellx-1,2*celly-1,2*cellz-1);
    FFTN22.resize(2*cellx-1,2*celly-1,2*cellz-1);
    FFTN33.resize(2*cellx-1,2*celly-1,2*cellz-1);
    FFTN12.resize(2*cellx-1,2*celly-1,2*cellz-1);
    FFTN13.resize(2*cellx-1,2*celly-1,2*cellz-1);
    FFTN23.resize(2*cellx-1,2*celly-1,2*cellz-1);

    std::cout << __LINE__ << std::endl;

    N11plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N11.ptr(), FFTN11.ptr(), FFTW_MEASURE);
    N22plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N22.ptr(), FFTN22.ptr(), FFTW_MEASURE);
    N33plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N33.ptr(), FFTN33.ptr(), FFTW_MEASURE);
    N12plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N12.ptr(), FFTN12.ptr(), FFTW_MEASURE);
    N13plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N13.ptr(), FFTN13.ptr(), FFTW_MEASURE);
    N23plan = fftw_plan_dft_r2c_3d(2*cellx-1,2*celly-1,2*cellz-1, N23.ptr(), FFTN23.ptr(), FFTW_MEASURE);

    std::cout << __LINE__ << std::endl;


    FFTN11.IFill(0);
    FFTN22.IFill(0);
    FFTN33.IFill(0);
    FFTN12.IFill(0);
    FFTN13.IFill(0);
    FFTN23.IFill(0);
}

long double f(double x1,double y1,double z1){

    long double x = std::abs(x1);
    long double y = std::abs(y1);
    long double z = std::abs(z1);

    long double fres = 
        0.5*y*(z*z-x*x)*asinh(y/(sqrt(x*x+z*z)+eps)) 
        + 0.5*z*(y*y-x*x)*asinh(z/(sqrt(x*x+y*y)+eps)) 
        - x*y*z*atan(y*z/(x*sqrt(x*z+y*y+z*z)+eps)) 
        + 1.0/6.0*(2*x*x-y*y-z*z)*sqrt(x*x+y*y+z*z);
    
    return fres;

}

long double g(long double x,long double y,long double z1){

    long double z = std::abs(z1);

    long double gres = x*y*z*asinh(z/(sqrt(x*x+y*y)+eps))                 
        + y/6.0*(3.0*z*z-y*y)*asinh(x/(sqrt(y*y+z*z)+eps)) 
        + x/6.0*(3.0*z*z-x*x)*asinh(y/(sqrt(x*x+z*z)+eps)) 
        - z*z*z/6.0*atan(x*y/(z*sqrt(x*x+y*y+z*z)+eps))    
        - z*y*y*0.5*atan(x*z/(y*sqrt(x*x+y*y+z*z)+eps))    
        - z*x*x*0.5*atan(y*z/(x*sqrt(x*x+y*y+z*z)+eps))    
        - x*y*sqrt(x*x+y*y+z*z)/3.0;

    return gres;
}



int main(){

    N11.resize(2*cellx-1,2*celly-1,2*cellz-1);
    N22.resize(2*cellx-1,2*celly-1,2*cellz-1);
    N33.resize(2*cellx-1,2*celly-1,2*cellz-1);
    N12.resize(2*cellx-1,2*celly-1,2*cellz-1);
    N13.resize(2*cellx-1,2*celly-1,2*cellz-1);
    N23.resize(2*cellx-1,2*celly-1,2*cellz-1);

    plans();

    N11.IFill(0.0);
    N22.IFill(0.0);
    N33.IFill(0.0);
    N12.IFill(0.0);
    N13.IFill(0.0);
    N23.IFill(0.0);

    int x,y,z;
    double exponent;

    for (int i1 = 0; i1 < cellx; i1++){
        for (int i2 = 0; i2 < celly; i2++){
            for (int i3 = 0; i3 < cellz; i3++){
                for (int j1 = 0; j1 < cellx; j1++){
                    for (int j2 = 0; j2 < celly; j2++){
                        for (int j3 = 0; j3 < cellz; j3++){


                            x = modfunc(cellx,i1-j1);
                            y = modfunc(celly,i2-j2);
                            z = modfunc(cellz,i3-j3);

                            N11(x,y,z) = 0.0;
                            N22(x,y,z) = 0.0;
                            N33(x,y,z) = 0.0;
                            N12(x,y,z) = 0.0;
                            N13(x,y,z) = 0.0;
                            N23(x,y,z) = 0.0;


                            for (int k1 = 0; k1 < 2; k1++){
                                for (int k2 = 0; k2 < 2; k2++){
                                    for (int k3 = 0; k3 < 2; k3++){
                                        for (int l1 = 0; l1 < 2; l1++){
                                            for (int l2 = 0; l2 < 2; l2++){
                                                for (int l3 = 0; l3 < 2; l3++){

                                                    exponent = k1+k2+k3+l1+l2+l3;

                                                    N11(x,y,z) += pow(-1,exponent) * f((i1-j1+k1-l1)*dr1, (i2-j2+k2-l2) * dr2, (i3-j3+k3-l3) * dr3);
                                                    N22(x,y,z) += pow(-1,exponent) * f((i2-j2+k1-l1)*dr1, (i3-j3+k2-l2) * dr2, (i1-j1+k3-l3) * dr3);
                                                    N33(x,y,z) += pow(-1,exponent) * f((i3-j3+k1-l1)*dr1, (i1-j1+k2-l2) * dr2, (i2-j2+k3-l3) * dr3);
                                                    N12(x,y,z) += pow(-1,exponent) * g((i1-i2+k1-l1)*dr1, (i2-j2+k2-l2) * dr2, (i3-j3+k3-l3) * dr3);
                                                    N13(x,y,z) += pow(-1,exponent) * g((i1-j3+k1-l1)*dr1, (i3-j3+k2-l2) * dr2, (i2-j2+k3-l3) * dr3);
                                                    N23(x,y,z) += pow(-1,exponent) * g((i2-i3+k1-l1)*dr1, (i3-j3+k2-l2) * dr2, (i1-j1+k3-l3) * dr3);
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            N11(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3);
                            N22(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3);
                            N33(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3);
                            N12(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3);
                            N13(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3);
                            N23(x,y,z) *= 1.0/(4*M_PI*dr1*dr2*dr3); 

                            std::cout << i1 << " " << i2 << " " << i3 << " " << j1 << " " << j2 << " " << j3 << " ";
                            std::cout << x << " " << y << " " << z << " ";
                            std::cout << N11(x,y,z) << " " << N22(x,y,z) << " " << N33(x,y,z) << " ";
                            std::cout << N12(x,y,z) << " " << N13(x,y,z) << " " << N23(x,y,z) << "\n";

                        }
                    }
                }
    
            }    
        }   
    }
    std::cout << __LINE__ << std::endl;
    fftw_execute(N11plan);
    std::cout << __LINE__ << std::endl;
    fftw_execute(N22plan);
    fftw_execute(N33plan);
    fftw_execute(N12plan);
    fftw_execute(N13plan);
    fftw_execute(N23plan);

}