#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>
#include <iomanip>
#include "../inc/mathfuncs.h"
#include "../inc/fftw3.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/array3d.h"
#include "../inc/array4d.h"


double J0;

double gamma1 = 1.76e11;
double mu_b = 9.2740e-24;
double mu_s = 2.22  * mu_b;
double a1 = 0.287e-9; //iron
double Nk = 1000;
double hbar = 6.583e-16 / 0.001; // eV s 
int Lx = 100;
int Ly = 100;
int Lz = 100;
int dx = 2;
int dy = 2;
int dz = 2;
int IzC = (Lz * dz)/2 + 1;
int zdimc = (Lz/2) + 1;

// testing
fftw_plan planFFTJzz;
Array3D<double> Jzz;
Array3D<fftw_complex> FFTJzz;

double M_S = (mu_s * 2)/ (a1 * a1 * a1);

double normsize = Lx * Ly * Lz;

std::complex<double> I(0,1);
std::complex<double> Jkx;
std::complex<double> Jky;
std::complex<double> Jkz;
std::complex<double> NxC;
std::complex<double> NyC;
std::complex<double> NzC;
std::complex<double> jC;
std::complex<double> kx;
std::complex<double> ky;
std::complex<double> kz;

double xval, yval;

int main(int argc, char* argv[]){
    

    const double Magnetisation = (atof(argv[1])); 

    Jzz.resize(Lx* dx, Ly * dy , Lz * dz);
    Jzz.IFill(0);

    FFTJzz.resize(Lx*dx, Ly*dy, IzC);
    FFTJzz.IFill(0);

    planFFTJzz = fftw_plan_dft_r2c_3d(Lx* dx, Ly * dy , Lz * dz, Jzz.ptr(), FFTJzz.ptr(), FFTW_MEASURE);
    
    // std::ifstream input("../Mn2Au/mn2au-lmto-jij_new/rsj-gga-af-p3.mn2au");
    // std::ifstream input("../../SimpleCrystal_3D/input/SC_test1.dat");
    std::ifstream input("../../Fe/BCC_testlswt.dat");
    double a, b, c, d, e, f;
    std::vector<double> Nx;
    std::vector<double> Ny;
    std::vector<double> Nz;
    std::vector<double> Jij;
    std::vector<int> ib;
    std::vector<int> jb;


    if (!input.is_open()) {
        std::cout << "Error opening file" << std::endl;
    }

    while (input >> a >> b >> c >> d >> e >> f)
    {
        // if (a == 1) { 
            ib.push_back(a);
            jb.push_back(b);
            Nx.push_back(c);
            Ny.push_back(d);
            Nz.push_back(e);
            Jij.push_back(f);
        // }
    }

    input.close();
    
    int length = Jij.size();

    // FOR JEROMES FILES
    // for (int i = 0; i < length; i++){
    //     if (jb[i] == 3){
    //        Nz[i] =  Nz[i] / 2.5639422;
    //     }
    //     if (jb[i] == 4){
    //         if (Nz[i] = 0.4235633){ Nz[i] = 1.0} 
    //         if (Nz[i] = -2.140379){ Nz[i] = -4.0}
    //         if (Nz[i] = 2.9875055){ Nz[i] = 7.0}
    //         if (Nz[i] = 2.9875055){ Nz[i] = 7.0}
    //     }
    //     if (jb[i] == 5){

            
    //     }
    //     if (jb[i] == 6){

            
    //     }
    // }

    for (int i = 0; i < length; i++){
        Jzz(modfunc(Lx * dx, dx * Nx[i]), modfunc(Ly * dy, dy * Ny[i]), modfunc(Lz * dz, dz * Nz[i])) = Jij[i];
        // Jzz(modfunc(Lx, dx * Nx[i]), modfunc(Ly, dy * Ny[i]), modfunc(Lz, dz * Nz[i])) = Jij[i];
        std::cout << dx * Nx[i] << " ";
        std::cout << dy * Ny[i] << " ";
        std::cout << dz * Nz[i] << std::endl;
    } 

    fftw_execute(planFFTJzz);


    // for (int x = 0; x < Lx* dx; x++){
    //     for (int y = 0; y < Ly * dy; y++){
    //         for (int z = 0; z < IzC; z++){
    //             FFTJzz(x,y,z)[0] /= normsize;
    //             FFTJzz(x,y,z)[1] /= normsize;
    //         }
    //     }
    // }   
    
    // for (int i = 0; i < length; i++){
    //     Jij[i] = Jij[i] *  2.179872e-21 * 2;
    // }

    //  sign convention =====================================================================================
    // for (int i = 0; i < length; i++){
    //     if (ib[i] == jb[i]){
    //         Jij[i] = Jij[i];
    //     }
    //     else if (((ib[i] == 3) && (jb[i] == 5)) || ((ib[i] == 5) && (jb[i] == 3))) {
    //         Jij[i] = Jij[i];
    //     }
    //     else if (((ib[i] == 4) && (jb[i] == 6)) || ((ib[i] == 6) && (jb[i] == 4))) {
    //         Jij[i] = Jij[i];
    //     }
    //     else  {
    //         Jij[i] = -1 * Jij[i];
    //     }   
    // }
    // =====================================================================================================

    //  calculate J(0) =====================================================================================
    for (int i = 0; i < length; i++){
        J0 += Jij[i];
    }

    std::cout << "J0 = " << J0 << std::endl;
    // ======================================================================================================

    std::stringstream sstr;
    sstr << "output/bcc_fm_lswt_" << Magnetisation << ".txt";
    std::ofstream myfile;
    myfile.open(sstr.str());

    // calculate J(k) =====================================================================================
    for (int z = 0; z < IzC; z+=2){    

        std::cout << z << std::endl;    

        yval = Magnetisation * (gamma1 / (dz * 2 * M_PI * mu_s)) * (J0 - FFTJzz(0,0,z)[0]);
        xval = z * ((2 * M_PI) / 0.287e-9)/ IzC;
        myfile << xval << " " << yval << "\n";

    }
    for (int z = 0; z < IzC/2; z+=2){        
        yval = ((Magnetisation * gamma1) / (dz * 2 * M_PI * mu_s)) * (J0 - FFTJzz(0,z,IzC-1-z)[0]);
        xval = z * ((2 * M_PI) / 0.287e-9)/ IzC;
        myfile << xval << " " << yval << "\n";

    }
    for (int z = 0; z < IzC/2; z+=2){        
        yval = ((Magnetisation * gamma1) / ( dz * 2 * M_PI * mu_s)) * (J0 - FFTJzz(0,IzC/2-z,IzC/2-z)[0]);
        xval = z * ((2 * M_PI) / 0.287e-9)/ IzC;
        myfile << xval << " " << yval << "\n";

    }
    for (int z = 0; z < IzC/2; z+=2){        
        yval = ((Magnetisation * gamma1) / ( dz * 2 * M_PI * mu_s)) * (J0 - FFTJzz(z,z,z)[0]);
        xval = z * ((2 * M_PI) / 0.287e-9)/ IzC;
        myfile << xval << " " << yval << "\n";

    }
    // ====================================================================================================

    
    // calculate J(k) =====================================================================================
    // for (int x = 0; x < 1; x++){
    //     for (int y = 0; y < 1; y++){
    //         for (int z = 0; z < Nk; z++){        
    
    //             kx = std::complex<double>(x * (2*M_PI/a1) / Nk,0);
    //             ky = std::complex<double>(y * (2*M_PI/a1) / Nk,0);
    //             kz = std::complex<double>(z * (2*M_PI/a1) / Nk,0);

    //             // kx = std::complex<double>(Nk - (x * (M_PI/a1) / Nk),0);
    //             // ky = std::complex<double>(Nk - (z * (M_PI/a1) / Nk),0);
    //             // kz = std::complex<double>(Nk - (z * (M_PI/a1) / Nk),0);

    //             Jkx = std::complex<double>(0,0);
    //             Jky = std::complex<double>(0,0);
    //             Jkz = std::complex<double>(0,0);

    //             for (int i = 0; i < length; i++){

    //                 jC = std::complex<double>(Jij[i],0);
    //                 NzC = std::complex<double>(Nx[i] * a1,0);
    //                 NyC = std::complex<double>(Ny[i] * a1,0);
    //                 NzC = std::complex<double>(Nz[i] * a1,0);
    //                 Jkz += jC * std::exp(I * two * M_PI * (kx * NxC + ky * NyC + kz * NzC));
    //             }
                
    //             std::cout << J0 - Jkz.real() << std::endl;

    //             val = (gamma1 / (2 * M_PI * mu_s)) * (J0 * J0 - Jkz.real() * Jkz.real() - Jkz.imag() * Jkz.imag());

    //             // std::cout << M_S << std::endl;

    //             myfile << kx.real() << " " << ky.real() << " " << kz.real() << " ";
    //             // myfile << - (J0 - Jkz.real()) * (M_S / (2 * M_PI)) << "\n";
    //             myfile << val << "\n";
    //             // myfile << (gamma1 / mu_s) * sqrt(J0*J0 - Jkz.real() * Jkz.real()) << "\n";

    //         }
    //     }
    // }
    // ======================================================================================================

    myfile.close();

    return 0;
}
