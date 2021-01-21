
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>
#include <iomanip>
#include "inc/array.h"
#include "inc/array2d.h"
#include "inc/array3d.h"
#include "inc/array4d.h"


double J0;

double gamma1 = 1.76e11;
double mu_b = 9.2740e-24;
double mu_s = 2 * mu_b;
double a1 = 2.87e-10; //iron
double Nk = 10000;

std::complex<double> I(0,1);
std::complex<double> Jkx;
std::complex<double> Jky;
std::complex<double> Jkz;
std::complex<double> NxC;
std::complex<double> NyC;
std::complex<double> NzC;
std::complex<double> jC;
std::complex<double> k;

int main(){

    // std::ifstream input("../Mn2Au/mn2au-lmto-jij_new/rsj-gga-af-p3.mn2au");
    std::ifstream input("SC_test1.dat");
    double a, b, c, d, e, f;
    std::vector<double> Nx;
    std::vector<double> Ny;
    std::vector<double> Nz;
    std::vector<double> Jij;
    std::vector<int> ib;
    std::vector<int> jb;


    while (input >> a >> b >> c >> d >> e >> f)
    {
        if (a == 1) { 
            ib.push_back(a);
            jb.push_back(b);
            Nx.push_back(c);
            Ny.push_back(d);
            Nz.push_back(e);
            Jij.push_back(f);
        }
    }
    input.close();
    
    int length = Jij.size();

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
    // ======================================================================================================

    std::stringstream sstr;
    sstr << "outputfile.txt";
    std::ofstream myfile;
    myfile.open(sstr.str());

    // calculate J(k) =====================================================================================
    for (int j= 0; j < Nk; j++){


        k = std::complex<double>(j * 1e6,0);
        Jkx = std::complex<double>(0,0);
        Jky = std::complex<double>(0,0);
        Jkz = std::complex<double>(0,0);

        for (int i = 0; i < length; i++){

            jC = std::complex<double>(Jij[i],0);
            NzC = std::complex<double>(Nx[i] * a1,0);
            NyC = std::complex<double>(Ny[i] * a1,0);
            NzC = std::complex<double>(Nz[i] * a1,0);
            Jkx += jC * std::exp(I * k * NxC);
            Jky += jC * std::exp(I * k * NyC);
            Jkz += jC * std::exp(I * k * NzC);

        }

        myfile << k.real() << " ";
        // myfile << (gamma1 / mu_s) * (J0 - Jkx.real()) << " ";   
        // myfile << (gamma1 / mu_s) * (J0 - Jky.real()) << " ";
        // myfile << (gamma1 / mu_s) * (J0 - Jkz.real()) << "\n";
        myfile << (gamma1 / mu_s) * sqrt(J0*J0 - Jkz.real() * Jkz.real()) << "\n";

        // Plot along brillion zone path


    }
    // ======================================================================================================

    myfile.close();


    return 0;
}