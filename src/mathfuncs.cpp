#include <cmath>
#include <fstream>
#include "../inc/mathfuncs.h"
#include "../inc/params1.h"
#include "../inc/util.h"
#include "../inc/NeighbourList.h"
#include "../inc/array.h"
#include "../inc/array2d.h"

int modfunc(int L, int x){
        int y = ((x % L) + L) % L;
        return y;
    }

void CrossP(double vect_A[], double vect_B[], double crossP[]){
    crossP[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
    crossP[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
    crossP[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}

void DotP(double vect_A1[], double vect_B1[], double dotP){
    dotP = vect_A1[0] * vect_B1[0] + vect_A1[1] * vect_B1[1] + vect_A1[2] * vect_B1[2]; 
}

void Inverse3x3(double m[3][3], double minv[3][3]){
    double det = m[0][ 0] * (m[1][ 1] * m[2][ 2] - m[2][ 1] * m[1][ 2]) -
                    m[0][ 1] * (m[1][ 0] * m[2][ 2] - m[1][ 2] * m[2][ 0]) +
                    m[0][ 2] * (m[1][ 0] * m[2][ 1] - m[1][ 1] * m[2][ 0]);

    double invdet = 1 / det;

    minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
    minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
    minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

}


//rotation matrix
double xu[3] = {1.0,0.0,0.0};
double C;
double D[3];
double X[3][3];
double Xsqr[3][3];
double R[3][3];
double Mnewb[3];
double norm;
double pref;
double Mnew[3];
double Mu[3];
Array2D<double> Id;


void IdentityMatrix(){

    Id.resize(3,3);
    Id.IFill(0);
    Id(0,0) = 1;
    Id(1,1) = 1;
    Id(2,2) = 1;

}

void Rotation(){

            //normalise M
            Mnew[0] = util::M(0,0) / params::NmomentsSubLat;
            Mnew[1] = util::M(0,1) / params::NmomentsSubLat;
            Mnew[2] = util::M(0,2) / params::NmomentsSubLat;
            norm = sqrt(Mnew[0] * Mnew[0] + Mnew[1] * Mnew[1] + Mnew[2] * Mnew[2]);

            Mu[0] = Mnew[0] / norm;
            Mu[1] = Mnew[1] / norm;
            Mu[2] = Mnew[2] / norm;

            C = Mu[0] * xu[0] + Mu[1] * xu[1] + Mu[2] * xu[2]; 
            CrossP(Mu, xu, D);

            X[0][0] = 0;
            X[0][1] = -1 * D[2];
            X[0][2] = D[1];
            X[1][0] = D[2];
            X[1][1] = 0;
            X[1][2] = -1 * D[0];
            X[2][0] = -1 * D[1];
            X[2][1] = D[0];
            X[2][2] = 0;
            norm = sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]); 
            pref = (1 - C) / (norm * norm);

            for (int w = 0; w < 3; w++){
                for (int e = 0; e < 3; e++){
                    Xsqr[w][e] = 0;
                }
            }

            for (int w = 0; w < 3; w++){
                for (int e = 0; e < 3; e++){
                    for (int r = 0; r < 3; r++){
                        Xsqr[w][e] += X[w][r] * X[r][e];
                    }
                }
            }

            for (int w = 0; w < 3; w++){
                for (int e = 0; e < 3; e++){
                    R[w][e] = Id(w,e) + X[w][e] + ( (1 - C) / (norm*norm) ) * Xsqr[w][e];
                }
            }

            for (int w = 0; w < 3; w++){
                Mnewb[w] = 0;
            }

            for (int a = 0; a < params::Nspins; a++){     
                    neigh::Sx1d[a] = R[0][0] * neigh::Sx1d[a] + R[0][1] * neigh::Sy1d[a] + R[0][2] * neigh::Sz1d[a];
                    neigh::Sy1d[a] = R[1][0] * neigh::Sx1d[a] + R[1][1] * neigh::Sy1d[a] + R[1][2] * neigh::Sz1d[a];
                    neigh::Sz1d[a] = R[2][0] * neigh::Sx1d[a] + R[2][1] * neigh::Sy1d[a] + R[2][2] * neigh::Sz1d[a];
            }
}

