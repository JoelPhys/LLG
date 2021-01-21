#include <cmath>
#include <fstream>
#include "../inc/mathfuncs.h"
#include "../inc/array.h"

int modfunc(int L, int x){
        int y = ((x % L) + L) % L;
        return y;
    }

void CrossP(double vect_A[], double vect_B[], double crossP[]){
    crossP[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
    crossP[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
    crossP[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}

void Inverse3x3(double m[3][3], double minv[3][3]){
    double det = m[0][ 0] * (m[1][ 1] * m[2][ 2] - m[2][ 1] * m[1][ 2]) -
                    m[0][ 1] * (m[1][ 0] * m[2][ 2] - m[1][ 2] * m[2][ 0]) +
                    m[0][ 2] * (m[1][ 0] * m[2][ 1] - m[1][ 1] * m[2][ 0]);

    double invdet = 1 / det;

    minv[0][ 0] = (m[1][ 1] * m[2][ 2] - m[2][ 1] * m[1][ 2]) * invdet;
    minv[0][ 1] = (m[0][ 2] * m[2][ 1] - m[0][ 1] * m[2][ 2]) * invdet;
    minv[0][ 2] = (m[0][ 1] * m[1][ 2] - m[0][ 2] * m[1][ 1]) * invdet;
    minv[1][ 0] = (m[1][ 2] * m[2][ 0] - m[1][ 0] * m[2][ 2]) * invdet;
    minv[1][ 1] = (m[0][ 0] * m[2][ 2] - m[0][ 2] * m[2][ 0]) * invdet;
    minv[1][ 2] = (m[1][ 0] * m[0][ 2] - m[0][ 0] * m[1][ 2]) * invdet;
    minv[2][ 0] = (m[1][ 0] * m[2][ 1] - m[2][ 0] * m[1][ 1]) * invdet;
    minv[2][ 1] = (m[2][ 0] * m[0][ 1] - m[0][ 0] * m[2][ 1]) * invdet;
    minv[2][ 2] = (m[0][ 0] * m[1][ 1] - m[1][ 0] * m[0][ 1]) * invdet;

}