#ifndef _MATHFUNCS_H_
#define _MATHFUNCS_H_

#include "../inc/array2d.h"

// global variables
extern Array2D<double> R;

// Functions
void IdentityMatrix();
int modfunc(int L, int x);
void Inverse3x3(double m[][3], double minv[][3]);
void DotP(double vect_A1[], double vect_B1[], double dotP);
void CrossP(double vect_A[], double vect_B[], double crossP[]);

#endif
