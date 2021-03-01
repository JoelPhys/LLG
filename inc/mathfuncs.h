#ifndef __MATHFUNCS_H__
#define __MATHFUNCS_H__

#include "../inc/array2d.h"

extern Array2D<double> R;

int modfunc(int L, int x);
void CrossP(double vect_A[], double vect_B[], double crossP[]);
void DotP(double vect_A1[], double vect_B1[], double dotP);
void Inverse3x3(double m[][3], double minv[][3]);
void IdentityMatrix();
void Rotation();

#endif
