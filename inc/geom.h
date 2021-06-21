#ifndef __GEOM_H__
#define __GEOM_H__
#include "array4d.h"
#include "array3d.h"

    namespace geom {

        extern Array4D<double> latticeX;
        extern Array4D<double> latticeY;
        extern Array4D<double> latticeZ;
        extern Array3D<double> Sx;
        extern Array3D<double> Sy; 
        extern Array3D<double> Sz;
        extern Array3D<int> Scount;
        extern Array4D<int> LatCount;

        extern int Ix, Iy, Iz, IzC;
        extern int latXsize, latYsize, latZsize, latZsizeS;

        void CreateLattice();     
        void CountDistinct(); 
        void CreateIntLattice();           
        void InitSpins();
    }

#endif