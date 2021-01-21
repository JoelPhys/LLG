#ifndef __GEOM_H__
#define __GEOM_H__
#include "array4d.h"

    namespace geom {

        extern Array4D<double> latticeX;
        extern Array4D<double> latticeY;
        extern Array4D<double> latticeZ;
        extern Array4D<int> LatCount;

        void CreateLattice();                 
    }

#endif