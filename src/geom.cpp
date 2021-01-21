#include "../inc/params1.h"
#include "../inc/array4d.h"
#include "../inc/geom.h"

namespace geom {


    Array4D<double> latticeX;
    Array4D<double> latticeY;
    Array4D<double> latticeZ;
    Array4D<int> LatCount;

    void CreateLattice(){
        int counter = 0;
        int Cx = 0;
        int Cy = 0;
        int Cz = 0;

        latticeX.resize(params::Lx,params::Ly,params::Lz,params::Nq);
        latticeY.resize(params::Lx,params::Ly,params::Lz,params::Nq);
        latticeZ.resize(params::Lx,params::Ly,params::Lz,params::Nq);
        LatCount.resize(params::Lx,params::Ly,params::Lz,params::Nq);

        latticeX.IFill(0);
        latticeY.IFill(0);
        latticeZ.IFill(0);
        LatCount.IFill(0);

        for (int x = 0; x < params::Lx; ++x){ 
            Cy = 0;          
            for (int y = 0; y < params::Ly; ++y){     
                Cz = 0;     
                for (int z = 0; z < params::Lz; z++){   
                    for (int q = 0; q < params::Nq; q++){   
                        latticeX(x,y,z,q) = params::sites[q][0] + params::Plat[0][0]*Cx + params::Plat[1][0]*Cy + params::Plat[2][0]*Cz;
                        latticeY(x,y,z,q) = params::sites[q][1] + params::Plat[0][1]*Cx + params::Plat[1][1]*Cy + params::Plat[2][1]*Cz;
                        latticeZ(x,y,z,q) = params::sites[q][2] + params::Plat[0][2]*Cx + params::Plat[1][2]*Cy + params::Plat[2][2]*Cz;
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


}