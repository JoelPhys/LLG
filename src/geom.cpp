#include "../inc/params1.h"
#include "../inc/array4d.h"
#include "../inc/geom.h"
#include <vector>
#include <iostream>
#include <unordered_set>

namespace geom {


    Array4D<double> latticeX;
    Array4D<double> latticeY;
    Array4D<double> latticeZ;
    Array4D<int> LatCount;
    Array3D<double> Sx;
    Array3D<double> Sy; 
    Array3D<double> Sz;
    Array3D<int> Scount;

    int Ix, Iy, Iz, IzC;

    int latXsize, latYsize, latZsize, latZsizeS;

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

    void CountDistinct(){
        int size = params::sites.size();
        
        for (int c = 0; c < 3; c++){
            int res[3] = {0,0,0};

            std::unordered_set<double> s;



            for (int cc = 0; cc < size; cc++){
                if (s.find(params::sites[cc][c]) == s.end()){
                    s.insert(params::sites[cc][c]);
                    res[c]++;
                }
            }
        }
    }

    void CreateIntLattice(){

        Ix = params::Idx * params::Lx;
        Iy = params::Idy * params::Ly;
        Iz = params::Idz * params::Lz;
        IzC = Iz/2 + 1;

        latXsize = Ix + params::Idx;
        latYsize = Iy + params::Idy;
        latZsize = Iz + params::Idz;
        latZsizeS = latZsize/2 + 1;


        Sx.resize(Ix, Iy, Iz);
        Sy.resize(Ix, Iy, Iz);
        Sz.resize(Ix, Iy, Iz);
        Scount.resize(Ix, Iy, Iz);
        Sx.IFill(0);
        Sy.IFill(0);
        Sz.IFill(0);
        Scount.IFill(0);
        int counter = 0;

        for (int l = 0; l < Ix; l += params::Idx){
            for (int m = 0; m < Iy; m += params::Idy){
                for (int n = 0; n < Iz; n += params::Idz){
                    for (int q = 0; q < params::Nq; q++){

                        // std::cout << l + params::Isites[q][0] << " " << m + params::Isites[q][1] << " " << n + params::Isites[q][2] << std::endl;
                
                        Sx(l + params::Isites[q][0], m + params::Isites[q][1], n + params::Isites[q][2]) = 1;
                        Sy(l + params::Isites[q][0], m + params::Isites[q][1], n + params::Isites[q][2]) = 1;
                        Sz(l + params::Isites[q][0], m + params::Isites[q][1], n + params::Isites[q][2]) = 1;
                        Scount(l + params::Isites[q][0], m + params::Isites[q][1], n + params::Isites[q][2]) = counter;
                        counter++;                    
                    }
                }
            }
        }

        // for (int l = 0; l < 8; l++){
        //     for (int m = 0; m < 24; m++){

        //         std::cout << Sz(1,l,m) <<  " ";
        //     }
        //     std::cout << std::endl;
        // }

    }

}