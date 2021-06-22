#include "../inc/params1.h"
#include "../inc/array4d.h"
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/geom.h"
#include "../inc/NeighbourList.h"
#include "../inc/mathfuncs.h"
#include <vector>
#include <iostream>
#include <unordered_set>

namespace geom {

    Array<int> lw;
    Array<int> rw;
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

    void InitSpins(){

        Array4D<double> Sx4, Sy4, Sz4;
		Sx4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
		Sy4.resize(params::Lx, params::Ly, params::Lz, params::Nq);
		Sz4.resize(params::Lx, params::Ly, params::Lz, params::Nq);

		Sx4.IFill(0);
		Sy4.IFill(0);
		Sz4.IFill(0);


		int count1d = 0;

		if (params::afmflag != "NiO"){
		for (int x = 0; x < params::Lx; x++){
			for (int y = 0; y < params::Ly; y++){
				for (int z = 0; z < params::Lz; z++){
					for (int q = 0; q < params::Nq; q++){
						Sx4(x,y,z,q) = params::initm[q][0];
						Sy4(x,y,z,q) = params::initm[q][1];
						Sz4(x,y,z,q) = params::initm[q][2];

                        // testing domain walls
                        if (x < params::Lx/2){
                            Sx4(x,y,z,q) = params::initm[q][0];
                            Sy4(x,y,z,q) = params::initm[q][1];
                            Sz4(x,y,z,q) = params::initm[q][2];                       
                        }
                        else if (x > params::Lx/2){
                            Sx4(x,y,z,q) = -1 * params::initm[q][0];
                            Sy4(x,y,z,q) = -1 * params::initm[q][1];
                            Sz4(x,y,z,q) = -1 * params::initm[q][2];                         
                        }
                        else if (x == params::Lx/2){
                            Sx4(x,y,z,q) = 0.0;
                            Sz4(x,y,z,q) = 1.0;                        
                        }

						neigh::Sx1d(count1d + q) = Sx4(x,y,z,q);
						neigh::Sy1d(count1d + q) = Sy4(x,y,z,q);
						neigh::Sz1d(count1d + q) = Sz4(x,y,z,q);




					}	
					count1d += params::Nq; 
				}
			}
		}
		}
		else {
		for (int a = 0; a < params::Nspins; a++){
			if (a / params::Nq % 2 == 0){
				if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
					neigh::Sz1d(a) = 1.0; 
				}
				else if (modfunc(params::Nq,a) == 2) {
					neigh::Sz1d(a) = -1.0;
				}
				else {
					std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
					exit(0);
				}
			}
			else if (a / params::Nq % 2 == 1) {
				if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
					neigh::Sz1d(a) = -1.0; 
				}
				else if (modfunc(params::Nq,a) == 2) {
					neigh::Sz1d(a) = 1.0;
				}
				else {
					std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
					exit(0);
				}
			}
		}
		}

    }

    void InitDomainWall(){


		lw.resize(params::Ly*params::Lz);
		rw.resize(params::Ly*params::Lz);
		lw.IFill(0);
		rw.IFill(0);
    
        int inc = 0;

		// find all sites when x is 0
		for (int j = 0; j < params::Ly; j++){
			for (int k = 0; k < params::Lz; k++){
                lw(inc) = geom::LatCount(0,j,k,2);
                rw(inc) = geom::LatCount(params::Lx-1,j,k,2);
                inc++;
			}
		}

        std::cout << "Initialised Domain Wall \n";

	}

}