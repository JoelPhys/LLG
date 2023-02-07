// cpp header files
#include <vector>
#include <iomanip>
#include <iostream>
#include <unordered_set>

// my header files
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/config.h"
#include "../inc/array4d.h"
#include "../inc/array3d.h"
#include "../inc/defines.h"
#include "../inc/mathfuncs.h"

namespace geom {

    Array<int> lw;
    Array<int> rw;

   
    Array3D<double> Sx, Sy, Sz;
    Array4D<double> latticeX;
    Array4D<double> latticeY;
    Array4D<double> latticeZ;
    Array4D<int> LatCount;
    Array3D<int> Scount;
    Array<int> xlayer;
    Array<int> ylayer;
    Array<int> zlayer;
    Array<int> block;

    int Ix, Iy, Iz, IzC, nblocks;

    //testing for hedgehog
    Array<double> surfx;
    Array<double> surfy;
    Array<double> surfz;

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
        xlayer.resize(params::Nspins);
        ylayer.resize(params::Nspins);
        zlayer.resize(params::Nspins);

        latticeX.IFill(0);
        latticeY.IFill(0);
        latticeZ.IFill(0);
        LatCount.IFill(0);
        xlayer.IFill(0);
        ylayer.IFill(0);
        zlayer.IFill(0);


        // Discretise into micromagnetic cells
        int bx;
        int by;
        int bz;
        
		params::cfgmissing("Util.blockx");	
		params::cfgmissing("Util.blocky");
		params::cfgmissing("Util.blockz");		
		
		bx = params::cfg.lookup("Util.blockx");
        by = params::cfg.lookup("Util.blocky");
        bz = params::cfg.lookup("Util.blockz");
        INFO_OUT("blockx:", bx);
        INFO_OUT("blocky:", by);
        INFO_OUT("blockz:", bz);

        int blockx, blocky, blockz;
        block.resize(params::Nspins);
        int nblockx = params::Lx/bx;
        int nblocky = params::Ly/by;
        int nblockz = params::Lz/bz;
        nblocks = (nblockx)*(nblocky)*(nblockz);

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
                        xlayer(counter) = x;
                        ylayer(counter) = y;
                        zlayer(counter) = z;

                        // For micromagnetic cells
                        blockx = x / bx;
                        blocky = y / by;
                        blockz = z / bz;
                        block(counter) =  blockz + blocky*nblockz + blockx*nblocky*nblockz;
                        
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
    }

//TODO: Move this to spins file
void initdw(){

		// lw.resize(params::Ly*params::Lz*params::Nq);
        
        // testing for hedgehogs - points on surface of cubic lattice 
		//int nsurface = (2*(params::Lx*params::Ly) + 2*(params::Lx*(params::Ly-2)) + 2*((params::Lx-2)*(params::Ly-2)))*params::Nq;
        //lw.resize(nsurface);
        //surfx.resize(nsurface);
        //surfy.resize(nsurface);
        //surfz.resize(nsurface);
		
		lw.resize(params::Ly*params::Lz*params::Nq);
		rw.resize(params::Ly*params::Lz*params::Nq);
		lw.IFill(0);
		rw.IFill(0);
        
		int inc = 0;
        int test = 1;
		DEBUGGER;
		// initialise domain wall along x-direction.
		// set spins at either half of domain wall to anti-parallel alignment.
        for (int i = 0; i < params::Lx; i++){
            for (int j = 0; j < params::Ly; j++){
                for (int k = 0; k < params::Lz; k++){
                    for (int q = 0; q < params::Nq; q++){

						// up spins in first half of lattice
						if (i < params::Lx/2){
							spins::sx1d(geom::LatCount(i,j,k,q)) = params::initm[q][0];
                    	    spins::sy1d(geom::LatCount(i,j,k,q)) = params::initm[q][1];
                    	    spins::sz1d(geom::LatCount(i,j,k,q)) = params::initm[q][2];                       
                    	}
						// down spins in second half of lattice
                    	else if (i > params::Lx/2){
							spins::sx1d(geom::LatCount(i,j,k,q)) = -1 * params::initm[q][0];
                    	    spins::sy1d(geom::LatCount(i,j,k,q)) = -1 * params::initm[q][1];
                    	    spins::sz1d(geom::LatCount(i,j,k,q)) = -1 * params::initm[q][2];                         
                    	}
						// spins in middle set to another orientation so that domain wall forms in middle of sample
                    	else if (i == params::Lx/2){
							spins::sx1d(geom::LatCount(i,j,k,q)) = 0.5;
                    	    spins::sz1d(geom::LatCount(i,j,k,q)) = 0.5;                        
                    	}

						// for sublattice 1	
                        //lw(inc) = geom::LatCount(0,j,k,q);
                        //rw(inc) = geom::LatCount(params::Lx-1,j,k,q);
                        //
						//// for sublattice 2
						//rw(inc) = geom::LatCount(0,j,k,q);
                        //lw(inc) = geom::LatCount(params::Lx-1,j,k,q);

                    }
                }
            }
		}
        
		std::cout << "Initialised Domain Wall \n";

	}

}
