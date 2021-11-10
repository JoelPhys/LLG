
// cpp header files
#include <cmath>
#include <vector>
#include <iostream>


// my header files
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/array2d.h"
#include "../inc/defines.h"
#include "../inc/libconfig.h++"

namespace defects {

    int centrex;
    int centrey;
    int centrez;
    int radius;

    std::vector<int> list;

    void init(){

        if (params::cfg.exists("Defects") == 1){
            centrex = params::cfg.lookup("Defects.centrex");
            centrey = params::cfg.lookup("Defects.centrey");
            centrez = params::cfg.lookup("Defects.centrez");
            radius = params::cfg.lookup("Defects.radius");

            INFO_OUT("Spherical Defect:", "true");
            INFO_OUT("Defect Location:", "[" << centrex << ", " << centrey << ", " << centrez << "] [unit cells]");
            INFO_OUT("Spherical Defect radius:", radius << " [unit cells]");




            int minx = centrex - radius;
            int miny = centrey - radius;
            int minz = centrez - radius;
            int maxx = centrex + radius;
            int maxy = centrey + radius;
            int maxz = centrez + radius;
            int index;

            for (int x = minx; x <= maxx; x++){
                for (int y = miny; y <= maxy; y++){
                    for (int z = minz; z <= maxz; z++){
                        for (int q = 0; q < params::Nq; q++){
            
                            // loop through sites that are within simulation box
                            if ((minx >= 0) && (miny >= 0) && (minz >= 0) && (maxx < params::Lx) && (maxy < params::Ly) && (maxz < params::Lz)){
                                
                                // circle equation
                                if ((x - centrex)*(x - centrex) + (y - centrey)*(y - centrey) + (z - centrez)*(z - centrez) <= radius*radius){

                                    // find 1d index of spin and add it to the defect list.
                                    index = geom::LatCount(x,y,z,q);
                                    list.push_back(index);

                                }
                            }
                        }
                    }
                }
            }
     }
    }
    
    void populate(){
        int count = 0;
        int index = 0;
        if (params::cfg.exists("Defects") == 1){
            if (list.size() != 0){
                for (int x = 0; x < params::Lx; x++){
                    for (int y = 0; y < params::Ly; y++){
                        for (int z = 0; z < params::Lz; z++){
                            for (int q = 0; q < params::Nq; q++){

                                    index = geom::LatCount(x,y,z,q);
                                    if (index == list[count]){
                                        spins::sx1d(index) = 0.0;
                                        spins::sy1d(index) = 0.0;
                                        spins::sz1d(index) = 0.0;
                                        count++;
                                    }
                            }
                        }
                    }
                }
            }
        }
    }

}