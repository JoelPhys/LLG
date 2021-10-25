
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

namespace defects {

    int centrex = 10;
    int centrey = 10;
    int centrez = 10;
    int radius = 5;

    std::vector<int> list;

    void init(){

    int minx = centrex - radius;
    int miny = centrex - radius;
    int minz = centrex - radius;
    int maxx = centrex + radius;
    int maxy = centrex + radius;
    int maxz = centrex + radius;
    int index;

        for (int x = minx; x <= maxx; x++){
            for (int y = miny; y <= maxy; y++){
                for (int z = minz; z <= maxz; z++){
                    for (int q = 0; q < params::Nq; q++){
        
                        if ((x - centrex)*(x - centrex) + (y - centrey)*(y - centrey) + (z - centrez)*(z - centrez) <= radius*radius){

                            index = geom::LatCount(x,y,z,q);
                            list.push_back(index);

                        }
                    }
                }
            }
        }

    }
    
    void populate(){
        int count = 0;
        int index = 0;
        
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