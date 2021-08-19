#ifndef _SPINS_H_
#define _SPINS_H_

#include "../inc/array.h"

namespace spins {

    // Global Variables
    extern Array<double> sx1d;
    extern Array<double> sy1d;
    extern Array<double> sz1d;

    // Testing for hedgehog
    extern Array<double> surfx;
    extern Array<double> surfy;
    extern Array<double> surfz;

    extern Array<int> lw;
    extern Array<int> rw;
    extern Array<int> zlayer;

    // Functions
    void init();
    void populate();

}


#endif