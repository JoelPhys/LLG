#ifndef __SPINWAVES_H__
#define __SPINWAVES_H__

    #include <iostream>
    #include <sstream>
    #include <cmath>
    #include <fstream>
    #include <vector>
    #include <iomanip>
    #include <random>
    #include <algorithm>
    #include "fftw3.h"

    namespace spinwaves {

        extern int icount;
       

        //intiialise arrays
        void initialiseFFT();
        void FFTspace();
        void FFTtime();

    }

#endif
