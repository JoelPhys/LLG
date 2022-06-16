#ifndef __SPINWAVES_H__
#define __SPINWAVES_H__

    #include <iostream>
    #include <sstream>

    namespace spinwaves {

        // global variables
	    extern int start;
        extern int icount;
        extern double dt_spinwaves;        
        extern std::ofstream file_spnwvs;

        //intiialise arrays
        void init();
        void FFTspace();
        void FFTtime();

    }

#endif
