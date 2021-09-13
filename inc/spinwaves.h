#ifndef __SPINWAVES_H__
#define __SPINWAVES_H__

    #include <iostream>
    #include <sstream>

    namespace spinwaves {

        extern int icount;
        extern std::ofstream file_spnwvs;
        
       

        //intiialise arrays
        void init();
        void FFTspace();
        void FFTtime();

    }

#endif
