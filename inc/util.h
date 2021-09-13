#ifndef _UTIL_H_
#define _UTIL_H_

// cpp header files
#include <cmath>
#include <fstream>

// my header files
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/array2d.h"
#include "../inc/neighbourlist.h"

namespace util {

    void SortSublat();
    void init();
    void ResetMag();
    void MagLength();
    void OutputMagToFile(int i);
    void OutputDWtoFile(int i);
    void OutputMagToTerm(int i);
    void InitMagFile(double temp);
    void InitDWFile(double temp);
    void CloseMagFile();
    void SumMag(int i);
    void OutputSumMag();
    void startclock();
    void endclock();
    void readexternalspins(std::string);
}

#endif