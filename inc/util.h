#include <cmath>
#include <fstream>
#include "../inc/config.h"
#include "../inc/neighbourlist.h"
#include "../inc/array.h"
#include "../inc/array2d.h"


namespace util {

    extern Array2D<double> M;
    extern Array<double> Mt;
    extern Array<double> MmaG;
    extern Array<double> MdivMs;
    extern Array<double> MdivMsSum;
    extern Array2D<double> Msum;
    extern Array2D<double> MsumSQR;

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
    

}