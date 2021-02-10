#include <cmath>
#include <fstream>
#include "../inc/params1.h"
#include "../inc/NeighbourList.h"
#include "../inc/array.h"
#include "../inc/array2d.h"


namespace util {

    extern Array2D<double> M;
    extern Array<double> MmaG;
    extern Array<double> MdivMs;
    extern Array<double> MdivMsSum;
    extern Array2D<double> Msum;
    extern Array2D<double> MsumSQR;

    void SortSublat();
    void InitUtil();
    void ResetMag();
    void MagLength();
    void OutputMagToFile(int i);
    void OutputMagToTerm(int i);
    void InitOutputFile(double temp);
    void CloseMagFile();
    void SumMag(int i);
    void OutputSumMag();
    

}