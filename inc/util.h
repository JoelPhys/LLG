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

	// outputs to files and terminal
    void OutputMagToFile(int i);
    void OutputFldToFile(int i);
    void OutputDWtoFile(int i);
    void OutputMagToTerm(int i);
    void OutputSpinToFile(int i);
    
	// initialise files
	void InitMagFile(double temp);
    void InitFldFile(double temp);
    void InitDWFile(double temp);
    void InitSpinFile(double temp);
    
	void CloseMagFile();
    void SumMag(int i);
    void OutputSumMag();
	
	// outputting lattice
	void OutputLatticetoFile(double temp);
    void OutputLatticeAverageOverX();
    void OutputLatticeAverageOverY();
    void OutputLatticeAverageOverZ();
    void OutputLatticeAverageOverQ();
    
	// clock
	void startclock();
    void endclock();

	// read initial spin config from a file
    void readexternalspins(std::string);
}

#endif
