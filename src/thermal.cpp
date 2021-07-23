#include <fstream>
#include <vector>
#include <iostream>
#include "../inc/params1.h"
#include "../inc/libconfig.h++"
#include "../inc/util.h"


namespace thermal {


    // Define global variables
    std::vector<double> ttime;     // time
    std::vector<double> tz;        // z plane
    std::vector<double> te;        // electron temp
    std::vector<double> pp;        // pump power
    std::vector<double> tp;        // phonon temp

    int nsteps = 0;
    int nz = 0;


    void ReadThermalFile(){
		// Define variables
        std::cout << "READING THERMAL FILE " << std::endl;
		std::ifstream input(params::temp_filename);
		int a, b; 
        double c, d, e;
        int count = 0;

		// Read file
        if (!input){
			//Output Error if file can't be read. Quit program.
            std::cout << "ERROR: Could not open file " << params::temp_filename << std::endl;
            exit(0);
        }
		else {
			while (input >> a >> b >> c >> d >> e)
            {
				ttime.push_back(a);
				tz.push_back(b);
				pp.push_back(c);
				te.push_back(d);
				tp.push_back(e);
                count++;
            }
        	std::cout << "temp input file has been read" << std::endl;
		}

        // number of timesteps
        int nsteps = ttime.back();
        int nz = tz.back();
	}

}