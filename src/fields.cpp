#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sstream>
#include "../inc/mathfuncs.h"
#include "../inc/array2d.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/error.h"

namespace fields {


	// Applied field settings
	std::string Type;
	bool Hfield;
	Array<double> H_appx;
	Array<double> H_appy;
	Array<double> H_appz;

	double start_time;
	double end_time;
	double height;
	double freq;

	double std_dev;
	double centre_pos;

	void readfields(){

		H_appx.resize(params::Nmoments);
		H_appy.resize(params::Nmoments);
		H_appz.resize(params::Nmoments);

		H_appx.IFill(0);
		H_appy.IFill(0);
		H_appz.IFill(0);

		Type = params::cfg.lookup("ExternalField.Type").c_str();
		libconfig::Setting& setting1 = params::cfg.lookup("ExternalField");

		if (Type == "Uniform") {
			std::cout.width(75); std::cout << std::left << "Field type:"; std::cout << Type << std::endl;
			for (int a = 0; a < params::Nmoments; a++){
				H_appx(a) = setting1["Field"][0];
				H_appy(a) = setting1["Field"][1];
				H_appz(a) = setting1["Field"][2];
			}
			std::cout.width(75); std::cout << std::left << "Field values:"; std::cout << "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)" << std::endl;              
		}
		else if (Type == "Split") {
			std::cout.width(75); std::cout << std::left <<"Field type:"; std::cout << Type << std::endl;
			for (int a = 0; a < params::Nmoments; a++){
				if ((modfunc(params::Nq,a) == 0) || (modfunc(params::Nq,a) == 2)) {
					H_appx(a) = setting1["Field"][0];
					H_appy(a) = setting1["Field"][1];
					H_appz(a) = setting1["Field"][2];
				}
				else if ((modfunc(params::Nq,a) == 1) || (modfunc(params::Nq,a) == 3)) {
					H_appx(a) = -1.0 * static_cast<double>(setting1["Field"][0]);
					H_appy(a) = -1.0 * static_cast<double>(setting1["Field"][1]);
					H_appz(a) = -1.0 * static_cast<double>(setting1["Field"][2]);
				}
				else std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
			}
			std::cout.width(75); std::cout << std::left << "Field values:"; std::cout << "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)" << std::endl;
		}
		else if (Type == "Square_Pulse"){
			std::cout.width(75); std::cout << std::left << "Field type:"; std::cout << Type << std::endl;
			height = params::cfg.lookup("ExternalField.height");
			start_time = params::cfg.lookup("ExternalField.start_time");
			end_time = params::cfg.lookup("ExternalField.end_time");
			std::cout.width(75); std::cout << std::left << "Start time of pulse = "; std::cout << start_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "End time of pulse = "; std::cout << end_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "Magniture of pulse = "; std::cout << height << " (T) \n";
		}
		else if (Type == "Gaussian_Pulse"){
			std::cout << "Field type = "; std::cout << Type << std::endl;
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			std::cout.width(75); std::cout << std::left << "Central Position of Pulse = "; std::cout << start_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "Standard Deviation of Pulse = "; std::cout << end_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "Magniture of pulse = "; std::cout << height << " (T) \n";
		}
		else if (Type == "Multi_Cycle_Pulse"){
			std::cout << "Field type = "; std::cout << Type << std::endl;
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			freq = params::cfg.lookup("ExternalField.frequency");
			std::cout.width(75); std::cout << std::left << "Central Position of Pulse = " ; std::cout<< start_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "Standard Deviation of Pulse = "; std::cout << end_time << " timesteps \n";
			std::cout.width(75); std::cout << std::left << "Magniture of pulse = "; std::cout << height << " (T) \n";
			std::cout.width(75); std::cout << std::left << "Frequency of pulse = "; std::cout << height << " (Hz) \n";
		}
		else {	
			std::cout << "ERROR: Unknown Field Type." << std::endl;
			exit(0);
		} 

	}
}
