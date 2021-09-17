// cpp header files
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

// my header files
#include "../inc/mathfuncs.h"
#include "../inc/array2d.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/defines.h"

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
	double gauss;

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
			INFO_OUT("Field type:", Type);
			for (int a = 0; a < params::Nmoments; a++){
				H_appx(a) = setting1["Field"][0];
				H_appy(a) = setting1["Field"][1];
				H_appz(a) = setting1["Field"][2];
			}
			INFO_OUT("Field values:", "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)");              
		}
		else if (Type == "Split") {
			INFO_OUT("Field type:", Type);
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
			INFO_OUT("Field values:", "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)");
		}
		else if (Type == "Square_Pulse"){
			INFO_OUT("Field type:", Type);
			height = params::cfg.lookup("ExternalField.height");
			start_time = params::cfg.lookup("ExternalField.start_time");
			end_time = params::cfg.lookup("ExternalField.end_time");
			INFO_OUT("Start time of pulse = ", start_time << " timesteps");
			INFO_OUT("End time of pulse = ", end_time << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " (T)");
		}
		else if (Type == "Gaussian_Pulse"){
			INFO_OUT("Field type = ", Type);
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			INFO_OUT("Central Position of Pulse = ", start_time << " timesteps");
			INFO_OUT("Standard Deviation of Pulse = ", end_time << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " (T)");
		}
		else if (Type == "Multi_Cycle_Pulse"){
			INFO_OUT("Field type = ", Type);
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			freq = params::cfg.lookup("ExternalField.frequency");
			INFO_OUT("Central Position of Pulse = ", start_time << " timesteps");
			INFO_OUT("Standard Deviation of Pulse = ", end_time << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " (T)");
			INFO_OUT("Frequency of pulse = ", height << " (Hz)");
		}
		else {	
			std::cout << "ERROR: Unknown Field Type." << std::endl;
			exit(0);
		} 

	}

	void square_pulse(double time){


		if ((time >= start_time) && (time < end_time)){
			for (int i = 0; i < params::Nspins; i++){
				H_appx(i) = height;
				H_appy(i) = 0.0;
				H_appz(i) = 0.0;  
			}
		}	
		else {
			for (int i = 0; i < params::Nspins; i++){
				H_appx(i) = 0.0;
				H_appy(i) = 0.0;
				H_appz(i) = 0.0; 
			}
		}

	}

	void gaussian_pulse(double time){

		gauss = height * exp(-1 * (((time - centre_pos) * (time - centre_pos))/(2 * std_dev * std_dev)));

		for (int i = 0; i < params::Nspins; i++){

			H_appx(i) = gauss;
			H_appy(i) = 0.0;
			H_appz(i) = 0.0; 

		}

	}

	void multi_cycle_pulse(double time){

		gauss = height * exp(-1 * (((time - centre_pos) * (time - centre_pos))/(2 * std_dev * std_dev))) * sin(2*M_PI*freq*(time - centre_pos));

		for (int i = 0; i < params::Nspins; i++){

			H_appx(i) = gauss;
			H_appy(i) = 0.0;
			H_appz(i) = 0.0; 

		}


	}

	void calculate(double time){
		if (Type == "Uniform"){
			(time);
		}
		else if (Type == "split"){
			(time);
		}
		else if (Type == "Square_Pulse"){
			square_pulse(time);
		}
		else if (Type == "Gaussian_Pulse"){
			gaussian_pulse(time);
		}
		else if (Type == "Multi_Cycle_Pulse"){
			multi_cycle_pulse(time);
		}
	}

}
