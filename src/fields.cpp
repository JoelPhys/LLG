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
#include <cstring>
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
	std::string type;
	bool Hfield;
	Array<double> H_appx;
	Array<double> H_appy;
	Array<double> H_appz;

	double start_time = 0.0;
	double end_time = 0.0;
	double height;
	double freq;
	double gauss;
	double cuniform[3];
	int sublatsites;

	double std_dev;
	double centre_pos;

	void readfields(){

		H_appx.resize(params::Nmoments);
		H_appy.resize(params::Nmoments);
		H_appz.resize(params::Nmoments);

		H_appx.IFill(0);
		H_appy.IFill(0);
		H_appz.IFill(0);

		type = params::cfg.lookup("ExternalField.Type").c_str();
		libconfig::Setting& setting1 = params::cfg.lookup("ExternalField");

		TITLE("EXTERNAL FIELD");

		if (type == "Uniform") {
			INFO_OUT("Field type:", type);
			for (int a = 0; a < params::Nmoments; a++){
				H_appx(a) = setting1["Field"][0];
				H_appy(a) = setting1["Field"][1];
				H_appz(a) = setting1["Field"][2];
			}
			cuniform[0] = static_cast<double>(setting1["Field"][0]);
			cuniform[1] = static_cast<double>(setting1["Field"][1]);
			cuniform[2] = static_cast<double>(setting1["Field"][2]);
			INFO_OUT("Field values:", "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)");              
		}
		else if (type == "Uniform_Staggered") {
			INFO_OUT("Field type:", type);
			for (int a = 0; a < params::Nmoments; a++){
				
				sublatsites = params::sublat_sites[a % params::Nq];

				if (sublatsites == 0){
					H_appx(a) = setting1["Field"][0];
					H_appy(a) = setting1["Field"][1];
					H_appz(a) = setting1["Field"][2];
				}
				else if (sublatsites == 1){
					H_appx(a) = -1.0 * static_cast<double>(setting1["Field"][0]);
					H_appy(a) = -1.0 * static_cast<double>(setting1["Field"][1]);
					H_appz(a) = -1.0 * static_cast<double>(setting1["Field"][2]);
				}
				else std::cout << "WARNING: unasigned modulo value  = " << modfunc(params::Nq,a) << std::endl;
			}
			cuniform[0] = static_cast<double>(setting1["Field"][0]);
			cuniform[1] = static_cast<double>(setting1["Field"][1]);
			cuniform[2] = static_cast<double>(setting1["Field"][2]);
			INFO_OUT("Field values:", "[" << static_cast<double>(setting1["Field"][0]) << " , " << static_cast<double>(setting1["Field"][1]) << " , " << static_cast<double>(setting1["Field"][2]) << "] (T)");
		}
		else if ((type == "Square_Pulse") || (type == "Square_Pulse_Staggered")){
			INFO_OUT("Field type:", type);
			height = params::cfg.lookup("ExternalField.height");
			start_time = params::cfg.lookup("ExternalField.start_time");
			end_time = params::cfg.lookup("ExternalField.end_time");
			INFO_OUT("Start time of pulse = ", start_time << " timesteps");
			INFO_OUT("End time of pulse = ", end_time << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " (T)");
		}
		else if (type == "Gaussian_Pulse"){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre_pos");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			INFO_OUT("Central Position of Pulse = ", centre_pos << " timesteps");
			INFO_OUT("Standard Deviation of Pulse = ", std_dev << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " (T)");
		}
		else if ((type == "Multi_Cycle_Pulse") || (type == "Multi_Cycle_Pulse_Staggered")){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre_pos");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			freq = params::cfg.lookup("ExternalField.freq");
			INFO_OUT("Central Position of Pulse = ", centre_pos << " [s]");
			INFO_OUT("Standard Deviation of Pulse = ", std_dev << " [s]");
			INFO_OUT("Magniture of pulse = ", height << " [T]");
			INFO_OUT("Frequency of pulse = ", freq << " [Hz]");
		}
		else if ((type == "Sine_Pulse") || (type == "Sine_Pulse_Staggered")){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			freq = params::cfg.lookup("ExternalField.freq");
			INFO_OUT("Magniture of pulse = ", height << " [T]");
			INFO_OUT("Frequency of pulse = ", freq << " [Hz]");
		}
		else {	
			std::cout << "ERROR: Unknown Field type." << std::endl;
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
	
	void multi_cycle_pulse_staggered(double time){

		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev))) * sin(2.0*M_PI*freq*(time - centre_pos));

		for (int i = 0; i < params::Nspins; i++){
			
			int sublatsites = params::sublat_sites[i % params::Nq];
			
			if (sublatsites == 0){
				H_appx[i] = cos(0.25*M_PI)*gauss;
				H_appy[i] = cos(0.25*M_PI)*gauss;
				H_appz[i] = 0.0;  
			}
			else if (sublatsites == 1){
				H_appx[i] = -1.0 *cos(0.25*M_PI)* gauss;
				H_appy[i] = -1.0 *cos(0.25*M_PI)* gauss;
				H_appz[i] = 0.0;  
			}
		}


	}


	void calculate(double time){
		if (type == "Uniform"){
			(time);
		}
		else if (type == "split"){
			(time);
		}
		else if (type == "Square_Pulse"){
			square_pulse(time);
		}
		else if (type == "Gaussian_Pulse"){
			gaussian_pulse(time);
		}
		else if (type == "Multi_Cycle_Pulse"){
			multi_cycle_pulse(time);
		}
		else if (type == "Multi_Cycle_Pulse_Staggered"){
			multi_cycle_pulse_staggered(time);
		}
		else {
			std::cout << "ERROR. Field type not programmed in src/fields.cpp. Exiting." << std::endl;
			exit(0);
		}
	}

}
