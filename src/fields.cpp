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
#include <string>

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
	double heightac;
	double freq;
	double gauss;
	double gauss1;
	double gauss2;
	double kpoint;
	double kstep;
	double cuniform[3];
	int sublatsites;
	double std_dev;
	double centre_pos;
	int npump;

	// field direction and magnitude
	double direc_mag;
	double direc[3];

	// sublattice staggered order
	std::vector<int> sublat_stag;

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

		// Calculate direction of field	
		params::cfgmissing("ExternalField.direction");	
		direc[0] = static_cast<double>(setting1["direction"][0]);
		direc[1] = static_cast<double>(setting1["direction"][1]);
		direc[2] = static_cast<double>(setting1["direction"][2]);
		INFO_OUT("Field direction:", "[" << direc[0]<< " , " << direc[1]<< " , " << direc[2]<< "] [arb.]");              
		
		// Calculate magnitude of field along each direction	
		direc_mag = 1.0 / sqrt((direc[0]*direc[0]) + (direc[1]*direc[1]) + (direc[2]*direc[2]));
		INFO_OUT("Normalisation value of field along x,y,z:",direc_mag);	

		// If the field is staggered, assign the field to the desired sublattice
		int index_str, position;
		std::string sublatstr;
		sublat_stag.resize(params::Nsublat);

		// Check setting exists in cfg file
		params::cfgmissing("ExternalField.SublatStagger");
		
		// Check if sublattice staggering is the same length as the number of sublattices
		if (setting1["SublatStagger"].getLength() != params::Nsublat){ 
			std::cout << "ERROR: Sublattice Staggering is not the same length as number of sublattices. \n Exiting." << std::endl; 
			exit(0);
		}	

		for (int sublat = 0; sublat < params::Nsublat; sublat++){
			
			// read sublat staggering from config file
			sublat_stag[sublat] = static_cast<int>(setting1["SublatStagger"][sublat]);
			
			// print the field staggering to stdout
			sublatstr = "Sublattice " + std::to_string(sublat) + " being multiplied by:";
			INFO_OUT(sublatstr, sublat_stag[sublat]);
					
			// check if the staggering is 0,1,-1. If is isn't, exit the program.	
			if ((sublat_stag[sublat] != 0) && (sublat_stag[sublat] != 1) && (sublat_stag[sublat] != -1)){
				std::cout << "ERROR: Unexpected value for staggered field. \n";
				std::cout << "sublat_stag[" << sublat << "] = "  << sublat_stag[sublat] << "\n";
				std::cout << "Exiting. \n";
				exit(0);
			}
		}	

		// Loop over the possible types of field
		if (type == "Uniform") {
			INFO_OUT("Field type:", type);
			cuniform[0] = static_cast<double>(setting1["Field"][0]);
			cuniform[1] = static_cast<double>(setting1["Field"][1]);
			cuniform[2] = static_cast<double>(setting1["Field"][2]);
				for (int a = 0; a < params::Nmoments; a++){
				sublatsites = params::sublat_sites[a % params::Nq];
				H_appx(a) =sublat_stag[sublatsites]*cuniform[0];
				H_appy(a) =sublat_stag[sublatsites]*cuniform[1];
				H_appz(a) =sublat_stag[sublatsites]*cuniform[2];
			}
			
			INFO_OUT("Field values:", "[" << cuniform[0] << " , " << cuniform[1] << " , " << cuniform[2] << "] [T]");
		}
		else if ((type == "Square_Pulse")){
			INFO_OUT("Field type:", type);
			height = params::cfg.lookup("ExternalField.height");
			
			params::cfgmissing("ExternalField.start_time");
			params::cfgmissing("ExternalField.end_time");
			start_time = params::cfg.lookup("ExternalField.start_time");
			end_time = params::cfg.lookup("ExternalField.end_time");
			INFO_OUT("Start time of pulse = ", start_time << " timesteps");
			INFO_OUT("End time of pulse = ", end_time << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " [T]");

			//////////////////////////////////////////////////////////////////////////////////////////
			// TODO: create a way that two fields of different types can be simulated at the same time
			INFO_OUT("Field type = ", type);
			params::cfgmissing("ExternalField.heightac");
			heightac = params::cfg.lookup("ExternalField.heightac");
			freq = params::cfg.lookup("ExternalField.freq");
			kpoint = params::cfg.lookup("ExternalField.kpoint");
			INFO_OUT("Magniture of pulse:", height << " [T]");
			INFO_OUT("Frequency of pulse:", freq << " [Hz]");
			INFO_OUT("kpoint of pulse:", kpoint << " [pi]");
	
			// check if pumping upto a certain atom
			if (params::cfg.exists("ExternalField.layer")){
				npump = params::cfg.lookup("ExternalField.layer");
			}			
			else {
				npump = params::Nspins;
			}

			// Print npump
			INFO_OUT("pumping spins up to:", "N = " << npump);
			//////////////////////////////////////////////////////////////////////////////////////////


		}
		else if (type == "Gaussian_Pulse"){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			centre_pos = params::cfg.lookup("ExternalField.centre_pos");
			std_dev = params::cfg.lookup("ExternalField.std_dev");
			INFO_OUT("Central Position of Pulse = ", centre_pos << " timesteps");
			INFO_OUT("Standard Deviation of Pulse = ", std_dev << " timesteps");
			INFO_OUT("Magniture of pulse = ", height << " [T]");
		}
		else if ((type == "Multi_Cycle_Pulse")){
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
		else if ((type == "Sine_Pulse")){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			freq = params::cfg.lookup("ExternalField.freq");
			INFO_OUT("Magniture of pulse = ", height << " [T]");
			INFO_OUT("Frequency of pulse = ", freq << " [Hz]");
		}
		else if ((type == "Sine_Pulse_Linear") || (type == "Sine_Pulse_Circular")){
			INFO_OUT("Field type = ", type);
			height = params::cfg.lookup("ExternalField.height");
			freq = params::cfg.lookup("ExternalField.freq");
			kpoint = params::cfg.lookup("ExternalField.kpoint");
			INFO_OUT("Magniture of pulse:", height << " [T]");
			INFO_OUT("Frequency of pulse:", freq << " [Hz]");
			INFO_OUT("kpoint of pulse:", kpoint << " [pi]");
	
			// check if pumping upto a certain atom
			if (params::cfg.exists("ExternalField.layer")){
				npump = params::cfg.lookup("ExternalField.layer");
			}			
			else {
				npump = params::Nspins;
			}

			// Print npump
			INFO_OUT("pumping spins up to:", "N = " << npump);

		}
		else {	
			std::cout << "ERROR: Unknown Field type." << std::endl;
			exit(0);
		}

		

	}

	void square_pulse(double time){


		if ((time >= start_time) && (time < end_time)){
			for (int i = 0; i < params::Nspins; i++){
				sublatsites = params::sublat_sites[i % params::Nq];
				H_appx(i) = sublat_stag[sublatsites]*direc_mag*direc[0]*height;
				H_appy(i) = sublat_stag[sublatsites]*direc_mag*direc[1]*height;
				H_appz(i) = sublat_stag[sublatsites]*direc_mag*direc[2]*height;  
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

		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev)));

		for (int i = 0; i < params::Nspins; i++){
			sublatsites = params::sublat_sites[i % params::Nq];
			H_appx(i) = sublat_stag[sublatsites]*direc_mag*direc[0]*gauss;
			H_appy(i) = sublat_stag[sublatsites]*direc_mag*direc[1]*gauss;
			H_appz(i) = sublat_stag[sublatsites]*direc_mag*direc[2]*gauss;  

		}

	}
	
	void multi_cycle_pulse(double time){

		gauss = height * exp(-1.0 * (((time - centre_pos) * (time - centre_pos))/(2.0 * std_dev * std_dev))) * sin(2.0*M_PI*freq*(time - centre_pos));

		for (int i = 0; i < params::Nspins; i++){
			
			sublatsites = params::sublat_sites[i % params::Nq];
			
			H_appx[i] = sublat_stag[sublatsites]*direc_mag*direc[0]*gauss;
			H_appy[i] = sublat_stag[sublatsites]*direc_mag*direc[1]*gauss;
			H_appz[i] = sublat_stag[sublatsites]*direc_mag*direc[2]*gauss;  
			
		}


	}

	// TODO: Fix the polarisation of the fields
	void sine_pulse_circular(double time){


		for (int i = 0; i < params::Nspins; i++){
		
            sublatsites = params::sublat_sites[i % params::Nq];
			gauss  = heightac * sin(kpoint * M_PI * static_cast<double>(i) + 2.0*M_PI*freq*time);
            gauss2 = heightac * cos(kpoint * M_PI * static_cast<double>(i) + 2.0*M_PI*freq*time);
			//gauss1 = 0.0;
			//gauss2 = 0.0;	   
			//
			//for (int k = 0; k < params::Nspins; k++){
			//	kstep = static_cast<double>(k)/static_cast<double>(params::Nspins);
			//	gauss1 += sublat_stag[sublatsites] * sin( kstep * M_PI * i + 2.0*M_PI*freq*time);
			//	gauss2 += sublat_stag[sublatsites] * cos( kstep * M_PI * i + 2.0*M_PI*freq*time);
			//}

			H_appx[i] += sublat_stag[sublatsites] * gauss;
            H_appy[i] += sublat_stag[sublatsites] * gauss2;
            H_appz[i] += sublat_stag[sublatsites] * 0.0; 	
			
			//H_appx[i] = height * gauss;
            //H_appy[i] = height * gauss2;
            //H_appz[i] = 0.0; 	
		
		}

	}

	void sine_pulse_linear(double time){


		for (int i = 0; i < npump; i++){
		
            sublatsites = params::sublat_sites[i % params::Nq];
            gauss = height * sin(kpoint * M_PI * static_cast<double>(i) + 2.0*M_PI*freq*time);
            H_appx[i] = sublat_stag[sublatsites] * direc_mag*direc[0]*gauss;
            H_appy[i] = sublat_stag[sublatsites] * direc_mag*direc[1]*gauss;
            H_appz[i] = sublat_stag[sublatsites] * direc_mag*direc[2]*gauss;
		
		}

	}



	void calculate(double time){
		if (type == "Square_Pulse"){
			square_pulse(time);
			
			//////////////////////////////////////////////////////////////////////////////////////////
			// TODO: create a way that two fields of different types can be simulated at the same time
			if ((time < start_time) || (time >= end_time)){
				sine_pulse_circular(time);	
			}
			//////////////////////////////////////////////////////////////////////////////////////////
			
		}
		else if (type == "Gaussian_Pulse"){
			gaussian_pulse(time);
		}
		else if (type == "Multi_Cycle_Pulse"){
			multi_cycle_pulse(time);
		}
		else if (type == "Sine_Pulse_Linear"){
			sine_pulse_linear(time);	
		}
		else if (type == "Sine_Pulse_Circular"){
			sine_pulse_circular(time);	
		}
		else if (type == "Uniform"){
		}
		else {
			std::cout << "ERROR. Field type not programmed in src/fields.cpp. Exiting." << std::endl;
			exit(0);
		}
	}

}
