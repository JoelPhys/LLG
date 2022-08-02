// cpp header files
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// My Header Files
#include "../inc/defines.h"
#include "../inc/config.h"
#include "../inc/libconfig.h++"
#include "../inc/util.h"


namespace thermal {

	// Temperature
    std::string temptype;
    double ttm_start;
    double temp_gradient;

    // Define global variables
    std::vector<double> ttime;     // time
    std::vector<double> tz;        // z plane
    std::vector<double> Te;        // electron temp
    std::vector<double> P_it;        // pump power
    std::vector<double> Tp;        // phonon temp
    
	//two temperature model variables
    double gamma_e;           
    double Cp;                
    double kappa_0;           
    double delta;             
    double Gep;               
    double P_0;               
    double t0;                
    double tau;               
	double oneOvrdzdz;
    double oneOvr2dz;
	double Tpp1 = 0.0;
	double Tep1 = 0.0;

	// Output temperature profile to file
	std::ofstream tempfile;
	
    void initthermal(double temp){
		
		// Resize arrays
		Te.resize(params::Lz,temp);   
        P_it.resize(params::Lz,0.0); 
        Tp.resize(params::Lz,temp);   
		
		// Temparature
		params::cfgmissing("Temperature.method");					
		params::cfgmissing("Temperature.ttm_start");				
		temptype = params::cfg.lookup("Temperature.method").c_str();
		ttm_start = params::cfg.lookup("Temperature.ttm_start");
		
		if ((temptype == "uniform_gradient") && (params::cfg.exists("Temperature.gradient") == 0)){
			std::cout << "ERROR: Method is uniform_gradient but no gradient value has been provided" << std::endl;
			std::cout << "EXITING SIMULATION" << std::endl;
			exit(0);
		}
		else if ((temptype == "uniform_gradient")){
			double temp_gradient_km = params::cfg.lookup("Temperature.gradient");
			INFO_OUT("Temperature Gradient:", temp_gradient_km << "[K/m]");
			temp_gradient = temp_gradient_km * params::a1;
			INFO_OUT("Temperature Gradient:", temp_gradient << "[K / unit cell]");
		}

		INFO_OUT("Temperature method: ", temptype);
		INFO_OUT("two temperature model start time: ", ttm_start);
	
		if (temptype == "ttm"){
			
			// Check two-temperature model variables are in cfg file
			params::cfgmissing("Temperature.gamma_e");   
        	params::cfgmissing("Temperature.Cp");        
        	params::cfgmissing("Temperature.kappa_0");  
        	params::cfgmissing("Temperature.delta");  
        	params::cfgmissing("Temperature.Gep");     
        	params::cfgmissing("Temperature.P_0");    
        	params::cfgmissing("Temperature.t0");   
        	params::cfgmissing("Temperature.tau"); 

			// Read two-temperature model variables	
			gamma_e = params::cfg.lookup("Temperature.gamma_e");   
        	Cp      = params::cfg.lookup("Temperature.Cp");        
        	kappa_0 = params::cfg.lookup("Temperature.kappa_0");  
        	delta   = params::cfg.lookup("Temperature.delta");  
        	Gep     = params::cfg.lookup("Temperature.Gep");     
        	P_0     = params::cfg.lookup("Temperature.P_0");    
        	t0      = params::cfg.lookup("Temperature.t0");   
        	tau		= params::cfg.lookup("Temperature.tau"); 
        	
			// Reduced variables - c1 is lattice constant in z-direction
			oneOvrdzdz = 1./(params::c1*params::c1);
        	oneOvr2dz  = 1./(2.0*params::c1);

			// Output ttm variables to log
			INFO_OUT("ttm variable gamma_e: ", gamma_e << " [J/m^3/K^2]");   
            INFO_OUT("ttm variable Cp: ", Cp << " [J/m^3/K]");        
            INFO_OUT("ttm variable kappa_0: ", kappa_0 << " [J/m/K/s]");  
            INFO_OUT("ttm variable delta: ", delta << " [m]");  
            INFO_OUT("ttm variable Gep: ", Gep << " [J/m^3/s/K]");     
            INFO_OUT("ttm variable P_0: ", P_0 << " [J/m^3/s]");    
            INFO_OUT("ttm variable t0: ", t0 << " [s]");   
            INFO_OUT("ttm variable tau: ", tau << " [s]"); 


			// Sort out temp profile file
			std::stringstream sstr;
			sstr << params::filepath << "temp_file_tsteps_" << params::Nt << "_T_" << std::setw(4) << std::setfill('0') << temp << ".out";
			tempfile.open(sstr.str());

		}
	}
	
	void ttm(double time){
	
        // loop over layes
		for (int i = 0; i < params::Lz; i++){
				
			//P_it[i]=P_0*exp(-((time-t0)/tau)*((time-t0)/tau));
			//Tep1=Te[0] + (params::dt/(gamma_e*Te[0]))*(Gep*(Tp[0]-Te[0]) + P_it[0] + kappa_0*( (Te[0]/Tp[0]) * 2.0*(Te[1]-Te[0])*oneOvrdzdz));
			//Tpp1=Tp[0]+(params::dt*Gep/Cp)*(Te[0]-Tp[0]);
			  P_it[i]=P_0*exp(-((time-t0)/tau)*((time-t0)/tau));
              Tep1=Te[0] + (params::dt/(gamma_e*Te[0]))*(Gep*(Tp[0]-Te[0]) + P_it[0] + kappa_0*( (Te[0]/Tp[0]) * 2.0*(Te[1]-Te[0])*oneOvrdzdz));
              Tpp1=Tp[0]+(params::dt*Gep/Cp)*(Te[0]-Tp[0]);

			// }
			// if (i == Nz-1)
			// {
			//     z=staticast<double>(Nz-1)*dz;
			//     P_it[Nz-1]=P_0*exp(-((time-t0)/tau)*((time-t0)/tau))*exp(-z/delta);
			//     Tep1=Te[Nz-1]+(dt/(gamma_e*Te[Nz-1]))*(Gep*(Tp[Nz-1]-Te[Nz-1])+P_it[Nz-1]+kappa_0*( (Te[Nz-1]/Tp[Nz-1]) * 2.0*(Te[Nz-2]-Te[Nz-1])*oneOvrdzdz));
			//     Tpp1=Tp[Nz-1]+(dt*Gep/Cp)*(Te[Nz-1]-Tp[Nz-1]);
			// }
			// if ((1 <= i) && (i < Nz-1))
			// {
			//     z=staticast<double>(i)*dz;
			//     P_it[i]=P_0*exp(-((time-t0)/tau)*((time-t0)/tau))*exp(-z/delta);
			//     Tep1=Te[i] + (dt/(gamma_e*Te[i]))*(Gep*(Tp[i]-Te[i]) + P_it[i]+kappa_0*( (Te[i]/Tp[i]) * (Te[i+1]-2.0*Te[i]+Te[i-1])*oneOvrdzdz+(Tp[i]*((Te[i+1]-Te[i-1])*oneOvr2dz) - Te[i]*(Tp[i+1]-Tp[i-1])*oneOvr2dz)/(Tp[i]*Tp[i])*((Te[i+1]-Te[i-1])*oneOvr2dz)));
			//     Tpp1=Tp[i]+(dt*Gep/Cp)*(Te[i]-Tp[i]);
			// }

			//update the values of Te[i] and Tp[i]
			Te[i]=Tep1;
			Tp[i]=Tpp1;

    	}

	}

	void cputemperature(double time){
		
		if (temptype == "ttm"){
			ttm(time);
		}
		else if (temptype == "constant"){

		}
		else if (temptype == "uniform_gradient") {

		}
		else {
			std::cout << "ERROR: Unknown temptype. Exiting." << std::endl;
			exit(0);
		}
	}


	void ttmtofile(double time){

		// Output to temp file
		tempfile << time << " " << P_it[0] << " " << Te[0] << " " << Tp[0] << "\n";	

	}

	void closettmfile(){
		// close temp file
		tempfile << std::flush;
		tempfile.close();
	}

    //void ReadThermalFile(){
	//	// Define variables
    //    std::cout << "READING THERMAL FILE " << std::endl;
	//	std::ifstream input(params::temp_filename);
	//	int a, b; 
    //    double c, d, e;
    //    int count = 0;

	//	// Read file
    //    if (!input){
	//		//Output Error if file can't be read. Quit program.
    //        std::cout << "ERROR: Could not open file " << params::temp_filename << std::endl;
    //        exit(0);
    //    }
	//	else {
	//		while (input >> a >> b >> c >> d >> e)
    //        {
	//			ttime.push_back(a);
	//			tz.push_back(b);
	//			pp.push_back(c);
	//			te.push_back(d);
	//			tp.push_back(e);
    //            count++;
    //        }
    //    	std::cout << "temp input file has been read" << std::endl;
	//	}

    //    // number of timesteps
    //    int nsteps = ttime.back();
    //    int nz = tz.back();
	//}

}
