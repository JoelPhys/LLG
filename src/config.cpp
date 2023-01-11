// cpp header files
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>
#include <libconfig.h++>

// my header files
#include "../inc/array.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/array2d.h"
#include "../inc/defines.h"
#include "../inc/mathfuncs.h"

namespace params {

	libconfig::Config cfg;

	std::string simtype;

	double k_B, mu_b, gamma;
	double dt, Nt, dtau, half_dtau;
	int relaxtime, outputstep;

	std::vector<double> lambda;
	std::vector<double> lambdaPrime;
	std::vector<double> thermal_const;
	std::vector<double> mu_s;
	std::vector<double> INVmu_s;

	// Uniaxial Anisotropy
	std::vector<double> dxu, dyu, dzu;
	std::vector<double> dzup, dxup, dyup;

	//Cubic Anisotropy
	double dzc, dzcp;

	int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat;
	int Idx, Idy, Idz; // For integer lattice
	double a1, b1, c1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;
	double angle;
	std::vector<int> sublat_sites;
	std::vector<int> NmomentsSubLat;


	// Array for finding first instance of sublat within unit clel
	std::vector<int> uniquesublat;

	//Boundary Conditions
	std::string xbound;
	std::string ybound;
	std::string zbound;

	// Util variables
	std::string afmflag;
	std::string format;
	std::string filepath;
	std::string filepath_sw;
	bool OutputToTerminal;

	// Random Number Seed
	int seed;

	// Jij SETTINGS
	std::string Jij_filename;
	std::string Jij_units;
	std::string changesign;
	bool JijCutoff, Jijhalf;
	double Jij_min;
	int ibtoq;

	// Lattive Vectors
	double Plat[3][3];
	double PlatINV[3][3];

	// Site positions 
	std::vector< std::vector<double> > sites;
	std::vector< std::vector<int> > Isites;
	std::vector< std::vector<double> > initm;

	// Output Lattice
	bool OutputLattice = false;		
	int OutputLatticeStep = 10000000;
	int OutputLatticeStart = 10000000;	
	std::string OutputLatticeFilepath;		
	std::string OutputLatticeAverageOver = "false";

	void banner(){
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "                                       █████           █████████    █████████  ██████████                     " << std::endl;
		std::cout << "                                      ░░███           ███░░░░░███  ███░░░░░███░░███░░░░███                    " << std::endl;
		std::cout << "                                       ░███          ░███    ░███ ░███    ░░░  ░███   ░░███                   " << std::endl;
		std::cout << "                                       ░███   █████  ░███████████ ░░█████████  ░███    ░███                   " << std::endl;
		std::cout << "                                       ░███  ░░░░░   ░███░░░░░███  ░░░░░░░░███ ░███    ░███                   " << std::endl;
		std::cout << "                                 ███   ░███          ░███    ░███  ███    ░███ ░███    ███                    " << std::endl;
		std::cout << "                                ░░████████           █████   █████░░█████████  ██████████                     " << std::endl;
		std::cout << "                                 ░░░░░░░░           ░░░░░   ░░░░░  ░░░░░░░░░  ░░░░░░░░░░                      " << std::endl;
		std::cout << "                                                                                                              " << std::endl;
		std::cout << "                                             Joel's Atomistic Spin Dynamics                                   " << std::endl;
		std::cout << "                                                   														   " << std::endl;	
		std::cout << "                                                   														   " << std::endl;		
		std::cout << "                                                   														   " << std::endl;	
	}

	//=======================================================================================================
	// Check config file exists =============================================================================
	//=======================================================================================================

	void cfgmissing(std::string in){
		if (!cfg.exists(in)){
			std::cout << "ERROR: Missing config file setting " << in << std::endl;
			std::cout << "EXITING SIMULATION" << std::endl;
			exit(0);
		}
	}

	//=======================================================================================================
	//Intialise Config File =================================================================================
	//=======================================================================================================

	void intitialiseConfig(const char* cfg_filename){

		TITLE("COMPILATION INFO");
		INFO_OUT("CPU Compiler: ", CPUCOMP);
#ifdef CUDA 
		INFO_OUT("NVCC Compiler: ", GPUCOMP); 
#endif
		INFO_OUT("Compile Date and Time: ", __DATE__ << " " << __TIME__);
		INFO_OUT("Compiled on Machine: ", HOSTNAME);
		if(GITDIRTY!="0")
		{
			INFO_OUT("WARNING: Your git build is dirty. Recompile with Git SHA:", GIT_SHA1 << ", Dirty");
		}
		else
		{
			INFO_OUT("Git SHA: ", GIT_SHA1 << ", Clean");
		}

		if(!cfg_filename)
		{
			error::Info(__FILE__,__LINE__);
			error::Message("You must give a config file, exiting");
		}
		try
		{
			cfg.readFile(cfg_filename);
		}
		catch(const libconfig::FileIOException &fioex)
		{
			error::Info(__FILE__,__LINE__);
			error::Message("I/O error while reading config file");
		}
		catch(const libconfig::ParseException &pex)
		{
			error::Info(__FILE__,__LINE__);
			std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
			exit(EXIT_FAILURE);
		}
		cfg.setAutoConvert(true);
	}

	//=======================================================================================================
	// Read Parameters ======================================================================================
	//=======================================================================================================

	void readparams(){

		// Output Simulation Start time
		time_t now = time(0);
		std::string simtim = strtok(ctime(&now), "\n");
		INFO_OUT("Time of Simulation:", simtim);

		// Simulation Type
		cfgmissing("SimulationType");
		simtype = cfg.lookup("SimulationType").c_str();

		// Global Constants
		cfgmissing("PhysicalConsts.BoltzmannConstant"); 
		cfgmissing("PhysicalConsts.BohrMagneton");      
		cfgmissing("PhysicalConsts.GyromagneticRatio"); 
		k_B = cfg.lookup("PhysicalConsts.BoltzmannConstant");
		mu_b = cfg.lookup("PhysicalConsts.BohrMagneton");
		gamma = cfg.lookup("PhysicalConsts.GyromagneticRatio");

		// Time parameters
		cfgmissing("Time.SizeOfStep");					
		cfgmissing("Time.NumberOfSteps");			    
		cfgmissing("Time.RelaxationTime");			    
		cfgmissing("Time.OutputStep");					
		dt = cfg.lookup("Time.SizeOfStep");
		Nt = cfg.lookup("Time.NumberOfSteps");
		relaxtime = cfg.lookup("Time.RelaxationTime");
		outputstep = cfg.lookup("Time.OutputStep");

		// Reduced Time variables
		dtau = gamma * dt;
		half_dtau = 0.5 * dtau;   


		// Cubic Anisotropy
		cfgmissing("Cubic_Anisotropy.d_c");				
		dzc = cfg.lookup("Cubic_Anisotropy.d_c");			
		//dzcp = 2 * ( dzc / mu_s[0] );

		// system dimensions
		cfgmissing("Geom.UnitCellsInX");				
		cfgmissing("Geom.UnitCellsInY");				
		cfgmissing("Geom.UnitCellsInZ");				
		cfgmissing("Geom.NumberOfSites");				
		cfgmissing("Geom.NumberOfSublat");
		Lx = cfg.lookup("Geom.UnitCellsInX");
		Ly = cfg.lookup("Geom.UnitCellsInY");
		Lz = cfg.lookup("Geom.UnitCellsInZ");
		Nq = cfg.lookup("Geom.NumberOfSites");
		Nsublat = cfg.lookup("Geom.NumberOfSublat");				
		ax = 2;
		ay = 2;
		az = 2;

		// Angle of sublattice rotation
		angle = cfg.lookup("angle");
		angle *= M_PI / 180.0;

		// System geometry
		Nspins = Nq*Lx*Ly*Lz;
		Nmoments = (Nq*Lx*Ly*Lz); 
		NsitesINV_S = 1/(Lx*Ly*Lz);
		xdim = ax*Lx;
		ydim = ay*Ly;
		zdim = az*Lz;
		NsitesINV = 1/(xdim*ydim*zdim);
		zdimC = zdim/2+1;

		//Boundary Conditions
		cfgmissing("Geom.BoundaryConditionsX");				
		cfgmissing("Geom.BoundaryConditionsY");				
		cfgmissing("Geom.BoundaryConditionsZ");				
		xbound = cfg.lookup("Geom.BoundaryConditionsX").c_str(); 
		ybound = cfg.lookup("Geom.BoundaryConditionsY").c_str(); 
		zbound = cfg.lookup("Geom.BoundaryConditionsZ").c_str(); 
		
		// Lattice parameters
		cfgmissing("LatticeVectors.a");
		cfgmissing("LatticeVectors.b");
		cfgmissing("LatticeVectors.c");	
		a1 = cfg.lookup("LatticeVectors.a");
		b1 = cfg.lookup("LatticeVectors.b");
		c1 = cfg.lookup("LatticeVectors.c");

		//Read integer lattice spacing
		Idx = cfg.lookup("IntegerSites.Idx");
		Idy = cfg.lookup("IntegerSites.Idy");
		Idz = cfg.lookup("IntegerSites.Idz");

		// Print key parameters to log file
		TITLE("MATERIAL CONSTANTS");
		INFO_OUT("Lattice Parameter, a:", a1*1e9 << " [nm]");
		INFO_OUT("Lattice Parameter, b:", b1*1e9 << " [nm]");
		INFO_OUT("Lattice Parameter, c:", c1*1e9 << " [nm]");
		INFO_OUT("Timestep:", dt << " (s)");
		INFO_OUT("Number of timesteps:", Nt);
		INFO_OUT("Outputting every ", outputstep << " timesteps");
		INFO_OUT("Outputting to terminal: ", OutputToTerminal);
		INFO_OUT("Number of unit cells:", Lx << "x" << Ly << "x" << Lz);
		INFO_OUT("Number of sites in unit cell:", Nq);
		INFO_OUT("Number of sublattices:", Nsublat);
		INFO_OUT("Boundary Conditions:", "[" << xbound << " , " << ybound << " , " << zbound << "]");

		// resize std::vectors 
		lambdaPrime.resize(Nq);
		thermal_const.resize(Nq);
		INVmu_s.resize(Nq);
		NmomentsSubLat.resize(Nsublat);
		sites.resize(Nq);
		Isites.resize(Nq);
		initm.resize(Nq);
		lambda.resize(Nq);
		mu_s.resize(Nq);
		sublat_sites.resize(Nq);
		dxu.resize(Nq);		
		dyu.resize(Nq);
		dzu.resize(Nq);
		dxup.resize(Nq);		
		dyup.resize(Nq);
		dzup.resize(Nq);
	


		//=======================================================================================================
		// Read Lattice Vectors and constants ===================================================================
		//=======================================================================================================

		for (int v = 0; v < 3; v++){

			std::stringstream sstr1;
			sstr1 << "LatVec" << v;
			std::string str1 = sstr1.str();

			libconfig::Setting& setting = cfg.lookup("LatticeVectors");  
			cfgmissing("LatticeVectors." + str1);
			Plat[v][0] = setting[str1.c_str()][0];
			Plat[v][1] = setting[str1.c_str()][1];
			Plat[v][2] = setting[str1.c_str()][2];
			INFO_OUT("Lattice Vectors:", std::fixed << std::setprecision(5) << Plat[v][0] << " " << Plat[v][1] << " " << Plat[v][2]);
		}

		// GET MATERIAL STUFF EITHER FROM A SEPERATE UNIT CELL FILE OR FROM CONFIG
		if (cfg.exists("MaterialConsts.file")){

			// get filename
			std::string matfile_name = cfg.lookup("MaterialConsts.file");
			std::ifstream matfile(matfile_name);

			INFO_OUT("Unit cell File: ", matfile_name);

			//Check file can be opened
			if (!matfile.is_open()){
				std::cout << "ERROR: Could not open mat file. Exiting." << std::endl;
				exit(0);
			}
			
			// variables for reading matfile
			double f1,f2,f3,f4,f5,f9,f10,f11,f13,f14,f15;
			int f6,f7,f8,f12;
			int q =0;

			// Create table in log file
			std::cout << std::setw(8) << "lambda";
            std::cout << std::setw(8) << "mu_s";
            std::cout << std::setw(8) << "posx";
            std::cout << std::setw(8) << "posy";
            std::cout << std::setw(8) << "posz";
            std::cout << std::setw(8) << "iposx";
            std::cout << std::setw(8) << "iposy";
            std::cout << std::setw(8) << "iposz";
            std::cout << std::setw(8) << "mx";
            std::cout << std::setw(8) << "my";
            std::cout << std::setw(8) << "mz";
            std::cout << std::setw(8) << "sublat";
            std::cout << std::setw(12) << "dx";
            std::cout << std::setw(12) << "dy";
            std::cout << std::setw(12) << "dy";
            std::cout << std::endl;

			while (matfile >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10 >> f11 >> f12 >> f13 >> f14 >> f15)
            {

				// Print values to log file		
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)  << f1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f2;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f3;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f4;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f5;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)  << f6;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)  << f7;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)  << f8;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f9;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f10;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8) << f11;
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)  << f12;
                std::cout << std::scientific << std::setprecision(4) << std::setw(12) << f13;
                std::cout << std::scientific << std::setprecision(4) << std::setw(12) << f14;
                std::cout << std::scientific << std::setprecision(4) << std::setw(12) << f15;
                std::cout << std::endl;
	
				//resize 2D std::vectors
				sites[q].resize(3);
				Isites[q].resize(3);
				initm[q].resize(3);
			
				// Read Values from file
				lambda[q] = f1;
				mu_s[q] = f2;
				sites[q][0] = f3;
				sites[q][1] = f4;
				sites[q][2] = f5;
				Isites[q][0] = f6;
				Isites[q][1] = f7;
				Isites[q][2] = f8;
				initm[q][0] = f9;
				initm[q][1] = f10;
				initm[q][2] = f11;
				sublat_sites[q] = f12;
				dxu[q] = f13; 
                dyu[q] = f14; 
                dzu[q] = f15; 

				// count number of sites within sublattices
				NmomentsSubLat[f12] += Lx*Ly*Lz;

				// Sort out magnetic moments
				mu_s[q] *= mu_b;
				INVmu_s[q] = 1.0 / mu_s[q];

				// Thermal parameters
				thermal_const[q] = sqrt( (2.0 * lambda[q] * k_B)  / (mu_s[q] * dtau) );
				lambdaPrime[q] = 1.0 / (1.0+(lambda[q]*lambda[q]));
	
				// Uniaxial Anisotropy
				dxup[q] = 2*dxu[q] / mu_s[q];
				dyup[q] = 2*dyu[q] / mu_s[q];
				dzup[q] = 2*dzu[q] / mu_s[q];

				// Increment through sites in unit cell
				q++;

            }

			// Check mat file length is same as number of sites specified in cfg file
			if (q != Nq){
				std::cout << "ERROR: Number of lines in file does not match Nq in cfg file. Exiting" << std::endl;
				std::cout << q << " " << Nq << std::endl;
				exit(0);
			}
		}
		else {

			//=======================================================================================================
			//Read Site positions ===================================================================================
			//=======================================================================================================

			for (int s = 0; s < Nq; s++){
				sites[s].resize(3);

				std::stringstream sstr;
				sstr << "Site" << s;
				std::string str = sstr.str();

				libconfig::Setting& setting = cfg.lookup("Sites");  
				sites[s][0] = setting[str.c_str()][0];
				sites[s][1] = setting[str.c_str()][1];
				sites[s][2] = setting[str.c_str()][2];
			}

			//=======================================================================================================
			// Read integer site position ===========================================================================
			//=======================================================================================================

			for (int s = 0; s < Nq; s++){
				Isites[s].resize(3);

				std::stringstream sstr;
				sstr << "ISite" << s;
				std::string str = sstr.str();

				libconfig::Setting& setting = cfg.lookup("IntegerSites"); 
				cfgmissing("IntegerSites." + str);
				Isites[s][0] = setting[str.c_str()][0];
				Isites[s][1] = setting[str.c_str()][1];
				Isites[s][2] = setting[str.c_str()][2];
			}
			
			//=======================================================================================================
			//Read Magnetic Moments =================================================================================
			//=======================================================================================================

			std::string mutext;
			for (int i = 0; i < Nq; i++){

				// Read from config file
				libconfig::Setting& setting = cfg.lookup("MaterialConsts");  
				mu_s[i]   = setting["mu_s"][i];

				// Print information to log file
				mutext = "Magnetic Moment for site " + std::to_string(i) + ":";
				INFO_OUT(mutext, mu_s[i] << " [mu_b]")

				// Convert to SI units
				mu_s[i] *= mu_b;
				INVmu_s[i] = 1.0 / mu_s[i];
			}

			//=======================================================================================================
			//Read Damping Consts ===================================================================================
			//=======================================================================================================

			cfgmissing("MaterialConsts.lambda");	

			std::string lmtext;
			for (int i = 0; i < Nq; i++){

				// Read from config file
				libconfig::Setting& setting = cfg.lookup("MaterialConsts");  
				lambda[i] = setting["lambda"][i];

				// Print information to log file
				lmtext = "Damping Constant for site  " + std::to_string(i) + ":";
				INFO_OUT(lmtext, lambda[i] << " [arb.]");

				// Thermal parameters
				thermal_const[i] = sqrt( (2.0 * lambda[i] * k_B)  / (mu_s[i] * dtau) );
				lambdaPrime[i] = 1.0 / (1.0+(lambda[i]*lambda[i]));

			}


			//=======================================================================================================
			// Read Initial Magnetisation Vectors ===================================================================
			//=======================================================================================================

			for (int v = 0; v < Nq; v++){

				initm[v].resize(3);
				std::stringstream sstr2;
				sstr2 << "initm" << v;
				std::string str2 = sstr2.str();

				libconfig::Setting& setting = cfg.lookup("InitialMagnetisation");   
				initm[v][0] = setting[str2.c_str()][0];
				initm[v][1] = setting[str2.c_str()][1];
				initm[v][2] = setting[str2.c_str()][2];
				INFO_OUT("Initial Magnestaion Vectors:", std::fixed << std::setprecision(5) << initm[v][0] << " " << initm[v][1] << " " << initm[v][2]);
			}
			std::cout.unsetf(std::ios_base::fixed);

			//=======================================================================================================
			// Read Sublat vector ===================================================================================
			//=======================================================================================================
			cfgmissing("Geom.sublat_sites");
			for (int v = 0; v < Nq; v++){
				libconfig::Setting& setting = cfg.lookup("Geom");   
				sublat_sites[v] = setting["sublat_sites"][v];

				//count number of sites within each sublattice
				NmomentsSubLat[sublat_sites[v]] += 1*Lx*Ly*Lz;

			}
			
			// Output number of magnetic moments in each sublattice.
			for (int v = 0; v < Nsublat; v++){
				lmtext = "Number of Moments in Sublattice " + std::to_string(v) + ":";   		
				INFO_OUT(lmtext, NmomentsSubLat[v]);
			}

			//=======================================================================================================
			// Uniaxial Anisotropy ==================================================================================
			//=======================================================================================================
			cfgmissing("Uniaxial_Anisotropy.d_x");			
			cfgmissing("Uniaxial_Anisotropy.d_y");			
			cfgmissing("Uniaxial_Anisotropy.d_z");	
			cfgmissing("Cubic_Anisotropy.d_c");
			libconfig::Setting& setting = cfg.lookup("Uniaxial_Anisotropy");   
			
			if (setting["d_x"].getLength() != Nq){ std::cout << "ERROR: ani not same length as number of sites. Exiting." << std::endl; exit(0);}
			if (setting["d_y"].getLength() != Nq){ std::cout << "ERROR: ani not same length as number of sites. Exiting." << std::endl; exit(0);}
			if (setting["d_z"].getLength() != Nq){ std::cout << "ERROR: ani not same length as number of sites. Exiting." << std::endl; exit(0);}
			
			for (int v = 0; v < Nq; v++){	
					
				dxu[v] = setting["d_x"][v];
				dyu[v] = setting["d_y"][v];
				dzu[v] = setting["d_z"][v];			
				dxup[v] = 2 * ( dxu[v] / mu_s[v] );
				dyup[v] = 2 * ( dyu[v] / mu_s[v] );	
				dzup[v] = 2 * ( dzu[v] / mu_s[v] );

				//Cubic Anisotropy
				dzc = cfg.lookup("Cubic_Anisotropy.d_c");
				dzcp = 2 * ( dzc / mu_s[0] );

				// Print Anisotropy to terminal
				lmtext = "Uniaxial Anisotropy Energy for site " + std::to_string(v) + ":";
				INFO_OUT(lmtext, "[" << dxu[v] << " , " << dyu[v] << " , " << dzu[v] << "] [J]");
				lmtext = "Uniaxial Anisotropy Field for site " + std::to_string(v) + ":";
				INFO_OUT(lmtext, "[" << dxup[v] << " , " << dyup[v] << " , " << dzup[v] << "] [T]");
			}

			INFO_OUT("Cubic Anisotropy:", dzc << " [J]");

		}
		
		//=======================================================================================================
		// Find the index of first instance of each sublat in unit cell =========================================
		//=======================================================================================================

		for (int j =0 ; j < Nsublat; j++){
                for (int i = 0; i < sublat_sites.size(); i++){ 
					if (j == sublat_sites[i]){
           			     uniquesublat.push_back(i);
           			     std::cout << j << " " << i << std::endl;
           			     break;
           			 }
 			   	}	
		}
		
		//=======================================================================================================
		// Check filepaths exist ================================================================================
		//=======================================================================================================

		cfgmissing("Util.filepath");
		filepath = cfg.lookup("Util.filepath").c_str();      

		struct stat buffer;
		if (stat(filepath.c_str(), &buffer) != 0) {
				std::cout << "ERROR: Magnetisation output directory does not exist!" << std::endl;
			exit(0);
		}

		//=======================================================================================================
		// Output of Lattice ====================================================================================
		//=======================================================================================================

		cfgmissing("Util.OutputLattice");		
		OutputLattice = cfg.lookup("Util.OutputLattice");
		if (OutputLattice == true){
			cfgmissing("Util.OutputLatticeStep");	
			cfgmissing("Util.OutputLatticeFilepath");			
			cfgmissing("Util.OutputLatticeStart");			
			OutputLatticeStep = cfg.lookup("Util.OutputLatticeStep");  
			OutputLatticeStart = cfg.lookup("Util.OutputLatticeStart");  
			OutputLatticeFilepath = cfg.lookup("Util.OutputLatticeFilepath").c_str();  

			if (cfg.exists("Util.OutputLatticeAverageOver")){
				std::vector<std::string> averagelist{"false","x","y","z","q"};
				OutputLatticeAverageOver = cfg.lookup("Util.OutputLatticeAverageOver").c_str();
				if (std::find(std::begin(averagelist), std::end(averagelist), OutputLatticeAverageOver) != std::end(averagelist)){
					INFO_OUT("Averaging Lattice output over spin lattice component:",OutputLatticeAverageOver);
				}	
				else {
					std::cout << "ERROR: Unknown Util.OutputLatticeAverageOver specified in config file. \n";
					std::cout << "Util.OutputLatticeAverageOver = " << OutputLatticeAverageOver << std::endl;
					std::cout << "Exiting." << std::endl;
					exit(0);
				}
			}

	
			struct stat buffer;
			if (stat(OutputLatticeFilepath.c_str(), &buffer) != 0) {
					std::cout << "ERROR: Lattice output \"" <<  OutputLatticeFilepath << "\" directory does not exist!" << std::endl;
				exit(0);
			}
		}
	
		//=======================================================================================================
		// Random Number Seed ===================================================================================
		//=======================================================================================================
		if (cfg.exists("Util.Seed")){		
				seed = cfg.lookup("Util.Seed");		
		}	
		else {
			std::cout << "WARNING: Random Number Seed not specified in config file." << std::endl;
			std::cout << "Using std::random_device to generate seed." << std::endl;
			std::random_device device;
			seed = device();
		}	
		
		INFO_OUT("Random Number Seed:", seed);

		//=======================================================================================================
		// Rest of parameters ===================================================================================
		//=======================================================================================================

		cfgmissing("Util.OutputToTerminal");			
		cfgmissing("Util.afmflag");  				
		cfgmissing("Exchange.Format");  			
		cfgmissing("Exchange.InputFile");			
		cfgmissing("Exchange.Units");   			
		cfgmissing("Exchange.Cutoff");    			
		cfgmissing("Exchange.ChangeSign");  		
		cfgmissing("Exchange.Double_Jij");    		
		cfgmissing("Exchange.CutoffEnergy");    	
		cfgmissing("Exchange.ibtoq");  	
		OutputToTerminal = cfg.lookup("Util.OutputToTerminal");
		afmflag = cfg.lookup("Util.afmflag").c_str();  
		format = cfg.lookup("Exchange.Format").c_str();  
		Jij_filename = cfg.lookup("Exchange.InputFile").c_str();
		Jij_units = cfg.lookup("Exchange.Units").c_str();   
		JijCutoff = cfg.lookup("Exchange.Cutoff");    
		changesign = cfg.lookup("Exchange.ChangeSign").c_str();  
		Jijhalf = cfg.lookup("Exchange.Double_Jij");    
		Jij_min = cfg.lookup("Exchange.CutoffEnergy");    
		ibtoq = cfg.lookup("Exchange.ibtoq");  			

	}
}
