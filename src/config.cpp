#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <sys/stat.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include "../inc/mathfuncs.h"
#include "../inc/array2d.h"
#include "../inc/array.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/defines.h"

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
	double dxu, dyu, dzu;
	double dzup, dxup, dyup;

	//Cubic Anisotropy
	double dzc, dzcp;

	int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat;
	int Idx, Idy, Idz; // For integer lattice
	double a1, b1, c1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;
	int xdimS, ydimS, zdimS, start;
	double dt_spinwaves;
	double sg_spinwaves;
	double angle;
	std::vector<int> sublat_sites;
	std::vector<int> NmomentsSubLat;


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
	std::vector< std::vector<double> > Isites;
	std::vector< std::vector<double> > initm;

	//Temperature
	std::string temptype;
	double ttm_start;
	double temp_gradient;

	// Output Lattice
	bool OutputLattice = false;		
	int OutputLatticeStep = 10000000;
	int OutputLatticeStart = 10000000;	
	std::string OutputLatticeFilepath;		

	void banner(){
		std::cout << std::endl;
		std::cout << "                                            ███    █████████    █████████  ██████████                     " << std::endl;
		std::cout << "                                           ░░░    ███░░░░░███  ███░░░░░███░░███░░░░███                    " << std::endl;
		std::cout << "                                           █████ ░███    ░███ ░███    ░░░  ░███   ░░███                   " << std::endl;
		std::cout << "                                          ░░███  ░███████████ ░░█████████  ░███    ░███                   " << std::endl;
		std::cout << "                                           ░███  ░███░░░░░███  ░░░░░░░░███ ░███    ░███                   " << std::endl;
		std::cout << "                                           ░███  ░███    ░███  ███    ░███ ░███    ███                    " << std::endl;
		std::cout << "                                           ░███  █████   █████░░█████████  ██████████                     " << std::endl;
		std::cout << "                                           ░███ ░░░░░   ░░░░░  ░░░░░░░░░  ░░░░░░░░░░                      " << std::endl;
		std::cout << "                                       ███ ░███                                                           " << std::endl;
		std::cout << "                                      ░░██████                                                            " << std::endl;
		std::cout << "                                       ░░░░░░                                                             " << std::endl;		
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

		// Gilbert Damping
		cfgmissing("MaterialConsts.lambda");				
		// lambda = cfg.lookup("MaterialConsts.lambda");

		// Magnetic Moment
		// cfgmissing("MaterialConsts.mu_s");
		// mu_s = cfg.lookup("MaterialConsts.mu_s");


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

		//Temparature
		cfgmissing("Temperature.method");					
		cfgmissing("Temperature.ttm_start");				
		temptype = cfg.lookup("Temperature.method").c_str();
		ttm_start = cfg.lookup("Temperature.ttm_start");
		if ((temptype == "uniform_gradient") && (cfg.exists("Temperature.gradient") == 0)){
			std::cout << "ERROR: Method is uniform_gradient but no gradient value has been provided" << std::endl;
			std::cout << "EXITING SIMULATION" << std::endl;
			exit(0);
		}
		else if ((temptype == "uniform_gradient")){
			double temp_gradient_km = cfg.lookup("Temperature.gradient");
			INFO_OUT("Temperature Gradient:", temp_gradient_km << "[K/m]");
			temp_gradient = temp_gradient_km * a1;
			INFO_OUT("Temperature Gradient:", temp_gradient << "[K / unit cell]");
		}

		// Print key parameters to log file
		TITLE("MATERIAL CONSTANTS");
		// INFO_OUT("Damping constant:",lambda);
		// INFO_OUT("Magnetic Moment:", mu_s << " (mu_b)");
		INFO_OUT("Lattice Parameter, a:", a1 << " (m)");
		INFO_OUT("Lattice Parameter, b:", b1 << " (m)");
		INFO_OUT("Lattice Parameter, c:", c1 << " (m)");
		INFO_OUT("Timestep:", dt << " (s)");
		INFO_OUT("Number of timesteps:", Nt);
		INFO_OUT("Outputting every ", outputstep << " timesteps");
		INFO_OUT("Outputting to terminal: ", OutputToTerminal);
		INFO_OUT("Number of unit cells:", Lx << "x" << Ly << "x" << Lz);
		INFO_OUT("Number of sites in unit cell:", Nq);
		INFO_OUT("Number of sublattices:", Nsublat);
		INFO_OUT("Boundary Conditions:", "[" << xbound << " , " << ybound << " , " << zbound << "]");
		INFO_OUT("Temperature method: ", temptype);
		INFO_OUT("two temperature model start time: ", ttm_start);


		//=======================================================================================================
		//Read Magnetic Moments =================================================================================
		//=======================================================================================================

		mu_s.resize(Nq);
		INVmu_s.resize(Nq);
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
		// Uniaxial Anisotropy ==================================================================================
		//=======================================================================================================
		cfgmissing("Uniaxial_Anisotropy.d_x");			
		cfgmissing("Uniaxial_Anisotropy.d_y");			
		cfgmissing("Uniaxial_Anisotropy.d_z");
		dxu = cfg.lookup("Uniaxial_Anisotropy.d_x");
		dyu = cfg.lookup("Uniaxial_Anisotropy.d_y");
		dzu = cfg.lookup("Uniaxial_Anisotropy.d_z");			
		dxup = 2 * ( dxu / mu_s[0] );
		dyup = 2 * ( dyu / mu_s[0] );	
		dzup = 2 * ( dzu / mu_s[0] );

		INFO_OUT("Uniaxial Anisotropy Energy:", "[" << dxu << " , " << dyu << " , " << dzu << "] [J]");
		INFO_OUT("Uniaxial Anisotropy Field:", "[" << dxup << " , " << dyup << " , " << dzup << "] [J]");
		INFO_OUT("Cubic Anisotropy:", dzc << " [J]");



		//=======================================================================================================
		//Read Damping Consts ===================================================================================
		//=======================================================================================================

		lambda.resize(Nq);
		lambdaPrime.resize(Nq);
		thermal_const.resize(Nq);
		std::string lmtext;
		for (int i = 0; i < Nq; i++){

			// Read from config file
			libconfig::Setting& setting = cfg.lookup("MaterialConsts");  
			lambda[i] = setting["lambda"][i];

			// Print information to log file
			lmtext = "Damping Constant for site  " + std::to_string(i) + ":";
			INFO_OUT(lmtext, lambda[i] << " [arb.]");

			// Thermal parameters
			thermal_const[i] = sqrt( (2 * lambda[i] * k_B)  / (mu_s[i] * dtau) );
			lambdaPrime[i] = 1 / (1+(lambda[i]*lambda[i]));

		}

		//=======================================================================================================
		//Read Site positions ===================================================================================
		//=======================================================================================================

		sites.resize(Nq);
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

		Isites.resize(Nq);
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
		//Read integer lattice spacing
		Idx = cfg.lookup("IntegerSites.Idx");
		Idy = cfg.lookup("IntegerSites.Idy");
		Idz = cfg.lookup("IntegerSites.Idz");

		//=======================================================================================================
		// Read Lattice Vectors and constants ===================================================================
		//=======================================================================================================

		// lattice constants
		cfgmissing("LatticeVectors.a");
		cfgmissing("LatticeVectors.b");
		cfgmissing("LatticeVectors.c");	
		a1 = cfg.lookup("LatticeVectors.a");
		b1 = cfg.lookup("LatticeVectors.b");
		c1 = cfg.lookup("LatticeVectors.c");

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

		//=======================================================================================================
		// Read Initial Magnetisation Vectors ===================================================================
		//=======================================================================================================

		initm.resize(Nq);

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

		//=======================================================================================================
		// Read Sublat vector ===================================================================================
		//=======================================================================================================
		sublat_sites.resize(Nq);
		NmomentsSubLat.resize(Nsublat);

		for (int v = 0; v < Nq; v++){
			libconfig::Setting& setting = cfg.lookup("Geom");   
			sublat_sites[v] = setting["sublat_sites"][v];

			//count number of sites within each sublattice
			NmomentsSubLat[sublat_sites[v]] += 1*Lx*Ly*Lz;

		}

		std::cout << NmomentsSubLat[0] << " " << NmomentsSubLat[1] << " " << NmomentsSubLat[2] << std::endl;

		//Cubic Anisotropy
		cfgmissing("Cubic_Anisotropy.d_c");
		dzc = cfg.lookup("Cubic_Anisotropy.d_c");
		dzcp = 2 * ( dzc / mu_s[0] );

		//=======================================================================================================
		// Check filepaths exist ================================================================================
		//=======================================================================================================

		cfgmissing("Spinwaves.filepath");        		
		cfgmissing("Util.filepath");
		filepath = cfg.lookup("Util.filepath").c_str();      
		filepath_sw = cfg.lookup("Spinwaves.filepath").c_str();   

		struct stat buffer;
		if (stat(filepath.c_str(), &buffer) != 0) {
				std::cout << "ERROR: Output directory does not exist!" << std::endl;
			exit(0);
		}

		if (stat(filepath_sw.c_str(), &buffer) != 0) {
    		std::cout << "ERROR: Spinwaves output directory does not exist!" << std::endl;;
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
			OutputLatticeStep = cfg.lookup("Util.OutputLatticeStep");  
			OutputLatticeStart = cfg.lookup("Util.OutputLatticeStart");  
			OutputLatticeFilepath = cfg.lookup("Util.OutputLatticeFilepath").c_str();  
		}


		//=======================================================================================================
		// Rest of parameters ===================================================================================
		//=======================================================================================================

		cfgmissing("Util.OutputToTerminal");		
		cfgmissing("Spinwaves.StartTime");			
		cfgmissing("Util.afmflag");  				
		cfgmissing("Exchange.Format");  			
		cfgmissing("Spinwaves.TimeStep");
		cfgmissing("Spinwaves.smoothing");			
		cfgmissing("Exchange.InputFile");			
		cfgmissing("Exchange.Units");   			
		cfgmissing("Exchange.Cutoff");    			
		cfgmissing("Exchange.ChangeSign");  		
		cfgmissing("Exchange.Double_Jij");    		
		cfgmissing("Exchange.CutoffEnergy");    	
		cfgmissing("Exchange.ibtoq");  	
		OutputToTerminal = cfg.lookup("Util.OutputToTerminal");
		start = cfg.lookup("Spinwaves.StartTime");
		afmflag = cfg.lookup("Util.afmflag").c_str();  
		format = cfg.lookup("Exchange.Format").c_str();  
		dt_spinwaves = cfg.lookup("Spinwaves.TimeStep");
		sg_spinwaves = cfg.lookup("Spinwaves.smoothing");
		Jij_filename = cfg.lookup("Exchange.InputFile").c_str();
		Jij_units = cfg.lookup("Exchange.Units").c_str();   
		JijCutoff = cfg.lookup("Exchange.Cutoff");    
		changesign = cfg.lookup("Exchange.ChangeSign").c_str();  
		Jijhalf = cfg.lookup("Exchange.Double_Jij");    
		Jij_min = cfg.lookup("Exchange.CutoffEnergy");    
		ibtoq = cfg.lookup("Exchange.ibtoq");  			

		TITLE("EXCHANGE FILE INFO");
		INFO_OUT("Exchange filename: ", params::Jij_filename);        
		INFO_OUT("Exhchange Cutoff:", JijCutoff);
		INFO_OUT("Exhchange Energy Minimum:", Jij_min << " (" << Jij_units << ")");
		if (Jijhalf == true) {INFO_OUT("Have Jij values been doubled", "Yes");}
		else if (Jijhalf == false) {INFO_OUT("Have Jij values been doubled", "No");}


	}

}
