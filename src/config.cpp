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
	double lambda, lambdaPrime, mu_s, INVmu_s, thermal_const;
	int relaxtime, outputstep;


	// Uniaxial Anisotropy
	double dxu, dyu, dzu;
	double dzup, dxup, dyup;
	
	//Cubic Anisotropy
	double dxc, dyc, dzc;
	double dzcp, dxcp, dycp;

	int Lx, Ly, Lz, Nq, ax, ay, az, zdimC, Nspins, Nmoments, Nsublat, NmomentsSubLat;
	int Idx, Idy, Idz; // For integer lattice
	double a1, NsitesINV_S, xdim, ydim, zdim, NsitesINV;
	int xdimS, ydimS, zdimS, start;
	double dt_spinwaves;
	double angle;


	//Boundary Conditions
	std::string xbound;
	std::string ybound;
	std::string zbound;

	// Util variables
	std::string afmflag;
	std::string format;
	std::string filepath;
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

	//Intialise Config File ====================================================================================================================================//
	void intitialiseConfig(const char* cfg_filename){

		TITLE("COMPILATION INFO");
		std::cout.width(75); std::cout << std::left << "CPU Compiler: "; std::cout <<CPUCOMP << std::endl;
		std::cout.width(75); std::cout << std::left << "NVCC Compiler: "; std::cout <<GPUCOMP << std::endl;
		std::cout.width(75); std::cout << std::left << "Compile Date and Time: "; std::cout <<__DATE__ << " " << __TIME__ << std::endl;
		std::cout.width(75); std::cout << std::left << "Compiled on Machine: "; std::cout <<HOSTNAME << std::endl;
		if(GITDIRTY!="0")
        {
            std::cout.width(75); std::cout << std::left << "WARNING: Your git build is dirty. Recompile with Git SHA:"; std::cout << GIT_SHA1 << ", Dirty" << std::endl;
        }
        else
        {
            std::cout.width(75); std::cout << std::left << "Git SHA: "; std::cout <<GIT_SHA1 << ", Clean" << std::endl;
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
	//===========================================================================================================================================================//

	// Read Parameters ==========================================================================================================================================//
	void readparams(){

		// Simulation Type
		simtype = cfg.lookup("SimulationType").c_str();

		// Global Constants
		k_B = cfg.lookup("PhysicalConsts.BoltzmannConstant");
		mu_b = cfg.lookup("PhysicalConsts.BohrMagneton");
		gamma = cfg.lookup("PhysicalConsts.GyromagneticRatio");

		// Time parameters
		dt = cfg.lookup("Time.SizeOfStep");
		Nt = cfg.lookup("Time.NumberOfSteps");
		relaxtime = cfg.lookup("Time.RelaxationTime");
		outputstep = cfg.lookup("Time.OutputStep");
		
		// Reduced Time variables
		dtau = gamma * dt;
		half_dtau = 0.5 * dtau;   

		// Material Constants
		lambda = cfg.lookup("MaterialConsts.lambda");
		mu_s = cfg.lookup("MaterialConsts.mu_s");
		a1 = cfg.lookup("MaterialConsts.a");
    	mu_s *= mu_b;
    	INVmu_s = 1 / mu_s;

		// Uniaxial Anisotropy
		dxu = cfg.lookup("Uniaxial_Anisotropy.d_x");
		dyu = cfg.lookup("Uniaxial_Anisotropy.d_y");
		dzu = cfg.lookup("Uniaxial_Anisotropy.d_z");
		dxup = 2 * ( dxu / mu_s );
		dyup = 2 * ( dyu / mu_s );	
		dzup = 2 * ( dzu / mu_s );

		// Cubic Anisotropy
		dxc = cfg.lookup("Cubic_Anisotropy.d_x");
		dyc = cfg.lookup("Cubic_Anisotropy.d_y");
		dzc = cfg.lookup("Cubic_Anisotropy.d_z");
		dxcp = 2 * ( dxc / mu_s );
		dycp = 2 * ( dyc / mu_s );	
		dzcp = 2 * ( dzc / mu_s );

		// Thermal parameters
		thermal_const = sqrt( (2 * lambda * k_B)  / (mu_s * dtau) );
		lambdaPrime = 1 / (1+(lambda*lambda));


		// system dimensions
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
		NmomentsSubLat = Nspins / Nsublat;
		NsitesINV_S = 1/(Lx*Ly*Lz);
		xdim = ax*Lx;
		ydim = ay*Ly;
		zdim = az*Lz;
		NsitesINV = 1/(xdim*ydim*zdim);
		zdimC = zdim/2+1;

		//Boundary Conditions
		xbound = cfg.lookup("Geom.BoundaryConditionsX").c_str(); 
		ybound = cfg.lookup("Geom.BoundaryConditionsY").c_str(); 
		zbound = cfg.lookup("Geom.BoundaryConditionsZ").c_str(); 

		//Temparature
		temptype = cfg.lookup("Temperature.method").c_str();
		ttm_start = cfg.lookup("Temperature.ttm_start");

		// Print key parameters to log file
		TITLE("MATERIAL CONSTANTS");
		std::cout.width(75); std::cout << std::left << "Damping constant:"; std::cout <<lambda << std::endl;
		std::cout.width(75); std::cout << std::left << "Magnetic Moment:"<< mu_s << " (mu_b)"; std::cout <<std::endl;
		std::cout.width(75); std::cout << std::left << "Uniaxial Anisotropy:"; std::cout << "[" << dxu << " , " << dyu << " , " << dzu << "] (J)" << std::endl;
		std::cout.width(75); std::cout << std::left << "Cubic Anisotropy:"; std::cout << "[" << dxc << " , " << dyc << " , " << dzc << "] (J)" << std::endl;
		std::cout.width(75); std::cout << std::left << "Lattice Parameter:"; std::cout << a1 << " (m)" << std::endl;
		std::cout.width(75); std::cout << std::left << "Timestep:"; std::cout << dt << " (s)" << std::endl;
		std::cout.width(75); std::cout << std::left << "Number of timesteps:"; std::cout << Nt << std::endl;
		std::cout.width(75); std::cout << std::left << "Outputting every "; std::cout << outputstep << " timesteps" << std::endl;
		std::cout.width(75); std::cout << std::left << "Outputting to terminal: "; std::cout << OutputToTerminal << std::endl;
		std::cout.width(75); std::cout << std::left << "Number of unit cells:"; std::cout << Lx << "x" << Ly << "x" << Lz << std::endl;
		std::cout.width(75); std::cout << std::left << "Number of sites in unit cell:"; std::cout << Nq << std::endl;
		std::cout.width(75); std::cout << std::left << "Number of sublattices:"; std::cout << Nsublat << std::endl;
		std::cout.width(75); std::cout << std::left << "Boundary Conditions:"; std::cout << "[" << xbound << " , " << ybound << " , " << zbound << "]" << std::endl;
		std::cout.width(75); std::cout << std::left << "Temperature method: "; std::cout << temptype << std::endl;
		std::cout.width(75); std::cout << std::left << "two temperature model start time: "; std::cout << ttm_start << std::endl; 
		

		//Read Site positions ==============================================================================
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
		Isites.resize(Nq);
		for (int s = 0; s < Nq; s++){
			Isites[s].resize(3);

			std::stringstream sstr;
			sstr << "ISite" << s;
			std::string str = sstr.str();

			libconfig::Setting& setting = cfg.lookup("IntegerSites");  
			Isites[s][0] = setting[str.c_str()][0];
			Isites[s][1] = setting[str.c_str()][1];
			Isites[s][2] = setting[str.c_str()][2];
		}
		//Read integer lattice spacing
		Idx = cfg.lookup("IntegerSites.Idx");
		Idy = cfg.lookup("IntegerSites.Idy");
		Idz = cfg.lookup("IntegerSites.Idz");

		//=======================================================================================================

		// Read Lattice Vectors =================================================================================
		for (int v = 0; v < 3; v++){

			std::stringstream sstr1;
			sstr1 << "LatVec" << v;
			std::string str1 = sstr1.str();

			libconfig::Setting& setting = cfg.lookup("LatticeVectors");   
			Plat[v][0] = setting[str1.c_str()][0];
			Plat[v][1] = setting[str1.c_str()][1];
			Plat[v][2] = setting[str1.c_str()][2];
			std::cout.width(75); std::cout << std::left << "Lattice Vectors:";
			std::cout << std::fixed << std::setprecision(5) << Plat[v][0] << " " << Plat[v][1] << " " << Plat[v][2] << std::endl;
		}
		//=======================================================================================================

		// Read Initial Magnetisation Vectors ===================================================================
		initm.resize(Nq);
		
		for (int v = 0; v < Nq; v++){
			std::cout.width(75); std::cout << std::left << "Initial Magnestaion Vectors:";
			initm[v].resize(3);
			std::stringstream sstr2;
			sstr2 << "initm" << v;
			std::string str2 = sstr2.str();

			libconfig::Setting& setting = cfg.lookup("InitialMagnetisation");   
			initm[v][0] = setting[str2.c_str()][0];
			initm[v][1] = setting[str2.c_str()][1];
			initm[v][2] = setting[str2.c_str()][2];
			std::cout << std::fixed << std::setprecision(5) << initm[v][0] << " " << initm[v][1] << " " << initm[v][2] << std::endl;
		}
		//=======================================================================================================


		OutputToTerminal = cfg.lookup("Util.OutputToTerminal");
		start = cfg.lookup("Spinwaves.StartTime");
		afmflag = cfg.lookup("Util.afmflag").c_str();  
		format = cfg.lookup("Exchange.Format").c_str();  
		filepath = cfg.lookup("Util.filepath").c_str();        
		dt_spinwaves = cfg.lookup("Spinwaves.TimeStep");
		Jij_filename = cfg.lookup("Exchange.InputFile").c_str();
		Jij_units = cfg.lookup("Exchange.Units").c_str();   
		JijCutoff = cfg.lookup("Exchange.Cutoff");    
		changesign = cfg.lookup("Exchange.ChangeSign").c_str();  
		Jijhalf = cfg.lookup("Exchange.Double_Jij");    
		Jij_min = cfg.lookup("Exchange.CutoffEnergy");    
		ibtoq = cfg.lookup("Exchange.ibtoq");  

		TITLE("EXCHANGE FILE INFO");
		std::cout.width(75); std::cout << std::left << "Exchange filename: "; std::cout <<params::Jij_filename << std::endl;        
		std::cout.width(75); std::cout << std::left << "Exhchange Cutoff:"; std::cout <<JijCutoff << std::endl;
		std::cout.width(75); std::cout << std::left << "Exhchange Energy Minimum:"; std::cout <<Jij_min << " (" << Jij_units << ")" << std::endl;
		if (Jijhalf == true) {std::cout.width(75); std::cout << std::left << "Have Jij values been doubled"; std::cout << "Yes" << std::endl;}
		else if (Jijhalf == false) {std::cout.width(75); std::cout << std::left << "Have Jij values been doubled"; std::cout << "No" << std::endl;}
 

	}
	//========================================================================================================================================//


}
