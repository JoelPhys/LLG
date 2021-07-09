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
#include "../inc/params1.h"
#include "../inc/error.h"

namespace params {

	libconfig::Config cfg;

	double k_B, mu_b, gamma;
	double dt, Nt, dtau, half_dtau;   
	double lambda, lambdaPrime, mu_s, INVmu_s, thermal_const;
	

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

	int relaxtime;

	//Boundary Conditions
	std::string xbound;
	std::string ybound;
	std::string zbound;

	// Util variables
	std::string afmflag;
	std::string format;
	std::string filepath;

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

	//Intialise Config File ====================================================================================================================================//
	void intitialiseConfig(const char* cfg_filename){

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


		// Global Constants
		k_B = cfg.lookup("PhysicalConsts.BoltzmannConstant");
		mu_b = cfg.lookup("PhysicalConsts.BohrMagneton");
		gamma = cfg.lookup("PhysicalConsts.GyromagneticRatio");

		// Time parameters
		dt = cfg.lookup("Time.SizeOfStep");
		Nt = cfg.lookup("Time.NumberOfSteps");
		relaxtime = cfg.lookup("Time.RelaxationTime");
		
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

		// Print key parameters to log file
		std::cout << "Damping constant = " << lambda << std::endl;
		std::cout << "Magnetic Moment = "<< mu_s << " (mu_b)" << std::endl;
		std::cout << "Uniaxial Anisotropy = [" << dxu << " , " << dyu << " , " << dzu << "] (J)" << std::endl;
		std::cout << "Cubic Anisotropy = [" << dxc << " , " << dyc << " , " << dzc << "] (J)" << std::endl;
		std::cout << "Lattice Parameter = "<< a1 << " (m)" << std::endl;

		// system dimensions
		Lx = cfg.lookup("Geom.UnitCellsInX");
		Ly = cfg.lookup("Geom.UnitCellsInY");
		Lz = cfg.lookup("Geom.UnitCellsInZ");
		Nq = cfg.lookup("Geom.NumberOfSites");
		Nsublat = cfg.lookup("Geom.NumberOfSublat");
		ax = 2;
		ay = 2;
		az = 2;

		// Print system parameters to log file
		std::cout << "Number of unit cells = " << Lx << "x" << Ly << "x" << Lz << std::endl;
		std::cout << "Number of sites in unit cell = " << Nq << std::endl;
		std::cout << "Number of sublattices = " << Nsublat << std::endl;

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
		std::cout << "Boundary Conditions: [" << xbound << " , " << ybound << " , " << zbound << "]" << std::endl;

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
		std::cout << "Lattice Vectors = " << std::endl;
		for (int v = 0; v < 3; v++){

			std::stringstream sstr1;
			sstr1 << "LatVec" << v;
			std::string str1 = sstr1.str();

			libconfig::Setting& setting = cfg.lookup("LatticeVectors");   
			Plat[v][0] = setting[str1.c_str()][0];
			Plat[v][1] = setting[str1.c_str()][1];
			Plat[v][2] = setting[str1.c_str()][2];
			std::cout << Plat[v][0] << " " << Plat[v][1] << " " << Plat[v][2] << std::endl;
		}
		//=======================================================================================================

		// Read Initial Magnetisation Vectors ===================================================================
		initm.resize(Nq);
		std::cout << " Initial Magnestaion Vectors = " << std::endl;
		for (int v = 0; v < Nq; v++){
			
			initm[v].resize(3);
			std::stringstream sstr2;
			sstr2 << "initm" << v;
			std::string str2 = sstr2.str();

			libconfig::Setting& setting = cfg.lookup("InitialMagnetisation");   
			initm[v][0] = setting[str2.c_str()][0];
			initm[v][1] = setting[str2.c_str()][1];
			initm[v][2] = setting[str2.c_str()][2];
			std::cout << initm[v][0] << " " << initm[v][1] << " " << initm[v][2] << std::endl;
		}
		//=======================================================================================================

		start = cfg.lookup("Spinwaves.StartTime");
		afmflag = cfg.lookup("Util.afmflag").c_str();  
		format = cfg.lookup("Exchange.Format").c_str();  
		filepath = cfg.lookup("Util.filepath").c_str();        
		dt_spinwaves = cfg.lookup("Spinwaves.TimeStep");
		Jij_filename = cfg.lookup("Exchange.InputFile").c_str();
		std::cout << "Exchange filename: " << params::Jij_filename << std::endl;        
		Jij_units = cfg.lookup("Exchange.Units").c_str();   
		JijCutoff = cfg.lookup("Exchange.Cutoff");    
		std::cout << "Exhchange Cutoff = " << JijCutoff << std::endl;
		changesign = cfg.lookup("Exchange.ChangeSign").c_str();  
		Jijhalf = cfg.lookup("Exchange.Double_Jij");    
		Jij_min = cfg.lookup("Exchange.CutoffEnergy");    
		std::cout << "Exhchange Energy Minimum = " << Jij_min << " (" << Jij_units << ")" << std::endl;
		ibtoq = cfg.lookup("Exchange.ibtoq");    

	}
	//========================================================================================================================================//


}
