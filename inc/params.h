// #ifndef __PARAMS_H__
// #define __PARAMS_H__

// #include <vector>
// #include <string>
// #include <cmath>


// namespace params {

// 	// PHYSICAL CONSTANTS
// 	const double k_B = 1.3807e-23;
// 	const double mu_b = 9.2740e-24;
// 	const double gamma = 1.76e11;

// 	// SIMULATION PARAMETERS
// 	const double dt = 1e-16; 
// 	const double Nt = 10000; 
// 	const double dtau = gamma * dt;
// 	const double half_dtau = dtau * 0.5;    

// 	// MATERIAL CONSTANTS
// 	const double lambda = 1;
// 	const double lambdaPrime = 1 / (1+(lambda*lambda));
// 	const double mu_s = 3.8663 * mu_b;
// 	const double INVmu_s = 1 / mu_s;
// 	const double d_z = 0;
// 	const double thermal_const = sqrt( (2 * lambda * k_B)  / (mu_s * dtau) );
// 	const double d_z_prime = 2 * ( d_z / mu_s );


// 	// EXTERNAL FIELD
// 	const double H_app[3] = {0,0,0};

// 	// SYSTEM DIMENSIONS
// 	const int Lx = 30;
// 	const int Ly = 30;
// 	const int Lz = 30;
// 	const int Nq = 4;
// 	const int Nsublat = 2;
// 	const double a1 = 0.287e-9;
// 	const int ax = 2;
// 	const int ay = 2;
// 	const int az = 2;
// 	const int Nspins = Nq*Lx*Ly*Lz;
// 	const int Nmoments = (4*Lx*Ly*Lz);
// 	const int NmomentsSubLat = (2*Lx*Ly*Lz);
// 	const double NsitesINV_S = 1/(Lx*Ly*Lz);
// 	const double xdim = ax*Lx;
// 	const double ydim = ay*Ly;
// 	const double zdim = az*Lz;
// 	const double NsitesINV = 1/(xdim*ydim*zdim);
// 	const int zdimC = zdim/2+1;


// 	// playing with spinwaves
// 	int xdimS = Lx;
// 	int ydimS = Ly;
// 	int zdimS = Lz/2 + 1;
// 	double dt_spinwaves = 5e-15;
// 	int start = 50000;

// 	// Jij SETTINGS
// 	std::string Jij_filename = "../Mn2Au/modified_input/FullInteraction/rsj-gga-af-p3.mn2au";
// 	bool Jijhalf = true;
// 	std::string Jij_units = "mRy";
// 	bool JijCutoff = true;
// 	const double Jij_min = 0.01;
// 	const int ibtoq = -3;
// 	bool changesign = true;

// 	//atom sites
// 	double Plat[3][3] = {{1.0, 0.0, 0.0},
// 						 {0.0, 1.0, 0.0},
// 						 {0.0, 0.0, 2.563942}};

// 	double sites[Nq][3] = {{0.500000, 0.500000, 0.429204},
// 						  {0.000000, 0.000000, 0.853767},
// 						  {0.000000, 0.000000, 1.711175},
// 						  {0.500000, 0.500000, 2.134738}};

// 	double PlatINV[3][3];
// }


// #endif
