#include <fstream>
#include <vector>
#include <iostream>
#include "../inc/config.h"
#include "../inc/libconfig.h++"
#include "../inc/util.h"


namespace thermal {


    // Define global variables
    std::vector<double> ttime;     // time
    std::vector<double> tz;        // z plane
    std::vector<double> Te;        // electron temp
    std::vector<double> P_it;        // pump power
    std::vector<double> Tp;        // phonon temp

     //two temperature model variables
    double gamma_e=1e3;                //gamma_e defines the electron specific heat through, e = gamma_e * T_e. [J/m^3/K^2]
    double Cp=3e6;                       //Specific heat of phonons. [J/m^3/K]
    double kappa_0=11.0;                 //kappa_0 defines the thermal heat conductivity (kappa) through, kappa = kappa_0 * T_e/T_p [J/m/K/s]
    double delta=20.0e-9;                //Penetration depth of laser. [m]
    double Gep=10e17;                    //Electron-phonon coupling [ J/m^3/s/K]
    double P_0=2.5e21;                   //Pump fluence prefactor, P_0. P(z,t)=P_0*exp(-((t-t0)/tau)**2)*exp(-z/delta) [ J/m^3/s]
    double t0=5000e-15;                   //Pump temporal offset [s]
    double tau=4500e-16;                   //Pump temporal full width half max [s]
    double Nz=100;                       //number of unit cells in z-direction (assumed uniform heating perpendicular [unit cells in z]
    double dz=0.3e-9;                    //lattice constant (or difference between planes) [m]

    // Arrays

    // int nsteps = 0;
    // int nz = 0;

    void ttm(double time){

        // loop over layes

        P_it[i]=P_0*exp(-((time-t0)/tau)*((time-t0)/tau));
            Tep1=Te[0] + (dt/(gamma_e*Te[0]))*(Gep*(Tp[0]-Te[0]) + P_it[0] + kappa_0*( (Te[0]/Tp[0]) * 2.0*(Te[1]-Te[0])*oneOvrdzdz));
            Tpp1=Tp[0]+(dt*Gep/Cp)*(Te[0]-Tp[0]);
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


    void ReadThermalFile(){
		// Define variables
        std::cout << "READING THERMAL FILE " << std::endl;
		std::ifstream input(params::temp_filename);
		int a, b; 
        double c, d, e;
        int count = 0;

		// Read file
        if (!input){
			//Output Error if file can't be read. Quit program.
            std::cout << "ERROR: Could not open file " << params::temp_filename << std::endl;
            exit(0);
        }
		else {
			while (input >> a >> b >> c >> d >> e)
            {
				ttime.push_back(a);
				tz.push_back(b);
				pp.push_back(c);
				te.push_back(d);
				tp.push_back(e);
                count++;
            }
        	std::cout << "temp input file has been read" << std::endl;
		}

        // number of timesteps
        int nsteps = ttime.back();
        int nz = tz.back();
	}

}