#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>



int main(){


	double Natoms = 108000;
	double mu_b = 9.2740e-24;
	double mu_s = 3.8663 * mu_b;
	double k_B = 1.3807e-23;
	int stemp = 0;
	int etemp = 1401;
	int step = 100;
	int tstep = 10;
	double avgstart = 100000;
	double avgend = 200000;
	std::ofstream output;

	std::stringstream sstr;
	sstr << "/home/sr4871/PhD/materials/mn2au/AFMR/ASD/attempt2susfm.dat";
	output.open(sstr.str());


	for (int T = stemp; T < etemp; T+=step){

		//omit missing files
		// if (T != 600){


			std::stringstream sstr_eq;
			sstr_eq << "/home/sr4871/PhD/materials/mn2au/AFMR/ASD/attempt2/mag_tsteps_200000_T_" << std::setfill('0') << std::setw(4) << T << ".txt";
			int count  = 0;
			std::cout << "Reading file: " << sstr_eq.str() << "\n";
			std::ifstream input(sstr_eq.str());
			if (!input){
				std::cout << "ERROR: Could not open file" << std::endl;
				exit(0);
			}


			double a, b, c, d, e;
			double sumx = 0;
			double sumy = 0;
			double sumz = 0;
			double summ = 0;
			double sumxsqr = 0;
			double sumysqr = 0;
			double sumzsqr = 0;
			double summsqr = 0;
			std::string line;

			while (std::getline(input,line))
			{
				std::istringstream ss(line);
				ss >> a >> b >> c >> d >> e;
				if (a >= avgstart && a < avgend){
					sumx += b;
					sumy += c;
					sumz += d;
					summ += e;
					sumxsqr += b*b;
					sumysqr += c*c;
					sumzsqr += d*d;
					summsqr += e*e;
					count++;
				}
			}

			double x = 0;
			double y = 0;
			double z = 0;
			double m = 0;
			double xsqr = 0;
			double ysqr = 0;
			double zsqr = 0;
			double msqr = 0;

			x = sumx / count;
			y = sumy / count;
			z = sumz / count;
			m = summ / count;
			xsqr = sumxsqr / count;
			ysqr = sumysqr / count;
			zsqr = sumzsqr / count;
			msqr = summsqr / count;

			double susX = 0;
			double susY = 0;
			double susZ = 0;
			double susM = 0;

			susX = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (xsqr - x*x);
			susY = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (ysqr - y*y);
			susZ = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (zsqr - z*z);
			susM = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (msqr - m*m);

			output << T << " " << m << " " << susX << " " << susY << " " << susZ << " " << susM << std::endl;


		// }
	}

	output.close();
	return 0;
}
