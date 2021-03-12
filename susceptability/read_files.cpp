#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>



int main(){


	double Natoms = 54000;
	double mu_b = 9.2740e-24;
	double mu_s = 3.8664 * mu_b;
	double k_B = 1.3807e-23;
	int stemp = 600;
	int etemp = 605;
	int step = 5;
	std::ofstream output;

	std::stringstream sstr;
	sstr << "/Users/Hirst/Documents/PhD/LLG_code/LLB-master/src/sus_from_mag.txt";
	output.open(sstr.str());


	for (int T = stemp; T < etemp; T+=step){

		// if (T != 40){
		std::stringstream sstr_eq;
		sstr_eq << "/Users/Hirst/Documents/PhD/LLG_code/LLB-master/src/mag_tsteps_1e+06_T_" << T << ".txt";


		std::cout << "Reading file: " << sstr_eq.str() << std::endl;
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
		double avgstart = 50000;
		double avgend = 1000000;
		std::string line;

		while (std::getline(input,line))
		{


			std::istringstream ss(line);
			ss >> a >> b >> c >> d >> e;
			if (a >= avgstart){
				sumx += b;
				sumy += c;
				sumz += d;
				summ += e;
				sumxsqr += b*b;
				sumysqr += c*c;
				sumzsqr += d*d;
				summsqr += e*e;
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

		x = sumx / (avgend-avgstart);
		y = sumy / (avgend-avgstart);
		z = sumz / (avgend-avgstart);
		m = summ / (avgend-avgstart);
		std::cout << m << std::endl;
		xsqr = sumxsqr / (avgend-avgstart);
		ysqr = sumysqr / (avgend-avgstart);
		zsqr = sumzsqr / (avgend-avgstart);
		msqr = summsqr / (avgend-avgstart);

		double susX = 0;
		double susY = 0;
		double susZ = 0;
		double susM = 0;

		susX = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (xsqr - x*x);
		susY = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (ysqr - y*y);
		susZ = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (zsqr - z*z);
		susM = ((Natoms * mu_s) / (k_B * static_cast<double>(T))) * (msqr - m*m);

		output << T << " " << susX << " " << susY << " " << susZ << " " << susM << std::endl;


		// }
	}

	output.close();
	return 0;
}
