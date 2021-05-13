#include <cmath>
#include <fstream>
#include "../inc/mathfuncs.h"
#include "../inc/params1.h"
#include "../inc/util.h"
#include "../inc/NeighbourList.h"
#include "../inc/array.h"
#include "../inc/array2d.h"

int modfunc(int L, int x){
	int y = ((x % L) + L) % L;
	return y;
}

void CrossP(double vect_A[], double vect_B[], double crossP[]){
	crossP[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
	crossP[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
	crossP[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}

void DotP(double vect_A1[], double vect_B1[], double dotP){
	dotP = vect_A1[0] * vect_B1[0] + vect_A1[1] * vect_B1[1] + vect_A1[2] * vect_B1[2]; 
}

void Inverse3x3(double m[3][3], double minv[3][3]){
	double det = m[0][ 0] * (m[1][ 1] * m[2][ 2] - m[2][ 1] * m[1][ 2]) -
		m[0][ 1] * (m[1][ 0] * m[2][ 2] - m[1][ 2] * m[2][ 0]) +
		m[0][ 2] * (m[1][ 0] * m[2][ 1] - m[1][ 1] * m[2][ 0]);

	double invdet = 1 / det;

	minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
	minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
	minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
	minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
	minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
	minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
	minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
	minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
	minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

}


//rotation matrix
double xu[3] = {0.0,0.0,1.0};
double C;
double D[3];
double Stest[3];
double X[3][3];
double Xsqr[3][3];
double Mnewb[3];
double norm;
double norm1;
double pref;
double Mnew[3];
double Mu[3];
Array2D<double> Id;
Array2D<double> R;


void IdentityMatrix(){

	Id.resize(3,3);
	Id.IFill(0);
	Id(0,0) = 1;
	Id(1,1) = 1;
	Id(2,2) = 1;

	R.resize(3,3);
	R.IFill(0);
	R(0,0) = 1;
	R(1,1) = 1;
	R(2,2) = 1;

}


void Rotation(){

	//convert degrees to radians
	C = M_PI / 180.0;
	double angle = C * 90;
	double x,y,z;
	//loop through all spins
	for (int i = 0; i < params::Nspins; i++){

		//Find spins in the second sublattice
		if ((modfunc(params::Nq,i) == 0) || (modfunc(params::Nq,i) == 3) || (modfunc(params::Nq,i) == 5) || (modfunc(params::Nq,i) == 6)) {
			
			// Asign to temporary variables to avoid overwriting each component
			x = neigh::Sx1d[i]; 
			y = neigh::Sy1d[i] * cos(angle) - neigh::Sz1d[i] * sin(angle);
			z = neigh::Sy1d[i] * sin(angle) + neigh::Sz1d[i] * cos(angle);
			
			// assign each component to the magnetisation vector
			neigh::Sx1d[i] = x;
			neigh::Sy1d[i] = y;
			neigh::Sz1d[i] = z;
			
			std::cout << neigh::Sx1d[i] << " "; 
			std::cout << neigh::Sy1d[i] << " ";
			std::cout << neigh::Sz1d[i] << std::endl;
		}

	}


	}








	// void Rotation(){

	//             //normalise M
	//             Mnew[0] = util::M(0,0) / params::NmomentsSubLat;
	//             Mnew[1] = util::M(0,1) / params::NmomentsSubLat;
	//             Mnew[2] = util::M(0,2) / params::NmomentsSubLat;

	//             norm = sqrt(Mnew[0] * Mnew[0] + Mnew[1] * Mnew[1] + Mnew[2] * Mnew[2]);

	//             Mu[0] = Mnew[0] / norm;
	//             Mu[1] = Mnew[1] / norm;
	//             Mu[2] = Mnew[2] / norm;

	//             C = Mu[0] * xu[0] + Mu[1] * xu[1] + Mu[2] * xu[2]; 

	//             D[0] = Mu[1] * xu[2] - Mu[2] * xu[1]; 
	//             D[1] = Mu[2] * xu[0] - Mu[0] * xu[2]; 
	//             D[2] = Mu[0] * xu[1] - Mu[1] * xu[0]; 

	//             X[0][0] = 0;
	//             X[0][1] = -1 * D[2];
	//             X[0][2] = D[1];
	//             X[1][0] = D[2];
	//             X[1][1] = 0;
	//             X[1][2] = -1 * D[0];
	//             X[2][0] = -1 * D[1];
	//             X[2][1] = D[0];
	//             X[2][2] = 0;

	//             norm1 = sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]); 

	//             for (int w = 0; w < 3; w++){
	//                 for (int e = 0; e < 3; e++){
	//                     Xsqr[w][e] = 0;
	//                 }
	//             }

	//             for (int w = 0; w < 3; w++){
	//                 for (int e = 0; e < 3; e++){
	//                     for (int r = 0; r < 3; r++){
	//                         Xsqr[w][e] += X[w][r] * X[r][e];
	//                     }
	//                 }
	//             }

	//             for (int w = 0; w < 3; w++){
	//                 for (int e = 0; e < 3; e++){
	//                     R(w,e) = Id(w,e) + X[w][e] + ( (1 - C) / (norm1*norm1) ) * Xsqr[w][e];
	//                 }
	//             }

	//             for (int a = 0; a < params::Nspins; a++){     
	//                     Stest[0] = R(0,0) * neigh::Sx1d[a] + R(0,1) * neigh::Sy1d[a] + R(0,2) * neigh::Sz1d[a];
	//                     Stest[1] = R(1,0) * neigh::Sx1d[a] + R(1,1) * neigh::Sy1d[a] + R(1,2) * neigh::Sz1d[a];
	//                     Stest[2] = R(2,0) * neigh::Sx1d[a] + R(2,1) * neigh::Sy1d[a] + R(2,2) * neigh::Sz1d[a];

	//                     neigh::Sx1d[a] = Stest[0];
	//                     neigh::Sy1d[a] = Stest[1];
	//                     neigh::Sz1d[a] = Stest[2];

	//             }

	//             // std::cout << Stest[0] << " ";
	//             // std::cout << R(1][0] * Mnew[0] + R(1][1] * Mnew[1] + R(1][2] * Mnew[2] << " ";
	//             // std::cout << R(2][0] * Mnew[0] + R(2][1] * Mnew[1] + R(2][2] * Mnew[2] << std::endl;

	// }

