import numpy as np 
import matplotlib.pyplot as plt
import sys
import math
import os

count  = 0
min = 0
max = 1800
steps = 5
Nsteps = 360
Avg_start = 10000

T = np.zeros(Nsteps)
X = np.zeros(Nsteps)
Y = np.zeros(Nsteps)
Z = np.zeros(Nsteps)
Xsqr = np.zeros(Nsteps)
Ysqr = np.zeros(Nsteps)
Zsqr = np.zeros(Nsteps)
M = np.zeros(Nsteps)
Msum = np.zeros(Nsteps)
Msqr = np.zeros(Nsteps)
SusM = np.zeros(Nsteps)
TransX = np.zeros(Nsteps)
ParraX = np.zeros(Nsteps)
SusX = np.zeros(Nsteps)
SusY = np.zeros(Nsteps)
SusZ = np.zeros(Nsteps)
SusT = np.zeros(Nsteps)

Natoms = 54000
mu_b = 9.2740e-24
mu_s = 3.8663 * mu_b
k_B = 1.3807e-23

for i in range(min, max, steps):

    # omit missing files 
    if (i != 40):
        string = "/Volumes/ExternalHDD/Joel/Mn2Au_AFM_GGA/mag_tsteps_gga_1e+06_T_" + str(i) + ".txt"

        data = np.loadtxt(string, usecols =(0,1,2,3,4), max_rows=1000000)

        length = 1000000
        Num_x = data[Avg_start:length,1]
        Num_y = data[Avg_start:length,2]
        Num_z = data[Avg_start:length,3]
        Num_x_sqr = data[Avg_start:length,1] * data[Avg_start:length,1]
        Num_y_sqr = data[Avg_start:length,2] * data[Avg_start:length,2]
        Num_z_sqr = data[Avg_start:length,3] * data[Avg_start:length,3]
        X[count] = np.average(Num_x)
        Y[count] = np.average(Num_y)
        Z[count] = np.average(Num_z)
        Xsqr[count] = np.average(Num_x_sqr)
        Ysqr[count] = np.average(Num_y_sqr)
        Zsqr[count] = np.average(Num_z_sqr)

        Num_M = data[Avg_start:length,4]
        Num_M_sqr = data[Avg_start:length,4] * data[Avg_start:length,4]

        print(X[count], Xsqr[count])


        M[count] = np.average(Num_M)
        Msqr[count] = np.average(Num_M_sqr)
        T[count] = i

        if T[count] == 0:
            ParraX[count] = 0
            TransX[count] = 0
            SusX[count] = 0
            SusY[count] = 0
            SusZ[count] = 0
        else:	
            SusM[count] = ((Natoms * mu_s) / (k_B * T[count])) * (Msqr[count] - M[count]*M[count])
            SusX[count] = ((Natoms * mu_s) / (k_B * T[count])) * (Xsqr[count] - X[count]*X[count])
            SusY[count] = ((Natoms * mu_s) / (k_B * T[count])) * (Ysqr[count] - Y[count]*Y[count])
            SusZ[count] = ((Natoms * mu_s) / (k_B * T[count])) * (Zsqr[count] - Z[count]*Z[count])
            SusT[count] = ((Natoms * mu_s) / (2 * k_B * T[count])) * (Zsqr[count] + Ysqr[count] - Y[count]*Y[count] - Z[count]*Z[count])



        print(T[count], SusX[count], SusY[count], SusZ[count])
        count = count+1

# np.savetxt("/Users/Hirst/flucttest.txt", np.c_[T, SusM], fmt='%s')
np.savetxt("/Users/Hirst/Documents/testing/MvT2.txt", np.c_[T, X, Y, Z, M], fmt='%s')
np.savetxt("/Users/Hirst/Documents/testing/Sus2.txt", np.c_[T, SusX, SusY, SusZ], fmt='%s')
