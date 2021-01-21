import numpy as np 
import matplotlib.pyplot as plt
import sys
import math
import os

count  = 0
min = 1255
max = 1750
steps = 5
Nsteps = 99
testingT = np.arange(min,max,steps)
testingX = np.zeros(Nsteps)
testingY = np.zeros(Nsteps)
testingZ = np.zeros(Nsteps)
testingXsqr = np.zeros(Nsteps)
testingYsqr = np.zeros(Nsteps)
testingZsqr = np.zeros(Nsteps)
testingM = np.zeros(Nsteps)


testingXsqr = np.zeros(Nsteps)
testingYsqr = np.zeros(Nsteps)
testingZsqr = np.zeros(Nsteps)
TransX = np.zeros(Nsteps)
ParraX = np.zeros(Nsteps)
SusX = np.zeros(Nsteps)
SusY = np.zeros(Nsteps)
SusZ = np.zeros(Nsteps)

Natoms = 108000
mu_b = 9.2740e-24
mu_s = 3.8663 * mu_b
k_B = 1.3807e-23

for i in range(min, max, steps):
    string = "/Users/Hirst/Documents/PhD/LLG_code/Mn2Au/results/AFM/run2-(38x38x38)/rsj-gga-af-p3-v2/txtfiles/af_gga_Nsteps_100000_T_" + str(i) + ".txt"
    data = np.loadtxt(string, usecols =(0,1,2,3,4))

    length = len(data[:,1])
    Num_x = data[length//2:length,1]
    Num_y = data[length//2:length,2]
    Num_z = data[length//2:length,3]
    Num_x_sqr = data[length//2:length,1] * data[length//2:length,1]
    Num_y_sqr = data[length//2:length,2] * data[length//2:length,2]
    Num_z_sqr = data[length//2:length,3] * data[length//2:length,3]
    Num_M = data[length//2:length,4]
    testingX[count] = np.average(Num_x)
    testingY[count] = np.average(Num_y)
    testingZ[count] = np.average(Num_z)
    testingXsqr[count] = np.average(Num_x_sqr)
    testingYsqr[count] = np.average(Num_y_sqr)
    testingZsqr[count] = np.average(Num_z_sqr)
    testingM[count] = np.average(Num_M)

    if testingT[count] == 0:
        ParraX[count] = 0
        TransX[count] = 0
        SusX[count] = 0
        SusY[count] = 0
        SusZ[count] = 0
    else:	
        ParraX[count] = ((Natoms * mu_b) / (k_B * testingT[count])) * (testingZsqr[count] - testingZ[count]*testingZ[count])
        TransX[count] = ((Natoms * mu_b) / (2 * k_B * testingT[count])) * (testingXsqr[count] + testingYsqr[count] - testingX[count]*testingX[count] - testingY[count]*testingY[count])
        SusX[count] = ((Natoms * mu_b) / (k_B * testingT[count])) * (testingXsqr[count] - testingX[count]*testingX[count])
        SusY[count] = ((Natoms * mu_b) / (k_B * testingT[count])) * (testingYsqr[count] - testingY[count]*testingY[count])
        SusZ[count] = ((Natoms * mu_b) / (k_B * testingT[count])) * (testingZsqr[count] - testingZ[count]*testingZ[count])

    print(i)
    print(testingT[count])
    print(testingZ[count])
    print(testingZsqr[count])
    print(testingZ[count]*testingZ[count])
    print(testingZsqr[count] - testingZ[count]*testingZ[count])
    print(ParraX[count])
    count = count+1


# np.savetxt("/Users/Hirst/Documents/PhD/LLG_code/Mn2Au/results/AFM/run3 (30x30x30) 4 nearest neighbours/lda/Mag_v_temp.txt", np.c_[testingT, testingX, testingY, testingZ, testingM], fmt='%s')
np.savetxt('/Users/Hirst/Documents/PhD/LLG_code/Mn2Au/results/AFM/run2-(38x38x38)/rsj-gga-af-p3-v2/TEST_Susceptability-1250to1800.txt', np.c_[testingT, ParraX, TransX], fmt='%s')
