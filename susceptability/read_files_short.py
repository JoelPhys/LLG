import matplotlib.pyplot as plt
import numpy as np
import os

#Simulation values
Natoms = 54000
mu_b = 9.2740e-24
mu_s = 3.8663 * mu_b
k_B = 1.3807e-23

# Values for Loop
min = 0
max = 795
Nsteps = 159
steps = 5
count = 0

#Arrays
SusX = np.zeros(Nsteps)
SusY = np.zeros(Nsteps)
SusZ = np.zeros(Nsteps)
Temp = np.zeros(Nsteps)
Mmag = np.zeros(Nsteps)

for i in range(min, max, steps):
    string = "/Users/Hirst/Documents/PhD/LLG_code/Mn2Au/results/AFM/run5/txtfiles/mag_tsteps_lda_1e+06_T_" + str(i) + ".txt"
    # data = np.loadtxt(string, usecols =(0,1))
    f = open(string,"r")
    lines  = f.readlines()
    Nt = float(lines[2])
    Mmag1 = float(lines[0].split()[6])

    x = float(lines[0].split()[0])
    x_sqr = float(lines[1].split()[0])
    y = float(lines[0].split()[1])
    y_sqr = float(lines[1].split()[1])
    z = float(lines[0].split()[2])
    z_sqr = float(lines[1].split()[2])


    xval = x / Nt
    xval_sqr = x_sqr / Nt
    yval = y / Nt
    yval_sqr = y_sqr / Nt
    zval = z / Nt
    zval_sqr = z_sqr / Nt

    Mmag[count] = Mmag1 / Nt

    Temp[count] = float(i)

    print(Temp[count], y, y_sqr, Nt)

    if Temp[count] == 0:
        SusX[count] = 0
        SusY[count] = 0
        SusZ[count] = 0
    else:
        SusX[count] = ((Natoms * mu_b) / (k_B * Temp[count])) * (xval_sqr - xval*xval)
        SusY[count] = ((Natoms * mu_b) / (k_B * Temp[count])) * (yval_sqr - yval*yval)
        SusZ[count] = ((Natoms * mu_b) / (k_B * Temp[count])) * (zval_sqr - zval*zval)
    
    count = count+1


# plt.plot(Temp,Mmag)
# plt.plot(Temp[40:],SusX[40:])
plt.plot(Temp[1:],SusY[1:])
# plt.plot(Temp[1:],SusZ[1:])
plt.ylabel('Reduced Susceptibility (1 / T)')
plt.xlabel('Temperature (K)')
plt.title('$M_z$')


plt.show()