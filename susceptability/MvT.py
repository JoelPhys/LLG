import numpy as np 
import matplotlib.pyplot as plt
import sys
import math

data_lda1 = np.loadtxt('/Users/Hirst/Documents/PhD/LLG_code/SimpleCrystal_3D/Mag_v_temp.txt', usecols =(0,1,2,3,4))
T_vals_lda1 = data_lda1[:,0]
X_vals_lda1 = data_lda1[:,1]
Y_vals_lda1 = data_lda1[:,2]
Z_vals_lda1 = data_lda1[:,3]
M_vals_lda1 = data_lda1[:,4]

# data_lda2 = np.loadtxt('AFM/rsj-lda-af-p3-v2/af_lda_Mag_v_temp_1255to1750.txt', usecols =(0,1,2,3,4))
# T_vals_lda2 = data_lda2[:,0]
# X_vals_lda2 = data_lda2[:,1]
# Y_vals_lda2 = data_lda2[:,2]
# Z_vals_lda2 = data_lda2[:,3]
# M_vals_lda2 = data_lda2[:,4]

# data_gga1 = np.loadtxt('AFM/rsj-gga-af-p3-v2/af_gga_Mag_v_temp_0to1250.txt', usecols =(0,1,2,3,4))
# T_vals_gga1 = data_gga1[:,0]
# X_vals_gga1 = data_gga1[:,1]
# Y_vals_gga1 = data_gga1[:,2]
# Z_vals_gga1 = data_gga1[:,3]
# M_vals_gga1 = data_gga1[:,4]

# data_gga2 = np.loadtxt('AFM/rsj-gga-af-p3-v2/af_gga_Mag_v_temp_1255to1750.txt', usecols =(0,1,2,3,4))
# T_vals_gga2 = data_gga2[:,0]
# X_vals_gga2 = data_gga2[:,1]
# Y_vals_gga2 = data_gga2[:,2]
# Z_vals_gga2 = data_gga2[:,3]
# M_vals_gga2 = data_gga2[:,4]

data_lda_sus1 = np.loadtxt('/Users/Hirst/Documents/PhD/LLG_code/SimpleCrystal_3D/Susceptability.txt', usecols =(0,1,2))
testingT1_lda = data_lda_sus1[:,0]
TransX1_lda = data_lda_sus1[:,1]
ParraX1_lda = data_lda_sus1[:,2]

# data_lda_sus2 = np.loadtxt('AFM/rsj-lda-af-p3-v2/Susceptability-1255to1750.txt', usecols =(0,1,2))
# testingT2_lda = data_lda_sus2[:,0]
# TransX2_lda = data_lda_sus2[:,1]
# ParraX2_lda = data_lda_sus2[:,2]

# data_gga_sus1 = np.loadtxt('AFM/rsj-gga-af-p3-v2/Susceptability-0to1250.txt', usecols =(0,1,2))
# testingT1_gga = data_gga_sus1[:,0]
# TransX1_gga = data_gga_sus1[:,1]
# ParraX1_gga = data_gga_sus1[:,2]

# data_gga_sus2 = np.loadtxt('AFM/rsj-gga-af-p3-v2/Susceptability-1255to1750.txt', usecols =(0,1,2))
# testingT2_gga = data_gga_sus2[:,0]
# TransX2_gga = data_gga_sus2[:,1]
# ParraX2_gga = data_gga_sus2[:,2]


plt.figure(1)
plt.grid(True)
# plt.scatter(T_vals_gga1, M_vals_gga1, facecolors='none', edgecolors='b', label='GGA')
# plt.scatter(T_vals_gga2, M_vals_gga2, facecolors='none', edgecolors='b', label='GGA')
plt.scatter(T_vals_lda1, M_vals_lda1, facecolors='none', edgecolors='r', marker='s', label='LDA')
# plt.scatter(T_vals_lda2, M_vals_lda2, facecolors='none', edgecolors='r', marker='s', label='LDA')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('$M / M_s$')
plt.title('Magnetisation v Temperature AFM configurations')


fig, (ax1, ax2) = plt.subplots(2)
ax1.grid(True)
ax2.grid(True)
fig.tight_layout()
ax1.scatter(testingT1_lda, ParraX1_lda, facecolors='none', edgecolors='r', marker='s')
ax1.scatter(testingT1_lda, TransX1_lda, facecolors='none', edgecolors='b')
# ax1.scatter(testingT2_lda, ParraX2_lda, facecolors='none', edgecolors='r', marker='s')
# ax1.scatter(testingT2_lda, TransX2_lda, facecolors='none', edgecolors='b')
# ax1.set_yscale('log')
ax1.set_title('LDA', y=1.0,pad=-14)
ax1.legend([r"$ \tilde \chi_{⊥}$",r"$ \tilde \chi_{||}$"])
ax1.set_ylabel(r'$ \tilde \chi$ (1 / T)')
# ax2.set_title("GGA")
# ax2.scatter(testingT1_gga, ParraX1_gga, facecolors='none', edgecolors='r', marker='s')
# ax2.scatter(testingT1_gga, TransX1_gga, facecolors='none', edgecolors='b')
# ax2.scatter(testingT2_gga, ParraX2_gga, facecolors='none', edgecolors='r', marker='s')
# ax2.scatter(testingT2_gga, TransX2_gga, facecolors='none', edgecolors='b')
# ax2.legend([r"$ \tilde \chi_{||}$",r"$ \tilde \chi_{⊥}$"])
# ax2.set_xlabel('Temp (K)')
# ax2.set_ylabel(r'$ \tilde \chi$ (1 / T)')
plt.show()