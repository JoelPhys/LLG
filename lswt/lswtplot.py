import matplotlib.pyplot as plt 
import numpy as np 



data1 = np.loadtxt("sc_fm_lswt.txt", usecols =(0,1))
data2 = np.loadtxt("sc_af_lswt.txt", usecols =(0,1))

k1 = data1[:,0]
freq1 = data1[:,1]

k2 = data2[:,0]
freq2 = data2[:,1]

fig, (ax1, ax2) = plt.subplots(2)
ax1.grid(True)
ax2.grid(True)
fig.tight_layout()
ax1.scatter(k1, freq1, facecolors='none', edgecolors='r', marker='s')
ax1.set_title('FM', y=1.0,pad=-14)
ax1.set_ylabel('Frequency (Hz)')
ax2.set_title("AFM")
ax2.scatter(k2, freq2, facecolors='none', edgecolors='r', marker='s')
ax2.set_xlabel('wave vector (k)')
ax2.set_ylabel('Frequency (Hz)')
plt.show()