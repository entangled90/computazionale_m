#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
f='pressione.dat'
L_bin = 1000
t,dati = np.loadtxt(f,unpack=True)
N_bin = int(len(dati)/L_bin)
dati_binnati = [dati[i*L_bin:(i+1)*L_bin].mean() for i in range(N_bin)]
dati_np = np.asarray(dati_binnati,dtype='float64')
print ("media: %lf +- %lf"%(dati_np.mean(), dati_np.std()/math.sqrt(N_bin)))

plt.plot(dati_np)
plt.show()