#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys


f=sys.argv[1]
L_bin = 2500
t,dati = np.loadtxt(f,unpack=True)
N_bin = int(len(dati)/L_bin)
dati_binnati = [dati[i*L_bin:(i+1)*L_bin].mean() for i in range(N_bin)]
dati_np = np.asarray(dati_binnati,dtype='float64')
print ("media: %lf +- %lf"%(dati_np.mean(), dati_np.std()/math.sqrt(N_bin)))
plt.suptitle('Temperatura in funzione del tempo')
plt.ylabel(r'$T$',fontsize=18)
plt.xlabel("t",fontsize=18)
plt.plot(dati,'+')
print (len(dati))
plt.show()