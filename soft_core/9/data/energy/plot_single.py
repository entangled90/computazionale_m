#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.special  as sc
import math

x=[]
y=[]
ers=[]
L_bin=2000

filename = sys.argv[1]
N = int(filename[-7:-4])
x,y = np.loadtxt(filename, comments='#',unpack=True)

N_bin = int(len(y)/L_bin)
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
dati_binnati = [y[i*L_bin:(i+1)*L_bin].mean() for i in range(N_bin)]
bin_np = np.asarray(dati_binnati,dtype='float64')
print(bin_np.mean())
print(bin_np.std()/math.sqrt(len(bin_np)))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
plt.suptitle(r'Energia interna $U$ in funzione del tempo',fontsize='18')
plt.xlabel(r'$t$',fontsize='18')
plt.ylabel(r'$<\frac{U}{N}>$',fontsize='18')
plt.errorbar(xnp,ynp, fmt='+',label='%d'%(N))
plt.legend(loc='upper right')

plt.show()
