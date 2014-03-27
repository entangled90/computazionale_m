#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from numpy import vectorize


filename1 = sys.argv[1]
x,y,ers = np.loadtxt(filename1, comments='#',unpack=True)
#temp = sorted(temp, key=itemgetter(0))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
ersnp= np.asarray(ers,dtype='float64')
plt.suptitle(r'Tempo medio di collisione in funzione di $\eta$',fontsize=18)
plt.grid(True)
plt.xlabel(r'$\eta$',fontsize=18)
plt.ylabel(r'$t_c$', fontsize=18)
plt.errorbar(xnp,ynp,yerr=ersnp,fmt='|')
plt.show()
