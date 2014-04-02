#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.special  as sc


x=[]
y=[]
plt.xlabel(r'$t$',fontsize='18')
plt.ylabel(r'$<\frac{U}{N}>$',fontsize='18')
plt.suptitle(r'Energia interna $U$ in funzione del tempo',fontsize='18')

filename = sys.argv[1:]
for f in filename:
	N = int(f[-7:-4])
	x,y = np.loadtxt(f, comments='#',unpack=True)
	#Ordina i valori x,y secondo la x:
	#temp = sorted(temp, key=itemgetter(0))
	xnp = np.asarray(x,dtype='float64')
	ynp = np.asarray(y,dtype='float64')
	plt.errorbar(xnp,ynp, fmt='+',label='%d'%(N))
plt.legend(loc='upper right')

plt.show()
