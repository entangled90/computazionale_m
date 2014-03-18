#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from numpy import vectorize


def mag_fun(x):
	BETAC = 0.4406868
	if x >= BETAC:
		return 1 - 1/((np.sinh(2*x))**4)
	else:
		return 0

vmag_func = vectorize(mag_fun)
x=[]
y=[]
filename = sys.argv[1]
x,y,ers = np.loadtxt(filename, comments='#',unpack=True)
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
ersnp= np.asarray(ers,dtype='float64')
x_lin = np.linspace(0,1,100)
plt.suptitle('Magnetizzazione in funzione di Beta')
plt.xlabel(r'$\beta$')
plt.ylabel('<|M|>')
plt.plot(x_lin,vmag_func(x_lin))
plt.errorbar(xnp,ynp,yerr=ersnp, fmt='o')
plt.show()
