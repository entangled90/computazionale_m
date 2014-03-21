#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from numpy import vectorize


def mag_fun(x):
	BETAC = 0.4406868
	result=[]
	for i in x:
		if i > BETAC:
			result.append((1 - np.sinh(2*i)**(-4))**(0.125))
		else:
			result.append(0)
	return np.asarray(result,dtype='float64')
vmag_fun = vectorize(mag_fun)
x=[]
y=[]
filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3]
x,y,ers = np.loadtxt(filename1, comments='#',unpack=True)
x2,y2,ers2 = np.loadtxt(filename2, comments='#',unpack=True)
x3,y3,ers3 = np.loadtxt(filename3, comments='#',unpack=True)
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
ersnp= np.asarray(ers,dtype='float64')
xnp2 = np.asarray(x2,dtype='float64')
ynp2 = np.asarray(y2,dtype='float64')
ersnp2= np.asarray(ers2,dtype='float64')
xnp3 = np.asarray(x3,dtype='float64')
ynp3 = np.asarray(y3,dtype='float64')
ersnp3= np.asarray(ers3,dtype='float64')
x_lin = np.linspace(0,1,100)
plt.suptitle('Magnetizzazione in funzione di Beta')
plt.xlabel(r'$\beta$')
plt.ylabel('<|M|>')
plt.plot(x_lin,mag_fun(x_lin),label='Onsager')
plt.errorbar(xnp,ynp,yerr=ersnp,fmt='|', label='50')
plt.errorbar(xnp2,ynp2,yerr=ersnp2,fmt='|', label='75')
plt.errorbar(xnp3,ynp3,yerr=ersnp3,fmt='|', label='100')
plt.legend(loc='upper left')

plt.show()
