#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.special  as sc

def en_fun(x):
	return - 1/(np.tanh(2*x))*( 1 + 2/np.pi*(2*(np.tanh(2*x)**2)-1)*sc.ellipk((2*np.sinh(2*x))**2/(np.cosh(2*x)**4)))



BETAC = 0.4406868
x=[]
y=[]
ers=[]
filename = sys.argv[1]
x,y,ers = np.loadtxt(filename, comments='#',unpack=True)
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
ersnp= np.asarray(ers,dtype='float64')
x_lin = np.linspace(0,1,100)
plt.suptitle('Energia in funzione di Beta')
plt.xlabel(r'$\beta$')
plt.ylabel('<H>')
plt.plot(x_lin,en_fun(x_lin),label='Onsager')
plt.errorbar(xnp,ynp,yerr=ersnp, fmt='o',label='dati')
plt.legend(loc='upper right')

plt.show()
