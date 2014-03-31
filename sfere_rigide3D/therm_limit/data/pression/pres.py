#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from numpy import vectorize
from scipy.optimize import curve_fit
import math

def retta(x,*p):
	return p[0]*x + p[1]


filename1 = 'press-N.dat'
x,y,ers = np.loadtxt(filename1, comments='#',unpack=True)
#temp = sorted(temp, key=itemgetter(0))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
ersnp= np.asarray(ers,dtype='float64')
guess=[1,1]
popt,pcov = curve_fit(retta,xnp,ynp,sigma=ersnp,p0=guess)
print(popt)
print(pcov)
x= np.linspace(0,1/32,100)
plt.suptitle(r'Pressione a $\eta=%.1lf$, limite termodinamico'%(0.3),fontsize=18)
plt.grid(True)
plt.xlabel(r'$\frac{1}{N}$',fontsize=18)
plt.ylabel(r'$\frac{P V}{n k_b T} - 1$', fontsize=18)
plt.errorbar(xnp,ynp,yerr=ersnp,fmt='+',label='Dati')
plt.plot(x,retta(x,*popt), label='Fit')
plt.text(1/50,5.98,r'Fit con: $f(x)=m x + q$''\n' r'$m=%lf \pm %lf$' '\n' r'$q=%lf \pm %lf$'%(popt[0],math.sqrt(pcov[0][0]),popt[1],math.sqrt(pcov[1][1])),bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})

plt.show()
