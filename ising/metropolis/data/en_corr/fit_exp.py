#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os


def fit_fun(x,*p):
	return p[0]*(np.exp(x*p[1]))

np.seterr(all='log')

f = sys.argv[1]
beta = []
tau =  []
l_ord = []
temp = np.loadtxt(f,dtype='float64')
	#= np.loadtxt('en_')
xs= []
ys= []
for i in range(len(temp)):
	t = temp[i]
	if t[0]<7:
		xs.append(t[0])
		ys.append(t[1])
x_points = np.asarray(xs,dtype='float64')
y_points = np.asarray(ys,dtype='float64')
#x_points = (-1)*(x_points + (-BetaC))/BetaC
guess = [1,-0.01]
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
print('l\'expcritico è : %f'%(-1/popt[1]))
x = np.linspace(np.amin(x_points),np.amax(x_points),100)
fig = plt.figure()
fig.suptitle('Tempo di Correlazione (magnetizzazione)')
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
plt.xlabel(r'$\beta$',fontsize='15')
plt.ylabel(r'$\tau_{corr}$',fontsize='15')
#ax.text(BetaMin+0.05,20,r'Funzione di fit: $C(t) = A \tau^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$' '\n' r'$\; \beta_c=%lf\pm %lf$'%(-1/popt[1],ExpErr,popt[2],BetaErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(x_points, y_points ,fmt='o', ecolor='r',label='original data')
ax.plot(x,fit_fun(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()