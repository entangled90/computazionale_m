#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os
import math
import random


def fitfunc(x,*p): 
	return p[0]*np.exp(-p[1]*x)
numOfBin=100
np.seterr(all='warn')
file_temp = glob.glob('0.*.dat')	
#print(file_temp)
eta = []
tau =  []
files = []
for f in file_temp:
	eta_temp= float(f[:8])
#	print(eta_temp) 
	if(os.stat(f)[6]!= 0):
		files.append(f)
l_ord=[]
for f in files:
	eta_temp= float(f[:8])
	data = np.loadtxt(f,dtype='float64')
	hist , bin_edges = np.histogram(data, bins= numOfBin, density=True)
	step = float((bin_edges[-1] - bin_edges[0])/(2*numOfBin))
	x_data = bin_edges[:-1]+step
	guess = [ 100,0.001]
	popt,pcov = curve_fit(fitfunc,x_data,hist, p0=guess )
	eta.append(eta_temp)
	tau.append(1/popt[1]) #, math.sqrt(pcov[1][1])])


out_file = open('pdf_tau.dat',"a")
for i in range(len(eta)):
	out_file.write('%.8lf\t%.14e\n'%(eta[i],tau[i]))
out_file.close()

fig = plt.figure()
fig.suptitle(r'$\tau_{coll}$ in funzione di $\eta$')
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
#midplot = np.amin(tau)/2.0  + np.amax(tau)/2.0
plt.xlabel(r'$\eta$',fontsize='15')
plt.ylabel(r'$\tau_{coll}$',fontsize='15')
#ax.text(np.amin(eta)+0.0001,midplot,r'Funzione di fit: $C(t) = A \xi^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$'%(-popt[1],ExpErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(eta, tau,fmt='|',label='original data')
#ax.plot(x,fit_fun_corr_fitted(x,*pPlot),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()