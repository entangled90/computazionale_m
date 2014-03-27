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



eta,tau=np.loadtxt('pdf_tau.dat',unpack=True)
fig = plt.figure()
fig.suptitle(r'$\tau_{fit}$ in funzione di $\eta$')
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
#midplot = np.amin(tau)/2.0  + np.amax(tau)/2.0
plt.xlabel(r'$\eta$',fontsize='15')
plt.ylabel(r'$\tau_{fit}$',fontsize='15')
ax.text(0.55,0.07,r'Funzione di fit:' '\n' r' $P(t_{coll}) = A e^{-\frac{t_{coll}}{\tau}}$', fontsize='15',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(eta, tau,fmt='o',label=r'$\tau$ da fit')
#ax.plot(x,fit_fun_corr_fitted(x,*pPlot),label ='fitted curve')
plt.legend(loc='upper center')
plt.show()