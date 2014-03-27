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



eta,mfp,ers=np.loadtxt('eta-mfp.dat',unpack=True)
fig = plt.figure()
plt.suptitle(r'Libero cammino medio in funzione di $\eta$',fontsize=17)
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
#midplot = np.amin(mfp)/2.0  + np.amax(mfp)/2.0
plt.xlabel(r'$\eta$',fontsize='15')
plt.ylabel(r'$l_c$',fontsize='15')
#ax.text(0.55,0.0017,r'Funzione di fit:' '\n' r' $P(t_{coll}) = A e^{-\frac{t_{coll}}{\mfp}}$', fontsize='15',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(eta, mfp,yerr=ers,fmt='+',label=r'$l_c$')
#ax.plot(x,fit_fun_corr_fitted(x,*pPlot),label ='fitted curve')
plt.legend(loc='upper right')
plt.show()