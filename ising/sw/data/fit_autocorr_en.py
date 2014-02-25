#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys

filename = sys.argv[1]
np.seterr(all='warn')
lungh=20
y_points = np.loadtxt(filename,dtype='float64',usecols=[1])
		#= np.loadtxt('en_')
y_points = y_points[ :lungh]
x_points = range(len(y_points))
def fit_fun(x,*p):
	return p[0]*np.exp(x*p[1])

x = np.linspace(0,lungh,100)
guess = [1,-5]
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )

AErr = math.sqrt(pcov[0][0])
tauErr = math.sqrt(pcov[1][1])

fig = plt.figure()
fig.suptitle("Autocorrelazione dell'energia")
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
ax.text(10.2,0.4,r'Funzione di fit: $C(t) = A e^{-t/ \tau }$' '\n' r'$A=%lf \pm %lf$' '\n' r'$\; \tau=%lf\pm %lf$'%(popt[0],AErr,-1/popt[1],tauErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.plot(x_points,y_points,'ko',label='Original data')
plt.plot(x,fit_fun(x,*popt),label ='fitted curve')
plt.legend()
plt.show()