#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys


def fit_fun(x,*p):
	return p[0]*(np.power(np.absolute((x+(-p[2]))/x),p[1]))



filename = sys.argv[1]
np.seterr(all='warn')
lungh=20
BetaC = 0.43
BetaMin = 0.3
temp = np.loadtxt(filename,dtype='float64')
		#= np.loadtxt('en_')
xs= []
ys= []
for i in range(len(temp)):
	t = temp[i]
	if BetaMin<t[0]< BetaC:
		xs.append(t[0])
		ys.append(t[1])
x_points = np.asarray(xs,dtype='float64')
y_points = np.asarray(ys,dtype='float64')

#x_points = (-1)*(x_points + (-BetaC))/BetaC
guess = [1,-1,BetaC]
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
print(popt)
print(pcov)
AErr = math.sqrt(pcov[0][0])
tauErr = math.sqrt(pcov[1][1])
x = np.linspace(BetaMin,BetaC,100)
#AErr = 1
#tauErr =1
fig = plt.figure()
fig.suptitle("Esponente critico di %s"%(filename))
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
ax.text(BetaMin,np.amax(y_points)/2.0+np.amin(y_points)/2.0,r'Funzione di fit: $C(t) = A x^{-\gamma} }$' '\n' r'$A=%lf \pm %lf$' '\n' r'$\; \gamma=%lf\pm %lf$' ' \n' r'$\beta = %lf $'%(popt[0],AErr,-popt[1],tauErr,popt[2]) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.plot(x_points,y_points,'ko',label='Original data')
plt.plot(x,fit_fun(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()