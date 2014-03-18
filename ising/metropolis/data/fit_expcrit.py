#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys

filename = sys.argv[1]
np.seterr(all='warn')
lungh=20
BetaC = 0.41
temp = np.loadtxt(filename,dtype='float64')
		#= np.loadtxt('en_')
xs= []
ys= []
ers=[]
for i in range(len(temp)):
	t = temp[i]
	if 0.2<t[0]< BetaC:
		xs.append(t[0])
		ys.append(t[1])
		ers.append(t[2])
x_points = np.asarray(xs,dtype='float64')
y_points = np.asarray(ys,dtype='float64')
ers_points = np.asarray(ers,dtype='float64')
def fit_fun(x,*p):
	return p[0]*(np.power(np.fabs((x+(-p[2]))/p[2]),p[1]))

#x_points = (-1)*(x_points + (-BetaC))/BetaC
guess = [1,-1,0.44]
x = np.linspace(0.2,0.415,100)
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
print(popt)
print(pcov)
#AErr = math.sqrt(pcov[0][0])
gammaerr = math.sqrt(pcov[1][1])
fig = plt.figure()
fig.suptitle(r'Esponente critico di $\chi$ ')
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
ax.text(0.25,20,r'Funzione di fit: $C(t) = A t^{-\gamma}$' '\n' r'$\; \gamma=%lf\pm %lf$'%(-popt[1],gammaerr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.errorbar(x_points,y_points,yerr=ers_points, fmt='o',label='Original data')
plt.plot(x,fit_fun(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()