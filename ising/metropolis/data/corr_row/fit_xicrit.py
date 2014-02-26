#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob

np.seterr(all='warn')
BetaC = 0.43

files = glob.glob("corr_row*.dat")

print (files)



temp = np.loadtxt(filename,dtype='float64')
		#= np.loadtxt('en_')
xs= []
ys= []
for i in range(len(temp)):
	t = temp[i]
	if 0.1<t[0]< BetaC:
		xs.append(t[0])
		ys.append(t[1])
x_points = np.asarray(xs,dtype='float64')
y_points = np.asarray(ys,dtype='float64')

def fit_fun(x,*p):
	return p[0]*(np.power((-1)*(x+(-p[2]))/p[2],p[1]))

#x_points = (-1)*(x_points + (-BetaC))/BetaC
guess = [1,-1,BetaC]
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
print(popt)
print(pcov)
#AErr = math.sqrt(pcov[0][0])
#tauErr = math.sqrt(pcov[1][1])
x = np.linspace(0.01,0.43,100)
AErr = 1
tauErr =1
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