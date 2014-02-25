#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys


lungh=20

filename = sys.argv[1]


y_points = np.loadtxt(filename,usecols=[1])
		#= np.loadtxt('en_')
y_points = y_points[ :lungh]
x_points = range(len(y_points))
def fit_fun(x,*p):
	return p[0]*np.exp(x*p[1])
x = np.linspace(0,lungh,100)
guess = [1,-5]
popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )

print(popt)
print(pcov)

fig = plt.figure()
fig.suptitle('Autocorrelazione della magnetizzazione')
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
ax.text(10.2,0.4,r'Funzione di fit: $f(x) = A e^{-x/ \tau }$' '\n' r'$A=%lf \; \tau=%lf$'%(popt[0],-1/popt[1]), bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.plot(x_points,y_points,'ko',label='Original data')
plt.plot(x,fit_fun(x,*popt),label ='fitted curve')
plt.legend()
plt.show()