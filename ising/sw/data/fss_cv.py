#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os

np.seterr(all='warn')

def par(x,*p):
	return p[0]*x**2 + p[1]*x + p[2] 
xmin=0.43
xmax=0.445

def fit_par(filename):
	guess = [-100,100,-5]
	x_points,y_points = load_file(filename,0.35,0.5)
	popt,pcov = curve_fit(par,x_points,y_points, p0=guess )
	print(popt)
	print(pcov)
	print("Beta critico è %lf"%(-popt[1]/(2.0*popt[0])))
	return ( -popt[1]/(2.0*popt[0]),popt)


def load_file(filename, betamin, betamax):
	lungh=20
	temp = np.loadtxt(filename,dtype='float64')
	xs= []
	ys= []
	for i in range(len(temp)):
		t = temp[i]
		if betamin<t[0]< betamax:
			xs.append(t[0])
			ys.append(t[1])
	x_points = np.asarray(xs,dtype='float64')
	y_points = np.asarray(ys,dtype='float64')
	return (x_points,y_points)


N = int(sys.argv[1])
nu_corr = 1.0
f = "cv%s.dat"%(N)
data = np.loadtxt(f,dtype='float64')
data = sorted(data, key=itemgetter(0))
#fit_exp(f)
BETA_CRIT,popt = fit_par(f)
xs= []
ys= []
ers=[]


for i in range(len(data)):
	t = data[i]
	xs.append(t[0])
	ys.append(t[1])
	ers.append(t[2])

x = np.asarray(xs,dtype='float64')
y = np.asarray(ys,dtype='float64')
er = np.asarray(ers,dtype='float64')

del xs[:]
del ys[:]
del ers[:]

plt.errorbar(x,y,yerr=er,fmt='ro')
x_lin= np.linspace(xmin,xmax,100)
plt.plot(x_lin,par(x_lin,*popt))
plt.show()

x = (x + -BETA_CRIT)/x
x *= N**(1/nu_corr)
y /= math.log(N)
er/= math.log(N)

out_file = open('./cvfssN%d.dat'%(N),"w")
for i in range(len(x)):
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(x[i],y[i],er[i]))
out_file.close()
