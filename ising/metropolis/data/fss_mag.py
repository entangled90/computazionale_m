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
xmin=0.43
xmax=0.44
def par(x,*p):
	return p[0]*x**2 + p[1]*x + p[2] 

def fit_par(filename):
	guess = [-100,100,0]
	x_points,y_points = load_file(filename,xmin,xmax)
	popt,pcov = curve_fit(par,x_points,y_points, p0=guess )
	print(popt)
	print(pcov)
	print("Beta critico è %lf"%(-popt[1]/(2.0*popt[0])))
	return ( -popt[1]/(2.0*popt[0]),popt)

def load_file(filename, betamin, betamax):
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
mag_file = "mag_mean%s.dat"%(N)
mag_data = np.loadtxt(mag_file,dtype='float64',usecols=(0,1,2))
mag_data = sorted(mag_data, key=itemgetter(0))
gamma = 0.125
#fit_exp(mag_file)
BETA_CRIT = 0.44
#,popt= fit_par(mag_file)
xs= []
ys= []
ers=[]
for i in range(len(mag_data)):
	t = mag_data[i]
	xs.append(t[0])
	ys.append(t[1])
	ers.append(t[2])
x = np.asarray(xs,dtype='float64')
y = np.asarray(ys,dtype='float64')
er = np.asarray(ers,dtype='float64')
del xs[:]
del ys[:]
del ers[:]

#x_lin = np.linspace(xmin,xmax,100)
#plt.plot(x_lin,par(x_lin,*popt))
plt.errorbar(x,y,yerr=er)
plt.show()
x = (x + -BETA_CRIT)/x
x *= N**(1/nu_corr)
y /= N**(-gamma/nu_corr)
er/= N**(-gamma/nu_corr)

out_file = open('./magfssN%d.dat'%(N),"w")
for i in range(len(x)):
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(x[i],y[i],er[i]))
out_file.close()


