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
xmin=0.35
xmax=0.5	

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
data = np.loadtxt(f,dtype='float64',usecols=(0,1))
data = sorted(data, key=itemgetter(0))
#fit_exp(f)
BETA_CRIT,popt = fit_par(f)
xs= []
ys= []
for i in range(len(data)):
	t = data[i]
	xs.append(t[0])
	ys.append(t[1])
x_chi = np.asarray(xs,dtype='float64')
y_chi = np.asarray(ys,dtype='float64')

del xs[:]
del ys[:]
plt.plot(x_chi,y_chi,'ro')
x = np.linspace(xmin,xmax,100)
plt.plot(x,par(x,*popt))
plt.show()

x_chi = (x_chi + -BETA_CRIT)/x_chi
x_chi *= N**(1/nu_corr)
y_chi /= math.log(N)

out_file = open('./cvfssN%d.dat'%(N),"w")
for i in range(len(x_chi)):
	out_file.write('%.8lf\t%.14e\n'%(x_chi[i],y_chi[i]))
out_file.close()
