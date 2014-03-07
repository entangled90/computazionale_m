#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os


def fit_fun(x,*p):
	return p[0]*(np.exp(x*p[1]))




np.seterr(all='warn')
N = sys.argv[1]
file_temp = glob.glob('corr_row_N%s*.dat'%(N))	
beta = []
xi =  []
files = []
l_ord = []
for f in file_temp:
	beta_temp= float(f[-13:-4])
#	print(beta_temp) 
	if  (os.stat(f)[6]!= 0):
		files.append(f)
for f in files:
	beta_temp= float(f[-13:-4])
	temp = np.loadtxt(f,dtype='float64')
		#= np.loadtxt('en_')
	xs= []
	ys= []
	for i in range(len(temp)):
		t = temp[i]
		if t[0]<int(N)/5:
			xs.append(t[0])
			ys.append(t[1])
	x_points = np.asarray(xs,dtype='float64')
	y_points = np.asarray(ys,dtype='float64')
	#x_points = (-1)*(x_points + (-BetaC))/BetaC
	guess = [1,-1]
	popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
	if not (any( math.isnan(i) for i in popt) or any(math.isinf(i) for i in popt)):
		l_ord.append([beta_temp,-1/popt[1], math.sqrt(pcov[1][1])])
error = []
l_ord = sorted(l_ord, key=itemgetter(0))
for l in l_ord:
	beta.append(l[0])
	xi.append(l[1])
	error.append(l[2])
beta_np = np.asarray(beta,dtype='float64')
xi_np = np.asarray(xi,dtype='float64')
error_np = np.asarray(error,dtype='float64')


out_file = open('../xi_corrN%s.dat'%(N),"w")
for i in range(len(beta)):
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(beta_np[i],xi_np[i],error[i]))
out_file.close()
