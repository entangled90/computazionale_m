#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter


def fit_fun(x,*p):
	return p[0]*(np.exp(x*p[1]))

np.seterr(all='warn')
BetaC = 0.4
N = sys.argv[1]
files = glob.glob('corr_row_N%s*.dat'%(N))	
beta = []
xi =  []

l_ord = []
for f in files:
	beta_temp= float(f[-9:-4])
#	print(beta_temp) 
	if beta_temp > BetaC:
		files.remove(f)

for f in files:
	beta_temp= float(f[-9:-4])
	temp = np.loadtxt(f,dtype='float64')
		#= np.loadtxt('en_')
	xs= []
	ys= []
	for i in range(len(temp)):
		t = temp[i]
		if t[0]<20:
			xs.append(t[0])
			ys.append(t[1])
	x_points = np.asarray(xs,dtype='float64')
	y_points = np.asarray(ys,dtype='float64')
	#x_points = (-1)*(x_points + (-BetaC))/BetaC
	guess = [1,-1]
	popt,pcov = curve_fit(fit_fun,x_points,y_points, p0=guess )
	if not (any( math.isnan(i) for i in popt) or any(math.isinf(i) for i in popt)):
		l_ord.append([beta_temp,-1/popt[1]])

l_ord = sorted(l_ord, key=itemgetter(0))
for l in l_ord:
	beta.append(l[0])
	xi.append(l[1])
beta_np = np.asarray(beta,dtype='float64')
xi_np = np.asarray(xi,dtype='float64')
#print (beta_np)
beta_np = (beta_np + -1*BetaC)/BetaC
#print (beta_np)

#AErr = math.sqrt(pcov[0][0])
#tauErr = math.sqrt(pcov[1][1])
#x = np.linspace(np.amin(beta_np),np.amax(beta_np),100)
AErr = 1
tauErr =1
fig = plt.figure()
fig.suptitle("Lunghezza di Correlazione")
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
#ax.text(10.2,0.4,r'Funzione di fit: $C(t) = A e^{-t/ \tau }$' '\n' r'$A=%lf \pm %lf$' '\n' r'$\; \tau=%lf\pm %lf$'%(popt[0],AErr,-1/popt[1],tauErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.plot(beta_np, xi_np ,'ro',label='Original data')
plt.plot(beta,xi,label ='fitted curve')
plt.legend()
plt.show()