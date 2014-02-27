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
def fit_fun_corr(x,*p):
	return p[0]*(-1*(x + -1*p[2])/p[2])**p[1]


np.seterr(all='warn')
BetaC = 0.44
BetaMin = 0.3
N = sys.argv[1]
file_temp = glob.glob('corr_row_N%s*.dat'%(N))	
beta = []
xi =  []
files = []
l_ord = []
for f in file_temp:
	beta_temp= float(f[-9:-4])
#	print(beta_temp) 
	if BetaMin <beta_temp < BetaC:
		files.append(f)

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
#print (beta_np)
out_file = open('xi_corrN%s.dat'%(N),"w")
for i in range(len(beta)):
	out_file.write('%.14e\t%.14e\t%.14e\n'%(beta_np[i],xi_np[i],error[i]))
out_file.close()


guess = [1, -1 , 0.44]
popt,pcov = curve_fit(fit_fun_corr,beta_np,xi_np, sigma=error, p0=guess )

print(popt,pcov)
#beta_np = (beta_np + -1*popt[2])/popt[2]
#print (beta_np)
Aerr = math.sqrt(pcov[0][0])
ExpErr= math.sqrt(pcov[1][1])
BetaErr = math.sqrt(pcov[2][2]) 
x = np.linspace(np.amin(beta_np),np.amax(beta_np),100)
fig = plt.figure()
fig.suptitle("Lunghezza di Correlazione")
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
plt.xlabel(r'$\beta$',fontsize='15')
plt.ylabel(r'$\xi_{corr}$',fontsize='15')
ax.text(BetaMin,20,r'Funzione di fit: $C(t) = A \xi^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$' '\n' r'$\; \beta_c=%lf\pm %lf$'%(-1/popt[1],ExpErr,popt[2],BetaErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(beta_np, xi_np ,yerr=error_np, fmt='o',ecolor='r',label='original data')
ax.plot(x,fit_fun_corr(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()