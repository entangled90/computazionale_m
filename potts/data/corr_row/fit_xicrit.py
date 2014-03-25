#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os
import math


def fit_fun(x,*p):
	return p[0]*(np.exp(x*p[1]))
def fit_fun_corr(x,*p):
#	print(x)
	return p[0]*(np.fabs((x + -1*BETAES)/BETAES))**p[1]

def fit_fun_corr_fitted(x,*p):
	return p[0]*(np.fabs(x))**p[1]


np.seterr(all='warn')
BetaC = 0.995
BetaMin = 0.98
BETAES= math.log(1+math.sqrt(3))-0.0001
N = sys.argv[1]
file_temp = glob.glob('corr_row_N%s*.dat'%(N))	
#print(file_temp)
beta = []
xi =  []
files = []
l_ord = []
for f in file_temp:
	beta_temp= float(f[-14:-4])
	print(beta_temp) 
	if (BetaMin <beta_temp < BetaC) & (os.stat(f)[6]!= 0):
		files.append(f)

for f in files:
	beta_temp= float(f[-14:-4])
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


guess = [1, -0.83]
popt,pcov = curve_fit(fit_fun_corr,beta_np,xi_np, sigma=error, p0=guess )
file_list = glob.glob('../xi_corrN%s_NU*.dat'%(N))
for f in file_list:
	os.remove(f)
out_file = open('../xi_corrN%s_NU%.8lf.dat'%(N,-popt[1]),"w")
for i in range(len(beta)):
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(beta_np[i],xi_np[i],error[i]))
out_file.close()

print(beta_np)
print(popt,pcov)
#RISCALO BETA!
beta_np = (beta_np + -1*BETAES)/beta_np

#print (beta_np)
Aerr = math.sqrt(pcov[0][0])
ExpErr= math.sqrt(pcov[1][1])
#BetaErr = math.sqrt(pcov[2][2])
x = np.linspace(np.amin(beta_np),np.amax(beta_np),100)
fig = plt.figure()
fig.suptitle('Lunghezza di Correlazione N=%s'%(N))
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
midplot = np.amin(xi_np)/2.0  + np.amax(xi_np)/2.0
pPlot=popt[:2]
plt.xlabel(r'$\beta$',fontsize='15')
plt.ylabel(r'$\xi_{corr}$',fontsize='15')
ax.text(np.amin(beta_np)+0.0001,midplot,r'Funzione di fit: $C(t) = A \xi^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$'%(-popt[1],ExpErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(beta_np, xi_np ,yerr=error_np,fmt='|',ecolor='r',label='original data')
ax.plot(x,fit_fun_corr_fitted(x,*pPlot),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()