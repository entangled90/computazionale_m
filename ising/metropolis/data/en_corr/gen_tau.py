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
BetaC = 0.4406868
BetaMin = 0.3
N = sys.argv[1]
files = glob.glob('en_autocorrN%s_*.dat'%(N))	
lungh_max = 40 #prima che si appiattica al punto critico
tau_temp = 0.0
beta = []
tau =  []
l_ord = []
for f in files:
	if (os.stat(f)[6]!= 0):
		beta_temp= float(f[-13:-4])
		print(beta_temp)
		temp = np.loadtxt(f,dtype='float64')
			#= np.loadtxt('en_')
		xs= []
		ys= []
		for i in range(len(temp)):
			t = temp[i]
			if t[0]<lungh_max:
				xs.append(t[0])
				ys.append(t[1])

		tau_temp = 0.5
		for i in range(len(ys)):
			tau_temp+= ys[i]
		l_ord.append([beta_temp,tau_temp])

l_ord = sorted(l_ord, key=itemgetter(0))
for l in l_ord:
	beta.append(l[0])
	tau.append(l[1])
#	error.append(l[2])

beta_np = np.asarray(beta,dtype='float64')
tau_np = np.asarray(tau,dtype='float64')
#error_np = np.asarray(error,dtype='float64')
#print (beta_np)
out_file = open('tau_corrN%s.dat'%(N),"w")
for i in range(len(beta)):
	out_file.write('%.14e\t%.14e\n'%(beta_np[i],tau_np[i]))
out_file.close()
beta_np = (beta_np + (- 0.4406868))/beta_np
x = np.linspace(BetaMin,np.amax(beta_np),100)
fig = plt.figure()
fig.suptitle('Tempo di Correlazione N=%s'%(N))
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
plt.xlabel(r'$\beta$',fontsize='15')
plt.ylabel(r'$\tau_{corr}$',fontsize='15')
#ax.text(BetaMin+0.05,20,r'Funzione di fit: $C(t) = A \tau^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$' '\n' r'$\; \beta_c=%lf\pm %lf$'%(-1/popt[1],ExpErr,popt[2],BetaErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.plot(beta_np, tau_np ,'o',label='original data')
#ax.plot(x,fit_fun_corr(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()
