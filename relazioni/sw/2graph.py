#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys


file1= sys.argv[1]
file2 = sys.argv[2]

x1,y1,err1 = np.loadtxt(file1, unpack=True)

x2,y2,err2 = np.loadtxt(file2, unpack=True)


fig = plt.figure()
fig.suptitle('Tempo di autocorrelazione')
ax = fig.add_subplot(111)
ax.grid(True)
#parameter_plot = [popt[0],popt[1],0]
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
plt.xlabel(r'$\beta$',fontsize='15')
plt.ylabel(r'$\tau_{corr}$',fontsize='15')
#ax.text(BetaMin+0.05,20,r'Funzione di fit: $C(t) = A \tau^{-\nu}$' '\n' r'$\nu=%lf \pm %lf$' '\n' r'$\; \beta_c=%lf\pm %lf$'%(-1/popt[1],ExpErr,popt[2],BetaErr) , bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
ax.errorbar(x1,y1,yerr=err1, fmt='o' ,label='Metropolis')
ax.errorbar(x2,y2,yerr=err2, fmt='o' , label="Swendsen-Wang")
#ax.plot(x,fit_fun_corr(x,*popt),label ='fitted curve')
plt.legend(loc='upper left')
plt.show()