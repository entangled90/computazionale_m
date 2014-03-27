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

#Definisce la funzione per il fit
def fitfunc(x,*p): 
	return p[0]*np.exp(-p[1]*x)
# read data from a text file. One number per line
filename = sys.argv[1]
eta = float (filename[ :7])
numOfBin = float(100)
data = np.loadtxt(filename,dtype='float64')
print('La media è:%lf'%(data.mean()))

#Density normalizza l'istogramma 
hist , bin_edges = np.histogram(data, bins= numOfBin, density=False)
#  init dei parametri del fit:


step = float((bin_edges[-1] - bin_edges[0])/numOfBin)
x_data = bin_edges[:-1]+step
out_file = open('hist_pdftc%.8lf.dat'%(eta),"w")
print(x_data)
print(hist)
for i in range( len(hist) ):
	out_file.write('%.8lf\t%.14e\n'%(x_data[i],hist[i]))
out_file.close()

guess = [ 100, 1]
popt,pcov  = curve_fit(fitfunc, x_data, hist,p0=guess)

print(popt)
print(pcov)
tau = 1/popt[1]
#print("x è "+ str(len(x_data)) + "y è "+ str(len(hist)))
# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)
y = fitfunc(x_data,*popt)
fig = plt.figure()
fig.suptitle(r'Pdf dei tempi di collisione per $\eta = %lf$'%(eta))
ax = fig.add_subplot(111)
ax.grid(True)
width = 0.7 * (bin_edges[1] - bin_edges[0])
center = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(center, hist, align='center', width=width,color='g')
#plot
plt.plot(x_data,y,'r--',linewidth=2)
ax.annotate ('Funzione di fit: ' '\n' r'$P(t) = e^{-\frac{t}{\tau}}$'  '\n' r'$ \tau = %.4e $'%(tau), 
	xy=(1,1),xycoords='axes fraction', ha='right',va='top', xytext=(-10, -10), textcoords='offset points',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
plt.xlabel('Time collision')
plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.show()