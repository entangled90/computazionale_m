#!/usr/bin/python

from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
from pylab import *
from scipy.optimize import curve_fit
import numpy as np

#Definisce la funzione per il fit
def fitfunc(x,*p):
	return p[2]*exp(-(x-p[0])**2/p[1])
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
# read data from a text file. One number per line
filename = sys.argv[1]
numOfBin = 30
data = np.loadtxt(filename, dtype='float64')
#Density normalizza l'istogramma 
hist , bin_edges = np.histogram(data, bins= numOfBin, density=False)
#  init dei parametri del fit:
step = (float(bin_edges[-1] - bin_edges[0])/numOfBin)
x_data = bin_edges[:-1]+step
#print(hist)
guess = [0,1,1]
#print( hist)
#popt,pcov  = curve_fit(fitfunc ,x_data, hist, p0=guess)
#print(popt)
#print(pcov)
#print("x è "+ str(len(x_data)) + "y è "+ str(len(hist)))
# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)
#x= np.linspace(bin_edges[0],bin_edges[-1],100)
fig=plt.figure()
width = 0.7 * (bin_edges[1] - bin_edges[0])
center = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(center, hist, align='center', width=width,color='g')
#plt.plot(x,fitfunc(x,*popt), 'r')
fig.suptitle(r'Distribuzione di probabilità a $\beta = 1$')
plt.xlabel(r'Magnetizzazione $M$')
plt.ylabel(r'Densità di probabilità (normalizzata)')
#plt.xlabel('Time collision')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.show()