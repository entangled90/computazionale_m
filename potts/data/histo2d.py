#!/usr/bin/python

from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
from pylab import *
from scipy.optimize import curve_fit
import numpy as np


filename = sys.argv[1]
numOfBin = 140
x,y = np.loadtxt(filename, dtype='float64', unpack=True)
#Density normalizza l'istogramma 
H,xedges,yedges= np.histogram2d(x,y, bins= numOfBin)
fig=plt.figure()
H = np.rot90(H)
H = np.flipud(H)

#plt.plot(x,fitfunc(x,*popt), 'r')
Hmasked = np.ma.masked_where(H==0,H)
fig.suptitle(r'Distribuzione di probabilità')
plt.xlabel(r'Parte Reale $M$')
plt.ylabel(r'Parte Immaginaria')
#plt.xlabel('Time collision')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.pcolormesh(xedges,yedges,Hmasked)

cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')

plt.grid(True)
plt.show()