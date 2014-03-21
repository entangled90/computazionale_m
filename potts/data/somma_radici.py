#!/usr/bin/python

from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
from pylab import *
import numpy as np
import itertools as itools



L=30
N=L*L
n_max=3.0
x=[]
y=[]
l = itools.combinations_with_replacement('012', N )
for perm in l:
	p = np.asarray(perm, dtype='int16')
	real = np.sum( np.cos(p*2*np.pi/n_max))/N
	im = np.sum( np.sin(p*2*np.pi/n_max))/N
	x.append(real)
	y.append(im)

numOfBin=300
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