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
	return (4*math.pi)*((1/(p[0]*2*math.pi))**(1.5))*x*x*np.exp(-(x**2)/p[0]/2)
# read data from a text file. One number per line
filename = sys.argv[1]
numOfBin = float(60)
data = np.loadtxt(filename,dtype='float64')
#Density normalizza l'istogramma 
hist , bin_edges = np.histogram(data, bins= numOfBin, density=True)
#  init dei parametri del fit:


step = float((bin_edges[-1] - bin_edges[0])/(2*numOfBin))
x_data = bin_edges[:-1]+step
guess=[1]
popt , pcov  = curve_fit(fitfunc, x_data, hist, p0=guess )

print(popt)
print(pcov)
#print("x è "+ str(len(x_data)) + "y è "+ str(len(hist)))
# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)
y = fitfunc(x_data,*popt)
fig = plt.figure()
fig.suptitle(r'Distribuzione di probabilità per la velocità')
ax = fig.add_subplot(111)
ax.grid(True)
width = 0.7 * (bin_edges[1] - bin_edges[0])
center = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(center, hist, align='center', width=width,color='g')
#plot
plt.plot(x_data,y,'r--',linewidth=2)
plt.text(2.6,0.5, r'$ T = %lf\pm %lf$'%(popt[0],math.sqrt(pcov[0][0])),bbox={'facecolor':'green', 'alpha':0.5, 'pad':10} )
plt.xlabel(r'$|v|$')
plt.ylabel(r'P($|v|$)')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.show()