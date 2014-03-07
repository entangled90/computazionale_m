#!/usr/bin/python

from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
from pylab import *
from scipy.optimize import leastsq
import numpy as np

#Definisce la funzione per il fit
fitfunc  = lambda p, x: p[0]*exp(-p[1]*x)
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
# read data from a text file. One number per line
filename = sys.argv[1]
numOfBin = float(20)
data = []
for item in open(filename,'r'):
    item = item.strip()
    if item != '':
        try:
            data.append(float(item))
        except ValueError:
            pass

#Density normalizza l'istogramma 
hist , bin_edges = np.histogram(data, bins= numOfBin, density=True)
#  init dei parametri del fit:
step = float((bin_edges[-1] - bin_edges[0])/numOfBin)
x_data = bin_edges[:-1]+step
#print(hist)
init = [1,0.1]
#out   = leastsq( errfunc, init, args=(x_data, hist))
#coeff_fit = out[0]
#print("x è "+ str(len(x_data)) + "y è "+ str(len(hist)))
# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)
#y = fitfunc(coeff_fit,x_data)
#l1=plt.plot(x_data,y,'r--',linewidth=2)


width = 0.7 * (bin_edges[1] - bin_edges[0])
center = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(center, hist, align='center', width=width,color='g')
#plot
plt.xlabel('mfp')
plt.ylabel('Occurrence')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.show()