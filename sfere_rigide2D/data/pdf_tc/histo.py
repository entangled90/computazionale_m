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
	return p[0]*exp(-p[1]*x)
# read data from a text file. One number per line
filename = sys.argv[1]
eta = float (filename[ :7])
numOfBin = float(100)
data = np.loadtxt(filename,dtype='float64')
#Density normalizza l'istogramma 
hist , bin_edges = np.histogram(data, bins= numOfBin, density=False)
#  init dei parametri del fit:


step = float((bin_edges[-1] - bin_edges[0])/numOfBin)
x_data = bin_edges[:-1]+step
out_file = open('hist_pdftc%.8lf.dat'%(eta),"w")


for i in range( len(hist) ):
	out_file.write('%.8lf\t%.14e\n'%(x_data[i],hist[i]))
out_file.close()

guess = [ 100,-1]
popt,pcov  = curve_fit(fitfunc, x_data, hist, p0=guess )

print(popt)
print(pcov)
#tau = 1/popt[1]
tau=100
#print("x è "+ str(len(x_data)) + "y è "+ str(len(hist)))
# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)
#y = fitfunc(coeff,x_data)
fig = plt.figure()
fig.suptitle(r'Pdf dei tempi di collisione per $\eta = %lf$'%(eta))
ax = fig.add_subplot(111)
ax.grid(True)
width = 0.7 * (bin_edges[1] - bin_edges[0])
center = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(center, hist, align='center', width=width,color='g')
#plot
#plt.plot(x_data,y,'r--',linewidth=2)
plt.text(0.0008,7000,'Funzione di fit: ' '\n' r'$P(t) = e^{-\frac{t}{\tau}}$'  '\n' r'$ \tau = %.4e $'%(tau))
plt.xlabel('Time collision')
plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)
plt.show()