#!/usr/bin/python

from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
from pylab import *
from scipy.optimize import curve_fit
import numpy as np



f = 'press-eta.dat'

eta, p, ers = np.loadtxt(f,unpack=True)

p = eta*p

plt.errorbar(eta,p,fmt='+',yerr=ers)
plt.show()