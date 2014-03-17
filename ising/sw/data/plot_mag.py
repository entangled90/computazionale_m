#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

BETAC = 0.4406868
x=[]
y=[]
filename = sys.argv[1]
temp = np.loadtxt(filename, comments='#')
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
for vector in temp:
	x.append(float(vector[0]))
	y.append(float(abs(vector[1])))
xnp = np.asarray(x,dtype='float64')
ynp = np.asarray(y,dtype='float64')
xnp = (xnp + (- BETAC))/xnp
plt.plot(xnp,ynp,'ro')
plt.show()
