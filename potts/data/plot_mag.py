#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

x=[]
y=[]
filename = sys.argv[1]
temp = np.loadtxt(filename, comments='#')
#Ordina i valori x,y secondo la x:
#temp = sorted(temp, key=itemgetter(0))
for vector in temp:
	x.append(float(vector[0]))
	y.append(float(abs(vector[1])))
plt.plot(x,y,'ro')
plt.show()