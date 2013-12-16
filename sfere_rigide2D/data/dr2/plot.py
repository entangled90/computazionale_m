#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
x=[]
y=[]
filename = sys.argv[1]
#Carica il file fottendosene di "#"
temp = np.loadtxt(filename, comments='#')
#Ordina i valori x,y secondo la x:
sorted(temp, key=itemgetter(0))
for vector in temp:
	x.append(float(vector[0]))
	y.append(float(vector[1]))
plt.plot(x,y,'b')
plt.show()