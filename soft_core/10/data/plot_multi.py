#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.special  as sc


x=[]
y=[]
plt.xlabel(r'$t$',fontsize='18')
#plt.ylabel(r'$<\frac{U}{N}>$',fontsize='18')
plt.suptitle(r'Confronto tra varie osservabili',fontsize='18')

filename = sys.argv[1:]
i=0
labels=['Energia interna', 'Pressione','Temperatura']
for f in filename:
	x,y = np.loadtxt(f, comments='#',unpack=True)
	#Ordina i valori x,y secondo la x:
	#temp = sorted(temp, key=itemgetter(0))
	xnp = np.asarray(x,dtype='float64')
	ynp = np.asarray(y,dtype='float64')
	plt.errorbar(xnp,ynp, fmt='+',label=labels[i])
	i+=1
plt.legend(loc='center right')

plt.show()
