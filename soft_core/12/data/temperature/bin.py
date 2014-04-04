#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import glob 
import os
from scipy.optimize import curve_fit


def fit(x,*p):
	return p[0]*x + p[1]

L_bin = 20000
files = glob.glob('temperature*.dat')	
n=[]
m=[]
d=[]
for f in files:
	N=int(f[-8:-4])
	print(N)
	t,dati = np.loadtxt(f,unpack=True)
	N_bin = int(len(dati)/L_bin)
	dati_binnati = [dati[i*L_bin:(i+1)*L_bin].mean() for i in range(N_bin)]
	dati_np = np.asarray(dati_binnati,dtype='float64')
	media=dati_np.mean()
	dev=dati_np.std()/math.sqrt(N_bin)
	m.append(media)
	n.append(N)
	d.append(dev)
	print ("%d --- media: %lf +- %lf"%(N,media,dev ))
print("fine ciclo")
print(n[2])
out_file = open('binned_energy',"w")
for i in range(len(n)):
	out_file.write('%d\t%.14e\t%.14e\n'%(1/n[i],m[i],d[i]))
out_file.close()

for i in range(len(n)):
	n[i] = 1/n[i]

m_np=np.asarray(m)
n_np=np.asarray(n)
guess=[0,-3]
popt,pcov = curve_fit(fit,n_np,m_np,p0=guess)
print(popt , pcov)
x=np.linspace(0,0.012,1000)

fig = plt.figure()
fig.suptitle("Temperatura in funzione del numero di particelle")
plt.grid()
plt.xlabel(r'$\frac{1}{N}$')
plt.ylabel(r'$<T>$')
plt.errorbar(n,m,fmt='+',yerr=dev)
plt.plot(x,fit(x,*popt))

plt.show()