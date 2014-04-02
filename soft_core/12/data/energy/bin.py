#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import glob 
import os



L_bin = 1000
files = glob.glob('energy*.dat')	
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
	out_file.write('%d\t%.14e\t%.14e\n'%(n[i],m[i],d[i]))
out_file.close()

fig = plt.figure()
fig.suptitle("Energia interna in funzione del numero di particelle")
plt.grid()
plt.xlabel("Numero di particelle")
plt.ylabel(r'$\frac{U}{N \epsilon}$')
plt.errorbar(n,m,fmt='+',yerr=dev)
plt.show()