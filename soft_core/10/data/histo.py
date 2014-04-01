#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import glob

def force (r):
	t=1/r;
	return (24*(2*(t)**13 -(t)**7) -0.039);

f = sys.argv[1]
data=np.loadtxt(f,dtype='float64')
numOfBin=100
print("File caricato")
hist,bin_edges=np.histogram(data,bins=numOfBin, density=True)
hist*=100
step = (bin_edges[0]+bin_edges[-1])/(2*numOfBin)
x=bin_edges[:-1]+step
y=force(data)
press_vera= y*data/108/1.22
press_hist=force(x)*x*hist/108/1.22
valore_pres=press_vera.mean()
print("Valore vero:%lf"%(valore_pres))
valore_pres=press_hist.mean()

print("Valore hist:%lf"%(valore_pres))

plt.grid()
plt.plot(x,force(x),'-',label='force')
plt.plot(x,hist,label='hist')
plt.plot(x,press_hist,label='press')
plt.legend()
plt.show()