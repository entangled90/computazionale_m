#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os
import math
import random
import matplotlib.cm as cm
import matplotlib.colors as colors


NUM_COLORS = 2

cm = plt.get_cmap('Paired')

file_list = sys.argv[1:]
fig = plt.figure()
#plt.suptitle(r'Pdf per i tempi di collisione in funzione di $\eta$',fontsize=16)
ax = fig.add_subplot(111)
ax.grid(True)
plt.xlabel(r'$t_c$',fontsize='15')
plt.ylabel(r'$P(t_c)$',fontsize='15')
ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
i=0
Markers=['x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.']
col=['r','b']
#labels for cfr between fit and mean
labels=['fit','calcolato']
for f in file_list:
	#eta_temp=float(f[-13:-4])
	eta,mfp=np.loadtxt(f,unpack=True,usecols=(0,1))
	width = 0.7 * (eta[1] - eta[0])
	color = cm(1.*i/NUM_COLORS) 
	#plt.bar(eta, mfp, width=width,label=r'$\eta=%.5lf$'%(eta_temp), alpha=(0.8),color=color,linewidth=0)
	plt.plot(eta,mfp,Markers[i], label=labels[i],color=col[i])
	i+=1
	mean=0

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='upper right') #bbox_to_anchor=(1, 0.5), fontsize=12)
plt.show()