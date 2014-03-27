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


file_list = sys.argv[1:]
fig = plt.figure()
plt.suptitle(r'$\Delta r^2(t)$ in funzione di $\eta$',fontsize=17)
ax = fig.add_subplot(111)
ax.grid(True)
plt.xlabel(r'$t$',fontsize='15')
plt.ylabel(r'$\Delta r^2(t)$',fontsize='15')
for f in file_list:
	eta_temp=float(f[:8])
	eta,mfp,ers=np.loadtxt(f,unpack=True)
	ax.errorbar(eta, mfp,yerr=ers,fmt='+',label=r'$\eta=%.5lf$'%(eta_temp))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
plt.show()