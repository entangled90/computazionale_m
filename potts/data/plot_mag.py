#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter


filename = sys.argv[1]
x,y,ers = np.loadtxt(filename, comments='#',unpack=True)

fig = plt.figure()
fig.suptitle("Magnetizzazione modello di Potts")
plt.xlabel(r'$\beta$')
plt.ylabel(r'<|M|>')
plt.errorbar(x,y,yerr=ers,fmt='|')
plt.show()
