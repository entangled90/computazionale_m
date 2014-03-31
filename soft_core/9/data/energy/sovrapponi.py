#!/usr/bin/python


import math
import numpy as np
import matplotlib.pyplot as plt



f1='energy108.dat'
f2='energy508.dat'


x1,y1 = np.loadtxt(f1,unpack=True)
x2,y2 = np.loadtxt(f2,unpack=True)




y1_scarto = y1- y1.mean()
y2_scarto = y2-y2.mean()


plt.plot(x1,y1_scarto,label='108')
plt.plot(x1,y2_scarto,label='508')
plt.show()