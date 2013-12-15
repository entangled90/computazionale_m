#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]
x = []
y = []
for line in open(filename):
	values = line.strip('\n').split()
	x.append(float(values[0]))
	y.append(float(values[1]))


plt.plot(x,y,'b')
plt.show()
