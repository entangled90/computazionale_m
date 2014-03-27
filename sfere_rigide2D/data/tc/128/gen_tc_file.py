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



files=[]

file_temp = glob.glob('tc*.dat')
for f in file_temp:
	if (os.stat(f)[6]!= 0):
		files.append(f)

out_file = open('eta-tc.dat',"a")
for f in files:
	print(f[3:-4])
	eta_temp = float(f[3:-4])
	eta , tc_tmp = np.loadtxt(f, dtype='float64',unpack=True)
	print(tc_tmp)
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(eta_temp,tc_tmp.mean(),tc_tmp.std()/math.sqrt(len(tc_tmp))))
	os.rename(f,'checked'+f)
out_file.close()

