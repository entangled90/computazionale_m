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

file_temp = glob.glob('pression*.dat')	
for f in file_temp:
	if (os.stat(f)[6]!= 0):
		files.append(f)
out_file = open('press-N.dat',"a")
i=0
for f in files:
	print(f[-7:-4])
	N = int(f[-7:-4])
	N_vec,press_tmp = np.loadtxt(f, dtype='float64',unpack=True)
	print(press_tmp)
	out_file.write('%lf\t%.14e\t%.14e\n'%(float(1/float(N)),press_tmp.mean(),press_tmp.std()/math.sqrt(len(press_tmp))))
	os.rename(f,'checked'+f)
out_file.close()

