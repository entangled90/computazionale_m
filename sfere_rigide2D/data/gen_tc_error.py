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



N = int(sys.argv[1])

files=[]

file_temp = glob.glob('%s/tc.dat'%(N))	
for f in file_temp:
	if (os.stat(f)[6]!= 0):
		files.append(f)
b_num = 5
out_file = open('%d/press-eta.dat'%(N),"a")
i=0
for f in files:
	print(f[-12:-4])
	eta_temp = float(f[-12:-4])
	press_tmp = np.loadtxt(f, dtype='float64')
	print(press_tmp)
	out_file.write('%.8lf\t%.14e\t%.14e\n'%(eta_temp,press_tmp.mean(),press_tmp.std()/math.sqrt(len(press_tmp))))
	os.rename(f,f[:4]+'checked'+f[4:])
out_file.close()
