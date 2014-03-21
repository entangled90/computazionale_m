#!/usr/bin/python

import sys
import numpy as np
import glob
import os


N = sys.argv[1]

file_temp = glob.glob('corr_row_N%s*.dat'%(N))	
files=[]
for f in file_temp:
	if (os.stat(f)[6]!= 0):
		files.append(f)


for f in files:
	x,y,yer= np.loadtxt(f,unpack=True,dtype='float64')
	yer_np = np.asarray(yer)
	yer_np /= 10.0

	out_file = open('ok'+str(f),"w")
	for i in range(len(x)):
		out_file.write('%.8lf\t%.14e\t%.14e\n'%(x[i],y[i],yer_np[i]))
	out_file.close()
