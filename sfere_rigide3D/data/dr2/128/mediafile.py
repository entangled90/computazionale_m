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


file_list = glob.glob('dr2_*_*.dat')
temps=[]
for f in file_list:
	tmp=float(f[-12:-4])
	if tmp not in temps:
		temps.append(tmp)

print (temps)
for temp in temps:
	time_np= np.zeros(148,dtype='float64')
	mean = np.zeros(148,dtype='float64')
	var = np.zeros(148,dtype='float64')
	print('dr2_*_%.6lf.dat'%(temp))
	file_temp = glob.glob('dr2_*_%.6lf.dat'%(temp))
	print(file_temp)
	for f in file_temp:
		time , dr2 = np.loadtxt(f,dtype='float64',unpack=True,usecols=(0,1))
		mean += dr2
		var += dr2*dr2
		time_np = time

	mean /= (len(file_temp))
	var /= (len(file_temp))
	var -= mean*mean
	var = np.sqrt(var/(len(file_temp)))
	np.savetxt('mean_dr2_%1.3lf.dat'%(temp), np.c_[time_np,mean,var], fmt='%1.14e')
