#!/usr/bin/python

from subprocess import call
import time
import sys
import datetime 
import multiprocessing
import os
import numpy as np 
import random 


def loop (eta_values,index):
	start = time.time()
	f = open("log/core" + index + ".log", "w")
	for l in eta_values:
		print ("Starting program with particles=" + str(l)+ " with core ",str(index))
		call(["./main",str(l)],stdout=f)
		print ("Program with ",str(l)," particles finished")
	print("time elapsed for core"+index+": ", str(datetime.timedelta(seconds=(time.time() - start))))


iteration = int(sys.argv[2])
n_core = int(sys.argv[1])
# 32 son le iterazione fra 16 e 512 con passo di 16
step = int ((512 - 32)/ iteration)
iteration_per_core = int ( iteration/n_core)
allvalues =[]
allvalues.extend(range(32,512,step))
temp=[]
temp_index= 0
eta_cores=[]
for i in range(n_core):
	temp = []
	for j in range(iteration_per_core):
		temp_index = int(random.random()*len(allvalues))
		temp.append(allvalues[temp_index])
		allvalues.pop(temp_index)
	temp = sorted(temp)
	eta_cores.append(temp)

jobs=[]
for i in range(n_core):
	arg = [eta_cores[i],str(i)]
	p = multiprocessing.Process(target=loop, args=arg)
	jobs.append(p)
	p.start()
print("Script exited successfully")
