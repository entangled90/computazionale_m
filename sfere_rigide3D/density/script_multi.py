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
		print ("Starting sfere2D with eta=" + str(l))
		call(["./main",str(l)],stdout=f)
	print("time elapsed for core"+index+": ", str(datetime.timedelta(seconds=(time.time() - start))))

eta=float(sys.argv[1])
eta_max=float(sys.argv[2])
n_iteration=int(sys.argv[3])
n_core = int(sys.argv[4])
iteration_per_core= int(n_iteration/n_core)
step = (eta_max - eta)/n_iteration
allvalues = np.arange(eta,eta_max,step).tolist()
if ( n_iteration % n_core != 0):
	print("Iterazioni multiple del numero di core selezionati(",n_core,")!")
	sys.exit(1)
temp=[]
temp_index= 0
eta_cores=[]
for i in range(n_core):
	temp = []
	for j in range(iteration_per_core):
		temp_index = int(random.random()*len(allvalues))
		temp.append(allvalues[temp_index])
		allvalues.pop(temp_index)
	eta_cores.append(temp)

jobs=[]
for i in range(n_core):
	arg = [eta_cores[i],str(i)]
	p = multiprocessing.Process(target=loop, args=arg)
	jobs.append(p)
	p.start()
print("Script exited successfully")
