#!/usr/bin/python

from subprocess import call
import time
import sys
import datetime 
import multiprocessing
import os
import numpy as np 
import random 

####
# INPUT: x , x_max , numero di iterazioni

def loop (iter,index):
	start = time.time()
	f = open("log/core" + index + ".log", "w")
	for l in iter:
		print ("  Starting ising sw with Beta*J="+ str(l[0])+' ' + str(l[1]))
		command = []
		command.append("./ising")
		for i in l:
			command.append(i)
		call(command,stdout=f ,shell=False)
	print("time elapsed for core"+index+": ", str(datetime.timedelta(seconds=(time.time() - start))))

x=float(sys.argv[1])
x_max=float(sys.argv[2])
n_iteration=int(sys.argv[3])
n_lattice = int(sys.argv[4])
n_core = int(sys.argv[5])
iteration_per_core= int(n_iteration/n_core)
step = (x_max - x)/n_iteration
#Valori di Beta
allvalues = np.arange(x,x_max,step).tolist()
if ( n_iteration % n_core != 0):
	print("Iterazioni multiple del numero di core selezionati(",n_core,")!")
	sys.exit(1)
eta_cores = []
temp_index=0
for i in range(n_core):
	temp = []
	for j in range(iteration_per_core):
		temp_index = int(random.random()*len(allvalues))
		temp.append([str(allvalues[temp_index]),str(n_lattice)])
		allvalues.pop(temp_index)
	eta_cores.append(temp)
jobs=[]
for i in range(n_core):
	arg = [eta_cores[i],str(i)]
	p = multiprocessing.Process(target=loop, args=arg)
	jobs.append(p)
	p.start()
print("Script exited successfully")
