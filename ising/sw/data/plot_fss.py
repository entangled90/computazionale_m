#!/usr/bin/python

import sys
import glob
from operator import itemgetter
import os
import numpy as np
import os.path
import matplotlib.pyplot as plt
N_MAX = 200

osservabile = str(sys.argv[1]) # chi o cv per es
file_list=[]
for n in range(N_MAX):
	if os.path.isfile(osservabile+"fssN%d.dat"%(n)):
		file_list.append([n,osservabile+"fssN%d.dat"%(n)])

all_datas=[]
for i in file_list:
	f = i[1]
	print(f)
	x,y,err= np.loadtxt(f,dtype='float64',unpack=True)
	f = str(f)
	N= i[0]
	all_datas.append([x,y,err,N])
x_min = -10
x_max = 10
fig = plt.figure()
fig.suptitle("Studio di finite size scaling di %s"%(osservabile))
ax = fig.add_subplot(111)
ax.grid(True)
#Mette le griglie su asse x
#plt.xticks([i for i in range(0,lungh)])
ax.text(x_min,10, 'Legenda',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
for d in all_datas:
	plt.errorbar(d[0],d[1], fmt='|',yerr=d[2],label ='N = %d'%(int(d[3])))
plt.xlabel('Lx')
plt.legend(loc='upper left')
plt.show()
