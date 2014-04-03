#!/usr/bin/python
import numpy as np
import sys


f = sys.argv[1]
t,dr2,ers = np.loadtxt(f,unpack=True)

out_file=open(f,"w")
for i in range(len(t)):
	out_file.write("%e\t%e\t%e\n"%(t[i]/100,dr2[i],ers[i]))

out_file.close()